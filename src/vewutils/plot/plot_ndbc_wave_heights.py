"""Plot modeled significant wave heights against NDBC buoy observations.

Model values come from a SWAN global output file written by ADCIRC+SWAN
(``swan_HS.63.nc`` by default), interpolated to each buoy location inside its
containing mesh element, and observations come from NDBC via
:mod:`vewutils.plot.get_obswave`.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

DEFAULT_WAVE_FILE = 'swan_HS.63.nc'
DEFAULT_WAVE_VARIABLE = 'swan_HS'
DEFAULT_WAVE_MAP_HTML = '0000_wave_height_map.html'
DEFAULT_WAVE_FILENAME_PATTERN = '{index:04d}_NDBC_{station_id}_{name}.png'
# Values at or below this are SWAN/ADCIRC fill values (e.g. -99999 when dry).
FILL_THRESHOLD = -9999.0

# (n_node, n_elem, x0, y0, xN, yN) -> mesh lookup arrays, so the analysis
# lookback files (all on one mesh) are only searched once per process.
_mesh_memory: dict[tuple, dict[str, np.ndarray]] = {}


# ---------------------------------------------------------------------------
# Model extraction
# ---------------------------------------------------------------------------

def _mesh_key(nc) -> tuple:
    x = nc.variables['x']
    y = nc.variables['y']
    n_node = x.shape[0]
    n_elem = nc.variables['element'].shape[0]
    return (n_node, n_elem, float(x[0]), float(y[0]),
            float(x[n_node - 1]), float(y[n_node - 1]))


def _load_mesh(nc) -> dict[str, np.ndarray]:
    key = _mesh_key(nc)
    mesh = _mesh_memory.get(key)
    if mesh is None:
        x = np.asarray(nc.variables['x'][:], dtype=float)
        y = np.asarray(nc.variables['y'][:], dtype=float)
        elem = np.asarray(nc.variables['element'][:], dtype=np.int64) - 1
        ex = x[elem]
        ey = y[elem]
        mesh = {
            'x': x,
            'y': y,
            'element': elem,
            'exmin': ex.min(axis=1),
            'exmax': ex.max(axis=1),
            'eymin': ey.min(axis=1),
            'eymax': ey.max(axis=1),
        }
        _mesh_memory[key] = mesh
    return mesh


def locate_points_in_mesh(
        mesh: dict[str, np.ndarray],
        lon: float,
        lat: float) -> tuple[np.ndarray, np.ndarray] | None:
    """Return (3 node indices, 3 barycentric weights) for a point, or None.

    Element bounding boxes are prefiltered first, so this is cheap even on
    large meshes. Like :func:`vewutils.plot.get_f63wl_at.get_f63wl_at`, the
    first containing element wins.
    """
    candidates = np.nonzero(
        (mesh['exmin'] <= lon) & (lon <= mesh['exmax'])
        & (mesh['eymin'] <= lat) & (lat <= mesh['eymax'])
    )[0]
    x = mesh['x']
    y = mesh['y']
    for e in candidates:
        n1, n2, n3 = mesh['element'][e]
        x1, y1, x2, y2, x3, y3 = x[n1], y[n1], x[n2], y[n2], x[n3], y[n3]
        denominator = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
        if denominator == 0:
            continue
        a = ((y2 - y3) * (lon - x3) + (x3 - x2) * (lat - y3)) / denominator
        b = ((y3 - y1) * (lon - x3) + (x1 - x3) * (lat - y3)) / denominator
        c = 1.0 - a - b
        if a >= 0.0 and b >= 0.0 and c >= 0.0:
            return np.array([n1, n2, n3]), np.array([a, b, c])
    return None


def mesh_bbox(nc_path: str | Path) -> tuple[float, float, float, float]:
    """(xmin, xmax, ymin, ymax) of the mesh nodes in a *.63.nc file."""
    from netCDF4 import Dataset

    with Dataset(nc_path) as nc:
        mesh = _load_mesh(nc)
    return (float(mesh['x'].min()), float(mesh['x'].max()),
            float(mesh['y'].min()), float(mesh['y'].max()))


def _read_nodes_over_time(var, nodes: np.ndarray) -> np.ndarray:
    """Read var[:, nodes] as a float array, choosing a chunk-friendly strategy."""
    if var.shape[0] == 0:
        return np.empty((0, len(nodes)), dtype=float)
    chunking = var.chunking()
    node_chunked = (
        isinstance(chunking, list) and len(chunking) == 2
        and chunking[1] <= 10 and chunking[0] > 100
    )
    if node_chunked:
        data = var[:, nodes]
    elif isinstance(chunking, list):
        # Time-chunked (the usual ADCIRC NetCDF4 layout): every chunk holds a
        # whole time step, so read one time step at a time and pick the
        # nodes out, touching each chunk once for all buoys together.
        data = np.ma.vstack([var[i, :][nodes] for i in range(var.shape[0])])
    else:
        # Contiguous (classic NetCDF3): strided per-node reads, as in
        # get_f63wl_at.
        data = np.ma.column_stack([var[:, int(n)] for n in nodes])
    data = np.ma.filled(np.ma.asarray(data, dtype=float), np.nan)
    data[data <= FILL_THRESHOLD] = np.nan
    return data


def _extraction_cache_path(cache_dir, nc_path, variable, points) -> Path:
    stat = os.stat(nc_path)
    key = json.dumps([
        os.path.abspath(nc_path), stat.st_size, int(stat.st_mtime), variable,
        sorted((sid, round(lon, 6), round(lat, 6)) for sid, lon, lat in points),
    ])
    digest = hashlib.sha1(key.encode('utf-8')).hexdigest()[:20]
    return Path(cache_dir) / 'model_point_series' / f'{digest}.json'


def extract_series_at_points(
        nc_path: str | Path,
        variable: str,
        points: list[tuple[str, float, float]],
        *,
        cache_dir: str | Path | None = None) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Interpolate a (time, node) variable to several points at once.

    Parameters
    ----------
    points : list of (point_id, lon, lat)

    Returns
    -------
    times : numpy datetime64 array
    values : dict point_id -> float array (all NaN when the point is outside
        the mesh or the variable is missing there)

    When cache_dir is given, results are cached as JSON keyed by the file's
    path, size and mtime, so an analysis file reused by later forecast cycles
    is only read once.
    """
    from netCDF4 import Dataset

    from vewutils.plot.plot_f61_hydrographs import read_netcdf_time_series

    cache_path = None
    if cache_dir:
        cache_path = _extraction_cache_path(cache_dir, nc_path, variable, points)
        if cache_path.is_file():
            try:
                with open(cache_path, encoding='utf-8') as f:
                    cached = json.load(f)
                times = np.array(cached['time'], dtype='datetime64[ns]')
                values = {
                    pid: np.array([np.nan if v is None else v for v in vals], dtype=float)
                    for pid, vals in cached['values'].items()
                }
                return times, values
            except (OSError, ValueError, KeyError):
                pass

    times = read_netcdf_time_series(nc_path).to_numpy().astype('datetime64[ns]')
    values: dict[str, np.ndarray] = {}
    with Dataset(nc_path) as nc:
        if variable not in nc.variables:
            raise KeyError(f'{nc_path}: variable {variable!r} not found')
        mesh = _load_mesh(nc)
        located = {}
        for pid, lon, lat in points:
            hit = locate_points_in_mesh(mesh, lon, lat)
            if hit is None:
                values[pid] = np.full(len(times), np.nan)
            else:
                located[pid] = hit
        if located:
            nodes = np.unique(np.concatenate([n for n, _ in located.values()]))
            column = {int(n): i for i, n in enumerate(nodes)}
            data = _read_nodes_over_time(nc.variables[variable], nodes)
            for pid, (elem_nodes, weights) in located.items():
                cols = [column[int(n)] for n in elem_nodes]
                values[pid] = data[:, cols] @ weights

    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, 'w', encoding='utf-8') as f:
            json.dump({
                'file': os.path.abspath(nc_path),
                'variable': variable,
                'time': [str(t) for t in times],
                'values': {
                    pid: [None if not np.isfinite(v) else float(v) for v in vals]
                    for pid, vals in values.items()
                },
            }, f)
    return times, values


def concat_series(pieces: list[tuple[np.ndarray, np.ndarray]]) -> tuple[np.ndarray, np.ndarray]:
    """Concatenate (time, value) pieces in order; later files win on overlap."""
    pieces = [(t, v) for t, v in pieces if len(t)]
    if not pieces:
        return np.array([], dtype='datetime64[ns]'), np.array([], dtype=float)
    times = np.concatenate([t for t, _ in pieces])
    vals = np.concatenate([v for _, v in pieces])
    # Stable reverse so the last occurrence of a duplicated time is kept.
    rev_times = times[::-1]
    _, first_in_rev = np.unique(rev_times, return_index=True)
    keep = len(times) - 1 - first_in_rev
    keep.sort()
    return times[keep], vals[keep]


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def comparison_stats(
        obs_time: np.ndarray,
        obs_hs: np.ndarray,
        model_time: np.ndarray,
        model_hs: np.ndarray) -> dict[str, float] | None:
    """Bias, RMSE and N of the model interpolated to observation times."""
    valid = np.isfinite(model_hs)
    if len(obs_time) == 0 or valid.sum() < 2:
        return None
    mt = model_time[valid].astype('datetime64[s]').astype(np.int64)
    mv = model_hs[valid]
    ot = obs_time.astype('datetime64[s]').astype(np.int64)
    inside = (ot >= mt[0]) & (ot <= mt[-1])
    if not inside.any():
        return None
    interp = np.interp(ot[inside], mt, mv)
    diff = interp - obs_hs[inside]
    diff = diff[np.isfinite(diff)]
    if len(diff) == 0:
        return None
    return {
        'bias': float(diff.mean()),
        'rmse': float(np.sqrt((diff ** 2).mean())),
        'n': int(len(diff)),
    }


def plot_wave_height_at_buoy(
        fig, ax,
        station: dict[str, Any],
        obs_time: np.ndarray,
        obs_hs: np.ndarray,
        model_series: list[tuple[np.ndarray, np.ndarray]],
        date_start: datetime,
        date_end: datetime,
        *,
        labels: list[str] | None = None,
        colors: list[str] | None = None,
        nowcast_forecast_style: bool = False,
        plot_in_foot: bool = False) -> dict[str, float] | None:
    """Draw observed and modeled Hs at one buoy; return comparison stats."""
    from matplotlib.dates import DateFormatter

    from vewutils.plot.plot_hydrograph_at_station import _default_line_color

    scale = 3.28084 if plot_in_foot else 1.0

    if len(obs_time):
        ax.plot(obs_time, obs_hs * scale, '.', markersize=4,
                color='k', label='Obs.')

    n_series = len(model_series)
    for i, (mtime, mhs) in enumerate(model_series):
        if nowcast_forecast_style:
            color = 'b'
            linestyle = '--' if i == n_series - 1 and n_series > 1 else '-'
        else:
            color = colors[i] if colors and i < len(colors) else _default_line_color(i)
            linestyle = '-'
        if labels and i < len(labels):
            label = labels[i]
        elif nowcast_forecast_style and n_series > 1:
            label = 'Forecast' if i == n_series - 1 else 'Analysis'
        else:
            label = 'Model' if n_series == 1 else f'Model {i + 1}'
        ax.plot(mtime, mhs * scale, linestyle, color=color, label=label)

    all_time, all_hs = concat_series(model_series)
    stats = comparison_stats(obs_time, obs_hs, all_time, all_hs)
    if stats is not None:
        unit = 'ft' if plot_in_foot else 'm'
        ax.text(
            0.01, 0.98,
            f"Bias {stats['bias'] * scale:+.2f} {unit}   "
            f"RMSE {stats['rmse'] * scale:.2f} {unit}   N {stats['n']}",
            transform=ax.transAxes, va='top', ha='left', fontsize=9,
            bbox={'facecolor': 'white', 'alpha': 0.8, 'edgecolor': 'none'},
        )

    ax.xaxis.set_major_formatter(DateFormatter('%m-%d %H:%M'))
    ax.grid()
    ax.set_title(f"NDBC {station['station_id']}: {station['name']}")
    ax.set_ylabel(
        'Significant Wave Height (ft)' if plot_in_foot
        else 'Significant Wave Height (m)'
    )
    ax.set_xlim([date_start, date_end])
    ax.set_ylim(bottom=0)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc='upper right')
    fig.tight_layout()
    return stats


def plot_wave_heights_at_ndbc_buoys(
        model_files: list[list[str]],
        output_dir: str | Path,
        *,
        variable: str = DEFAULT_WAVE_VARIABLE,
        date_start: datetime | None = None,
        date_end: datetime | None = None,
        stations: list[dict[str, Any]] | None = None,
        station_ids: list[str] | None = None,
        bboxes: list[tuple[float, float, float, float]] | None = None,
        station_types=None,
        exclude_programs=None,
        cache_dir: str | Path | None = None,
        labels: list[str] | None = None,
        colors: list[str] | None = None,
        nowcast_forecast_style: bool = False,
        require_observations: bool = True,
        plot_in_foot: bool = False,
        filename_pattern: str = DEFAULT_WAVE_FILENAME_PATTERN,
        max_stations: int | None = None,
        skip_existing: bool = False,
        skip_on_error: bool = False,
        cache_only: bool = False,
        figsize: tuple[float, float] = (12.0, 5.0)) -> tuple[list[Path], list[dict[str, Any]]]:
    """Plot modeled vs NDBC-observed Hs for buoys inside the model mesh.

    Parameters
    ----------
    model_files : list of list of str
        One inner list per plotted model series; files in an inner list are
        concatenated in order (e.g. ``[[analysis files...], [forecast]]``).
    stations : list of dict, optional
        NDBC station list (see :func:`vewutils.plot.get_obswave.get_ndbc_active_stations`);
        fetched when omitted.
    station_ids : list of str, optional
        Explicit NDBC ids to plot, in order, instead of every buoy in the mesh.
    bboxes : list of (xmin, xmax, ymin, ymax), optional
        Keep only buoys inside any of these boxes (a union).
    require_observations : bool
        Skip buoys with no observed Hs in the window (default True), so
        stations that don't measure waves don't fill the report.
    cache_only : bool
        Read the NDBC station list and observations only from cache_dir
        (filled by ``vewutils plot refresh-obs-cache``); never download.

    Returns
    -------
    written, station_records
        Paths written this run, and per-station records (keys ``station``,
        ``index``, ``output_path``, ``status``) in the format
        :func:`vewutils.plot.plot_f61_hydrographs.generate_hydrograph_map_html`
        expects.
    """
    import matplotlib.pyplot as plt

    from vewutils.plot.get_obswave import (
        DEFAULT_NDBC_EXCLUDE_PROGRAMS,
        DEFAULT_NDBC_STATION_TYPES,
        NDBCRequestError,
        filter_ndbc_stations,
        get_ndbc_active_stations,
        get_obswave,
        obs_series_to_numpy,
    )
    from vewutils.plot.plot_f61_hydrographs import format_output_path

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    model_files = [list(series) for series in model_files if series]
    if not model_files:
        raise ValueError('No model wave files given')
    first_file = model_files[0][0]

    if stations is None:
        stations = get_ndbc_active_stations(cache_dir=cache_dir, cache_only=cache_only)
    candidates = filter_ndbc_stations(
        stations,
        bbox=mesh_bbox(first_file),
        station_ids=station_ids,
        station_types=(
            DEFAULT_NDBC_STATION_TYPES if station_types is None else station_types
        ),
        exclude_programs=(
            DEFAULT_NDBC_EXCLUDE_PROGRAMS if exclude_programs is None
            else exclude_programs
        ),
    )
    if bboxes:
        candidates = [
            st for st in candidates
            if any(
                xmin <= st['lon'] <= xmax and ymin <= st['lat'] <= ymax
                for xmin, xmax, ymin, ymax in bboxes
            )
        ]

    points = [(st['station_id'], st['lon'], st['lat']) for st in candidates]
    print(
        f'Extracting {variable} at {len(points)} NDBC station(s) from '
        f'{sum(len(s) for s in model_files)} model file(s)'
    )
    series_by_station: dict[str, list[tuple[np.ndarray, np.ndarray]]] = {
        st['station_id']: [] for st in candidates
    }
    for series_files in model_files:
        pieces: dict[str, list] = {sid: [] for sid, _, _ in points}
        for nc_path in series_files:
            if not points:
                break
            times, values = extract_series_at_points(
                nc_path, variable, points, cache_dir=cache_dir
            )
            for sid in pieces:
                pieces[sid].append((times, values[sid]))
        for sid in pieces:
            series_by_station[sid].append(concat_series(pieces[sid]))

    # A buoy is in the domain when the model has a value there at some time.
    in_domain = [
        st for st in candidates
        if any(np.isfinite(v).any() for _, v in series_by_station[st['station_id']])
    ]
    if max_stations is not None:
        in_domain = in_domain[:max_stations]
    print(f'{len(in_domain)} NDBC station(s) inside the wet model domain')

    all_times = np.concatenate([
        t for series in series_by_station.values() for t, _ in series if len(t)
    ]) if in_domain else np.array([], dtype='datetime64[ns]')
    if date_start is None:
        date_start = (
            all_times.min().astype('datetime64[s]').astype(datetime)
            if len(all_times) else None
        )
    if date_end is None:
        date_end = (
            all_times.max().astype('datetime64[s]').astype(datetime)
            if len(all_times) else None
        )

    written: list[Path] = []
    records: list[dict[str, Any]] = []
    for index, station in enumerate(in_domain, start=1):
        sid = station['station_id']
        output_path = format_output_path(
            filename_pattern, index, 'NDBC', sid, station['name'], output_dir,
        )
        # generate_hydrograph_map_html colors markers by owner label.
        record_station = {**station, 'owner': 'NDBC'}
        if skip_existing and output_path.is_file():
            print(f'Skipping existing wave height figure: {output_path.name}')
            records.append({
                'station': record_station, 'index': index,
                'output_path': output_path, 'status': 'skipped_existing',
            })
            continue

        try:
            _, _, _, obs_time, obs_hs = get_obswave(
                sid, date_start, date_end, station=station, cache_dir=cache_dir,
                cache_only=cache_only,
            )
        except NDBCRequestError as exc:
            print(f'Warning: NDBC {sid}: {exc}', file=sys.stderr)
            obs_time = obs_hs = None
        obs_time, obs_hs = obs_series_to_numpy(obs_time, obs_hs)
        if require_observations and len(obs_time) == 0:
            print(f'No observed wave heights at NDBC {sid} ({station["name"]}); skipping')
            continue

        print(f'Processing NDBC {sid} ({station["name"]}) -> {output_path.name}')
        try:
            fig, ax = plt.subplots(figsize=figsize)
            plot_wave_height_at_buoy(
                fig, ax, station, obs_time, obs_hs,
                series_by_station[sid], date_start, date_end,
                labels=labels, colors=colors,
                nowcast_forecast_style=nowcast_forecast_style,
                plot_in_foot=plot_in_foot,
            )
            fig.savefig(output_path)
            plt.close(fig)
            written.append(output_path)
            records.append({
                'station': record_station, 'index': index,
                'output_path': output_path, 'status': 'plotted',
            })
        except Exception as exc:
            plt.close('all')
            if not skip_on_error:
                raise
            print(f'Warning: skipped NDBC station {sid}: {exc}', file=sys.stderr)
            records.append({
                'station': record_station, 'index': index,
                'output_path': output_path, 'status': 'failed',
            })

    return written, records
