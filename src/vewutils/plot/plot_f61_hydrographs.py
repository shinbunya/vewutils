"""Plot hydrographs for stations listed in an ADCIRC elev_stat.151 file."""

from __future__ import annotations

import os
import re
import sys
from datetime import datetime
from glob import glob
from pathlib import Path
from typing import Any

# elev_stat.151 owner label -> observation source for :func:`plot_hydrograph_at_station`
ELEV_STAT_OWNER_TO_OBS_OWNER = {
    'NOAA/NOS': 'NOAA',
    'USGS': 'USGS',
    'NCEM': 'CONTRAIL',
}

DEFAULT_ELEV_STAT_OWNERS = tuple(ELEV_STAT_OWNER_TO_OBS_OWNER.keys())


def read_elev_stat_stations(
        elev_stat_path: str | Path,
        owners: list[str] | tuple[str, ...] | None = None) -> list[dict[str, Any]]:
    """Read station entries from an ADCIRC elev_stat.151 file.

    Each line has the form::

        lon lat ! station_id ! owner ! name

    Parameters
    ----------
    elev_stat_path : path-like
        Path to elev_stat.151.
    owners : sequence of str, optional
        If given, only stations whose owner label is in this list are returned.
        Owner labels are matched exactly (e.g. ``NOAA/NOS``, not ``NOAA``).

    Returns
    -------
    list of dict
        Keys: ``lon``, ``lat``, ``station_id``, ``owner``, ``name``, ``line_no``.
    """
    elev_stat_path = Path(elev_stat_path)
    owner_filter = set(owners) if owners is not None else None
    stations: list[dict[str, Any]] = []

    with open(elev_stat_path, encoding='utf-8') as f:
        header = f.readline().strip()
        if not header:
            raise ValueError(f'{elev_stat_path}: empty file')
        try:
            n_expected = int(header)
        except ValueError as exc:
            raise ValueError(
                f'{elev_stat_path}: first line must be station count, got {header!r}'
            ) from exc

        for line_no, line in enumerate(f, start=2):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            station = _parse_elev_stat_line(line, line_no)
            if owner_filter is not None and station['owner'] not in owner_filter:
                continue
            stations.append(station)

    if len(stations) != n_expected and owner_filter is None:
        print(
            f'Warning: {elev_stat_path}: header count is {n_expected} but '
            f'{len(stations)} station lines were parsed',
            file=sys.stderr,
        )
    return stations


def _parse_elev_stat_line(line: str, line_no: int) -> dict[str, Any]:
    parts = [part.strip() for part in line.split('!')]
    if len(parts) < 4:
        raise ValueError(
            f'Line {line_no}: expected "lon lat ! id ! owner ! name", got: {line!r}'
        )
    coord_tokens = parts[0].split()
    if len(coord_tokens) < 2:
        raise ValueError(
            f'Line {line_no}: could not parse longitude and latitude from {parts[0]!r}'
        )
    return {
        'lon': float(coord_tokens[0]),
        'lat': float(coord_tokens[1]),
        'station_id': parts[1],
        'owner': parts[2],
        'name': parts[3],
        'line_no': line_no,
    }


def load_station_ids_file(path: str | Path) -> list[str]:
    """Load station ids from a text file (one id per line; ``#`` comments allowed)."""
    ids: list[str] = []
    with open(path, encoding='utf-8') as f:
        for line_no, line in enumerate(f, start=1):
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            for token in line.split():
                ids.append(token)
    if not ids:
        raise ValueError(f'{path}: no station ids found')
    return ids


def select_elev_stat_stations_by_ids(
        stations: list[dict[str, Any]],
        station_ids: list[str] | tuple[str, ...]) -> tuple[list[dict[str, Any]], list[str]]:
    """Select elev_stat rows by fort.61 station id (case-insensitive).

    Parameters
    ----------
    stations : list of dict
        Stations from :func:`read_elev_stat_stations`.
    station_ids : sequence of str
        Station ids to keep, in the order requested. If the elev_stat file
        contains multiple rows with the same id, all matching rows are included
        for each requested id.

    Returns
    -------
    selected, missing
        ``selected`` is the filtered station list; ``missing`` lists requested
        ids that were not found among ``stations``.
    """
    by_id: dict[str, list[dict[str, Any]]] = {}
    for station in stations:
        key = station['station_id'].strip().upper()
        by_id.setdefault(key, []).append(station)

    selected: list[dict[str, Any]] = []
    missing: list[str] = []
    for station_id in station_ids:
        key = station_id.strip().upper()
        matches = by_id.get(key)
        if not matches:
            missing.append(station_id)
            continue
        selected.extend(matches)
    return selected, missing


def map_elev_stat_owner_to_obs_owner(owner: str) -> str:
    """Map elev_stat.151 owner label to :func:`plot_hydrograph_at_station` owner."""
    try:
        return ELEV_STAT_OWNER_TO_OBS_OWNER[owner]
    except KeyError as exc:
        valid = ', '.join(sorted(ELEV_STAT_OWNER_TO_OBS_OWNER))
        raise ValueError(
            f'Unsupported elev_stat owner {owner!r}. '
            f'Supported mappings: {valid}'
        ) from exc


def sanitize_filename_component(text: str) -> str:
    """Make a string safe for use in output file names."""
    text = text.strip().replace(' ', '_')
    return re.sub(r'[^A-Za-z0-9._-]+', '_', text)


def format_output_path(
        pattern: str,
        index: int,
        owner: str,
        station_id: str,
        name: str,
        output_dir: Path) -> Path:
    safe_owner = sanitize_filename_component(owner)
    safe_id = sanitize_filename_component(station_id)
    safe_name = sanitize_filename_component(name)
    filename = pattern.format(
        index=index,
        owner=safe_owner,
        station_id=safe_id,
        id=safe_id,
        name=safe_name,
    )
    return output_dir / filename


def expand_f61or63_patterns(patterns: list[str]) -> list[str]:
    """Expand glob patterns to a flat list of NetCDF paths."""
    files: list[str] = []
    for pattern in patterns:
        files.extend(glob(pattern))
    if not files:
        raise ValueError('No files matched --f61or63files patterns')
    return files


def flatten_f61or63_file_list(f61or63files) -> list[str]:
    """Flatten nested file lists produced when ``f61or63concat`` is used."""
    flat: list[str] = []
    for item in f61or63files:
        if isinstance(item, list):
            flat.extend(flatten_f61or63_file_list(item))
        else:
            flat.append(item)
    return flat


def _read_netcdf_time_series(path: str | Path):
    """Return pandas DatetimeIndex for the time dimension of a fort.61/63 NetCDF file."""
    import numpy as np
    import pandas as pd
    import xarray as xr

    with xr.open_dataset(path) as ds:
        if 'time' not in ds:
            raise ValueError(f'{path}: no time variable')
        tvar = ds['time']
        values = tvar.values
        if np.issubdtype(np.asarray(values).dtype, np.number):
            units = str(tvar.attrs.get('units', ''))
            if 'since' in units:
                origin = units.split('since', 1)[1].strip()
                return pd.to_datetime(values, unit='s', origin=origin)
            base_date = tvar.attrs.get('base_date')
            if base_date is not None:
                origin = pd.Timestamp(str(base_date))
                return origin + pd.to_timedelta(values, unit='s')
            raise ValueError(
                f'{path}: numeric time values without decodable units/base_date'
            )
        return pd.to_datetime(values.astype('datetime64[ms]'))


def infer_date_range_from_f61or63_files(
        f61or63files: list[str]) -> tuple[datetime, datetime]:
    """Infer inclusive plot date range from fort.61.nc / fort.63.nc time axes.

    Uses the earliest calendar date among all first timestamps and the latest
    calendar date among all last timestamps across the given files.
    """
    paths = flatten_f61or63_file_list(f61or63files)
    if not paths:
        raise ValueError('No fort.61/63 files provided for date inference')

    t_min = None
    t_max = None
    for path in paths:
        times = _read_netcdf_time_series(path)
        if len(times) == 0:
            raise ValueError(f'{path}: time dimension is empty')
        file_min = times.min().to_pydatetime()
        file_max = times.max().to_pydatetime()
        t_min = file_min if t_min is None else min(t_min, file_min)
        t_max = file_max if t_max is None else max(t_max, file_max)

    date_start = datetime(t_min.year, t_min.month, t_min.day)
    date_end = datetime(t_max.year, t_max.month, t_max.day)
    return date_start, date_end


def resolve_plot_date_range(
        date_start_str: str | None,
        date_end_str: str | None,
        f61or63files: list[str]) -> tuple[datetime, datetime]:
    """Resolve plot bounds from CLI strings and/or fort.61/63 file metadata."""
    inferred_start, inferred_end = None, None
    if date_start_str is None or date_end_str is None:
        inferred_start, inferred_end = infer_date_range_from_f61or63_files(
            f61or63files
        )

    if date_start_str is None:
        if inferred_start is None:
            raise ValueError(
                'date-start is required when it cannot be inferred from '
                '--f61or63files'
            )
        date_start = inferred_start
    else:
        date_start = datetime.strptime(date_start_str, '%Y-%m-%d')

    if date_end_str is None:
        if inferred_end is None:
            raise ValueError(
                'date-end is required when it cannot be inferred from '
                '--f61or63files'
            )
        date_end = inferred_end
    else:
        date_end = datetime.strptime(date_end_str, '%Y-%m-%d')

    if date_start > date_end:
        raise ValueError(
            f'date-start {date_start.date()} is after date-end {date_end.date()}'
        )
    return date_start, date_end


def _prepare_f61or63_lists(
        args, expanded_files: list[str] | None = None
) -> tuple[list, list | None, list | None]:
    """Expand glob patterns and parse start times (mirrors plot_hydrograph_at_station.main)."""
    f61or63files: list = (
        list(expanded_files)
        if expanded_files is not None
        else expand_f61or63_patterns(args.f61or63files)
    )

    if args.f61or63concat:
        f61or63files = [f61or63files]

    if args.f61or63starts:
        f61or63starts = [
            datetime.strptime(f61or63start, '%Y-%m-%d') if f61or63start else None
            for f61or63start in args.f61or63starts
        ]
    else:
        f61or63starts = None

    f63files_fallback = None
    if getattr(args, 'f63files_fallback', None):
        f63files_fallback = []
        for f63file in args.f63files_fallback:
            f63files_fallback.extend(glob(f63file))
        if args.f61or63concat:
            f63files_fallback = [f63files_fallback]

    return f61or63files, f61or63starts, f63files_fallback


def _contrail_options(args) -> dict[str, str] | None:
    username = getattr(args, 'username', None) or os.environ.get(
        getattr(args, 'contrail_username_env', 'CONTRAIL_USERNAME')
    )
    password = getattr(args, 'password', None) or os.environ.get(
        getattr(args, 'contrail_password_env', 'CONTRAIL_PASSWORD')
    )
    if not username or not password:
        return None
    return {
        'username': username,
        'password': password,
        'sensor_type': getattr(args, 'sensor_type', 'water_elevation'),
    }


def plot_f61_hydrographs_from_elev_stat(
        elev_stat_path: str | Path,
        output_dir: str | Path,
        date_start: datetime,
        date_end: datetime,
        f61or63files,
        *,
        owners: list[str] | tuple[str, ...] | None = None,
        station_datum: str | None = 'NAVD',
        f61or63starts=None,
        f61or63labels=None,
        f61or63colors=None,
        f63files_fallback=None,
        plot_movingaverage: bool = False,
        adjust_datum_by_mean_error_period_days: float = 0,
        cache_dir: str | Path | None = None,
        contrail_options: dict[str, str] | None = None,
        contrail_station_id_type: str | None = 'f61',
        filename_pattern: str = '{index:04d}_{owner}_{station_id}_{name}.png',
        station_ids: list[str] | tuple[str, ...] | None = None,
        skip_on_error: bool = False,
        skip_existing: bool = False,
        max_stations: int | None = None,
        plot_in_foot: bool = False,
        figsize: tuple[float, float] = (10.0, 6.0)) -> tuple[list[Path], int]:
    """Plot hydrographs for selected elev_stat.151 stations.

    Parameters
    ----------
    station_ids : sequence of str, optional
        If given, only stations whose elev_stat id matches one of these values
        (case-insensitive) are plotted, in the order given.
    skip_existing : bool, optional
        If True, skip stations whose output figure file already exists.
    figsize : tuple of float, optional
        Figure size ``(width, height)`` in inches (default: ``(10.0, 6.0)``).

    Returns
    -------
    written, skipped_existing
        ``written`` lists paths of figures saved this run; ``skipped_existing``
        is how many stations were skipped because the output file existed.
    """
    import matplotlib.pyplot as plt

    from vewutils.plot.plot_hydrograph_at_station import plot_hydrograph_at_station

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stations = read_elev_stat_stations(elev_stat_path, owners=owners)
    if station_ids:
        stations, missing_ids = select_elev_stat_stations_by_ids(
            stations, station_ids
        )
        if missing_ids:
            print(
                'Warning: station id(s) not found in elev_stat (after owner '
                f'filter): {", ".join(missing_ids)}',
                file=sys.stderr,
            )
        if not stations:
            raise ValueError(
                'No stations to plot: none of the requested --station-ids '
                'match the elev_stat file'
            )
    if max_stations is not None:
        stations = stations[:max_stations]

    written: list[Path] = []
    skipped_existing = 0
    n_contrail = sum(
        1 for st in stations
        if map_elev_stat_owner_to_obs_owner(st['owner']) == 'CONTRAIL'
    )
    if n_contrail > 0 and not contrail_options:
        raise ValueError(
            'CONTRAIL credentials are required for NCEM stations. '
            'Set CONTRAIL_USERNAME and CONTRAIL_PASSWORD or pass --username and --password.'
        )

    for index, station in enumerate(stations, start=1):
        owner_label = station['owner']
        obs_owner = map_elev_stat_owner_to_obs_owner(owner_label)
        station_id = station['station_id']
        output_path = format_output_path(
            filename_pattern,
            index,
            owner_label,
            station_id,
            station['name'],
            output_dir,
        )
        if skip_existing and output_path.is_file():
            skipped_existing += 1
            print(
                f'Skipping {index}/{len(stations)}: {owner_label} '
                f'{station_id} ({station["name"]}) -> {output_path.name} exists'
            )
            continue

        print(
            f'Processing {index}/{len(stations)}: {owner_label} '
            f'{station_id} ({station["name"]}) -> {output_path.name}'
        )

        options = contrail_options if obs_owner == 'CONTRAIL' else None
        station_id_type = (
            contrail_station_id_type if obs_owner == 'CONTRAIL' else None
        )

        try:
            fig, ax = plt.subplots(figsize=figsize)
            plot_hydrograph_at_station(
                fig,
                ax,
                obs_owner,
                station_id,
                station['lon'],
                station['lat'],
                station_datum,
                date_start,
                date_end,
                f61or63files,
                f61or63starts,
                f61or63labels,
                f61or63colors,
                f63files_fallback=f63files_fallback,
                plot_movingaverage=plot_movingaverage,
                plot_in_foot=plot_in_foot,
                options=options,
                cache_dir=cache_dir,
                station_id_type=station_id_type,
                adjust_datum_by_mean_error_period_days=(
                    adjust_datum_by_mean_error_period_days
                ),
            )
            fig.savefig(output_path)
            plt.close(fig)
            written.append(output_path)
        except Exception as exc:
            plt.close('all')
            if skip_on_error:
                print(
                    f'Warning: skipped station {station_id} ({owner_label}): {exc}',
                    file=sys.stderr,
                )
                continue
            raise

    return written, skipped_existing


def get_parser():
    import argparse

    parser = argparse.ArgumentParser(
        add_help=False,
        description=(
            'Plot observed and modeled hydrographs for NOAA/NOS, USGS, and NCEM '
            'stations listed in an elev_stat.151 file. NCEM uses CONTRAIL '
            'observations; NOAA/NOS maps to NOAA.'
        ),
    )
    parser.add_argument(
        '--elev-stat',
        required=True,
        help='Path to elev_stat.151 (first line is station count)',
    )
    parser.add_argument(
        '--owners',
        nargs='+',
        default=list(DEFAULT_ELEV_STAT_OWNERS),
        help=(
            'elev_stat owner labels to include (default: %(default)s). '
            'NCEM is fetched via CONTRAIL.'
        ),
    )
    parser.add_argument(
        '--station-ids',
        nargs='+',
        default=None,
        metavar='ID',
        help=(
            'Plot only these elev_stat / fort.61 station ids (case-insensitive). '
            'Order is preserved. If omitted, all stations matching --owners '
            'are plotted.'
        ),
    )
    parser.add_argument(
        '--station-ids-file',
        default=None,
        help=(
            'File with station ids to plot (one per line; # comments allowed). '
            'Combined with --station-ids when both are given.'
        ),
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Directory for per-station PNG files',
    )
    parser.add_argument(
        '--filename-pattern',
        default='{index:04d}_{owner}_{station_id}_{name}.png',
        help=(
            'Output filename pattern per station. Fields: index, owner, '
            'station_id (alias id), name'
        ),
    )
    parser.add_argument(
        '--skip-on-error',
        action='store_true',
        help='Continue with remaining stations if one plot fails',
    )
    parser.add_argument(
        '--skip-existing',
        action='store_true',
        help=(
            'Skip stations whose output figure already exists in --output-dir '
            '(useful for resuming interrupted runs)'
        ),
    )
    parser.add_argument(
        '--max-stations',
        type=int,
        default=None,
        help='Plot at most this many stations (for testing)',
    )
    parser.add_argument(
        '--station-datum',
        type=str,
        default='NAVD',
        help='Datum for observations (e.g. NAVD, MSL)',
    )
    parser.add_argument(
        '--date-start',
        type=str,
        default=None,
        help=(
            'Plot/observation start date (YYYY-MM-DD). Default: first calendar '
            'date in --f61or63files time axis'
        ),
    )
    parser.add_argument(
        '--date-end',
        type=str,
        default=None,
        help=(
            'Plot/observation end date (YYYY-MM-DD). Default: last calendar '
            'date in --f61or63files time axis'
        ),
    )
    parser.add_argument(
        '--f61or63files',
        type=str,
        nargs='+',
        required=True,
        help='fort.61.nc or fort.63.nc file(s); glob patterns allowed',
    )
    parser.add_argument(
        '--f61or63starts',
        type=str,
        nargs='+',
        default=None,
        help='Optional start date (YYYY-MM-DD) per file',
    )
    parser.add_argument(
        '--f61or63labels',
        type=str,
        nargs='+',
        default=None,
        help='Legend labels per model file',
    )
    parser.add_argument(
        '--f61or63colors',
        type=str,
        nargs='+',
        default=None,
        help='Line colors per model file',
    )
    parser.add_argument(
        '--f61or63concat',
        action='store_true',
        help='Concatenate all matched f61/f63 files into one series',
    )
    parser.add_argument(
        '--f63files-fallback',
        type=str,
        nargs='+',
        default=None,
        help='fort.63.nc fallback when a station is missing from fort.61.nc',
    )
    parser.add_argument(
        '--plot-movingaverage',
        action='store_true',
        help='Plot moving average overlays',
    )
    parser.add_argument(
        '--adjust-datum-by-mean-error',
        type=float,
        default=0,
        metavar='period_in_days',
        help='Adjust model datum using mean obs-model error over first N days',
    )
    parser.add_argument(
        '--cache-dir',
        help='Cache directory for observation downloads',
    )
    parser.add_argument(
        '--plot-in-foot',
        action='store_true',
        help='Plot water levels in feet',
    )
    parser.add_argument(
        '--fig-width',
        type=float,
        default=12.0,
        help='Figure width in inches (default: 10.0)',
    )
    parser.add_argument(
        '--fig-height',
        type=float,
        default=5.0,
        help='Figure height in inches (default: 6.0)',
    )

    contrail_group = parser.add_argument_group('CONTRAIL / NCEM options')
    contrail_group.add_argument(
        '--username',
        help='CONTRAIL username (default: CONTRAIL_USERNAME env var)',
    )
    contrail_group.add_argument(
        '--password',
        help='CONTRAIL password (default: CONTRAIL_PASSWORD env var)',
    )
    contrail_group.add_argument(
        '--contrail-username-env',
        default='CONTRAIL_USERNAME',
        help='Environment variable for CONTRAIL username (default: CONTRAIL_USERNAME)',
    )
    contrail_group.add_argument(
        '--contrail-password-env',
        default='CONTRAIL_PASSWORD',
        help='Environment variable for CONTRAIL password (default: CONTRAIL_PASSWORD)',
    )
    contrail_group.add_argument(
        '--sensor-type',
        choices=['water_elevation', 'stream_elevation', 'stage'],
        default='water_elevation',
        help='CONTRAIL sensor type (default: water_elevation)',
    )
    contrail_group.add_argument(
        '--station-id-type',
        choices=['auto', 'contrail', 'f61'],
        default='f61',
        help=(
            'How to interpret elev_stat station ids for CONTRAIL/NCEM '
            '(default: f61, i.e. fort.61 station code)'
        ),
    )
    return parser


def main(args=None):
    import matplotlib.pyplot as plt

    if args is None:
        args = get_parser().parse_args()

    station_ids: list[str] = []
    if args.station_ids:
        station_ids.extend(args.station_ids)
    if args.station_ids_file:
        station_ids.extend(load_station_ids_file(args.station_ids_file))
    station_ids = station_ids or None

    expanded_files = expand_f61or63_patterns(args.f61or63files)
    date_start, date_end = resolve_plot_date_range(
        args.date_start, args.date_end, expanded_files
    )
    if args.date_start is None or args.date_end is None:
        print(
            f'Inferred date range from model file(s): '
            f'{date_start.date()} to {date_end.date()}',
            file=sys.stderr,
        )

    f61or63files, f61or63starts, f63files_fallback = _prepare_f61or63_lists(
        args, expanded_files=expanded_files
    )

    contrail_options = _contrail_options(args)

    written, skipped_existing = plot_f61_hydrographs_from_elev_stat(
        args.elev_stat,
        args.output_dir,
        date_start,
        date_end,
        f61or63files,
        owners=args.owners,
        station_datum=args.station_datum,
        f61or63starts=f61or63starts,
        f61or63labels=args.f61or63labels,
        f61or63colors=args.f61or63colors,
        f63files_fallback=f63files_fallback,
        plot_movingaverage=args.plot_movingaverage,
        adjust_datum_by_mean_error_period_days=args.adjust_datum_by_mean_error,
        cache_dir=args.cache_dir,
        contrail_options=contrail_options,
        contrail_station_id_type=args.station_id_type,
        filename_pattern=args.filename_pattern,
        station_ids=station_ids,
        skip_on_error=args.skip_on_error,
        skip_existing=args.skip_existing,
        max_stations=args.max_stations,
        plot_in_foot=args.plot_in_foot,
        figsize=(args.fig_width, args.fig_height),
    )
    plt.close('all')
    print(f'Wrote {len(written)} figure(s) to {args.output_dir}')
    if skipped_existing:
        print(f'Skipped {skipped_existing} station(s) with existing figure(s)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
