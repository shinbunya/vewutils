"""Generate HTML reports with hydrograph maps and contour figure galleries."""

from __future__ import annotations

import argparse
import html
import os
import re
import sys
try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

from vewutils.plot.plot_f61_hydrographs import (
    DEFAULT_ELEV_STAT_OWNERS,
    DEFAULT_HYDROGRAPH_MAP_HTML,
    DEFAULT_MAP_THUMB_SCALE,
    generate_hydrograph_map_html,
    load_station_ids_file,
    plot_f61_hydrographs_from_elev_stat,
    resolve_plot_date_range,
)
from vewutils.plot.plot_max_ele_2d import plot_max_ele_2d

DEFAULT_LOOKBACK_DAYS = 5
DEFAULT_MODE = 'standard'
DEFAULT_REPORT_HTML = 'index.html'
DEFAULT_THUMB_WIDTH = 300
DEFAULT_SAMPLE_CONFIG_NAME = 'sample_report_config.toml'
ARCHIVE_PATH_RE = re.compile(r'/archive/(\d{8})/hour_(\d{2})/')

SAMPLE_REPORT_CONFIG_TOML = '''\
# Sample vewutils generate-report configuration.
#
# Fill in the paths marked "REQUIRED" for your simulation, then run:
#   vewutils plot generate-report --config this_file.toml

[report]
# REQUIRED. Directory containing the simulation's fort.61.nc / maxele NetCDF files.
# For forecast mode this must contain an "/archive/YYYYMMDD/hour_HH/" segment.
data_dir = "/path/to/simulation/archive/20260101/hour_00"

# Where the report (HTML, hydrograph figures, contour figures) is written.
# Defaults to "<data_dir>/report" if omitted.
# output_dir = "/path/to/output/report"

title = "ADCIRC Simulation Report"

# "standard" uses hydrographs.f61or63files (or fort.61.nc in data_dir) as-is.
# "forecast" parses data_dir as an archive path and concatenates lookback_days
# of prior analysis fort.61.nc files with the forecast fort.61.nc.
mode = "standard"

# Only used when mode = "forecast".
lookback_days = 5

# Hover-thumbnail scale for markers on the hydrograph station map.
map_thumb_scale = 1.0


[hydrographs]
# REQUIRED. Path to the ADCIRC elev_stat.151 station list.
elev_stat = "/path/to/elev_stat.151"

# Restrict which elev_stat owner labels are plotted (default: NOAA/NOS, USGS, NCEM).
# owners = ["NOAA/NOS", "USGS"]

# Optional explicit list/file of station ids to plot instead of all stations.
# station_ids = ["8651370", "8652587"]
# station_ids_file = "/path/to/station_ids.txt"

# Vertical datum used for plotted/observed water levels (default: "NAVD").
station_datum = "NAVD"

# fort.61/63 file(s) to read station output from. Ignored in forecast mode
# (which discovers files automatically). Defaults to "fort.61.nc" in data_dir
# if omitted.
# f61or63files = ["fort.61.nc"]

plot_movingaverage = false

# Directory used to cache downloaded observation data between runs.
# cache_dir = "/path/to/cache"

# Output filename pattern for per-station hydrograph PNGs.
filename_pattern = "{index:04d}_{owner}_{station_id}_{name}.png"

# CONTRAIL/NCEM observation credentials are read from --username/--password,
# hydrographs.username/password below, or the CONTRAIL_USERNAME/CONTRAIL_PASSWORD
# environment variables (or the env var names given by contrail_username_env /
# contrail_password_env). Leave unset if you have no NCEM/CONTRAIL stations.
# username = "..."
# password = "..."
# sensor_type = "water_elevation"


[contours]
# Defaults applied to every [[fields]] entry below unless overridden per-field
# or per-extent.
cmap = "viridis"
levels = 20
dpi = 300
figsizex = 12.0
figsizey = 10.0

# Thumbnail width (px) for the contour gallery tiles.
thumb_width = 300

draw_shorelines = false
drawmesh = false

# Optional storm track overlay (KMZ best-track file) drawn on every contour plot.
# track_file = "/path/to/track.kmz"
# track_color = "red"
# track_linewidth = 2.0
# track_markersize = 5.0
# track_annotate_datetime = false
# track_annotate_category = false


# Reusable extent (map zoom / colorbar range) presets, referenced by id from
# [[fields]] via extent_presets = ["region_a"] or [[fields.extents]] preset = "region_a".
[[extent_presets]]
id = "region_a"
label = "Region A"
xmin = -80.5
xmax = -79.5
ymin = 32.0
ymax = 33.0
vmin = 0.0
vmax = 3.0
cbar_label = "Water Surface Elevation (m, NAVD88)"
cbar_increment = 0.5


# Each [[fields]] entry is one NetCDF variable to render as contour figures.
# It must reference at least one extent, either via extent_presets and/or
# inline [[fields.extents]] tables.
[[fields]]
id = "maxele"
label = "Maximum Water Elevation"
file = "maxele.63.nc"          # relative to report.data_dir, or an absolute path
variable = "zeta_max"
extent_presets = ["region_a"]

[[fields]]
id = "maxwvel"
label = "Maximum Wind Speed"
file = "maxwvel.63.nc"
variable = "wind_max"

  [[fields.extents]]
  id = "full_domain"
  label = "Full Domain"
  vmin = 0.0
  vmax = 40.0
  cbar_label = "Wind Speed (m/s)"
'''

# Keys copied from an extent preset or inline extent table into resolved extents.
EXTENT_OPTION_KEYS = (
    'label', 'title', 'xmin', 'xmax', 'ymin', 'ymax',
    'vmin', 'vmax', 'cbar_label', 'cbar_increment', 'cbar_ticks_increment',
)


@dataclass(frozen=True)
class ArchiveContext:
    data_dir: Path
    archive_root: Path
    target_date: datetime
    target_hour: int


@dataclass(frozen=True)
class ContourRecord:
    field_id: str
    field_label: str
    extent_id: str
    extent_label: str
    image_path: Path
    thumb_path: Path
    title: str


def _parse_extent_presets(raw: dict[str, Any], path: Path) -> dict[str, dict[str, Any]]:
    """Parse [[extent_presets]] into a map of preset id -> options."""
    presets_raw = raw.get('extent_presets')
    if presets_raw is None:
        return {}
    if not isinstance(presets_raw, list):
        raise ValueError(f'{path}: [[extent_presets]] must be an array of tables')

    presets: dict[str, dict[str, Any]] = {}
    for index, preset in enumerate(presets_raw):
        if not isinstance(preset, dict):
            raise ValueError(f'{path}: extent_presets[{index}] must be a table')
        preset_id = preset.get('id')
        if not preset_id:
            raise ValueError(f'{path}: extent_presets[{index}].id is required')
        if preset_id in presets:
            raise ValueError(
                f'{path}: duplicate extent_presets id {preset_id!r}'
            )
        presets[preset_id] = dict(preset)
    return presets


def _merge_extent_options(
        base: dict[str, Any],
        overrides: dict[str, Any]) -> dict[str, Any]:
    """Merge preset and override tables; overrides win."""
    merged = dict(base)
    for key in EXTENT_OPTION_KEYS:
        if key in overrides and overrides[key] is not None:
            merged[key] = overrides[key]
    return merged


def _resolve_field_extents(
        field: dict[str, Any],
        presets: dict[str, dict[str, Any]],
        path: Path,
        field_index: int) -> list[dict[str, Any]]:
    """Resolve inline extents and/or extent_presets references for one field."""
    resolved: list[dict[str, Any]] = []

    preset_refs = field.get('extent_presets')
    if preset_refs is not None:
        if not isinstance(preset_refs, list) or not preset_refs:
            raise ValueError(
                f'{path}: fields[{field_index}].extent_presets must be a '
                'non-empty array of preset ids'
            )
        for ref_index, preset_id in enumerate(preset_refs):
            if not isinstance(preset_id, str):
                raise ValueError(
                    f'{path}: fields[{field_index}].extent_presets[{ref_index}] '
                    'must be a preset id string'
                )
            try:
                preset = presets[preset_id]
            except KeyError as exc:
                valid = ', '.join(sorted(presets))
                raise ValueError(
                    f'{path}: fields[{field_index}] references unknown '
                    f'extent preset {preset_id!r}. '
                    f'Available presets: {valid or "(none)"}'
                ) from exc
            resolved.append(_merge_extent_options(preset, {'id': preset_id}))

    inline_extents = field.get('extents')
    if inline_extents is not None:
        if not isinstance(inline_extents, list):
            raise ValueError(f'{path}: fields[{field_index}].extents must be an array')

        for extent_index, extent in enumerate(inline_extents):
            if not isinstance(extent, dict):
                raise ValueError(
                    f'{path}: fields[{field_index}].extents[{extent_index}] '
                    'must be a table'
                )
            preset_ref = extent.get('preset')
            if preset_ref:
                try:
                    preset = presets[preset_ref]
                except KeyError as exc:
                    valid = ', '.join(sorted(presets))
                    raise ValueError(
                        f'{path}: fields[{field_index}].extents[{extent_index}] '
                        f'references unknown preset {preset_ref!r}. '
                        f'Available presets: {valid or "(none)"}'
                    ) from exc
                extent_id = extent.get('id', preset_ref)
                merged = _merge_extent_options(preset, extent)
                merged['id'] = extent_id
            else:
                extent_id = extent.get('id')
                if not extent_id:
                    raise ValueError(
                        f'{path}: fields[{field_index}].extents[{extent_index}] '
                        'requires id (or preset referencing a preset with id)'
                    )
                merged = dict(extent)
            resolved.append(merged)

    if not resolved:
        raise ValueError(
            f'{path}: fields[{field_index}] must specify extent_presets and/or '
            '[[fields.extents]] entries'
        )
    return resolved


def write_sample_config(path: str | Path, *, force: bool = False) -> Path:
    """Write a sample report TOML config to path; refuse to overwrite unless force."""
    path = Path(path)
    if path.exists() and not force:
        raise FileExistsError(f'{path} already exists; use --force to overwrite')
    path.write_text(SAMPLE_REPORT_CONFIG_TOML, encoding='utf-8')
    return path


def load_report_config(path: str | Path) -> dict[str, Any]:
    """Load and validate a report TOML configuration file."""
    path = Path(path)
    with open(path, 'rb') as f:
        raw = tomllib.load(f)

    report = raw.get('report')
    if not isinstance(report, dict):
        raise ValueError(f'{path}: [report] section is required')

    data_dir = report.get('data_dir')
    if not data_dir:
        raise ValueError(f'{path}: report.data_dir is required')
    data_dir = Path(data_dir).resolve()
    if not data_dir.is_dir():
        raise ValueError(f'{path}: report.data_dir is not a directory: {data_dir}')

    output_dir = report.get('output_dir')
    if output_dir:
        output_dir = Path(output_dir).resolve()
    else:
        output_dir = data_dir / 'report'

    hydrographs = raw.get('hydrographs')
    if hydrographs is not None and not isinstance(hydrographs, dict):
        raise ValueError(f'{path}: [hydrographs] must be a table')
    hydrographs = dict(hydrographs or {})
    elev_stat = hydrographs.get('elev_stat')
    if not elev_stat:
        raise ValueError(f'{path}: hydrographs.elev_stat is required')
    elev_stat_path = Path(elev_stat)
    if not elev_stat_path.is_file():
        raise ValueError(f'{path}: hydrographs.elev_stat not found: {elev_stat_path}')

    fields = raw.get('fields')
    if not fields:
        raise ValueError(f'{path}: at least one [[fields]] entry is required')
    if not isinstance(fields, list):
        raise ValueError(f'{path}: [[fields]] must be an array of tables')

    extent_presets = _parse_extent_presets(raw, path)

    validated_fields: list[dict[str, Any]] = []
    for index, field in enumerate(fields):
        if not isinstance(field, dict):
            raise ValueError(f'{path}: fields[{index}] must be a table')
        field_id = field.get('id')
        if not field_id:
            raise ValueError(f'{path}: fields[{index}].id is required')
        netcdf_file = field.get('file')
        if not netcdf_file:
            raise ValueError(f'{path}: fields[{index}].file is required')
        variable = field.get('variable')
        if not variable:
            raise ValueError(f'{path}: fields[{index}].variable is required')

        validated_extents = _resolve_field_extents(
            field, extent_presets, path, index
        )

        validated_fields.append({
            **field,
            'extents': validated_extents,
        })

    contours = raw.get('contours')
    if contours is not None and not isinstance(contours, dict):
        raise ValueError(f'{path}: [contours] must be a table')

    return {
        'report': {
            **report,
            'data_dir': data_dir,
            'output_dir': output_dir,
        },
        'hydrographs': hydrographs,
        'contours': dict(contours or {}),
        'extent_presets': extent_presets,
        'fields': validated_fields,
    }


def parse_archive_context(data_dir: str | Path) -> ArchiveContext:
    """Parse archive date/hour and archive root from a simulation data directory."""
    data_dir = Path(data_dir).resolve()
    match = ARCHIVE_PATH_RE.search(str(data_dir))
    if not match:
        raise ValueError(
            f'Could not parse archive date/hour from data_dir: {data_dir}. '
            'Expected path containing /archive/YYYYMMDD/hour_HH/'
        )

    archive_root = data_dir
    while archive_root.name != 'archive' and archive_root.parent != archive_root:
        archive_root = archive_root.parent
    if archive_root.name != 'archive':
        raise ValueError(f'Could not locate archive root for data_dir: {data_dir}')

    target_date = datetime.strptime(match.group(1), '%Y%m%d')
    target_hour = int(match.group(2))
    return ArchiveContext(
        data_dir=data_dir,
        archive_root=archive_root,
        target_date=target_date,
        target_hour=target_hour,
    )


def discover_forecast_f61_files(
        ctx: ArchiveContext,
        lookback_days: int = DEFAULT_LOOKBACK_DAYS) -> list[str]:
    """Collect analysis fort.61.nc files for lookback window plus forecast file."""
    if lookback_days < 0:
        raise ValueError(f'lookback_days must be non-negative, got {lookback_days}')

    start_date = ctx.target_date - timedelta(days=lookback_days)
    end_date = ctx.target_date
    target_ymd = ctx.target_date.strftime('%Y%m%d')
    files: list[str] = []

    current = start_date
    while current <= end_date:
        ymd = current.strftime('%Y%m%d')
        day_dir = ctx.archive_root / ymd
        if not day_dir.is_dir():
            print(
                f'Warning: missing archive day directory: {day_dir}',
                file=sys.stderr,
            )
            current += timedelta(days=1)
            continue

        hour_dirs = sorted(day_dir.glob('hour_*'))
        for hour_dir in hour_dirs:
            hour = int(hour_dir.name.split('_', 1)[1])
            if ymd == target_ymd and hour > ctx.target_hour:
                continue
            f61_path = hour_dir / 'adcirc' / 'analysis' / 'fort.61.nc'
            if f61_path.is_file():
                files.append(str(f61_path))
            else:
                print(
                    f'Warning: missing analysis fort.61.nc: {f61_path}',
                    file=sys.stderr,
                )
        current += timedelta(days=1)

    forecast_f61 = ctx.data_dir / 'fort.61.nc'
    if forecast_f61.is_file():
        files.append(str(forecast_f61))
    else:
        raise ValueError(f'Forecast fort.61.nc not found: {forecast_f61}')

    return files


def _resolve_netcdf_path(data_dir: Path, file_spec: str) -> Path:
    path = Path(file_spec)
    if not path.is_absolute():
        path = data_dir / path
    return path.resolve()


def _contrail_options_from_config(
        hydrographs: dict[str, Any],
        args) -> dict[str, str] | None:
    username = getattr(args, 'username', None) or hydrographs.get('username')
    if not username:
        username = os.environ.get(
            hydrographs.get('contrail_username_env', 'CONTRAIL_USERNAME')
        )
    password = getattr(args, 'password', None) or hydrographs.get('password')
    if not password:
        password = os.environ.get(
            hydrographs.get('contrail_password_env', 'CONTRAIL_PASSWORD')
        )
    if not username or not password:
        return None
    return {
        'username': username,
        'password': password,
        'sensor_type': hydrographs.get('sensor_type', 'water_elevation'),
    }


def _resolve_f61or63files(
        config: dict[str, Any],
        *,
        mode: str,
        lookback_days: int) -> list[str]:
    data_dir = config['report']['data_dir']
    hydrographs = config['hydrographs']

    if mode == 'forecast':
        ctx = parse_archive_context(data_dir)
        return discover_forecast_f61_files(ctx, lookback_days=lookback_days)

    explicit = hydrographs.get('f61or63files')
    if explicit:
        files: list[str] = []
        for pattern in explicit:
            path = Path(pattern)
            if path.is_absolute():
                if path.is_file():
                    files.append(str(path))
                else:
                    raise ValueError(f'f61or63files entry not found: {path}')
            else:
                resolved = (data_dir / pattern).resolve()
                if resolved.is_file():
                    files.append(str(resolved))
                else:
                    raise ValueError(f'f61or63files entry not found: {resolved}')
        return files

    default_f61 = data_dir / 'fort.61.nc'
    if not default_f61.is_file():
        raise ValueError(
            f'No hydrographs.f61or63files configured and {default_f61} not found'
        )
    return [str(default_f61)]


def generate_hydrographs(
        config: dict[str, Any],
        output_dir: Path,
        *,
        mode: str,
        lookback_days: int,
        skip_existing: bool,
        skip_on_error: bool,
        contrail_options: dict[str, str] | None) -> tuple[list[Path], list[dict[str, Any]]]:
    """Plot station hydrographs and return written paths and station records."""
    hydrographs = config['hydrographs']
    hydrograph_dir = output_dir / 'hydrographs'
    hydrograph_dir.mkdir(parents=True, exist_ok=True)

    f61or63files = _resolve_f61or63files(
        config,
        mode=mode,
        lookback_days=lookback_days,
    )
    print(
        f'Using {len(f61or63files)} fort.61/63 file(s) for hydrographs '
        f'(mode={mode})',
        file=sys.stderr,
    )

    date_start_str = hydrographs.get('date_start')
    date_end_str = hydrographs.get('date_end')
    date_start, date_end = resolve_plot_date_range(
        date_start_str,
        date_end_str,
        f61or63files,
    )
    if date_start_str is None or date_end_str is None:
        print(
            f'Inferred hydrograph date range: '
            f'{date_start.date()} to {date_end.date()}',
            file=sys.stderr,
        )

    owners = hydrographs.get('owners', list(DEFAULT_ELEV_STAT_OWNERS))
    station_ids = hydrographs.get('station_ids')
    if station_ids is not None and not isinstance(station_ids, list):
        raise ValueError('hydrographs.station_ids must be a list of station ids')

    station_ids_file = hydrographs.get('station_ids_file')
    if station_ids_file:
        loaded_ids = load_station_ids_file(station_ids_file)
        station_ids = list(station_ids or []) + loaded_ids
    station_ids = station_ids or None

    f61or63concat = hydrographs.get('f61or63concat', mode == 'forecast')
    if f61or63concat:
        f61or63files = [f61or63files]

    nowcast_forecast_style = hydrographs.get(
        'nowcast_forecast_style',
        mode == 'forecast',
    )
    connect = hydrographs.get('connect', mode == 'forecast')

    written, _skipped, station_records = plot_f61_hydrographs_from_elev_stat(
        hydrographs['elev_stat'],
        hydrograph_dir,
        date_start,
        date_end,
        f61or63files,
        owners=owners,
        station_datum=hydrographs.get('station_datum', 'NAVD'),
        f61or63starts=hydrographs.get('f61or63starts'),
        f61or63labels=hydrographs.get('f61or63labels'),
        f61or63colors=hydrographs.get('f61or63colors'),
        f63files_fallback=hydrographs.get('f63files_fallback'),
        plot_movingaverage=hydrographs.get('plot_movingaverage', False),
        adjust_datum_by_mean_error_period_days=hydrographs.get(
            'adjust_datum_by_mean_error',
            0,
        ),
        cache_dir=hydrographs.get('cache_dir'),
        contrail_options=contrail_options,
        contrail_station_id_type=hydrographs.get('station_id_type', 'f61'),
        filename_pattern=hydrographs.get(
            'filename_pattern',
            '{index:04d}_{owner}_{station_id}_{name}.png',
        ),
        station_ids=station_ids,
        skip_on_error=skip_on_error,
        skip_existing=skip_existing,
        max_stations=hydrographs.get('max_stations'),
        plot_in_foot=hydrographs.get('plot_in_foot', False),
        connect=connect,
        nowcast_forecast_style=nowcast_forecast_style,
        figsize=(
            hydrographs.get('fig_width', 12.0),
            hydrographs.get('fig_height', 5.0),
        ),
    )
    return written, station_records


def _field_plot_kwargs(
        field: dict[str, Any],
        extent: dict[str, Any],
        contours_cfg: dict[str, Any],
        netcdf_path: Path) -> dict[str, Any]:
    title = extent.get('title')
    if not title:
        field_label = field.get('label', field['id'])
        extent_label = extent.get('label', extent['id'])
        title = f'{field_label} — {extent_label}'

    track_file = field.get('track_file', contours_cfg.get('track_file'))
    if track_file and not Path(track_file).is_absolute():
        track_file = str(Path(track_file).resolve())

    kwargs: dict[str, Any] = {
        'maxele_file': str(netcdf_path),
        'variable': field['variable'],
        'title': title,
        'cmap': field.get('cmap', contours_cfg.get('cmap', 'viridis')),
        'drawmesh': field.get('drawmesh', contours_cfg.get('drawmesh', False)),
        'levels': field.get('levels', contours_cfg.get('levels', 20)),
        'draw_shorelines': field.get(
            'draw_shorelines',
            contours_cfg.get('draw_shorelines', False),
        ),
        'track_file': track_file,
        'track_color': field.get(
            'track_color',
            contours_cfg.get('track_color', 'red'),
        ),
        'track_linewidth': field.get(
            'track_linewidth',
            contours_cfg.get('track_linewidth', 2.0),
        ),
        'track_markersize': field.get(
            'track_markersize',
            contours_cfg.get('track_markersize', 5.0),
        ),
        'track_annotate_datetime': field.get(
            'track_annotate_datetime',
            contours_cfg.get('track_annotate_datetime', False),
        ),
        'track_annotate_category': field.get(
            'track_annotate_category',
            contours_cfg.get('track_annotate_category', False),
        ),
        'track_annotate_category_inside_circle': field.get(
            'track_annotate_category_inside_circle',
            contours_cfg.get('track_annotate_category_inside_circle', False),
        ),
        'track_annotate_non_hurricane_inside_circle': field.get(
            'track_annotate_non_hurricane_inside_circle',
            contours_cfg.get('track_annotate_non_hurricane_inside_circle', False),
        ),
        'track_annotate_fontsize': field.get(
            'track_annotate_fontsize',
            contours_cfg.get('track_annotate_fontsize', 8.0),
        ),
    }

    for key in (
        'vmin', 'vmax', 'cbar_label', 'cbar_increment', 'cbar_ticks_increment',
        'xmin', 'xmax', 'ymin', 'ymax',
    ):
        value = extent.get(key, field.get(key, contours_cfg.get(key)))
        if value is not None:
            kwargs[key] = value

    return kwargs


def generate_contour_figures(
        config: dict[str, Any],
        output_dir: Path,
        *,
        skip_existing: bool) -> list[ContourRecord]:
    """Generate contour PNGs and thumbnails for all configured fields/extents."""
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt

    data_dir = config['report']['data_dir']
    contours_cfg = config['contours']
    contour_dir = output_dir / 'contours'
    contour_dir.mkdir(parents=True, exist_ok=True)

    figsizex = contours_cfg.get('figsizex', 12.0)
    figsizey = contours_cfg.get('figsizey', 10.0)
    dpi = contours_cfg.get('dpi', 300)
    thumb_width = contours_cfg.get('thumb_width', DEFAULT_THUMB_WIDTH)

    records: list[ContourRecord] = []

    for field in config['fields']:
        netcdf_path = _resolve_netcdf_path(data_dir, field['file'])
        if not netcdf_path.is_file():
            raise ValueError(f'Contour NetCDF file not found: {netcdf_path}')

        field_label = field.get('label', field['id'])
        for extent in field['extents']:
            extent_id = extent['id']
            extent_label = extent.get('label', extent_id)
            image_name = f"{field['id']}_{extent_id}.png"
            thumb_name = f"{field['id']}_{extent_id}_thumb.png"
            image_path = contour_dir / image_name
            thumb_path = contour_dir / thumb_name

            if skip_existing and image_path.is_file() and thumb_path.is_file():
                print(f'Skipping existing contour figure: {image_path.name}')
                records.append(ContourRecord(
                    field_id=field['id'],
                    field_label=field_label,
                    extent_id=extent_id,
                    extent_label=extent_label,
                    image_path=image_path,
                    thumb_path=thumb_path,
                    title=extent.get('title', f'{field_label} — {extent_label}'),
                ))
                continue

            plot_kwargs = _field_plot_kwargs(
                field,
                extent,
                contours_cfg,
                netcdf_path,
            )
            print(f'Creating contour figure: {image_path.name}')

            fig, ax = plt.subplots(figsize=(figsizex, figsizey))
            success = plot_max_ele_2d(fig, ax, **plot_kwargs)
            if not success:
                plt.close(fig)
                raise RuntimeError(
                    f'Failed to create contour plot for {field["id"]}/{extent_id}'
                )

            fig.savefig(image_path, dpi=dpi, bbox_inches='tight')
            plt.close(fig)

            img_array = mpimg.imread(image_path)
            height, width = img_array.shape[:2]
            step = max(1, width // thumb_width)
            thumb_array = img_array[::step, ::step]
            mpimg.imsave(thumb_path, thumb_array)

            records.append(ContourRecord(
                field_id=field['id'],
                field_label=field_label,
                extent_id=extent_id,
                extent_label=extent_label,
                image_path=image_path,
                thumb_path=thumb_path,
                title=plot_kwargs['title'],
            ))

    plt.close('all')
    return records


def assemble_report_html(
        config: dict[str, Any],
        station_records: list[dict[str, Any]],
        contour_records: list[ContourRecord],
        output_dir: Path,
        *,
        mode: str,
        lookback_days: int,
        map_thumb_scale: float = DEFAULT_MAP_THUMB_SCALE) -> Path:
    """Write hydrograph map HTML and the main report index page."""
    report_cfg = config['report']
    title = report_cfg.get('title', 'ADCIRC Simulation Report')
    hydrograph_dir = output_dir / 'hydrographs'
    map_path = hydrograph_dir / DEFAULT_HYDROGRAPH_MAP_HTML

    generate_hydrograph_map_html(
        map_path,
        station_records,
        title=f'{title} — Hydrograph Stations',
        thumb_scale=map_thumb_scale,
    )

    grouped: dict[str, list[ContourRecord]] = {}
    for record in contour_records:
        grouped.setdefault(record.field_id, []).append(record)

    field_order = [field['id'] for field in config['fields']]
    data_dir = report_cfg['data_dir']
    subtitle_parts = [
        f'Data directory: {html.escape(str(data_dir))}',
        f'Mode: {html.escape(mode)}',
    ]
    if mode == 'forecast':
        subtitle_parts.append(f'Lookback days: {lookback_days}')

    sections: list[str] = []
    for field_id in field_order:
        field_records = grouped.get(field_id, [])
        if not field_records:
            continue
        field_label = html.escape(field_records[0].field_label)
        tiles: list[str] = []
        for record in field_records:
            rel_image = html.escape(record.image_path.relative_to(output_dir).as_posix())
            rel_thumb = html.escape(record.thumb_path.relative_to(output_dir).as_posix())
            extent_label = html.escape(record.extent_label)
            plot_title = html.escape(record.title)
            tiles.append(
                f'<figure class="contour-tile">'
                f'<a href="{rel_image}" target="_blank" rel="noopener">'
                f'<img src="{rel_thumb}" alt="{plot_title}">'
                f'</a>'
                f'<figcaption>{extent_label}</figcaption>'
                f'</figure>'
            )
        sections.append(
            f'<section class="field-section">'
            f'<h2>{field_label}</h2>'
            f'<div class="contour-grid">{"".join(tiles)}</div>'
            f'</section>'
        )

    rel_map = html.escape(map_path.relative_to(output_dir).as_posix())
    index_path = output_dir / DEFAULT_REPORT_HTML
    index_path.write_text(
        f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    body {{
      font-family: Arial, Helvetica, sans-serif;
      margin: 0;
      padding: 1.5rem;
      color: #222;
      background: #fafafa;
    }}
    h1, h2 {{
      margin: 0 0 0.75rem 0;
    }}
    .subtitle {{
      color: #555;
      margin-bottom: 1.5rem;
    }}
    .map-frame {{
      width: 100%;
      height: 700px;
      border: 1px solid #ccc;
      background: #fff;
      margin-bottom: 2rem;
    }}
    .field-section {{
      margin-bottom: 2rem;
      background: #fff;
      border: 1px solid #ddd;
      border-radius: 6px;
      padding: 1rem;
    }}
    .contour-grid {{
      display: flex;
      flex-wrap: wrap;
      gap: 1rem;
    }}
    .contour-tile {{
      margin: 0;
      width: 320px;
    }}
    .contour-tile img {{
      width: 100%;
      height: auto;
      border: 1px solid #ccc;
      background: #fff;
    }}
    .contour-tile figcaption {{
      margin-top: 0.5rem;
      text-align: center;
      font-size: 0.95rem;
    }}
  </style>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <p class="subtitle">{'<br>'.join(subtitle_parts)}</p>
  <section>
    <h2>Station Hydrographs</h2>
    <iframe class="map-frame" src="{rel_map}" title="Hydrograph station map"></iframe>
  </section>
  {''.join(sections)}
</body>
</html>
""",
        encoding='utf-8',
    )
    return index_path


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description=(
            'Generate an HTML report with station hydrograph map and contour '
            'figure galleries from a TOML configuration file.'
        ),
    )
    parser.add_argument(
        '--config',
        help='Path to report TOML configuration file',
    )
    parser.add_argument(
        '--write-sample-config',
        nargs='?',
        const=DEFAULT_SAMPLE_CONFIG_NAME,
        metavar='PATH',
        help=(
            'Write a sample TOML config file to PATH (default: '
            f'{DEFAULT_SAMPLE_CONFIG_NAME}) and exit, without running the report'
        ),
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite PATH if it already exists (used with --write-sample-config)',
    )
    parser.add_argument(
        '--mode',
        choices=['standard', 'forecast'],
        default=None,
        help=(
            'Report mode. forecast: build hydrographs from analysis lookback '
            f'plus forecast fort.61.nc (default from TOML or {DEFAULT_MODE})'
        ),
    )
    parser.add_argument(
        '--lookback-days',
        type=int,
        default=None,
        help=(
            'Days of analysis fort.61.nc lookback for forecast mode '
            f'(default: {DEFAULT_LOOKBACK_DAYS})'
        ),
    )
    parser.add_argument(
        '--skip-existing',
        action='store_true',
        help='Skip hydrograph and contour figures that already exist',
    )
    parser.add_argument(
        '--skip-on-error',
        action='store_true',
        help='Continue hydrograph plotting if an individual station fails',
    )
    parser.add_argument(
        '--map-thumb-scale',
        type=float,
        default=None,
        help='Hover thumbnail scale for hydrograph map markers',
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
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()

    if args.write_sample_config:
        out_path = write_sample_config(args.write_sample_config, force=args.force)
        print(f'Wrote sample config to {out_path}')
        return 0

    if not args.config:
        get_parser().error('--config is required (or use --write-sample-config)')

    config = load_report_config(args.config)
    report_cfg = config['report']
    output_dir = Path(report_cfg['output_dir'])
    output_dir.mkdir(parents=True, exist_ok=True)

    mode = args.mode or report_cfg.get('mode', DEFAULT_MODE)
    lookback_days = (
        args.lookback_days
        if args.lookback_days is not None
        else report_cfg.get('lookback_days', DEFAULT_LOOKBACK_DAYS)
    )
    map_thumb_scale = (
        args.map_thumb_scale
        if args.map_thumb_scale is not None
        else report_cfg.get('map_thumb_scale', DEFAULT_MAP_THUMB_SCALE)
    )

    contrail_options = _contrail_options_from_config(config['hydrographs'], args)

    written, station_records = generate_hydrographs(
        config,
        output_dir,
        mode=mode,
        lookback_days=lookback_days,
        skip_existing=args.skip_existing,
        skip_on_error=args.skip_on_error,
        contrail_options=contrail_options,
    )
    print(f'Wrote {len(written)} hydrograph figure(s) to {output_dir / "hydrographs"}')

    contour_records = generate_contour_figures(
        config,
        output_dir,
        skip_existing=args.skip_existing,
    )
    print(f'Wrote {len(contour_records)} contour figure(s) to {output_dir / "contours"}')

    index_path = assemble_report_html(
        config,
        station_records,
        contour_records,
        output_dir,
        mode=mode,
        lookback_days=lookback_days,
        map_thumb_scale=map_thumb_scale,
    )
    print(f'Wrote report to {index_path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
