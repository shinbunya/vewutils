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
    filter_stations_by_extent,
    generate_hydrograph_map_html,
    load_station_ids_file,
    plot_f61_hydrographs_from_elev_stat,
    read_elev_stat_stations,
    read_netcdf_time_series,
    resolve_plot_date_range,
)
from vewutils.plot.plot_max_ele_2d import plot_max_ele_2d
from vewutils.plot.plot_solution_2d import plot_solutions_2d

DEFAULT_LOOKBACK_DAYS = 5
DEFAULT_MODE = 'standard'
DEFAULT_REPORT_HTML = 'index.html'
DEFAULT_THUMB_WIDTH = 300
DEFAULT_SAMPLE_CONFIG_NAME = 'sample_report_config.toml'
# Matches the forecast cycle date/hour (e.g. "gfs_2026-07-19-00_NCSCv3") embedded in
# a forecast fort.61.nc's rundes/description attribute.
CYCLE_DATE_HOUR_RE = re.compile(r'(\d{4})-(\d{2})-(\d{2})-(\d{2})')

SAMPLE_CONFIG_HEADER = '''\
# Sample vewutils generate-report configuration.
#
# Fill in the paths marked "REQUIRED" for your simulation, then run:
#   vewutils plot generate-report --config this_file.toml

'''

# [report] section content is mode-dependent: forecast mode needs data_dir_archive
# and lookback_days, while standard mode doesn't use them.
SAMPLE_REPORT_SECTIONS = {
    'standard': '''\
[report]
# REQUIRED. Directory containing the simulation's fort.61.nc / maxele NetCDF files.
data_dir = "/path/to/simulation_cycle/ADCIRC/simulation"

# Where the report (HTML, hydrograph figures, contour figures) is written.
# Defaults to "<data_dir>/report" if omitted.
# output_dir = "/path/to/output/report"

title = "ADCIRC Simulation Report"

# "standard" uses hydrographs.f61or63files (or fort.61.nc in data_dir) as-is.
# "forecast" concatenates lookback_days of prior analysis fort.61.nc files (found
# under data_dir_archive) with data_dir's own forecast fort.61.nc. The forecast
# cycle date/hour is read from that fort.61.nc's rundes/description attribute.
mode = "standard"

# Hover-thumbnail scale for markers on the hydrograph station map.
map_thumb_scale = 1.0

# Skip hydrograph/contour figures that already exist, and continue past a
# failed station instead of aborting. Both are also settable per-run via
# --skip-existing/--skip-on-error, which take precedence over these.
# skip_existing = false
# skip_on_error = false
''',
    'forecast': '''\
[report]
# REQUIRED. Directory containing the forecast simulation's fort.61.nc / maxele
# NetCDF files. This need not live under the archive tree.
data_dir = "/path/to/simulation_cycle/ADCIRC/simulation"

# Where the report (HTML, hydrograph figures, contour figures) is written.
# Defaults to "<data_dir>/report" if omitted.
# output_dir = "/path/to/output/report"

title = "ADCIRC Simulation Report"

# "forecast" concatenates lookback_days of prior analysis fort.61.nc files (found
# under data_dir_archive) with data_dir's own forecast fort.61.nc. The forecast
# cycle date/hour is read from that fort.61.nc's rundes/description attribute.
mode = "forecast"

# REQUIRED for mode = "forecast". Root of the archive tree
# ("<data_dir_archive>/YYYYMMDD/hour_HH/adcirc/analysis/fort.61.nc") used to look
# back for analysis cycles.
data_dir_archive = "/path/to/archive"

# Days of analysis fort.61.nc lookback to concatenate before the forecast cycle.
lookback_days = 5

# Hover-thumbnail scale for markers on the hydrograph station map.
map_thumb_scale = 1.0

# Skip hydrograph/contour figures that already exist, and continue past a
# failed station instead of aborting. Both are also settable per-run via
# --skip-existing/--skip-on-error, which take precedence over these.
# skip_existing = false
# skip_on_error = false
''',
}

SAMPLE_CONFIG_REST = '''\

[hydrographs]
# REQUIRED. Path to the ADCIRC elev_stat.151 station list.
elev_stat = "/path/to/elev_stat.151"

# Restrict which elev_stat owner labels are plotted (default: NOAA/NOS, USGS, NCEM).
# owners = ["NOAA/NOS", "USGS"]

# Optional explicit list/file of station ids to plot instead of all stations.
# station_ids = ["8651370", "8652587"]
# station_ids_file = "/path/to/station_ids.txt"

# Restrict to stations whose fort.61.nc (x, y) falls within these extent_presets
# id(s) bounding box; stations outside all of them are excluded. A single id or
# a list is accepted; a list is a union (station kept if inside any one of
# them). Combines with station_ids above (intersection) if both are set.
# extent_presets = ["region_a"]

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

# Optional: a `timesteps` list treats the field as a time-series *.63.nc file
# (plotted with plot_solutions_2d, one figure per extent per time step) instead
# of a single-snapshot maxele-style file. 0-based indices; -1 means the last
# time step. Each entry gets its own figure, titled/labeled "time step N" (or
# "last time step" for -1).
# [[fields]]
# id = "dynamic_water_level_correction"
# label = "Dynamic Water Level Correction"
# file = "dynamicWaterlevelCorrection.63.nc"
# variable = "dynamicWaterlevelCorrection"
# cmap = "bwr"
# timesteps = [0, -1]          # first and last time steps
# extent_presets = ["region_a"]


# Optional. Only used by `vewutils plot watch-reports`, which scans for
# newly-completed cycles under [watch].pattern and runs generate-report on
# each one automatically (meant to be invoked periodically by cron). See
# `vewutils plot watch-reports --help`.
# [watch]
# pattern = "/path/to/archive/20??????/hour_??/adcirc/forecast/forecast_base"
# max_age_days = 3        # ignore cycles older than this many days from now
# quiet_seconds = 60       # require this much idle time since the newest file
# max_attempts = 5         # give up on a cycle after this many failed attempts
# backoff_minutes = [0, 30, 60, 120, 240]   # wait between successive retries

# Alternative to max_age_days: an absolute cutoff (fixed watermark, doesn't
# slide with time) instead of a rolling window. Takes precedence over
# max_age_days when set.
# min_cycle_date = "2026-07-20"   # YYYY-MM-DD
# min_cycle_hour = 0               # 0-23, default 0

# Optional. If set, watch-reports uploads each finished report/ directory via
# SFTP after a successful run. Auth is via SSH private key only (no
# passwords). remote_path_template's YYYY/YYYYMMDDHH tokens are substituted
# with the cycle's own date/hour (not the upload time).
# [sftp]
# host = "sftp.example.org"
# port = 22
# username = "..."
# key_path = "~/.ssh/id_rsa"
# remote_path_template = "/remote/base/YYYY/YYYYMMDDHH/"
'''


def sample_report_config_toml(mode: str = DEFAULT_MODE) -> str:
    """Build sample report TOML text reflecting the given mode's [report] fields."""
    try:
        report_section = SAMPLE_REPORT_SECTIONS[mode]
    except KeyError:
        valid = ', '.join(sorted(SAMPLE_REPORT_SECTIONS))
        raise ValueError(f'Unknown mode {mode!r}; expected one of: {valid}') from None
    return SAMPLE_CONFIG_HEADER + report_section + SAMPLE_CONFIG_REST


# Keys copied from an extent preset or inline extent table into resolved extents.
EXTENT_OPTION_KEYS = (
    'label', 'title', 'xmin', 'xmax', 'ymin', 'ymax',
    'vmin', 'vmax', 'cbar_label', 'cbar_increment', 'cbar_ticks_increment',
    'figsizex', 'figsizey',
)


@dataclass(frozen=True)
class ArchiveContext:
    data_dir: Path
    archive_root: Path
    target_date: datetime
    target_hour: int


@dataclass(frozen=True)
class CycleSummary:
    """Best-effort cycle/time-window info shown in the report subtitle."""
    cycle_datetime: datetime | None = None
    analysis_window: tuple[datetime, datetime] | None = None
    forecast_window: tuple[datetime, datetime] | None = None


@dataclass(frozen=True)
class ContourRecord:
    field_id: str
    field_label: str
    extent_id: str
    extent_label: str
    image_path: Path
    thumb_path: Path
    title: str
    # Grouping key for the HTML gallery: records sharing one extent_group_id
    # are kept together as a row, with a break before the next group. Equal to
    # extent_id for a plain field (each extent its own group). For a timesteps
    # field, it's the time step's slug instead, so one row holds that time
    # step's tiles across all extents.
    extent_group_id: str = ''

    def __post_init__(self):
        if not self.extent_group_id:
            object.__setattr__(self, 'extent_group_id', self.extent_id)


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


def write_sample_config(
        path: str | Path,
        *,
        mode: str = DEFAULT_MODE,
        force: bool = False) -> Path:
    """Write a sample report TOML config reflecting mode to path.

    Refuses to overwrite an existing file unless force is set.
    """
    path = Path(path)
    if path.exists() and not force:
        raise FileExistsError(f'{path} already exists; use --force to overwrite')
    path.write_text(sample_report_config_toml(mode), encoding='utf-8')
    return path


def load_report_config(
        path: str | Path,
        *,
        data_dir_override: str | Path | None = None) -> dict[str, Any]:
    """Load and validate a report TOML configuration file.

    data_dir_override, if given, replaces report.data_dir (e.g. from a
    watcher applying one template config to many cycle directories) before
    validation, so it goes through the same resolve/is_dir() check as the
    TOML value and output_dir's default (data_dir / 'report') re-derives
    against it automatically.
    """
    path = Path(path)
    with open(path, 'rb') as f:
        raw = tomllib.load(f)

    report = raw.get('report')
    if not isinstance(report, dict):
        raise ValueError(f'{path}: [report] section is required')

    data_dir = data_dir_override or report.get('data_dir')
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

    data_dir_archive = report.get('data_dir_archive')
    if data_dir_archive:
        data_dir_archive = Path(data_dir_archive).resolve()

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

    # [sftp] and [watch] are only used by `watch-reports`; generate-report
    # itself ignores them, so they're just passed through unvalidated here.
    sftp = raw.get('sftp')
    if sftp is not None and not isinstance(sftp, dict):
        raise ValueError(f'{path}: [sftp] must be a table')

    watch = raw.get('watch')
    if watch is not None and not isinstance(watch, dict):
        raise ValueError(f'{path}: [watch] must be a table')

    return {
        'report': {
            **report,
            'data_dir': data_dir,
            'output_dir': output_dir,
            'data_dir_archive': data_dir_archive,
        },
        'hydrographs': hydrographs,
        'contours': dict(contours or {}),
        'extent_presets': extent_presets,
        'fields': validated_fields,
        'sftp': dict(sftp or {}),
        'watch': dict(watch or {}),
    }


def _read_forecast_cycle(f61_path: Path) -> tuple[datetime, int]:
    """Extract forecast cycle date/hour from a fort.61.nc rundes/description attribute."""
    import xarray as xr

    with xr.open_dataset(f61_path) as ds:
        rundes = ds.attrs.get('rundes') or ds.attrs.get('description')

    if not rundes:
        raise ValueError(
            f'{f61_path}: no rundes/description attribute to determine forecast cycle'
        )
    match = CYCLE_DATE_HOUR_RE.search(str(rundes))
    if not match:
        raise ValueError(
            f'{f61_path}: could not parse cycle date/hour (YYYY-MM-DD-HH) from '
            f'rundes/description {rundes!r}'
        )
    year, month, day, hour = match.groups()
    return datetime(int(year), int(month), int(day)), int(hour)


def _time_bounds_across_files(paths: list[str]) -> tuple[datetime, datetime] | None:
    """Earliest/latest exact timestamps across the given files' time axes.

    Best-effort: files that fail to open or have no/empty time axis are
    skipped rather than raising. Returns None if none of the files yielded a
    usable time axis.
    """
    t_min = t_max = None
    for path in paths:
        try:
            times = read_netcdf_time_series(path)
        except (ValueError, OSError, KeyError):
            continue
        if len(times) == 0:
            continue
        file_min = times.min().to_pydatetime()
        file_max = times.max().to_pydatetime()
        t_min = file_min if t_min is None else min(t_min, file_min)
        t_max = file_max if t_max is None else max(t_max, file_max)
    if t_min is None:
        return None
    return t_min, t_max


def compute_cycle_summary(
        config: dict[str, Any],
        *,
        mode: str,
        f61or63files: list[str]) -> CycleSummary:
    """Best-effort cycle date/hour and analysis/forecast time windows.

    f61or63files is the flat (pre-concat) file list already resolved for
    hydrograph generation, so this doesn't re-scan the archive tree; in
    forecast mode the forecast file is always last (see
    discover_forecast_f61_files), everything before it is the analysis
    lookback.
    """
    data_dir = config['report']['data_dir']
    cycle_datetime = None
    forecast_f61 = data_dir / 'fort.61.nc'
    if forecast_f61.is_file():
        try:
            target_date, target_hour = _read_forecast_cycle(forecast_f61)
            cycle_datetime = target_date.replace(hour=target_hour)
        except (ValueError, OSError):
            cycle_datetime = None

    analysis_window = None
    forecast_window = None
    if mode == 'forecast' and f61or63files:
        analysis_window = _time_bounds_across_files(f61or63files[:-1])
        forecast_window = _time_bounds_across_files(f61or63files[-1:])

    return CycleSummary(
        cycle_datetime=cycle_datetime,
        analysis_window=analysis_window,
        forecast_window=forecast_window,
    )


def resolve_forecast_archive_context(
        data_dir: str | Path,
        archive_root: str | Path) -> ArchiveContext:
    """Build a forecast ArchiveContext from a simulation data_dir and archive root.

    The forecast cycle date/hour is read from the rundes/description attribute of
    data_dir's fort.61.nc, since a forecast simulation directory need not live under
    the archive tree the way analysis cycles do.
    """
    data_dir = Path(data_dir).resolve()
    archive_root = Path(archive_root).resolve()
    if not archive_root.is_dir():
        raise ValueError(f'report.data_dir_archive is not a directory: {archive_root}')

    forecast_f61 = data_dir / 'fort.61.nc'
    if not forecast_f61.is_file():
        raise ValueError(f'Forecast fort.61.nc not found: {forecast_f61}')

    target_date, target_hour = _read_forecast_cycle(forecast_f61)
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
        archive_root = config['report'].get('data_dir_archive')
        if not archive_root:
            raise ValueError(
                'report.data_dir_archive is required when mode is "forecast"'
            )
        ctx = resolve_forecast_archive_context(data_dir, archive_root)
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
        contrail_options: dict[str, str] | None) -> tuple[list[Path], list[dict[str, Any]], CycleSummary]:
    """Plot station hydrographs and return written paths, station records, and cycle summary."""
    hydrographs = config['hydrographs']
    hydrograph_dir = output_dir / 'hydrographs'
    hydrograph_dir.mkdir(parents=True, exist_ok=True)

    f61or63files = _resolve_f61or63files(
        config,
        mode=mode,
        lookback_days=lookback_days,
    )
    cycle_summary = compute_cycle_summary(config, mode=mode, f61or63files=f61or63files)
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

    extent_preset_ids = hydrographs.get('extent_presets')
    if extent_preset_ids:
        if isinstance(extent_preset_ids, str):
            extent_preset_ids = [extent_preset_ids]
        elif not isinstance(extent_preset_ids, list):
            raise ValueError(
                'hydrographs.extent_presets must be a string or a list of '
                'extent preset ids'
            )

        extent_presets = config['extent_presets']
        all_stations = read_elev_stat_stations(hydrographs['elev_stat'], owners=owners)
        coord_source = f61or63files[0]
        extent_station_ids: list[str] = []
        seen_station_ids: set[str] = set()
        for extent_preset_id in extent_preset_ids:
            try:
                extent = extent_presets[extent_preset_id]
            except KeyError:
                valid = ', '.join(sorted(extent_presets))
                raise ValueError(
                    f'hydrographs.extent_presets references unknown extent preset '
                    f'{extent_preset_id!r}. Available presets: {valid or "(none)"}'
                ) from None
            # Union: a station is kept if it falls within any listed preset.
            for station in filter_stations_by_extent(all_stations, extent, coord_source):
                if station['station_id'] not in seen_station_ids:
                    seen_station_ids.add(station['station_id'])
                    extent_station_ids.append(station['station_id'])

        station_ids = (
            [sid for sid in station_ids if sid in extent_station_ids]
            if station_ids
            else extent_station_ids
        )

    station_ids = station_ids or None

    f61or63concat = hydrographs.get('f61or63concat', mode == 'forecast')
    if f61or63concat:
        if mode == 'forecast' and len(f61or63files) > 1:
            # Keep analysis (nowcast) cycles as one connected series, separate
            # from the forecast file, so nowcast_forecast_style draws the
            # nowcast cycles solid and only the forecast series dashed.
            f61or63files = [f61or63files[:-1], f61or63files[-1]]
        else:
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
    return written, station_records, cycle_summary


def _field_plot_kwargs(
        field: dict[str, Any],
        extent: dict[str, Any],
        contours_cfg: dict[str, Any]) -> dict[str, Any]:
    """Build plot kwargs shared by plot_max_ele_2d and plot_solutions_2d.

    Does not include the file path (plot_max_ele_2d's maxele_file vs
    plot_solutions_2d's solution_file/timestep), which the caller adds since the
    two functions take it differently.
    """
    title = extent.get('title')
    if not title:
        field_label = field.get('label', field['id'])
        extent_label = extent.get('label', extent['id'])
        title = f'{field_label} — {extent_label}'

    track_file = field.get('track_file', contours_cfg.get('track_file'))
    if track_file and not Path(track_file).is_absolute():
        track_file = str(Path(track_file).resolve())

    kwargs: dict[str, Any] = {
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


def _write_thumbnail(image_path: Path, thumb_path: Path, thumb_width: int) -> None:
    import matplotlib.image as mpimg

    img_array = mpimg.imread(image_path)
    width = img_array.shape[1]
    step = max(1, width // thumb_width)
    thumb_array = img_array[::step, ::step]
    mpimg.imsave(thumb_path, thumb_array)


def _emit_contour_figure(
        *,
        image_path: Path,
        thumb_path: Path,
        field_id: str,
        field_label: str,
        extent_id: str,
        extent_label: str,
        title: str,
        figsize: tuple[float, float],
        dpi: int,
        thumb_width: int,
        skip_existing: bool,
        plot_fn,
        extent_group_id: str = '') -> ContourRecord:
    """Render one contour PNG + thumbnail via plot_fn(fig, ax), or reuse existing.

    plot_fn must return True on success, like plot_max_ele_2d/plot_solutions_2d.
    """
    import matplotlib.pyplot as plt

    if skip_existing and image_path.is_file() and thumb_path.is_file():
        print(f'Skipping existing contour figure: {image_path.name}')
    else:
        print(f'Creating contour figure: {image_path.name}')
        fig, ax = plt.subplots(figsize=figsize)
        success = plot_fn(fig, ax)
        if not success:
            plt.close(fig)
            raise RuntimeError(
                f'Failed to create contour plot for {field_id}/{extent_id}'
            )
        fig.savefig(image_path, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        _write_thumbnail(image_path, thumb_path, thumb_width)

    return ContourRecord(
        field_id=field_id,
        field_label=field_label,
        extent_id=extent_id,
        extent_label=extent_label,
        image_path=image_path,
        thumb_path=thumb_path,
        title=title,
        extent_group_id=extent_group_id,
    )


def generate_contour_figures(
        config: dict[str, Any],
        output_dir: Path,
        *,
        skip_existing: bool) -> list[ContourRecord]:
    """Generate contour PNGs and thumbnails for all configured fields/extents.

    A field with a `timesteps` list (0-based indices; -1 for the last time step)
    is treated as a time-series *.63.nc file and plotted with plot_solutions_2d
    (one figure per extent per time step) instead of plot_max_ele_2d.
    """
    import matplotlib.pyplot as plt

    data_dir = config['report']['data_dir']
    contours_cfg = config['contours']
    contour_dir = output_dir / 'contours'
    contour_dir.mkdir(parents=True, exist_ok=True)

    dpi = contours_cfg.get('dpi', 300)
    thumb_width = contours_cfg.get('thumb_width', DEFAULT_THUMB_WIDTH)

    records: list[ContourRecord] = []

    for field in config['fields']:
        netcdf_path = _resolve_netcdf_path(data_dir, field['file'])
        if not netcdf_path.is_file():
            raise ValueError(f'Contour NetCDF file not found: {netcdf_path}')

        timesteps = field.get('timesteps')
        if timesteps is not None and (
            not isinstance(timesteps, list)
            or not timesteps
            or not all(isinstance(t, int) for t in timesteps)
        ):
            raise ValueError(
                f"fields[{field['id']!r}].timesteps must be a non-empty array of "
                'integers (0-based; -1 for the last time step)'
            )

        field_label = field.get('label', field['id'])

        if timesteps:
            # Timestep is the outer loop so each time step's tiles (across all
            # extents) land together as one row, separated from other time
            # steps' rows -- rather than grouping by extent.
            for timestep in timesteps:
                step_slug = 'last' if timestep == -1 else f't{timestep}'
                step_phrase = 'last time step' if timestep == -1 else f'time step {timestep}'

                for extent in field['extents']:
                    extent_id = extent['id']
                    extent_label = extent.get('label', extent_id)
                    plot_kwargs = _field_plot_kwargs(field, extent, contours_cfg)
                    figsizex = extent.get(
                        'figsizex', field.get('figsizex', contours_cfg.get('figsizex', 12.0))
                    )
                    figsizey = extent.get(
                        'figsizey', field.get('figsizey', contours_cfg.get('figsizey', 10.0))
                    )

                    solution_kwargs = dict(plot_kwargs)
                    solution_kwargs.pop('track_annotate_category_inside_circle', None)
                    solution_kwargs.pop('track_annotate_non_hurricane_inside_circle', None)
                    base_title = solution_kwargs.pop('title')

                    step_id = f'{extent_id}_{step_slug}'
                    step_title = f'{base_title} ({step_phrase})'
                    records.append(_emit_contour_figure(
                        image_path=contour_dir / f"{field['id']}_{step_id}.png",
                        thumb_path=contour_dir / f"{field['id']}_{step_id}_thumb.png",
                        field_id=field['id'],
                        field_label=field_label,
                        extent_id=step_id,
                        extent_label=f'{extent_label} ({step_phrase})',
                        extent_group_id=step_slug,
                        title=step_title,
                        figsize=(figsizex, figsizey),
                        dpi=dpi,
                        thumb_width=thumb_width,
                        skip_existing=skip_existing,
                        plot_fn=lambda fig, ax, ts=timestep, kw=solution_kwargs, title=step_title: (
                            plot_solutions_2d(
                                fig, ax, str(netcdf_path), ts, title=title, **kw
                            )
                        ),
                    ))
        else:
            for extent in field['extents']:
                extent_id = extent['id']
                extent_label = extent.get('label', extent_id)
                plot_kwargs = _field_plot_kwargs(field, extent, contours_cfg)
                figsizex = extent.get(
                    'figsizex', field.get('figsizex', contours_cfg.get('figsizex', 12.0))
                )
                figsizey = extent.get(
                    'figsizey', field.get('figsizey', contours_cfg.get('figsizey', 10.0))
                )
                records.append(_emit_contour_figure(
                    image_path=contour_dir / f"{field['id']}_{extent_id}.png",
                    thumb_path=contour_dir / f"{field['id']}_{extent_id}_thumb.png",
                    field_id=field['id'],
                    field_label=field_label,
                    extent_id=extent_id,
                    extent_label=extent_label,
                    title=plot_kwargs['title'],
                    figsize=(figsizex, figsizey),
                    dpi=dpi,
                    thumb_width=thumb_width,
                    skip_existing=skip_existing,
                    plot_fn=lambda fig, ax, kw=plot_kwargs: plot_max_ele_2d(
                        fig, ax, maxele_file=str(netcdf_path), **kw
                    ),
                ))

    plt.close('all')
    return records


def _html_id(raw: str) -> str:
    """Sanitize an arbitrary string (e.g. a field id) for use as an HTML id."""
    return re.sub(r'[^a-zA-Z0-9_-]+', '-', raw).strip('-') or 'section'


def _format_window(window: tuple[datetime, datetime]) -> str:
    start, end = window
    return f"{start.strftime('%Y-%m-%d %H:%M')} to {end.strftime('%Y-%m-%d %H:%M')} UTC"


def assemble_report_html(
        config: dict[str, Any],
        station_records: list[dict[str, Any]],
        contour_records: list[ContourRecord],
        output_dir: Path,
        *,
        cycle_summary: CycleSummary,
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
    generated_at = datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')
    subtitle_parts = []
    if cycle_summary.cycle_datetime:
        subtitle_parts.append(
            f"Cycle: {html.escape(cycle_summary.cycle_datetime.strftime('%Y-%m-%d %H'))} UTC"
        )
    if cycle_summary.analysis_window:
        subtitle_parts.append(f'Analysis: {html.escape(_format_window(cycle_summary.analysis_window))}')
    if cycle_summary.forecast_window:
        subtitle_parts.append(f'Forecast: {html.escape(_format_window(cycle_summary.forecast_window))}')
    subtitle_parts.append(f'Generated: {html.escape(generated_at)}')

    nav_links: list[str] = ['<a href="#section-hydrographs">Station Hydrographs</a>']
    sections: list[str] = []
    for field_id in field_order:
        field_records = grouped.get(field_id, [])
        if not field_records:
            continue
        field_label = html.escape(field_records[0].field_label)
        section_id = f'section-{_html_id(field_id)}'
        nav_links.append(f'<a href="#{section_id}">{field_label}</a>')
        # Multiple time steps produce several tiles per extent (same
        # extent_group_id); break between extent groups so they don't run
        # together in the wrapping grid.
        group_ids = [record.extent_group_id for record in field_records]
        grouped_by_extent = len(set(group_ids)) < len(group_ids)
        tiles: list[str] = []
        previous_group_id: str | None = None
        for record in field_records:
            if (
                grouped_by_extent
                and previous_group_id is not None
                and record.extent_group_id != previous_group_id
            ):
                tiles.append('<div class="contour-break"></div>')
            previous_group_id = record.extent_group_id

            rel_image = html.escape(record.image_path.relative_to(output_dir).as_posix())
            rel_thumb = html.escape(record.thumb_path.relative_to(output_dir).as_posix())
            extent_label = html.escape(record.extent_label)
            plot_title = html.escape(record.title)
            tiles.append(
                f'<figure class="contour-tile">'
                f'<img src="{rel_thumb}" alt="{plot_title}" class="lightbox-trigger" '
                f'data-full="{rel_image}" data-caption="{plot_title}" '
                f'tabindex="0" role="button">'
                f'<figcaption>{extent_label}</figcaption>'
                f'</figure>'
            )
        sections.append(
            f'<details class="field-section" id="{section_id}" open>'
            f'<summary><span class="chevron" aria-hidden="true"></span>'
            f'<h2>{field_label}</h2></summary>'
            f'<div class="contour-grid">{"".join(tiles)}</div>'
            f'</details>'
        )

    rel_map = html.escape(map_path.relative_to(output_dir).as_posix())
    index_path = output_dir / DEFAULT_REPORT_HTML
    index_path.write_text(
        f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <link rel="icon" href="data:image/svg+xml,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><text y='.9em' font-size='90'>🌊</text></svg>">
  <style>
    :root {{
      --accent: #1b4965;
      --accent-light: #e8f1f5;
      --border: #dde3e7;
      --bg: #f5f7f8;
      --text: #1f2933;
      --text-muted: #5b6b74;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      margin: 0;
      color: var(--text);
      background: var(--bg);
      line-height: 1.5;
    }}
    a {{ color: var(--accent); }}
    .site-header {{
      position: sticky;
      top: 0;
      z-index: 10;
      background: #fff;
      border-bottom: 1px solid var(--border);
      box-shadow: 0 1px 4px rgba(0,0,0,0.06);
    }}
    .site-header nav {{
      display: flex;
      flex-wrap: wrap;
      gap: 0.25rem 1rem;
      padding: 0.6rem 1.5rem;
      max-width: 1500px;
      margin: 0 auto;
    }}
    .site-header nav a {{
      text-decoration: none;
      color: var(--text-muted);
      font-size: 0.9rem;
      font-weight: 600;
      padding: 0.25rem 0;
      border-bottom: 2px solid transparent;
    }}
    .site-header nav a:hover {{
      color: var(--accent);
      border-bottom-color: var(--accent);
    }}
    .page {{
      max-width: 1500px;
      margin: 0 auto;
      padding: 1.5rem;
    }}
    .hero h1 {{
      margin: 0 0 0.4rem 0;
      font-size: 1.6rem;
      color: var(--accent);
    }}
    .subtitle {{
      color: var(--text-muted);
      margin: 0 0 1.5rem 0;
      font-size: 0.9rem;
    }}
    .map-frame {{
      width: 100%;
      height: 700px;
      border: 1px solid var(--border);
      background: #fff;
      border-radius: 0 0 8px 8px;
    }}
    .field-section {{
      margin-bottom: 1.5rem;
      background: #fff;
      border: 1px solid var(--border);
      border-radius: 8px;
      overflow: hidden;
    }}
    .field-section summary {{
      display: flex;
      align-items: center;
      gap: 0.6rem;
      cursor: pointer;
      padding: 0.9rem 1.25rem;
      background: var(--accent-light);
      list-style: none;
      user-select: none;
    }}
    .field-section summary::-webkit-details-marker {{ display: none; }}
    .field-section summary h2 {{
      margin: 0;
      font-size: 1.05rem;
      color: var(--accent);
    }}
    .field-section .chevron {{
      width: 0.6rem;
      height: 0.6rem;
      border-right: 2px solid var(--accent);
      border-bottom: 2px solid var(--accent);
      transform: rotate(-45deg);
      transition: transform 0.15s ease;
      flex: none;
    }}
    .field-section[open] > summary .chevron {{
      transform: rotate(45deg);
    }}
    .field-section > .contour-grid,
    .field-section > .map-frame {{
      padding: 1.25rem;
    }}
    .contour-grid {{
      display: flex;
      flex-wrap: wrap;
      gap: 1rem;
    }}
    .contour-break {{
      flex-basis: 100%;
      height: 0;
    }}
    .contour-tile {{
      margin: 0;
    }}
    .contour-tile img {{
      height: 260px;
      width: auto;
      max-width: 100%;
      border: 1px solid var(--border);
      background: #fff;
      cursor: zoom-in;
      display: block;
    }}
    .contour-tile figcaption {{
      margin-top: 0.5rem;
      text-align: center;
      font-size: 0.95rem;
    }}
    footer {{
      max-width: 1500px;
      margin: 0 auto;
      padding: 1rem 1.5rem 2.5rem;
      color: var(--text-muted);
      font-size: 0.85rem;
      display: flex;
      justify-content: space-between;
      gap: 1rem;
      flex-wrap: wrap;
    }}
    .lightbox {{
      position: fixed;
      inset: 0;
      background: rgba(10, 15, 20, 0.9);
      display: none;
      align-items: center;
      justify-content: center;
      z-index: 1000;
      padding: 2rem;
    }}
    .lightbox.open {{ display: flex; }}
    .lightbox img {{
      max-width: 100%;
      max-height: 78vh;
      box-shadow: 0 4px 24px rgba(0,0,0,0.5);
    }}
    .lightbox-content {{
      display: flex;
      flex-direction: column;
      align-items: center;
      gap: 0.75rem;
      max-width: 100%;
    }}
    .lightbox-caption {{
      color: #f0f0f0;
      font-size: 0.95rem;
      text-align: center;
    }}
    .lightbox-caption a {{
      color: #9fd3ff;
      margin-left: 0.75rem;
    }}
    .lightbox-close, .lightbox-prev, .lightbox-next {{
      position: fixed;
      background: rgba(255,255,255,0.1);
      color: #fff;
      border: none;
      cursor: pointer;
      font-size: 1.5rem;
      line-height: 1;
      padding: 0.5rem 0.9rem;
      border-radius: 6px;
    }}
    .lightbox-close:hover, .lightbox-prev:hover, .lightbox-next:hover {{
      background: rgba(255,255,255,0.25);
    }}
    .lightbox-close {{ top: 1rem; right: 1.5rem; }}
    .lightbox-prev {{ left: 1.5rem; top: 50%; transform: translateY(-50%); }}
    .lightbox-next {{ right: 1.5rem; top: 50%; transform: translateY(-50%); }}
  </style>
</head>
<body>
  <header class="site-header">
    <nav>{''.join(nav_links)}</nav>
  </header>
  <div class="page">
    <div class="hero">
      <h1>{html.escape(title)}</h1>
      <p class="subtitle">{' &nbsp;&middot;&nbsp; '.join(subtitle_parts)}</p>
    </div>
    <details class="field-section" id="section-hydrographs" open>
      <summary><span class="chevron" aria-hidden="true"></span><h2>Station Hydrographs</h2></summary>
      <iframe class="map-frame" src="{rel_map}" title="Hydrograph station map"></iframe>
    </details>
    {''.join(sections)}
  </div>
  <footer>
    <span>Generated {html.escape(generated_at)}</span>
    <a href="#top">Back to top</a>
  </footer>

  <div class="lightbox" id="lightbox">
    <button class="lightbox-close" id="lightbox-close" aria-label="Close">&times;</button>
    <button class="lightbox-prev" id="lightbox-prev" aria-label="Previous">&#8249;</button>
    <button class="lightbox-next" id="lightbox-next" aria-label="Next">&#8250;</button>
    <div class="lightbox-content">
      <img id="lightbox-img" src="" alt="">
      <div class="lightbox-caption">
        <span id="lightbox-caption-text"></span>
        <a id="lightbox-full-link" href="#" target="_blank" rel="noopener">Open full size ⤢</a>
      </div>
    </div>
  </div>

  <script>
    (function() {{
      // Keep collapsed/expanded state across reloads (the report regenerates
      // periodically, so this persists a viewer's preference between runs).
      document.querySelectorAll('details.field-section').forEach(function(details) {{
        var key = 'report-section:' + details.id;
        var stored = localStorage.getItem(key);
        if (stored !== null) {{
          details.open = stored === 'true';
        }}
        details.addEventListener('toggle', function() {{
          localStorage.setItem(key, details.open);
        }});
      }});

      // Ensure the target section is expanded when navigated to via a nav
      // link or a direct link with a #fragment.
      function openTargetSection() {{
        if (!location.hash) return;
        var target = document.querySelector(location.hash);
        if (target && target.tagName === 'DETAILS') {{
          target.open = true;
        }}
      }}
      window.addEventListener('hashchange', openTargetSection);
      openTargetSection();

      // Lightbox
      var triggers = Array.prototype.slice.call(document.querySelectorAll('.lightbox-trigger'));
      var lightbox = document.getElementById('lightbox');
      var lightboxImg = document.getElementById('lightbox-img');
      var lightboxCaption = document.getElementById('lightbox-caption-text');
      var lightboxFullLink = document.getElementById('lightbox-full-link');
      var currentIndex = -1;

      function show(index) {{
        if (triggers.length === 0) return;
        currentIndex = (index + triggers.length) % triggers.length;
        var el = triggers[currentIndex];
        lightboxImg.src = el.getAttribute('data-full');
        lightboxImg.alt = el.getAttribute('data-caption') || '';
        lightboxCaption.textContent = el.getAttribute('data-caption') || '';
        lightboxFullLink.href = el.getAttribute('data-full');
        lightbox.classList.add('open');
      }}

      function close() {{
        lightbox.classList.remove('open');
        lightboxImg.src = '';
      }}

      triggers.forEach(function(el, index) {{
        el.addEventListener('click', function() {{ show(index); }});
        el.addEventListener('keydown', function(evt) {{
          if (evt.key === 'Enter' || evt.key === ' ') {{
            evt.preventDefault();
            show(index);
          }}
        }});
      }});

      document.getElementById('lightbox-close').addEventListener('click', close);
      document.getElementById('lightbox-prev').addEventListener('click', function() {{ show(currentIndex - 1); }});
      document.getElementById('lightbox-next').addEventListener('click', function() {{ show(currentIndex + 1); }});
      lightbox.addEventListener('click', function(evt) {{
        if (evt.target === lightbox) close();
      }});
      document.addEventListener('keydown', function(evt) {{
        if (!lightbox.classList.contains('open')) return;
        if (evt.key === 'Escape') close();
        else if (evt.key === 'ArrowLeft') show(currentIndex - 1);
        else if (evt.key === 'ArrowRight') show(currentIndex + 1);
      }});
    }})();
  </script>
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
        '--data-dir',
        help=(
            'Override report.data_dir from the TOML (e.g. to apply one '
            'template config to a different simulation directory). '
            'output_dir still defaults to <data_dir>/report unless set.'
        ),
    )
    parser.add_argument(
        '--write-sample-config',
        nargs='?',
        const=DEFAULT_SAMPLE_CONFIG_NAME,
        metavar='PATH',
        help=(
            'Write a sample TOML config file to PATH (default: '
            f'{DEFAULT_SAMPLE_CONFIG_NAME}) and exit, without running the report. '
            'The [report] section reflects --mode (standard or forecast)'
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
        default=None,
        help=(
            'Skip hydrograph and contour figures that already exist '
            '(overrides report.skip_existing; default: false)'
        ),
    )
    parser.add_argument(
        '--skip-on-error',
        action='store_true',
        default=None,
        help=(
            'Continue hydrograph plotting if an individual station fails '
            '(overrides report.skip_on_error; default: false)'
        ),
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
        out_path = write_sample_config(
            args.write_sample_config,
            mode=args.mode or DEFAULT_MODE,
            force=args.force,
        )
        print(f'Wrote sample config to {out_path} (mode={args.mode or DEFAULT_MODE})')
        return 0

    if not args.config:
        get_parser().error('--config is required (or use --write-sample-config)')

    config = load_report_config(args.config, data_dir_override=args.data_dir)
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
    skip_existing = (
        args.skip_existing if args.skip_existing is not None
        else report_cfg.get('skip_existing', False)
    )
    skip_on_error = (
        args.skip_on_error if args.skip_on_error is not None
        else report_cfg.get('skip_on_error', False)
    )

    contrail_options = _contrail_options_from_config(config['hydrographs'], args)

    written, station_records, cycle_summary = generate_hydrographs(
        config,
        output_dir,
        mode=mode,
        lookback_days=lookback_days,
        skip_existing=skip_existing,
        skip_on_error=skip_on_error,
        contrail_options=contrail_options,
    )
    print(f'Wrote {len(written)} hydrograph figure(s) to {output_dir / "hydrographs"}')

    contour_records = generate_contour_figures(
        config,
        output_dir,
        skip_existing=skip_existing,
    )
    print(f'Wrote {len(contour_records)} contour figure(s) to {output_dir / "contours"}')

    index_path = assemble_report_html(
        config,
        station_records,
        contour_records,
        output_dir,
        cycle_summary=cycle_summary,
        map_thumb_scale=map_thumb_scale,
    )
    print(f'Wrote report to {index_path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
