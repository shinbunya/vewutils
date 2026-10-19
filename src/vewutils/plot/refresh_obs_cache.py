"""Fill the hydrograph observation cache from a report configuration.

Fully spanned days already in the cache are reused. The current incomplete
UTC day is fetched and overwritten on every run. Report generation can then
use ``hydrographs.observations = "cache"`` and read only those files.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timedelta

import pytz

from vewutils.plot.generate_report import (
    DEFAULT_LOOKBACK_DAYS,
    _contrail_options_from_config,
    load_report_config,
)
from vewutils.plot.get_obswl import get_obswl
from vewutils.plot.plot_f61_hydrographs import (
    DEFAULT_ELEV_STAT_OWNERS,
    observation_id_for_station,
    read_elev_stat_stations,
    resolve_hydrograph_station_ids,
    resolve_plot_date_range,
    select_elev_stat_stations_by_ids,
)


def observation_refresh_window(config, *, days=None):
    """Return the UTC window the cache refresh should cover.

    When the report sets both ``hydrographs.date_start`` and ``date_end``,
    those bounds are used. Otherwise the window is ``lookback_days`` plus one
    extra day, ending now, so the earliest analysis calendar day is included.
    ``days`` overrides ``report.lookback_days``.
    """
    hydrographs = config['hydrographs']
    date_start = hydrographs.get('date_start')
    date_end = hydrographs.get('date_end')
    if (date_start is None) != (date_end is None):
        raise ValueError(
            'hydrographs.date_start and hydrographs.date_end must both be set '
            'or both be omitted'
        )
    if date_start is not None and date_end is not None:
        return resolve_plot_date_range(date_start, date_end, [])

    lookback = (
        days if days is not None
        else config['report'].get('lookback_days', DEFAULT_LOOKBACK_DAYS)
    )
    span = int(lookback) + 1
    if span < 1:
        raise ValueError('lookback days must be at least 0')
    now = datetime.now(pytz.UTC)
    start = (now - timedelta(days=span)).replace(
        hour=0, minute=0, second=0, microsecond=0,
    )
    return start, now


def refresh_obs_cache(config, *, days=None, contrail_options=None) -> int:
    """Fetch observations for the report's stations and write the daily cache.

    Returns the number of stations whose fetch failed. A failed fetch does not
    replace an existing cache file for the open day.
    """
    hydrographs = config['hydrographs']
    cache_dir = hydrographs.get('cache_dir')
    if not cache_dir:
        raise ValueError(
            'hydrographs.cache_dir is required to refresh the observation cache'
        )

    owners = hydrographs.get('owners', list(DEFAULT_ELEV_STAT_OWNERS))
    station_ids = resolve_hydrograph_station_ids(
        hydrographs['elev_stat'],
        owners=owners,
        station_ids=hydrographs.get('station_ids'),
        station_ids_file=hydrographs.get('station_ids_file'),
        extent_presets=config.get('extent_presets'),
        extent_preset_ids=hydrographs.get('extent_presets'),
    )
    stations = read_elev_stat_stations(hydrographs['elev_stat'], owners=owners)
    if station_ids:
        stations, missing = select_elev_stat_stations_by_ids(stations, station_ids)
        if missing:
            print(
                'Warning: station id(s) not found in elev_stat (after owner '
                f'filter): {", ".join(missing)}',
                file=sys.stderr,
            )
        if not stations:
            raise ValueError('No stations selected for the observation cache')

    date_start, date_end = observation_refresh_window(config, days=days)
    print(
        f'Refreshing observations for {len(stations)} station(s) '
        f'from {date_start} to {date_end} into {cache_dir}',
        file=sys.stderr,
    )

    datum = hydrographs.get('station_datum', 'NAVD')
    station_id_type = hydrographs.get('station_id_type', 'f61')
    overrides = hydrographs.get('contrail_id_overrides')
    n_contrail = sum(
        1 for station in stations
        if observation_id_for_station(station, overrides)[0] == 'CONTRAIL'
    )
    if n_contrail and not contrail_options:
        raise ValueError(
            'CONTRAIL credentials are required for NCEM stations. '
            'Set CONTRAIL_USERNAME and CONTRAIL_PASSWORD or pass '
            '--username and --password.'
        )

    failures = 0
    for index, station in enumerate(stations, start=1):
        obs_owner, obs_station_id = observation_id_for_station(station, overrides)
        print(
            f'Fetching {index}/{len(stations)}: {station["owner"]} '
            f'{station["station_id"]} ({station["name"]})',
            file=sys.stderr,
        )
        options = None
        if obs_owner == 'CONTRAIL':
            options = dict(contrail_options or {})
            options['station_id_type'] = station_id_type
        try:
            _name, _lon, _lat, obs_time, obs_wl = get_obswl(
                obs_owner,
                obs_station_id,
                date_start,
                date_end,
                datum,
                options,
                cache_dir=cache_dir,
                cache_open_day=True,
            )
        except Exception as exc:
            failures += 1
            print(
                f'Warning: could not refresh {obs_owner} station '
                f'{obs_station_id} ({exc})',
                file=sys.stderr,
            )
            continue
        if obs_time is None or obs_wl is None:
            failures += 1
            print(
                f'Warning: fetch failed for {obs_owner} station {obs_station_id}; '
                'leaving any existing open-day cache in place',
                file=sys.stderr,
            )
    print(
        f'Finished observation cache refresh ({failures} station(s) failed)',
        file=sys.stderr,
    )
    return failures


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description=(
            'Download water-level observations for the stations in a report '
            'config and write them to hydrographs.cache_dir. Fully spanned '
            'days already cached are reused. The current incomplete UTC day '
            'is fetched and overwritten on every run.'
        ),
    )
    parser.add_argument(
        '--config',
        required=True,
        help='Report TOML configuration file',
    )
    parser.add_argument(
        '--days',
        type=int,
        default=None,
        help=(
            'Override report.lookback_days. One extra day is included so the '
            'earliest analysis calendar day is covered. Ignored when '
            'hydrographs.date_start and date_end are both set.'
        ),
    )
    parser.add_argument('--username', help='CONTRAIL username')
    parser.add_argument('--password', help='CONTRAIL password')
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    config = load_report_config(args.config)
    contrail_options = _contrail_options_from_config(config['hydrographs'], args)
    failures = refresh_obs_cache(
        config,
        days=args.days,
        contrail_options=contrail_options,
    )
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
