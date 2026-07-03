# Default line colors for model series; cycled when more series than colors.
DEFAULT_LINE_COLORS = ['b', 'r', 'y', 'm', 'c', 'g']


def _default_line_color(index):
    """Return the color for model series ``index``, cycling through the palette."""
    return DEFAULT_LINE_COLORS[index % len(DEFAULT_LINE_COLORS)]


def _coerce_observation_series(obs_time, obs_wl):
    """Normalize observation arrays; return empty series when data are missing."""
    import numpy as np
    import pandas as pd

    if obs_time is None or obs_wl is None:
        return np.array([]), np.array([], dtype=float), False

    if isinstance(obs_wl, pd.Series):
        obs_wl = obs_wl.to_numpy(dtype=float, na_value=np.nan)
    else:
        obs_wl = np.asarray(obs_wl, dtype=float)

    if obs_wl.size == 0:
        return np.array([]), np.array([], dtype=float), False

    try:
        obs_time = pd.to_datetime(obs_time, utc=True)
    except (TypeError, ValueError):
        return np.array([]), np.array([], dtype=float), False

    if isinstance(obs_time, pd.DatetimeIndex):
        if obs_time.tz is not None:
            obs_time = obs_time.tz_convert('UTC').tz_localize(None)
        obs_time = obs_time.to_numpy()
    else:
        if obs_time.dt.tz is not None:
            obs_time = obs_time.dt.tz_convert('UTC').dt.tz_localize(None)
        obs_time = obs_time.to_numpy()

    if obs_time.size == 0:
        return np.array([]), np.array([], dtype=float), False

    valid = np.isfinite(obs_wl)
    if not np.any(valid):
        return np.array([]), np.array([], dtype=float), False

    return obs_time, obs_wl, True


def plot_hydrograph_at_station(
        fig, ax,
        station_owner, station_id, station_lon, station_lat, station_datum,
        date_start, date_end, 
        f61or63files=None, f61or63starts=None, f61or63labels=None, f61or63colors=None,
        station_id_f61=None,
        # Backward compatibility aliases
        f63files=None, f63starts=None, f63labels=None, f63colors=None,
        f63files_fallback=None, legend_loc='best', legend_loc_rect=None, 
        plot_movingaverage=False, plot_movingaverage_position='backward', 
        plot_in_foot=False, movingaverage_window=0, options=None, cache_dir=None,
        station_id_type=None, adjust_datum_by_mean_error_period_days=0,
        connect=False):
    """Plot observed and modeled water levels at a station.

    Parameters
    ----------
    station_owner : str or None
        Data source for observations: 'NOAA', 'USGS', 'CONTRAIL', 'SECOORA', or None
        to skip observed data.
    station_id_type : str, optional
        For CONTRAIL only: how to interpret ``station_id`` when it is not
        ``contrail_site_id/f61_code``. One of ``'auto'``, ``'contrail'``, ``'f61'``,
        or ``None`` (legacy: same id for observations and fort.61). Combined ids such
        as ``1205/EGHN7`` are always split. See
        :func:`vewutils.plot.get_obswl.resolve_contrail_station_ids`.
    options : dict, optional
        Passed to :func:`vewutils.plot.get_obswl.get_obswl`. For CONTRAIL, must
        include ``username`` and ``password``; may include ``sensor_type``
        (water_elevation, stream_elevation, or stage).
    cache_dir : str or Path, optional
        Directory for caching observation data and the CONTRAIL station list
        (via :func:`vewutils.plot.get_obswl.get_obswl`).
    connect : bool, optional
        When multiple files are concatenated into a single model series, prepend
        the last ``(time, value)`` of the previous file to the next file so the
        plotted line bridges the gap between consecutive files instead of
        breaking. Only affects files concatenated within the same series
        (default: ``False``).
    """
    import os
    import sys
    import requests
    from datetime import datetime, timedelta
    import matplotlib.pyplot as plt
    import pytz
    from matplotlib.dates import DateFormatter
    import pandas as pd
    import numpy as np
    from vewutils.plot.get_obswl import get_obswl, resolve_contrail_station_ids
    from vewutils.plot.get_adcwl import get_adcwl
    
    # Handle backward compatibility: f63* parameters map to f61or63*
    if f63files is not None:
        f61or63files = f63files
    if f63starts is not None:
        f61or63starts = f63starts
    if f63labels is not None:
        f61or63labels = f63labels
    if f63colors is not None:
        f61or63colors = f63colors
    
    # Validate that we have files to process
    if f61or63files is None:
        raise ValueError("Either f61or63files or f63files must be provided")

    display_station_id = station_id
    contrail_site_id = station_id
    contrail_obs_unavailable = False

    # Resolve CONTRAIL site id vs fort.61 station code when needed
    if station_owner == 'CONTRAIL' and station_id is not None:
        opts = options or {}
        needs_resolve = (
            '/' in str(station_id)
            or station_id_type is not None
            or (
                station_id_f61 is not None
                and str(station_id_f61) != str(station_id)
            )
        )
        if needs_resolve:
            try:
                resolved = resolve_contrail_station_ids(
                    station_id,
                    station_id_type=station_id_type,
                    username=opts.get('username'),
                    password=opts.get('password'),
                    cache_dir=cache_dir,
                    f61_station_id=station_id_f61,
                )
                contrail_site_id = resolved['contrail_site_id']
                station_id_f61 = resolved['f61_station_id']
                display_station_id = resolved['display_station_id']
            except Exception as exc:
                contrail_obs_unavailable = True
                station_id_f61 = station_id_f61 or station_id
                display_station_id = station_id
                print(
                    f'Warning: could not resolve CONTRAIL id for fort.61 station '
                    f'{station_id!r} ({exc}); plotting model only.',
                    file=sys.stderr,
                )
        elif station_id_f61 is None:
            station_id_f61 = station_id
    elif station_id_f61 is None:
        station_id_f61 = station_id
        
    # Set unit conversion factor
    if plot_in_foot:
        m2ft = 1.0/0.3048
    else:
        m2ft = 1.0
        
    # Initialization of station_name
    station_name = None
    has_obs_data = False
    obs_time = np.array([])
    obs_wl = np.array([], dtype=float)

    # Get the observed water level data
    if station_owner is not None and not contrail_obs_unavailable:
        if station_owner == 'CONTRAIL':
            opts = options or {}
            if not opts.get('username') or not opts.get('password'):
                raise ValueError(
                    "CONTRAIL requires 'username' and 'password' in options")
        obs_station_id = (
            contrail_site_id if station_owner == 'CONTRAIL' else station_id
        )
        obswl_options = options
        obs_fetch_failed = False
        try:
            station_name, station_lon_, station_lat_, obs_time, obs_wl = \
                get_obswl(
                    station_owner, obs_station_id, date_start, date_end,
                    station_datum, obswl_options, cache_dir=cache_dir)
        except Exception as exc:
            obs_fetch_failed = True
            print(
                f'Warning: could not retrieve observation data for '
                f'{station_owner} station {obs_station_id} ({exc}); '
                f'plotting model only.',
                file=sys.stderr,
            )
            station_lon_ = None
            station_lat_ = None
            obs_time = None
            obs_wl = None
        if station_lon is None and station_lon_ is not None:
            station_lon = station_lon_
            station_lat = station_lat_
        if station_id is None:
            raise ValueError('Station ID is required if station_owner is not NONE')
        if station_lon is None:
            raise ValueError('Station longitude is required if station_owner is not NONE')

        obs_time, obs_wl, has_obs_data = _coerce_observation_series(obs_time, obs_wl)
        if not has_obs_data and not obs_fetch_failed:
            print(
                f'Warning: no observation data for {station_owner} station '
                f'{obs_station_id}; plotting model only.',
                file=sys.stderr,
            )
    
    # Get the forecast water level data
    f61or63_times = []
    f61or63_wls = []
    for i in range(len(f61or63files)):
        if isinstance(f61or63files[i], list):
            f61or63file = f61or63files[i]
            f61or63_start = [None for _ in range(len(f61or63file))]
            if f61or63starts:
                if isinstance(f61or63starts, list):
                    if len(f61or63starts) == len(f61or63file):
                        f61or63_start = f61or63starts[i]
                    else:
                        f61or63_start = f61or63starts[0]
                else:
                    f61or63_start = f61or63starts
        else:
            f61or63file = [f61or63files[i]]
            f61or63_start = [None for _ in range(len(f61or63file))]
            if f61or63starts:
                if isinstance(f61or63starts, list):
                    if len(f61or63starts) == len(f61or63file):
                        f61or63_start = f61or63starts
                    else:
                        f61or63_start = [f61or63starts for _ in range(len(f61or63file))]
                else:
                    f61or63_start = [f61or63starts for _ in range(len(f61or63file))]
        
        f61or63_time = []
        f61or63_wl = []
        print('Processing ADCIRC files', end='')
        for j in range(len(f61or63file)):
            f61or63filej = f61or63file[j]
            f61or63_startj = f61or63_start[j]
            # print('Processing {}...'.format(f61or63filej))
            print('.', end='')
            
            # Determine fallback file if available
            fallback_filej = None
            if f63files_fallback is not None and i < len(f63files_fallback):
                fallback_file = f63files_fallback[i]
                if isinstance(fallback_file, list) and j < len(fallback_file):
                    fallback_filej = fallback_file[j]
                elif isinstance(fallback_file, str):
                    fallback_filej = fallback_file
            
            # Read data using get_adcwl (with automatic fallback if needed)
            f61or63_timej, f61or63_wlj = get_adcwl(
                f61or63filej, 
                station_name=station_id_f61,  # Use station_id for fort.61.nc files
                stx=station_lon, 
                sty=station_lat,
                fallback_file=fallback_filej  # Fallback fort.63.nc file if station not found
            )
            
            if f61or63_startj:
                tdj = f61or63_startj - f61or63_timej[0]
                f61or63_timej = [tdj + t for t in f61or63_timej]
            if f61or63_timej is not None:
                timej_list = (
                    f61or63_timej.tolist()
                    if hasattr(f61or63_timej, 'tolist')
                    else list(f61or63_timej)
                )
                wlj_list = (
                    f61or63_wlj.tolist()
                    if hasattr(f61or63_wlj, 'tolist')
                    else list(f61or63_wlj)
                )
                # Bridge the gap between consecutive files in the same series
                if connect and f61or63_time and timej_list:
                    last_t = f61or63_time[-1]
                    last_v = f61or63_wl[-1]
                    if last_v is not None and not (
                        isinstance(last_v, float) and np.isnan(last_v)
                    ):
                        timej_list = [last_t] + timej_list
                        wlj_list = [last_v] + wlj_list
                f61or63_time.extend(timej_list)
                f61or63_wl.extend(wlj_list)
        f61or63_times.append(f61or63_time)
        f61or63_wl = [wl if abs(wl) <= 100 else np.nan for wl in f61or63_wl]
        f61or63_wls.append(f61or63_wl)
        print(' Done.')

    # Datum adjustment by mean error over first period_in_days: interpolate obs at solution times, compute mean error, subtract from solution
    if adjust_datum_by_mean_error_period_days > 0 and has_obs_data and obs_time.size > 0:
        from matplotlib.dates import date2num
        obs_num = date2num(obs_time)
        obs_wl_valid = np.asarray(obs_wl)
        sort_idx = np.argsort(obs_num)
        obs_num = obs_num[sort_idx]
        obs_wl_valid = obs_wl_valid[sort_idx]
        for i in range(len(f61or63_times)):
            times_i = np.asarray(f61or63_times[i])
            wls_i = np.asarray(f61or63_wls[i], dtype=float)
            t0 = times_i[0]
            cutoff = t0 + timedelta(days=adjust_datum_by_mean_error_period_days)
            mask = times_i <= cutoff
            if not np.any(mask):
                continue
            sol_times_in_window = times_i[mask]
            sol_wl_in_window = wls_i[mask]
            sol_num = date2num(sol_times_in_window)
            obs_interp = np.interp(sol_num, obs_num, obs_wl_valid)
            errors = sol_wl_in_window - obs_interp
            valid = np.isfinite(errors)
            if np.any(valid):
                mean_error = np.mean(errors[valid])
                f61or63_wls[i] = [w - mean_error for w in f61or63_wls[i]]

    # Calculate moving average over a 2-day window
    if plot_movingaverage:
        window_size = timedelta(days=2)
        
        if has_obs_data:
            obs_time_ma = []
            obs_wl_ma = []
            for i in range(len(obs_time)):
                if plot_movingaverage_position == 'backward':
                    start_time = obs_time[i] - window_size
                    end_time = obs_time[i]
                elif plot_movingaverage_position == 'center':
                    start_time = obs_time[i] - window_size / 2
                    end_time = obs_time[i] + window_size / 2
                elif plot_movingaverage_position == 'forward':
                    start_time = obs_time[i]
                    end_time = obs_time[i] + window_size
                else:
                    raise ValueError('Invalid plot_movingaverage_position: {}'.format(plot_movingaverage_position))
                start_i = 0
                end_i = 0
                for start_i in reversed(range(i)):
                    if obs_time[start_i] < start_time:
                        start_i = min(i, start_i + 1)
                        break
                for end_i in range(i, len(obs_time)):
                    if obs_time[end_i] > end_time:
                        end_i = max(i, end_i - 1)
                        break
                wl_window = [obs_wl[k] for k in range(start_i, end_i+1)]
                if wl_window:
                    wl_avg = sum(wl_window) / len(wl_window)
                    obs_time_ma.append(obs_time[i])
                    obs_wl_ma.append(wl_avg)
                
        f61or63_times_ma = []
        f61or63_wls_ma = []
        for i in range(len(f61or63_times)):
            times = f61or63_times[i]
            wls = f61or63_wls[i]
            times_ma = []
            wls_ma = []
            for j in range(len(times)):
                start_time = times[j] - window_size / 2
                end_time = times[j] + window_size / 2
                # wl_window = [wls[k] for k in range(len(times)) if start_time <= times[k] <= end_time]
                start_i = 0
                end_i = 0
                for start_i in reversed(range(j)):
                    if times[start_i] < start_time:
                        start_i = min(j, start_i + 1)
                        break
                for end_i in range(j, len(times)):
                    if times[end_i] > end_time:
                        end_i = max(j, end_i - 1)
                        break
                wl_window = [wls[k] for k in range(start_i, end_i+1)]
                if wl_window:
                    wl_avg = sum(wl_window) / len(wl_window)
                    times_ma.append(times[j])
                    wls_ma.append(wl_avg)
            f61or63_times_ma.append(times_ma)
            f61or63_wls_ma.append(wls_ma)
        
    # Plot the observed data
    if has_obs_data:
        if options is not None and 'datum_adjustment_to_observation' in options:
            obs_wl = obs_wl + options['datum_adjustment_to_observation']
        valid = np.isfinite(obs_wl)
        obs_time_plot = obs_time[valid]
        obs_wl_plot = obs_wl[valid]
        if obs_time_plot.size > 0:
            ax.plot(
                obs_time_plot, obs_wl_plot * m2ft, '-',
                color=[0.5, 0.5, 0.5], label='Obs.',
            )
        
    # Plot the forecast data
    for i in range(len(f61or63files)):
        if f61or63colors is not None:
            if isinstance(f61or63colors, list):
                f61or63color = f61or63colors[i % len(f61or63colors)]
            else:
                f61or63color = f61or63colors
        else:
            f61or63color = _default_line_color(i)
        # Apply a 3-point moving maximum to f61or63_wls[i]
        f61or63_wls_ma3 = []
        wls = f61or63_wls[i]
        if movingaverage_window > 0:
            for j in range(len(wls)):
                # Get indices for the 3-point window centered at j
                start = max(0, j-1)
                end = min(len(wls)-1, j+1)
                window = wls[start:end+1]
                # Only consider non-nan values
                window_valid = [v for v in window if v is not None and not (isinstance(v, float) and np.isnan(v))]
                if window_valid:
                    ma_val = np.nanmean(window_valid)
                else:
                    ma_val = np.nan
                f61or63_wls_ma3.append(ma_val)
        else:
            f61or63_wls_ma3 = f61or63_wls[i]

        f61or63label = None
        if f61or63labels:
            if isinstance(f61or63labels, list):
                if len(f61or63labels) == len(f61or63files):
                    f61or63label = f61or63labels[i]
                else:
                    f61or63label = [f61or63labels[0] for _ in range(len(f61or63files))]
            else:
                f61or63label = f61or63labels

        ax.plot(f61or63_times[i], np.array(f61or63_wls_ma3) * m2ft, '-', color=f61or63color, label=f61or63label)
        # ax.plot(f61or63_times[i], f61or63_wls[i], fmt, label=f61or63labels[i])
    
    # Plot moving averages
    if plot_movingaverage:
        fmt = '--'
        
        if has_obs_data:
            ax.plot(obs_time_ma, np.array(obs_wl_ma) * m2ft, fmt, color='gray', label='Obs. (2d MA)')
            
        for i in range(len(f61or63files)):
            if f61or63colors is not None:
                if isinstance(f61or63colors, list):
                    f61or63color = f61or63colors[i % len(f61or63colors)]
                else:
                    f61or63color = f61or63colors
            else:
                f61or63color = _default_line_color(i)
            ax.plot(f61or63_times_ma[i], np.array(f61or63_wls_ma[i]) * m2ft, fmt, color=f61or63color, label=f'{f61or63labels[i]} (2d MA)')
            
    ax.xaxis.set_major_formatter(DateFormatter('%m-%d %H:%M'))
    fig.autofmt_xdate()
    ax.grid()
    if station_owner is not None:
        title_id = display_station_id if station_owner == 'CONTRAIL' else station_id
        title_name = station_name or station_id or title_id
        ax.set_title('{} {}: {}'.format(station_owner, title_id, title_name))
    else:
        ax.set_title('{}'.format(display_station_id))
    ax.set_ylabel('Water Level (ft)' if plot_in_foot else 'Water Level (m)')
    ax.set_xlim([date_start, date_end])
    # if min(obs_wl) > 0 and min(f61or63_wls[0]) > 0:
    #     ax.set_ylim([0, 1.05*max(max(obs_wl), max(f61or63_wls[0]))])
    
    # Handle legend placement - support "outside" locations
    if legend_loc == 'outside right':
        fig.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
        rect = legend_loc_rect if legend_loc_rect is not None else [0, 0, 0.9, 1]
        fig.tight_layout(rect=rect)
    elif legend_loc == 'outside left':
        fig.legend(bbox_to_anchor=(0.0, 1.0), loc='upper right')
        rect = legend_loc_rect if legend_loc_rect is not None else [0.1, 0, 1, 1]
        fig.tight_layout(rect=rect)
    elif legend_loc == 'outside top':
        fig.legend(bbox_to_anchor=(0.5, 1.0), loc='lower center')
        rect = legend_loc_rect if legend_loc_rect is not None else [0, 0, 1, 0.9]
        fig.tight_layout(rect=rect)
    elif legend_loc == 'outside bottom':
        fig.legend(bbox_to_anchor=(0.5, 0.0), loc='upper center')
        rect = legend_loc_rect if legend_loc_rect is not None else [0, 0.1, 1, 1]
        fig.tight_layout(rect=rect)
    else:
        # Use ax.legend() for 'best' location since fig.legend() doesn't support it
        if legend_loc == 'best' or legend_loc == 0:
            ax.legend(loc='best')
        else:
            fig.legend(loc=legend_loc)
        if legend_loc_rect is not None:
            fig.tight_layout(rect=legend_loc_rect)
        else:
            fig.tight_layout()

    return station_name

def get_parser():
    import argparse
    parser = argparse.ArgumentParser(add_help=False, description='Plot hydrograph at a station.')
    parser.add_argument(
        '--station-owner', type=str, required=True,
        choices=['NOAA', 'USGS', 'CONTRAIL', 'SECOORA', 'NONE'],
        help='Station owner: NOAA, USGS, CONTRAIL, SECOORA, or NONE. '
             'Observation data will not be plotted if station_owner is NONE.')
    parser.add_argument('--station-id', type=str, required=False, default=None, help='Station ID')
    parser.add_argument('--station-lon', type=float, required=False, default=None, help='Station longitude')
    parser.add_argument('--station-lat', type=float, required=False, default=None, help='Station latitude')
    parser.add_argument('--station-datum', type=str, required=False, default=None, help='Station datum: MSL or NAVD')
    parser.add_argument('--date-start', type=str, required=True, help='Start date (YYYY-MM-DD)')
    parser.add_argument('--date-end', type=str, required=True, help='End date (YYYY-MM-DD)')
    parser.add_argument('--f61or63files', type=str, nargs='+', required=True, help='List of fort.61.nc or fort.63.nc files')
    parser.add_argument('--f61or63starts', type=str, nargs='+', required=False, default=None, help='List of file start times')
    parser.add_argument('--f61or63labels', type=str, nargs='+', required=False, default=None, help='List of file labels')
    parser.add_argument('--f61or63colors', type=str, nargs='+', required=False, default=None, help='List of file colors')
    parser.add_argument('--f61or63concat', action='store_true', help='Concatenate files')
    parser.add_argument('--connect', action='store_true',
                        help='Bridge the line between consecutive concatenated '
                             'files by carrying over the last time/value to the '
                             'start of the next file (used with --f61or63concat)')
    parser.add_argument('--f63files-fallback', type=str, nargs='+', required=False, default=None, help='List of fort.63.nc fallback files for when station not found in fort.61.nc')
    parser.add_argument('--plot-movingaverage', action='store_true', help='Plot moving average')
    parser.add_argument('--adjust-datum-by-mean-error', type=float, default=0, metavar='period_in_days',
                        help='Adjust solution datum by subtracting the mean instantaneous error over the first period_in_days days (obs interpolated at solution times). Default 0 = no adjustment.')
    parser.add_argument('--outputfile', type=str, required=True, help='Output figure file name')
    parser.add_argument(
        '--fig-width',
        type=float,
        default=10.0,
        help='Figure width in inches (default: 10.0)',
    )
    parser.add_argument(
        '--fig-height',
        type=float,
        default=6.0,
        help='Figure height in inches (default: 6.0)',
    )
    parser.add_argument(
        '--cache-dir',
        help='Directory for caching downloaded observation data (passed to get_obswl)')

    contrail_group = parser.add_argument_group('CONTRAIL options')
    contrail_group.add_argument(
        '--username',
        help='Username for CONTRAIL authentication (required when station-owner is CONTRAIL)')
    contrail_group.add_argument(
        '--password',
        help='Password for CONTRAIL authentication (required when station-owner is CONTRAIL)')
    contrail_group.add_argument(
        '--sensor-type',
        choices=['water_elevation', 'stream_elevation', 'stage'],
        default='water_elevation',
        help='Sensor type for CONTRAIL (default: water_elevation)')
    contrail_group.add_argument(
        '--station-id-type',
        choices=['auto', 'contrail', 'f61'],
        default=None,
        help=(
            'CONTRAIL station id interpretation: contrail=integer site id (e.g. 1205), '
            'f61=fort.61 code (e.g. EGHN7), auto=detect. Omit for legacy behavior. '
            'Combined form 1205/EGHN7 always splits contrail/f61 ids.'
        ),
    )
    return parser

def main(args=None):
    import sys
    from datetime import datetime
    import matplotlib.pyplot as plt
    from glob import glob
    if args is None:
        args = get_parser().parse_args()
    date_start = datetime.strptime(args.date_start, '%Y-%m-%d')
    date_end = datetime.strptime(args.date_end, '%Y-%m-%d')

    station_owner = None if args.station_owner == 'NONE' else args.station_owner

    options = None
    if station_owner == 'CONTRAIL':
        if not args.username or not args.password:
            print("Error: CONTRAIL requires --username and --password", file=sys.stderr)
            sys.exit(1)
        options = {
            'username': args.username,
            'password': args.password,
            'sensor_type': args.sensor_type,
        }
    
    f61or63files = []
    for f61or63file in args.f61or63files:
        f61or63files.extend(glob(f61or63file))
    if args.f61or63concat:
        f61or63files = [f61or63files]
    
    if args.f61or63starts:
        f61or63starts = [datetime.strptime(f61or63start, '%Y-%m-%d') if f61or63start else None for f61or63start in args.f61or63starts]
    else:
        f61or63starts = None
    
    # Process fallback files if provided
    f63files_fallback = None
    if hasattr(args, 'f63files_fallback') and args.f63files_fallback:
        f63files_fallback = []
        for f63file in args.f63files_fallback:
            f63files_fallback.extend(glob(f63file))
        if args.f61or63concat:
            f63files_fallback = [f63files_fallback]
    
    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
    plot_hydrograph_at_station(
        fig, ax,
        station_owner, args.station_id, args.station_lon, args.station_lat, args.station_datum,
        date_start, date_end,
        f61or63files, f61or63starts, args.f61or63labels, args.f61or63colors,
        f63files_fallback=f63files_fallback,
        plot_movingaverage=args.plot_movingaverage,
        options=options,
        cache_dir=args.cache_dir,
        station_id_type=getattr(args, 'station_id_type', None),
        adjust_datum_by_mean_error_period_days=args.adjust_datum_by_mean_error,
        connect=getattr(args, 'connect', False))
    plt.savefig(args.outputfile)

if __name__ == '__main__':
    main()