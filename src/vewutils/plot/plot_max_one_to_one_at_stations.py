def plot_max_one_to_one_at_stations(
        fig, ax,
        station_owners, station_ids, station_lons, station_lats, station_datums,
        date_start, date_end, 
        f63files, f63starts, f63labels, movingaverage_window=0, f63colors=None, plot_in_foot=False, options_list=None,
        xmin=None, xmax=None, ymin=None, ymax=None, show_station_names=False, show_stats=True):
    import os
    import sys
    import requests
    from datetime import datetime, timedelta
    import matplotlib.pyplot as plt
    import pytz
    from matplotlib.dates import DateFormatter
    import pandas as pd
    import numpy as np
    from vewutils.plot.get_obswl import get_obswl
    from vewutils.plot.get_f63wl_at import get_f63wl_at
    
    # Set unit conversion factor
    if plot_in_foot:
        m2ft = 1.0/0.3048
    else:
        m2ft = 1.0
    
    # Ensure all station parameters are lists
    if not isinstance(station_owners, list):
        station_owners = [station_owners]
    if not isinstance(station_ids, list):
        station_ids = [station_ids]
    if not isinstance(station_lons, list):
        station_lons = [station_lons]
    if not isinstance(station_lats, list):
        station_lats = [station_lats]
    if not isinstance(station_datums, list):
        station_datums = [station_datums]
    
    # Ensure all lists have the same length
    n_stations = len(station_owners)
    if len(station_ids) != n_stations:
        raise ValueError('All station parameter lists must have the same length')
    if len(station_lons) != n_stations:
        raise ValueError('All station parameter lists must have the same length')
    if len(station_lats) != n_stations:
        raise ValueError('All station parameter lists must have the same length')
    if len(station_datums) != n_stations:
        raise ValueError('All station parameter lists must have the same length')
    
    # Collect all observed and predicted maximum values organized by forecast
    all_obs_max = [[] for _ in range(len(f63files))]
    all_pred_max = [[] for _ in range(len(f63files))]
    all_station_names = [[] for _ in range(len(f63files))]
    station_names = []
    
    # Process each station
    for station_idx in range(n_stations):
        station_owner = station_owners[station_idx]
        station_id = station_ids[station_idx]
        station_lon = station_lons[station_idx]
        station_lat = station_lats[station_idx]
        station_datum = station_datums[station_idx]
        options = options_list[station_idx]
        
        print(f'Processing station {station_idx + 1}/{n_stations}: {station_id}')
        
        # Get the observed water level data
        if station_owner is not None and station_owner != 'NONE':
            station_name, station_lon_, station_lat_, obs_time, obs_wl = \
                get_obswl(station_owner, station_id, date_start, date_end, station_datum, options)
            if station_lon is None:
                station_lon = station_lon_
                station_lat = station_lat_
        else:
            if station_id is None:
                raise ValueError('Station ID is required if station_owner is NONE')
            if station_lon is None:
                raise ValueError('Station longitude is required if station_owner is NONE')
            station_name = station_id
            obs_time = []
            obs_wl = []
        
        # Convert obs_time and obs_wl to numpy arrays
        obs_time = np.array(obs_time)
        obs_wl = np.array(obs_wl)
        
        # Get the forecast water level data for this station
        f63_times = []
        f63_wls = []
        for i in range(len(f63files)):
            if isinstance(f63files[i], list):
                f63file = f63files[i]
                f63_start = f63starts[i]
            else:
                f63file = [f63files[i]]
                f63_start = [f63starts[i]]
            
            f63_time = []
            f63_wl = []
            for j in range(len(f63file)):
                f63filej = f63file[j]
                f63_startj = f63_start[j]
                f63_timej, f63_wlj = get_f63wl_at(f63filej, station_lon, station_lat)
                if f63_startj:
                    tdj = f63_startj - f63_timej[0]
                    f63_timej = [tdj + t for t in f63_timej]
                if f63_timej is not None:
                    f63_time.extend(f63_timej.tolist())
                    f63_wl.extend(f63_wlj.tolist())
            f63_times.append(f63_time)
            f63_wl = [wl if abs(wl) <= 100 else np.nan for wl in f63_wl]
            f63_wls.append(f63_wl)
        
        # Apply moving average if specified
        if movingaverage_window > 0:
            # Apply moving average to observed data
            if len(obs_wl) > 0:
                obs_wl_ma = []
                for i in range(len(obs_wl)):
                    start = max(0, i - movingaverage_window // 2)
                    end = min(len(obs_wl) - 1, i + movingaverage_window // 2)
                    window = obs_wl[start:end+1]
                    window_valid = [v for v in window if not np.isnan(v)]
                    if window_valid:
                        ma_val = np.nanmean(window_valid)
                    else:
                        ma_val = np.nan
                    obs_wl_ma.append(ma_val)
                obs_wl_processed = np.array(obs_wl_ma)
            else:
                obs_wl_processed = obs_wl
            
            # Apply moving average to forecast data
            f63_wls_processed = []
            for i in range(len(f63_wls)):
                wls = f63_wls[i]
                f63_wls_ma = []
                for j in range(len(wls)):
                    start = max(0, j - movingaverage_window // 2)
                    end = min(len(wls) - 1, j + movingaverage_window // 2)
                    window = wls[start:end+1]
                    window_valid = [v for v in window if v is not None and not (isinstance(v, float) and np.isnan(v))]
                    if window_valid:
                        ma_val = np.nanmean(window_valid)
                    else:
                        ma_val = np.nan
                    f63_wls_ma.append(ma_val)
                f63_wls_processed.append(f63_wls_ma)
        else:
            obs_wl_processed = obs_wl
            f63_wls_processed = f63_wls
        
        # Apply datum adjustment to observations if specified
        if station_owner is not None and station_owner != 'NONE':
            if options is not None and 'datum_adjustment_to_observation' in options:
                obs_wl_processed = obs_wl_processed + options['datum_adjustment_to_observation']
        
        # Calculate maximum values for this station
        if len(obs_wl_processed) > 0:
            obs_max = np.nanmax(obs_wl_processed) * m2ft
        else:
            obs_max = np.nan
        
        # Store station name
        station_names.append(station_name)
        
        # For each forecast, calculate and store maximum values
        for i in range(len(f63files)):
            if len(f63_wls_processed[i]) > 0:
                f63_max = np.nanmax(f63_wls_processed[i]) * m2ft
            else:
                f63_max = np.nan
            
            # Only store valid values (not NaN)
            if not np.isnan(obs_max) and not np.isnan(f63_max):
                all_obs_max[i].append(obs_max)
                all_pred_max[i].append(f63_max)
                all_station_names[i].append(station_name)
    
    # Calculate bias and MAE for each forecast
    forecast_stats = []
    for i in range(len(f63files)):
        forecast_obs = all_obs_max[i]
        forecast_pred = all_pred_max[i]
        
        if forecast_obs and forecast_pred:
            # Calculate bias (mean error)
            errors = np.array(forecast_pred) - np.array(forecast_obs)
            bias = np.nanmean(errors)
            mae = np.nanmean(np.abs(errors))
            forecast_stats.append((bias, mae))
        else:
            forecast_stats.append((np.nan, np.nan))
    
    # Plot one-to-one for each forecast
    for i in range(len(f63files)):
        if f63colors is not None:
            f63color = f63colors[i]
        else:
            if i == 0:
                f63color = 'b'
            elif i == 1:
                f63color = 'r'
            elif i == 2:
                f63color = 'y'
            else:
                f63color = 'g'
        
        # Get values for this forecast
        forecast_obs = all_obs_max[i]
        forecast_pred = all_pred_max[i]
        
        # Plot the points for this forecast
        if forecast_obs and forecast_pred:
            ax.scatter(forecast_obs, forecast_pred, color=f63color, s=50, label=f63labels[i], alpha=0.7)
            
            # Add station names at marker locations if requested
            if show_station_names:
                forecast_station_names = all_station_names[i]
                for j, (x, y, name) in enumerate(zip(forecast_obs, forecast_pred, forecast_station_names)):
                    ax.annotate(name, (x, y), xytext=(5, 5), textcoords='offset points',
                               fontsize=8, alpha=0.8, color=f63color)
    
    # Plot 1:1 line
    all_obs_values = []
    all_pred_values = []
    for i in range(len(f63files)):
        all_obs_values.extend(all_obs_max[i])
        all_pred_values.extend(all_pred_max[i])
    
    if all_obs_values and all_pred_values:
        if xmin is not None and xmax is not None:
            min_val = xmin
            max_val = xmax
        else:
            min_val = min(min(all_obs_values), min(all_pred_values))
            max_val = max(max(all_obs_values), max(all_pred_values))
        if ymin is not None and ymax is not None:
            min_val = ymin
            max_val = ymax
        else:
            min_val = min(min(all_obs_values), min(all_pred_values))
            max_val = max(max(all_obs_values), max(all_pred_values))
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='1:1 Line')
    
    # Set labels and title
    ax.set_xlabel('Observed Maximum Water Level (ft)' if plot_in_foot else 'Observed Maximum Water Level (m)')
    ax.set_ylabel('Predicted Maximum Water Level (ft)' if plot_in_foot else 'Predicted Maximum Water Level (m)')
    
    if n_stations == 1:
        if station_owners[0] is not None and station_owners[0] != 'NONE':
            ax.set_title('{} {}: {} - Maximum Values'.format(station_owners[0], station_ids[0], station_names[0]))
        else:
            ax.set_title('{} - Maximum Values'.format(station_ids[0]))
    else:
        ax.set_title('Multiple Stations - Maximum Values ({} stations)'.format(n_stations))
    
    # Add grid and legend
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Draw text on ax of bias and MAE for each forecast (if requested)
    if show_stats:
        colors = f63colors if f63colors is not None else ['b', 'r', 'y', 'g']
        for i, (bias, mae) in enumerate(forecast_stats):
            if not np.isnan(bias) and not np.isnan(mae):
                textstr = f'Bias: {bias:.3f}\nMAE: {mae:.3f}'
                props = dict(boxstyle='round', facecolor=colors[i], alpha=0.5)
                ax.text(0.05, 0.95 - i*0.1, textstr, transform=ax.transAxes, fontsize=10,
                        verticalalignment='top', bbox=props)
    
    # Set equal aspect ratio for proper 1:1 visualization
    ax.set_aspect('equal', adjustable='box')
    
    # Set axis limits if specified
    if xmin is not None or xmax is not None:
        ax.set_xlim(left=xmin, right=xmax)
    if ymin is not None or ymax is not None:
        ax.set_ylim(bottom=ymin, top=ymax)
    
    return station_names

def get_parser():
    import argparse
    parser = argparse.ArgumentParser(add_help=False, description='Plot one-to-one maximum values at multiple stations.')
    parser.add_argument('--station-owners', type=str, nargs='+', required=True, help='List of station owners: NOAA, USGS, SECOORA, or NONE. Observation data will not be plotted if station_owner is NONE.')
    parser.add_argument('--station-ids', type=str, nargs='+', required=False, default=None, help='List of station IDs')
    parser.add_argument('--station-lons', type=float, nargs='+', required=False, default=None, help='List of station longitudes')
    parser.add_argument('--station-lats', type=float, nargs='+', required=False, default=None, help='List of station latitudes')
    parser.add_argument('--station-datums', type=str, nargs='+', required=False, default=None, help='List of station datums: MSL or NAVD')
    parser.add_argument('--date-start', type=str, required=True, help='Start date (YYYY-MM-DD)')
    parser.add_argument('--date-end', type=str, required=True, help='End date (YYYY-MM-DD)')
    parser.add_argument('--f63files', type=str, nargs='+', required=True, help='List of f63 files')
    parser.add_argument('--f63starts', type=str, nargs='+', required=True, help='List of f63 start times')
    parser.add_argument('--f63labels', type=str, nargs='+', required=True, help='List of f63 labels')
    parser.add_argument('--movingaverage-window', type=int, default=0, help='Moving average window size (0 to disable)')
    parser.add_argument('--f63colors', type=str, nargs='*', default=None, help='List of f63 colors')
    parser.add_argument('--plot-in-foot', action='store_true', help='Plot in feet instead of meters')
    parser.add_argument('--xmin', type=float, default=None, help='Minimum x-axis value (observed water level)')
    parser.add_argument('--xmax', type=float, default=None, help='Maximum x-axis value (observed water level)')
    parser.add_argument('--ymin', type=float, default=None, help='Minimum y-axis value (predicted water level)')
    parser.add_argument('--ymax', type=float, default=None, help='Maximum y-axis value (predicted water level)')
    parser.add_argument('--show-station-names', action='store_true', help='Show station names at marker locations')
    parser.add_argument('--no-stats', action='store_true', help='Hide bias and MAE statistics display')
    parser.add_argument('--outputfile', type=str, required=True, help='Output figure file name')
    return parser

def main(args=None):
    from datetime import datetime
    import matplotlib.pyplot as plt
    if args is None:
        args = get_parser().parse_args()
    date_start = datetime.strptime(args.date_start, '%Y-%m-%d')
    date_end = datetime.strptime(args.date_end, '%Y-%m-%d')
    f63starts = [datetime.strptime(f63start, '%Y-%m-%d') if f63start else None for f63start in args.f63starts]
    fig, ax = plt.subplots()
    plot_max_one_to_one_at_station(
        fig, ax,
        args.station_owners, args.station_ids, args.station_lons, args.station_lats, args.station_datums,
        date_start, date_end, 
        args.f63files, f63starts, args.f63labels, args.movingaverage_window, args.f63colors, args.plot_in_foot, None,
        args.xmin, args.xmax, args.ymin, args.ymax, args.show_station_names, not args.no_stats)
    plt.savefig(args.outputfile)

if __name__ == '__main__':
    main()
