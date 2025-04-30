def plot_hydrograph_at_station(
        fig, ax,
        station_owner, station_id, station_lon, station_lat, station_datum,
        date_start, date_end, 
        f63files, f63starts, f63labels, plot_movingaverage=False):
    import os
    import sys
    import requests
    from datetime import datetime, timedelta
    import matplotlib.pyplot as plt
    import pytz
    from matplotlib.dates import DateFormatter
    import pandas as pd
    import numpy as np
    from adcircutils.plot.get_obswl import get_obswl
    from adcircutils.plot.get_f63wl_at import get_f63wl_at
    
    # Get the observed water level data
    if station_owner is None:
        station_name, station_lon_, station_lat_, obs_time, obs_wl = \
            get_obswl(station_owner, station_id, date_start, date_end, station_datum)
        if station_lon is None:
            station_lon = station_lon_
            station_lat = station_lat_
    else:
        if station_id is None:
            raise ValueError('Station ID is required if station_owner is not NONE')
        if station_lon is None:
            raise ValueError('Station longitude is required if station_owner is not NONE')
            
    # Get the forecast water level data
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
        print('Processing f63 files', end='')
        for j in range(len(f63file)):
            f63filej = f63file[j]
            f63_startj = f63_start[j]
            # print('Processing {}...'.format(f63filej))
            print('.', end='')
            f63_timej, f63_wlj = get_f63wl_at(f63filej, station_lon, station_lat)
            if f63_startj:
                tdj = f63_startj - f63_timej[0]
                f63_timej = [tdj + t for t in f63_timej]
            f63_time.extend(f63_timej.tolist())
            f63_wl.extend(f63_wlj.tolist())
        f63_times.append(f63_time)
        f63_wl = [wl if abs(wl) <= 100 else np.nan for wl in f63_wl]
        f63_wls.append(f63_wl)
        print(' Done.')
        
    # Calculate moving average over a 2-day window
    if plot_movingaverage:
        window_size = timedelta(days=2)
        
        obs_time_ma = []
        obs_wl_ma = []
        for i in range(len(obs_time)):
            start_time = obs_time[i] - window_size / 2
            end_time = obs_time[i] + window_size / 2
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
                
        f63_times_ma = []
        f63_wls_ma = []
        for i in range(len(f63_times)):
            times = f63_times[i]
            wls = f63_wls[i]
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
            f63_times_ma.append(times_ma)
            f63_wls_ma.append(wls_ma)
        
    # Plot the observed data
    if station_owner is not None:
        ax.plot(obs_time, obs_wl, 'k-', label='Obs.')
        
    # Plot the forecast data
    for i in range(len(f63files)):
        if i == 0:
            fmt = 'b-'
        elif i == 1:
            fmt = 'r-'
        else:
            fmt = 'y-'
        ax.plot(f63_times[i], f63_wls[i], fmt, label=f63labels[i])
    
    # Plot moving averages
    if plot_movingaverage:
        fmt = '-'
        ax.plot(obs_time_ma, obs_wl_ma, fmt, color='gray', label='Obs. (2d MA)')
        for i in range(len(f63files)):
            if i == 0:
                color='lightblue'
            elif i == 1:
                color='lightcoral'
            else:
                color='lightgreen'
            ax.plot(f63_times_ma[i], f63_wls_ma[i], fmt, color=color, label=f'{f63labels[i]} (2d MA)')
            
    ax.xaxis.set_major_formatter(DateFormatter('%m-%d %H:%M'))
    fig.autofmt_xdate()
    ax.grid()
    if station_owner is not None:
        ax.set_title('{} {}: {}'.format(station_owner, station_id, station_name))
    else:
        ax.set_title('{}'.format(station_id))
    ax.set_ylabel('Water Level (m)')
    ax.set_xlim([date_start, date_end])
    # if min(obs_wl) > 0 and min(f63_wls[0]) > 0:
    #     ax.set_ylim([0, 1.05*max(max(obs_wl), max(f63_wls[0]))])
    ax.legend()

    return station_name

def main():
    from datetime import datetime
    import matplotlib.pyplot as plt
    import argparse
    parser = argparse.ArgumentParser(description='Plot hydrograph at a station.')
    parser.add_argument('--station_owner', type=str, required=True, help='Station owner: NOAA, USGS, SECOORA, or NONE. Observation data will not be plotted if station_owner is NONE.')
    parser.add_argument('--station_id', type=str, required=False, default=None, help='Station ID')
    parser.add_argument('--station_lon', type=float, required=False, default=None, help='Station longitude')
    parser.add_argument('--station_lat', type=float, required=False, default=None, help='Station latitude')
    parser.add_argument('--station_datum', type=str, required=False, default=None, help='Station datum: MSL or NAVD')
    parser.add_argument('--date_start', type=str, required=True, help='Start date (YYYY-MM-DD)')
    parser.add_argument('--date_end', type=str, required=True, help='End date (YYYY-MM-DD)')
    parser.add_argument('--f63files', type=str, nargs='+', required=True, help='List of f63 files')
    parser.add_argument('--f63starts', type=str, nargs='+', required=True, help='List of f63 start times')
    parser.add_argument('--f63labels', type=str, nargs='+', required=True, help='List of f63 labels')
    parser.add_argument('--plot_movingaverage', action='store_true', help='Plot moving average')
    parser.add_argument('--outputfile', type=str, required=True, help='Output figure file name')

    args = parser.parse_args()

    date_start = datetime.strptime(args.date_start, '%Y-%m-%d')
    date_end = datetime.strptime(args.date_end, '%Y-%m-%d')
    f63starts = [datetime.strptime(f63start, '%Y-%m-%d') if f63start else None for f63start in args.f63starts]

    fig, ax = plt.subplots()
    plot_hydrograph_at_station(
        fig, ax,
        args.station_owner, args.station_id, args.station_lon, args.station_lat, args.station_datum,
        date_start, date_end, 
        args.f63files, f63starts, args.f63labels, args.plot_movingaverage)
    plt.savefig(args.outputfile)

if __name__ == '__main__':
    main()