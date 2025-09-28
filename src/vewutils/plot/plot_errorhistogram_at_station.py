import os
import sys
import requests
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import pytz
from matplotlib.dates import DateFormatter
from vewutils.plot.get_obswl import get_obswl
from vewutils.plot.get_f63wl_at import get_f63wl_at
import pandas as pd
import argparse

def plot_errorhistogram_at_station(
        fig, ax,
        station_owner, station_id, station_lon, station_lat, station_datum,
        date_start, date_end, 
        f63files, f63starts, f63labels, plot_movingaverage=False, f63colors=None, options=None):
    # Get the observed water level data
    station_name, station_lon_, station_lat_, obs_time, obs_wl = \
        get_obswl(station_owner, station_id, date_start, date_end, station_datum, options)
    if station_lon is None:
        station_lon = station_lon_
        station_lat = station_lat_
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
            
            # Filter f63_timej and f63_wlj to only include times between date_start and date_end
            f63_timej = np.array(f63_timej)
            f63_wlj = np.array(f63_wlj)
            mask = (f63_timej >= date_start) & (f63_timej <= date_end)
            f63_timej = f63_timej[mask]
            f63_wlj = f63_wlj[mask]
            
            if f63_startj:
                tdj = f63_startj - f63_timej[0]
                f63_timej = [tdj + t for t in f63_timej]
            f63_time.extend(f63_timej.tolist())
            f63_wl.extend(f63_wlj.tolist())
        f63_times.append(f63_time)
        f63_wl = [wl if abs(wl) <= 100 else np.nan for wl in f63_wl]
        f63_wls.append(f63_wl)
        print('done.')
            
    # Interpolate obs_wl at each member of f63_times
    obs_wls_interp = []
    for times in f63_times:
        obs_df = pd.DataFrame({'time': obs_time, 'wl': obs_wl})
        obs_df.set_index('time', inplace=True)
        obs_df.index = obs_df.index.tz_localize(None)  # Remove timezone information
        times = pd.to_datetime(times).tz_localize(None)  # Remove timezone information
        obs_df = obs_df.reindex(times, method='nearest')
        obs_df = obs_df.interpolate(method='linear')
        wl_interp = obs_df['wl'].values
        obs_wls_interp.append(wl_interp)
        
    # Calculate errors f63_wls - obs_wls_interp for each member
    f63_errs = []
    for i in range(len(f63_wls)):
        error = np.array(f63_wls[i]) - np.array(obs_wls_interp[i])
        f63_errs.append(error)

    # Calculate the mean and the standard deviation of each member of f63_errs
    f63_errs_mean_std = []
    for error in f63_errs:
        mean = np.nanmean(error)
        std_dev = np.nanstd(error)
        f63_errs_mean_std.append((mean, std_dev))

    # Plot the error histogram
    if f63colors is not None:
        colors = f63colors
    else:
        colors = ['b', 'r', 'y', 'g']
    ax.hist(f63_errs, bins=20, label=f63labels, color=[colors[i] for i in range(len(f63_errs))])
    ax.set_title('{} {}: {}'.format(station_owner, station_id, station_name))
    ax.set_xlabel('Error (m)')
    ax.set_ylabel('Frequency')
    ax.legend()
    
    # Draw text on ax of mean and std_dev for each member of f63_errs
    for i, (mean, std_dev) in enumerate(f63_errs_mean_std):
        textstr = f'Mean: {mean:.2f}\nStd Dev: {std_dev:.2f}'
        props = dict(boxstyle='round', facecolor=colors[i], alpha=0.5)
        ax.text(0.05, 0.95 - i*0.1, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
                
    return station_name

def get_parser():
    parser = argparse.ArgumentParser(add_help=False, description='Plot error histogram at a station.')
    parser.add_argument('--station-owner', type=str, required=True, help='Station owner: NOAA, USGS, SECOORA, or NONE. Observation data will not be plotted if station-owner is NONE.')
    parser.add_argument('--station-id', type=str, required=False, default=None, help='Station ID')
    parser.add_argument('--station-lon', type=float, required=False, default=None, help='Station longitude')
    parser.add_argument('--station-lat', type=float, required=False, default=None, help='Station latitude')
    parser.add_argument('--station-datum', type=str, required=False, default=None, help='Station datum: MSL or NAVD')
    parser.add_argument('--date-start', type=str, required=True, help='Start date (YYYY-MM-DD)')
    parser.add_argument('--date-end', type=str, required=True, help='End date (YYYY-MM-DD)')
    parser.add_argument('--f63files', type=str, nargs='+', required=True, help='List of f63 files')
    parser.add_argument('--f63starts', type=str, nargs='+', required=True, help='List of f63 start times (YYYY-MM-DD)')
    parser.add_argument('--f63labels', type=str, nargs='+', required=True, help='List of f63 labels')
    parser.add_argument('--plot-movingaverage', action='store_true', help='Plot moving average')
    parser.add_argument('--f63colors', type=str, nargs='+', required=False, default=None, help='List of f63 colors')
    parser.add_argument('--outputfile', type=str, required=True, help='Output figure file name')
    return parser

def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    from datetime import datetime
    import matplotlib.pyplot as plt
    date_start = datetime.strptime(args.date_start, '%Y-%m-%d')
    date_end = datetime.strptime(args.date_end, '%Y-%m-%d')
    f63starts = [datetime.strptime(f63start, '%Y-%m-%d') if f63start else None for f63start in args.f63starts]
    fig, ax = plt.subplots()
    plot_errorhistogram_at_station(
        fig, ax,
        args.station_owner, args.station_id, args.station_lon, args.station_lat, args.station_datum,
        date_start, date_end,
        args.f63files, f63starts, args.f63labels, args.plot_movingaverage
    )
    plt.savefig(args.outputfile)

if __name__ == '__main__':
    main()