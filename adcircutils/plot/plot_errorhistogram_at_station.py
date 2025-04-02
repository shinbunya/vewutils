import os
import sys
import requests
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import pytz
from matplotlib.dates import DateFormatter
sys.path.append('/projects/storm_surge/COMT/sbunya/projects/SGS/work/EC120_Ian/0_src')
from get_obswl import get_obswl
from get_f63wl_at import get_f63wl_at
import pandas as pd

def plot_errorhistogram_at_station(
        fig, ax,
        station_owner, station_id, station_lon, station_lat, station_datum,
        date_start, date_end, 
        f63files, f63starts, f63labels, plot_movingaverage=False):
    # Get the observed water level data
    station_name, station_lon_, station_lat_, obs_time, obs_wl = \
        get_obswl(station_owner, station_id, date_start, date_end, station_datum)
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
    colors = ['b', 'r', 'y']
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