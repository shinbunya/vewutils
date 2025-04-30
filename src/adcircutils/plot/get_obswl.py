# %%
import requests
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import pytz
import dataretrieval.nwis as nwis
from erddapy import ERDDAP

def get_obswl(station_owner, station_id, date_start, date_end, datum):
    ft2m = 0.3048
    tzutc = pytz.timezone('UTC')

    if station_owner == 'NOAA':
        obs_time = []
        obs_wl = []
        date_start_i = date_start
        while date_start_i <= date_end:
            date_start_str = date_start_i.strftime('%Y%m%d')
            date_end_i = date_start_i + timedelta(days=30)
            if date_end_i > date_end:
                date_end_i = date_end
            date_end_str = date_end_i.strftime('%Y%m%d')
            date_start_str = date_start_i.strftime('%Y%m%d')
            date_end_str = date_end_i.strftime('%Y%m%d')
            obs_url = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=water_level&application=NOS.COOPS.TAC.WL&begin_date={:s}&end_date={:s}&datum={:s}&station={:s}&time_zone=GMT&units=metric&format=json'\
                .format(date_start_str, date_end_str, datum, station_id)
            print(obs_url)
            response = requests.get(obs_url)
            if response.status_code == 200:
                obs_data = response.json()
            else:
                print(f"Failed to retrieve data: {response.status_code}")
                return None, None, None, None, None
            obs_time.extend([datetime.strptime(obs_data['data'][i]['t'], '%Y-%m-%d %H:%M') for i in range(len(obs_data['data']))])
            obs_wl.extend([float(obs_data['data'][i]['v']) if obs_data['data'][i]['v'] else np.nan for i in range(len(obs_data['data']))])
            date_start_i += timedelta(days=31)
            
        station_name = obs_data['metadata']['name']
        station_lon = float(obs_data['metadata']['lon'])
        station_lat = float(obs_data['metadata']['lat'])
    elif station_owner == 'USGS':
        date_start_str = date_start.strftime('%Y-%m-%dT%H:%M')
        date_end_str = date_end.strftime('%Y-%m-%dT%H:%M')
        dfst = nwis.get_record(sites=station_id, service='site')
        dfiv = nwis.get_record(sites=station_id, service='iv', start=date_start_str, end=date_end_str)
        station_name = dfst['station_nm'][0]
        station_lon = dfst['dec_long_va'][0]
        station_lat = dfst['dec_lat_va'][0]
        obs_time = pd.to_datetime(dfiv.index)
        print('alt_va = ', dfst['alt_va'][0])
        if '00065' in dfiv.columns:
            obs_wl = (dfiv['00065'] + dfst['alt_va'][0])*ft2m
        elif '62620' in dfiv.columns:
            obs_wl = dfiv['62620']*ft2m
        else:
            raise KeyError('No valid column found in dfiv')
        # station_name = obs_data['value']['timeSeries'][0]['sourceInfo']['siteName']
        # station_lon = float(obs_data['value']['timeSeries'][0]['sourceInfo']['geoLocation']['geogLocation']['longitude'])
        # station_lat = float(obs_data['value']['timeSeries'][0]['sourceInfo']['geoLocation']['geogLocation']['latitude'])
        # obs_time = [datetime.strptime(obs_data['value']['timeSeries'][0]['values'][0]['value'][i]['dateTime'], '%Y-%m-%dT%H:%M:%S.%f%z').astimezone(tzutc) for i in range(len(obs_data['value']['timeSeries'][0]['values'][0]['value']))]
        # obs_wl = [float(obs_data['value']['timeSeries'][0]['values'][0]['value'][i]['value'])*ft2m for i in range(len(obs_data['value']['timeSeries'][0]['values'][0]['value']))]
    elif station_owner == 'SECOORA':
        if datum != 'NAVD':
            raise ValueError('SECOORA only supports NAVD datum')
        
        e = ERDDAP(
            server='https://erddap.secoora.org/erddap',
            protocol='tabledap'
        )
        
        e.response = 'csv'
        e.dataset_id = station_id
        e.variables = [
            'time',
            'latitude',
            'longitude',
            'short_name',
            'sea_surface_height_above_sea_level_geoid_navd88_surveyed_navd88'
        ]
        e.constraints = {
            'time>=': date_start_str,
            'time<=': date_end_str
        }
        obs_data = e.to_pandas(parse_dates=True)
        station_name = obs_data['short_name'][0]
        station_lon = obs_data['longitude'][0]
        station_lat = obs_data['latitude'][0]
        obs_time = obs_data['time']
        obs_wl = obs_data['sea_surface_height_above_sea_level_geoid_navd88_surveyed_navd88']
        
    else:
        raise ValueError('Invalid station owner')
    
    return station_name, station_lon, station_lat, obs_time, obs_wl