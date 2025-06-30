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
        
        date_start_str = date_start.strftime('%Y-%m-%dT%H:%M')
        date_end_str = date_end.strftime('%Y-%m-%dT%H:%M')

        # Create ERDDAP client
        e = ERDDAP(
            server='https://erddap.secoora.org/erddap',
            protocol='tabledap'
        )
        
        # Get the dataset metadata to find available variables
        try:
            # Get the full metadata for the dataset
            metadata_url = f'https://erddap.secoora.org/erddap/info/{station_id}/index.json'
            metadata_response = requests.get(metadata_url)
            
            if metadata_response.status_code == 200:
                metadata_json = metadata_response.json()
                
                # Find the station name from metadata
                station_name = station_id
                for attr in metadata_json['table']['rows']:
                    if attr[0] == 'attribute' and attr[1] == 'station' and attr[2] == 'long_name':
                        station_name = attr[4]
                        break
                
                # Find water level variable name
                water_level_var = None
                for attr in metadata_json['table']['rows']:
                    if attr[0] == 'variable':
                        var_name = attr[1]
                        if 'water_level' in var_name.lower() or 'sea_surface_height' in var_name.lower() or 'water_surface_height' in var_name.lower():
                            water_level_var = var_name
                            break
                
                if not water_level_var:
                    raise ValueError(f"Couldn't find water level variable for station {station_id}")
                
                print(f"Found water level variable: {water_level_var}")
            else:
                raise ValueError(f"Failed to retrieve metadata for station {station_id}")
        except Exception as e:
            print(f"Error retrieving metadata: {str(e)}")
            raise
        
        # Now get the actual water level data with the found variable
        e.response = 'csv'
        e.dataset_id = station_id
        
        # First try to get the station's fixed position
        e.variables = ['time', 'station']
        e.constraints = {'time>=': date_start_str, 'time<=': date_start_str}
        
        try:
            # Get station info first to get lat/lon
            info_url = f'https://erddap.secoora.org/erddap/info/{station_id}/index.json'
            info_response = requests.get(info_url)
            
            if info_response.status_code == 200:
                info_json = info_response.json()
                
                # Find the fixed station latitude and longitude from global attributes
                station_lon = None
                station_lat = None
                
                for attr in info_json['table']['rows']:
                    if attr[0] == 'attribute' and attr[1] == 'NC_GLOBAL' and attr[2] == 'geospatial_lon_min':
                        station_lon = float(attr[4])
                    if attr[0] == 'attribute' and attr[1] == 'NC_GLOBAL' and attr[2] == 'geospatial_lat_min':
                        station_lat = float(attr[4])
                
                if station_lon is None or station_lat is None:
                    # Try looking for longitude and latitude variables instead
                    for attr in info_json['table']['rows']:
                        if attr[0] == 'variable' and 'longitude' in attr[1].lower():
                            lon_var = attr[1]
                            lat_var = None
                            # Look for corresponding latitude variable
                            for attr2 in info_json['table']['rows']:
                                if attr2[0] == 'variable' and 'latitude' in attr2[1].lower():
                                    lat_var = attr2[1]
                                    break
                            
                            if lat_var:
                                # Now get both variables
                                e.variables = ['time', lon_var, lat_var, water_level_var]
                                break
            
            # If we still don't have lat/lon variables, use defaults
            if not ('lon_var' in locals() and 'lat_var' in locals()):
                e.variables = ['time', water_level_var]
                
        except Exception as e:
            print(f"Error getting station position: {str(e)}")
            # Default to basic variables
            e.variables = ['time', water_level_var]
        
        # Get the actual water level data
        e.constraints = {
            'time>=': date_start_str,
            'time<=': date_end_str
        }
        
        try:
            obs_data = e.to_pandas(parse_dates=True)
            print(f"Data columns: {obs_data.columns.tolist()}")
            
            # Handle column names with units in parentheses
            time_col = next((col for col in obs_data.columns if col.startswith('time')), None)
            if not time_col:
                raise ValueError("Cannot find time column in data")
            
            water_level_col = next((col for col in obs_data.columns if water_level_var in col), None)
            if not water_level_col:
                raise ValueError(f"Cannot find {water_level_var} column in data")
            
            obs_time = obs_data[time_col]
            # Convert time values to datetime objects if they're not already
            if obs_time.dtype == 'object' or isinstance(obs_time.iloc[0], str):
                obs_time = pd.to_datetime(obs_time)
            
            obs_wl = obs_data[water_level_col]
            
            # If we have lat/lon in the data, use it
            if 'lon_var' in locals() and any(lon_var in col for col in obs_data.columns):
                lon_col = next(col for col in obs_data.columns if lon_var in col)
                lat_col = next(col for col in obs_data.columns if lat_var in col)
                station_lon = obs_data[lon_col][0]
                station_lat = obs_data[lat_col][0]
            
            # If we still don't have lat/lon, use hardcoded values from the script
            if station_lon is None or station_lat is None:
                # Use the values specified in the plot_hydrographs_at_stations.py
                # Those should be passed to this function
                pass
                
            print(f"Station position: {station_lon}, {station_lat}")
            
        except Exception as e:
            print(f"Error processing SECOORA data: {str(e)}")
            raise
    else:
        raise ValueError('Invalid station owner')
    
    return station_name, station_lon, station_lat, obs_time, obs_wl