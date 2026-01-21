# %%
import requests
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import pytz
import dataretrieval.nwis as nwis
from erddapy import ERDDAP
from bs4 import BeautifulSoup
import logging
from urllib.parse import urlencode
import json
import os
import re
from pathlib import Path

def _sanitize_station_id(station_id):
    """
    Sanitize station ID for use in filenames.
    
    Parameters:
    -----------
    station_id : str
        Station identifier
    
    Returns:
    --------
    str
        Sanitized station ID (lowercase, special chars replaced with underscores)
    """
    # Replace non-alphanumeric characters (except hyphens and underscores) with underscores
    sanitized = re.sub(r'[^a-zA-Z0-9_-]', '_', str(station_id))
    # Convert to lowercase
    sanitized = sanitized.lower()
    # Collapse multiple consecutive underscores into one
    sanitized = re.sub(r'_+', '_', sanitized)
    # Remove leading/trailing underscores
    sanitized = sanitized.strip('_')
    return sanitized

def _get_cache_filename(station_owner, station_id, date_start, date_end, datum):
    """
    Generate cache filename based on naming convention.
    
    Format: {owner}_{sanitized_station}_{datum}_{start_date}_{end_date}.json
    
    Parameters:
    -----------
    station_owner : str
        Source of the data ('NOAA', 'USGS', 'CONTRAIL', 'SECOORA')
    station_id : str
        Station identifier
    date_start : datetime
        Start date for data retrieval
    date_end : datetime
        End date for data retrieval
    datum : str
        Datum for water level measurements
    
    Returns:
    --------
    str
        Cache filename
    """
    sanitized_station = _sanitize_station_id(station_id)
    sanitized_datum = _sanitize_station_id(datum)  # Reuse sanitization function for datum
    
    # Format dates as YYYYMMDDTHHMMSS (ISO format without colons/hyphens)
    date_start_naive = date_start.replace(tzinfo=None) if date_start.tzinfo else date_start
    date_end_naive = date_end.replace(tzinfo=None) if date_end.tzinfo else date_end
    start_str = date_start_naive.strftime('%Y%m%dT%H%M%S')
    end_str = date_end_naive.strftime('%Y%m%dT%H%M%S')
    
    filename = f"{station_owner.upper()}_{sanitized_station}_{sanitized_datum}_{start_str}_{end_str}.json"
    return filename

def _load_cache(cache_path):
    """
    Load data from cache file.
    
    Parameters:
    -----------
    cache_path : str or Path
        Path to cache file
    
    Returns:
    --------
    tuple or None
        (station_name, station_lon, station_lat, obs_time, obs_wl) if cache exists, None otherwise
    """
    try:
        with open(cache_path, 'r') as f:
            cache_data = json.load(f)
        
        # Reconstruct pandas Series from JSON
        obs_time = pd.Series(pd.to_datetime(cache_data['obs_time']))
        # Ensure timezone-aware (UTC) like the original data
        if obs_time.dt.tz is None:
            obs_time = obs_time.dt.tz_localize('UTC')
        else:
            obs_time = obs_time.dt.tz_convert('UTC')
        
        obs_wl = pd.Series(cache_data['obs_wl'])
        
        return (
            cache_data['station_name'],
            cache_data['station_lon'],
            cache_data['station_lat'],
            obs_time,
            obs_wl
        )
    except (FileNotFoundError, json.JSONDecodeError, KeyError) as e:
        return None

def _save_cache(cache_path, station_name, station_lon, station_lat, obs_time, obs_wl, datum_offset=None):
    """
    Save data to cache file.
    
    Parameters:
    -----------
    cache_path : str or Path
        Path to cache file
    station_name : str
        Station name
    station_lon : float
        Station longitude
    station_lat : float
        Station latitude
    obs_time : pd.Series
        Observation times
    obs_wl : pd.Series
        Observation water levels
    """
    # Ensure cache directory exists
    cache_dir = os.path.dirname(cache_path)
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
    
    # Convert pandas Series to JSON-serializable format
    cache_data = {
        'station_name': station_name,
        'station_lon': float(station_lon) if station_lon is not None else None,
        'station_lat': float(station_lat) if station_lat is not None else None,
        'datum_offset': float(datum_offset) if datum_offset is not None else None,
        'obs_time': [t.isoformat() if hasattr(t, 'isoformat') else str(t) for t in obs_time.tolist()],
        'obs_wl': [float(wl) if not pd.isna(wl) else None for wl in obs_wl.tolist()]
    }
    
    with open(cache_path, 'w') as f:
        json.dump(cache_data, f, indent=2)

def _get_noaa_data(station_id, date_start, date_end, datum, **kwargs):
    """Local method to retrieve NOAA water level data"""
    import pytz
    ft2m = 0.3048
    tzutc = pytz.timezone('UTC')
    
    # Check cache if cache_dir is provided
    cache_dir = kwargs.get('cache_dir')
    if cache_dir:
        cache_filename = _get_cache_filename('NOAA', station_id, date_start, date_end, datum)
        cache_path = os.path.join(cache_dir, cache_filename)
        cached_data = _load_cache(cache_path)
        if cached_data is not None:
            print(f"Loading NOAA data from cache: {cache_path}")
            return cached_data
    
    obs_time = []
    obs_wl = []
    date_start_i = date_start
    while date_start_i <= date_end:
        # Handle timezone-aware datetime objects
        date_start_naive = date_start_i.replace(tzinfo=None) if date_start_i.tzinfo else date_start_i
        date_start_str = date_start_naive.strftime('%Y%m%d')
        date_end_i = date_start_i + timedelta(days=30)
        if date_end_i > date_end:
            date_end_i = date_end
        date_end_naive = date_end_i.replace(tzinfo=None) if date_end_i.tzinfo else date_end_i
        date_end_str = date_end_naive.strftime('%Y%m%d')
        obs_url = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=water_level&application=NOS.COOPS.TAC.WL&begin_date={:s}&end_date={:s}&datum={:s}&station={:s}&time_zone=GMT&units=metric&format=json'\
            .format(date_start_str, date_end_str, datum, station_id)
        print(obs_url)
        response = requests.get(obs_url)
        if response.status_code == 200:
            obs_data = response.json()
        else:
            print(f"Failed to retrieve data: {response.status_code}")
            return None, None, None, None, None
        # Parse times and make them UTC-aware
        times_parsed = [datetime.strptime(obs_data['data'][i]['t'], '%Y-%m-%d %H:%M') for i in range(len(obs_data['data']))]
        times_utc = [tzutc.localize(t) for t in times_parsed]
        obs_time.extend(times_utc)
        obs_wl.extend([float(obs_data['data'][i]['v']) if obs_data['data'][i]['v'] else np.nan for i in range(len(obs_data['data']))])
        date_start_i += timedelta(days=31)
        
    station_name = obs_data['metadata']['name']
    station_lon = float(obs_data['metadata']['lon'])
    station_lat = float(obs_data['metadata']['lat'])
    
    # Convert time list to pandas Series with UTC timezone
    obs_time = pd.Series(obs_time)
    obs_wl = pd.Series(obs_wl)
    
    # Save to cache if cache_dir is provided
    if cache_dir:
        cache_filename = _get_cache_filename('NOAA', station_id, date_start, date_end, datum)
        cache_path = os.path.join(cache_dir, cache_filename)
        print(f"Saving NOAA data to cache: {cache_path}")
        _save_cache(cache_path, station_name, station_lon, station_lat, obs_time, obs_wl)
    
    return station_name, station_lon, station_lat, obs_time, obs_wl

def _get_usgs_data(station_id, date_start, date_end, datum, **kwargs):
    """Local method to retrieve USGS water level data"""
    import pytz
    ft2m = 0.3048
    tzutc = pytz.timezone('UTC')
    
    # Check cache if cache_dir is provided
    cache_dir = kwargs.get('cache_dir')
    if cache_dir:
        cache_filename = _get_cache_filename('USGS', station_id, date_start, date_end, datum)
        cache_path = os.path.join(cache_dir, cache_filename)
        cached_data = _load_cache(cache_path)
        if cached_data is not None:
            print(f"Loading USGS data from cache: {cache_path}")
            return cached_data
    
    # Handle timezone-aware datetime objects
    date_start_naive = date_start.replace(tzinfo=None) if date_start.tzinfo else date_start
    date_end_naive = date_end.replace(tzinfo=None) if date_end.tzinfo else date_end
    date_start_str = date_start_naive.strftime('%Y-%m-%dT%H:%M')
    date_end_str = date_end_naive.strftime('%Y-%m-%dT%H:%M')
    dfst = nwis.get_record(sites=station_id, service='site')
    dfiv = nwis.get_record(sites=station_id, service='iv', start=date_start_str, end=date_end_str)
    station_name = dfst['station_nm'][0]
    station_lon = dfst['dec_long_va'][0]
    station_lat = dfst['dec_lat_va'][0]
    
    # Convert time to UTC timezone-aware
    obs_time = pd.to_datetime(dfiv.index)
    if obs_time.tz is None:
        # If timezone-naive, assume UTC
        obs_time = obs_time.tz_localize('UTC')
    else:
        # Convert to UTC if it has a different timezone
        obs_time = obs_time.tz_convert('UTC')
    
    print('alt_va = ', dfst['alt_va'][0])
    if '00065' in dfiv.columns:
        obs_wl = (dfiv['00065'] + dfst['alt_va'][0])*ft2m
    elif '62620' in dfiv.columns:
        obs_wl = dfiv['62620']*ft2m
    elif '62623' in dfiv.columns:
        obs_wl = dfiv['62623']*ft2m
    elif '00062' in dfiv.columns:
        obs_wl = dfiv['00062']*ft2m
    else:
        print(f"Available columns: {dfiv.columns}")
        raise KeyError('No valid column found in dfiv')
    
    # Save to cache if cache_dir is provided
    if cache_dir:
        cache_filename = _get_cache_filename('USGS', station_id, date_start, date_end, datum)
        cache_path = os.path.join(cache_dir, cache_filename)
        print(f"Saving USGS data to cache: {cache_path}")
        _save_cache(cache_path, station_name, station_lon, station_lat, obs_time, obs_wl)
    
    return station_name, station_lon, station_lat, obs_time, obs_wl

def _get_contrail_metadata(station_id, session, username, password):
    """Retrieve CONTRAIL station metadata including device IDs and coordinates"""
    metadata_url = f"https://contrail.nc.gov/site/?site_id={station_id}"
    
    try:
        # Try to access metadata page
        response = session.get(metadata_url, allow_redirects=True)
        
        # Check if we were redirected to login page
        if 'login' in response.url.lower():
            # Need to authenticate first
            soup = BeautifulSoup(response.text, 'html.parser')
            form = soup.find('form')
            
            if not form:
                raise ValueError("No login form found for metadata access")
            
            # Extract form data
            form_data = {}
            for input_tag in form.find_all('input'):
                name = input_tag.get('name')
                value = input_tag.get('value', '')
                input_type = input_tag.get('type', 'text')
                
                if name:
                    if input_type == 'hidden':
                        form_data[name] = value
                    elif name == 'username':
                        form_data[name] = username
                    elif name == 'password':
                        form_data[name] = password
            
            form_data['login'] = 'login'
            
            # Submit credentials
            login_headers = {
                'Referer': response.url,
                'Origin': 'https://contrail.nc.gov',
                'Content-Type': 'application/x-www-form-urlencoded',
            }
            
            auth_response = session.post(response.url, 
                                       data=form_data, 
                                       headers=login_headers,
                                       allow_redirects=True)
            
            # Now try to access metadata page again
            response = session.get(metadata_url, allow_redirects=True)
        
        if response.status_code != 200:
            raise ValueError(f"Failed to retrieve metadata: HTTP {response.status_code}")
        
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Extract sensor information
        sensors = {}
        sensor_links = soup.find_all('a', href=lambda x: x and 'sensor/' in x and 'device_id=' in x)
        
        for link in sensor_links:
            # Extract device_id from href
            href = link.get('href')
            device_id_match = href.split('device_id=')[1].split('&')[0] if 'device_id=' in href else None
            
            if device_id_match:
                device_id = int(device_id_match)
                sensor_text = link.get_text().strip().lower()
                
                # Map sensor types to device IDs
                if 'water elevation' in sensor_text:
                    sensors['water_elevation'] = device_id
                elif 'stage' in sensor_text:
                    sensors['stage'] = device_id
                elif 'stream elevation' in sensor_text:
                    sensors['stream_elevation'] = device_id
        
        # Extract coordinates from the map link
        station_lon = None
        station_lat = None
        coord_links = soup.find_all('a', href=lambda x: x and 'map/' in x and 'find_site_id=' in x)
        
        for link in coord_links:
            coord_text = link.get_text().strip()
            if ',' in coord_text:
                try:
                    lat_str, lon_str = coord_text.split(',')
                    station_lat = float(lat_str.strip())
                    station_lon = float(lon_str.strip())
                    break
                except ValueError:
                    continue
        
        # Extract station name from title
        station_name = "Unknown Station"
        title_element = soup.find('h3', id='title')
        if title_element:
            # Remove icon and small elements to get clean name
            for elem in title_element.find_all(['i', 'small']):
                elem.decompose()
            station_name = title_element.get_text().strip()
        
        return {
            'sensors': sensors,
            'station_lon': station_lon,
            'station_lat': station_lat,
            'station_name': station_name
        }
        
    except Exception as e:
        raise ValueError(f"Failed to parse CONTRAIL metadata: {e}")

def _get_vdatum_offset(station_lon, station_lat, source_datum='NAVD88', target_datum='LMSL', region='contiguous'):
    """
    Get datum offset from VDATUM API for converting between vertical datums.
    
    Parameters:
    -----------
    station_lon : float
        Station longitude
    station_lat : float
        Station latitude
    source_datum : str
        Source vertical reference frame (default: 'NAVD88')
    target_datum : str
        Target vertical reference frame (default: 'MSL')
    region : str
        VDATUM region (default: 'contiguous')
    
    Returns:
    --------
    float
        Datum offset in meters (target_datum = source_datum - offset)
        Returns None if API call fails
    """
    vdatum_url = 'https://vdatum.noaa.gov/vdatumweb/api/convert'
    
    params = {
        's_x': station_lon,
        's_y': station_lat,
        's_z': 0.0,  # Use zero to get the datum offset
        's_v_frame': source_datum,
        't_v_frame': target_datum,
        'region': region
    }
    
    try:
        response = requests.get(vdatum_url, params=params)
        print(f"VDATUME API URL: {response.url}")
        
        if response.status_code == 200:
            data = response.json()
            # The offset is the difference: t_z - s_z
            # Since s_z = 0, offset = t_z
            t_z = float(data.get('t_z', 0.0))
            
            if t_z == -999999.0:
                print("Warning: VDATUM API returned -999999.0. Datum offset is set to 0.0")
                t_z = 0.0
                
            return t_z
        else:
            print(f"Warning: VDATUM API returned status code {response.status_code}")
            return None
            
    except requests.exceptions.RequestException as e:
        print(f"Warning: Failed to retrieve datum offset from VDATUM API: {e}")
        return None
    except (KeyError, ValueError, TypeError) as e:
        print(f"Warning: Failed to parse VDATUM API response: {e}")
        return None

def _get_contrail_data(station_id, date_start, date_end, datum, **kwargs):
    """Local method to retrieve Contrail water level data"""
    # Extract credentials from options
    username = kwargs.get('username')
    password = kwargs.get('password')
    sensor_type = kwargs.get('sensor_type', 'auto')  # Default sensor type
    
    if not username or not password:
        raise ValueError("Contrail requires 'username' and 'password' in options")
    
    # Check cache if cache_dir is provided
    cache_dir = kwargs.get('cache_dir')
    if cache_dir:
        cache_filename = _get_cache_filename('CONTRAIL', station_id, date_start, date_end, datum)
        cache_path = os.path.join(cache_dir, cache_filename)
        cached_data = _load_cache(cache_path)
        if cached_data is not None:
            print(f"Loading CONTRAIL data from cache: {cache_path}")
            return cached_data
    
    # CONTRAIL datum validation
    # Note: CONTRAIL data typically comes in local reference datum (often NAVD88 for NC)
    # MSL conversion will be performed via VDATUM API when requested
    supported_datums = ['NAVD88', 'NAVD', 'MSL']
    if datum not in supported_datums:
        print(f"Warning: CONTRAIL data is returned in its native datum (typically NAVD88 for North Carolina).")
        print(f"Requested datum '{datum}' may not match the actual datum of the data.")
        print(f"Supported datums: {supported_datums}")
        if datum.upper() == 'MSL':
            print(f"MSL conversion will be attempted via VDATUM API.")
        else:
            print(f"Proceeding with data retrieval but datum transformation is NOT applied.")
    
    # Create session with realistic headers
    session = requests.Session()
    session.headers.update({
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Accept-Language': 'en-US,en;q=0.5',
        'Accept-Encoding': 'gzip, deflate',
        'DNT': '1',
        'Connection': 'keep-alive',
        'Upgrade-Insecure-Requests': '1',
    })
    
    # Get station metadata to find device_id and coordinates
    print(f"Retrieving CONTRAIL metadata for station {station_id}...")
    metadata = _get_contrail_metadata(station_id, session, username, password)
    
    # Map sensor type to device_id
    available_sensors = list(metadata['sensors'].keys())
    if sensor_type == 'auto':
        if 'water_elevation' in available_sensors:
            sensor_type = 'water_elevation'
        elif 'stream_elevation' in available_sensors:
            sensor_type = 'stream_elevation'
        else:
            raise ValueError(f"Automatic sensor type determination failed. Available sensors: {available_sensors}")
    if sensor_type not in metadata['sensors']:
        raise ValueError(f"Sensor type '{sensor_type}' not found. Available sensors: {available_sensors}")
    
    device_id = metadata['sensors'][sensor_type]
    station_name = metadata['station_name']
    station_lon = metadata['station_lon'] or kwargs.get('station_lon', -80.0)
    station_lat = metadata['station_lat'] or kwargs.get('station_lat', 35.0)
    
    print(f"Found sensor '{sensor_type}' with device_id={device_id}")
    print(f"Station: {station_name} at ({station_lat}, {station_lon})")
    
    # Build the export URL
    # Handle timezone-aware datetime objects
    date_start_naive = date_start.replace(tzinfo=None) if date_start.tzinfo else date_start
    date_end_naive = date_end.replace(tzinfo=None) if date_end.tzinfo else date_end
    date_start_str = date_start_naive.strftime('%Y-%m-%d %H:%M:%S')
    date_end_str = date_end_naive.strftime('%Y-%m-%d %H:%M:%S')
    
    source_tz = 'US/Eastern' # CONTRAIL data is in US/Eastern timezone
    
    url_params = {
        'site_id': station_id,
        'device_id': device_id,
        'hours': '',
        'data_start': date_start_str,
        'data_end': date_end_str,
        'tz': source_tz,
        'format_datetime': '%Y-%m-%d %H:%i:%S',
        'mime': 'txt',
        'delimiter': 'comma'
    }
    
    export_url = f"https://contrail.nc.gov/export/file/?{urlencode(url_params)}"
    print(f"Contrail export URL: {export_url}")
    
    try:
        # Step 1: Access the data URL (will redirect to login)
        response = session.get(export_url, allow_redirects=True)
        
        if 'login' in response.url.lower():
            # Step 2: Parse the login form
            soup = BeautifulSoup(response.text, 'html.parser')
            form = soup.find('form')
            
            if not form:
                raise ValueError("No login form found on redirected page")
            
            # Extract form data
            form_data = {}
            for input_tag in form.find_all('input'):
                name = input_tag.get('name')
                value = input_tag.get('value', '')
                input_type = input_tag.get('type', 'text')
                
                if name:
                    if input_type == 'hidden':
                        form_data[name] = value
                    elif name == 'username':
                        form_data[name] = username
                    elif name == 'password':
                        form_data[name] = password
            
            form_data['login'] = 'login'
            
            # Step 3: Submit credentials
            login_headers = {
                'Referer': response.url,
                'Origin': 'https://contrail.nc.gov',
                'Content-Type': 'application/x-www-form-urlencoded',
            }
            
            auth_response = session.post(response.url, 
                                       data=form_data, 
                                       headers=login_headers,
                                       allow_redirects=True)
            
            # Step 4: Check if we got data
            if export_url in auth_response.url or 'export/file' in auth_response.url:
                if not auth_response.text.strip().startswith('<'):
                    data = auth_response.text
                else:
                    # Try multipart form submission
                    files = {}
                    for key, value in form_data.items():
                        files[key] = (None, value)
                    
                    multipart_headers = {
                        'Referer': response.url,
                        'Origin': 'https://contrail.nc.gov',
                    }
                    
                    multipart_response = session.post(response.url, 
                                                   files=files, 
                                                   headers=multipart_headers,
                                                   allow_redirects=True)
                    
                    if (export_url in multipart_response.url or 'export/file' in multipart_response.url) and \
                       not multipart_response.text.strip().startswith('<'):
                        data = multipart_response.text
                    else:
                        raise ValueError("Failed to retrieve data after authentication")
            else:
                raise ValueError("Authentication failed - unexpected redirect")
        else:
            # Data returned without authentication
            if not response.text.strip().startswith('<'):
                data = response.text
            else:
                raise ValueError("No data received and no authentication required")
        
        # Parse the CSV data
        from io import StringIO
        df = pd.read_csv(StringIO(data))
        
        print(f"CONTRAIL CSV columns: {df.columns.tolist()}")
        print(f"CONTRAIL data shape: {df.shape}")
        
        # CONTRAIL specific column format: Reading,Receive,Value,Unit,Data Quality
        if 'Reading' not in df.columns or 'Value' not in df.columns:
            raise ValueError(f"Expected CONTRAIL columns 'Reading' and 'Value' not found. Available columns: {df.columns.tolist()}")
        
        # Convert time column (Reading) to datetime with UTC timezone
        import pytz
        obs_time_raw = pd.to_datetime(df['Reading'])
        obs_time = obs_time_raw.copy()
        
        # Ensure time is UTC timezone-aware
        if obs_time.dt.tz is None:
            # CONTRAIL data requested in UTC, so localize as UTC
            # Handle ambiguous times during DST transitions by marking them as NaT
            # (they will be filtered out later in the valid_data check)
            # During DST fall-back (e.g., Nov 2, 2025 2:00 AM -> 1:00 AM),
            # times between 1:00-2:00 AM occur twice and are ambiguous
            obs_time = obs_time.dt.tz_localize(source_tz, ambiguous='NaT')
            
            # Check if any ambiguous times were marked as NaT
            ambiguous_count = obs_time.isna().sum() - obs_time_raw.isna().sum()
            if ambiguous_count > 0:
                print(f"Warning: {ambiguous_count} ambiguous time(s) during DST transition marked as NaT and will be excluded")
            
            obs_time = obs_time.dt.tz_convert('UTC')
        else:
            # Convert to UTC if it has a different timezone
            obs_time = obs_time.dt.tz_convert('UTC')
        
        # Get water level values from Value column
        obs_wl = pd.to_numeric(df['Value'], errors='coerce')
        
        # Check for and remove NaN values
        valid_data = ~(obs_time.isna() | obs_wl.isna())
        if not valid_data.any():
            raise ValueError("No valid data found after parsing")
        
        obs_time = obs_time[valid_data]
        obs_wl = obs_wl[valid_data]
        
        print(f"Valid data points: {len(obs_time)} (after removing NaN values)")
        
        # Check if we have Unit information for conversion
        if 'Unit' in df.columns:
            unit = df['Unit'].iloc[0] if len(df) > 0 else 'unknown'
            print(f"CONTRAIL data unit: {unit}")
            
            # Convert from feet to meters if needed (CONTRAIL typically uses feet)
            if unit.lower() in ['ft', 'feet']:
                ft2m = 0.3048
                obs_wl = obs_wl * ft2m
                print(f"Converted from feet to meters (factor: {ft2m})")
        
        # Apply datum conversion if MSL is requested
        datum_offset = None  # Initialize to None
        if datum.upper() == 'MSL':
            # Validate that station coordinates are available
            if station_lon is None or station_lat is None:
                print(f"Warning: Station coordinates not available. Cannot convert to MSL.")
                print(f"Returning data in native NAVD88 datum.")
            else:
                print(f"Converting from NAVD88 to local MSL using VDATUM API...")
                datum_offset = _get_vdatum_offset(station_lon, station_lat, source_datum='NAVD88', target_datum='LMSL')
                
                if datum_offset is not None:
                    # Apply offset: MSL = NAVD88 + offset
                    obs_wl = obs_wl + datum_offset
                    print(f"Applied datum offset: {datum_offset:.4f} m (converted to local MSL)")
                else:
                    print(f"Warning: Failed to retrieve datum offset from VDATUM API.")
                    print(f"Returning data in native NAVD88 datum.")
        
        # CONTRAIL data is in descending time order - reverse to make it ascending
        # (consistent with other data sources)
        obs_time = obs_time[::-1].reset_index(drop=True)
        obs_wl = obs_wl[::-1].reset_index(drop=True)
        
        print(f"Data time range: {obs_time.iloc[0]} to {obs_time.iloc[-1]}")
        print(f"Water level range: {obs_wl.min():.3f} to {obs_wl.max():.3f} m")
        
        # Save to cache if cache_dir is provided
        if cache_dir:
            cache_filename = _get_cache_filename('CONTRAIL', station_id, date_start, date_end, datum)
            cache_path = os.path.join(cache_dir, cache_filename)
            print(f"Saving CONTRAIL data to cache: {cache_path}")
            _save_cache(cache_path, station_name, station_lon, station_lat, obs_time, obs_wl, datum_offset=datum_offset)
        
        return station_name, station_lon, station_lat, obs_time, obs_wl
        
    except Exception as e:
        print(f"Error retrieving Contrail data: {e}")
        return None, None, None, None, None

def _get_secoora_data(station_id, date_start, date_end, datum, **kwargs):
    """Local method to retrieve SECOORA water level data"""
    if datum != 'NAVD':
        raise ValueError('SECOORA only supports NAVD datum')
    
    # Check cache if cache_dir is provided
    cache_dir = kwargs.get('cache_dir')
    if cache_dir:
        cache_filename = _get_cache_filename('SECOORA', station_id, date_start, date_end, datum)
        cache_path = os.path.join(cache_dir, cache_filename)
        cached_data = _load_cache(cache_path)
        if cached_data is not None:
            print(f"Loading SECOORA data from cache: {cache_path}")
            return cached_data
    
    # Handle timezone-aware datetime objects
    date_start_naive = date_start.replace(tzinfo=None) if date_start.tzinfo else date_start
    date_end_naive = date_end.replace(tzinfo=None) if date_end.tzinfo else date_end
    date_start_str = date_start_naive.strftime('%Y-%m-%dT%H:%M')
    date_end_str = date_end_naive.strftime('%Y-%m-%dT%H:%M')

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
            water_level_vars = []
            for attr in metadata_json['table']['rows']:
                if attr[0] == 'variable':
                    var_name = attr[1]
                    if 'water_surface_above' in var_name.lower() or 'sea_surface_height' in var_name.lower():
                        water_level_vars.append(var_name)
            
            if len(water_level_vars) == 0:
                for attr in metadata_json['table']['rows']:
                    if attr[0] == 'variable':
                        print(f"Available variables: {attr[1]}")
                raise ValueError(f"Couldn't find water level variable for station {station_id}")
            
            if len(water_level_vars) == 1:
                water_level_var = water_level_vars[0]
            else:
                water_level_var = water_level_vars[0]
                for var in water_level_vars:
                    if 'navd' in var.lower():
                        water_level_var = var
                        break
            
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
        
        # Ensure time is UTC timezone-aware
        if obs_time.dt.tz is None:
            # SECOORA data is typically in UTC, so localize as UTC
            obs_time = obs_time.dt.tz_localize('UTC')
        else:
            # Convert to UTC if it has a different timezone
            obs_time = obs_time.dt.tz_convert('UTC')
        
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
        
        # Save to cache if cache_dir is provided
        if cache_dir:
            cache_filename = _get_cache_filename('SECOORA', station_id, date_start, date_end, datum)
            cache_path = os.path.join(cache_dir, cache_filename)
            print(f"Saving SECOORA data to cache: {cache_path}")
            _save_cache(cache_path, station_name, station_lon, station_lat, obs_time, obs_wl)
        
        return station_name, station_lon, station_lat, obs_time, obs_wl
        
    except Exception as e:
        print(f"Error processing SECOORA data: {str(e)}")
        return None, None, None, None, None

def get_obswl(station_owner, station_id, date_start, date_end, datum, options=None, cache_dir=None):
    """
    Get observed water level data from various sources.
    
    Parameters:
    -----------
    station_owner : str
        Source of the data ('NOAA', 'USGS', 'CONTRAIL', 'SECOORA')
    station_id : str
        Station identifier
    date_start : datetime
        Start date for data retrieval (assumed to be in UTC)
    date_end : datetime
        End date for data retrieval (assumed to be in UTC)
    datum : str
        Datum for water level measurements
    options : dict, optional
        Additional options for specific data sources:
        - For CONTRAIL: 'username', 'password', 'sensor_type' (water_elevation, stream_elevation, stage)
          Note: CONTRAIL only supports NAVD88/NAVD datum; other datums will generate a warning
          Station coordinates are automatically extracted from metadata
        - For other sources: additional parameters as needed
    cache_dir : str or Path, optional
        Directory path for caching downloaded observation data in JSON format.
        If provided, data will be saved to and loaded from cache files.
        Cache files are named: {owner}_{sanitized_station}_{start_date}_{end_date}.json
    
    Returns:
    --------
    tuple
        (station_name, station_lon, station_lat, obs_time, obs_wl)
        obs_time will be timezone-aware pandas Series in UTC
    """
    if options is None:
        options = {}
    
    # Add cache_dir to options if provided
    if cache_dir is not None:
        options['cache_dir'] = cache_dir
    
    # Dispatch to appropriate local method
    if station_owner == 'NOAA':
        return _get_noaa_data(station_id, date_start, date_end, datum, **options)
    elif station_owner == 'USGS':
        return _get_usgs_data(station_id, date_start, date_end, datum, **options)
    elif station_owner == 'CONTRAIL':
        return _get_contrail_data(station_id, date_start, date_end, datum, **options)
    elif station_owner == 'SECOORA':
        return _get_secoora_data(station_id, date_start, date_end, datum, **options)
    else:
        raise ValueError(f'Invalid station owner: {station_owner}. Valid options are: NOAA, USGS, CONTRAIL, SECOORA')

def get_parser():
    """Get argument parser for get_obswl command"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Retrieve observed water level data from various sources',
        add_help=False
    )
    
    # Required arguments
    parser.add_argument(
        'station_owner',
        choices=['NOAA', 'USGS', 'CONTRAIL', 'SECOORA'],
        help='Data source (NOAA, USGS, CONTRAIL, SECOORA)'
    )
    parser.add_argument(
        'station_id',
        help='Station identifier'
    )
    parser.add_argument(
        'date_start',
        help='Start date (YYYY-MM-DD or YYYY-MM-DD HH:MM:SS)'
    )
    parser.add_argument(
        'date_end',
        help='End date (YYYY-MM-DD or YYYY-MM-DD HH:MM:SS)'
    )
    parser.add_argument(
        'datum',
        nargs='?',
        help='Datum for water level measurements (e.g., MLLW, NAVD, MSL). If not specified, uses default for each source: NOAA=MLLW, USGS=NAVD88, CONTRAIL=NAVD88, SECOORA=NAVD'
    )
    
    # Optional arguments
    parser.add_argument(
        '-o', '--output',
        help='Output file path (CSV format). If not specified, prints to stdout'
    )
    parser.add_argument(
        '--format',
        choices=['csv', 'json', 'summary'],
        default='csv',
        help='Output format (default: csv)'
    )
    parser.add_argument(
        '--cache-dir',
        help='Directory path for caching downloaded observation data in JSON format'
    )
    
    # CONTRAIL specific options
    contrail_group = parser.add_argument_group('CONTRAIL options')
    contrail_group.add_argument(
        '--username',
        help='Username for CONTRAIL authentication'
    )
    contrail_group.add_argument(
        '--password',
        help='Password for CONTRAIL authentication'
    )
    contrail_group.add_argument(
        '--sensor-type',
        choices=['water_elevation', 'stream_elevation', 'stage'],
        default='water_elevation',
        help='Sensor type for CONTRAIL (default: water_elevation)'
    )
    
    return parser

def main(args):
    """Main function for CLI"""
    import json
    import sys
    from datetime import datetime
    
    # Parse date strings
    def parse_date(date_str):
        """Parse date string in various formats and return as UTC datetime"""
        import pytz
        formats = [
            '%Y-%m-%d',
            '%Y-%m-%d %H:%M:%S',
            '%Y-%m-%dT%H:%M:%S'
        ]
        for fmt in formats:
            try:
                dt = datetime.strptime(date_str, fmt)
                # Assume input dates are in UTC
                return pytz.UTC.localize(dt)
            except ValueError:
                continue
        raise ValueError(f"Could not parse date: {date_str}. Expected format: YYYY-MM-DD or YYYY-MM-DD HH:MM:SS")
    
    # Get default datum for each data source
    def get_default_datum(station_owner):
        """Get default datum for each data source"""
        defaults = {
            'NOAA': 'MLLW',      # NOAA typically uses MLLW for tidal stations
            'USGS': 'NAVD88',    # USGS typically uses NAVD88 for stream gauges
            'CONTRAIL': 'NAVD88', # CONTRAIL (North Carolina) uses NAVD88
            'SECOORA': 'NAVD'    # SECOORA only supports NAVD
        }
        return defaults.get(station_owner, 'NAVD88')
    
    try:
        date_start = parse_date(args.date_start)
        date_end = parse_date(args.date_end)
        print(f"Date range (UTC): {date_start} to {date_end}", file=sys.stderr)
    except ValueError as e:
        print(f"Error parsing dates: {e}", file=sys.stderr)
        return 1
    
    # Handle optional datum argument
    datum = args.datum
    if datum is None:
        datum = get_default_datum(args.station_owner)
        print(f"Using default datum for {args.station_owner}: {datum}", file=sys.stderr)
    
    # Build options dictionary
    options = {}
    if args.station_owner == 'CONTRAIL':
        if not args.username or not args.password:
            print("Error: CONTRAIL requires --username and --password", file=sys.stderr)
            return 1
        options['username'] = args.username
        options['password'] = args.password
        options['sensor_type'] = args.sensor_type
    
    # Retrieve data
    try:
        print(f"Retrieving data from {args.station_owner} for station {args.station_id}...", file=sys.stderr)
        station_name, station_lon, station_lat, obs_time, obs_wl = get_obswl(
            args.station_owner,
            args.station_id,
            date_start,
            date_end,
            datum,
            options,
            cache_dir=args.cache_dir
        )
        
        if station_name is None:
            print("Error: No data retrieved", file=sys.stderr)
            return 1
        
        print(f"Retrieved {len(obs_time)} data points for {station_name}", file=sys.stderr)
        
        # Format output
        if args.format == 'summary':
            output = {
                'station_name': station_name,
                'station_lon': station_lon,
                'station_lat': station_lat,
                'data_points': len(obs_time),
                'time_range': {
                    'start': obs_time.min().isoformat() if hasattr(obs_time, 'min') else str(obs_time[0]),
                    'end': obs_time.max().isoformat() if hasattr(obs_time, 'max') else str(obs_time[-1])
                },
                'water_level_range': {
                    'min': float(obs_wl.min()) if hasattr(obs_wl, 'min') else float(min(obs_wl)),
                    'max': float(obs_wl.max()) if hasattr(obs_wl, 'max') else float(max(obs_wl))
                }
            }
            output_str = json.dumps(output, indent=2)
        elif args.format == 'json':
            # Convert to JSON format
            data = []
            for i in range(len(obs_time)):
                time_val = obs_time.iloc[i] if hasattr(obs_time, 'iloc') else obs_time[i]
                wl_val = obs_wl.iloc[i] if hasattr(obs_wl, 'iloc') else obs_wl[i]
                data.append({
                    'time': time_val.isoformat() if hasattr(time_val, 'isoformat') else str(time_val),
                    'water_level': float(wl_val) if not pd.isna(wl_val) else None
                })
            output = {
                'station_name': station_name,
                'station_lon': station_lon,
                'station_lat': station_lat,
                'datum': datum,
                'data': data
            }
            output_str = json.dumps(output, indent=2)
        else:  # CSV format
            output_lines = ['time,water_level']
            for i in range(len(obs_time)):
                time_val = obs_time.iloc[i] if hasattr(obs_time, 'iloc') else obs_time[i]
                wl_val = obs_wl.iloc[i] if hasattr(obs_wl, 'iloc') else obs_wl[i]
                time_str = time_val.isoformat() if hasattr(time_val, 'isoformat') else str(time_val)
                wl_str = f"{wl_val:.6f}" if not pd.isna(wl_val) else ""
                output_lines.append(f"{time_str},{wl_str}")
            output_str = '\n'.join(output_lines)
        
        # Write output
        if args.output:
            with open(args.output, 'w') as f:
                f.write(output_str)
            print(f"Data written to {args.output}", file=sys.stderr)
        else:
            print(output_str)
        
        return 0
        
    except Exception as e:
        print(f"Error retrieving data: {e}", file=sys.stderr)
        return 1

if __name__ == '__main__':
    import sys
    import argparse
    
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(main(args))