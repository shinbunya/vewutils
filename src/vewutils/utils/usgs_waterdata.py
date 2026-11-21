"""Client for the USGS Water Data OGC API (api.waterdata.usgs.gov).

Replaces the legacy NWIS WaterServices (waterservices.usgs.gov), whose
servers have become slow and occasionally accept a request without ever
answering it. Every request here has a timeout so a stalled server fails
the call instead of hanging the caller.

An API key from https://api.waterdata.usgs.gov/signup/ raises the request
rate limit; set it in the API_USGS_PAT environment variable (the name the
dataretrieval package also uses).
"""
import os
import time

import pandas as pd
import requests

BASE_URL = 'https://api.waterdata.usgs.gov/ogcapi/v1'
REQUEST_TIMEOUT_SEC = (15, 60)  # (connect, read)
PAGE_LIMIT = 10000
MAX_ATTEMPTS = 3

_session = None


def _get_session():
    global _session
    if _session is None:
        _session = requests.Session()
        api_key = os.environ.get('API_USGS_PAT')
        if api_key:
            _session.headers['X-Api-Key'] = api_key
    return _session


def monitoring_location_id(site_id):
    """Return the 'USGS-<number>' id for a site number.

    Site numbers shorter than 8 digits had their leading zero dropped
    somewhere upstream (e.g. 2136350 for 02136350), so they are padded.
    """
    site = str(site_id).strip()
    if site.upper().startswith('USGS-'):
        site = site[5:]
    if site.isdigit() and len(site) < 8:
        site = site.zfill(8)
    return 'USGS-' + site


def _get_json(url, params=None):
    """GET a URL and return its JSON, retrying timeouts and server errors.

    Returns None for a 404 (unknown id). Raises on anything else that keeps
    failing.
    """
    session = _get_session()
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            response = session.get(url, params=params, timeout=REQUEST_TIMEOUT_SEC)
            if response.status_code == 404:
                return None
            if response.status_code == 429 or response.status_code >= 500:
                raise requests.HTTPError(
                    f'{response.status_code} from {response.url}', response=response
                )
            response.raise_for_status()
            return response.json()
        except (requests.Timeout, requests.ConnectionError, requests.HTTPError) as exc:
            status = getattr(getattr(exc, 'response', None), 'status_code', None)
            if status is not None and status != 429 and status < 500:
                raise
            if attempt == MAX_ATTEMPTS:
                raise
            retry_after = 5 * attempt
            if status == 429:
                try:
                    retry_after = int(exc.response.headers.get('Retry-After', 60))
                except ValueError:
                    retry_after = 60
            print(f'USGS request failed ({exc}); retrying in {retry_after} s')
            time.sleep(retry_after)


def get_monitoring_location(site_id):
    """Return the monitoring-location properties for a site, or None if unknown.

    Notable keys: monitoring_location_name, altitude (gage datum, ft),
    vertical_datum, plus the point geometry as 'longitude' and 'latitude'.
    """
    data = _get_json(
        f'{BASE_URL}/collections/monitoring-locations/items/{monitoring_location_id(site_id)}',
        params={'f': 'json'},
    )
    if data is None:
        return None
    props = dict(data.get('properties', {}))
    coords = (data.get('geometry') or {}).get('coordinates') or [None, None]
    props['longitude'], props['latitude'] = coords[0], coords[1]
    return props


def _primary_time_series_ids(site_id, parameter_code):
    data = _get_json(
        f'{BASE_URL}/collections/time-series-metadata/items',
        params={
            'monitoring_location_id': monitoring_location_id(site_id),
            'parameter_code': parameter_code,
            'computation_period_identifier': 'Points',
            'skipGeometry': 'true',
            'f': 'json',
        },
    )
    if not data:
        return set()
    return {
        f['id']
        for f in data.get('features', [])
        if f.get('properties', {}).get('primary') == 'Primary'
    }


def get_continuous(site_id, parameter_codes, start, end):
    """Return instantaneous values of one or more parameters.

    Returns a dict mapping each parameter code that has data in
    [start, end] to a pandas Series indexed by UTC time, holding float
    values in the parameter's native unit (ft, ft^3/s, ...). Codes without
    data are left out. When the site has several series of a parameter
    (e.g. two sensors), the one USGS marks primary is used.
    """
    if isinstance(parameter_codes, str):
        parameter_codes = [parameter_codes]
    start_utc = pd.Timestamp(start)
    end_utc = pd.Timestamp(end)
    start_utc = start_utc.tz_localize('UTC') if start_utc.tzinfo is None else start_utc.tz_convert('UTC')
    end_utc = end_utc.tz_localize('UTC') if end_utc.tzinfo is None else end_utc.tz_convert('UTC')

    url = f'{BASE_URL}/collections/continuous/items'
    params = {
        'monitoring_location_id': monitoring_location_id(site_id),
        'parameter_code': ','.join(parameter_codes),
        'time': f"{start_utc.strftime('%Y-%m-%dT%H:%M:%SZ')}/{end_utc.strftime('%Y-%m-%dT%H:%M:%SZ')}",
        'properties': 'time,value,parameter_code,time_series_id',
        'skipGeometry': 'true',
        'limit': PAGE_LIMIT,
        'f': 'json',
    }
    records = []
    while url:
        data = _get_json(url, params=params)
        if not data:
            break
        for feature in data.get('features', []):
            records.append(feature['properties'])
        # The 'next' link carries every query parameter plus a cursor
        url = next(
            (link['href'] for link in data.get('links', []) if link.get('rel') == 'next'),
            None,
        )
        params = None

    result = {}
    if not records:
        return result
    df = pd.DataFrame.from_records(records)
    for code, group in df.groupby('parameter_code'):
        series_ids = group['time_series_id'].unique()
        if len(series_ids) > 1:
            primary = _primary_time_series_ids(site_id, code) & set(series_ids)
            if primary:
                chosen = sorted(primary)[0]
            else:
                chosen = group['time_series_id'].value_counts().index[0]
            group = group[group['time_series_id'] == chosen]
        values = pd.to_numeric(group['value'], errors='coerce')
        index = pd.DatetimeIndex(pd.to_datetime(group['time'], utc=True))
        result[code] = pd.Series(values.to_numpy(), index=index, dtype=float).sort_index()
    return result
