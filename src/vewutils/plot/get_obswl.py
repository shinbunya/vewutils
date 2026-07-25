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
import time
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

CONTRAIL_LIST_URL = 'https://contrail.nc.gov/list/'
CONTRAIL_STATION_LIST_CACHE = 'contrail_station_list.json'

# None of these sources have observations for the future; a window that
# reaches past "now" routinely comes back as an error-shaped response
# instead of the partial (historical) data it should still return (see
# _fetch_noaa_range/_fetch_usgs_range). Clip the requested end time to a
# few minutes behind the current time before ever making a request --
# confirmed empirically that both NOAA and USGS return real data for a
# window ending this close to the present.
OBS_FUTURE_MARGIN = timedelta(minutes=5)


def _clip_to_present(date_end):
    """Clip date_end so it never reaches past OBS_FUTURE_MARGIN before now."""
    now = datetime.now(pytz.UTC)
    cutoff = now if date_end.tzinfo else now.replace(tzinfo=None)
    cutoff = cutoff - OBS_FUTURE_MARGIN
    return min(date_end, cutoff)

def _create_contrail_session():
    """Create a requests session with browser-like headers for CONTRAIL."""
    session = requests.Session()
    session.headers.update({
        'User-Agent': (
            'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 '
            '(KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        ),
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Accept-Language': 'en-US,en;q=0.5',
        'Accept-Encoding': 'gzip, deflate',
        'DNT': '1',
        'Connection': 'keep-alive',
        'Upgrade-Insecure-Requests': '1',
    })
    return session

def _contrail_login(session, login_response, username, password):
    """Submit CONTRAIL login form from a login page response."""
    soup = BeautifulSoup(login_response.text, 'html.parser')
    form = soup.find('form')
    if not form:
        raise ValueError("No login form found on CONTRAIL login page")

    form_data = {}
    for input_tag in form.find_all('input'):
        name = input_tag.get('name')
        value = input_tag.get('value', '')
        input_type = input_tag.get('type', 'text')
        if not name:
            continue
        if input_type == 'hidden':
            form_data[name] = value
        elif name == 'username':
            form_data[name] = username
        elif name == 'password':
            form_data[name] = password

    form_data['login'] = 'login'
    login_headers = {
        'Referer': login_response.url,
        'Origin': 'https://contrail.nc.gov',
        'Content-Type': 'application/x-www-form-urlencoded',
    }
    session.post(
        login_response.url,
        data=form_data,
        headers=login_headers,
        allow_redirects=True,
    )

def _contrail_fetch_url(session, url, username, password):
    """GET a CONTRAIL URL, authenticating via login redirect when required."""
    response = session.get(url, allow_redirects=True)
    if 'login' in response.url.lower():
        _contrail_login(session, response, username, password)
        response = session.get(url, allow_redirects=True)
    if response.status_code != 200:
        raise ValueError(f"Failed to retrieve {url}: HTTP {response.status_code}")
    if 'login' in response.url.lower():
        raise ValueError(f"CONTRAIL authentication failed for {url}")
    return response

def _parse_contrail_station_list_html(html):
    """
    Parse CONTRAIL /list/ HTML into station records.

    Each record has site_id (str), code (str or None), and name (str or None).
    """
    soup = BeautifulSoup(html, 'html.parser')
    stations = []
    for heading in soup.find_all('h4', class_='list-group-item-heading'):
        link = heading.find('a', href=lambda x: x and 'site_id=' in x)
        if not link:
            continue
        href = link.get('href', '')
        site_match = re.search(r'site_id=(\d+)', href)
        if not site_match:
            continue
        site_id = site_match.group(1)
        name = link.get_text(strip=True)
        code = None
        small = heading.find('small')
        if small:
            # CONTRAIL's fort.61-style codes routinely contain underscores
            # (e.g. "DE_03", "BF_01") or hyphens (e.g. "RUTHE-029"); the
            # previous [A-Za-z0-9]+ class silently dropped every one of
            # those, leaving `code` None even though CONTRAIL does publish
            # a code for the station.
            code_match = re.search(r'\(([\w.-]+)\)', small.get_text())
            if code_match:
                code = code_match.group(1).upper()
        stations.append({'site_id': site_id, 'code': code, 'name': name})
    return stations

def _build_contrail_station_map(stations):
    """Build lookup dicts from parsed CONTRAIL station list records."""
    by_site_id = {}
    by_code = {}
    for station in stations:
        site_id = station['site_id']
        by_site_id[site_id] = station
        code = station.get('code')
        if code:
            by_code[code] = station
    return {
        'fetched_at': datetime.utcnow().isoformat() + 'Z',
        'stations': stations,
        'by_site_id': by_site_id,
        'by_code': by_code,
    }

def _get_contrail_station_list_cache_path(cache_dir):
    return os.path.join(cache_dir, CONTRAIL_STATION_LIST_CACHE)

def _load_contrail_station_list_cache(cache_path):
    if not os.path.isfile(cache_path):
        return None
    try:
        with open(cache_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        if 'by_site_id' in data and 'by_code' in data:
            print(f"Loading CONTRAIL station list from cache: {cache_path}")
            return data
    except (json.JSONDecodeError, OSError) as e:
        print(f"Warning: Could not load CONTRAIL station list cache: {e}")
    return None

def _save_contrail_station_list_cache(cache_path, station_map):
    cache_dir = os.path.dirname(cache_path)
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
    with open(cache_path, 'w', encoding='utf-8') as f:
        json.dump(station_map, f, indent=2)
    print(f"Saved CONTRAIL station list to cache: {cache_path}")

def get_contrail_station_map(username, password, cache_dir=None, force_refresh=False):
    """
    Return CONTRAIL station metadata keyed by site_id and fort.61-style code.

    Downloads https://contrail.nc.gov/list/ when no valid cache is available.
    """
    if not username or not password:
        raise ValueError("CONTRAIL station list requires 'username' and 'password'")

    cache_path = None
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        cache_path = _get_contrail_station_list_cache_path(cache_dir)
        if not force_refresh:
            cached = _load_contrail_station_list_cache(cache_path)
            if cached is not None:
                return cached

    print("Retrieving CONTRAIL station list...")
    session = _create_contrail_session()
    response = _contrail_fetch_url(session, CONTRAIL_LIST_URL, username, password)
    stations = _parse_contrail_station_list_html(response.text)
    if not stations:
        raise ValueError(
            "No CONTRAIL stations parsed from list page; "
            "the page layout may have changed or authentication failed"
        )
    station_map = _build_contrail_station_map(stations)
    print(f"Found {len(stations)} CONTRAIL stations ({len(station_map['by_code'])} with codes)")
    if cache_path:
        _save_contrail_station_list_cache(cache_path, station_map)
    return station_map

def resolve_contrail_station_ids(
        station_id,
        station_id_type=None,
        username=None,
        password=None,
        cache_dir=None,
        f61_station_id=None,
        station_map=None):
    """
    Resolve CONTRAIL site_id (integer string) and fort.61 station code.

    Parameters
    ----------
    station_id : str
        Station identifier. May be ``contrail_site_id/f61_code`` (e.g. ``1205/EGHN7``).
    station_id_type : str, optional
        How to interpret ``station_id`` when it does not contain ``/``:

        - ``None`` (legacy): use ``station_id`` for both CONTRAIL and fort.61 unless
          lookup is required (see ``auto``).
        - ``'auto'``: numeric ids are CONTRAIL site ids; otherwise treat as fort.61 code.
        - ``'contrail'``: ``station_id`` is a CONTRAIL site id (lookup fort.61 code).
        - ``'f61'``: ``station_id`` is a fort.61 / ADCIRC station code (lookup site id).
    f61_station_id : str, optional
        Explicit fort.61 station code; overrides lookup for the model side only.
    station_map : dict, optional
        Pre-loaded map from :func:`get_contrail_station_map`.

    Returns
    -------
    dict
        ``contrail_site_id``, ``f61_station_id``, and ``display_station_id`` (for titles).
    """
    if not station_id:
        raise ValueError("station_id is required for CONTRAIL")

    station_id = str(station_id).strip()
    if f61_station_id is not None:
        f61_station_id = str(f61_station_id).strip()

    if '/' in station_id:
        contrail_part, f61_part = station_id.split('/', 1)
        contrail_site_id = contrail_part.strip()
        f61_from_slash = f61_part.strip()
        if not contrail_site_id or not f61_from_slash:
            raise ValueError(
                "CONTRAIL combined station_id must be 'contrail_site_id/f61_code', "
                f"got {station_id!r}"
            )
        f61_station_id = f61_station_id or f61_from_slash
        return {
            'contrail_site_id': contrail_site_id,
            'f61_station_id': f61_station_id,
            'display_station_id': f61_station_id,
        }

    id_type = (station_id_type or 'legacy').lower()
    if id_type == 'legacy':
        if f61_station_id is not None:
            return {
                'contrail_site_id': station_id,
                'f61_station_id': f61_station_id,
                'display_station_id': f61_station_id,
            }
        return {
            'contrail_site_id': station_id,
            'f61_station_id': station_id,
            'display_station_id': station_id,
        }

    needs_map = id_type in ('auto', 'contrail', 'f61')
    if needs_map and station_map is None:
        station_map = get_contrail_station_map(username, password, cache_dir=cache_dir)

    if id_type == 'auto':
        id_type = 'contrail' if station_id.isdigit() else 'f61'

    if id_type == 'contrail':
        contrail_site_id = station_id
        entry = station_map['by_site_id'].get(contrail_site_id)
        if entry and entry.get('code'):
            resolved_f61 = entry['code']
        elif f61_station_id is not None:
            resolved_f61 = f61_station_id
        else:
            resolved_f61 = station_id
            print(
                f"Warning: No fort.61 code found for CONTRAIL site {contrail_site_id}; "
                f"using {resolved_f61!r} for fort.61 lookup"
            )
    elif id_type == 'f61':
        code = station_id.upper()
        entry = station_map['by_code'].get(code)
        if not entry:
            raise ValueError(
                f"Unknown CONTRAIL fort.61 station code {code!r}. "
                "Refresh the station list cache or verify the code on "
                "https://contrail.nc.gov/list/"
            )
        contrail_site_id = entry['site_id']
        resolved_f61 = f61_station_id or code
    else:
        raise ValueError(
            f"Invalid station_id_type: {station_id_type!r}. "
            "Valid values: None, 'auto', 'contrail', 'f61'"
        )

    if f61_station_id is not None:
        resolved_f61 = f61_station_id

    return {
        'contrail_site_id': contrail_site_id,
        'f61_station_id': resolved_f61,
        'display_station_id': resolved_f61,
    }

def _get_daily_cache_filename(station_owner, station_id, datum, day):
    """
    Generate the cache filename for one calendar day's worth of data.

    Format: {owner}_{sanitized_station}_{datum}_{YYYYMMDD}.json

    Parameters:
    -----------
    station_owner : str
        Source of the data ('NOAA', 'USGS', 'CONTRAIL', 'SECOORA')
    station_id : str
        Station identifier
    datum : str
        Datum for water level measurements
    day : date
        The calendar day (UTC) this cache file holds

    Returns:
    --------
    str
        Cache filename
    """
    sanitized_station = _sanitize_station_id(station_id)
    sanitized_datum = _sanitize_station_id(datum)
    day_str = day.strftime('%Y%m%d')
    return f"{station_owner.upper()}_{sanitized_station}_{sanitized_datum}_{day_str}.json"

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

def _get_data_with_daily_cache(station_owner, station_id, date_start, date_end, datum, fetch_range_fn, cache_dir=None):
    """
    Fetch water level data for [date_start, date_end] using a per-day on-disk cache.

    Cache files are stored one per calendar day (UTC), keyed by
    (station_owner, station_id, datum, day). A day is only read from or written
    to cache when [date_start, date_end] fully spans it (from that day's
    midnight through the next); a day only partially covered by the requested
    window (typically the first or last day of the range) is always fetched
    fresh and never cached. Data with missing/NaN points within an otherwise
    fully-spanned day is still cached as-is -- only the query coverage matters,
    not the completeness of what came back.

    Contiguous runs of uncached, fully-spanned days are fetched together in one
    fetch_range_fn call (then split and cached per day), so a request doesn't
    turn into one network call per day. A short pause follows each such live
    call to avoid hammering these APIs with back-to-back requests across the
    many stations in a report.

    A day (or the boundary fragment) that comes back with no data -- most
    commonly because the window reaches past what the source has published
    yet -- is treated as "no data for this piece" rather than failing the
    whole request: whatever other days were already gathered (from cache or
    earlier fetches) are kept instead of being discarded.

    Parameters:
    -----------
    station_owner : str
        Source of the data ('NOAA', 'USGS', 'CONTRAIL', 'SECOORA'); used only
        for the cache filename.
    station_id : str
        Station identifier; used only for the cache filename.
    date_start, date_end : datetime
        Requested window (assumed UTC).
    datum : str
        Datum for water level measurements; used only for the cache filename.
    fetch_range_fn : callable
        fetch_range_fn(sub_date_start, sub_date_end) -> (station_name,
        station_lon, station_lat, obs_time, obs_wl), matching the return
        contract of the source-specific fetchers (a None station_name signals
        no data available for that sub-range, not necessarily an error).
    cache_dir : str or Path, optional
        Directory for per-day cache files. If not given, caching is skipped
        entirely and fetch_range_fn is called once for the whole window.

    Returns:
    --------
    tuple
        (station_name, station_lon, station_lat, obs_time, obs_wl)
    """
    if not cache_dir:
        return fetch_range_fn(date_start, date_end)

    date_start_naive = date_start.replace(tzinfo=None) if getattr(date_start, 'tzinfo', None) else date_start
    date_end_naive = date_end.replace(tzinfo=None) if getattr(date_end, 'tzinfo', None) else date_end

    # For each calendar day touched by the request, work out its window, whether
    # the request fully spans it, and (if so) whether it's already cached.
    day_infos = []
    day = date_start_naive.date()
    while day <= date_end_naive.date():
        day_start = datetime.combine(day, datetime.min.time())
        if day_start >= date_end_naive:
            # No overlap with the requested window (date_end lands exactly on
            # this day's start); nothing to fetch or cache for it.
            break
        day_end = day_start + timedelta(days=1)
        fully_spanned = date_start_naive <= day_start and date_end_naive >= day_end
        cache_path = None
        cached = None
        if fully_spanned:
            cache_filename = _get_daily_cache_filename(station_owner, station_id, datum, day)
            cache_path = os.path.join(cache_dir, cache_filename)
            cached = _load_cache(cache_path)
            if cached is not None:
                print(f"Loading {station_owner} data from cache: {cache_path}")
        day_infos.append({
            'start': day_start, 'end': day_end,
            'fully_spanned': fully_spanned,
            'cache_path': cache_path, 'cached': cached,
        })
        day += timedelta(days=1)

    station_name = station_lon = station_lat = None
    time_pieces = []
    wl_pieces = []

    i = 0
    n = len(day_infos)
    while i < n:
        info = day_infos[i]

        if info['cached'] is not None:
            name, lon, lat, day_time, day_wl = info['cached']
            time_pieces.append(day_time)
            wl_pieces.append(day_wl)
            i += 1

        elif not info['fully_spanned']:
            # Boundary day only partially requested: always fetch, never cache.
            frag_start = max(date_start_naive, info['start'])
            frag_end = min(date_end_naive, info['end'])
            name, lon, lat, f_time, f_wl = fetch_range_fn(frag_start, frag_end)
            time.sleep(1)
            if name is not None:
                time_pieces.append(f_time)
                wl_pieces.append(f_wl)
            # else: no data available for this boundary fragment (e.g. it
            # reaches past what's been published yet) -- keep whatever
            # pieces were already gathered instead of discarding them.
            i += 1

        else:
            # Coalesce a contiguous run of uncached, fully-spanned days into one
            # fetch, then split and cache each day's slice individually.
            j = i
            while j < n and day_infos[j]['fully_spanned'] and day_infos[j]['cached'] is None:
                j += 1
            range_start = day_infos[i]['start']
            range_end = day_infos[j - 1]['end']
            name, lon, lat, f_time, f_wl = fetch_range_fn(range_start, range_end)
            time.sleep(1)
            if name is not None:
                f_time = pd.Series(f_time).reset_index(drop=True)
                f_wl = pd.Series(f_wl).reset_index(drop=True)
                for k in range(i, j):
                    k_info = day_infos[k]
                    day_start_utc = pd.Timestamp(k_info['start'], tz='UTC')
                    day_end_utc = pd.Timestamp(k_info['end'], tz='UTC')
                    mask = (f_time >= day_start_utc) & (f_time < day_end_utc)
                    day_time = f_time[mask].reset_index(drop=True)
                    day_wl = f_wl[mask].reset_index(drop=True)
                    print(f"Saving {station_owner} data to cache: {k_info['cache_path']}")
                    _save_cache(k_info['cache_path'], name, lon, lat, day_time, day_wl)
                    time_pieces.append(day_time)
                    wl_pieces.append(day_wl)
            # else: no data available for this whole coalesced batch -- skip
            # caching it (so it's retried on a later run) and keep whatever
            # pieces were already gathered instead of discarding them.
            i = j

        # Prefer the latest chunk's metadata over the earliest: days are
        # processed oldest-first, and the newest (boundary) day is always
        # freshly fetched live rather than read from a possibly much older
        # cache entry, so it reflects the current sensor/name most reliably.
        if name is not None:
            station_name = name
        if lon is not None:
            station_lon = lon
        if lat is not None:
            station_lat = lat

    obs_time = pd.concat(time_pieces, ignore_index=True) if time_pieces else pd.Series([], dtype='datetime64[ns, UTC]')
    obs_wl = pd.concat(wl_pieces, ignore_index=True) if wl_pieces else pd.Series([], dtype=float)
    return station_name, station_lon, station_lat, obs_time, obs_wl

def _fetch_noaa_range(station_id, date_start, date_end, datum, **kwargs):
    """Fetch NOAA water level data for [date_start, date_end]; no caching.

    Returns a None station_name (with empty time/wl series) if no chunk of
    the window yields data -- e.g. a sub-range that reaches past the most
    recently published reading gets a 200 response shaped like
    {"error": {"message": ...}} rather than an outright HTTP failure -- so
    callers can treat it as "no data" rather than a hard error.
    """
    tzutc = pytz.timezone('UTC')

    obs_time = []
    obs_wl = []
    station_name = station_lon = station_lat = None
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
        if response.status_code != 200:
            print(f"Failed to retrieve data: {response.status_code}")
            date_start_i += timedelta(days=31)
            continue
        obs_data = response.json()
        if 'data' not in obs_data:
            # NOAA responds 200 with an {"error": {...}} body (e.g. "No data
            # was found") when this chunk has nothing -- most commonly the
            # tail of a window reaching past the latest published reading.
            message = obs_data.get('error', {}).get('message', 'unknown error')
            print(f"No NOAA data for station {station_id} {date_start_str}-{date_end_str}: {message}")
            date_start_i += timedelta(days=31)
            continue
        # Parse times and make them UTC-aware
        times_parsed = [datetime.strptime(obs_data['data'][i]['t'], '%Y-%m-%d %H:%M') for i in range(len(obs_data['data']))]
        times_utc = [tzutc.localize(t) for t in times_parsed]
        obs_time.extend(times_utc)
        obs_wl.extend([float(obs_data['data'][i]['v']) if obs_data['data'][i]['v'] else np.nan for i in range(len(obs_data['data']))])
        station_name = obs_data['metadata']['name']
        station_lon = float(obs_data['metadata']['lon'])
        station_lat = float(obs_data['metadata']['lat'])
        date_start_i += timedelta(days=31)

    # Convert time list to pandas Series with UTC timezone
    obs_time = pd.Series(obs_time)
    obs_wl = pd.Series(obs_wl)

    return station_name, station_lon, station_lat, obs_time, obs_wl

def _get_noaa_data(station_id, date_start, date_end, datum, **kwargs):
    """Local method to retrieve NOAA water level data (per-day cached)."""
    cache_dir = kwargs.get('cache_dir')
    return _get_data_with_daily_cache(
        'NOAA', station_id, date_start, date_end, datum,
        fetch_range_fn=lambda s, e: _fetch_noaa_range(station_id, s, e, datum, **kwargs),
        cache_dir=cache_dir,
    )

def _fetch_usgs_range(station_id, date_start, date_end, datum, **kwargs):
    """Fetch USGS water level data for [date_start, date_end]; no caching.

    Returns a None station_name (with empty time/wl series) if the interval
    ('iv') service has no data for this window -- e.g. it reaches past what's
    been published yet, which can come back as an empty dataframe or, for a
    narrow enough window, even fail to parse as JSON at all -- so callers can
    treat it as "no data" rather than a hard error. Station metadata ('site'
    service) is not time-windowed, so a failure there still propagates as a
    genuine error (e.g. an invalid station id).
    """
    ft2m = 0.3048

    # Handle timezone-aware datetime objects
    date_start_naive = date_start.replace(tzinfo=None) if date_start.tzinfo else date_start
    date_end_naive = date_end.replace(tzinfo=None) if date_end.tzinfo else date_end
    date_start_str = date_start_naive.strftime('%Y-%m-%dT%H:%M')
    date_end_str = date_end_naive.strftime('%Y-%m-%dT%H:%M')
    dfst = nwis.get_record(sites=station_id, service='site')
    station_name = dfst['station_nm'][0]
    station_lon = dfst['dec_long_va'][0]
    station_lat = dfst['dec_lat_va'][0]

    empty_time = pd.Series([], dtype='datetime64[ns, UTC]')
    empty_wl = pd.Series([], dtype=float)

    try:
        dfiv = nwis.get_record(sites=station_id, service='iv', start=date_start_str, end=date_end_str)
    except Exception as exc:
        print(f"No USGS interval data for station {station_id} {date_start_str} to {date_end_str}: {exc}")
        return None, station_lon, station_lat, empty_time, empty_wl

    if dfiv.empty:
        print(f"No USGS interval data for station {station_id} {date_start_str} to {date_end_str} (empty response)")
        return None, station_lon, station_lat, empty_time, empty_wl

    # Convert time to UTC timezone-aware. pd.to_datetime() on an already
    # DatetimeIndex input (dfiv.index) returns a DatetimeIndex, not a
    # Series -- wrap it so callers can pd.concat() it together with the
    # Series obs_time from other sources/cached days without erroring.
    obs_time = pd.Series(pd.to_datetime(dfiv.index)).reset_index(drop=True)
    if obs_time.dt.tz is None:
        # If timezone-naive, assume UTC
        obs_time = obs_time.dt.tz_localize('UTC')
    else:
        # Convert to UTC if it has a different timezone
        obs_time = obs_time.dt.tz_convert('UTC')

    def _pick_usgs_param_column(df, base_code):
        # Prefer exact legacy names, then tolerate suffixed variants
        # such as "00065_primary sensor" from newer NWIS responses.
        if base_code in df.columns:
            return base_code
        pattern = re.compile(rf"^{re.escape(base_code)}($|[_\\s])")
        for col in df.columns:
            if pattern.match(str(col)):
                return col
        return None

    print('alt_va = ', dfst['alt_va'][0])
    col_00065 = _pick_usgs_param_column(dfiv, '00065')
    col_62620 = _pick_usgs_param_column(dfiv, '62620')
    col_62623 = _pick_usgs_param_column(dfiv, '62623')
    col_00062 = _pick_usgs_param_column(dfiv, '00062')

    if col_00065:
        obs_wl = (dfiv[col_00065] + dfst['alt_va'][0]) * ft2m
    elif col_62620:
        obs_wl = dfiv[col_62620] * ft2m
    elif col_62623:
        obs_wl = dfiv[col_62623] * ft2m
    elif col_00062:
        obs_wl = dfiv[col_00062] * ft2m
    else:
        print(f"No recognized water-level column for station {station_id} in this window. Available columns: {list(dfiv.columns)}")
        return None, station_lon, station_lat, empty_time, empty_wl

    return station_name, station_lon, station_lat, obs_time, obs_wl

def _get_usgs_data(station_id, date_start, date_end, datum, **kwargs):
    """Local method to retrieve USGS water level data (per-day cached)."""
    cache_dir = kwargs.get('cache_dir')
    return _get_data_with_daily_cache(
        'USGS', station_id, date_start, date_end, datum,
        fetch_range_fn=lambda s, e: _fetch_usgs_range(station_id, s, e, datum, **kwargs),
        cache_dir=cache_dir,
    )

def _get_contrail_metadata(station_id, session, username, password):
    """Retrieve CONTRAIL station metadata including device IDs and coordinates"""
    metadata_url = f"https://contrail.nc.gov/site/?site_id={station_id}"
    
    try:
        response = _contrail_fetch_url(session, metadata_url, username, password)
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

def _fetch_contrail_range(station_id, date_start, date_end, datum, **kwargs):
    """Fetch CONTRAIL water level data for [date_start, date_end]; no caching."""
    username = kwargs.get('username')
    password = kwargs.get('password')
    sensor_type = kwargs.get('sensor_type', 'auto')  # Default sensor type

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
    
    session = _create_contrail_session()
    
    # Get station metadata to find device_id and coordinates
    print(f"Retrieving CONTRAIL metadata for station {station_id}...")
    metadata = _get_contrail_metadata(station_id, session, username, password)
    
    # Map sensor type to device_id. 'auto' prefers water_elevation (a true
    # water-surface elevation), then stream_elevation (also elevation-like),
    # and only falls back to stage (often a raw gauge height on an arbitrary
    # local reference, not necessarily comparable to the model's datum) as a
    # last resort.
    available_sensors = list(metadata['sensors'].keys())
    if sensor_type == 'auto':
        for candidate in ('water_elevation', 'stream_elevation', 'stage'):
            if candidate in available_sensors:
                sensor_type = candidate
                break
        else:
            raise ValueError(f"Automatic sensor type determination failed. Available sensors: {available_sensors}")
    if sensor_type not in metadata['sensors']:
        raise ValueError(f"Sensor type '{sensor_type}' not found. Available sensors: {available_sensors}")

    device_id = metadata['sensors'][sensor_type]
    station_name = metadata['station_name']
    if station_name:
        # Surface which sensor was used (esp. relevant when 'auto' resolved
        # it) since it ends up in the plotted hydrograph's title.
        station_name = f"{station_name} [sensor: {sensor_type}]"
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

        return station_name, station_lon, station_lat, obs_time, obs_wl

    except Exception as e:
        print(f"Error retrieving Contrail data: {e}")
        return None, None, None, None, None

def _get_contrail_data(station_id, date_start, date_end, datum, **kwargs):
    """Local method to retrieve Contrail water level data (per-day cached)."""
    username = kwargs.get('username')
    password = kwargs.get('password')
    if not username or not password:
        raise ValueError("Contrail requires 'username' and 'password' in options")

    cache_dir = kwargs.get('cache_dir')
    return _get_data_with_daily_cache(
        'CONTRAIL', station_id, date_start, date_end, datum,
        fetch_range_fn=lambda s, e: _fetch_contrail_range(station_id, s, e, datum, **kwargs),
        cache_dir=cache_dir,
    )

def _fetch_secoora_range(station_id, date_start, date_end, datum, **kwargs):
    """Fetch SECOORA water level data for [date_start, date_end]; no caching."""
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

        return station_name, station_lon, station_lat, obs_time, obs_wl

    except Exception as e:
        print(f"Error processing SECOORA data: {str(e)}")
        return None, None, None, None, None

def _get_secoora_data(station_id, date_start, date_end, datum, **kwargs):
    """Local method to retrieve SECOORA water level data (per-day cached)."""
    if datum != 'NAVD':
        raise ValueError('SECOORA only supports NAVD datum')

    cache_dir = kwargs.get('cache_dir')
    return _get_data_with_daily_cache(
        'SECOORA', station_id, date_start, date_end, datum,
        fetch_range_fn=lambda s, e: _fetch_secoora_range(station_id, s, e, datum, **kwargs),
        cache_dir=cache_dir,
    )

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
        - For CONTRAIL: 'username', 'password', 'sensor_type' (water_elevation, stream_elevation, stage),
          and optionally 'station_id_type' ('auto', 'contrail', 'f61') passed to
          :func:`resolve_contrail_station_ids` before download. Omit or use legacy behavior
          (same id for observation and fort.61) when not set.
          Note: CONTRAIL only supports NAVD88/NAVD datum; other datums will generate a warning
          Station coordinates are automatically extracted from metadata
        - For other sources: additional parameters as needed
    cache_dir : str or Path, optional
        Directory path for caching downloaded observation data in JSON format.
        If provided, one cache file is kept per calendar day (UTC):
        {owner}_{sanitized_station}_{sanitized_datum}_{YYYYMMDD}.json. A day is
        only read from/written to cache when the requested window fully spans
        it; a day only partially covered (typically the first or last day of
        the window) is always fetched fresh and never cached.
    
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

    # Never ask for observations in the future -- these sources have none,
    # and requesting a window that reaches past "now" is what causes the
    # error-shaped ("no data") responses handled deeper in the fetchers.
    date_end = _clip_to_present(date_end)
    if date_start > date_end:
        empty_time = pd.Series([], dtype='datetime64[ns, UTC]')
        empty_wl = pd.Series([], dtype=float)
        return None, None, None, empty_time, empty_wl

    # Dispatch to appropriate local method
    if station_owner == 'NOAA':
        return _get_noaa_data(station_id, date_start, date_end, datum, **options)
    elif station_owner == 'USGS':
        return _get_usgs_data(station_id, date_start, date_end, datum, **options)
    elif station_owner == 'CONTRAIL':
        opts = dict(options)
        station_id_type = opts.pop('station_id_type', None)
        needs_resolve = (
            '/' in str(station_id)
            or (
                station_id_type is not None
                and str(station_id_type).lower() not in ('legacy', '')
            )
        )
        if needs_resolve:
            resolved = resolve_contrail_station_ids(
                station_id,
                station_id_type=station_id_type,
                username=opts.get('username'),
                password=opts.get('password'),
                cache_dir=opts.get('cache_dir'),
            )
            station_id = resolved['contrail_site_id']
        return _get_contrail_data(station_id, date_start, date_end, datum, **opts)
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
        choices=['auto', 'water_elevation', 'stream_elevation', 'stage'],
        default='auto',
        help="Sensor type for CONTRAIL. 'auto' (default) prefers water_elevation, "
             'then stream_elevation, then stage.'
    )
    contrail_group.add_argument(
        '--station-id-type',
        choices=['auto', 'contrail', 'f61'],
        default=None,
        help=(
            'How to interpret station_id for CONTRAIL: '
            'contrail=integer site id, f61=fort.61/ADCIRC code (e.g. EGHN7), '
            'auto=detect (default when omitted: legacy single-id behavior). '
            'Combined ids (e.g. 1205/EGHN7) are always split.'
        ),
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
        if getattr(args, 'station_id_type', None):
            options['station_id_type'] = args.station_id_type
    
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