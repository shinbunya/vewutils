"""Fetch observed significant wave heights from NDBC buoys.

Observations come from NDBC's standard meteorological ("stdmet") text files:

* ``data/realtime2/<id>.txt`` for roughly the last 45 days,
* ``data/stdmet/<Mon>/<id><m><year>.txt.gz`` for earlier months of the
  current year,
* ``data/historical/stdmet/<id>h<year>.txt.gz`` for earlier years.

Only ``requests`` (already a vewutils dependency) and the standard library
are used, so this works anywhere generate-report already runs. Downloads go
through :func:`vewutils.plot.get_obswl._get_data_with_daily_cache`, so they
share the same per-day JSON cache layout as the water level observations.
"""

from __future__ import annotations

import gzip
import json
import os
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import pytz
import requests

from vewutils.plot.get_obswl import (
    _as_utc,
    _atomic_json_dump,
    _clip_to_present,
    _get_data_with_daily_cache,
)

NDBC_BASE_URL = 'https://www.ndbc.noaa.gov'
NDBC_ACTIVE_STATIONS_URL = f'{NDBC_BASE_URL}/activestations.xml'
NDBC_REALTIME_URL = NDBC_BASE_URL + '/data/realtime2/{station}.txt'
NDBC_MONTHLY_URL = NDBC_BASE_URL + '/data/stdmet/{mon}/{station}{mcode}{year}.txt.gz'
NDBC_HISTORICAL_URL = NDBC_BASE_URL + '/data/historical/stdmet/{station}h{year}.txt.gz'

NDBC_STATION_LIST_CACHE = 'ndbc_active_stations.json'
NDBC_STATION_LIST_MAX_AGE = timedelta(days=7)
# realtime2 files hold the last 45 days; stay a little inside that.
NDBC_REALTIME_SPAN = timedelta(days=44)
NDBC_REQUEST_TIMEOUT_SEC = 60
NDBC_REQUEST_PAUSE_SEC = 1
# Cache key "datum" slot for wave heights (there's no datum for Hs).
NDBC_HS_CACHE_TAG = 'WVHT'

# activestations.xml station types. DART (tsunami) and TAO (equatorial
# Pacific) stations don't report waves.
DEFAULT_NDBC_STATION_TYPES = ('buoy', 'fixed', 'other', 'oilrig')
# NOS/CO-OPS water level stations are listed too but report no waves.
DEFAULT_NDBC_EXCLUDE_PROGRAMS = ('NOS/CO-OPS',)

_MONTH_ABBR = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
_MONTH_CODE = '123456789abc'

# url -> parsed DataFrame (or None for 404), reused within one process so a
# report touching the same realtime2 file for several sub-ranges downloads it
# once.
_download_memory: dict[str, pd.DataFrame | None] = {}


class NDBCRequestError(RuntimeError):
    """An NDBC request failed for a reason other than "no such file"."""


def _http_get(url: str) -> requests.Response | None:
    """GET url; None for a 404, NDBCRequestError for any other failure."""
    print(url)
    try:
        response = requests.get(url, timeout=NDBC_REQUEST_TIMEOUT_SEC)
    except requests.RequestException as exc:
        raise NDBCRequestError(f'{url}: {exc}') from exc
    finally:
        time.sleep(NDBC_REQUEST_PAUSE_SEC)
    if response.status_code == 404:
        return None
    if response.status_code != 200:
        raise NDBCRequestError(f'{url}: HTTP {response.status_code}')
    return response


# ---------------------------------------------------------------------------
# Station list
# ---------------------------------------------------------------------------

def parse_ndbc_active_stations_xml(xml_text: str | bytes) -> list[dict]:
    """Parse NDBC activestations.xml into a list of station dicts.

    Keys: ``station_id``, ``lon``, ``lat``, ``name``, ``owner``, ``program``,
    ``type``, ``met``.
    """
    root = ET.fromstring(xml_text)
    stations = []
    for node in root.iter('station'):
        try:
            lat = float(node.get('lat'))
            lon = float(node.get('lon'))
        except (TypeError, ValueError):
            continue
        stations.append({
            'station_id': (node.get('id') or '').strip().upper(),
            'lon': lon,
            'lat': lat,
            'name': (node.get('name') or '').strip(),
            'owner': (node.get('owner') or '').strip(),
            'program': (node.get('pgm') or '').strip(),
            'type': (node.get('type') or '').strip().lower(),
            'met': (node.get('met') or '').strip().lower(),
        })
    return stations


def get_ndbc_active_stations(
        cache_dir=None, force_refresh=False, cache_only=False) -> list[dict]:
    """Return NDBC's active station list, cached for up to a week in cache_dir.

    ``cache_only`` returns the cached list whatever its age and never
    downloads; a missing list raises ValueError.
    """
    cache_path = (
        os.path.join(cache_dir, NDBC_STATION_LIST_CACHE) if cache_dir else None
    )
    if cache_path and not force_refresh and os.path.isfile(cache_path):
        try:
            with open(cache_path, encoding='utf-8') as f:
                cached = json.load(f)
            fetched_at = datetime.fromisoformat(cached['fetched_at'])
            if (
                cache_only
                or datetime.now(pytz.UTC) - fetched_at < NDBC_STATION_LIST_MAX_AGE
            ):
                return cached['stations']
        except (OSError, ValueError, KeyError, TypeError):
            pass
    if cache_only:
        raise ValueError(
            'NDBC station list is not in the cache. '
            'Run `vewutils plot refresh-obs-cache` before a cache-only report.'
        )

    response = _http_get(NDBC_ACTIVE_STATIONS_URL)
    if response is None:
        raise NDBCRequestError(f'{NDBC_ACTIVE_STATIONS_URL}: HTTP 404')
    stations = parse_ndbc_active_stations_xml(response.content)

    if cache_path:
        _atomic_json_dump(cache_path, {
            'fetched_at': datetime.now(pytz.UTC).isoformat(),
            'stations': stations,
        })
    return stations


def filter_ndbc_stations(
        stations: list[dict],
        *,
        bbox: tuple[float, float, float, float] | None = None,
        station_ids=None,
        station_types=DEFAULT_NDBC_STATION_TYPES,
        exclude_programs=DEFAULT_NDBC_EXCLUDE_PROGRAMS) -> list[dict]:
    """Filter an NDBC station list by bounding box, ids, type and program.

    bbox is ``(xmin, xmax, ymin, ymax)`` in degrees. station_ids, when given,
    overrides the type/program filters (an explicitly requested id is kept
    whatever its type) and sets the output order.
    """
    def in_bbox(st):
        if bbox is None:
            return True
        xmin, xmax, ymin, ymax = bbox
        return xmin <= st['lon'] <= xmax and ymin <= st['lat'] <= ymax

    if station_ids:
        by_id = {st['station_id']: st for st in stations}
        selected = []
        for sid in station_ids:
            st = by_id.get(str(sid).strip().upper())
            if st is None:
                print(
                    f'Warning: NDBC station {sid} not in the active station list',
                    file=sys.stderr,
                )
            elif in_bbox(st):
                selected.append(st)
        return selected

    types = {t.lower() for t in station_types} if station_types else None
    excluded = set(exclude_programs or ())
    return [
        st for st in stations
        if in_bbox(st)
        and (types is None or st['type'] in types)
        and st['program'] not in excluded
    ]


# ---------------------------------------------------------------------------
# stdmet files
# ---------------------------------------------------------------------------

def parse_ndbc_stdmet_text(text: str) -> pd.DataFrame:
    """Parse an NDBC stdmet text file into a DataFrame with ``time`` and ``wvht``.

    Handles realtime2 files (``#YY MM DD hh mm ...`` with ``MM`` for missing)
    and historical files (``YY``/``YYYY`` header, optional minute column,
    ``99.00`` for missing WVHT). ``time`` is UTC-aware; missing heights are
    dropped.
    """
    header = None
    rows = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith('#') or not stripped[0].isdigit():
            tokens = stripped.lstrip('#').split()
            # First comment/header line names the columns; the second holds units.
            if header is None and tokens and tokens[0].upper() in ('YY', 'YYYY'):
                header = [t.upper() for t in tokens]
            continue
        rows.append(stripped.split())

    empty = pd.DataFrame({
        'time': pd.Series([], dtype='datetime64[ns, UTC]'),
        'wvht': pd.Series([], dtype=float),
    })
    if header is None or 'WVHT' not in header or not rows:
        return empty

    ncol = len(header)
    rows = [r[:ncol] for r in rows if len(r) >= ncol]
    if not rows:
        return empty
    # Columns are addressed by position: "MM" appears twice (month and minute).
    table = np.array(rows, dtype=object)

    def column(index):
        return pd.to_numeric(pd.Series(table[:, index]), errors='coerce')

    year = column(0)
    year = year.where(year >= 100, year + 1900)  # pre-1999 two-digit years
    has_minute = len(header) > 4 and header[4] == 'MM'
    parts = pd.DataFrame({
        'year': year,
        'month': column(1),
        'day': column(2),
        'hour': column(3),
        'minute': column(4) if has_minute else 0,
    })
    valid_time = parts.notna().all(axis=1)
    times = pd.Series(pd.NaT, index=parts.index, dtype='datetime64[ns, UTC]')
    times[valid_time] = pd.to_datetime(
        parts[valid_time].astype(int), errors='coerce', utc=True
    )
    wvht = column(header.index('WVHT'))
    wvht = wvht.where(wvht < 99.0)

    out = pd.DataFrame({'time': times, 'wvht': wvht.astype(float)})
    out = out.dropna().sort_values('time').drop_duplicates('time', keep='last')
    return out.reset_index(drop=True)


def _download_stdmet(url: str, gzipped: bool) -> pd.DataFrame | None:
    if url in _download_memory:
        return _download_memory[url]
    response = _http_get(url)
    if response is None:
        df = None
    else:
        content = response.content
        if gzipped:
            try:
                content = gzip.decompress(content)
            except OSError:
                pass  # server already decompressed it
        df = parse_ndbc_stdmet_text(content.decode('utf-8', errors='replace'))
    _download_memory[url] = df
    return df


def _fetch_ndbc_hs_range(station_id, date_start, date_end, *, station=None, now=None):
    """Fetch NDBC WVHT for [date_start, date_end]; no on-disk caching.

    Returns ``(name, lon, lat, obs_time, obs_hs)`` like the water level
    fetchers in :mod:`vewutils.plot.get_obswl`.
    """
    now = now or datetime.now(pytz.UTC)
    date_start = _as_utc(date_start)
    date_end = _as_utc(date_end)
    sid = str(station_id).lower()
    realtime_start = now - NDBC_REALTIME_SPAN

    pieces: list[pd.DataFrame] = []
    covered_until = date_start
    if date_start < realtime_start:
        # Older than realtime2: one file per past year, one per month this year.
        months = pd.period_range(
            date_start.replace(tzinfo=None),
            min(date_end, realtime_start).replace(tzinfo=None),
            freq='M',
        )
        done_years = set()
        for period in months:
            if period.year < now.year:
                if period.year in done_years:
                    continue
                done_years.add(period.year)
                url = NDBC_HISTORICAL_URL.format(station=sid, year=period.year)
            else:
                url = NDBC_MONTHLY_URL.format(
                    mon=_MONTH_ABBR[period.month - 1],
                    station=sid,
                    mcode=_MONTH_CODE[period.month - 1],
                    year=period.year,
                )
            df = _download_stdmet(url, gzipped=True)
            if df is not None and len(df):
                pieces.append(df)
        if pieces:
            covered_until = max(df['time'].max() for df in pieces)

    if date_end >= realtime_start or covered_until < date_end:
        df = _download_stdmet(
            NDBC_REALTIME_URL.format(station=str(station_id).upper()), gzipped=False
        )
        if df is not None and len(df):
            pieces.append(df)

    if pieces:
        data = pd.concat(pieces, ignore_index=True)
        data = data.sort_values('time').drop_duplicates('time', keep='last')
        data = data[(data['time'] >= date_start) & (data['time'] <= date_end)]
    else:
        data = parse_ndbc_stdmet_text('')

    station = station or {}
    name = station.get('name') or str(station_id).upper()
    return (
        name,
        station.get('lon'),
        station.get('lat'),
        data['time'].reset_index(drop=True),
        data['wvht'].reset_index(drop=True),
    )


def get_obswave(
        station_id, date_start, date_end, *, station=None, cache_dir=None,
        cache_only=False, cache_open_day=False):
    """Observed significant wave height (m) at an NDBC station.

    Parameters
    ----------
    station_id : str
        NDBC station id (e.g. ``41025``).
    date_start, date_end : datetime
        Window to fetch (naive values are UTC).
    station : dict, optional
        Station metadata (``name``, ``lon``, ``lat``) to return alongside the
        data, e.g. one entry of :func:`get_ndbc_active_stations`.
    cache_dir : str or Path, optional
        Per-day JSON cache directory shared with the water level observations
        (files are named ``NDBC_<id>_wvht_<YYYYMMDD>.json``).
    cache_only, cache_open_day : bool, optional
        As in :func:`vewutils.plot.get_obswl.get_obswl`: read only the cache,
        or also store the current incomplete UTC day (marked partial).

    Returns
    -------
    tuple
        ``(station_name, station_lon, station_lat, obs_time, obs_hs)`` with
        ``obs_time`` a UTC-aware pandas Series.
    """
    date_end = _clip_to_present(_as_utc(date_end))
    date_start = _as_utc(date_start)
    if date_start > date_end:
        return None, None, None, pd.Series([], dtype='datetime64[ns, UTC]'), pd.Series([], dtype=float)

    return _get_data_with_daily_cache(
        'NDBC',
        str(station_id).upper(),
        date_start,
        date_end,
        NDBC_HS_CACHE_TAG,
        fetch_range_fn=lambda s, e: _fetch_ndbc_hs_range(station_id, s, e, station=station),
        cache_dir=str(cache_dir) if cache_dir else None,
        # _http_get already pauses after every request.
        pause_after_fetch=False,
        cache_only=cache_only,
        cache_open_day=cache_open_day,
    )


def obs_series_to_numpy(obs_time, obs_hs):
    """UTC-aware Series -> (naive datetime64 array, float array), finite only."""
    if obs_time is None or obs_hs is None or len(obs_time) == 0:
        return np.array([], dtype='datetime64[ns]'), np.array([], dtype=float)
    times = pd.to_datetime(pd.Series(obs_time), utc=True).dt.tz_localize(None).to_numpy()
    values = pd.to_numeric(pd.Series(obs_hs), errors='coerce').to_numpy(dtype=float)
    valid = np.isfinite(values)
    return times[valid], values[valid]
