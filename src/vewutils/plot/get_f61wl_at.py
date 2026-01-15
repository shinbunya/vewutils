# %%
import os
import pickle
import re

import numpy as np
import xarray as xr
import pandas as pd


# %%
class StationNotFoundError(Exception):
    """Raised when a station is not found in fort.61.nc file"""
    pass


# %%
def get_f61wl_at(f61file, station_name):
    """
    Extract water level time series at a given station_name from an ADCIRC f61
    NetCDF file. The zeta variable is already indexed by station, so no
    interpolation is needed.
    Results are cached on disk for subsequent calls.
    
    Raises
    ------
    StationNotFoundError
        If the specified station_name is not found in the file or an error occurs
        while reading the file.
    """
    # Sanitize station name for cache filename (uppercase, replace special chars)
    sanitized_name = re.sub(r'[^a-zA-Z0-9_]', '_', station_name.strip().upper())
    cachefile = f61file.replace(
        ".nc",
        f"_{sanitized_name}_cache.pkl",
    )

    # ------------------------------------------------------------------
    # 1. Cache short‑circuit
    # ------------------------------------------------------------------
    if os.path.exists(cachefile):
        with open(cachefile, "rb") as f:
            f61_time, f61_wl = pickle.load(f)
        return f61_time, f61_wl

    # ------------------------------------------------------------------
    # 2. Read station names and find matching station
    # ------------------------------------------------------------------
    try:
        with xr.open_dataset(f61file) as ds:
            # Read station names - they are stored as fixed-length byte strings (S50)
            station_names_raw = ds["station_name"].values  # numpy.bytes_ objects
            station_names = []
            for i in range(len(station_names_raw)):
                # Each element is a numpy.bytes_ object, decode it directly
                name = station_names_raw[i].decode("utf-8", errors='ignore').strip().upper()
                station_names.append(name)
            
            # Time axis
            f61_time = pd.to_datetime(
                ds["time"].values.astype("datetime64[ms]")
            ).to_pydatetime()
            nt = len(f61_time)
    except Exception as e:
        print(f"Error reading fort.61.nc file '{f61file}': {str(e)}")
        raise StationNotFoundError(
            f"Error reading fort.61.nc file '{f61file}': {str(e)}"
        ) from e

    # ------------------------------------------------------------------
    # 3. Match station name (case-insensitive, trimmed)
    # ------------------------------------------------------------------
    search_name = station_name.strip().upper()
    station_idx = None
    
    for i, name in enumerate(station_names):
        if name == search_name:
            station_idx = i
            break

    # ------------------------------------------------------------------
    # 4. Station not found - raise exception
    # ------------------------------------------------------------------
    if station_idx is None:
        # Check for partial matches to help diagnose
        partial_matches = [s for s in station_names if search_name in s or s in search_name]
        available_stations = station_names[:10] if len(station_names) > 10 else station_names
        error_msg = f"Station '{station_name}' not found in file '{f61file}'. "
        if partial_matches:
            error_msg += f"Partial matches found: {partial_matches[:5]}. "
        error_msg += f"Available stations (first {len(available_stations)}): {available_stations}"
        raise StationNotFoundError(error_msg)

    # ------------------------------------------------------------------
    # 5. Read zeta data for the matched station using xarray
    # ------------------------------------------------------------------
    try:
        with xr.open_dataset(f61file) as ds:
            # Read zeta for the specific station (time, station) -> (time,)
            f61_wl = ds["zeta"][:, station_idx].values
            print(f"f61_wl.min(), f61_wl.max() = {f61_wl.min(), f61_wl.max()}")
            
            # Handle fill values
            fill_value = ds["zeta"].attrs.get("_FillValue", -99999.0)
            f61_wl = np.where(f61_wl == fill_value, np.nan, f61_wl)
    except Exception as e:
        raise StationNotFoundError(
            f"Error reading zeta data for station '{station_name}' from file '{f61file}': {str(e)}"
        ) from e

    # ------------------------------------------------------------------
    # 7. Cache results
    # ------------------------------------------------------------------
    with open(cachefile, "wb") as f:
        pickle.dump([f61_time, f61_wl], f)

    return f61_time, f61_wl

