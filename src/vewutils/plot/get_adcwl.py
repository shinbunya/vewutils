# %%
import os

from vewutils.plot.get_f61wl_at import get_f61wl_at, StationNotFoundError
from vewutils.plot.get_f63wl_at import get_f63wl_at


# %%
def get_adcwl(file, station_name=None, stx=None, sty=None, fallback_file=None):
    """
    Wrapper function to extract water level time series from ADCIRC NetCDF files.
    Auto-detects file type based on filename pattern and routes to the appropriate
    extraction function.
    
    Parameters
    ----------
    file : str
        Path to ADCIRC NetCDF file (fort.61.nc or fort.63.nc)
    station_name : str, optional
        Station name/ID for fort.61.nc files. Required for fort.61.nc files.
        Station names are matched case-insensitively with trimmed whitespace.
    stx : float, optional
        Longitude for fort.63.nc files or for fort.61.nc files with full mesh data.
    sty : float, optional
        Latitude for fort.63.nc files or for fort.61.nc files with full mesh data.
    fallback_file : str, optional
        Fallback file path (typically fort.63.nc) to use if station_name is not found
        in fort.61.nc file. Requires stx and sty to be specified.
    
    Returns
    -------
    time : array-like
        Array of datetime objects representing model time
    wl : array-like
        Array of water level values (meters)
    
    Raises
    ------
    StationNotFoundError
        If station_name is not found in fort.61.nc file and no valid fallback is available
    
    Examples
    --------
    >>> # For fort.61.nc files (station-based)
    >>> time, wl = get_adcwl('fort.61.nc', station_name='CTIN7')
    >>> 
    >>> # For fort.63.nc files (coordinate-based)
    >>> time, wl = get_adcwl('fort.63.nc', stx=-76.5, sty=35.2)
    >>> 
    >>> # For fort.61.nc with fallback to fort.63.nc
    >>> time, wl = get_adcwl('fort.61.nc', station_name='8651370', 
    ...                       stx=-75.74, sty=36.18, fallback_file='fort.63.nc')
    """
    # Determine file type from filename
    filename = os.path.basename(file)
    
    if 'fort.61' in filename or 'f61' in filename.lower() or 'fort61' in filename.lower():
        # fort.61.nc file - requires station_name (or coordinates if it has mesh data)
        if station_name is not None:
            # Try station-based lookup first
            try:
                return get_f61wl_at(file, station_name)
            except StationNotFoundError:
                # If station not found and fallback is available, use it
                if fallback_file is not None and stx is not None and sty is not None:
                    print(f'\nWarning: Station {station_name} not found in {file}')
                    print(f'Falling back to {fallback_file} with coordinates ({stx}, {sty})')
                    return get_f63wl_at(fallback_file, stx, sty)
                else:
                    # No fallback available, re-raise the exception
                    raise
        elif stx is not None and sty is not None:
            # Try coordinate-based lookup (works if file has full mesh data)
            return get_f63wl_at(file, stx, sty)
        else:
            raise ValueError(
                f"Either station_name or (stx, sty) coordinates are required for fort.61.nc files. "
                f"Got: file={file}, station_name={station_name}, stx={stx}, sty={sty}"
            )
    
    elif 'fort.63' in filename or 'f63' in filename.lower() or 'fort63' in filename.lower():
        # fort.63.nc file - requires coordinates
        if stx is None or sty is None:
            raise ValueError(
                f"stx and sty coordinates are required for fort.63.nc files. "
                f"Got: file={file}, stx={stx}, sty={sty}"
            )
        return get_f63wl_at(file, stx, sty)
    
    else:
        raise ValueError(
            f"Unrecognized file type. Expected fort.61.nc or fort.63.nc. "
            f"Got: {filename}"
        )

