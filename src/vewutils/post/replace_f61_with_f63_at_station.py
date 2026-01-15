#!/usr/bin/env python3
"""
Module for replacing incorrect water level data for a station in fort.61.nc files
with data extracted from fort.63.nc at corrected coordinates.

Usage:
    python replace_f61_with_f63_at_station.py file1.fort.61.nc file2.fort.61.nc ... STATION_NAME x y

The program automatically finds corresponding fort.63.nc files in the same directory
as each fort.61.nc file and generates output files with '_corrected' suffix.
"""

import argparse
import numpy as np
import netCDF4 as nc
from netCDF4 import num2date
from os import path
from typing import Tuple
import xarray as xr
import pandas as pd
import re


def decode_station_names(station_name_array: np.ndarray) -> list:
    """
    Decode 2D char array to list of station name strings.
    
    Args:
        station_name_array: numpy array with shape (num_stations, namelen)
    
    Returns:
        List of station names as strings
    """
    if hasattr(station_name_array, 'tobytes'):
        # netCDF4 chartostring is the standard way
        try:
            return [name.strip() for name in nc.chartostring(station_name_array)]
        except:
            pass
    
    # Fallback: manual decoding
    station_names = []
    for i in range(station_name_array.shape[0]):
        # Convert char array to bytes, decode to string, strip whitespace/nulls
        name_bytes = station_name_array[i, :].tobytes()
        name = name_bytes.decode('utf-8', errors='ignore').strip('\x00').strip()
        station_names.append(name)
    
    return station_names


def find_station_index(f61file: str, station_name: str) -> int:
    """
    Find the index of a station by name in fort.61.nc file.
    
    Args:
        f61file: Path to fort.61.nc file
        station_name: Name of the station to find (case-insensitive)
    
    Returns:
        Station index
    
    Raises:
        ValueError: If station is not found
    """
    with nc.Dataset(f61file) as src:
        station_names = decode_station_names(src['station_name'][:])
    
    search_name = station_name.strip().upper()
    
    for i, name in enumerate(station_names):
        if name.upper() == search_name:
            return i
    
    # Station not found - provide helpful error message
    partial_matches = [s for s in station_names if search_name in s.upper() or s.upper() in search_name]
    available_stations = station_names[:10] if len(station_names) > 10 else station_names
    error_msg = f"Station '{station_name}' not found in file '{f61file}'. "
    if partial_matches:
        error_msg += f"Partial matches found: {partial_matches[:5]}. "
    error_msg += f"Available stations (first {len(available_stations)}): {available_stations}"
    raise ValueError(error_msg)


def extract_f63_water_level(f63file: str, stx: float, sty: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract water level time series from fort.63.nc at given coordinates.
    
    Args:
        f63file: Path to fort.63.nc file
        stx: X coordinate
        sty: Y coordinate
    
    Returns:
        Tuple of (time_array, water_level_array)
        If point is outside domain, water_level_array will be all NaN
    """
    # Helper: point-in-triangle test (vectorized over many elements)
    def is_point_in_triangle_multi(px, py, x1, y1, x2, y2, x3, y3):
        denominator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3))
        a = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / denominator
        b = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / denominator
        c = 1.0 - a - b
        return np.all(
            [
                np.all([0.0 <= a, a <= 1.0], axis=0),
                np.all([0.0 <= b, b <= 1.0], axis=0),
                np.all([0.0 <= c, c <= 1.0], axis=0),
            ],
            axis=0,
        )
    
    # Use xarray for mesh + time (cheap compared to zeta)
    with xr.open_dataset(f63file) as ds:
        # Mesh
        adc_x = ds["x"].values
        adc_y = ds["y"].values
        adc_e = ds["element"].values - 1  # 0-based
        
        # Find containing element (this can be slow for large meshes)
        num_elements = len(adc_e)
        if num_elements > 100000:
            # For very large meshes, this operation can take time
            pass  # Could add progress feedback here if needed
        isinside = is_point_in_triangle_multi(
            stx,
            sty,
            adc_x[adc_e[:, 0]],
            adc_y[adc_e[:, 0]],
            adc_x[adc_e[:, 1]],
            adc_y[adc_e[:, 1]],
            adc_x[adc_e[:, 2]],
            adc_y[adc_e[:, 2]],
        )
        
        # Time axis - check if it has units attribute for proper conversion
        time_var = ds["time"]
        time_raw = time_var.values
        
        if hasattr(time_var, 'units'):
            # Use netCDF4's num2date for proper time conversion
            f63_time_cf = num2date(time_raw, units=time_var.units, calendar=getattr(time_var, 'calendar', 'standard'))
            # Convert cftime objects to pandas datetime
            try:
                # Try to convert directly - works for Python datetime objects
                f63_time = pd.to_datetime(f63_time_cf)
            except (TypeError, ValueError):
                # If that fails, it's likely cftime objects - convert via datetime property
                try:
                    f63_time = pd.to_datetime([dt.datetime for dt in f63_time_cf])
                except (AttributeError, TypeError):
                    # Last resort: convert to string then to datetime
                    f63_time = pd.to_datetime([str(dt) for dt in f63_time_cf])
        else:
            # Fallback to datetime64 conversion
            print(f"  WARNING: No time units attribute found, using datetime64 conversion")
            f63_time = pd.to_datetime(time_raw.astype("datetime64[ms]")).to_pydatetime()
        
        nt = len(f63_time)
    
    # Point not inside any triangle
    if not np.any(isinside):
        print(f"ERROR: Point ({stx}, {sty}) is not inside any triangle in {f63file}")
        if stx < np.min(adc_x) or stx > np.max(adc_x) or sty < np.min(adc_y) or sty > np.max(adc_y):
            print(f"ERROR:   Point is outside mesh bounding box!")
            print(f"ERROR:   Mesh extent: X=[{np.min(adc_x):.6f}, {np.max(adc_x):.6f}], Y=[{np.min(adc_y):.6f}, {np.max(adc_y):.6f}]")
            print(f"ERROR:   This suggests a coordinate system mismatch.")
            print(f"ERROR:   The point coordinates ({stx:.6f}, {sty:.6f}) may be in a different")
            print(f"ERROR:   coordinate system (e.g., WGS84 longitude/latitude) than the mesh")
            print(f"ERROR:   coordinates (e.g., projected coordinates).")
            print(f"ERROR:   Please verify the coordinate system of both the point and the mesh.")
        else:
            print(f"ERROR:   Point is within mesh bounding box but not inside any element")
            print(f"ERROR:   This may indicate a mesh issue or coordinate system mismatch")
        f63_wl = np.full(nt, np.nan)
        return f63_time, f63_wl
    
    # Element geometry + barycentric weights (only 3 nodes)
    elemn = adc_e[isinside][0]
    elemx = adc_x[elemn]
    elemy = adc_y[elemn]
    
    x1, y1 = elemx[0], elemy[0]
    x2, y2 = elemx[1], elemy[1]
    x3, y3 = elemx[2], elemy[2]
    
    denominator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3))
    a = ((y2 - y3) * (stx - x3) + (x3 - x2) * (sty - y3)) / denominator
    b = ((y3 - y1) * (stx - x3) + (x1 - x3) * (sty - y3)) / denominator
    c = 1.0 - a - b
    
    if abs(a + b + c - 1.0) > 1e-6:
        print(f"WARNING: Barycentric weights don't sum to 1.0 (sum={a + b + c:.6f})")
    
    # Heavy part: read zeta with netCDF4 directly (time, node)
    # This can be slow for large files, so ensure proper resource management
    with nc.Dataset(f63file, mode="r") as nc_file:
        zeta_var = nc_file.variables["zeta"]  # dims typically ('time', 'node')
        
        # Detect chunking to choose optimal read strategy
        chunking = zeta_var.chunking()
        if chunking and len(chunking) == 2:
            time_chunk, node_chunk = chunking
            # If chunked by node (node_chunk is small, e.g., 1), fancy indexing is faster
            # If chunked by time (time_chunk is small, e.g., 1), individual reads are faster
            is_node_chunked = node_chunk <= 10 and time_chunk > 100
        else:
            # Default: assume bad chunking (chunked by time)
            is_node_chunked = False
        
        # Check for fill values in zeta
        fill_value = zeta_var._FillValue if hasattr(zeta_var, '_FillValue') else None
        
        if is_node_chunked:
            # Optimized chunking: read all 3 nodes at once (faster for node-chunked files)
            elemv_all = zeta_var[:, elemn]
        else:
            # Bad chunking: read each node separately (faster for time-chunked files)
            # When file is chunked [1, 506054], reading nodes individually is much faster
            # than fancy indexing because it benefits from internal caching
            node1_data = zeta_var[:, elemn[0]]
            node2_data = zeta_var[:, elemn[1]]
            node3_data = zeta_var[:, elemn[2]]
            elemv_all = np.column_stack([node1_data, node2_data, node3_data])
        
        # Convert fill values to NaN before interpolation
        if fill_value is not None:
            fill_count_total = np.sum(elemv_all == fill_value)
            if fill_count_total > 0:
                print(f"  Converting {fill_count_total} fill values to NaN...")
                elemv_all = np.where(elemv_all == fill_value, np.nan, elemv_all)
    
    # Vectorized interpolation over time
    f63_wl = a * elemv_all[:, 0] + b * elemv_all[:, 1] + c * elemv_all[:, 2]
    
    if np.any(np.isnan(f63_wl)):
        nan_count = np.sum(np.isnan(f63_wl))
        print(f"  WARNING: {nan_count} NaN values in interpolated water level data")
    
    return f63_time, f63_wl


def find_corresponding_f63_file(f61file: str) -> str:
    """
    Find the corresponding fort.63.nc file in the same directory as the fort.61.nc file.
    
    Args:
        f61file: Path to fort.61.nc file
    
    Returns:
        Path to corresponding fort.63.nc file
    
    Raises:
        FileNotFoundError: If corresponding fort.63.nc file is not found
    """
    f61_dir = path.dirname(path.abspath(f61file))
    f61_basename = path.basename(f61file)
    
    # Try various patterns to find fort.63.nc file
    # Pattern 1: Replace "61" with "63" in filename
    f63_basename1 = f61_basename.replace("61", "63")
    f63_file1 = path.join(f61_dir, f63_basename1)
    if path.exists(f63_file1):
        return f63_file1
    
    # Pattern 2: Replace "f61" with "f63" (case insensitive)
    f63_basename2 = re.sub(r'f61', 'f63', f61_basename, flags=re.IGNORECASE)
    if f63_basename2 != f61_basename:
        f63_file2 = path.join(f61_dir, f63_basename2)
        if path.exists(f63_file2):
            return f63_file2
    
    # Pattern 3: Replace "fort.61" with "fort.63"
    f63_basename3 = f61_basename.replace("fort.61", "fort.63")
    if f63_basename3 != f61_basename:
        f63_file3 = path.join(f61_dir, f63_basename3)
        if path.exists(f63_file3):
            return f63_file3
    
    # Pattern 4: Look for "fort.63.nc" in the same directory
    f63_file4 = path.join(f61_dir, "fort.63.nc")
    if path.exists(f63_file4):
        return f63_file4
    
    # If none found, raise error with list of attempted files
    tried_files = [f63_file1]
    if f63_basename2 != f61_basename:
        tried_files.append(f63_file2)
    if f63_basename3 != f61_basename:
        tried_files.append(f63_file3)
    tried_files.append(f63_file4)
    
    raise FileNotFoundError(
        f"Could not find corresponding fort.63.nc file for {f61file}. "
        f"Tried: {', '.join(tried_files)}"
    )


def generate_output_filename(f61file: str) -> str:
    """
    Generate output filename by adding '_corrected' suffix before the extension.
    
    Args:
        f61file: Path to input fort.61.nc file
    
    Returns:
        Path to output file with '_corrected' suffix
    """
    f61_dir = path.dirname(f61file)
    f61_basename = path.basename(f61file)
    
    # Split filename and extension
    if f61_basename.endswith('.nc'):
        name_without_ext = f61_basename[:-3]
        output_basename = f"{name_without_ext}_corrected.nc"
    else:
        # If no .nc extension, just append _corrected
        output_basename = f"{f61_basename}_corrected"
    
    return path.join(f61_dir, output_basename)


def interpolate_to_f61_time(
    f61_time: np.ndarray,
    f63_time: np.ndarray,
    f63_wl: np.ndarray
) -> np.ndarray:
    """
    Interpolate fort.63.nc water level data to match fort.61.nc time axis.
    
    Args:
        f61_time: Time array from fort.61.nc (target time axis)
        f63_time: Time array from fort.63.nc (source time axis)
        f63_wl: Water level array from fort.63.nc
    
    Returns:
        Interpolated water level array matching f61_time length
        Values outside f63_time range will be NaN
    """
    if np.any(np.isnan(f63_wl)):
        nan_count = np.sum(np.isnan(f63_wl))
        print(f"  WARNING: Input f63_wl contains {nan_count} NaN values")
    
    # Convert to pandas datetime (remove timezone if present)
    f61_time_pd = pd.to_datetime(f61_time)
    if hasattr(f61_time_pd, 'tz') and f61_time_pd.tz is not None:
        f61_time_pd = f61_time_pd.tz_localize(None)
    
    f63_time_pd = pd.to_datetime(f63_time)
    if hasattr(f63_time_pd, 'tz') and f63_time_pd.tz is not None:
        f63_time_pd = f63_time_pd.tz_localize(None)
    
    # Create DataFrame with f63 data
    f63_df = pd.DataFrame({
        'time': f63_time_pd,
        'wl': f63_wl
    })
    f63_df.set_index('time', inplace=True)
    
    # Reindex to f61_time (this will create NaN for times outside f63 range)
    f63_df_reindexed = f63_df.reindex(f61_time_pd)
    
    # Check overlap
    overlap_mask = (f61_time_pd >= f63_time_pd.min()) & (f61_time_pd <= f63_time_pd.max())
    num_overlap = np.sum(overlap_mask)
    
    # Interpolate linearly (NaN values outside range will remain NaN - no extrapolation)
    f63_df_interp = f63_df_reindexed.interpolate(method='linear')
    
    # Extract interpolated values
    f63_wl_interp = f63_df_interp['wl'].values
    
    if np.any(np.isnan(f63_wl_interp)):
        nan_count = np.sum(np.isnan(f63_wl_interp))
        print(f"  WARNING: {nan_count} NaN values in interpolated data")
    
    # Count how many values are outside the f63 time range
    f63_time_min = f63_time_pd.min()
    f63_time_max = f63_time_pd.max()
    outside_range = (f61_time_pd < f63_time_min) | (f61_time_pd > f63_time_max)
    num_outside = np.sum(outside_range)
    
    if num_outside > 0:
        print(f"  Warning: {num_outside} time steps in fort.61.nc are outside fort.63.nc time range")
        print(f"    fort.61.nc range: {f61_time_pd.min()} to {f61_time_pd.max()}")
        print(f"    fort.63.nc range: {f63_time_min} to {f63_time_max}")
        print(f"    These will be set to NaN (converted to fill value)")
    
    return f63_wl_interp


def replace_f61_station_with_f63(
    f61file: str,
    f63file: str,
    station_name: str,
    correct_x: float,
    correct_y: float,
    output_file: str
) -> None:
    """
    Replace water level data for a station in fort.61.nc with data from fort.63.nc.
    
    Args:
        f61file: Path to input fort.61.nc file (with wrong data)
        f63file: Path to fort.63.nc file (source of correct data)
        station_name: Name of the station to correct
        correct_x: Correct x coordinate
        correct_y: Correct y coordinate
        output_file: Path to output fort.61.nc file
    
    Raises:
        ValueError: If input and output files are the same, or other validation errors
        FileNotFoundError: If input files don't exist
    """
    # Validate input and output files are different
    if path.abspath(f61file) == path.abspath(output_file):
        raise ValueError(
            f"Input and output files must be different. "
            f"Both point to: {path.abspath(f61file)}"
        )
    
    # Validate input files exist
    if not path.exists(f61file):
        raise FileNotFoundError(f"Input fort.61.nc file not found: {f61file}")
    if not path.exists(f63file):
        raise FileNotFoundError(f"Input fort.63.nc file not found: {f63file}")
    
    print(f"Finding station '{station_name}' in {f61file}...")
    station_idx = find_station_index(f61file, station_name)
    print(f"  Found station at index {station_idx}")
    
    print(f"Extracting water level from {f63file} at coordinates ({correct_x}, {correct_y})...")
    try:
        f63_time, f63_wl = extract_f63_water_level(f63file, correct_x, correct_y)
        print(f"  Extracted {len(f63_wl)} time steps")
    except Exception as e:
        print(f"  ERROR: Failed to extract water level from {f63file}: {str(e)}")
        raise
    
    # Read fort.61.nc to get time axis
    with nc.Dataset(f61file) as src:
        # Read time - check if it has units attribute for proper conversion
        time_var = src['time']
        time_raw = time_var[:]
        
        # Check for time units and convert properly
        if hasattr(time_var, 'units'):
            # Use netCDF4's num2date for proper time conversion
            f61_time_cf = num2date(time_raw, units=time_var.units, calendar=getattr(time_var, 'calendar', 'standard'))
            # Convert cftime objects to pandas datetime
            try:
                # Try to convert directly - works for Python datetime objects
                f61_time = pd.to_datetime(f61_time_cf)
            except (TypeError, ValueError):
                # If that fails, it's likely cftime objects - convert via datetime property
                try:
                    f61_time = pd.to_datetime([dt.datetime for dt in f61_time_cf])
                except (AttributeError, TypeError):
                    # Last resort: convert to string then to datetime
                    f61_time = pd.to_datetime([str(dt) for dt in f61_time_cf])
        else:
            # Fallback to datetime64 conversion (may not work correctly)
            print(f"  WARNING: No time units attribute found, using datetime64 conversion")
            f61_time = pd.to_datetime(time_raw.astype("datetime64[ms]")).to_pydatetime()
        
        fill_value = src['zeta']._FillValue if hasattr(src['zeta'], '_FillValue') else -99999.0
    
    # Interpolate f63 data to match f61 time axis
    print("Interpolating fort.63.nc data to match fort.61.nc time axis...")
    f63_wl_interp = interpolate_to_f61_time(f61_time, f63_time, f63_wl)
    print(f"  Interpolated to {len(f63_wl_interp)} time steps (matching fort.61.nc)")
    
    # Create output file by copying input structure
    print(f"Creating output file {output_file}...")
    with nc.Dataset(f61file) as src, nc.Dataset(output_file, "w") as dst:
        # Copy global attributes
        dst.setncatts(src.__dict__)
        
        # Copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        
        # Copy all variables with their attributes (but not data yet)
        for name, variable in src.variables.items():
            dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name].setncatts(src[name].__dict__)
        
        # Copy all data from input file
        for name, variable in src.variables.items():
            dst[name][:] = src[name][:]
        
        # Replace zeta data for the target station
        print(f"Replacing zeta data for station {station_idx}...")
        zeta_data = dst['zeta'][:]
        zeta_data[:, station_idx] = f63_wl_interp
        
        # Handle NaN values - convert to fill value if needed
        if np.any(np.isnan(zeta_data[:, station_idx])):
            nan_mask = np.isnan(zeta_data[:, station_idx])
            zeta_data[nan_mask, station_idx] = fill_value
            print(f"  Warning: {np.sum(nan_mask)} NaN values converted to fill value {fill_value}")
        
        dst['zeta'][:] = zeta_data
        
        # Update x and y coordinates for the target station
        print(f"Updating coordinates for station {station_idx} to ({correct_x}, {correct_y})...")
        dst['x'][station_idx] = correct_x
        dst['y'][station_idx] = correct_y
    
    print(f"Successfully created corrected file: {output_file}")


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Replace incorrect water level data for a station in fort.61.nc "
                    "files with data extracted from corresponding fort.63.nc files "
                    "at corrected coordinates. Corresponding fort.63.nc files are "
                    "automatically found in the same directory as each fort.61.nc file."
    )
    parser.add_argument(
        "fort61_files",
        nargs='+',
        help="Paths to input fort.61.nc files (with wrong station data)"
    )
    parser.add_argument(
        "station_name",
        help="Name of the station to correct"
    )
    parser.add_argument(
        "correct_x",
        type=float,
        help="Correct x coordinate"
    )
    parser.add_argument(
        "correct_y",
        type=float,
        help="Correct y coordinate"
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip processing if output file already exists"
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    
    num_files = len(args.fort61_files)
    print(f"Replacing station data in {num_files} fort.61.nc file(s)...")
    print(f"  Station: {args.station_name}")
    print(f"  Correct coordinates: ({args.correct_x}, {args.correct_y})")
    if args.skip_existing:
        print(f"  Mode: Skip existing output files")
    if num_files > 1:
        print(f"  Progress will be reported for each file")
    print()
    
    success_count = 0
    error_count = 0
    skipped_count = 0
    
    for i, f61file in enumerate(args.fort61_files, 1):
        if num_files > 1:
            print(f"[{i}/{num_files}] Processing {f61file}...")
        else:
            print(f"Processing {f61file}...")
        
        try:
            # Generate output filename first to check if it exists
            output_file = generate_output_filename(f61file)
            
            # Check if output file already exists
            if args.skip_existing and path.exists(output_file):
                skipped_count += 1
                if num_files > 1:
                    print(f"  [{i}/{num_files}] ⊘ Skipped (output file already exists: {output_file})")
                else:
                    print(f"  ⊘ Skipped (output file already exists: {output_file})")
            else:
                # Find corresponding fort.63.nc file
                if num_files > 1:
                    print(f"  [{i}/{num_files}] Finding corresponding fort.63.nc file...")
                f63file = find_corresponding_f63_file(f61file)
                print(f"  Found corresponding fort.63.nc: {f63file}")
                print(f"  Output file: {output_file}")
                
                # Process the file
                replace_f61_station_with_f63(
                    f61file,
                    f63file,
                    args.station_name,
                    args.correct_x,
                    args.correct_y,
                    output_file
                )
                
                success_count += 1
                if num_files > 1:
                    print(f"  [{i}/{num_files}] ✓ Successfully processed {f61file}")
                else:
                    print(f"  ✓ Successfully processed {f61file}")
            
        except Exception as e:
            error_count += 1
            if num_files > 1:
                print(f"  [{i}/{num_files}] ✗ Error processing {f61file}: {str(e)}")
            else:
                print(f"  ✗ Error processing {f61file}: {str(e)}")
        
        if num_files > 1:
            print(f"  [{i}/{num_files}] Progress: {success_count} succeeded, {error_count} failed, {skipped_count} skipped so far")
        print()
    
    # Summary
    print("=" * 60)
    if args.skip_existing:
        print(f"Processing complete: {success_count}/{num_files} succeeded, {error_count}/{num_files} failed, {skipped_count}/{num_files} skipped")
    else:
        print(f"Processing complete: {success_count}/{num_files} succeeded, {error_count}/{num_files} failed")
    
    if error_count > 0:
        return 1
    return 0


if __name__ == "__main__":
    main()

