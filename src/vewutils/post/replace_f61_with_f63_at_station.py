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
        
        # Find containing element
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
        
        # Time axis
        f63_time = pd.to_datetime(
            ds["time"].values.astype("datetime64[ms]")
        ).to_pydatetime()
        nt = len(f63_time)
    
    # Point not inside any triangle
    if not np.any(isinside):
        print(f"Warning: Point ({stx}, {sty}) is not inside any triangle in {f63file}")
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
    
    # Heavy part: read zeta with netCDF4 directly (time, node)
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
    
    # Vectorized interpolation over time
    f63_wl = a * elemv_all[:, 0] + b * elemv_all[:, 1] + c * elemv_all[:, 2]
    
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


def verify_time_match(f61_time: np.ndarray, f63_time: np.ndarray) -> None:
    """
    Verify that time axes from fort.61.nc and fort.63.nc match.
    
    Args:
        f61_time: Time array from fort.61.nc
        f63_time: Time array from fort.63.nc
    
    Raises:
        ValueError: If time axes don't match
    """
    if len(f61_time) != len(f63_time):
        raise ValueError(
            f"Time axes have different lengths: fort.61.nc has {len(f61_time)} timesteps, "
            f"fort.63.nc has {len(f63_time)} timesteps"
        )
    
    # Convert to comparable format (numpy datetime64)
    f61_time64 = pd.to_datetime(f61_time).values.astype("datetime64[ms]")
    f63_time64 = pd.to_datetime(f63_time).values.astype("datetime64[ms]")
    
    if not np.allclose(
        (f61_time64 - f63_time64).astype('timedelta64[ms]').astype(float),
        0.0,
        atol=1.0  # 1 ms tolerance
    ):
        raise ValueError(
            "Time axes do not match between fort.61.nc and fort.63.nc files. "
            "The files must have the same time steps."
        )


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
    f63_time, f63_wl = extract_f63_water_level(f63file, correct_x, correct_y)
    print(f"  Extracted {len(f63_wl)} time steps")
    
    # Read fort.61.nc to get time axis and verify match
    with nc.Dataset(f61file) as src:
        f61_time = src['time'][:]
        fill_value = src['zeta']._FillValue if hasattr(src['zeta'], '_FillValue') else -99999.0
    
    # Verify time axes match
    print("Verifying time axes match...")
    verify_time_match(f61_time, f63_time)
    print("  Time axes match")
    
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
        zeta_data[:, station_idx] = f63_wl
        
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
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    
    print(f"Replacing station data in {len(args.fort61_files)} fort.61.nc file(s)...")
    print(f"  Station: {args.station_name}")
    print(f"  Correct coordinates: ({args.correct_x}, {args.correct_y})")
    print()
    
    success_count = 0
    error_count = 0
    
    for i, f61file in enumerate(args.fort61_files, 1):
        print(f"[{i}/{len(args.fort61_files)}] Processing {f61file}...")
        
        try:
            # Find corresponding fort.63.nc file
            f63file = find_corresponding_f63_file(f61file)
            print(f"  Found corresponding fort.63.nc: {f63file}")
            
            # Generate output filename
            output_file = generate_output_filename(f61file)
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
            print(f"  ✓ Successfully processed {f61file}")
            
        except Exception as e:
            error_count += 1
            print(f"  ✗ Error processing {f61file}: {str(e)}")
        
        print()
    
    # Summary
    print("=" * 60)
    print(f"Processing complete: {success_count} succeeded, {error_count} failed")
    
    if error_count > 0:
        return 1
    return 0


if __name__ == "__main__":
    main()

