#!/usr/bin/env python3
"""
Module for concatenating time series data from multiple ADCIRC fort.61.nc files.

Usage:
    python concat_fort61.py file1.fort.61.nc file2.fort.61.nc file3.fort.61.nc output.fort.61.nc

The input files should be provided in ascending chronological order.
"""

import argparse
import numpy as np
import netCDF4 as nc
from os import path
from typing import List, Dict, Tuple


def decode_station_names(station_name_array: np.ndarray) -> List[str]:
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


def build_station_registry(ncfnames: List[str]) -> Dict:
    """
    Pre-scan all files to build unified station registry.
    
    Returns dictionary with:
        - unified_stations: ordered list of all unique station names (first appearance order)
        - station_data: dict of {station_name: (x, y, first_file_idx)}
        - file_station_names: list of station name lists per file
        - num_stations_per_file: list of station counts per file
    """
    unified_stations = []  # Ordered list of unique stations
    station_data = {}  # {station_name: (x, y, first_file_idx)}
    file_station_names = []
    num_stations_per_file = []
    seen_stations = set()
    
    for file_idx, ncfname in enumerate(ncfnames):
        with nc.Dataset(ncfname) as src:
            station_names = decode_station_names(src['station_name'][:])
            x_coords = src['x'][:]
            y_coords = src['y'][:]
            
            file_station_names.append(station_names)
            num_stations_per_file.append(len(station_names))
            
            # Add new stations in order of first appearance
            for i, station_name in enumerate(station_names):
                if station_name not in seen_stations:
                    unified_stations.append(station_name)
                    station_data[station_name] = (
                        float(x_coords[i]),
                        float(y_coords[i]),
                        file_idx
                    )
                    seen_stations.add(station_name)
    
    return {
        'unified_stations': unified_stations,
        'station_data': station_data,
        'file_station_names': file_station_names,
        'num_stations_per_file': num_stations_per_file
    }


def create_output_structure(template_file: str, ncfout: str, station_registry: Dict) -> None:
    """
    Create output file structure with unified station list.
    
    Args:
        template_file: Path to first input file (used as template)
        ncfout: Path to output file
        station_registry: Station registry from build_station_registry()
    """
    unified_stations = station_registry['unified_stations']
    station_data = station_registry['station_data']
    num_unified = len(unified_stations)
    
    with nc.Dataset(template_file) as src, nc.Dataset(ncfout, "w") as dst:
        # Copy global attributes
        dst.setncatts(src.__dict__)
        
        # Create dimensions (update station dimension to unified count)
        for name, dimension in src.dimensions.items():
            if name == 'station':
                dst.createDimension('station', num_unified)
            else:
                dst.createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else None))
        
        # Create all variables
        for name, variable in src.variables.items():
            dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name].setncatts(src[name].__dict__)
        
        # Populate station metadata
        # station_name
        max_name_len = dst.dimensions['namelen'].size
        
        # Create properly formatted char array
        station_name_array = np.empty((num_unified, max_name_len), dtype='c')
        
        for i, station_name in enumerate(unified_stations):
            # Pad station name with spaces to max_name_len and encode
            padded_name = station_name.ljust(max_name_len).encode('utf-8')[:max_name_len]
            # Ensure we have exactly max_name_len bytes
            if len(padded_name) < max_name_len:
                padded_name = padded_name + b' ' * (max_name_len - len(padded_name))
            # Convert to numpy char array and assign
            char_array = np.frombuffer(padded_name, dtype='S1')
            station_name_array[i, :] = char_array
        
        dst['station_name'][:] = station_name_array
        
        # x and y coordinates
        for i, station_name in enumerate(unified_stations):
            x_coord, y_coord, _ = station_data[station_name]
            dst['x'][i] = x_coord
            dst['y'][i] = y_coord


def process_timeseries_data(
    ncfnames: List[str],
    ncfout: str,
    station_registry: Dict,
    fill_value: float = -99999.0
) -> None:
    """
    Read, remap, and write time series data to output file.
    
    Args:
        ncfnames: List of input file paths
        ncfout: Output file path (must already exist with structure)
        station_registry: Station registry from build_station_registry()
        fill_value: Value to use for missing stations
    """
    unified_stations = station_registry['unified_stations']
    file_station_names = station_registry['file_station_names']
    num_unified = len(unified_stations)
    
    # Build mapping from unified station name to index
    unified_station_to_idx = {name: i for i, name in enumerate(unified_stations)}
    
    time_arrays = []
    zeta_arrays = []
    
    for file_idx, ncfname in enumerate(ncfnames):
        with nc.Dataset(ncfname) as src:
            # Read time and zeta data
            time_data = src['time'][:]
            zeta_data = src['zeta'][:]  # shape: (num_times, num_stations_in_file)
            
            num_times = zeta_data.shape[0]
            file_stations = file_station_names[file_idx]
            
            # Create remapped array filled with fill_value
            zeta_remapped = np.full((num_times, num_unified), fill_value, dtype=zeta_data.dtype)
            
            # Map data from file stations to unified stations
            for file_station_idx, station_name in enumerate(file_stations):
                unified_idx = unified_station_to_idx[station_name]
                zeta_remapped[:, unified_idx] = zeta_data[:, file_station_idx]
            
            time_arrays.append(time_data)
            zeta_arrays.append(zeta_remapped)
    
    # Concatenate along time dimension
    concatenated_time = np.concatenate(time_arrays, axis=0)
    concatenated_zeta = np.concatenate(zeta_arrays, axis=0)
    
    # Write to output file
    with nc.Dataset(ncfout, "r+") as dst:
        dst['time'][:] = concatenated_time
        dst['zeta'][:] = concatenated_zeta


def _concatenate_fast(ncfnames: List[str], ncfout: str) -> None:
    """
    Fast concatenation assuming identical station lists.
    
    Raises:
        ValueError: If station counts differ between files
    """
    # Pre-check: verify all files have same number of stations
    station_counts = []
    for ncfname in ncfnames:
        with nc.Dataset(ncfname) as src:
            station_counts.append(len(src.dimensions['station']))
    
    if len(set(station_counts)) > 1:
        unique_counts = sorted(set(station_counts))
        raise ValueError(
            f"Station counts differ across files: {unique_counts}. "
            f"Use --mode robust to handle varying station lists."
        )
    
    # Variables that depend on time dimension
    time_dependent_vars = ['zeta']
    
    # Use first file as template
    with nc.Dataset(ncfnames[0]) as src1, nc.Dataset(ncfout, "w") as dst:
        
        # Copy global attributes all at once via dictionary
        dst.setncatts(src1.__dict__)
        
        # Copy dimensions
        for name, dimension in src1.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        
        # Copy all variables with their attributes (but not data yet)
        for name, variable in src1.variables.items():
            x = dst.createVariable(name, variable.datatype, variable.dimensions)
            # Copy variable attributes all at once via dictionary
            dst[name].setncatts(src1[name].__dict__)
        
        # Copy station-related variables from first file (these don't change)
        for name, variable in src1.variables.items():
            if name not in time_dependent_vars and name != 'time':
                dst[name][:] = src1[name][:]
        
        # Concatenate time-dependent data
        time_arrays = []
        zeta_arrays = []
        
        for ncfname in ncfnames:
            with nc.Dataset(ncfname) as src:
                time_arrays.append(src['time'][:])
                zeta_arrays.append(src['zeta'][:])
        
        # Concatenate along time dimension (axis 0)
        concatenated_time = np.concatenate(time_arrays, axis=0)
        concatenated_zeta = np.concatenate(zeta_arrays, axis=0)
        
        # Write concatenated data to output file
        dst['time'][:] = concatenated_time
        dst['zeta'][:] = concatenated_zeta


def _concatenate_robust(ncfnames: List[str], ncfout: str) -> None:
    """
    Robust concatenation that matches stations by name.
    Handles cases where station lists differ across files.
    """
    # Phase 1: Pre-scan all files to build station registry
    print("Phase 1/3: Scanning station metadata...")
    station_registry = build_station_registry(ncfnames)
    
    unified_stations = station_registry['unified_stations']
    num_unified = len(unified_stations)
    min_count = min(station_registry['num_stations_per_file'])
    max_count = max(station_registry['num_stations_per_file'])
    
    print(f"  Found {num_unified} unique stations across all files")
    print(f"  Station counts per file range: {min_count} to {max_count}")
    
    # Phase 2: Create output file with unified station list
    print("Phase 2/3: Creating output file structure...")
    create_output_structure(ncfnames[0], ncfout, station_registry)
    print(f"  Output file created with {num_unified} stations")
    
    # Phase 3: Process and write time series data
    print("Phase 3/3: Processing time series data...")
    
    # Get fill value from first file
    with nc.Dataset(ncfnames[0]) as src:
        fill_value = src['zeta']._FillValue if hasattr(src['zeta'], '_FillValue') else -99999.0
    
    process_timeseries_data(ncfnames, ncfout, station_registry, fill_value)
    
    print(f"  Time series data concatenated successfully")


def concatenate_fort61(ncfnames: List[str], ncfout: str, mode: str = 'fast') -> None:
    """
    Concatenate time series data from multiple fort.61.nc files.
    
    Args:
        ncfnames: List of paths to fort.61.nc files in ascending time order
        ncfout: Path to output file
        mode: 'fast' (assumes identical stations) or 'robust' (matches by name)
    
    Raises:
        FileExistsError: If the output file already exists
    """
    # Check if output file already exists
    if path.exists(ncfout):
        raise FileExistsError(f"Output file already exists: {ncfout}")
    
    if mode == 'fast':
        _concatenate_fast(ncfnames, ncfout)
    elif mode == 'robust':
        _concatenate_robust(ncfnames, ncfout)
    else:
        raise ValueError(f"Unknown mode: {mode}. Must be 'fast' or 'robust'.")


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False, 
        description="Concatenate time series data from multiple ADCIRC fort.61.nc files"
    )
    parser.add_argument(
        "fort61_files",
        nargs='+',
        help="Paths to fort.61.nc files in ascending time order (requires at least 2 files)"
    )
    parser.add_argument(
        "output",
        help="Path to output fort.61.nc file with concatenated time series"
    )
    parser.add_argument(
        "--mode",
        choices=['fast', 'robust'],
        default='fast',
        help="Concatenation mode: 'fast' assumes identical station lists (default), "
             "'robust' matches stations by name and handles varying station lists"
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    
    # Validate at least 2 files are provided
    if len(args.fort61_files) < 2:
        parser = get_parser()
        parser.error("At least 2 fort.61.nc files are required")
    
    print(f"Concatenating {len(args.fort61_files)} fort.61.nc files (mode: {args.mode})...")
    for i, fort61_file in enumerate(args.fort61_files, 1):
        print(f"File {i}: {fort61_file}")
    print(f"Output: {args.output}")
    
    concatenate_fort61(
        args.fort61_files,
        args.output,
        args.mode
    )
    
    print("Concatenation completed successfully.")
    return 0


if __name__ == "__main__":
    main()

