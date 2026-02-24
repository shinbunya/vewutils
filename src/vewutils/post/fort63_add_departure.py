#!/usr/bin/env python3
"""
Module for adding departure field to ADCIRC fort.63.nc files.
"""

import argparse
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
from os import path
from typing import List, Optional
from datetime import datetime


def add_departure(
    ncfname: str, 
    ncfout: str, 
    reference_file: str,
    reference_time: str,
    reference_variable: str = 'zeta',
    toexclude: List[str] = None
) -> None:
    """
    Add departure field to a fort.63.nc file.
    
    The departure is defined as the water surface elevation at each time step
    minus the reference water level. For dry nodes (depth < 0), the 
    reference value is the elevation (-depth).
    
    Args:
        ncfname: Path to input fort.63.nc file
        ncfout: Path to output file
        reference_file: Path to reference solution file (e.g., fort.63.nc)
        reference_time: Reference time in format "YYYY-MM-DD HH:mm:ss"
        reference_variable: Variable name to use from reference file (default: 'zeta')
        toexclude: List of variables to exclude from copying
    """
    if toexclude is None:
        toexclude = []
    
    # Check if reference file exists
    if not path.exists(reference_file):
        raise FileNotFoundError(f"Reference file {reference_file} not found.")
    
    # Parse reference time
    try:
        # Try parsing with pandas first
        try:
            ref_time = pd.to_datetime(reference_time)
        except:
            # Try datetime parsing
            ref_time = datetime.strptime(reference_time, "%Y-%m-%d %H:%M:%S")
            ref_time = pd.to_datetime(ref_time)
    except Exception as e:
        raise ValueError(
            f"Cannot parse reference_time '{reference_time}'. "
            f"Use format 'YYYY-MM-DD HH:mm:ss'. Error: {e}"
        )
    
    # Read reference file to find matching time step
    print(f"Reading reference data from {reference_file}")
    try:
        ds_ref = xr.open_dataset(reference_file)
    except Exception as e:
        raise ValueError(f"Error reading {reference_file}: {e}")
    
    # Check if reference variable exists
    if reference_variable not in ds_ref.variables:
        available_vars = list(ds_ref.variables.keys())
        ds_ref.close()
        raise ValueError(
            f"Variable '{reference_variable}' not found in {reference_file}. "
            f"Available variables: {available_vars}"
        )
    
    # Get time values from reference dataset
    # Use xarray's time selection which handles cftime objects automatically
    try:
        # Try using xarray's time selection with pandas Timestamp (handles cftime automatically)
        ds_ref_selected = ds_ref.sel(time=ref_time, method='nearest')
        # Get the selected time value for reporting
        selected_time_value = ds_ref_selected.time.values
        if hasattr(selected_time_value, 'item'):
            selected_time_value = selected_time_value.item()
        
        # Extract reference variable data directly from selected dataset
        # This handles all time formats including cftime automatically
        ref_var_selected = ds_ref_selected[reference_variable]
        if len(ref_var_selected.dims) >= 1:
            # Remove time dimension, get the data
            zeta_ref = ref_var_selected.values
            # If it's still multi-dimensional, flatten or take first slice
            if zeta_ref.ndim > 1:
                # Should be (node,) or (element, 1) after time selection
                if zeta_ref.shape[-1] == 1:
                    zeta_ref = zeta_ref[..., 0]
                else:
                    zeta_ref = zeta_ref.flatten() if zeta_ref.size == len(zeta_ref) else zeta_ref
        else:
            zeta_ref = ref_var_selected.values
        
        # Find the index for reporting purposes
        ref_time_values = ds_ref['time'].values
        try:
            # Try to find index by comparing values
            if isinstance(selected_time_value, np.datetime64):
                # For datetime64, use direct comparison
                ref_timestep = np.where(ref_time_values == selected_time_value)[0]
                if len(ref_timestep) > 0:
                    ref_timestep = ref_timestep[0]
                else:
                    ref_timestep = 0
            else:
                # For cftime or other types, try to find match
                matches = [i for i, t in enumerate(ref_time_values) 
                          if str(t) == str(selected_time_value) or t == selected_time_value]
                ref_timestep = matches[0] if matches else 0
        except Exception:
            ref_timestep = 0
            
    except Exception as e:
        # Fallback: manual time comparison (same as maxele_add_departure.py)
        print(f"Warning: Could not use xarray time selection. Using manual time comparison. Error: {e}")
        ref_time_values = ds_ref['time'].values
        
        # Convert to pandas Timestamp for comparison
        if isinstance(ref_time_values[0], np.datetime64):
            ref_time_pd = pd.to_datetime(ref_time)
            # Find closest time step using vectorized operations
            ref_time_values_pd = pd.to_datetime(ref_time_values)
            time_diffs = np.abs(ref_time_values_pd - ref_time_pd)
            ref_timestep = np.argmin(time_diffs)
        else:
            # Try to convert and find matching time
            print(f"Warning: Unexpected time format in reference file. Attempting to find matching time step.")
            try:
                # Check if we have cftime objects
                if hasattr(ref_time_values[0], 'datetime'):
                    # Convert cftime objects to Python datetime, then to pandas
                    ref_time_values_pd = pd.to_datetime([dt.datetime for dt in ref_time_values])
                    ref_time_pd = pd.to_datetime(ref_time)
                else:
                    # Try pandas conversion directly
                    ref_time_values_pd = pd.to_datetime(ref_time_values)
                    ref_time_pd = pd.to_datetime(ref_time)
                
                time_diffs = np.abs(ref_time_values_pd - ref_time_pd)
                ref_timestep = np.argmin(time_diffs)
            except Exception as e2:
                print(f"Warning: Could not convert time values. Using first time step. Error: {e2}")
                ref_timestep = 0
        
        # Extract reference variable data using integer index
        ref_shape = ds_ref[reference_variable].shape
        if len(ref_shape) == 2:
            # Variables with shape (time, node) - defined on nodes
            zeta_ref = ds_ref[reference_variable][ref_timestep, :].values
        elif len(ref_shape) == 3 and ref_shape[2] == 1:
            # Variables with shape (time, element, 1) - defined on elements
            zeta_ref = ds_ref[reference_variable][ref_timestep, :, 0].values
        else:
            # For other cases, try to get the data at the specified time step
            zeta_ref = ds_ref[reference_variable][ref_timestep, :].values
    
    print(f"Using reference time step {ref_timestep} (target: {reference_time})")
    print(f"Using reference variable '{reference_variable}' from reference file")
    
    # Close reference dataset
    ds_ref.close()
    
    # Now process the fort.63.nc file
    with nc.Dataset(ncfname, "r") as src, nc.Dataset(ncfout, "w") as dst:
        # Check if required variables exist
        if 'zeta' not in src.variables:
            raise ValueError(f"'zeta' variable not found in {ncfname}")
        
        # Check if depth exists (may not be present in all fort.63.nc files)
        has_depth = 'depth' in src.variables
        
        # Get dimensions
        zeta_var = src.variables['zeta']
        zeta_dims = zeta_var.dimensions
        
        # Verify zeta has expected dimensions (time, node)
        if len(zeta_dims) != 2:
            raise ValueError(
                f"Expected 'zeta' to have 2 dimensions (time, node), "
                f"but found {len(zeta_dims)}: {zeta_dims}"
            )
        
        # Get dimension sizes
        time_dim_name = zeta_dims[0]
        node_dim_name = zeta_dims[1]
        num_times = len(src.dimensions[time_dim_name])
        num_nodes = len(src.dimensions[node_dim_name])
        
        print(f"Input file dimensions: {num_times} time steps, {num_nodes} nodes")
        
        # Check dimensions match
        if len(zeta_ref) != num_nodes:
            raise ValueError(
                f"Dimension mismatch for departure computation:\n"
                f"  zeta: {num_nodes} nodes, Reference: {len(zeta_ref)} nodes"
            )
        
        # Get depth if available
        if has_depth:
            depth = src.variables['depth'][:]
            if len(depth) != num_nodes:
                raise ValueError(
                    f"Dimension mismatch: depth has {len(depth)} nodes, "
                    f"but zeta has {num_nodes} nodes"
                )
        else:
            print("Warning: 'depth' variable not found. Assuming all nodes are wet (depth >= 0).")
            print("         Reference values will use zeta_ref for all nodes.")
            depth = None
        
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            if name not in toexclude:
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)
                dst[name][:] = src[name][:]
        
        # Compute reference values: use zeta_ref for wet nodes, use elevation (-depth) for dry nodes
        # A node is dry if depth < 0
        if depth is not None:
            reference_values = np.where(depth < 0, -depth, zeta_ref)
        else:
            # No depth available, use zeta_ref for all nodes
            reference_values = zeta_ref
        
        # Compute departure for each time step: zeta[time, node] - reference[node]
        print(f"Computing departure field from 'zeta' and reference values for {num_times} time steps...")
        
        # Read zeta data time step by time step to avoid loading entire array into memory
        # Create departure variable
        dst_name = 'departure'
        x = dst.createVariable(dst_name, zeta_var.datatype, zeta_var.dimensions)
        attr = zeta_var.__dict__.copy()
        attr['long_name'] = 'water surface elevation departure from reference condition'
        attr['standard_name'] = 'departure'
        dst[dst_name].setncatts(attr)
        
        # Process time steps in chunks to balance memory usage and performance
        chunk_size = min(1000, num_times)  # Process up to 1000 time steps at a time
        for t_start in range(0, num_times, chunk_size):
            t_end = min(t_start + chunk_size, num_times)
            print(f"  Processing time steps {t_start} to {t_end-1}...")
            
            # Read chunk of zeta data
            zeta_chunk = zeta_var[t_start:t_end, :]
            
            # Compute departure for this chunk
            # Broadcasting: zeta_chunk is (chunk_size, num_nodes), reference_values is (num_nodes,)
            # Result is (chunk_size, num_nodes)
            departure_chunk = zeta_chunk - reference_values[np.newaxis, :]
            
            # Write chunk to output
            dst[dst_name][t_start:t_end, :] = departure_chunk
        
        print(f"Departure field computed and written for all {num_times} time steps.")


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False, 
        description="Add departure field to ADCIRC fort.63.nc file"
    )
    parser.add_argument(
        "input_fort63",
        help="Path to input fort.63.nc file"
    )
    parser.add_argument(
        "output_fort63",
        help="Path to output fort.63.nc file with departure field"
    )
    parser.add_argument(
        "--reference-file",
        required=True,
        help="Path to reference solution file (e.g., fort.63.nc)"
    )
    parser.add_argument(
        "--reference-time",
        required=True,
        help='Reference time in format "YYYY-MM-DD HH:mm:ss"'
    )
    parser.add_argument(
        "--reference-variable",
        default="zeta",
        help="Variable name to use from the reference file (default: zeta)"
    )
    parser.add_argument(
        "-e", "--exclude",
        default="",
        help="Comma-separated list of variables to exclude (default: none)"
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    toexclude = args.exclude.split(",") if args.exclude else []
    print(f"Adding departure field to fort.63.nc file...")
    print(f"Input:  {args.input_fort63}")
    print(f"Output: {args.output_fort63}")
    print(f"Reference file: {args.reference_file}")
    print(f"Reference time: {args.reference_time}")
    print(f"Reference variable: {args.reference_variable}")
    print(f"Excluded variables: {toexclude}")
    try:
        add_departure(
            args.input_fort63,
            args.output_fort63,
            args.reference_file,
            args.reference_time,
            args.reference_variable,
            toexclude
        )
        print("Departure calculation completed successfully.")
        return 0
    except Exception as e:
        print(f"Error: {e}")
        return 1


if __name__ == "__main__":
    main()
