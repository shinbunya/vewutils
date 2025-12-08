#!/usr/bin/env python3
"""
Module for adding departure field to ADCIRC maxele.63.nc files.
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
    Add departure field to a maxele.63.nc file.
    
    The departure is defined as the maximum water surface elevation 
    minus the reference water level. For dry nodes (depth < 0), the 
    reference value is the elevation (-depth).
    
    Args:
        ncfname: Path to input maxele file
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
        # Fallback: manual time comparison (same as plot_solution_2d.py)
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
    
    # Now process the maxele file
    with nc.Dataset(ncfname, "r") as src, nc.Dataset(ncfout, "w") as dst:
        # Check if required variables exist
        if 'zeta_max' not in src.variables:
            raise ValueError(f"'zeta_max' variable not found in {ncfname}")
        if 'depth' not in src.variables:
            raise ValueError(f"'depth' variable not found in {ncfname} (required for departure computation)")
        
        # Check dimensions match
        zeta_max = src.variables['zeta_max'][:]
        depth = src.variables['depth'][:]
        
        if len(zeta_max) != len(zeta_ref) or len(zeta_max) != len(depth):
            raise ValueError(
                f"Dimension mismatch for departure computation:\n"
                f"  zeta_max: {len(zeta_max)}, Reference: {len(zeta_ref)}, Depth: {len(depth)}"
            )
        
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
        reference_values = np.where(depth < 0, -depth, zeta_ref)
        
        # Compute departure: zeta_max - reference
        print(f"Computing departure field from 'zeta_max' and reference values")
        departure = zeta_max - reference_values
        
        # add departure
        dst_name = 'departure'
        x = dst.createVariable(dst_name, src['zeta_max'].datatype, src['zeta_max'].dimensions)
        attr = src['zeta_max'].__dict__.copy()
        attr['long_name'] = 'maximum water surface elevation departure from reference condition'
        attr['standard_name'] = 'departure'
        dst[dst_name].setncatts(attr)
        dst[dst_name][:] = departure


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False, 
        description="Add departure field to ADCIRC maxele file"
    )
    parser.add_argument(
        "input_maxele",
        help="Path to input maxele file"
    )
    parser.add_argument(
        "output_maxele",
        help="Path to output maxele file with departure field"
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
    print(f"Adding departure field to maxele file...")
    print(f"Input:  {args.input_maxele}")
    print(f"Output: {args.output_maxele}")
    print(f"Reference file: {args.reference_file}")
    print(f"Reference time: {args.reference_time}")
    print(f"Reference variable: {args.reference_variable}")
    print(f"Excluded variables: {toexclude}")
    try:
        add_departure(
            args.input_maxele,
            args.output_maxele,
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

