#!/usr/bin/env python3
"""
Module for converting data assimilation error values to maxele.63.nc format.
"""

import argparse
import netCDF4 as nc
import os
from datetime import datetime
import numpy as np


def read_error_surface(error_file: str, num_nodes: int) -> np.ndarray:
    """
    Read error surface file and return array of error values for all nodes.
    
    Args:
        error_file: Path to error surface file
        num_nodes: Number of nodes in the mesh
        
    Returns:
        Array of error values, one per node (0-indexed array of size num_nodes)
    """
    with open(error_file, 'r') as f:
        lines = f.readlines()
    
    # Line 3 (index 2) contains the default value
    default_value = float(lines[2].strip())
    
    # Initialize array with default values (0-indexed, size num_nodes)
    error_values = np.full(num_nodes, default_value, dtype=np.float64)
    
    # Lines 5 onwards (index 4+) contain node,value pairs
    # Nodes in file are 1-indexed, convert to 0-indexed for array
    for line in lines[4:]:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split(',')
        if len(parts) == 2:
            node_num = int(parts[0])  # 1-indexed node number
            error_val = float(parts[1])
            # Convert to 0-indexed array index
            error_values[node_num - 1] = error_val
    
    return error_values


def create_error_surface_nc(maxele_file: str, error_file: str, output_file: str, no_wetdry_mask: bool = False, reference_time: str = None) -> None:
    """
    Create a maxele.63.nc format file with error surface values.
    
    Args:
        maxele_file: Path to input maxele.63.nc file (used as template)
        error_file: Path to error surface file
        output_file: Path to output netcdf file
        no_wetdry_mask: If True, do not copy -99999.0 values from original zeta_max (default: False, so masking is copied by default)
        reference_time: Reference time string in format YYYY-MM-DDTHH:mm:ss to overwrite time base_date
    """
    with nc.Dataset(maxele_file) as src, \
         nc.Dataset(output_file, "w") as dst:
        # Get number of nodes
        num_nodes = len(src.dimensions['node'])
        
        # Read error surface values
        print(f"Reading error surface file: {error_file}")
        error_values = -1.0 * read_error_surface(error_file, num_nodes) # Flip the sign to make it model - obs
        
        # Copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        
        # Copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        
        # Copy all variables except zeta_max and time_of_zeta_max
        # Also handle time variable separately if reference_time is specified
        toexclude = ['zeta_max', 'time_of_zeta_max']
        if reference_time:
            toexclude.append('time')
        
        for name, variable in src.variables.items():
            if name not in toexclude:
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                # Copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)
                dst[name][:] = src[name][:]
        
        # Handle time variable if reference_time is specified
        if reference_time:
            print(f"Setting reference time to: {reference_time}")
            # Parse the reference time
            try:
                ref_dt = datetime.strptime(reference_time, "%Y-%m-%dT%H:%M:%S")
                ref_time_str = ref_dt.strftime("%Y-%m-%d %H:%M:%S")
            except ValueError:
                raise ValueError(f"Invalid reference time format: {reference_time}. Expected YYYY-MM-DDTHH:mm:ss")
            
            # Create time variable with new reference time
            x = dst.createVariable('time', src['time'].datatype, src['time'].dimensions)
            attr = src['time'].__dict__.copy()
            attr['base_date'] = ref_time_str
            attr['units'] = f"seconds since {ref_time_str}"
            dst['time'].setncatts(attr)
            # Keep the same time values (they're relative to the reference time)
            dst['time'][:] = src['time'][:]
        
        # Create zeta_max with error values
        print(f"Setting zeta_max to error surface values...")
        x = dst.createVariable('zeta_max', src['zeta_max'].datatype, src['zeta_max'].dimensions)
        attr = src['zeta_max'].__dict__.copy()
        # Update attributes to reflect that this is error surface data
        attr['long_name'] = 'data assimilation error surface'
        attr['standard_name'] = 'data_assimilation_error_surface'
        dst['zeta_max'].setncatts(attr)
        
        # Copy wet/dry values by default (unless disabled)
        if not no_wetdry_mask:
            print(f"Copying wet/dry values (-99999.0) from original maxele file...")
            original_zeta_max = src['zeta_max'][:]
            
            # Handle masked arrays (netCDF4 uses masked arrays for fill values)
            if isinstance(original_zeta_max, np.ma.MaskedArray):
                # Get the mask - True means the value is masked (dry node)
                dry_mask = original_zeta_max.mask
            else:
                # If not masked array, check for fill value or NaN
                fill_value = src['zeta_max']._FillValue if hasattr(src['zeta_max'], '_FillValue') else -99999.0
                dry_mask = (original_zeta_max == fill_value) | np.isnan(original_zeta_max)
            
            num_dry_nodes = np.sum(dry_mask)
            print(f"Found {num_dry_nodes} dry nodes")
            if num_dry_nodes > 0:
                # Print the first 10 indexes of the dry_mask where the values are true
                dry_indices = np.where(dry_mask)[0]
                print(f"First 10 dry node indices: {dry_indices[:10]}")
                # Preserve -99999.0 values for dry nodes
                error_values[dry_mask] = -99999.0
        
        dst['zeta_max'][:] = error_values
        
        # Create time_of_zeta_max with all zeros
        print(f"Setting time_of_zeta_max to zeros...")
        x = dst.createVariable('time_of_zeta_max', src['time_of_zeta_max'].datatype, src['time_of_zeta_max'].dimensions)
        attr = src['time_of_zeta_max'].__dict__.copy()
        # Update time_of_zeta_max attributes if reference_time is specified
        if reference_time:
            ref_dt = datetime.strptime(reference_time, "%Y-%m-%dT%H:%M:%S")
            ref_time_str = ref_dt.strftime("%Y-%m-%d %H:%M:%S")
            # Update units to match the new reference time if it exists
            if 'units' in attr:
                # Extract the units format and update with new reference time
                attr['units'] = f"seconds since {ref_time_str}"
        dst['time_of_zeta_max'].setncatts(attr)
        dst['time_of_zeta_max'][:] = np.zeros(num_nodes, dtype=np.float64)


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Convert data assimilation error surface to maxele.63.nc format"
    )
    parser.add_argument(
        "maxele_file",
        nargs='?',
        default=None,
        help="Path to input maxele.63.nc file (used as template). Required if --dir is not specified."
    )
    parser.add_argument(
        "error_file",
        nargs='?',
        default=None,
        help="Path to error surface file (da_error_surface.dat.0 format). Required if --dir is not specified."
    )
    parser.add_argument(
        "output",
        nargs='?',
        default=None,
        help="Path to output netcdf file. Required if --dir is not specified."
    )
    parser.add_argument(
        "--dir",
        type=str,
        nargs='+',
        default=None,
        help="Directory path(s). When specified, uses 'maxele.63.nc', 'da_error_surface.dat.1', and 'da_error_surface_model_minus_obs.1.nc' in each directory. Can specify multiple directories."
    )
    parser.add_argument(
        "--no-wetdry-mask",
        action="store_true",
        dest="no_wetdry_mask",
        help="Do not copy -99999.0 values from original maxele.63.nc (masking is copied by default)"
    )
    parser.add_argument(
        "--reference-time",
        type=str,
        default=None,
        help="Reference time in format YYYY-MM-DDTHH:mm:ss to overwrite time and time_of_zeta_max base_date"
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    
    # Determine file paths based on --dir option or positional arguments
    if args.dir:
        # Process multiple directories
        directories = args.dir
        num_dirs = len(directories)
        
        print(f"Processing {num_dirs} directory(ies)...")
        if args.no_wetdry_mask:
            print(f"Wet/dry mask copying: disabled")
        if args.reference_time:
            print(f"Reference time: {args.reference_time}")
        print()
        
        success_count = 0
        error_count = 0
        
        for idx, dir_path in enumerate(directories, 1):
            print(f"[{idx}/{num_dirs}] Processing directory: {dir_path}")
            
            maxele_file = os.path.join(dir_path, 'maxele.63.nc')
            error_file = os.path.join(dir_path, 'da_error_surface.dat.1')
            output_file = os.path.join(dir_path, 'da_error_surface_model_minus_obs.1.maxele.63.nc')
            
            try:
                create_error_surface_nc(
                    maxele_file,
                    error_file,
                    output_file,
                    no_wetdry_mask=args.no_wetdry_mask,
                    reference_time=args.reference_time
                )
                print(f"  ✓ Successfully created: {output_file}")
                success_count += 1
            except Exception as e:
                print(f"  ✗ Error processing {dir_path}: {e}")
                error_count += 1
            print()
        
        # Summary
        print(f"Processing complete: {success_count} succeeded, {error_count} failed")
        return 0 if error_count == 0 else 1
    else:
        # Use positional arguments (single directory mode)
        if args.maxele_file is None or args.error_file is None or args.output is None:
            parser = get_parser()
            parser.error("Either --dir must be specified, or all three positional arguments (maxele_file, error_file, output) must be provided")
        maxele_file = args.maxele_file
        error_file = args.error_file
        output_file = args.output
        
        print(f"Creating error surface netcdf file...")
        print(f"Template maxele file: {maxele_file}")
        print(f"Error surface file: {error_file}")
        print(f"Output file: {output_file}")
        if args.no_wetdry_mask:
            print(f"Wet/dry mask copying: disabled")
        if args.reference_time:
            print(f"Reference time: {args.reference_time}")
        
        create_error_surface_nc(
            maxele_file,
            error_file,
            output_file,
            no_wetdry_mask=args.no_wetdry_mask,
            reference_time=args.reference_time
        )
        
        print("Error surface netcdf file created successfully.")
        return 0


if __name__ == "__main__":
    main()
