#!/usr/bin/env python3
"""
Module for replacing variable names in NetCDF files.
"""

import argparse
import xarray as xr
import os


def replace_variable(filename: str, old_var_name: str, new_var_name: str, output_filename: str) -> None:
    """
    Replace a variable name in a NetCDF file.
    
    Args:
        filename: Path to input NetCDF file
        old_var_name: Name of the variable to be replaced
        new_var_name: New name for the variable
        output_filename: Path to output file
    """
    print(f"Reading {filename}...")
    # Read the dataset
    data = xr.open_dataset(filename)
    
    # Check if the old variable exists
    if old_var_name not in data.variables:
        print(f"Error: Variable '{old_var_name}' not found in the dataset.")
        print(f"Available variables: {list(data.variables.keys())}")
        data.close()
        return
    
    # Check if the new variable name already exists
    if new_var_name in data.variables:
        print(f"Warning: Variable '{new_var_name}' already exists in the dataset.")
        print("The existing variable will be overwritten.")
    
    print(f"Replacing variable '{old_var_name}' with '{new_var_name}'...")
    
    # Create a new dataset with the renamed variable
    data_renamed = data.rename({old_var_name: new_var_name})
    
    print(f"Saving to {output_filename}...")
    # Create encoding that preserves original _FillValue and disables automatic scaling
    encoding = {}
    for var in data_renamed.variables:
        encoding[var] = {'zlib': False}
        # Copy the original _FillValue if it exists
        if hasattr(data[var] if var in data.variables else data_renamed[var], '_FillValue'):
            original_var = data[var] if var in data.variables else data_renamed[var]
            if original_var._FillValue is not None:
                encoding[var]['_FillValue'] = original_var._FillValue
            else:
                encoding[var]['_FillValue'] = None
        else:
            # If no _FillValue in original, don't add one
            encoding[var]['_FillValue'] = None
        
        # Handle time variable encoding
        if var == 'time':
            original_time_var = data['time'] if 'time' in data.variables else data_renamed['time']
            # Check for units in attrs dictionary first
            if 'units' in original_time_var.attrs:
                encoding[var]['units'] = original_time_var.attrs['units']
                print(f"Inheriting time units from source file: {encoding[var]['units']}")
            # Check for base_date and create units from it
            elif 'base_date' in original_time_var.attrs:
                base_date = original_time_var.attrs['base_date']
                encoding[var]['units'] = f'seconds since {base_date}'
                print(f"Creating time units from base_date: {encoding[var]['units']}")
            elif hasattr(original_time_var, 'units'):
                encoding[var]['units'] = original_time_var.units
                print(f"Inheriting time units from source file (via hasattr): {encoding[var]['units']}")
    
    data_renamed.to_netcdf(
        output_filename,
        format='NETCDF4_CLASSIC',
        engine='netcdf4',
        encoding=encoding
    )
    
    data.close()
    data_renamed.close()
    
    print(f"Successfully created {output_filename}")
    print(f"Variable '{old_var_name}' has been renamed to '{new_var_name}'")


def get_parser():
    parser = argparse.ArgumentParser(add_help=False, description="Replace a variable name in a NetCDF file")
    parser.add_argument(
        "filename",
        help="Path to input NetCDF file"
    )
    parser.add_argument(
        "old_var_name",
        help="Name of the variable to be replaced"
    )
    parser.add_argument(
        "new_var_name",
        help="New name for the variable"
    )
    parser.add_argument(
        "-o", "--output",
        help="Path to output file (default: 'renamed_<input_filename>')"
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    
    print(f"Reading {args.filename} for variable information...")
    data = xr.open_dataset(args.filename)
    
    # Display available variables
    print(f"Available variables in the file: {list(data.variables.keys())}")
    
    # Validate input
    if args.old_var_name not in data.variables:
        print(f"Error: Variable '{args.old_var_name}' not found in the dataset.")
        data.close()
        return 1
    
    if args.new_var_name in data.variables and args.old_var_name != args.new_var_name:
        print(f"Warning: Variable '{args.new_var_name}' already exists and will be overwritten.")
    
    output_filename = args.output if args.output else f"renamed_{os.path.basename(args.filename)}"
    
    print(f"Replacing variable name in NetCDF file...")
    print(f"Input file: {args.filename}")
    print(f"Old variable name: {args.old_var_name}")
    print(f"New variable name: {args.new_var_name}")
    print(f"Output file: {output_filename}")
    
    data.close()
    
    replace_variable(
        args.filename,
        args.old_var_name,
        args.new_var_name,
        output_filename
    )
    
    print("Variable replacement completed successfully.")
    return 0


if __name__ == "__main__":
    main()
