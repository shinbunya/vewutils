#!/usr/bin/env python3
"""
NetCDF File Merger Script

This script reads a NetCDF file (File A), saves its data to a new NetCDF file (File C),
and copies specific variables from another NetCDF file (File B) to File C.

Variables to copy from File B: swan_rsx1, swan_rsy1, swan_rsx2, swan_rsy2
"""

import xarray as xr
import os
import sys
import argparse
from pathlib import Path


def merge_netcdf_files(file_a_path, file_b_path, file_c_path, variables_to_copy):
    """
    Merge NetCDF files by copying specific variables from File B to File C.
    
    Parameters:
    -----------
    file_a_path : str
        Path to the source NetCDF file (File A)
    file_b_path : str
        Path to the NetCDF file containing variables to copy (File B)
    file_c_path : str
        Path to the output NetCDF file (File C)
    variables_to_copy : list
        List of variable names to copy from File B to File C
    """
    
    print(f"Reading File A: {file_a_path}")
    if not os.path.exists(file_a_path):
        raise FileNotFoundError(f"File A not found: {file_a_path}")
    
    # Read File A
    ds_a = xr.open_dataset(file_a_path)
    print(f"File A loaded successfully. Variables: {list(ds_a.data_vars.keys())}")
    
    print(f"Reading File B: {file_b_path}")
    if not os.path.exists(file_b_path):
        raise FileNotFoundError(f"File B not found: {file_b_path}")
    
    # Read File B
    ds_b = xr.open_dataset(file_b_path)
    print(f"File B loaded successfully. Variables: {list(ds_b.data_vars.keys())}")
    
    # Check which variables exist in File B
    available_vars = []
    missing_vars = []
    
    for var in variables_to_copy:
        if var in ds_b.data_vars:
            available_vars.append(var)
            print(f"Found variable '{var}' in File B")
        else:
            missing_vars.append(var)
            print(f"Warning: Variable '{var}' not found in File B")
    
    if not available_vars:
        raise ValueError("None of the specified variables were found in File B")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(file_c_path)
    os.makedirs(output_dir, exist_ok=True)
    
    # Start with File A data
    ds_output = ds_a.copy()
    
    # Copy specified variables from File B (preserving time coordinate from File A)
    for var in available_vars:
        print(f"Copying variable '{var}' from File B to output")
        # Copy the variable data without its time coordinate to preserve File A's time
        var_data = ds_b[var].values
        # Create the variable with proper encoding to preserve _FillValue
        ds_output[var] = (ds_b[var].dims, var_data, ds_b[var].attrs)
        # Preserve the encoding attributes (including _FillValue)
        ds_output[var].encoding.update(ds_b[var].encoding)
    
    # Save the merged dataset to File C
    print(f"Saving merged data to File C: {file_c_path}")
    ds_output.to_netcdf(file_c_path)
    
    print(f"Successfully created File C with {len(ds_output.data_vars)} variables")
    print(f"Variables in output file: {list(ds_output.data_vars.keys())}")
    
    if missing_vars:
        print(f"Note: The following variables were not found in File B: {missing_vars}")
    
    # Close datasets to free memory
    ds_a.close()
    ds_b.close()
    ds_output.close()
    
    return True


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge NetCDF files by copying specific variables from File B to File C",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python netcdf_merge.py file_a.nc file_b.nc output.nc
  python netcdf_merge.py --file-a input1.nc --file-b input2.nc --file-c output.nc
  python netcdf_merge.py -a data1.nc -b data2.nc -c merged.nc
        """
    )
    
    # Add arguments for file paths
    parser.add_argument('file_a', nargs='?', 
                       help='Path to source NetCDF file (File A)')
    parser.add_argument('file_b', nargs='?',
                       help='Path to NetCDF file containing variables to copy (File B)')
    parser.add_argument('file_c', nargs='?',
                       help='Path to output NetCDF file (File C)')
    
    # Add optional arguments with flags
    parser.add_argument('-a', '--file-a', dest='file_a_flag',
                       help='Path to source NetCDF file (File A)')
    parser.add_argument('-b', '--file-b', dest='file_b_flag',
                       help='Path to NetCDF file containing variables to copy (File B)')
    parser.add_argument('-c', '--file-c', dest='file_c_flag',
                       help='Path to output NetCDF file (File C)')
    
    # Add optional argument for variables to copy
    parser.add_argument('-v', '--variables', nargs='+', 
                       default=['rsx1', 'rsy1', 'rsx2', 'rsy2', 'swan_rsx1', 'swan_rsy1', 'swan_rsx2', 'swan_rsy2'],
                       help='Variables to copy from File B (default: rsx1 rsy1 rsx2 rsy2 swan_rsx1 swan_rsy1 swan_rsx2 swan_rsy2)')
    
    return parser.parse_args()


def main():
    """Main function to execute the NetCDF merge operation."""
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Determine file paths (positional args take precedence over flag args)
    file_a = args.file_a or args.file_a_flag
    file_b = args.file_b or args.file_b_flag
    file_c = args.file_c or args.file_c_flag
    
    # Check if all required arguments are provided
    if not file_a or not file_b or not file_c:
        print("Error: All three file paths are required.")
        print("Usage: python netcdf_merge.py <file_a> <file_b> <file_c>")
        print("   or: python netcdf_merge.py -a <file_a> -b <file_b> -c <file_c>")
        print("   or: python netcdf_merge.py --file-a <file_a> --file-b <file_b> --file-c <file_c>")
        sys.exit(1)
    
    # Variables to copy from File B
    variables_to_copy = args.variables
    
    try:
        print("Starting NetCDF file merge operation...")
        print("=" * 50)
        print(f"File A (source): {file_a}")
        print(f"File B (variables source): {file_b}")
        print(f"File C (output): {file_c}")
        print(f"Variables to copy: {variables_to_copy}")
        print("=" * 50)
        
        success = merge_netcdf_files(file_a, file_b, file_c, variables_to_copy)
        
        if success:
            print("=" * 50)
            print("Operation completed successfully!")
        else:
            print("Operation failed!")
            sys.exit(1)
            
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
