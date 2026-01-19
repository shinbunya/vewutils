#!/usr/bin/env python3
"""
Module for converting maxele.63.nc-like files to fort.63.nc-like format.
Takes multiple maxele.63.nc-like files and converts a specified variable
into a time series format like fort.63.nc.
"""

import argparse
import netCDF4 as nc
import numpy as np
from typing import List


def convert_maxelelike_to_fort63like(
    maxele_files: List[str],
    variable_name: str,
    output_file: str
) -> None:
    """
    Convert multiple maxele.63.nc-like files to a single fort.63.nc-like file.
    
    Args:
        maxele_files: List of paths to maxele.63.nc-like files
        variable_name: Name of variable to extract (e.g., 'zeta_max')
        output_file: Path to output fort.63.nc-like file
    """
    if not maxele_files:
        raise ValueError("At least one maxele file must be provided")
    
    num_files = len(maxele_files)
    print(f"Processing {num_files} maxele.63.nc-like file(s)...")
    
    # First pass: collect time values and verify structure
    time_values = []
    variable_data_list = []
    template_file = None
    num_nodes = None
    
    for file_idx, maxele_file in enumerate(maxele_files, 1):
        print(f"[{file_idx}/{num_files}] Reading: {maxele_file}")
        
        with nc.Dataset(maxele_file) as src:
            # Verify variable exists
            if variable_name not in src.variables:
                raise ValueError(f"Variable '{variable_name}' not found in {maxele_file}")
            
            var = src.variables[variable_name]
            
            # Get time value
            if 'time' not in src.variables:
                raise ValueError(f"No 'time' variable found in {maxele_file}")
            
            time_var = src.variables['time']
            time_data = time_var[:]
            
            # Get the time value (usually single element, but handle multiple)
            if time_data.size == 0:
                raise ValueError(f"Time variable is empty in {maxele_file}")
            time_value = float(time_data[0])  # Use first time value
            time_values.append(time_value)
            
            # Read variable data (preserve masking)
            var_data = var[:]
            
            # Check if variable has _FillValue and create proper mask if needed
            # netCDF4 may return MaskedArray with all False masks even when fill values exist
            if hasattr(var, '_FillValue') and var._FillValue is not None:
                fill_val = var._FillValue
                # If it's a MaskedArray, check if we need to update the mask
                if isinstance(var_data, np.ma.MaskedArray):
                    # Check if any values in the data match the fill value
                    data_values = var_data.data if hasattr(var_data, 'data') else var_data
                    fill_mask = (data_values == fill_val) | np.isnan(data_values)
                    # Update mask to include fill values
                    if np.any(fill_mask):
                        var_data = np.ma.MaskedArray(data_values, mask=fill_mask | var_data.mask, fill_value=fill_val)
                    else:
                        # Ensure fill_value matches
                        var_data = np.ma.MaskedArray(var_data.data, mask=var_data.mask, fill_value=fill_val)
                else:
                    # Regular array - check for fill values and create mask
                    fill_mask = (var_data == fill_val) | np.isnan(var_data)
                    if np.any(fill_mask):
                        var_data = np.ma.MaskedArray(var_data, mask=fill_mask, fill_value=fill_val)
            
            # Verify node dimension matches
            if num_nodes is None:
                num_nodes = var_data.shape[0]
                template_file = maxele_file
            elif var_data.shape[0] != num_nodes:
                raise ValueError(
                    f"Node dimension mismatch: {maxele_file} has {var_data.shape[0]} nodes, "
                    f"but template has {num_nodes} nodes"
                )
            
            variable_data_list.append(var_data)
    
    print(f"\nFound {num_nodes} nodes across all files")
    print(f"Time values: {time_values}")
    
    # Sort by time if needed (optional - could add flag for this)
    # For now, keep original order
    
    # Create output file using first file as template
    print(f"\nCreating output file: {output_file}")
    print(f"Using template: {template_file}")
    
    with nc.Dataset(template_file) as src_template, \
         nc.Dataset(output_file, "w") as dst:
        
        # Copy global attributes
        dst.setncatts(src_template.__dict__)
        
        # Copy dimensions
        # time dimension should be UNLIMITED for fort.63-like files
        for name, dimension in src_template.dimensions.items():
            if name == 'time':
                # Create UNLIMITED time dimension
                dst.createDimension('time', None)
            else:
                dst.createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else None))
        
        # Copy all variables except the one we're converting, time, and time_of_zeta_max
        # (time_of_zeta_max is maxele-specific and shouldn't be in fort.63-like files)
        toexclude = [variable_name, 'time', 'time_of_zeta_max']
        for name, variable in src_template.variables.items():
            if name not in toexclude:
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name].setncatts(variable.__dict__)
                dst[name][:] = variable[:]
        
        # Create time variable with multiple time steps
        print("Creating time variable...")
        time_template = src_template.variables['time']
        x = dst.createVariable('time', time_template.datatype, ('time',))
        dst['time'].setncatts(time_template.__dict__)
        dst['time'][:] = np.array(time_values, dtype=time_template.datatype)
        
        # Create zeta variable (fort.63-like) from the specified variable
        print(f"Creating zeta variable from '{variable_name}'...")
        var_template = src_template.variables[variable_name]
        
        # Get fill value from source variable or use default
        fill_value = var_template._FillValue if hasattr(var_template, '_FillValue') else -99999.0
        
        # Create zeta with (time, node) dimensions
        x = dst.createVariable('zeta', var_template.datatype, ('time', 'node'), fill_value=fill_value)
        
        # Copy attributes from source variable, but update for zeta
        attr = var_template.__dict__.copy()
        # Update attributes to match fort.63 format
        attr['long_name'] = 'water surface elevation above geoid'
        attr['standard_name'] = 'sea_surface_height_above_geoid'
        attr['coordinates'] = 'time y x'
        attr['location'] = 'node'
        attr['mesh'] = 'adcirc_mesh'
        if 'units' not in attr:
            attr['units'] = 'm'
        # Ensure _FillValue is set
        attr['_FillValue'] = fill_value
        
        dst['zeta'].setncatts(attr)
        
        # Stack variable data from all files into time series
        # Shape: (num_files, num_nodes) -> (time, node)
        # Preserve masking if any of the arrays are masked
        if any(isinstance(arr, np.ma.MaskedArray) for arr in variable_data_list):
            # Convert all to masked arrays if any are masked, using consistent fill value
            masked_list = []
            for arr in variable_data_list:
                if isinstance(arr, np.ma.MaskedArray):
                    # Ensure fill value matches
                    if arr.fill_value != fill_value:
                        arr = np.ma.MaskedArray(arr.data, mask=arr.mask, fill_value=fill_value)
                    masked_list.append(arr)
                else:
                    # Convert regular array to masked array (no mask)
                    masked_list.append(np.ma.MaskedArray(arr, fill_value=fill_value))
            zeta_data = np.ma.stack(masked_list, axis=0)
        else:
            zeta_data = np.stack(variable_data_list, axis=0)
        
        # Report masking information for each time step
        if isinstance(zeta_data, np.ma.MaskedArray):
            print(f"\nMasking information per time step:")
            for time_idx in range(zeta_data.shape[0]):
                num_masked = np.sum(zeta_data.mask[time_idx, :])
                num_total = zeta_data.shape[1]
                percentage = (num_masked / num_total * 100) if num_total > 0 else 0
                print(f"  Time step {time_idx + 1}/{zeta_data.shape[0]}: {num_masked:,} masked nodes ({percentage:.2f}%)")
        else:
            print(f"\nNo masking in output data (all {zeta_data.shape[1]:,} nodes have values at all {zeta_data.shape[0]} time steps)")
        
        # Write to NetCDF (masked values will be automatically handled with _FillValue)
        dst['zeta'][:] = zeta_data
        
        print(f"Successfully created fort.63-like file with {num_files} time steps")


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Convert maxele.63.nc-like files to fort.63.nc-like format"
    )
    parser.add_argument(
        "maxele_files",
        nargs='+',
        help="One or more paths to maxele.63.nc-like files"
    )
    parser.add_argument(
        "--variable",
        type=str,
        required=True,
        dest="variable_name",
        help="Name of variable to extract (e.g., 'zeta_max', 'zeta_max_diff')"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        dest="output",
        help="Path to output fort.63.nc-like file"
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    
    print(f"Converting maxele.63.nc-like files to fort.63.nc-like format...")
    print(f"Variable: {args.variable_name}")
    print(f"Output: {args.output}")
    print(f"Input files: {len(args.maxele_files)}")
    
    convert_maxelelike_to_fort63like(
        args.maxele_files,
        args.variable_name,
        args.output
    )
    
    print("Conversion completed successfully.")
    return 0


if __name__ == "__main__":
    main()
