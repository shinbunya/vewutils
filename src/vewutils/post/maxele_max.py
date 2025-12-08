#!/usr/bin/env python3
"""
Module for computing the maximum values across multiple ADCIRC maxele.63.nc files.
"""

import argparse
import numpy as np
import netCDF4 as nc
from os import path
from typing import List


def add_max_of_max(ncfnames: List[str], ncfout: str) -> None:
    """
    Calculate max of max values across multiple maxele.63.nc files.
    
    For each node, finds the maximum zeta_max value across all input files
    and records the corresponding time_of_zeta_max from the file where that
    maximum occurred.
    
    Args:
        ncfnames: List of paths to maxele files
        ncfout: Path to output file
    
    Raises:
        FileExistsError: If the output file already exists
    """
    # Check if output file already exists
    if path.exists(ncfout):
        raise FileExistsError(f"Output file already exists: {ncfout}")
    
    # Use first file as template
    with nc.Dataset(ncfnames[0]) as src1, \
         nc.Dataset(ncfout, "w") as dst:
        # Copy global attributes all at once via dictionary
        dst.setncatts(src1.__dict__)
        
        # Copy dimensions
        for name, dimension in src1.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        
        # Variables that will be computed (don't copy data from first file)
        computed_vars = ['zeta_max', 'time_of_zeta_max']
        
        # Create all variables and copy data for non-computed variables
        for name, variable in src1.variables.items():
            x = dst.createVariable(name, variable.datatype, variable.dimensions)
            # Copy variable attributes all at once via dictionary
            dst[name].setncatts(src1[name].__dict__)
            
            # Only copy data for non-computed variables
            if name not in computed_vars:
                dst[name][:] = src1[name][:]

        # Read both zeta_max and time_of_zeta_max from all files
        zeta_max_arrays = []
        time_of_zeta_max_arrays = []
        
        for ncfname in ncfnames:
            with nc.Dataset(ncfname) as src:
                zeta_max_arrays.append(src['zeta_max'][:])
                time_of_zeta_max_arrays.append(src['time_of_zeta_max'][:])
        
        # Stack into 2D arrays: (num_files, num_nodes)
        zeta_max_stack = np.array(zeta_max_arrays)
        time_of_zeta_max_stack = np.array(time_of_zeta_max_arrays)
        
        # Find which file has max for each node
        max_indices = np.nanargmax(zeta_max_stack, axis=0)
        
        # Get the actual max values
        zeta_max_of_max = np.nanmax(zeta_max_stack, axis=0)
        
        # Extract corresponding time values using advanced indexing
        time_of_max = np.take_along_axis(
            time_of_zeta_max_stack,
            max_indices[np.newaxis, :],
            axis=0
        ).squeeze()
        
        # Write computed values to output
        dst['zeta_max'][:] = zeta_max_of_max
        dst['time_of_zeta_max'][:] = time_of_max


def get_parser():
    parser = argparse.ArgumentParser(add_help=False, description="Calculate maximum values across multiple ADCIRC maxele files")
    parser.add_argument(
        "maxele_files",
        nargs='+',
        help="Paths to maxele files (requires at least 2 files)"
    )
    parser.add_argument(
        "output",
        help="Path to output maxele file with max values"
    )
    return parser

def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    
    # Validate at least 2 files are provided
    if len(args.maxele_files) < 2:
        parser = get_parser()
        parser.error("At least 2 maxele files are required")
    
    print(f"Computing max across {len(args.maxele_files)} maxele files...")
    for i, maxele_file in enumerate(args.maxele_files, 1):
        print(f"File {i}: {maxele_file}")
    print(f"Output: {args.output}")
    
    add_max_of_max(args.maxele_files, args.output)
    
    print("Max calculation completed successfully.")
    return 0


if __name__ == "__main__":
    main() 