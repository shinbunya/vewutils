#!/usr/bin/env python3
"""
Module for computing differences between two ADCIRC maxele.63.nc files.
"""

import argparse
import netCDF4 as nc
from os import path
from typing import List


def add_diff(ncfname1: str, ncfname2: str, ncfout: str, toexclude: List[str]) -> None:
    """
    Calculate differences between two maxele.63.nc files.
    
    Args:
        ncfname1: Path to first maxele file
        ncfname2: Path to second maxele file
        ncfout: Path to output file
        toexclude: List of variables to exclude from copying
    """
    with nc.Dataset(ncfname1) as src1, \
         nc.Dataset(ncfname2) as src2, \
         nc.Dataset(ncfout, "w") as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src1.__dict__)
        # copy dimensions
        for name, dimension in src1.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        # copy all file data except for the excluded
        for name, variable in src1.variables.items():
            if name not in toexclude:
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src1[name].__dict__)
                dst[name][:] = src1[name][:]

        # add diff
        dst_name = 'zeta_max_diff'
        x = dst.createVariable(dst_name, src1['zeta_max'].datatype, src1['zeta_max'].dimensions)
        attr = src1['zeta_max'].__dict__
        attr['long_name'] = 'maximum water surface elevation difference'
        attr['standard_name'] = 'zeta_max_diff'
        dst[dst_name].setncatts(attr)
        zeta_max1 = src1['zeta_max'][:]
        zeta_max2 = src2['zeta_max'][:]
        zeta_max_diff = zeta_max1 - zeta_max2
        dst[dst_name][:] = zeta_max_diff


def main():
    """Parse command line arguments and run the diff calculation."""
    parser = argparse.ArgumentParser(
        description="Calculate differences between two ADCIRC maxele files"
    )
    parser.add_argument(
        "maxele1",
        help="Path to first maxele file"
    )
    parser.add_argument(
        "maxele2",
        help="Path to second maxele file"
    )
    parser.add_argument(
        "output",
        help="Path to output maxele file with diff values"
    )
    parser.add_argument(
        "-e", "--exclude",
        default="",
        help="Comma-separated list of variables to exclude (default: none)"
    )
    
    args = parser.parse_args()
    
    # Process excluded variables
    toexclude = args.exclude.split(",") if args.exclude else []
    
    print(f"Computing differences between maxele files...")
    print(f"File 1: {args.maxele1}")
    print(f"File 2: {args.maxele2}")
    print(f"Output: {args.output}")
    print(f"Excluded variables: {toexclude}")
    
    add_diff(
        args.maxele1,
        args.maxele2,
        args.output,
        toexclude
    )
    
    print("Difference calculation completed successfully.")
    
    return 0


if __name__ == "__main__":
    main() 