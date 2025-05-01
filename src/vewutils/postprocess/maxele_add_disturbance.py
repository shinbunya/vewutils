#!/usr/bin/env python3
"""
Module for adding disturbance field to ADCIRC maxele.63.nc files.
"""

import argparse
import netCDF4 as nc
from os import path
from typing import List


def add_disturbance(ncfname: str, ncfout: str, toexclude: List[str]) -> None:
    """
    Add disturbance field to a maxele.63.nc file.
    
    The disturbance is defined as the maximum water surface elevation 
    plus the negative of the depth where depth is negative (dry land).
    
    Args:
        ncfname: Path to input maxele file
        ncfout: Path to output file
        toexclude: List of variables to exclude from copying
    """
    with nc.Dataset(ncfname) as src, nc.Dataset(ncfout, "w") as dst:
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

        # add disturbance
        dst_name = 'disturbance'
        x = dst.createVariable(dst_name, src['zeta_max'].datatype, src['zeta_max'].dimensions)
        attr = src['zeta_max'].__dict__
        attr['long_name'] = 'maximum water surface elevation departure above initial condition'
        attr['standard_name'] = 'disturbance'
        dst[dst_name].setncatts(attr)
        zeta_max = src['zeta_max'][:]
        depth = src['depth'][:]
        disturbance = [zeta_max[i] + min(depth[i], 0.0) for i in range(len(zeta_max))]
        dst[dst_name][:] = disturbance


def main():
    """Parse command line arguments and run the disturbance calculation."""
    parser = argparse.ArgumentParser(
        description="Add disturbance field to ADCIRC maxele file"
    )
    parser.add_argument(
        "input_maxele",
        help="Path to input maxele file"
    )
    parser.add_argument(
        "output_maxele",
        help="Path to output maxele file with disturbance field"
    )
    parser.add_argument(
        "-e", "--exclude",
        default="",
        help="Comma-separated list of variables to exclude (default: none)"
    )
    
    args = parser.parse_args()
    
    # Process excluded variables
    toexclude = args.exclude.split(",") if args.exclude else []
    
    print(f"Adding disturbance field to maxele file...")
    print(f"Input:  {args.input_maxele}")
    print(f"Output: {args.output_maxele}")
    print(f"Excluded variables: {toexclude}")
    
    add_disturbance(
        args.input_maxele,
        args.output_maxele,
        toexclude
    )
    
    print("Disturbance calculation completed successfully.")
    
    return 0


if __name__ == "__main__":
    main() 