#!/usr/bin/env python3
"""
Module for adding attribution information to ADCIRC maxele.63.nc files.
"""

import argparse
import numpy as np
import netCDF4 as nc
from os import path
from typing import List


def add_attribution(ncfin_coastal: str, ncfin_runoff: str, ncfin_compound: str, 
                   ncfout: str, toexclude: List[str], threshold: float = 0.05) -> None:
    """
    Add attribution variables to a compound maxele file.
    
    Args:
        ncfin_coastal: Path to maxele file for coastal process
        ncfin_runoff: Path to maxele file for runoff process
        ncfin_compound: Path to maxele file for compound process
        ncfout: Path to output file
        toexclude: List of variables to exclude from copying
        threshold: Threshold for identifying compound effects
    """
    with nc.Dataset(ncfin_coastal) as src_coastal, \
         nc.Dataset(ncfin_runoff) as src_runoff, \
         nc.Dataset(ncfin_compound) as src_compound, \
         nc.Dataset(ncfout, "w") as dst:
        # copy global attributes all at once via dictionary
        src1 = src_compound
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
            if name == 'zeta_max':
                variable_zeta_max = variable
                
        # add values
        def add_variable(src, dst, dst_name, long_name, standard_name):
            x = dst.createVariable(dst_name, variable_zeta_max.datatype, 
                                  variable_zeta_max.dimensions, 
                                  fill_value=variable_zeta_max._FillValue)
            attr = src['zeta_max'].__dict__
            attr['long_name'] = long_name
            attr['standard_name'] = standard_name
            dst[dst_name].setncatts(attr)
        
        src1 = src_compound
        add_variable(src1, dst, 'zeta_max_coastal', 'maximum water surface elevation for coastal', 'zeta_max_coastal')
        add_variable(src1, dst, 'zeta_max_runoff', 'maximum water surface elevation for runoff', 'zeta_max_runoff')
        add_variable(src1, dst, 'zeta_max_compound', 'maximum water surface elevation for compound', 'zeta_max_compound')
        add_variable(src1, dst, 'zeta_max_individual', 'maximum water surface elevation for coastal and runoff', 'zeta_max_individual')
        add_variable(src1, dst, 'zeta_max_compound_minus_individual', 'maximum water surface elevation for compound minus individual', 'zeta_max_compound_minus_individual')
        add_variable(src1, dst, 'zeta_max_attribution', 'maximum water surface elevation attribution - 1: coastal, 2: coastal, compound, 3: runoff compound, 4: runoff', 'zeta_max_attribution')
        
        zeta_max_coastal = src_coastal['zeta_max'][:]
        zeta_max_runoff = src_runoff['zeta_max'][:]
        zeta_max_compound = src_compound['zeta_max'][:]

        zeta_max_coastal[np.where(zeta_max_coastal.mask)] = np.nan
        zeta_max_runoff[np.where(zeta_max_runoff.mask)] = np.nan
        zeta_max_compound[np.where(zeta_max_compound.mask)] = np.nan
        
        zeta_max_coastal[zeta_max_coastal < -99998.0] = np.nan
        zeta_max_runoff[zeta_max_runoff < -99998.0] = np.nan
        zeta_max_compound[zeta_max_compound < -99998.0] = np.nan

        zeta_max_individual = np.nanmax([zeta_max_coastal, zeta_max_runoff], axis=0)
        zeta_max_compound_minus_individual = zeta_max_compound - zeta_max_individual

        nn = len(zeta_max_coastal)
        iscoastal = np.array([1 if zeta_max_coastal[i] >= zeta_max_runoff[i] or 
                            (not np.isnan(zeta_max_coastal[i]) and np.isnan(zeta_max_runoff[i])) 
                            else 0 for i in range(nn)]).astype(np.float64)
        iscompound = np.array([1 if zeta_max_compound[i] - threshold >= max(zeta_max_coastal[i], zeta_max_runoff[i]) 
                             else 0 for i in range(nn)]).astype(np.float64)
        zeta_max_attribution = np.array([1.0 if iscoa == 1 and iscom == 0 
                                       else 2.0 if iscoa == 1 and iscom == 1 
                                       else 3.0 if iscoa == 0 and iscom == 1 
                                       else 4.0 for iscoa, iscom in zip(iscoastal, iscompound)])
        
        # Handle nan values in attribution
        zeta_max_attribution[np.all([np.isnan(zeta_max_coastal), iscoastal == 1.0], axis=0)] = -99999.0
        zeta_max_attribution[np.all([np.isnan(zeta_max_runoff), iscoastal == 0.0], axis=0)] = -99999.0
        
        dst['zeta_max_coastal'][:] = zeta_max_coastal
        dst['zeta_max_runoff'][:] = zeta_max_runoff
        dst['zeta_max_compound'][:] = zeta_max_compound
        dst['zeta_max_individual'][:] = zeta_max_individual
        dst['zeta_max_compound_minus_individual'][:] = zeta_max_compound_minus_individual
        dst['zeta_max_attribution'][:] = zeta_max_attribution


def main():
    """Parse command line arguments and run the attribution process."""
    parser = argparse.ArgumentParser(
        description="Add attribution information to ADCIRC maxele files"
    )
    parser.add_argument(
        "coastal_maxele",
        help="Path to maxele file for coastal process"
    )
    parser.add_argument(
        "runoff_maxele",
        help="Path to maxele file for runoff process"
    )
    parser.add_argument(
        "compound_maxele",
        help="Path to maxele file for compound process"
    )
    parser.add_argument(
        "output_maxele",
        help="Path to output maxele file with attribution"
    )
    parser.add_argument(
        "-e", "--exclude",
        default="zeta_max,time_of_zeta_max",
        help="Comma-separated list of variables to exclude (default: zeta_max,time_of_zeta_max)"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=0.05,
        help="Threshold for identifying compound effects (default: 0.05)"
    )
    
    args = parser.parse_args()
    
    # Process excluded variables
    toexclude = args.exclude.split(",") if args.exclude else []
    
    print(f"Adding attribution to maxele file...")
    print(f"Coastal:  {args.coastal_maxele}")
    print(f"Runoff:   {args.runoff_maxele}")
    print(f"Compound: {args.compound_maxele}")
    print(f"Output:   {args.output_maxele}")
    print(f"Excluded variables: {toexclude}")
    print(f"Threshold: {args.threshold}")
    
    add_attribution(
        args.coastal_maxele,
        args.runoff_maxele,
        args.compound_maxele,
        args.output_maxele,
        toexclude,
        args.threshold
    )
    
    print("Attribution completed successfully.")
    
    return 0


if __name__ == "__main__":
    main() 