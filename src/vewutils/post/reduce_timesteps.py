#!/usr/bin/env python3
"""
Module for reducing time steps in NetCDF files based on a given interval.
"""

import argparse
import xarray as xr
from datetime import datetime
import os
import numpy as np


def reduce_timesteps(filename: str, start_time, end_time, interval: int, output_filename: str) -> None:
    """
    Reduce the time steps in a NetCDF file.
    
    Args:
        filename: Path to input NetCDF file
        start_time: Start time for slicing
        end_time: End time for slicing
        interval: Interval for selecting time steps
        output_filename: Path to output file
    """
    print(f"Reading {filename}...")
    # First read the dataset with decode_cf=True to get proper time handling
    data_decoded = xr.open_dataset(filename)
    
    print(f"Reducing time steps (start: {start_time}, end: {end_time}, interval: {interval})...")
    
    # Use the decoded dataset for selection
    if isinstance(start_time, datetime):
        # Convert datetime objects to np.datetime64 if needed
        start_time = np.datetime64(start_time)
    if isinstance(end_time, datetime):
        end_time = np.datetime64(end_time)
        
    # First select time range
    data_time_slice = data_decoded.sel(time=slice(start_time, end_time))
    
    # Then select every nth time step
    time_values = data_time_slice.time.values
    selected_times = time_values[::interval]
    data_reduced = data_time_slice.sel(time=selected_times)
    print(f"Selected {len(selected_times)} time steps")
    print(f"Times in reduced data: {data_reduced.time.values}")
    
    print(f"Saving to {output_filename}...")
    # Create encoding that preserves original _FillValue and disables automatic scaling
    encoding = {}
    for var in data_reduced.variables:
        encoding[var] = {'zlib': False}
        # Copy the original _FillValue if it exists
        if hasattr(data_decoded[var], '_FillValue') and data_decoded[var]._FillValue is not None:
            encoding[var]['_FillValue'] = data_decoded[var]._FillValue
        else:
            # If no _FillValue in original, don't add one
            encoding[var]['_FillValue'] = None
        
        # Special handling for zeta: always set _FillValue to -99999.
        if var == 'zeta':
            encoding[var]['_FillValue'] = -99999.
        
        # Inherit time base_date and create proper units
        if var == 'time':
            # Check for units in attrs dictionary first
            if 'units' in data_decoded[var].attrs:
                encoding[var]['units'] = data_decoded[var].attrs['units']
                print(f"Inheriting time units from source file: {encoding[var]['units']}")
            # Check for base_date and create units from it
            elif 'base_date' in data_decoded[var].attrs:
                base_date = data_decoded[var].attrs['base_date']
                encoding[var]['units'] = f'seconds since {base_date}'
                print(f"Creating time units from base_date: {encoding[var]['units']}")
            elif hasattr(data_decoded[var], 'units'):
                encoding[var]['units'] = data_decoded[var].units
                print(f"Inheriting time units from source file (via hasattr): {encoding[var]['units']}")
            else:
                print("Warning: Source file does not have time units or base_date specified")
                print(f"Available time attributes: {data_decoded[var].attrs}")
    
    data_reduced.to_netcdf(
        output_filename,
        format='NETCDF4_CLASSIC',
        engine='netcdf4',
        encoding=encoding
    )
    
    data_decoded.close()
    data_reduced.close()
    
    print(f"Successfully created {output_filename}")


def get_parser():
    parser = argparse.ArgumentParser(add_help=False, description="Reduce time steps in a NetCDF file based on a given interval")
    parser.add_argument(
        "filename",
        help="Path to input NetCDF file"
    )
    parser.add_argument(
        "--start-time",
        help="Start time for slicing (format: YYYY-MM-DDTHH:MM:SS). If not specified, uses first time in file.",
        default=None
    )
    parser.add_argument(
        "--end-time",
        help="End time for slicing (format: YYYY-MM-DDTHH:MM:SS). If not specified, uses last time in file.",
        default=None
    )
    parser.add_argument(
        "interval",
        type=int,
        help="Interval for selecting time steps"
    )
    parser.add_argument(
        "-o", "--output",
        help="Path to output file (default: 'reduced_<input_filename>')"
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    print(f"Reading {args.filename} for time information...")
    data = xr.open_dataset(args.filename)
    first_time = data.time.values[0]
    last_time = data.time.values[-1]
    start_time = None
    end_time = None
    if args.start_time:
        try:
            start_time = np.datetime64(datetime.strptime(args.start_time, '%Y-%m-%dT%H:%M:%S'))
        except ValueError:
            print("Error: Start time must be in YYYY-MM-DDTHH:MM:SS format")
            return 1
    else:
        start_time = first_time
        print(f"Using first time from file: {start_time}")
    if args.end_time:
        try:
            end_time = np.datetime64(datetime.strptime(args.end_time, '%Y-%m-%dT%H:%M:%S'))
        except ValueError:
            print("Error: End time must be in YYYY-MM-DDTHH:MM:SS format")
            return 1
    else:
        end_time = last_time
        print(f"Using last time from file: {end_time}")
    output_filename = args.output if args.output else f"reduced_{os.path.basename(args.filename)}"
    print(f"Reducing time steps in NetCDF file...")
    print(f"Input file: {args.filename}")
    print(f"Start time: {start_time}")
    print(f"End time: {end_time}")
    print(f"Interval: {args.interval}")
    print(f"Output file: {output_filename}")
    data.close()
    reduce_timesteps(
        args.filename,
        start_time,
        end_time,
        args.interval,
        output_filename
    )
    print("Time step reduction completed successfully.")
    return 0


if __name__ == "__main__":
    main() 