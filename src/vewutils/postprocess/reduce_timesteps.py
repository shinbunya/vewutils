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
    
    data_reduced.to_netcdf(
        output_filename,
        format='NETCDF4_CLASSIC',
        engine='netcdf4',
        encoding=encoding
    )
    
    data_decoded.close()
    data_reduced.close()
    
    print(f"Successfully created {output_filename}")


def main():
    """Parse command line arguments and run the time step reduction."""
    parser = argparse.ArgumentParser(
        description="Reduce time steps in a NetCDF file based on a given interval"
    )
    parser.add_argument(
        "filename",
        help="Path to input NetCDF file"
    )
    parser.add_argument(
        "--start_time",
        help="Start time for slicing (format: YYYY-MM-DD). If not specified, uses first time in file.",
        default=None
    )
    parser.add_argument(
        "--end_time",
        help="End time for slicing (format: YYYY-MM-DD). If not specified, uses last time in file.",
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
    
    args = parser.parse_args()
    
    # Load dataset to get time information if needed
    print(f"Reading {args.filename} for time information...")
    # Always decode CF conventions for consistent datetime handling
    data = xr.open_dataset(args.filename)
    
    # Get first and last times from the dataset
    first_time = data.time.values[0]
    last_time = data.time.values[-1]
    
    # Convert string times to datetime if provided, otherwise use dataset times
    start_time = None
    end_time = None
    
    if args.start_time:
        try:
            start_time = np.datetime64(datetime.strptime(args.start_time, '%Y-%m-%d'))
        except ValueError:
            print("Error: Start time must be in YYYY-MM-DD format")
            return 1
    else:
        start_time = first_time
        print(f"Using first time from file: {start_time}")
    
    if args.end_time:
        try:
            end_time = np.datetime64(datetime.strptime(args.end_time, '%Y-%m-%d'))
        except ValueError:
            print("Error: End time must be in YYYY-MM-DD format")
            return 1
    else:
        end_time = last_time
        print(f"Using last time from file: {end_time}")
    
    # Set default output filename if not provided
    output_filename = args.output if args.output else f"reduced_{os.path.basename(args.filename)}"
    
    print(f"Reducing time steps in NetCDF file...")
    print(f"Input file: {args.filename}")
    print(f"Start time: {start_time}")
    print(f"End time: {end_time}")
    print(f"Interval: {args.interval}")
    print(f"Output file: {output_filename}")
    
    # Close the dataset opened for checking times
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