#!/usr/bin/env python3
"""
Module for rechunking NetCDF files to optimize for time-series access.
Uses NCO tools (nccopy) as primary method with Python fallback.
"""

import argparse
import os
import shutil
import subprocess
import sys

import numpy as np
import xarray as xr


def rechunk_with_nco(
    input_file: str,
    output_file: str,
    time_chunk_size: int,
    node_chunk_size: int = 1,
) -> bool:
    """
    Rechunk NetCDF file using NCO nccopy tool.
    
    Args:
        input_file: Path to input NetCDF file
        output_file: Path to output rechunked file
        time_chunk_size: Chunk size for time dimension (integer, use full dimension size for all)
        node_chunk_size: Chunk size for node dimension
        
    Returns:
        True if successful, False otherwise
    """
    nccopy_path = shutil.which("nccopy")
    if not nccopy_path:
        return False
    
    # Build chunk specification: time/<size>,node/<size>
    chunk_spec = f"time/{time_chunk_size},node/{node_chunk_size}"
    
    try:
        cmd = [nccopy_path, "-c", chunk_spec, input_file, output_file]
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
        print(f"Successfully rechunked using NCO nccopy")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running nccopy: {e}")
        print(f"stderr: {e.stderr}")
        return False
    except Exception as e:
        print(f"Unexpected error with nccopy: {e}")
        return False


def rechunk_with_python(
    input_file: str,
    output_file: str,
    time_chunk_size: int = -1,
    node_chunk_size: int = 1,
    preserve_filters: bool = True,
) -> None:
    """
    Rechunk NetCDF file using Python (xarray/netCDF4).
    
    Args:
        input_file: Path to input NetCDF file
        output_file: Path to output rechunked file
        time_chunk_size: Chunk size for time dimension (-1 means all time steps)
        node_chunk_size: Chunk size for node dimension
        preserve_filters: Whether to preserve compression filters from original file
    """
    print(f"Reading {input_file}...")
    ds = xr.open_dataset(input_file)
    
    print(f"Rechunking with Python (time={time_chunk_size}, node={node_chunk_size})...")
    
    # Determine actual time chunk size (handle -1 as "all" for that dimension)
    time_chunk = time_chunk_size
    if time_chunk == -1:
        # Use full dimension size
        if "time" in ds.sizes:
            time_chunk = ds.sizes["time"]
        elif "time" in ds.dims:
            # Fallback for older xarray versions
            time_chunk = ds.dims["time"] if isinstance(ds.dims["time"], int) else len(ds.dims["time"])
    
    # Build encoding with chunking specification
    encoding = {}
    for var_name, var in ds.variables.items():
        var_encoding = {}
        
        # Set chunking for variables that have time and node dimensions
        if "time" in var.dims and "node" in var.dims:
            var_encoding["chunksizes"] = (time_chunk, node_chunk_size)
        elif "time" in var.dims:
            # For variables with only time dimension, chunk by time
            var_encoding["chunksizes"] = (time_chunk,)
        elif "node" in var.dims:
            # For variables with only node dimension, chunk by node
            var_encoding["chunksizes"] = (node_chunk_size,)
        
        # Preserve compression filters if requested
        if preserve_filters:
            # Check if variable has filters in original file
            if hasattr(var, "encoding") and "zlib" in var.encoding:
                var_encoding["zlib"] = var.encoding["zlib"]
                if "complevel" in var.encoding:
                    var_encoding["complevel"] = var.encoding["complevel"]
                if "shuffle" in var.encoding:
                    var_encoding["shuffle"] = var.encoding["shuffle"]
            
            # Also check _FillValue
            if "_FillValue" in var.attrs:
                var_encoding["_FillValue"] = var.attrs["_FillValue"]
            elif hasattr(var, "encoding") and "_FillValue" in var.encoding:
                var_encoding["_FillValue"] = var.encoding["_FillValue"]
        
        encoding[var_name] = var_encoding
    
    print(f"Saving to {output_file}...")
    ds.to_netcdf(
        output_file,
        format="NETCDF4_CLASSIC",
        engine="netcdf4",
        encoding=encoding,
    )
    
    ds.close()
    print(f"Successfully created {output_file}")


def rechunk_netcdf(
    input_file: str,
    output_file: str,
    time_chunk_size: str = "ALL",
    node_chunk_size: int = 1,
    force_python: bool = False,
    preserve_filters: bool = True,
) -> None:
    """
    Rechunk NetCDF file to optimize for time-series access.
    Tries NCO first, falls back to Python if NCO is not available.
    
    Args:
        input_file: Path to input NetCDF file
        output_file: Path to output rechunked file
        time_chunk_size: Chunk size for time dimension ("ALL" for all time steps, or integer)
        node_chunk_size: Chunk size for node dimension (default: 1 for time-series optimization)
        force_python: Force use of Python implementation even if NCO is available
        preserve_filters: Preserve compression filters from original file
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    # First, we need to get the actual time dimension size for NCO
    # Read the file to get dimension info
    ds_temp = xr.open_dataset(input_file)
    actual_time_size = ds_temp.sizes.get("time", None)
    if actual_time_size is None and "time" in ds_temp.dims:
        # Fallback for older xarray versions
        actual_time_size = ds_temp.dims["time"] if isinstance(ds_temp.dims["time"], int) else len(ds_temp.dims["time"])
    ds_temp.close()
    
    # Convert time_chunk_size for both implementations
    if time_chunk_size == "ALL" or time_chunk_size == "all" or time_chunk_size == -1:
        python_time_chunk = -1
        # For NCO, use the actual dimension size
        if actual_time_size is None:
            raise ValueError("Could not determine time dimension size from file")
        nco_time_chunk = actual_time_size
    else:
        try:
            python_time_chunk = int(time_chunk_size)
            nco_time_chunk = int(time_chunk_size)
        except ValueError:
            raise ValueError(f"Invalid time_chunk_size: {time_chunk_size}. Use 'ALL' or an integer.")
    
    # Try NCO first unless forced to use Python
    if not force_python:
        print("Attempting to rechunk using NCO nccopy...")
        success = rechunk_with_nco(
            input_file,
            output_file,
            time_chunk_size=nco_time_chunk,
            node_chunk_size=node_chunk_size,
        )
        if success:
            return
    
    # Fall back to Python implementation
    print("Using Python implementation for rechunking...")
    rechunk_with_python(
        input_file,
        output_file,
        time_chunk_size=python_time_chunk,
        node_chunk_size=node_chunk_size,
        preserve_filters=preserve_filters,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Rechunk NetCDF file to optimize for time-series access",
    )
    parser.add_argument(
        "filename",
        help="Path to input NetCDF file",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to output rechunked file",
    )
    parser.add_argument(
        "--time-chunk-size",
        default="ALL",
        help='Chunk size for time dimension (default: "ALL" for all time steps, or specify integer)',
    )
    parser.add_argument(
        "--node-chunk-size",
        type=int,
        default=1,
        help="Chunk size for node dimension (default: 1 for time-series optimization)",
    )
    parser.add_argument(
        "--force-python",
        action="store_true",
        help="Force use of Python implementation even if NCO is available",
    )
    parser.add_argument(
        "--preserve-filters",
        action="store_true",
        default=True,
        help="Preserve compression filters from original file (default: True)",
    )
    parser.add_argument(
        "--no-preserve-filters",
        action="store_false",
        dest="preserve_filters",
        help="Do not preserve compression filters from original file",
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    
    print(f"Rechunking NetCDF file...")
    print(f"Input file: {args.filename}")
    print(f"Output file: {args.output}")
    print(f"Time chunk size: {args.time_chunk_size}")
    print(f"Node chunk size: {args.node_chunk_size}")
    print(f"Force Python: {args.force_python}")
    print(f"Preserve filters: {args.preserve_filters}")
    
    try:
        rechunk_netcdf(
            args.filename,
            args.output,
            time_chunk_size=args.time_chunk_size,
            node_chunk_size=args.node_chunk_size,
            force_python=args.force_python,
            preserve_filters=args.preserve_filters,
        )
        print("Rechunking completed successfully.")
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

