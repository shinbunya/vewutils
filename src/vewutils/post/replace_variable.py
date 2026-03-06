#!/usr/bin/env python3
"""
Module for replacing variable names in NetCDF files.
"""

import argparse
import xarray as xr
import os


def replace_variable(
    filename: str,
    renames: dict[str, str],
    output_filename: str,
) -> None:
    """
    Replace one or more variable names in a NetCDF file.

    Args:
        filename: Path to input NetCDF file
        renames: Mapping of old variable name -> new variable name (one or more)
        output_filename: Path to output file
    """
    if not renames:
        raise ValueError("renames must contain at least one old_name -> new_name mapping")

    print(f"Reading {filename}...")
    data = xr.open_dataset(filename)

    # Check that all old variables exist
    for old_name in renames:
        if old_name not in data.variables:
            print(f"Error: Variable '{old_name}' not found in the dataset.")
            print(f"Available variables: {list(data.variables.keys())}")
            data.close()
            return

    # Warn if any new name already exists (and is not just a no-op rename)
    for old_name, new_name in renames.items():
        if new_name in data.variables and old_name != new_name:
            print(f"Warning: Variable '{new_name}' already exists and will be overwritten.")

    for old_name, new_name in renames.items():
        print(f"Replacing variable '{old_name}' with '{new_name}'...")

    # Create a new dataset with all renames applied at once
    data_renamed = data.rename(renames)
    
    print(f"Saving to {output_filename}...")
    # Map current (possibly renamed) variable name back to original name in source dataset
    name_in_source = {new: old for old, new in renames.items()}
    for v in data_renamed.variables:
        if v not in name_in_source:
            name_in_source[v] = v

    # Create encoding that preserves original _FillValue and disables automatic scaling
    encoding = {}
    for var in data_renamed.variables:
        encoding[var] = {'zlib': False}
        src_name = name_in_source[var]
        original_var = data[src_name]
        if hasattr(original_var, '_FillValue') and original_var._FillValue is not None:
            encoding[var]['_FillValue'] = original_var._FillValue
        else:
            encoding[var]['_FillValue'] = None

        # Handle time variable encoding
        if var == 'time':
            original_time_var = data[name_in_source['time']]
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
    for old_name, new_name in renames.items():
        print(f"  '{old_name}' -> '{new_name}'")


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Replace one or more variable names in a NetCDF file.",
    )
    parser.add_argument(
        "filename",
        help="Path to input NetCDF file",
    )
    parser.add_argument(
        "old_var_name",
        nargs="?",
        default=None,
        help="Name of the variable to be replaced (for single rename; omit if using --old/--new)",
    )
    parser.add_argument(
        "new_var_name",
        nargs="?",
        default=None,
        help="New name for the variable (for single rename; omit if using --old/--new)",
    )
    parser.add_argument(
        "--old",
        nargs="+",
        metavar="VAR",
        help="Old variable names for multiple renames (use with --new; same order and length as --new)",
    )
    parser.add_argument(
        "--new",
        nargs="+",
        metavar="VAR",
        dest="new_vars",
        help="New variable names for multiple renames (use with --old; same order and length as --old)",
    )
    parser.add_argument(
        "-o", "--output",
        help="Path to output file (default: 'renamed_<input_filename>')",
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()

    # Build renames mapping: either (old_var_name, new_var_name) or (--old, --new)
    if args.old is not None and args.new_vars is not None:
        if len(args.old) != len(args.new_vars):
            print(
                f"Error: --old and --new must have the same number of names "
                f"(got {len(args.old)} and {len(args.new_vars)})."
            )
            return 1
        renames = dict(zip(args.old, args.new_vars))
    elif args.old_var_name is not None and args.new_var_name is not None:
        renames = {args.old_var_name: args.new_var_name}
    else:
        print(
            "Error: specify either two positionals (old_var_name new_var_name) "
            "or both --old and --new with the same number of arguments."
        )
        return 1

    print(f"Reading {args.filename} for variable information...")
    data = xr.open_dataset(args.filename)

    # Display available variables
    print(f"Available variables in the file: {list(data.variables.keys())}")

    # Validate input (replace_variable will also check, but we can fail early)
    for old_name in renames:
        if old_name not in data.variables:
            print(f"Error: Variable '{old_name}' not found in the dataset.")
            data.close()
            return 1

    output_filename = args.output if args.output else f"renamed_{os.path.basename(args.filename)}"

    print(f"Replacing variable name(s) in NetCDF file...")
    print(f"Input file: {args.filename}")
    for old_name, new_name in renames.items():
        print(f"  {old_name} -> {new_name}")
    print(f"Output file: {output_filename}")

    data.close()

    replace_variable(args.filename, renames, output_filename)

    print("Variable replacement completed successfully.")
    return 0


if __name__ == "__main__":
    main()
