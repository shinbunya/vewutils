#!/usr/bin/env python3
"""
Script to adjust VEW boundary elevations by lowering channel node elevations relative to bank nodes.
"""

import argparse
from adcircpy import AdcircMesh
from vewutils.mesh.vew_boundary_manipulator import VEWBoundaryManipulator

def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Adjust VEW boundary elevations by lowering channel node elevations"
    )
    parser.add_argument(
        "input_mesh",
        help="Path to the input mesh file"
    )
    parser.add_argument(
        "-o", "--output",
        default="adjusted_mesh.14",
        help="Path to save the output mesh file (default: adjusted_mesh.14)"
    )
    parser.add_argument(
        "-t", "--tolerance",
        type=float,
        default=0.001,
        help="Amount to lower channel node elevations below bank node elevations (in meters, default: 0.001)"
    )
    return parser

def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    # Read the mesh file
    mesh = AdcircMesh.open(args.input_mesh)
    # Adjust VEW elevations
    mesh = VEWBoundaryManipulator.lower_channel_elevations_above_banks(mesh, args.tolerance)
    # Save the modified mesh
    mesh.write(args.output, overwrite=True)
    print(f"Modified mesh saved to: {args.output}")

if __name__ == "__main__":
    main() 