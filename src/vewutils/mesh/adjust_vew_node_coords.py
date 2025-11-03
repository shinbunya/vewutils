#!/usr/bin/env python3
"""
Script to match coordinates of VEW boundary node pairs.
"""

import argparse
from adcircpy import AdcircMesh
from vewutils.mesh.vew_boundary_manipulator import VEWBoundaryManipulator

def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Match coordinates of VEW boundary node pairs"
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
        "-s", "--adopted-side",
        choices=['bank', 'channel'],
        default='bank',
        help="Which side's coordinates to adopt (default: bank). "
             "'bank' copies bank node coordinates to channel nodes, "
             "'channel' copies channel node coordinates to bank nodes."
    )
    return parser

def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    # Read the mesh file
    mesh = AdcircMesh.open(args.input_mesh)
    # Match node coordinates
    mesh = VEWBoundaryManipulator.match_vew_node_coordinates(mesh, args.adopted_side)
    # Save the modified mesh
    mesh.write(args.output, overwrite=True)
    print(f"Modified mesh saved to: {args.output}")

if __name__ == "__main__":
    main()

