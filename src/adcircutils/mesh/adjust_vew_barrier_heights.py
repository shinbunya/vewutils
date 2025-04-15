#!/usr/bin/env python3
"""
Script to ensure VEW boundary barrier heights are above bank elevations.
"""

import argparse
from adcircpy import AdcircMesh
from adcircutils.channelmodeling.vew_boundary_manipulator import VEWBoundaryManipulator

def main():
    """Main function to handle command line arguments and process the mesh."""
    parser = argparse.ArgumentParser(
        description="Ensure VEW boundary barrier heights are above bank elevations"
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
        help="Minimum amount that barrier heights should be above bank elevations (in meters, default: 0.001)"
    )
    
    args = parser.parse_args()
    
    # Read the mesh file
    mesh = AdcircMesh.open(args.input_mesh)
    
    # Adjust barrier heights
    mesh = VEWBoundaryManipulator.ensure_barrier_heights_above_banks(mesh, args.tolerance)
    
    # Save the modified mesh
    mesh.write(args.output, overwrite=True)
    print(f"Modified mesh saved to: {args.output}")

if __name__ == "__main__":
    main() 