#!/usr/bin/env python3
"""
Module for repairing corrupted VEW boundaries where channel and bank nodes are identical
or have incorrect elevation ordering.
"""

import argparse
import pandas as pd
import numpy as np
from adcircpy import AdcircMesh
from scipy.spatial import cKDTree
from typing import Dict, List, Tuple, Optional, Set


class VEWCorruptionError(Exception):
    """Custom exception for VEW boundary corruption that cannot be repaired."""
    pass


class VEWBoundaryRepairer:
    """Class for repairing corrupted VEW boundaries in ADCIRC meshes."""

    def __init__(self, mesh: AdcircMesh):
        """Initialize the VEW boundary repairer with a mesh."""
        self._mesh = mesh
        self._config = self._create_default_config()
        self._coordinate_map = None  # Lazy initialization

    @staticmethod
    def _create_default_config() -> Dict:
        """Create default configuration for VEW boundary repair."""
        return {
            'tolerance': 1e-6,       # Tolerance for identifying identical nodes
            'max_candidates': 10     # Maximum number of candidate nodes to consider
        }

    def with_tolerance(self, tolerance: float) -> 'VEWBoundaryRepairer':
        """Set the tolerance for identifying identical nodes."""
        self._config['tolerance'] = tolerance
        return self

    def with_max_candidates(self, max_candidates: int) -> 'VEWBoundaryRepairer':
        """Set the maximum number of candidate nodes to consider."""
        self._config['max_candidates'] = max_candidates
        return self

    def _build_coordinate_map(self) -> Dict[Tuple[float, float], List[int]]:
        """Build a mapping from coordinates to node IDs for efficient lookup."""
        if self._coordinate_map is not None:
            return self._coordinate_map
            
        print("Building coordinate lookup table...")
        coordinate_map = {}
        tolerance = self._config['tolerance']
        
        # Round coordinates to handle floating point precision
        decimal_places = max(0, int(abs(np.log10(tolerance))) + 1)
        
        for node_id, row in self._mesh.nodes.iterrows():
            # Round coordinates to specified precision
            rounded_x = round(row['x'], decimal_places)
            rounded_y = round(row['y'], decimal_places)
            coord_key = (rounded_x, rounded_y)
            
            if coord_key not in coordinate_map:
                coordinate_map[coord_key] = []
            coordinate_map[coord_key].append(node_id)
        
        # Filter to keep only coordinates with multiple nodes
        self._coordinate_map = {
            coord: nodes for coord, nodes in coordinate_map.items() 
            if len(nodes) > 1
        }
        
        print(f"Found {len(self._coordinate_map)} coordinate locations with duplicate nodes")
        return self._coordinate_map

    def _find_corrupted_vew_boundaries(self) -> List[Tuple[int, int, str]]:
        """Find VEW boundaries with corruption (identical nodes or incorrect ordering)."""
        corrupted_boundaries = []
        
        boundaries = self._mesh.boundaries.to_dict()
        if '64' not in boundaries:
            print("No VEW boundaries found in mesh")
            return corrupted_boundaries

        vew_boundaries = boundaries['64']
        
        for boundary_idx, boundary in enumerate(vew_boundaries):
            node_pairs = boundary['node_id']
            
            for pair_idx, (bank_node_str, channel_node_str) in enumerate(node_pairs):
                bank_node = int(bank_node_str)
                channel_node = int(channel_node_str)
                
                # Check for identical nodes
                if bank_node == channel_node:
                    corrupted_boundaries.append((boundary_idx, pair_idx, 'identical'))
                    print(f"Found identical node corruption: boundary {boundary_idx}, "
                          f"pair {pair_idx}, node {bank_node}")
                else:
                    # Check for incorrect elevation ordering
                    # Bank node should have higher elevation than channel node
                    bank_elevation = self._mesh.nodes.loc[bank_node, 'value_1']
                    channel_elevation = self._mesh.nodes.loc[channel_node, 'value_1']
                    
                    if bank_elevation <= channel_elevation:
                        corrupted_boundaries.append((boundary_idx, pair_idx, 'ordering'))
                        print(f"Found elevation ordering corruption: boundary {boundary_idx}, "
                              f"pair {pair_idx}, bank node {bank_node} (elev: {bank_elevation:.3f}) "
                              f"<= channel node {channel_node} (elev: {channel_elevation:.3f})")
        
        return corrupted_boundaries

    def _find_counterpart_node(self, node_id: int) -> Optional[int]:
        """Find the best counterpart node for a given node with identical coordinates."""
        # Build coordinate map if not already done
        coordinate_map = self._build_coordinate_map()
        
        # Get the coordinates of the target node
        target_node = self._mesh.nodes.loc[node_id]
        tolerance = self._config['tolerance']
        decimal_places = max(0, int(abs(np.log10(tolerance))) + 1)
        
        # Round coordinates to match the coordinate map
        rounded_x = round(target_node['x'], decimal_places)
        rounded_y = round(target_node['y'], decimal_places)
        coord_key = (rounded_x, rounded_y)
        
        # Check if there are other nodes at the same location
        if coord_key not in coordinate_map:
            return None
            
        candidate_nodes = [nid for nid in coordinate_map[coord_key] if nid != node_id]
        
        if not candidate_nodes:
            return None
        
        # If only one candidate, return it
        if len(candidate_nodes) == 1:
            return candidate_nodes[0]
            
        # If multiple candidates, find the one with most different elevation
        target_elevation = target_node['value_1']
        best_candidate = None
        max_elevation_diff = 0
        
        for candidate_id in candidate_nodes:
            candidate_elevation = self._mesh.nodes.loc[candidate_id, 'value_1']
            elevation_diff = abs(candidate_elevation - target_elevation)
            if elevation_diff > max_elevation_diff:
                max_elevation_diff = elevation_diff
                best_candidate = candidate_id
        
        return best_candidate

    def _determine_channel_bank_pair(self, node1_id: int, node2_id: int) -> Tuple[int, int]:
        """Determine which node is the bank and which is the channel based on elevation.
        
        Bank node should have higher elevation than channel node.
        Channel node should have lower elevation than bank node.
        """
        elevation1 = self._mesh.nodes.loc[node1_id, 'value_1']
        elevation2 = self._mesh.nodes.loc[node2_id, 'value_1']
        
        # Bank node should have higher elevation than channel node
        if elevation1 > elevation2:
            return node1_id, node2_id  # bank, channel
        else:
            return node2_id, node1_id  # bank, channel

    def repair_vew_boundaries(self) -> AdcircMesh:
        """Repair corrupted VEW boundaries in the mesh."""
        print("Starting VEW boundary repair...")
        
        # Find corrupted boundaries
        corrupted_boundaries = self._find_corrupted_vew_boundaries()
        
        if not corrupted_boundaries:
            print("No corrupted VEW boundaries found")
            return self._mesh
            
        print(f"Found {len(corrupted_boundaries)} corrupted VEW boundary pairs")
        
        # Create a copy of the mesh boundaries
        boundaries = self._mesh.boundaries.to_dict().copy()
        vew_boundaries = boundaries['64']
        
        # Track repairs made
        repairs_made = 0
        identical_repairs = 0
        ordering_repairs = 0
        
        # Process each corrupted boundary
        for boundary_idx, pair_idx, corruption_type in corrupted_boundaries:
            boundary = vew_boundaries[boundary_idx]
            bank_node_str, channel_node_str = boundary['node_id'][pair_idx]
            bank_node = int(bank_node_str)
            channel_node = int(channel_node_str)
            
            print(f"Repairing {corruption_type} corruption: boundary {boundary_idx}, "
                  f"pair {pair_idx}, nodes {bank_node}, {channel_node}")
            
            if corruption_type == 'identical':
                # Handle identical nodes
                corrupted_node = bank_node  # Same as channel_node
                counterpart_node = self._find_counterpart_node(corrupted_node)
                
                if counterpart_node is None:
                    raise VEWCorruptionError(
                        f"Could not find counterpart node for corrupted node {corrupted_node} "
                        f"in boundary {boundary_idx}, pair {pair_idx}. "
                        f"Try decreasing tolerance to find more precise matches."
                    )
                
                # Determine which is bank and which is channel based on elevation
                repaired_bank, repaired_channel = self._determine_channel_bank_pair(
                    corrupted_node, counterpart_node
                )
                identical_repairs += 1
                
            elif corruption_type == 'ordering':
                # Handle incorrect ordering - just swap them
                repaired_bank, repaired_channel = self._determine_channel_bank_pair(
                    bank_node, channel_node
                )
                ordering_repairs += 1
            
            print(f"  Repaired: bank node {repaired_bank} "
                  f"(elev: {self._mesh.nodes.loc[repaired_bank, 'value_1']:.3f}), "
                  f"channel node {repaired_channel} "
                  f"(elev: {self._mesh.nodes.loc[repaired_channel, 'value_1']:.3f})")
            
            # Update the boundary
            boundary['node_id'][pair_idx] = (str(repaired_bank), str(repaired_channel))
            repairs_made += 1
        
        print(f"Successfully repaired {repairs_made} VEW boundary pairs:")
        print(f"  - {identical_repairs} identical node repairs")
        print(f"  - {ordering_repairs} elevation ordering repairs")
        
        # Create new mesh with repaired boundaries
        repaired_mesh = AdcircMesh(
            nodes=self._mesh.nodes,
            elements=self._mesh.elements.elements,
            boundaries=boundaries
        )
        
        return repaired_mesh

    def validate_vew_boundaries(self) -> bool:
        """Validate that all VEW boundaries have distinct nodes with correct elevation ordering."""
        boundaries = self._mesh.boundaries.to_dict()
        if '64' not in boundaries:
            return True
            
        vew_boundaries = boundaries['64']
        validation_passed = True
        
        for boundary_idx, boundary in enumerate(vew_boundaries):
            node_pairs = boundary['node_id']
            
            for pair_idx, (bank_node_str, channel_node_str) in enumerate(node_pairs):
                bank_node = int(bank_node_str)
                channel_node = int(channel_node_str)
                
                # Check for identical nodes
                if bank_node == channel_node:
                    print(f"Validation failed: boundary {boundary_idx}, "
                          f"pair {pair_idx} has identical nodes {bank_node}")
                    validation_passed = False
                    continue
                
                # Check elevation ordering
                bank_elevation = self._mesh.nodes.loc[bank_node, 'value_1']
                channel_elevation = self._mesh.nodes.loc[channel_node, 'value_1']
                
                if bank_elevation <= channel_elevation:
                    print(f"Validation failed: boundary {boundary_idx}, pair {pair_idx} "
                          f"has incorrect elevation ordering - bank node {bank_node} "
                          f"(elev: {bank_elevation:.3f}) <= channel node {channel_node} "
                          f"(elev: {channel_elevation:.3f})")
                    validation_passed = False
        
        if validation_passed:
            print("VEW boundary validation passed")
        
        return validation_passed

    def get_repair_statistics(self) -> Dict:
        """Get statistics about VEW boundaries and potential repairs."""
        stats = {
            'total_vew_boundaries': 0,
            'total_node_pairs': 0,
            'identical_pairs': 0,
            'ordering_corrupted_pairs': 0,
            'repairable_identical_pairs': 0,
            'total_corrupted_pairs': 0
        }
        
        boundaries = self._mesh.boundaries.to_dict()
        if '64' not in boundaries:
            return stats
            
        vew_boundaries = boundaries['64']
        stats['total_vew_boundaries'] = len(vew_boundaries)
        
        for boundary in vew_boundaries:
            node_pairs = boundary['node_id']
            stats['total_node_pairs'] += len(node_pairs)
            
            for bank_node_str, channel_node_str in node_pairs:
                bank_node = int(bank_node_str)
                channel_node = int(channel_node_str)
                
                # Check for identical nodes
                if bank_node == channel_node:
                    stats['identical_pairs'] += 1
                    
                    # Check if repairable
                    counterpart = self._find_counterpart_node(bank_node)
                    if counterpart is not None:
                        stats['repairable_identical_pairs'] += 1
                else:
                    # Check for incorrect elevation ordering
                    bank_elevation = self._mesh.nodes.loc[bank_node, 'value_1']
                    channel_elevation = self._mesh.nodes.loc[channel_node, 'value_1']
                    
                    if bank_elevation <= channel_elevation:
                        stats['ordering_corrupted_pairs'] += 1
        
        stats['total_corrupted_pairs'] = stats['identical_pairs'] + stats['ordering_corrupted_pairs']
        
        return stats


def get_parser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        'input_mesh',
        help='Path to input mesh file with corrupted VEW boundaries'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output mesh file path (default: input_mesh with _repaired suffix)'
    )
    parser.add_argument(
        '-t', '--tolerance',
        type=float,
        default=1e-6,
        help='Tolerance for identifying identical nodes (default: 1e-6)'
    )
    parser.add_argument(
        '--max-candidates',
        type=int,
        default=10,
        help='Maximum number of candidate nodes to consider (default: 10)'
    )
    parser.add_argument(
        '--validate-only',
        action='store_true',
        help='Only validate VEW boundaries without repairing'
    )
    parser.add_argument(
        '--statistics',
        action='store_true',
        help='Show statistics about VEW boundaries'
    )
    parser.add_argument(
        '-d', '--description',
        default='repaired',
        help='Description for the output mesh (default: repaired)'
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    
    # Read the mesh file
    mesh = AdcircMesh.open(args.input_mesh)
    print(f"Successfully read mesh from: {args.input_mesh}")
    print(f"Number of nodes: {len(mesh.nodes)}")
    print(f"Number of elements: {len(mesh.elements.elements)}")
    
    # Create repairer
    repairer = (VEWBoundaryRepairer(mesh)
                .with_tolerance(args.tolerance)
                .with_max_candidates(args.max_candidates))
    
    # Show statistics if requested
    if args.statistics:
        print("\nVEW Boundary Statistics:")
        stats = repairer.get_repair_statistics()
        for key, value in stats.items():
            print(f"  {key.replace('_', ' ').title()}: {value}")
        print()
    
    # Validate only if requested
    if args.validate_only:
        print("Validating VEW boundaries...")
        is_valid = repairer.validate_vew_boundaries()
        if is_valid:
            print("All VEW boundaries are valid")
            return 0
        else:
            print("Found corrupted VEW boundaries")
            return 1
    
    # Repair the boundaries
    try:
        repaired_mesh = repairer.repair_vew_boundaries()
        
        # Validate the repaired mesh
        repairer_validation = VEWBoundaryRepairer(repaired_mesh)
        if not repairer_validation.validate_vew_boundaries():
            raise VEWCorruptionError("Repair process failed validation")
        
        # Determine output filename
        if args.output:
            output_path = args.output
        else:
            # Add _repaired suffix to input filename
            input_parts = args.input_mesh.rsplit('.', 1)
            if len(input_parts) == 2:
                output_path = f"{input_parts[0]}_repaired.{input_parts[1]}"
            else:
                output_path = f"{args.input_mesh}_repaired"
        
        # Set description and write
        repaired_mesh.description = args.description
        repaired_mesh.write(output_path, overwrite=True)
        print(f"\nRepaired mesh saved to: {output_path}")
        
        return 0
        
    except VEWCorruptionError as e:
        print(f"Error: {e}")
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    exit(main()) 