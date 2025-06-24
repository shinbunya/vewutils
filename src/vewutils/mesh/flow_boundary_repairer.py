#!/usr/bin/env python3
"""
Module for repairing corrupted flow boundaries where bank nodes are used instead of channel nodes.
"""

import argparse
import pandas as pd
import numpy as np
from adcircpy import AdcircMesh
from typing import Dict, List, Tuple, Optional, Set


class FlowCorruptionError(Exception):
    """Custom exception for flow boundary corruption that cannot be repaired."""
    pass


class FlowBoundaryRepairer:
    """Class for repairing corrupted flow boundaries in ADCIRC meshes."""

    def __init__(self, mesh: AdcircMesh):
        """Initialize the flow boundary repairer with a mesh."""
        self._mesh = mesh
        self._vew_node_map = None  # Lazy initialization

    def _build_vew_node_map(self) -> Dict[int, int]:
        """Build a mapping from bank nodes to channel nodes from VEW boundaries."""
        if self._vew_node_map is not None:
            return self._vew_node_map
            
        print("Building VEW node mapping (bank -> channel)...")
        vew_node_map = {}
        
        boundaries = self._mesh.boundaries.to_dict()
        if '64' not in boundaries:
            print("No VEW boundaries found in mesh")
            self._vew_node_map = {}
            return self._vew_node_map

        vew_boundaries = boundaries['64']
        
        for boundary_idx, boundary in enumerate(vew_boundaries):
            node_pairs = boundary['node_id']
            
            for pair_idx, (bank_node_str, channel_node_str) in enumerate(node_pairs):
                bank_node = int(bank_node_str)
                channel_node = int(channel_node_str)
                
                # Map bank node to channel node
                vew_node_map[bank_node] = channel_node
        
        self._vew_node_map = vew_node_map
        print(f"Found {len(vew_node_map)} VEW bank->channel node mappings")
        return self._vew_node_map

    def _get_flow_boundaries(self) -> Dict[str, List]:
        """Get all flow boundaries (IBTYPE ending with 2)."""
        boundaries = self._mesh.boundaries.to_dict()
        flow_boundaries = {}
        
        for ibtype, boundary_list in boundaries.items():
            if ibtype is not None and ibtype.endswith('2'):
                flow_boundaries[ibtype] = boundary_list
                
        return flow_boundaries

    def _find_corrupted_flow_boundaries(self) -> List[Tuple[str, int, int]]:
        """Find flow boundaries with corrupted nodes (bank nodes instead of channel nodes)."""
        corrupted_boundaries = []
        
        # Build VEW node mapping
        vew_node_map = self._build_vew_node_map()
        
        if not vew_node_map:
            print("No VEW boundaries found, no flow boundary corruption possible")
            return corrupted_boundaries
            
        # Get flow boundaries
        flow_boundaries = self._get_flow_boundaries()
        
        if not flow_boundaries:
            print("No flow boundaries found in mesh")
            return corrupted_boundaries
            
        print(f"Checking {len(flow_boundaries)} flow boundary types for corruption...")
        
        for ibtype, boundary_list in flow_boundaries.items():
            print(f"Checking IBTYPE {ibtype} for corruption...")
            for boundary_idx, boundary in enumerate(boundary_list):
                print(f"Checking boundary {boundary_idx} for corruption...")
                node_ids = boundary['node_id']
                
                for node_idx, node_id_str in enumerate(node_ids):
                    print(f"Checking node {node_idx} ({node_id_str}) for corruption...")
                    node_id = int(node_id_str)
                    
                    print(f"node_id: {node_id}")
                    print(f"this node is a bank node: {node_id in vew_node_map}")
                    print(f"this node is a channel node: {node_id in vew_node_map.values()}")
                    
                    # Check if this node is a bank node in VEW boundaries
                    if node_id in vew_node_map:
                        corrupted_boundaries.append((ibtype, boundary_idx, node_idx))
                        channel_node = vew_node_map[node_id]
                        print(f"Found corrupted flow boundary: IBTYPE {ibtype}, "
                              f"boundary {boundary_idx}, node index {node_idx}, "
                              f"bank node {node_id} should be channel node {channel_node}")
        
        return corrupted_boundaries

    def repair_flow_boundaries(self) -> AdcircMesh:
        """Repair corrupted flow boundaries in the mesh."""
        print("Starting flow boundary repair...")
        
        # Find corrupted boundaries
        corrupted_boundaries = self._find_corrupted_flow_boundaries()
        
        if not corrupted_boundaries:
            print("No corrupted flow boundaries found")
            return self._mesh
            
        print(f"Found {len(corrupted_boundaries)} corrupted flow boundary nodes")
        
        # Build VEW node mapping
        vew_node_map = self._build_vew_node_map()
        
        # Create a copy of the mesh boundaries
        boundaries = self._mesh.boundaries.to_dict().copy()
        
        # Track repairs made
        repairs_made = 0
        
        # Process each corrupted boundary
        for ibtype, boundary_idx, node_idx in corrupted_boundaries:
            boundary = boundaries[ibtype][boundary_idx]
            bank_node_str = boundary['node_id'][node_idx]
            bank_node = int(bank_node_str)
            channel_node = vew_node_map[bank_node]
            
            print(f"Repairing IBTYPE {ibtype}, boundary {boundary_idx}, "
                  f"node index {node_idx}: {bank_node} -> {channel_node}")
            
            # Update the boundary node
            boundary['node_id'][node_idx] = str(channel_node)
            repairs_made += 1
        
        print(f"Successfully repaired {repairs_made} flow boundary nodes")
        
        # Create new mesh with repaired boundaries
        repaired_mesh = AdcircMesh(
            nodes=self._mesh.nodes,
            elements=self._mesh.elements.elements,
            boundaries=boundaries
        )
        
        return repaired_mesh

    def validate_flow_boundaries(self) -> bool:
        """Validate that no flow boundary nodes are bank nodes from VEW boundaries."""
        # Build VEW node mapping
        vew_node_map = self._build_vew_node_map()
        
        if not vew_node_map:
            print("No VEW boundaries found, flow boundary validation passed")
            return True
            
        # Get flow boundaries
        flow_boundaries = self._get_flow_boundaries()
        
        if not flow_boundaries:
            print("No flow boundaries found, flow boundary validation passed")
            return True
            
        # Check for corrupted boundaries
        for ibtype, boundary_list in flow_boundaries.items():
            for boundary_idx, boundary in enumerate(boundary_list):
                node_ids = boundary['node_id']
                
                for node_idx, node_id_str in enumerate(node_ids):
                    node_id = int(node_id_str)
                    
                    if node_id in vew_node_map:
                        channel_node = vew_node_map[node_id]
                        print(f"Validation failed: IBTYPE {ibtype}, boundary {boundary_idx}, "
                              f"node index {node_idx} has bank node {node_id} "
                              f"instead of channel node {channel_node}")
                        return False
        
        print("Flow boundary validation passed")
        return True

    def get_repair_statistics(self) -> Dict:
        """Get statistics about flow boundaries and potential repairs."""
        stats = {
            'total_flow_boundary_types': 0,
            'total_flow_boundaries': 0,
            'total_flow_nodes': 0,
            'corrupted_nodes': 0,
            'vew_bank_nodes': 0
        }
        
        # Build VEW node mapping
        vew_node_map = self._build_vew_node_map()
        stats['vew_bank_nodes'] = len(vew_node_map)
        
        # Get flow boundaries
        flow_boundaries = self._get_flow_boundaries()
        stats['total_flow_boundary_types'] = len(flow_boundaries)
        
        if not flow_boundaries:
            return stats
            
        for ibtype, boundary_list in flow_boundaries.items():
            stats['total_flow_boundaries'] += len(boundary_list)
            
            for boundary in boundary_list:
                node_ids = boundary['node_id']
                stats['total_flow_nodes'] += len(node_ids)
                
                for node_id_str in node_ids:
                    node_id = int(node_id_str)
                    if node_id in vew_node_map:
                        stats['corrupted_nodes'] += 1
        
        return stats


def get_parser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        'input_mesh',
        help='Path to input mesh file with corrupted flow boundaries'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output mesh file path (default: input_mesh with _flow_repaired suffix)'
    )
    parser.add_argument(
        '--validate-only',
        action='store_true',
        help='Only validate flow boundaries without repairing'
    )
    parser.add_argument(
        '--statistics',
        action='store_true',
        help='Show statistics about flow boundaries'
    )
    parser.add_argument(
        '-d', '--description',
        default='flow_repaired',
        help='Description for the output mesh (default: flow_repaired)'
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
    repairer = FlowBoundaryRepairer(mesh)
    
    # Show statistics if requested
    if args.statistics:
        print("\nFlow Boundary Statistics:")
        stats = repairer.get_repair_statistics()
        for key, value in stats.items():
            print(f"  {key.replace('_', ' ').title()}: {value}")
        print()
    
    # Validate only if requested
    if args.validate_only:
        print("Validating flow boundaries...")
        is_valid = repairer.validate_flow_boundaries()
        if is_valid:
            print("All flow boundaries are valid")
            return 0
        else:
            print("Found corrupted flow boundaries")
            return 1
    
    # Repair the boundaries
    try:
        repaired_mesh = repairer.repair_flow_boundaries()
        
        # Validate the repaired mesh
        repairer_validation = FlowBoundaryRepairer(repaired_mesh)
        if not repairer_validation.validate_flow_boundaries():
            raise FlowCorruptionError("Repair process failed validation")
        
        # Determine output filename
        if args.output:
            output_path = args.output
        else:
            # Add _flow_repaired suffix to input filename
            input_parts = args.input_mesh.rsplit('.', 1)
            if len(input_parts) == 2:
                output_path = f"{input_parts[0]}_flow_repaired.{input_parts[1]}"
            else:
                output_path = f"{args.input_mesh}_flow_repaired"
        
        # Set description and write
        repaired_mesh.description = args.description
        repaired_mesh.write(output_path, overwrite=True)
        print(f"\nRepaired mesh saved to: {output_path}")
        
        return 0
        
    except FlowCorruptionError as e:
        print(f"Error: {e}")
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    exit(main()) 