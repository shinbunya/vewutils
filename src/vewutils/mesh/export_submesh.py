#!/usr/bin/env python3
"""
Module for exporting detached submeshes from ADCIRC meshes.
"""

import argparse
import pandas as pd
import numpy as np
from adcircpy import AdcircMesh
from typing import Set


class SubmeshExporter:
    """Class for exporting detached submeshes from ADCIRC meshes."""

    def __init__(self, mesh: AdcircMesh):
        """Initialize with an ADCIRC mesh."""
        self._mesh = mesh

    def _find_connected_elements(self, start_element_id: int) -> Set[int]:
        """
        Find all elements in the connected component containing the given element.
        
        Args:
            start_element_id: Element ID to start from (1-based)
            
        Returns:
            Set of element IDs in the connected component
        """
        # Validate that the element exists
        if start_element_id not in self._mesh.elements.elements.index:
            raise ValueError(f"Element ID {start_element_id} not found in mesh")
        
        # Use BFS to find all connected elements
        connected_elements = set()
        queue = [start_element_id]
        visited = set()
        
        # Get node_elements mapping (node_id -> set of element_ids)
        node_elements = self._mesh.node_elements
        
        while queue:
            current_element_id = queue.pop(0)
            
            if current_element_id in visited:
                continue
                
            visited.add(current_element_id)
            connected_elements.add(current_element_id)
            
            # Get nodes of current element
            element_row = self._mesh.elements.elements.loc[current_element_id]
            nodes = [
                int(element_row['node_1']),
                int(element_row['node_2']),
                int(element_row['node_3'])
            ]
            
            # Find all elements connected through shared nodes
            for node_id in nodes:
                if node_id in node_elements:
                    # Get all elements containing this node
                    # Convert to int to match element index type
                    connected_element_ids = [int(elem_id) for elem_id in node_elements[node_id]]
                    for elem_id in connected_element_ids:
                        if elem_id not in visited:
                            queue.append(elem_id)
        
        return connected_elements

    def export_submesh(self, element_id: int) -> AdcircMesh:
        """
        Export a submesh containing the specified element and all connected elements.
        
        Args:
            element_id: Element ID to include in the submesh (1-based)
            
        Returns:
            AdcircMesh containing only the connected component
        """
        print(f"Finding connected elements starting from element {element_id}...")
        connected_elements = self._find_connected_elements(element_id)
        print(f"Found {len(connected_elements)} connected elements")
        
        # Extract elements DataFrame
        print("Extracting elements...")
        elements_subset = self._mesh.elements.elements.loc[list(connected_elements)].copy()
        
        # Find unique nodes used by these elements
        print("Finding unique nodes...")
        node_cols = ['node_1', 'node_2', 'node_3']
        unique_nodes = np.unique(elements_subset[node_cols].values)
        unique_nodes = sorted(unique_nodes.tolist())
        print(f"Found {len(unique_nodes)} unique nodes")
        
        # Create node mapping: old_node_id -> new_node_id (1-based consecutive)
        print("Creating node mapping...")
        node_mapping = {old_id: new_id for new_id, old_id in enumerate(unique_nodes, start=1)}
        
        # Create new nodes DataFrame with renumbered indices
        print("Creating new nodes DataFrame...")
        nodes_subset = self._mesh.nodes.loc[unique_nodes].copy()
        nodes_subset.index = range(1, len(unique_nodes) + 1)
        
        # Renumber element node references
        print("Renumbering element node references...")
        for col in node_cols:
            elements_subset[col] = elements_subset[col].map(node_mapping)
        
        # Reset element indices to be 1-based consecutive
        print("Renumbering element indices...")
        elements_subset.index = range(1, len(elements_subset) + 1)
        
        # Create new mesh with empty boundaries
        print("Creating submesh...")
        submesh = AdcircMesh(
            nodes=nodes_subset,
            elements=elements_subset,
            boundaries={}
        )
        
        print(f"Submesh created: {len(nodes_subset)} nodes, {len(elements_subset)} elements")
        return submesh


def get_parser():
    """Create argument parser for export-submesh command."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        'input_mesh_file',
        help='Path to input mesh file'
    )
    parser.add_argument(
        '--include-element',
        type=int,
        required=True,
        help='Element ID to include in the submesh (1-based)'
    )
    parser.add_argument(
        '-o', '--output',
        default='submesh.grd',
        help='Output mesh file path (default: submesh.grd)'
    )
    parser.add_argument(
        '-d', '--description',
        default='submesh',
        help='Description for the output mesh (default: submesh)'
    )
    return parser


def main(args=None):
    """Main function for export-submesh command."""
    if args is None:
        args = get_parser().parse_args()
    
    # Read the mesh file
    print(f"Reading input mesh from: {args.input_mesh_file}")
    mesh = AdcircMesh.open(args.input_mesh_file)
    print(f"Successfully read mesh")
    print(f"Number of nodes: {len(mesh.nodes)}")
    print(f"Number of elements: {len(mesh.elements.elements)}")
    
    # Validate element_id exists
    if args.include_element not in mesh.elements.elements.index:
        print(f"Error: Element ID {args.include_element} not found in mesh")
        print(f"Valid element IDs range from {mesh.elements.elements.index.min()} to {mesh.elements.elements.index.max()}")
        return 1
    
    # Create exporter and export submesh
    print(f"\nExporting submesh containing element {args.include_element}...")
    exporter = SubmeshExporter(mesh)
    submesh = exporter.export_submesh(args.include_element)
    
    # Set description and write mesh
    submesh.description = args.description
    submesh.write(args.output, overwrite=True)
    print(f"\nSubmesh saved to: {args.output}")
    
    return 0


if __name__ == "__main__":
    main()

