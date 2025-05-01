#!/usr/bin/env python3
"""
Module for adding land boundaries (ibtype = 20) along unassigned boundary segments.
"""

import argparse
import numpy as np
from collections import defaultdict
from typing import Dict, List, Set, Tuple
from adcircpy import AdcircMesh


class LandBoundaryAdder:
    """Class for adding land boundaries along unassigned boundary segments."""

    def __init__(self, mesh: AdcircMesh):
        """Initialize with an ADCIRC mesh."""
        self._mesh = mesh
        # Compute edge counts, and pre-compute boundary edges and nodes
        print("Computing mesh topology...")
        self._edge_counts = self._find_edges(self._mesh.elements.elements)
        self._boundary_edges = {edge for edge, count in self._edge_counts.items() if count == 1}
        self._boundary_nodes = set()
        for edge in self._boundary_edges:
            self._boundary_nodes.update(edge)
        print(f"Found {len(self._boundary_edges)} boundary edges and {len(self._boundary_nodes)} boundary nodes")

    @staticmethod
    def _find_edges(elements_df) -> Dict:
        """
        Find all edges in the mesh and count how many elements each edge belongs to.
        
        Args:
            elements_df: DataFrame containing element definitions
            
        Returns:
            Dictionary mapping edges (as frozenset of node IDs) to the number of elements they belong to
        """
        edge_counts = defaultdict(int)
        
        # Use list comprehension for speed instead of loops
        nodes_array = elements_df.iloc[:, 1:4].astype(int).values
        for nodes in nodes_array:
            edges = [
                frozenset([nodes[0], nodes[1]]),
                frozenset([nodes[1], nodes[2]]),
                frozenset([nodes[2], nodes[0]])
            ]
            for edge in edges:
                edge_counts[edge] += 1
        
        return edge_counts

    def _get_boundary_nodes(self) -> Set[int]:
        """Get all nodes that are part of boundary edges."""
        return self._boundary_nodes

    def _get_assigned_boundary_nodes(self) -> Tuple[Set[int], Set[int]]:
        """
        Get all nodes that are already part of boundary conditions and their endpoints.
        
        Returns:
            Tuple of (set of all assigned nodes, set of endpoint nodes that are part of unassigned boundary edges)
        """
        assigned_nodes = set()
        assigned_edges = set()
        
        # Build node-to-assigned-edges mapping for faster lookups
        node_to_assigned_edges = defaultdict(set)
        
        # First pass: collect all assigned nodes and edges
        for ibtype, boundaries in self._mesh.boundaries.to_dict().items():
            for boundary in boundaries:
                node_ids = boundary['node_id']
                if isinstance(node_ids[0], tuple):  # Weir boundary
                    # For weir boundaries, we have two separate node strings
                    # First string: all first nodes of each pair
                    first_nodes = [int(pair[0]) for pair in node_ids]
                    assigned_nodes.update(first_nodes)
                    # Add edges between consecutive first nodes
                    for i in range(len(first_nodes) - 1):
                        edge = frozenset([first_nodes[i], first_nodes[i + 1]])
                        assigned_edges.add(edge)
                        node_to_assigned_edges[first_nodes[i]].add(edge)
                        node_to_assigned_edges[first_nodes[i + 1]].add(edge)
                    
                    # Second string: all second nodes of each pair
                    second_nodes = [int(pair[1]) for pair in node_ids]
                    assigned_nodes.update(second_nodes)
                    # Add edges between consecutive second nodes
                    for i in range(len(second_nodes) - 1):
                        edge = frozenset([second_nodes[i], second_nodes[i + 1]])
                        assigned_edges.add(edge)
                        node_to_assigned_edges[second_nodes[i]].add(edge)
                        node_to_assigned_edges[second_nodes[i + 1]].add(edge)
                else:  # Other boundary types
                    node_ids = [int(node_id) for node_id in node_ids]
                    assigned_nodes.update(node_ids)
                    # Add edges between consecutive nodes
                    for i in range(len(node_ids) - 1):
                        edge = frozenset([node_ids[i], node_ids[i + 1]])
                        assigned_edges.add(edge)
                        node_to_assigned_edges[node_ids[i]].add(edge)
                        node_to_assigned_edges[node_ids[i + 1]].add(edge)
        
        # Find unassigned boundary edges and their endpoints
        endpoint_nodes = set()
        unassigned_boundary_edges = self._boundary_edges - assigned_edges
        
        for edge in unassigned_boundary_edges:
            node1, node2 = edge
            # Check if both nodes have assigned edges (with pre-computed lookup)
            if not (node_to_assigned_edges[node1] and node_to_assigned_edges[node2]):
                endpoint_nodes.update(edge)
        
        return assigned_nodes, endpoint_nodes

    def _split_node_string(self, nodes: List[int], node_neighbors: Dict[int, List[int]]) -> List[List[int]]:
        """Split a string of nodes into segments based on domain boundaries."""
        # Create lookup dictionaries for faster access
        node_to_idx = {node: idx for idx, node in enumerate(nodes)}
        
        segments = []
        i = 0  # Start with the first node as in the original algorithm
        current_segment = [i]
        processed = set([i])
        
        while len(processed) < len(nodes):
            found = False
            current_node = nodes[i]
            
            # Try to extend segment forward (same as original algorithm)
            for neighbor in node_neighbors[current_node]:
                if neighbor in node_to_idx:  # Much faster than nodes.index()
                    j = node_to_idx[neighbor]
                    if j not in processed:
                        edge = frozenset([current_node, neighbor])
                        if edge in self._boundary_edges:  # Using pre-computed boundary edges
                            current_segment.append(j)
                            processed.add(j)
                            i = j
                            found = True
                            break

            # If can't extend forward, try extending backward (same as original algorithm)
            if not found and len(current_segment) > 1:
                current_segment = list(reversed(current_segment))
                i = current_segment[-1]
                current_node = nodes[i]
                for neighbor in node_neighbors[current_node]:
                    if neighbor in node_to_idx:  # Much faster than nodes.index()
                        j = node_to_idx[neighbor]
                        if j not in processed:
                            edge = frozenset([current_node, neighbor])
                            if edge in self._boundary_edges:  # Using pre-computed boundary edges
                                current_segment.append(j)
                                processed.add(j)
                                i = j
                                found = True
                                break

            # If can't extend in either direction, start new segment (same as original algorithm)
            if not found:
                if len(current_segment) > 1:
                    segments.append(current_segment)
                
                # Find next unprocessed node
                next_i = None
                for idx in range(len(nodes)):
                    if idx not in processed:
                        next_i = idx
                        break
                        
                if next_i is None:
                    break
                    
                i = next_i
                current_segment = [i]
                processed.add(i)
        
        # Add the last segment if it has more than one node (same as original algorithm)
        if len(current_segment) > 1:
            segments.append(current_segment)
        
        return segments

    def add_land_boundaries(self) -> AdcircMesh:
        """Add land boundaries (ibtype = 20) along unassigned boundary segments."""
        print("\nStarting to add land boundaries...")
        
        # Get all boundary nodes and those already assigned to boundaries
        print("Finding boundary nodes...")
        boundary_nodes = self._get_boundary_nodes()
        print(f"Found {len(boundary_nodes)} boundary nodes")
        
        print("Finding nodes with existing boundary conditions...")
        assigned_nodes, endpoint_nodes = self._get_assigned_boundary_nodes()
        print(f"Found {len(assigned_nodes)} nodes with existing boundary conditions")
        print(f"Found {len(endpoint_nodes)} endpoints of unassigned boundary edges")
        
        # Find unassigned boundary nodes, including endpoints of unassigned boundary edges
        unassigned_nodes = list((boundary_nodes - assigned_nodes) | endpoint_nodes)
        print(f"Found {len(unassigned_nodes)} unassigned boundary nodes (including relevant endpoints)")
        
        if not unassigned_nodes:
            print("No unassigned boundary nodes found - no land boundaries to add")
            return self._mesh
        
        # Split unassigned nodes into continuous segments
        print("Splitting nodes into continuous segments...")
        segments = self._split_node_string(unassigned_nodes, self._mesh.node_neighbors)
        print(f"Split into {len(segments)} continuous segments")
        
        # Create land boundaries for each segment
        print("Creating land boundary definitions...")
        land_boundaries = []
        for i, segment in enumerate(segments, 1):
            node_ids = [str(unassigned_nodes[j]) for j in segment]
            land_boundary = {
                'node_id': node_ids
            }
            land_boundaries.append(land_boundary)
            print(f"  Created land boundary {i} with {len(node_ids)} nodes")
        
        # Add land boundaries to mesh
        print("Adding land boundaries to mesh...")
        boundaries = self._mesh.boundaries.to_dict()
        if '20' in boundaries:
            print("  Appending to existing land boundaries...")
            boundaries['20'].extend(land_boundaries)
        else:
            print("  Creating new land boundaries section...")
            boundaries['20'] = land_boundaries
        
        # Create new mesh with added boundaries
        print("Creating new mesh with added boundaries...")
        new_mesh = AdcircMesh(
            nodes=self._mesh.nodes,
            elements=self._mesh.elements.elements,
            boundaries=boundaries
        )
        
        print(f"\nSuccessfully added {len(land_boundaries)} land boundary segments")
        print("Land boundary addition complete!")
        
        return new_mesh


def main():
    """Main function to handle command line arguments and process the mesh."""
    parser = argparse.ArgumentParser(
        description="Add land boundaries (ibtype = 20) along unassigned boundary segments"
    )
    parser.add_argument(
        "input_mesh",
        help="Path to the input mesh file"
    )
    parser.add_argument(
        '-o', '--output',
        default='mesh_with_land_boundaries.grd',
        help='Path to save the output mesh (default: mesh_with_land_boundaries.grd)'
    )
    parser.add_argument(
        '-d', '--description',
        default='mesh with land boundaries',
        help='Description for the output mesh (default: mesh with land boundaries)'
    )
    
    args = parser.parse_args()
    
    # Read the mesh file
    print("Reading input mesh...")
    mesh = AdcircMesh.open(args.input_mesh)
    print(f"Successfully read mesh from: {args.input_mesh}")
    
    # Print mesh information
    print("\nInput mesh info:")
    print(f"Number of nodes: {len(mesh.nodes)}")
    print(f"Number of elements: {len(mesh.elements.elements)}")
    
    # Add land boundaries
    adder = LandBoundaryAdder(mesh)
    new_mesh = adder.add_land_boundaries()
    
    print("\nOutput mesh info:")
    print(f"Number of nodes: {len(new_mesh.nodes)}")
    print(f"Number of elements: {len(new_mesh.elements.elements)}")
    
    # Save the new mesh
    print("\nSaving output mesh...")
    new_mesh.description = args.description
    new_mesh.write(args.output, overwrite=True)
    print(f"Successfully saved to: {args.output}")
    
    return 0


if __name__ == "__main__":
    main() 