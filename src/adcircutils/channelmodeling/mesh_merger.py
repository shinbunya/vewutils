#!/usr/bin/env python3
"""
Module for merging ADCIRC meshes with optional VEW boundary generation.
"""

import argparse
import pandas as pd
import numpy as np
from adcircpy import AdcircMesh
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Set


class MeshMerger:
    """Class for merging ADCIRC meshes with configurable VEW boundaries."""

    def __init__(self, channel_mesh: AdcircMesh, land_mesh: AdcircMesh):
        """
        Initialize the mesh merger with channel and land meshes.

        Args:
            channel_mesh: The channel mesh to be merged
            land_mesh: The land mesh to merge into
        """
        self._channel_mesh = channel_mesh
        self._land_mesh = land_mesh
        self._config = self._create_default_config()

    @staticmethod
    def _create_default_config() -> Dict:
        """Create default configuration for mesh merging."""
        return {
            'tolerance': 1e-6,
            'height_offset': 0.001,
            'subcritical_coefficient': 1.0,
            'supercritical_coefficient': 1.0,
            'vew_enabled': False
        }

    @staticmethod
    def _find_edges(elements_df: pd.DataFrame) -> Dict:
        """
        Find all edges in the mesh and count how many elements each edge belongs to.
        
        Args:
            elements_df: DataFrame containing element definitions
            
        Returns:
            Dictionary mapping edges (as frozenset of node IDs) to the number of elements they belong to
        """
        edge_counts = defaultdict(int)
        
        for _, element in elements_df.iterrows():
            nodes = element[1:4].astype(int).tolist()
            edges = [
                frozenset([nodes[0], nodes[1]]),
                frozenset([nodes[1], nodes[2]]),
                frozenset([nodes[2], nodes[0]])
            ]
            
            for edge in edges:
                edge_counts[edge] += 1
        
        return edge_counts

    @staticmethod
    def _split_node_string(nodes: List[int], 
                          node_neighbors: Dict[int, List[int]], 
                          elements_df: pd.DataFrame) -> List[List[int]]:
        """
        Split a string of nodes into segments based on domain boundaries and valid edges.
        
        Args:
            nodes: List of node IDs
            node_neighbors: Dictionary of node neighbors
            elements_df: DataFrame containing element definitions
            
        Returns:
            List of node segments
        """
        edge_counts = MeshMerger._find_edges(elements_df)
        
        segments = []
        i = 0
        current_segment = [i]
        processed = set([i])
        
        while len(processed) < len(nodes):
            found = False
            
            # Try to extend segment forward
            for jj in node_neighbors[nodes[i]]:
                if jj in nodes:
                    j = nodes.index(jj)
                    if j not in processed:
                        edge = frozenset([nodes[i], nodes[j]])
                        if edge in edge_counts and edge_counts[edge] == 1:
                            current_segment.append(j)
                            processed.add(j)
                            i = j
                            found = True
                            break

            # If can't extend forward, try extending backward
            if not found and len(current_segment) > 1:
                current_segment = list(reversed(current_segment))
                i = current_segment[-1]
                for jj in node_neighbors[nodes[i]]:
                    if jj in nodes:
                        j = nodes.index(jj)
                        if j not in processed:
                            edge = frozenset([nodes[i], nodes[j]])
                            if edge in edge_counts and edge_counts[edge] == 1:
                                current_segment.append(j)
                                processed.add(j)
                                i = j
                                found = True
                                break

            # If can't extend in either direction, start new segment
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
        
        # Add the last segment if it has more than one node
        if len(current_segment) > 1:
            segments.append(current_segment)
        
        return segments

    def _find_paired_nodes(self, nodes_df: pd.DataFrame) -> List[Tuple[int, int]]:
        """
        Find pairs of nodes that are at identical locations.
        
        Args:
            nodes_df: DataFrame containing node coordinates
            
        Returns:
            List of tuples containing paired node IDs
        """
        paired_nodes = []
        processed: Set[int] = set()
        
        coords = nodes_df[['x', 'y']].values
        
        for i in range(len(coords)):
            if i in processed:
                continue
                
            distances = np.sqrt(np.sum((coords - coords[i])**2, axis=1))
            matches = np.where(distances < self._config['tolerance'])[0]
            
            if len(matches) > 1:
                for j in matches:
                    if j != i and j not in processed:
                        paired_nodes.append((i+1, j+1))
                        processed.add(j)
                processed.add(i)
        
        return paired_nodes

    def _create_vew_boundaries(self, paired_nodes: List[Tuple[int, int]], 
                             node_mapping: Dict[int, int],
                             combined_elements: pd.DataFrame) -> List[Dict]:
        """Create VEW boundary definitions."""
        nodes_land = [node for node, _ in paired_nodes]
        nodes_channel = [node for _, node in paired_nodes]
        
        segments = self._split_node_string(
            nodes_land, 
            self._land_mesh.node_neighbors, 
            combined_elements
        )
        
        vewboundaries = []
        for segment in segments:
            segment_pairs = [(nodes_land[i], nodes_channel[i]) for i in segment]
            vew_node_id = [(str(node1), str(node2)) for node1, node2 in segment_pairs]
            
            vew_barrier_height = []
            for node1, _ in segment_pairs:
                elevation = (self._land_mesh.nodes.loc[node1, 'value_1'] + 
                           self._config['height_offset'])
                vew_barrier_height.append(elevation)
            
            vewboundary = {
                'node_id': vew_node_id,
                'barrier_height': vew_barrier_height,
                'subcritical_flow_coefficient': [self._config['subcritical_coefficient']] * len(segment_pairs),
                'supercritical_flow_coefficient': [self._config['supercritical_coefficient']] * len(segment_pairs)
            }
            
            vewboundaries.append(vewboundary)
            
        return vewboundaries

    # Builder methods
    def with_tolerance(self, tolerance: float) -> 'MeshMerger':
        """Set the tolerance for node matching."""
        self._config['tolerance'] = tolerance
        return self

    def with_height_offset(self, offset: float) -> 'MeshMerger':
        """Set the height offset for VEW boundaries."""
        self._config['height_offset'] = offset
        return self

    def with_flow_coefficients(self, subcritical: float, supercritical: float) -> 'MeshMerger':
        """Set the flow coefficients for VEW boundaries."""
        self._config['subcritical_coefficient'] = subcritical
        self._config['supercritical_coefficient'] = supercritical
        return self

    def with_vew_boundaries(self, enabled: bool = True) -> 'MeshMerger':
        """Enable or disable VEW boundaries."""
        self._config['vew_enabled'] = enabled
        return self

    def merge(self) -> AdcircMesh:
        """
        Merge the meshes and add VEW boundaries if enabled.
        
        Returns:
            AdcircMesh: The merged mesh
        """
        # Get the number of nodes in the land mesh
        land_node_count = len(self._land_mesh.nodes)
        
        # Create node mapping
        node_mapping = {old_id: new_id for old_id, new_id in 
                       zip(self._channel_mesh.nodes.index, 
                           range(land_node_count + 1, 
                                land_node_count + len(self._channel_mesh.nodes) + 1))}
        
        # Combine nodes
        channel_nodes = self._channel_mesh.nodes.copy()
        channel_nodes.index = [node_mapping[i] for i in channel_nodes.index]
        combined_nodes = pd.concat([self._land_mesh.nodes, channel_nodes])
        
        # Combine elements
        channel_elements = self._channel_mesh.elements.elements.copy()
        for col in channel_elements.columns[1:]:
            channel_elements[col] = channel_elements[col].map(node_mapping)
        combined_elements = pd.concat([self._land_mesh.elements.elements, channel_elements])
        
        # Create the merged mesh
        merged_mesh = AdcircMesh(nodes=combined_nodes, elements=combined_elements)
        
        # Add VEW boundaries if enabled
        if self._config['vew_enabled']:
            paired_nodes = self._find_paired_nodes(merged_mesh.nodes)
            
            if paired_nodes:
                vewboundaries = self._create_vew_boundaries(
                    paired_nodes, 
                    node_mapping,
                    combined_elements
                )
                
                boundaries = merged_mesh.boundaries.to_dict()
                if '64' in boundaries:
                    boundaries['64'].extend(vewboundaries)
                else:
                    boundaries['64'] = vewboundaries
                
                merged_mesh = AdcircMesh(
                    nodes=merged_mesh.nodes,
                    elements=merged_mesh.elements.elements,
                    boundaries=boundaries
                )
                
                print(f"Added {len(vewboundaries)} VEW boundary segments")
        
        return merged_mesh


def main():
    """Main function to handle command line arguments and process the meshes."""
    parser = argparse.ArgumentParser(
        description="Merge ADCIRC meshes with optional VEW boundary generation"
    )
    parser.add_argument(
        "channelmesh",
        help="Path to the channel mesh file"
    )
    parser.add_argument(
        "landmesh",
        help="Path to the land mesh file"
    )
    parser.add_argument(
        "output",
        help="Path to save the merged mesh"
    )
    parser.add_argument(
        "--enable-vew",
        action="store_true",
        help="Enable VEW boundary generation between channel and land meshes"
    )
    
    args = parser.parse_args()
    
    try:
        # Read the mesh files
        channel_mesh = AdcircMesh.open(args.channelmesh)
        land_mesh = AdcircMesh.open(args.landmesh)
        
        print(f"Successfully read channel mesh from: {args.channelmesh}")
        print(f"Successfully read land mesh from: {args.landmesh}")
        
        # Print mesh information
        print(f"\nChannel mesh info:")
        print(f"Number of nodes: {len(channel_mesh.nodes)}")
        print(f"Number of elements: {len(channel_mesh.elements.elements)}")
        
        print(f"\nLand mesh info:")
        print(f"Number of nodes: {len(land_mesh.nodes)}")
        print(f"Number of elements: {len(land_mesh.elements.elements)}")
        
        # Merge the meshes
        merger = (MeshMerger(channel_mesh, land_mesh)
                 .with_vew_boundaries(args.enable_vew)
                 .with_tolerance(1e-6)
                 .with_height_offset(0.001)
                 .with_flow_coefficients(1.0, 1.0))
        
        merged_mesh = merger.merge()
        
        print(f"\nMerged mesh info:")
        print(f"Number of nodes: {len(merged_mesh.nodes)}")
        print(f"Number of elements: {len(merged_mesh.elements.elements)}")
        
        # Save the merged mesh
        merged_mesh.write(args.output, overwrite=True)
        print(f"\nMerged mesh saved to: {args.output}")
        
    except Exception as e:
        print(f"Error processing meshes: {str(e)}")
        return 1
    
    return 0


if __name__ == "__main__":
    main()
