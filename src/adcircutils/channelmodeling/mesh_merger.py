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
from abc import ABC, abstractmethod


class MergeStrategy(ABC):
    """Abstract base class for mesh merging strategies."""
    
    @abstractmethod
    def merge(self, channel_mesh: AdcircMesh, land_mesh: AdcircMesh, config: Dict) -> AdcircMesh:
        """Merge two meshes according to the strategy."""
        pass


class VEWBoundaryStrategy(MergeStrategy):
    """Strategy that keeps duplicate nodes and adds VEW boundaries."""
    
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
        """Split a string of nodes into segments based on domain boundaries."""
        edge_counts = VEWBoundaryStrategy._find_edges(elements_df)
        
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

    def _find_paired_nodes(self, nodes_df: pd.DataFrame, tolerance: float) -> List[Tuple[int, int]]:
        """Find pairs of nodes that are at identical locations."""
        paired_nodes = []
        processed: Set[int] = set()
        
        coords = nodes_df[['x', 'y']].values
        
        for i in range(len(coords)):
            if i in processed:
                continue
                
            distances = np.sqrt(np.sum((coords - coords[i])**2, axis=1))
            matches = np.where(distances < tolerance)[0]
            
            if len(matches) > 1:
                for j in matches:
                    if j != i and j not in processed:
                        paired_nodes.append((i+1, j+1))
                        processed.add(j)
                processed.add(i)
        
        return paired_nodes

    def _create_vew_boundaries(self, 
                             paired_nodes: List[Tuple[int, int]], 
                             land_mesh: AdcircMesh,
                             node_mapping: Dict[int, int],
                             combined_elements: pd.DataFrame,
                             config: Dict) -> List[Dict]:
        """Create VEW boundary definitions."""
        nodes_land = [node for node, _ in paired_nodes]
        nodes_channel = [node for _, node in paired_nodes]
        
        segments = self._split_node_string(
            nodes_land, 
            land_mesh.node_neighbors, 
            combined_elements
        )
        
        vewboundaries = []
        for segment in segments:
            segment_pairs = [(nodes_land[i], nodes_channel[i]) for i in segment]
            vew_node_id = [(str(node1), str(node2)) for node1, node2 in segment_pairs]
            
            vew_barrier_height = []
            for node1, _ in segment_pairs:
                elevation = (land_mesh.nodes.loc[node1, 'value_1'] + 
                           config['height_offset'])
                vew_barrier_height.append(elevation)
            
            vewboundary = {
                'node_id': vew_node_id,
                'barrier_height': vew_barrier_height,
                'subcritical_flow_coefficient': [config['subcritical_coefficient']] * len(segment_pairs),
                'supercritical_flow_coefficient': [config['supercritical_coefficient']] * len(segment_pairs)
            }
            
            vewboundaries.append(vewboundary)
            
        return vewboundaries

    def merge(self, channel_mesh: AdcircMesh, land_mesh: AdcircMesh, config: Dict) -> AdcircMesh:
        """Merge meshes using VEW boundaries."""
        # Get the number of nodes in the land mesh
        land_node_count = len(land_mesh.nodes)
        
        # Create node mapping
        node_mapping = {old_id: new_id for old_id, new_id in 
                       zip(channel_mesh.nodes.index, 
                           range(land_node_count + 1, 
                                land_node_count + len(channel_mesh.nodes) + 1))}
        
        # Combine nodes
        channel_nodes = channel_mesh.nodes.copy()
        channel_nodes.index = [node_mapping[i] for i in channel_nodes.index]
        combined_nodes = pd.concat([land_mesh.nodes, channel_nodes])
        
        # Combine elements
        channel_elements = channel_mesh.elements.elements.copy()
        for col in channel_elements.columns[1:]:
            channel_elements[col] = channel_elements[col].map(node_mapping)
        combined_elements = pd.concat([land_mesh.elements.elements, channel_elements])
        
        # Create the merged mesh
        merged_mesh = AdcircMesh(nodes=combined_nodes, elements=combined_elements)
        
        # Add VEW boundaries
        paired_nodes = self._find_paired_nodes(merged_mesh.nodes, config['tolerance'])
        
        if paired_nodes:
            vewboundaries = self._create_vew_boundaries(
                paired_nodes, 
                land_mesh,
                node_mapping,
                combined_elements,
                config
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


class ConsolidatedNodesStrategy(MergeStrategy):
    """Strategy that consolidates duplicate nodes at boundaries."""
    
    def _find_matching_nodes(self, nodes_df1: pd.DataFrame, nodes_df2: pd.DataFrame, 
                           tolerance: float) -> Dict[int, int]:
        """
        Find matching nodes between two meshes within tolerance.
        
        Args:
            nodes_df1: Channel mesh nodes DataFrame
            nodes_df2: Land mesh nodes DataFrame
            tolerance: Distance tolerance for matching nodes
            
        Returns:
            Dictionary mapping channel mesh node IDs to land mesh node IDs
        """
        coords1 = nodes_df1[['x', 'y']].values
        coords2 = nodes_df2[['x', 'y']].values
        matching_nodes = {}
        
        for i, coord1 in enumerate(coords1):
            distances = np.sqrt(np.sum((coords2 - coord1)**2, axis=1))
            matches = np.where(distances < tolerance)[0]
            if len(matches) > 0:
                # Match to the closest node
                closest = matches[np.argmin(distances[matches])]
                matching_nodes[i + 1] = closest + 1  # Convert to 1-based indexing
        
        return matching_nodes

    def merge(self, channel_mesh: AdcircMesh, land_mesh: AdcircMesh, config: Dict) -> AdcircMesh:
        """Merge meshes by consolidating duplicate nodes."""
        # Find matching nodes
        matching_nodes = self._find_matching_nodes(
            channel_mesh.nodes, 
            land_mesh.nodes, 
            config['tolerance']
        )
        
        # Create node mapping for non-matching nodes
        land_node_count = len(land_mesh.nodes)
        node_mapping = {}
        next_node_id = land_node_count + 1
        
        for old_id in channel_mesh.nodes.index:
            if old_id in matching_nodes:
                # Use existing node ID from land mesh
                node_mapping[old_id] = matching_nodes[old_id]
            else:
                # Assign new node ID
                node_mapping[old_id] = next_node_id
                next_node_id += 1
        
        # Start with a copy of the land mesh nodes
        combined_nodes = land_mesh.nodes.copy()
        
        # Update values for matching nodes with channel mesh values
        for channel_id, land_id in matching_nodes.items():
            combined_nodes.loc[land_id, 'value_1'] = channel_mesh.nodes.loc[channel_id, 'value_1']
        
        # Add non-matching nodes from channel mesh
        new_nodes = []
        for old_id in channel_mesh.nodes.index:
            if old_id not in matching_nodes:
                node = channel_mesh.nodes.loc[old_id].copy()
                node.name = node_mapping[old_id]
                new_nodes.append(node)
        
        if new_nodes:
            new_nodes_df = pd.DataFrame(new_nodes)
            combined_nodes = pd.concat([combined_nodes, new_nodes_df])
        
        # Update element connectivity
        channel_elements = channel_mesh.elements.elements.copy()
        for col in channel_elements.columns[1:]:
            channel_elements[col] = channel_elements[col].map(node_mapping)
        
        combined_elements = pd.concat([land_mesh.elements.elements, channel_elements])
        
        merged_mesh = AdcircMesh(nodes=combined_nodes, elements=combined_elements)
        print(f"Consolidated {len(matching_nodes)} duplicate nodes")
        print(f"Updated elevations at {len(matching_nodes)} matching nodes with channel mesh values")
        
        return merged_mesh


class MeshMerger:
    """Class for merging ADCIRC meshes with configurable strategies."""

    def __init__(self, channel_mesh: AdcircMesh, land_mesh: AdcircMesh):
        """Initialize the mesh merger with channel and land meshes."""
        self._channel_mesh = channel_mesh
        self._land_mesh = land_mesh
        self._config = self._create_default_config()
        self._strategy = None

    @staticmethod
    def _create_default_config() -> Dict:
        """Create default configuration for mesh merging."""
        return {
            'tolerance': 1e-6,
            'height_offset': 0.001,
            'subcritical_coefficient': 1.0,
            'supercritical_coefficient': 1.0,
            'vew_enabled': False,
            'consolidate_nodes': False
        }

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
        self._config['consolidate_nodes'] = not enabled
        return self

    def with_consolidated_nodes(self, enabled: bool = True) -> 'MeshMerger':
        """Enable or disable node consolidation."""
        self._config['consolidate_nodes'] = enabled
        self._config['vew_enabled'] = not enabled
        return self

    def merge(self) -> AdcircMesh:
        """Merge the meshes using the selected strategy."""
        if self._config['vew_enabled']:
            strategy = VEWBoundaryStrategy()
        else:
            strategy = ConsolidatedNodesStrategy()
            
        return strategy.merge(self._channel_mesh, self._land_mesh, self._config)


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
    parser.add_argument(
        "--consolidate-nodes",
        action="store_true",
        help="Consolidate duplicate nodes at mesh boundaries"
    )
    
    args = parser.parse_args()
    
    if args.enable_vew and args.consolidate_nodes:
        print("Error: Cannot use both --enable-vew and --consolidate-nodes together")
        return 1
    
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
                 .with_tolerance(1e-6)
                 .with_height_offset(0.001)
                 .with_flow_coefficients(1.0, 1.0))
        
        if args.enable_vew:
            merger.with_vew_boundaries(True)
        elif args.consolidate_nodes:
            merger.with_consolidated_nodes(True)
        
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