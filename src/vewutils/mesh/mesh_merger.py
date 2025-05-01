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
from scipy.spatial import cKDTree
from vewutils.mesh.vew_boundary_manipulator import VEWBoundaryManipulator


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
                # Before starting a new segment, check if current segment forms an island
                if len(current_segment) > 2:
                    first_node = nodes[current_segment[0]]
                    last_node = nodes[current_segment[-1]]
                    edge = frozenset([first_node, last_node])
                    if edge in edge_counts and edge_counts[edge] == 1:
                        # This is an island - add the first node to close it
                        current_segment.append(current_segment[0])
                
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
        
        # Check the last segment for island closure
        if len(current_segment) > 2:
            first_node = nodes[current_segment[0]]
            last_node = nodes[current_segment[-1]]
            edge = frozenset([first_node, last_node])
            if edge in edge_counts and edge_counts[edge] == 1:
                # This is an island - add the first node to close it
                current_segment.append(current_segment[0])
        
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

    def _filter_boundary_nodes(self, node_ids: List[str], original_length: int, matching_nodes: Set[str]) -> Tuple[List[int], bool]:
        """Filter out nodes that are members of matching nodes from a boundary.
        
        Args:
            node_ids: List of node IDs in the boundary
            original_length: Original length of the boundary
            matching_nodes: Set of node IDs that match between the two meshes
            
        Returns:
            Tuple of (list of indices to keep, whether to keep the boundary)
        """
        # Find indices of nodes to keep (not in matching_nodes)
        keep_indices = [i for i, node_id in enumerate(node_ids) if node_id not in matching_nodes]
        
        # If all nodes are removed, don't keep the boundary
        if not keep_indices:
            return [], False
            
        # Add one node before and after the filtered section if possible
        if keep_indices[0] > 0:
            keep_indices.insert(0, keep_indices[0]-1)
        if keep_indices[-1] < original_length - 1:
            keep_indices.append(keep_indices[-1]+1)
            
        # If only one node remains, don't keep the boundary
        if len(keep_indices) == 1:
            return [], False
            
        # If the boundary is reduced to just its endpoints, don't keep it
        if len(keep_indices) == 2 and keep_indices[0] == 0 and keep_indices[-1] == original_length - 1:
            return [], False
            
        return keep_indices, True

    def merge(self, channel_mesh: AdcircMesh, land_mesh: AdcircMesh, config: Dict) -> AdcircMesh:
        """Merge meshes using VEW boundaries."""
        print("\nStarting mesh merge with VEW boundaries...")
        
        # Get the number of nodes in the land mesh
        land_node_count = len(land_mesh.nodes)
        
        # Create node mapping
        print("Creating node mapping...")
        node_mapping = {old_id: new_id for old_id, new_id in 
                       zip(channel_mesh.nodes.index, 
                           range(land_node_count + 1, 
                                land_node_count + len(channel_mesh.nodes) + 1))}
        
        # Combine nodes
        print("Combining nodes...")
        channel_nodes = channel_mesh.nodes.copy()
        channel_nodes.index = [node_mapping[i] for i in channel_nodes.index]
        combined_nodes = pd.concat([land_mesh.nodes, channel_nodes])
        print(f"Combined nodes: {len(combined_nodes)}")
        
        # Combine elements
        print("Combining elements...")
        channel_elements = channel_mesh.elements.elements.copy()
        for col in channel_elements.columns[1:]:
            channel_elements[col] = channel_elements[col].map(node_mapping)
        combined_elements = pd.concat([land_mesh.elements.elements, channel_elements])
        # Reset element IDs to be 1-based
        combined_elements.index = range(1, len(combined_elements) + 1)
        print(f"Combined elements: {len(combined_elements)}")
        
        # Find matching nodes
        manipulator = VEWBoundaryManipulator()
        matching_nodes = manipulator.find_matching_nodes(
            channel_mesh.nodes, 
            land_mesh.nodes, 
            config['tolerance']
        )
        
        # Update boundaries of channel mesh
        print("Updating boundary node numbers...")
        boundaries_new = channel_mesh.boundaries.to_dict().copy()
        for ibtype in boundaries_new.keys():
            boundaries_ibtype_new = boundaries_new[ibtype]
            
            for i in range(len(boundaries_ibtype_new)):
                node_ids = boundaries_ibtype_new[i]['node_id']
                if ibtype is None: # open boundary
                    node_ids1 = node_ids
                elif ibtype.endswith('4'): # weir boundary
                    node_ids1 = [node_id[0] for node_id in node_ids]
                    node_ids2 = [node_id[1] for node_id in node_ids]
                else: # other boundaries
                    node_ids1 = node_ids
                    
                boundary_new = boundaries_ibtype_new[i]
                
                if ibtype is None: # open boundary
                    new_node_ids = [str(node_mapping[nid]) for nid in node_ids]
                elif ibtype.endswith('4'): # weir boundary
                    new_node_ids = [(str(node_mapping[node_ids1[j]]), str(node_mapping[node_ids2[j]])) for j in range(len(node_ids))]
                else: # other boundaries
                    new_node_ids = [str(node_mapping[nid]) for nid in node_ids]
                    
                boundary_new['node_id'] = new_node_ids

        # Merge land mesh boundaries and channel mesh boundaries
        boundaries_land = land_mesh.boundaries.to_dict().copy()
        for ibtype in land_mesh.boundaries.to_dict().keys():
            for i in range(len(boundaries_land[ibtype])):
                if ibtype is None: # open boundary
                    boundaries_land[ibtype][i]['node_id'] = [str(nid) for nid in boundaries_land[ibtype][i]['node_id']]
                elif ibtype.endswith('4'): # weir boundary
                    boundaries_land[ibtype][i]['node_id'] = [(str(nid[0]), str(nid[1])) for nid in boundaries_land[ibtype][i]['node_id']]
                else: # other boundaries
                    boundaries_land[ibtype][i]['node_id'] = [str(nid) for nid in boundaries_land[ibtype][i]['node_id']]

            if ibtype in boundaries_new.keys():
                boundaries_new[ibtype] = boundaries_land[ibtype] + boundaries_new[ibtype]
            else:
                boundaries_new[ibtype] = boundaries_land[ibtype]

        # Remove nodes that are members of matching nodes from the merged boundaries
        print("Removing matching nodes from boundaries...")
        # Create a set of all matching nodes (both from land and channel mesh)
        # For channel mesh nodes (keys), use the mapped IDs
        matching_node_strings = set(str(node_mapping[nid]) for nid in matching_nodes.keys()) | set(str(nid) for nid in matching_nodes.values())
        for ibtype in boundaries_new.keys():
            boundaries_ibtype = boundaries_new[ibtype]
            new_boundaries = []
            
            for boundary in boundaries_ibtype:
                node_ids = boundary['node_id']
                original_length = len(node_ids)
                
                # Extract first side of nodes for all boundary types
                if isinstance(node_ids[0], tuple):
                    node_ids1 = [node_id[0] for node_id in node_ids]
                    node_ids2 = [node_id[1] for node_id in node_ids]
                else:
                    node_ids1 = node_ids
                
                # Filter nodes once for all boundary types
                keep_indices, should_keep = self._filter_boundary_nodes(node_ids1, original_length, matching_node_strings)
                
                if should_keep:
                    if ibtype is None: # open boundary
                        boundary['node_id'] = [node_ids[i] for i in keep_indices]
                    elif ibtype.endswith('4'): # weir boundary
                        boundary['node_id'] = [(node_ids1[i], node_ids2[i]) for i in keep_indices]
                        # Filter other weir parameters using the same indices
                        boundary['barrier_height'] = [boundary['barrier_height'][i] for i in keep_indices]
                        boundary['subcritical_flow_coefficient'] = [boundary['subcritical_flow_coefficient'][i] for i in keep_indices]
                        boundary['supercritical_flow_coefficient'] = [boundary['supercritical_flow_coefficient'][i] for i in keep_indices]
                    else: # other boundaries
                        boundary['node_id'] = [node_ids[i] for i in keep_indices]
                    new_boundaries.append(boundary)
            
            boundaries_new[ibtype] = new_boundaries

        print("Creating final merged mesh...")
        merged_mesh = AdcircMesh(nodes=combined_nodes, elements=combined_elements, boundaries=boundaries_new)
        
        # Add VEW boundaries
        print("Finding paired nodes...")
        paired_nodes = self._find_paired_nodes(merged_mesh.nodes, config['tolerance'])
        
        if paired_nodes:
            print(f"Found {len(paired_nodes)} paired nodes")
            print("Creating VEW boundaries...")
            vewboundaries = self._create_vew_boundaries(
                paired_nodes, 
                land_mesh,
                node_mapping,
                combined_elements,
                config
            )
            
            print("Adding VEW boundaries to mesh...")
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
        else:
            print("No paired nodes found for VEW boundaries")
            
 
        if '64' in boundaries:
            # Merge endpoints of non-looped VEW boundaries
            print("Merging endpoints of non-looped VEW boundaries...")

            # Create a list of endpoints of VEW boundaries that are not looped
            endpoint_pairs = {}
            vew_boundaries = boundaries['64']
            for boundary in vew_boundaries:
                node_ids = boundary['node_id']
                if len(node_ids) >= 2:
                    # Check if the boundary is non-looped
                    if node_ids[0] != node_ids[-1]:
                        # Create a list of node pairs
                        endpoint_pairs[int(node_ids[0][0])] = int(node_ids[0][1])
                        endpoint_pairs[int(node_ids[-1][0])] = int(node_ids[-1][1])
            endpoints_land = list(endpoint_pairs.keys())
            endpoints_land_str = [str(i) for i in endpoints_land]
            
            # Update elevations at the endpoint land nodes to the mean of the land and channel nodes
            if endpoint_pairs:
                new_nodes = merged_mesh.nodes.copy()
                new_elements = merged_mesh.elements.elements.copy()
                boundaries = merged_mesh.boundaries.to_dict().copy()
                
                for land_node, channel_node in endpoint_pairs.items():
                    new_nodes.loc[land_node, 'value_1'] = (new_nodes.loc[land_node, 'value_1'] + new_nodes.loc[channel_node, 'value_1']) * 0.5
                    
                # Drop the endpoint channel nodes
                new_nodes = new_nodes[~new_nodes.index.isin(endpoint_pairs.values())]
                
                # Create a node mapping
                node_mapping = {str(i_old): str(i+1) for i, i_old in enumerate(new_nodes.index)}
                for land_node, channel_node in endpoint_pairs.items():
                    node_mapping[str(channel_node)] = node_mapping[str(land_node)]
                node_mapping_int = {int(k): int(v) for k, v in node_mapping.items()}
                    
                # Update the index of the combined_nodes
                new_nodes.index = range(1, len(new_nodes) + 1)
                
                # Update the node indices in new_elements
                for col in new_elements.columns[1:]:
                    new_elements[col] = new_elements[col].map(node_mapping_int)
                
                # Update the boundaries
                new_boundaries = {}
                for ibtype in boundaries.keys():
                    if ibtype == '64':
                        boundaries_ibtype = boundaries[ibtype]
                        new_boundaries_ibtype = []
                        for boundary in boundaries_ibtype:
                            keep_indices = [i for i, nid in enumerate(boundary['node_id']) if nid[0] not in endpoints_land_str]
                            if len(keep_indices) >= 2:
                                # Get the filtered node pairs using keep_indices
                                filtered_node_pairs = [boundary['node_id'][i] for i in keep_indices]
                                filtered_barrier_heights = [boundary['barrier_height'][i] for i in keep_indices]
                                filtered_subcrit = [boundary['subcritical_flow_coefficient'][i] for i in keep_indices]
                                filtered_supercrit = [boundary['supercritical_flow_coefficient'][i] for i in keep_indices]
                                
                                # Create new boundary with mapped node IDs
                                boundary['node_id'] = [(node_mapping[nid[0]], node_mapping[nid[1]]) for nid in filtered_node_pairs]
                                boundary['barrier_height'] = filtered_barrier_heights
                                boundary['subcritical_flow_coefficient'] = filtered_subcrit
                                boundary['supercritical_flow_coefficient'] = filtered_supercrit
                                new_boundaries_ibtype.append(boundary)
                        if new_boundaries_ibtype:
                            new_boundaries[ibtype] = new_boundaries_ibtype
                    else:
                        boundaries_ibtype = boundaries[ibtype]
                        new_boundaries_ibtype = []
                        for boundary in boundaries_ibtype:
                            if isinstance(boundary['node_id'][0], tuple):
                                boundary['node_id'] = [(node_mapping[nid[0]], node_mapping[nid[1]]) for nid in boundary['node_id']]
                            else:
                                boundary['node_id'] = [node_mapping[nid] for nid in boundary['node_id']]
                            new_boundaries_ibtype.append(boundary)
                        new_boundaries[ibtype] = new_boundaries_ibtype
                
                merged_mesh = AdcircMesh(
                    nodes=new_nodes,
                    elements=new_elements,
                    boundaries=new_boundaries
                )
        
        # Reorder the boundaries to have the VEW boundaries last
        boundaries = merged_mesh.boundaries.to_dict()
        if '64' in boundaries:
            vew_boundaries = boundaries.pop('64')
            boundaries['64'] = vew_boundaries           
        
        merged_mesh = AdcircMesh(
            nodes=merged_mesh.nodes,
            elements=merged_mesh.elements.elements,
            boundaries=boundaries
        )
        
        print("Mesh merge complete!")
        return merged_mesh


class MergedNodesStrategy(MergeStrategy):
    """Strategy that merges duplicate nodes at boundaries."""
    
    def merge(self, channel_mesh: AdcircMesh, land_mesh: AdcircMesh, config: Dict) -> AdcircMesh:
        """Merge meshes by merging duplicate nodes."""
        print("\nStarting mesh merge with node merging...")
        
        # Find matching nodes
        manipulator = VEWBoundaryManipulator()
        matching_nodes = manipulator.find_matching_nodes(
            channel_mesh.nodes, 
            land_mesh.nodes, 
            config['tolerance']
        )
        
        # Create node mapping for non-matching nodes
        print("Creating node mapping...")
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
        print("Combining nodes...")
        combined_nodes = land_mesh.nodes.copy()
        
        # Update values for matching nodes based on configuration
        if config['use_channel_values']:
            source_mesh = channel_mesh
            source_ids = matching_nodes.keys()
            target_ids = matching_nodes.values()
            value_source = "channel"
        else:
            source_mesh = land_mesh
            source_ids = matching_nodes.values()
            target_ids = matching_nodes.values()
            value_source = "land"

        print(f"Updating values at matching nodes using {value_source} mesh values...")
        for source_id, target_id in zip(source_ids, target_ids):
            combined_nodes.loc[target_id, 'value_1'] = source_mesh.nodes.loc[source_id, 'value_1']
        
        # Add non-matching nodes from channel mesh
        print("Adding non-matching nodes...")
        new_nodes = []
        for old_id in channel_mesh.nodes.index:
            if old_id not in matching_nodes:
                node = channel_mesh.nodes.loc[old_id].copy()
                node.name = node_mapping[old_id]
                new_nodes.append(node)
        
        if new_nodes:
            new_nodes_df = pd.DataFrame(new_nodes)
            combined_nodes = pd.concat([combined_nodes, new_nodes_df])
        print(f"Combined nodes: {len(combined_nodes)}")
        
        # Update element connectivity
        print("Updating element connectivity...")
        channel_elements = channel_mesh.elements.elements.copy()
        for col in channel_elements.columns[1:]:
            channel_elements[col] = channel_elements[col].map(node_mapping)
        
        combined_elements = pd.concat([land_mesh.elements.elements, channel_elements])
        # Reset element IDs to be 1-based
        combined_elements.index = range(1, len(combined_elements) + 1)
        print(f"Combined elements: {len(combined_elements)}")
        
        # Update boundaries of channel mesh
        print("Updating boundary node numbers...")
        boundaries_new = channel_mesh.boundaries.to_dict().copy()
        for ibtype in boundaries_new.keys():
            boundaries_ibtype_new = boundaries_new[ibtype]
            
            for i in range(len(boundaries_ibtype_new)):
                node_ids = boundaries_ibtype_new[i]['node_id']
                if ibtype is None: # open boundary
                    node_ids1 = node_ids
                elif ibtype.endswith('4'): # weir boundary
                    node_ids1 = [node_id[0] for node_id in node_ids]
                    node_ids2 = [node_id[1] for node_id in node_ids]
                else: # other boundaries
                    node_ids1 = node_ids
                    
                boundary_new = boundaries_ibtype_new[i]
                
                if ibtype is None: # open boundary
                    new_node_ids = [str(node_mapping[nid]) for nid in node_ids]
                elif ibtype.endswith('4'): # weir boundary
                    new_node_ids = [(str(node_mapping[node_ids1[j]]), str(node_mapping[node_ids2[j]])) for j in range(len(node_ids))]
                else: # other boundaries
                    new_node_ids = [str(node_mapping[nid]) for nid in node_ids]
                    
                boundary_new['node_id'] = new_node_ids

        # Merge land mesh boundaries and channel mesh boundaries
        boundaries_land = land_mesh.boundaries.to_dict().copy()
        for ibtype in land_mesh.boundaries.to_dict().keys():
            for i in range(len(boundaries_land[ibtype])):
                if ibtype is None: # open boundary
                    boundaries_land[ibtype][i]['node_id'] = [str(nid) for nid in boundaries_land[ibtype][i]['node_id']]
                elif ibtype.endswith('4'): # weir boundary
                    boundaries_land[ibtype][i]['node_id'] = [(str(nid[0]), str(nid[1])) for nid in boundaries_land[ibtype][i]['node_id']]
                else: # other boundaries
                    boundaries_land[ibtype][i]['node_id'] = [str(nid) for nid in boundaries_land[ibtype][i]['node_id']]

            if ibtype in boundaries_new.keys():
                boundaries_new[ibtype] = boundaries_land[ibtype] + boundaries_new[ibtype]
            else:
                boundaries_new[ibtype] = boundaries_land[ibtype]

        # Find boundary edges in the final mesh
        print("Finding boundary edges in final mesh...")
        edge_counts = defaultdict(int)
        for _, element in combined_elements.iterrows():
            nodes = element[1:4].astype(int).tolist()
            edges = [
                frozenset([nodes[0], nodes[1]]),
                frozenset([nodes[1], nodes[2]]),
                frozenset([nodes[2], nodes[0]])
            ]
            for edge in edges:
                edge_counts[edge] += 1

        # Create set of boundary nodes (nodes that belong to edges with count 1)
        boundary_nodes = set()
        for edge, count in edge_counts.items():
            if count == 1:
                boundary_nodes.update(edge)

        # Remove non-boundary nodes from the merged boundaries
        print("Removing non-boundary nodes...")
        for ibtype in boundaries_new.keys():
            boundaries_ibtype = boundaries_new[ibtype]
            new_boundaries = []
            
            for boundary in boundaries_ibtype:
                node_ids = boundary['node_id']
                if ibtype is None: # open boundary
                    # Only keep nodes that are on the boundary
                    node_ids = [nid for nid in node_ids if int(nid) in boundary_nodes]
                    # Only keep node strings with at least 2 nodes
                    if len(node_ids) >= 2:
                        # Check if the string has been reduced to just endpoints
                        if len(node_ids) == 2 and len(boundary['node_id']) > 3:
                            # Skip this boundary as it's been reduced to just endpoints
                            continue
                        boundary['node_id'] = node_ids
                        new_boundaries.append(boundary)
                elif ibtype.endswith('4'): # weir boundary
                    node_ids1 = [int(node_id[0]) for node_id in node_ids]
                    node_ids2 = [int(node_id[1]) for node_id in node_ids]
                    # Only keep nodes that are on the boundary
                    keep_indices = [i for i, (n1, n2) in enumerate(zip(node_ids1, node_ids2)) 
                                  if n1 in boundary_nodes and n2 in boundary_nodes]
                    # Only keep weir segments where both sides have nodes
                    if len(keep_indices) >= 2:
                        # Check if either side has been reduced to just endpoints
                        if len(keep_indices) == 2 and len(boundary['node_id']) > 3:
                            # Skip this boundary as one side has been reduced to just endpoints
                            continue
                        boundary['node_id'] = [(str(node_ids1[i]), str(node_ids2[i])) for i in keep_indices]
                        # Filter other weir parameters using the same indices
                        boundary['barrier_height'] = [boundary['barrier_height'][i] for i in keep_indices]
                        boundary['subcritical_flow_coefficient'] = [boundary['subcritical_flow_coefficient'][i] for i in keep_indices]
                        boundary['supercritical_flow_coefficient'] = [boundary['supercritical_flow_coefficient'][i] for i in keep_indices]
                        new_boundaries.append(boundary)
                else: # other boundaries
                    # Only keep nodes that are on the boundary
                    node_ids = [nid for nid in node_ids if int(nid) in boundary_nodes]
                    # Only keep node strings with at least 2 nodes
                    if len(node_ids) >= 2:
                        # Check if the string has been reduced to just endpoints
                        if len(node_ids) == 2 and len(boundary['node_id']) > 3:
                            # Skip this boundary as it's been reduced to just endpoints
                            continue
                        boundary['node_id'] = node_ids
                        new_boundaries.append(boundary)
            
            boundaries_new[ibtype] = new_boundaries
        
        # Reorder the boundaries to have the VEW boundaries last
        if '64' in boundaries_new:
            vew_boundaries = boundaries_new.pop('64')
            boundaries_new['64'] = vew_boundaries           
        
        print("Creating final merged mesh...")
        merged_mesh = AdcircMesh(nodes=combined_nodes, elements=combined_elements, boundaries=boundaries_new)
        print(f"Merged {len(matching_nodes)} duplicate nodes")
        print(f"Used {value_source} mesh values at {len(matching_nodes)} matching nodes")
        print("Mesh merge complete!")
        
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
            'use_channel_values': True  # Default to using channel mesh values
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

    def with_channel_values(self, enabled: bool = True) -> 'MeshMerger':
        """Choose whether to use channel mesh values at matching nodes."""
        self._config['use_channel_values'] = enabled
        return self

    def merge(self, boundary_mode: str = "merge") -> AdcircMesh:
        """Merge the meshes using the selected strategy."""
        if boundary_mode == "vew":
            strategy = VEWBoundaryStrategy()
        else:  # "merge"
            strategy = MergedNodesStrategy()
            
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
        '-o', '--output',
        default='merged_mesh.grd',
        help='Path to save the merged mesh (default: merged_mesh.grd)'
    )
    parser.add_argument(
        '-d', '--description',
        default='merged',
        help='Description for the merged mesh (default: merged)'
    )
    parser.add_argument(
        "-b", "--boundary-mode",
        choices=["merge", "vew"],
        default="merge",
        help="How to handle mesh boundaries: 'merge' to merge nodes, 'vew' to use VEW boundaries (default: merge)"
    )
    parser.add_argument(
        "--use-land-values",
        action="store_true",
        help="Use land mesh values at matching nodes (default: use channel mesh values)"
    )
    
    args = parser.parse_args()
    
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
    
    if args.use_land_values:
        merger.with_channel_values(False)
    
    merged_mesh = merger.merge(args.boundary_mode)
    
    print(f"\nMerged mesh info:")
    print(f"Number of nodes: {len(merged_mesh.nodes)}")
    print(f"Number of elements: {len(merged_mesh.elements.elements)}")
    
    # Save the merged mesh
    merged_mesh.description = args.description
    merged_mesh.write(args.output, overwrite=True)
    print(f"\nMerged mesh saved to: {args.output}")
    
    return 0


if __name__ == "__main__":
    main()