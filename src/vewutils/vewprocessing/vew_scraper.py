#!/usr/bin/env python3
"""
Module for scraping VEW boundaries from ADCIRC meshes and saving them to YAML format.
"""

import argparse
import os
import yaml
import numpy as np
import pandas as pd
from adcircpy import AdcircMesh
from typing import Dict, List, Tuple


class VEWScraper:
    """Class for scraping VEW boundaries from ADCIRC meshes."""

    def __init__(self, mesh: AdcircMesh):
        """Initialize the VEW scraper with a mesh."""
        self._mesh = mesh
        self._map_elem_node_next = {0: 1, 1: 2, 2: 0}

    def _get_node_elements_for_nodes(self, mesh: AdcircMesh, node_ids: List[int]) -> Dict[int, set]:
        """Get elements connected to specific nodes only (optimized version)."""
        node_elements = {}
        triangles = mesh.elements.elements
        
        for node_id in node_ids:
            node_elements[node_id] = set()
        
        # Only iterate through triangles once and check for our target nodes
        for idx, row in triangles.iterrows():
            triangle_nodes = [row['node_1'], row['node_2'], row['node_3']]
            for node_id in node_ids:
                if node_id in triangle_nodes:
                    node_elements[node_id].add(idx)
        
        return node_elements

    def _get_node_neighbors_for_nodes(self, mesh: AdcircMesh, node_ids: List[int]) -> Dict[int, set]:
        """Get neighbors for specific nodes only (optimized version)."""
        node_neighbors = {}
        for node_id in node_ids:
            node_neighbors[node_id] = set()
        
        triangles = mesh.elements.elements
        
        # Only iterate through triangles once and build neighbors for our target nodes
        for idx, row in triangles.iterrows():
            triangle_nodes = [row['node_1'], row['node_2'], row['node_3']]
            
            # Check if any of our target nodes are in this triangle
            target_nodes_in_triangle = [n for n in node_ids if n in triangle_nodes]
            
            if target_nodes_in_triangle:
                # For each target node in this triangle, add the other nodes as neighbors
                for target_node in target_nodes_in_triangle:
                    for triangle_node in triangle_nodes:
                        if triangle_node != target_node:
                            node_neighbors[target_node].add(triangle_node)
        
        return node_neighbors

    def _remap_elements_vectorized(self, elements_df, map_node: Dict[int, int]):
        """Vectorized remapping of element node IDs using numpy operations."""
        # Convert to numpy array for faster operations
        elements_array = elements_df[['node_1', 'node_2', 'node_3']].values.copy()
        
        # Create arrays for old and new node IDs
        old_nodes = np.array(list(map_node.keys()))
        new_nodes = np.array(list(map_node.values()))
        
        # Create a mapping array that covers the full range of node IDs
        max_node_id = max(elements_array.max(), old_nodes.max(), new_nodes.max()) + 1
        node_map_array = np.arange(max_node_id)  # Default: no change
        
        # Apply the mapping
        node_map_array[old_nodes] = new_nodes
        
        # Only remap nodes that are within the mapping array bounds
        mask = elements_array < len(node_map_array)
        elements_array[mask] = node_map_array[elements_array[mask]]
        
        # Update the dataframe
        elements_df[['node_1', 'node_2', 'node_3']] = elements_array
        
        return elements_df

    def _remap_boundary_nodes_vectorized(self, boundaries: Dict, map_node: Dict[int, int]):
        """Vectorized remapping of boundary node IDs using numpy operations."""
        # Create numpy arrays for mapping
        old_nodes = np.array(list(map_node.keys()))
        new_nodes = np.array(list(map_node.values()))
        
        # Create a mapping array that covers the full range of node IDs
        max_node_id = max(max(map_node.keys()), max(map_node.values())) + 1
        node_map_array = np.arange(max_node_id)  # Default: no change
        node_map_array[old_nodes] = new_nodes
        
        # Process each boundary type
        for ibtype in boundaries.keys():
            for i in range(len(boundaries[ibtype])):
                node_ids = boundaries[ibtype][i]['node_id']
                
                if ibtype is None:
                    # Single node per entry
                    node_array = np.array(node_ids)
                    # Only remap nodes that are in our mapping
                    mask = node_array < len(node_map_array)
                    node_array[mask] = node_map_array[node_array[mask]]
                    boundaries[ibtype][i]['node_id'] = node_array.tolist()
                    
                elif ibtype.endswith('4'):
                    # Pairs of nodes (like VEW boundaries)
                    if node_ids:  # Check if not empty
                        node_pairs = np.array(node_ids)
                        # Remap both columns
                        for col in range(node_pairs.shape[1]):
                            mask = node_pairs[:, col] < len(node_map_array)
                            node_pairs[mask, col] = node_map_array[node_pairs[mask, col]]
                        boundaries[ibtype][i]['node_id'] = [tuple(row) for row in node_pairs]
                    
                else:
                    # Other boundary types with variable node lists
                    if isinstance(node_ids[0], (list, tuple)) and len(node_ids[0]) > 1:
                        # Multi-dimensional array
                        node_array = np.array(node_ids)
                        if node_array.ndim == 2:
                            for col in range(node_array.shape[1]):
                                mask = node_array[:, col] < len(node_map_array)
                                node_array[mask, col] = node_map_array[node_array[mask, col]]
                        boundaries[ibtype][i]['node_id'] = node_array.tolist()
                    else:
                        # Single node per entry
                        node_array = np.array(node_ids)
                        mask = node_array < len(node_map_array)
                        node_array[mask] = node_map_array[node_array[mask]]
                        boundaries[ibtype][i]['node_id'] = node_array.tolist()
        
        return boundaries

    def _strip_vewstring(self, mesh: AdcircMesh, ivewboundary: int) -> Tuple[AdcircMesh, List[Dict]]:
        """Strip a single VEW boundary from the mesh and return the modified mesh and VEW string."""
        print(f"    - Initializing boundary data...")
        boundaries = mesh.boundaries.to_dict().copy()
        vewboundaries = boundaries['64']
        
        print(f"    - Copying VEW boundary data...")
        vewboundary = vewboundaries[ivewboundary]

        print(f"    - Extracting bank and channel nodes...")
        bank_nodes = [int(vewboundary['node_id'][i][0]) for i in range(len(vewboundary['node_id']))]
        channel_nodes = [int(vewboundary['node_id'][i][1]) for i in range(len(vewboundary['node_id']))]
        print(f"    - Boundary has {len(bank_nodes)} node pairs")
        print(f"    - Debug: Initial bank_nodes: {bank_nodes}")
        print(f"    - Debug: Initial channel_nodes: {channel_nodes}")
        
        # Check for duplicates or issues
        bank_set = set(bank_nodes)
        channel_set = set(channel_nodes)
        overlap = bank_set.intersection(channel_set)
        if overlap:
            print(f"    - WARNING: Found {len(overlap)} nodes that appear in both bank and channel: {list(overlap)}")
        
        # Check if all nodes exist in the mesh
        for i, (bank_node, channel_node) in enumerate(zip(bank_nodes, channel_nodes)):
            if bank_node not in mesh.nodes.index:
                print(f"    - ERROR: bank_node {bank_node} (pair {i}) not found in mesh!")
            if channel_node not in mesh.nodes.index:
                print(f"    - ERROR: channel_node {channel_node} (pair {i}) not found in mesh!")

        print(f"    - Getting node elements for specific nodes only...")
        # Only get node elements for the two nodes we actually need
        target_nodes_for_elements = [bank_nodes[0], bank_nodes[1]]
        node_elements = self._get_node_elements_for_nodes(mesh, target_nodes_for_elements)
        
        print(f"    - Getting node neighbors for boundary endpoints...")
        # Only get node neighbors for the nodes we actually need
        target_nodes_for_neighbors = []
        if bank_nodes[0] != bank_nodes[-1]:  # Not a loop
            target_nodes_for_neighbors = [bank_nodes[0], bank_nodes[-1], channel_nodes[0], channel_nodes[-1]]
        node_neighbors = self._get_node_neighbors_for_nodes(mesh, target_nodes_for_neighbors)

        print(f"    - Copying mesh data structures...")
        # Making a new mesh without the vew boundary
        nodes = mesh.nodes.copy()
        nodes_new = mesh.nodes.copy()
        elements_new = mesh.elements.elements.copy()
        boundaries_new = boundaries.copy()

        print(f"    - Removing bank nodes from node table...")
        # STEP 1: Drop the bank nodes first
        nodes_new = nodes_new.drop(index=bank_nodes)
        
        print(f"    - Creating node mapping table for remaining nodes...")
        # STEP 2: Create mapping table AFTER dropping bank nodes
        # Map old node IDs to new sequential IDs (1, 2, 3, ...)
        map_node = {}
        idx = 1
        for old_node_id in nodes_new.index:
            map_node[old_node_id] = idx
            idx += 1
        
        print(f"    - Debug: bank_nodes = {bank_nodes}")
        print(f"    - Debug: channel_nodes = {channel_nodes}")
        print(f"    - Debug: len(bank_nodes) = {len(bank_nodes)}, len(channel_nodes) = {len(channel_nodes)}")
        print(f"    - Debug: nodes_new.index contains {len(nodes_new.index)} nodes")
        print(f"    - Debug: map_node contains {len(map_node)} entries")
        
        # Add bank node mappings to their corresponding channel nodes
        for i, bank_node in enumerate(bank_nodes):
            channel_node = channel_nodes[i]
            print(f"    - Debug: Mapping bank_node {bank_node} (index {i}) to channel_node {channel_node}")
            
            # Check if channel_node exists in map_node
            if channel_node not in map_node:
                print(f"    - ERROR: channel_node {channel_node} not found in map_node!")
                print(f"    - Debug: First 10 keys in map_node: {list(map_node.keys())[:10]}")
                print(f"    - Debug: Last 10 keys in map_node: {list(map_node.keys())[-10:]}")
                print(f"    - Debug: Channel node {channel_node} in original nodes: {channel_node in mesh.nodes.index}")
                print(f"    - Debug: Channel node {channel_node} in nodes_new: {channel_node in nodes_new.index}")
                print(f"    - Debug: Bank node {bank_node} in original nodes: {bank_node in mesh.nodes.index}")
                print(f"    - Debug: Bank node {bank_node} in nodes_new: {bank_node in nodes_new.index}")
                
                # Check if this channel node was accidentally removed as a bank node
                if channel_node in bank_nodes:
                    print(f"    - WARNING: channel_node {channel_node} is also in bank_nodes at indices: {[j for j, bn in enumerate(bank_nodes) if bn == channel_node]}")
                
                raise KeyError(f"Channel node {channel_node} not found in map_node dictionary. This suggests the channel node was removed during bank node removal.")
            
            map_node[bank_node] = map_node[channel_node]
        
        print(f"    - Reindexing remaining nodes...")
        # STEP 3: Reindex the nodes using the new mapping
        nodes_new.index = nodes_new.index.map(lambda x: map_node[x])

        print(f"    - Updating element node IDs (vectorized)...")
        # STEP 4: Now remap elements with the correct mapping
        elements_new = self._remap_elements_vectorized(elements_new, map_node)

        print(f"    - Removing VEW boundary from boundary table...")
        # Drop the targeted vew boundary
        boundaries_new['64'].pop(ivewboundary)

        print(f"    - Updating boundary node IDs (vectorized)...")
        # STEP 5: Remap boundaries with the correct mapping
        boundaries_new = self._remap_boundary_nodes_vectorized(boundaries_new, map_node)

        print(f"    - Creating new mesh object...")            
        new_mesh = AdcircMesh(nodes=nodes_new, elements=elements_new, boundaries=boundaries_new)

        print(f"    - Processing triangle connectivity...")
        triangles = mesh.elements.elements.iloc[:, 1:4]

        # Making a vewstring
        # Reverse node strings if it's not seeing the bank side on the right
        print(f"    - Determining node orientation...")
        first_node_elems = node_elements[bank_nodes[0]]
        second_node_elems = node_elements[bank_nodes[1]]
        common_elems = list(set(first_node_elems) & set(second_node_elems))
        if len(common_elems) != 1:
            raise ValueError("The number of common elements between the two bank nodes is not 1")
        elem_nodes = list(triangles.loc[common_elems[0]])
        index = elem_nodes.index(bank_nodes[1])
        if elem_nodes[self._map_elem_node_next[index]] != bank_nodes[0]:
            bank_nodes.reverse()
            channel_nodes.reverse()

        print(f"    - Processing boundary endpoints...")
        # If this node string is not a loop, find the common node between the bank and channel nodes at the ends
        if bank_nodes[0] != bank_nodes[-1]:
            nneigh_bank = node_neighbors[bank_nodes[0]]
            nneigh_channel = node_neighbors[channel_nodes[0]]
            common_node = list(set(nneigh_bank) & set(nneigh_channel))
            if len(common_node) != 1:
                raise ValueError("The number of common nodes between the bank and channel nodes is not 1")
            bank_nodes.insert(0, common_node[0])
            channel_nodes.insert(0, common_node[0])
            
            nneigh_bank = node_neighbors[bank_nodes[-1]]
            nneigh_channel = node_neighbors[channel_nodes[-1]]
            common_node = list(set(nneigh_bank) & set(nneigh_channel))
            if len(common_node) != 1:
                raise ValueError("The number of common nodes between the bank and channel nodes is not 1")
            bank_nodes.append(common_node[0])
            channel_nodes.append(common_node[0])
        
        print(f"    - Building VEW string data...")    
        # Create a vewstring
        vewstring = []
        for i in range(len(bank_nodes)):
            bank_node = bank_nodes[i]
            channel_node = channel_nodes[i]
            
            # Check that both nodes exist in their respective data structures
            if channel_node not in map_node:
                print(f"    - ERROR: channel_node {channel_node} not found in map_node during vewstring creation!")
                raise KeyError(f"Channel node {channel_node} not found in map_node during vewstring creation")
            if bank_node not in nodes.index:
                print(f"    - ERROR: bank_node {bank_node} not found in original nodes during vewstring creation!")
                raise KeyError(f"Bank node {bank_node} not found in original nodes during vewstring creation")
            if channel_node not in nodes.index:
                print(f"    - ERROR: channel_node {channel_node} not found in original nodes during vewstring creation!")
                raise KeyError(f"Channel node {channel_node} not found in original nodes during vewstring creation")
            
            vewstring.append({
                'node_id': int(map_node[channel_node]),
                'x': float(nodes.loc[bank_node, 'x']),
                'y': float(nodes.loc[channel_node, 'y']),
                'bank_elevation': float(nodes.loc[bank_node, 'value_1']),
                'bank_mannings_n': 0.03
            })
        
        print(f"    - Boundary processing complete")       
        return new_mesh, vewstring

    def strip_vewstrings(self, boundary_indices=None) -> Tuple[AdcircMesh, Dict]:
        """Strip specified VEW boundaries from the mesh and return the modified mesh and VEW strings.
        
        Args:
            boundary_indices: List of boundary indices to scrape. If None, scrapes all boundaries.
        """
        boundaries = self._mesh.boundaries.to_dict().copy()
        vewboundaries = boundaries['64']
        
        # If no specific boundaries specified, scrape all
        if boundary_indices is None:
            boundary_indices = list(range(len(vewboundaries)))
            print(f"Processing all {len(vewboundaries)} VEW boundaries...")
        else:
            print(f"Processing {len(boundary_indices)} specified VEW boundaries: {boundary_indices}")
        
        # Validate boundary indices
        for idx in boundary_indices:
            if idx < 0 or idx >= len(vewboundaries):
                raise ValueError(f"Boundary index {idx} is out of range. Available indices: 0-{len(vewboundaries)-1}")

        vewstrings = []
        mesh = self._mesh
        
        # Sort indices in descending order to avoid index shifting issues when removing boundaries
        sorted_indices = sorted(boundary_indices, reverse=True)
        for i, original_idx in enumerate(sorted_indices):
            print(f"  Processing boundary {original_idx + 1} ({i+1}/{len(boundary_indices)})...")
            # When processing in descending order, the current index is simply the original index
            # because we remove boundaries from higher indices first, leaving lower indices unchanged
            current_idx = original_idx
            mesh, vewstring = self._strip_vewstring(mesh, current_idx)
            vewstrings.append((original_idx, vewstring))
            print(f"    ✓ Extracted {len(vewstring)} nodes from boundary {original_idx}")
        
        # Sort back to original order and extract vewstrings
        vewstrings.sort(key=lambda x: x[0])
        vewstrings = [vew[1] for vew in vewstrings]
        vewstrings = {'vewstrings': vewstrings}

        print(f"✓ Successfully processed {len(boundary_indices)} VEW boundaries")
        return mesh, vewstrings


def get_parser():
    parser = argparse.ArgumentParser(add_help=False, description="Scrape VEW boundaries from an ADCIRC mesh and save them to YAML format")
    parser.add_argument(
        "input_mesh",
        help="Path to the input mesh file with VEW boundaries"
    )
    parser.add_argument(
        '-o', '--output-mesh',
        required=True,
        help='Path to save the mesh without VEW boundaries'
    )
    parser.add_argument(
        '-y', '--output-yaml',
        required=True,
        help='Path to save the extracted VEW strings in YAML format'
    )
    parser.add_argument(
        '-d', '--description',
        default='Generated by vew_scraper',
        help='Description for the output mesh (default: Generated by vew_scraper)'
    )
    parser.add_argument(
        '-b', '--boundaries',
        type=int,
        nargs='+',
        help='Indices of VEW boundaries to scrape (0-based). If not specified, all boundaries will be scraped. Example: -b 0 2 3'
    )
    return parser

def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    from adcircpy import AdcircMesh
    import yaml
    
    print("=== VEW Boundary Scraper ===")
    print(f"Input mesh: {args.input_mesh}")
    print(f"Output mesh: {args.output_mesh}")
    print(f"Output YAML: {args.output_yaml}")
    print()
    
    # Read the mesh file
    print("Loading mesh file...")
    mesh = AdcircMesh.open(args.input_mesh)
    print(f"✓ Mesh loaded successfully")
    print(f"  Nodes: {len(mesh.nodes)}")
    print(f"  Elements: {len(mesh.elements.elements)}")
    
    # Check for VEW boundaries
    boundaries = mesh.boundaries.to_dict()
    if '64' not in boundaries or not boundaries['64']:
        print("⚠ Warning: No VEW boundaries (type 64) found in the mesh")
        return 1
    
    vew_count = len(boundaries['64'])
    print(f"  VEW boundaries found: {vew_count}")
    print()
    
    # Strip VEW boundaries
    print("Extracting VEW boundaries...")
    scraper = VEWScraper(mesh)
    mesh_new, vewstrings = scraper.strip_vewstrings(boundary_indices=args.boundaries)
    mesh_new.description = args.description
    print()
    
    # Save the new mesh file
    print("Saving modified mesh...")
    mesh_new.write(args.output_mesh, overwrite=True)
    print(f"✓ Mesh saved to: {args.output_mesh}")
    print(f"  Final node count: {len(mesh_new.nodes)}")
    print(f"  Final element count: {len(mesh_new.elements.elements)}")
    print()
    
    # Write the YAML output to a file
    print("Saving VEW strings to YAML...")
    with open(args.output_yaml, 'w') as f:
        yaml.dump(vewstrings, f, sort_keys=False)
    print(f"✓ VEW strings saved to: {args.output_yaml}")
    
    total_vew_nodes = sum(len(vew) for vew in vewstrings['vewstrings'])
    print(f"  Total VEW nodes extracted: {total_vew_nodes}")
    print()
    print("=== Processing Complete ===")
    
    return 0

if __name__ == "__main__":
    main() 