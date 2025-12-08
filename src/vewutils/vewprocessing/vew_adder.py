#!/usr/bin/env python3
"""
Add VEW boundaries to an ADCIRC mesh based on VEW string definitions in a YAML file.

This program reads an ADCIRC mesh in fort.14 format and VEW string definitions from a YAML file,
then adds VEW boundaries to the mesh. The modified mesh is written out.

The program supports two modes:
1. Node ID mode: VEW strings specify node_id values that directly reference mesh nodes
2. Coordinate mode: VEW strings specify x,y coordinates that are matched to mesh nodes using spatial search
"""

import argparse
import yaml
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from adcircpy import AdcircMesh
import time
from scipy.spatial import cKDTree


class VEWBoundaryAdder:
    """Class for adding VEW boundaries to an ADCIRC mesh."""

    def __init__(self, mesh: AdcircMesh, coordinate_mode: bool = False, tolerance: float = 1e-6):
        """Initialize with an ADCIRC mesh.
        
        Args:
            mesh: ADCIRC mesh object
            coordinate_mode: If True, use coordinate-based node matching instead of node IDs
            tolerance: Distance tolerance for coordinate matching (in mesh units)
        """
        self._mesh = mesh
        self._map_elem_node_prev = {0: 2, 1: 0, 2: 1}
        self._map_elem_node_next = {0: 1, 1: 2, 2: 0}
        self._coordinate_mode = coordinate_mode
        self._tolerance = tolerance
        self._coordinate_tree = None
        self._node_id_array = None
        
        if coordinate_mode:
            self._build_coordinate_tree()

    def _build_coordinate_tree(self):
        """Build KDTree for fast coordinate-based node lookup."""
        print("Building coordinate search tree for fast node lookup...")
        start_time = time.time()
        
        # Extract coordinates and node IDs
        coordinates = np.column_stack([self._mesh.nodes['x'].values, self._mesh.nodes['y'].values])
        self._node_id_array = self._mesh.nodes.index.to_numpy()
        
        # Build KDTree for fast spatial search
        self._coordinate_tree = cKDTree(coordinates)
        
        build_time = time.time() - start_time
        print(f"✓ Coordinate search tree built in {build_time:.3f} seconds for {len(coordinates)} nodes")

    def _find_node_by_coordinates(self, x: float, y: float) -> int:
        """Find the closest mesh node to given coordinates within tolerance.
        
        Args:
            x: X coordinate
            y: Y coordinate
            
        Returns:
            Node ID of closest mesh node
            
        Raises:
            ValueError: If no node found within tolerance
        """
        if self._coordinate_tree is None:
            raise RuntimeError("Coordinate tree not built. Initialize with coordinate_mode=True.")
        
        # Query the tree for nearest neighbor
        query_point = np.array([x, y])
        distance, index = self._coordinate_tree.query(query_point)
        
        # Check if within tolerance
        if distance > self._tolerance:
            raise ValueError(f"No mesh node found within tolerance {self._tolerance} of coordinates ({x}, {y}). "
                           f"Closest node is at distance {distance:.6f}. "
                           f"Consider increasing tolerance or checking coordinate accuracy.")
        
        # Return the node ID
        node_id = self._node_id_array[index]
        return node_id

    def _process_vew_node_data(self, nodedata: Dict) -> Tuple[int, float]:
        """Process a single VEW node data entry to extract node ID and bank elevation.
        
        Args:
            nodedata: Dictionary containing either node_id or coordinates (x,y) plus bank_elevation
            
        Returns:
            Tuple of (node_id, bank_elevation)
        """
        if self._coordinate_mode:
            # Coordinate mode: use x,y coordinates to find node
            if 'x' not in nodedata or 'y' not in nodedata:
                raise ValueError("In coordinate mode, each VEW node must have 'x' and 'y' coordinates. "
                               f"Missing coordinates in node data: {nodedata}")
            
            x = float(nodedata['x'])
            y = float(nodedata['y'])
            node_id = self._find_node_by_coordinates(x, y)
            
        else:
            # Node ID mode: use node_id directly
            if 'node_id' not in nodedata:
                raise ValueError("In node ID mode, each VEW node must have 'node_id'. "
                               f"Missing node_id in node data: {nodedata}")
            
            node_id = nodedata['node_id']
            
            # Validate that the node exists in the mesh
            if node_id not in self._mesh.nodes.index:
                raise ValueError(f"Node ID {node_id} from VEW string not found in mesh. "
                               f"Please verify that the VEW string file contains valid node IDs that exist in the mesh.")
        
        # Get bank elevation
        bank_elevation = nodedata.get('bank_elevation', -99999.0)
        
        return node_id, bank_elevation

    def _get_node_identifier(self, nodedata: Dict) -> str:
        """Get a string identifier for a node (used for validation and error messages).
        
        Args:
            nodedata: Dictionary containing node information
            
        Returns:
            String identifier for the node
        """
        if self._coordinate_mode:
            x = nodedata.get('x', 'N/A')
            y = nodedata.get('y', 'N/A')
            return f"({x}, {y})"
        else:
            return str(nodedata.get('node_id', 'N/A'))

    def add_vew_string(self, vewstring: List[Dict]) -> AdcircMesh:
        """Add a single VEW string to the mesh."""
        start_time = time.time()
        print(f"Processing VEW string with {len(vewstring)} nodes...")

        # Create bank nodes and add them to the mesh
        node_start_time = time.time()
        nn = self._mesh.nodes.shape[0]
        map_node = {}
        id_new = []
        x_new = []
        y_new = []
        z_new = []

        # Validate minimum string length
        min_length = 4 if self._coordinate_mode else 4
        if len(vewstring) < min_length:
            raise ValueError(
                f"The length of a vewstring is {len(vewstring)}. It should be at least {min_length} nodes "
                f"for proper VEW boundary definition.")

        # For coordinate mode, we need to process nodes to get their IDs first
        if self._coordinate_mode:
            # Convert coordinates to node IDs
            processed_vewstring = []
            for i, nodedata in enumerate(vewstring):
                try:
                    node_id, bank_elevation = self._process_vew_node_data(nodedata)
                    processed_nodedata = {
                        'node_id': node_id,
                        'bank_elevation': bank_elevation
                    }
                    processed_vewstring.append(processed_nodedata)
                except ValueError as e:
                    node_identifier = self._get_node_identifier(nodedata)
                    raise ValueError(f"Error processing VEW node {i+1} {node_identifier}: {e}")
            
            vewstring = processed_vewstring

        # Check for closed string (first and last nodes are the same)
        first_node_id = vewstring[0]['node_id']
        last_node_id = vewstring[-1]['node_id']
        
        if len(vewstring) < 4 or (len(vewstring) < 3 and first_node_id == last_node_id):
            raise ValueError(
                "The length of a vewstring is {:d}. It should be greater than 2 for an open node string "
                "and should be greater than 2 for a closed node string.".format(len(vewstring)))

        # Add the second last node to the beginning and the second node to the end
        if first_node_id == last_node_id:
            vewstring.insert(0, vewstring[-2])
            vewstring.append(vewstring[2])

        # Add the bank nodes
        for nodedata in vewstring[1:-1]:
            node = nodedata['node_id']
            if node in map_node.keys():  # Skip if node already added (for closed strings)
                continue
            
            # Validate that the node exists in the mesh
            if node not in self._mesh.nodes.index:
                raise ValueError(f"Node ID {node} from VEW string not found in mesh. "
                               f"Please verify that the VEW string file contains valid node IDs that exist in the mesh.")
            
            nn += 1
            map_node.update({node: nn})
            id_new.append(nn)
            x_new.append(self._mesh.x[node])
            y_new.append(self._mesh.y[node])
            zval = nodedata['bank_elevation']
            if zval == -99999.0:
                zval = self._mesh.nodes.loc[node, 'value_1']
            z_new.append(zval)

        # Update nodes
        nodes_update_time = time.time()
        nodes_new = pd.concat([
            self._mesh.nodes,
            pd.DataFrame(index=id_new, data={'x': x_new, 'y': y_new, 'value_1': z_new})
        ])

        # Update elements
        elements_start_time = time.time()
        elements = self._mesh.elements.elements
        elements_new = elements.copy()
        nodestring = [vewstring[i]['node_id'] for i in range(len(vewstring))]
        node_elements = self._mesh.node_elements.copy()

        # Detach elements along the vewstring
        for i in range(1, len(nodestring)-1):
            node1 = nodestring[i-1]
            node2 = nodestring[i]
            node3 = nodestring[i+1]
            
            print(f"\n[DEBUG] Processing VEW segment {i}/{len(nodestring)-2}: node1={node1}, node2={node2}, node3={node3}")
            
            # Validate that all nodes exist in the mesh
            for node_idx, node in enumerate([node1, node2, node3], start=i-1):
                if node not in self._mesh.nodes.index:
                    raise ValueError(f"Node ID {node} at position {node_idx} in VEW string not found in mesh.")
            
            # Validate that node2 has associated elements
            if node2 not in node_elements:
                raise ValueError(f"Node ID {node2} has no associated elements in the mesh. "
                               f"This may indicate a problem with the mesh connectivity.")

            # Find elements on the right side of the line segment
            eids = [int(node) for node in node_elements[node2]]
            
            print(f"[DEBUG] Node {node2} has {len(eids)} associated elements: {eids[:10]}{'...' if len(eids) > 10 else ''}")
            
            if not eids:
                raise ValueError(f"Node ID {node2} has no elements associated with it. "
                               f"Cannot process VEW string through this node.")

            # Find elements counterclockwise from node1 to node3
            # Special case: if node1 == node3, there's nothing to traverse
            if node1 == node3:
                print(f"[DEBUG] WARNING: node1 ({node1}) equals node3 ({node3}) - skipping traversal (no elements to process)")
                continue
            
            # Special case: if node2 == node3, we're already at the target node
            # In this case, find all elements connected to node2 that contain node1
            # These are the elements on the "right side" of the line from node1 to node2
            if node2 == node3:
                print(f"[DEBUG] WARNING: node2 ({node2}) equals node3 ({node3}) - collecting elements directly (no traversal needed)")
                eids_right = []
                for eid in eids:
                    if eid not in elements.index:
                        continue
                    try:
                        elem_nodes = [int(e) for e in elements.loc[eid].to_list()[1:4]]
                    except Exception as e:
                        continue
                    
                    # Check if element contains both node1 and node2
                    if node1 in elem_nodes and node2 in elem_nodes:
                        # Verify correct orientation: node2 should be before node1 counterclockwise
                        node1_idx = elem_nodes.index(node1)
                        node2_idx = elem_nodes.index(node2)
                        # Check if node2 is the previous node before node1 (counterclockwise)
                        if elem_nodes[self._map_elem_node_prev[node1_idx]] == node2:
                            eids_right.append(eid)
                            print(f"[DEBUG]   Found element {eid} with nodes {elem_nodes} (node2 before node1 counterclockwise)")
                
                print(f"[DEBUG] ✓ Found {len(eids_right)} elements on right side: {eids_right}")
            else:
                # Normal case: traverse from node1 to node3 through node2
                curr_node = node1
                eids_right = []
                max_iterations = len(eids) * 10  # Safety limit to prevent infinite loops
                iteration_count = 0
                visited_nodes = set()  # Track visited nodes to detect cycles
                
                print(f"[DEBUG] Starting traversal from node1={node1} to node3={node3} through node2={node2}")
                print(f"[DEBUG] Initial visited_nodes: {visited_nodes}")
                
                while curr_node != node3:
                    iteration_count += 1
                    
                    print(f"[DEBUG] Iteration {iteration_count}: curr_node={curr_node}, visited_nodes={visited_nodes}")
                    
                    # Safety check: prevent infinite loops
                    if iteration_count > max_iterations:
                        print(f"[DEBUG] ERROR: Max iterations ({max_iterations}) exceeded")
                        print(f"[DEBUG] Final visited_nodes: {visited_nodes}")
                        print(f"[DEBUG] Final eids_right: {eids_right}")
                        raise RuntimeError(f"Infinite loop detected while processing VEW string. "
                                         f"Cannot traverse from node {node1} to node {node3} through node {node2}. "
                                         f"This may indicate corrupted mesh connectivity or invalid VEW string definition.")
                    
                    # Safety check: detect circular paths
                    if curr_node in visited_nodes:
                        print(f"[DEBUG] ERROR: Circular path detected!")
                        print(f"[DEBUG] Current node {curr_node} was already in visited_nodes: {visited_nodes}")
                        print(f"[DEBUG] Traversal path so far: {sorted(visited_nodes)}")
                        print(f"[DEBUG] Target node3: {node3}")
                        print(f"[DEBUG] Starting node1: {node1}")
                        print(f"[DEBUG] Central node2: {node2}")
                        print(f"[DEBUG] Elements found so far: {eids_right}")
                        
                        # Additional debug: check if node1 == node3
                        if node1 == node3:
                            print(f"[DEBUG] WARNING: node1 ({node1}) equals node3 ({node3}) - this might be the issue!")
                        
                        raise RuntimeError(f"Circular path detected while traversing from node {node1} to node {node3}. "
                                         f"Current node {curr_node} has been visited before. "
                                         f"This may indicate invalid mesh topology or VEW string definition.")
                    
                    visited_nodes.add(curr_node)
                    
                    found_element = False
                    elements_checked = []
                    for eid in eids:
                        # Safety check: validate element existence
                        if eid not in elements.index:
                            print(f"[DEBUG] Warning: Element ID {eid} not found in mesh elements. Skipping.")
                            continue
                            
                        try:
                            elem_nodes = [int(e) for e in elements.loc[eid].to_list()[1:4]]
                        except Exception as e:
                            print(f"[DEBUG] Warning: Error reading element {eid}: {e}. Skipping.")
                            continue
                        
                        elements_checked.append((eid, elem_nodes))
                        
                        if curr_node in elem_nodes:
                            curr_index = elem_nodes.index(curr_node)
                            prev_node = elem_nodes[self._map_elem_node_prev[curr_index]]
                            next_node = elem_nodes[self._map_elem_node_next[curr_index]]
                            print(f"[DEBUG]   Checking element {eid}: nodes={elem_nodes}, curr_node at index {curr_index}, prev_node={prev_node}, next_node={next_node}, node2={node2}")
                            
                            if prev_node == node2:
                                found_element = True
                                print(f"[DEBUG]   ✓ Found matching element {eid}: will move from {curr_node} to {next_node}")
                                break
                    
                    if not found_element:
                        print(f"[DEBUG] ERROR: No valid element found for curr_node={curr_node}")
                        print(f"[DEBUG] Elements checked: {elements_checked[:5]}{'...' if len(elements_checked) > 5 else ''}")
                        raise RuntimeError(f"Cannot find valid element containing current node {curr_node} "
                                         f"while traversing from node {node1} to node {node3}. "
                                         f"This may indicate disconnected mesh regions or invalid VEW string.")
                    
                    eids_right.append(eid)
                    curr_node = elem_nodes[self._map_elem_node_next[curr_index]]
                    print(f"[DEBUG]   Moving to next node: {curr_node}")
                
                print(f"[DEBUG] ✓ Traversal completed: reached node3={node3}")
                print(f"[DEBUG] Final visited_nodes: {sorted(visited_nodes)}")
                print(f"[DEBUG] Elements found: {len(eids_right)} elements: {eids_right[:10]}{'...' if len(eids_right) > 10 else ''}")

            for eid in eids_right:
                # Safety check before modifying elements
                if eid in elements.index:
                    elements_new.loc[eid, elements.columns[elements.loc[eid] == node2]] = map_node[node2]

        # Update boundaries
        boundary_start_time = time.time()
        boundaries = self._mesh.boundaries.to_dict()
        boundaries_new = boundaries.copy()

        # Add VEW boundary
        vew_node_id = [(str(map_node[nodestring[i]]), str(nodestring[i])) 
                       for i in range(1, len(nodestring)-1)]
        vew_barrier_height = [vewstring[i]['bank_elevation'] + 1e-3 
                            for i in range(1, len(nodestring)-1)]
        vew_subcritical_flow_coefficient = [1.0 for _ in range(1, len(nodestring)-1)]
        vew_supercritical_flow_coefficient = [1.0 for _ in range(1, len(nodestring)-1)]

        vewboundary = {
            'node_id': vew_node_id,
            'barrier_height': vew_barrier_height,
            'subcritical_flow_coefficient': vew_subcritical_flow_coefficient,
            'supercritical_flow_coefficient': vew_supercritical_flow_coefficient
        }

        vewboundaries = boundaries_new.get('64', [])
        vewboundaries.append(vewboundary)
        boundaries_new['64'] = vewboundaries

        total_time = time.time() - start_time

        # Return updated mesh
        return AdcircMesh(nodes=nodes_new, elements=elements_new, boundaries=boundaries_new)

    def add_vew_strings(self, vewstrings: List[List[Dict]]) -> AdcircMesh:
        """Add multiple VEW strings to the mesh.
        
        Args:
            vewstrings: List of VEW string definitions
            
        Returns:
            Modified AdcircMesh object
        """
        # In coordinate mode, we need to convert coordinates to node IDs first for validation
        if self._coordinate_mode:
            print("Converting coordinates to node IDs for validation...")
            processed_vewstrings = []
            for string_idx, vewstring in enumerate(vewstrings, 1):
                processed_vewstring = []
                for i, nodedata in enumerate(vewstring):
                    try:
                        node_id, bank_elevation = self._process_vew_node_data(nodedata)
                        processed_nodedata = {
                            'node_id': node_id,
                            'bank_elevation': bank_elevation
                        }
                        processed_vewstring.append(processed_nodedata)
                    except ValueError as e:
                        node_identifier = self._get_node_identifier(nodedata)
                        raise ValueError(f"Error processing VEW string {string_idx}, node {i+1} {node_identifier}: {e}")
                processed_vewstrings.append(processed_vewstring)
            
            # Use processed strings for validation
            validation_strings = processed_vewstrings
        else:
            validation_strings = vewstrings

        # Validate for consecutive identical nodes within each VEW string
        print("Validating VEW strings for consecutive identical nodes...")
        consecutive_errors = []
        for string_idx, vewstring in enumerate(validation_strings, 1):
            if not vewstring:
                continue
            nodestring = [vewstring[i]['node_id'] for i in range(len(vewstring))]
            consecutive_positions = []
            for j in range(1, len(nodestring)):
                if nodestring[j] == nodestring[j-1]:
                    consecutive_positions.append((j, nodestring[j]))
            if consecutive_positions:
                consecutive_errors.append((string_idx, consecutive_positions))
        if consecutive_errors:
            error_msg = "Invalid VEW strings: consecutive identical nodes detected.\n"
            for string_idx, positions in consecutive_errors:
                details = ", ".join([f"(positions {pos-1}-{pos}, node {node_id})" for pos, node_id in positions])
                error_msg += f"  • VEW string {string_idx}: {details}\n"
            error_msg += "\nConsecutive identical nodes create degenerate segments (e.g., ... A, A, ...)\n"
            error_msg += "which lead to traversal loops and ambiguous topology modifications.\n\n"
            error_msg += "SOLUTION: Remove duplicate adjacent entries in the VEW string so that no\n"
            error_msg += "two neighboring nodes are the same. If a closed string is desired, only the\n"
            error_msg += "first and last nodes should be identical, not adjacent interior nodes."
            raise ValueError(error_msg)

        # Validate for duplicate nodes between VEW strings
        print("Validating VEW strings for duplicate nodes...")
        all_processed_nodes = set()
        duplicate_conflicts = []
        
        for string_idx, vewstring in enumerate(validation_strings, 1):
            # Extract node IDs from this string (excluding padding nodes)
            nodestring = [vewstring[i]['node_id'] for i in range(len(vewstring))]
            
            # For closed strings, remove the duplicate end node
            if len(nodestring) > 1 and nodestring[0] == nodestring[-1]:
                nodestring = nodestring[:-1]  # Remove duplicate end node
            
            # Check for duplicates with previously processed strings
            for node_id in nodestring:
                if node_id in all_processed_nodes:
                    duplicate_conflicts.append((node_id, string_idx))
                
            # Add nodes from this string to the processed set
            all_processed_nodes.update(nodestring)
        
        # Report any duplicate conflicts
        if duplicate_conflicts:
            error_msg = "CRITICAL ERROR: Duplicate nodes found between different VEW strings:\n"
            node_to_strings = {}
            
            # Group conflicts by node
            for node_id, string_idx in duplicate_conflicts:
                if node_id not in node_to_strings:
                    node_to_strings[node_id] = []
                node_to_strings[node_id].append(string_idx)
            
            # Find which strings contain each duplicate node
            for node_id, conflicting_strings in node_to_strings.items():
                # Find all strings containing this node
                all_strings_with_node = []
                for string_idx, vewstring in enumerate(validation_strings, 1):
                    nodestring = [vewstring[i]['node_id'] for i in range(len(vewstring))]
                    if node_id in nodestring:
                        all_strings_with_node.append(string_idx)
                
                error_msg += f"  • Node {node_id} appears in VEW strings: {all_strings_with_node}\n"
            
            error_msg += "\nThis will cause infinite loops during processing because the mesh connectivity "
            error_msg += "around these nodes gets modified multiple times.\n\n"
            error_msg += "SOLUTION: Remove the duplicate nodes from all but one of the conflicting VEW strings."
            
            raise ValueError(error_msg)
        
        print(f"✓ All {len(vewstrings)} VEW strings validated - no duplicate nodes between strings")
        
        # Process each VEW string
        mesh = self._mesh
        for string_idx, vewstring in enumerate(vewstrings, 1):
            print(f"Processing VEW string {string_idx}/{len(vewstrings)}...")
            mesh = self.add_vew_string(vewstring)
            self._mesh = mesh  # Update internal mesh for next iteration
        return mesh


def add_vews_to_mesh(f14file: str, vewfile: str, output_f14: str = None, 
                     coordinate_mode: bool = False, tolerance: float = 1e-6) -> None:
    """Add VEW boundaries to an ADCIRC mesh based on VEW string definitions.
    
    Args:
        f14file: Path to input fort.14 file
        vewfile: Path to input YAML file containing VEW string definitions
        output_f14: Path to output fort.14 file (optional)
        coordinate_mode: If True, use coordinate-based node matching instead of node IDs
        tolerance: Distance tolerance for coordinate matching (in mesh units)
    """
    start_time = time.time()
    
    # Read mesh and VEW strings
    print(f"Reading mesh from {f14file}...")
    mesh_read_time = time.time()
    mesh = AdcircMesh.open(f14file)
    
    print(f"Reading VEW definitions from {vewfile}...")
    vew_read_time = time.time()
    with open(vewfile, 'r') as file:
        vewdata = yaml.safe_load(file)
    vewstrings = vewdata['vewstrings']
    print(f"Found {len(vewstrings)} VEW strings to process")
    
    # Print mode information
    if coordinate_mode:
        print(f"Using coordinate mode with tolerance: {tolerance}")
    else:
        print("Using node ID mode")

    # Add VEW boundaries
    processing_time = time.time()
    adder = VEWBoundaryAdder(mesh, coordinate_mode=coordinate_mode, tolerance=tolerance)
    mesh_new = adder.add_vew_strings(vewstrings)
    mesh_new.description = "Generated by vew_boundary_adder"

    # Write output
    if output_f14 is None:
        output_f14 = str(Path(f14file).with_suffix('')) + '_vew.grd'
    
    print(f"Writing output mesh to {output_f14}...")
    write_time = time.time()
    mesh_new.write(output_f14, overwrite=True)
    
    total_time = time.time() - start_time
    print(f"✓ Total processing time: {total_time:.3f} seconds")


def get_parser():
    parser = argparse.ArgumentParser(add_help=False, description="Add VEW boundaries to an ADCIRC mesh")
    parser.add_argument("f14file", help="Input fort.14 file")
    parser.add_argument("vewfile", help="Input YAML file containing VEW string definitions")
    parser.add_argument("-o", "--output", help="Output fort.14 file")
    parser.add_argument("-c", "--coordinate-mode", action="store_true", 
                       help="Use coordinate-based node matching instead of node IDs")
    parser.add_argument("-t", "--tolerance", type=float, default=1e-6,
                       help="Distance tolerance for coordinate matching (default: 1e-6)")
    return parser

def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    add_vews_to_mesh(args.f14file, args.vewfile, args.output, 
                     coordinate_mode=args.coordinate_mode, tolerance=args.tolerance)


if __name__ == "__main__":
    main() 