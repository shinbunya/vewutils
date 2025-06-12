#!/usr/bin/env python3
"""
Add VEW boundaries to an ADCIRC mesh based on VEW string definitions in a YAML file.

This program reads an ADCIRC mesh in fort.14 format and VEW string definitions from a YAML file,
then adds VEW boundaries to the mesh. The modified mesh is written out.
"""

import argparse
import yaml
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple
from adcircpy import AdcircMesh
import time


class VEWBoundaryAdder:
    """Class for adding VEW boundaries to an ADCIRC mesh."""

    def __init__(self, mesh: AdcircMesh):
        """Initialize with an ADCIRC mesh."""
        self._mesh = mesh
        self._map_elem_node_prev = {0: 2, 1: 0, 2: 1}
        self._map_elem_node_next = {0: 1, 1: 2, 2: 0}

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

        if len(vewstring) < 4 or (len(vewstring) < 3 and vewstring[0]['node_id'] == vewstring[-1]['node_id']):
            raise ValueError(
                "The length of a vewstring is {:d}. It should be greater than 2 for an open node string "
                "and should be greater than 2 for a closed node string.".format(len(vewstring)))

        # Add the second last node to the beginning and the second node to the end
        if vewstring[0]['node_id'] == vewstring[-1]['node_id']:
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
            
            if not eids:
                raise ValueError(f"Node ID {node2} has no elements associated with it. "
                               f"Cannot process VEW string through this node.")

            # Find elements counterclockwise from node1 to node3
            curr_node = node1
            eids_right = []
            max_iterations = len(eids) * 10  # Safety limit to prevent infinite loops
            iteration_count = 0
            visited_nodes = set()  # Track visited nodes to detect cycles
            
            while curr_node != node3:
                iteration_count += 1
                
                # Safety check: prevent infinite loops
                if iteration_count > max_iterations:
                    raise RuntimeError(f"Infinite loop detected while processing VEW string. "
                                     f"Cannot traverse from node {node1} to node {node3} through node {node2}. "
                                     f"This may indicate corrupted mesh connectivity or invalid VEW string definition.")
                
                # Safety check: detect circular paths
                if curr_node in visited_nodes:
                    raise RuntimeError(f"Circular path detected while traversing from node {node1} to node {node3}. "
                                     f"Current node {curr_node} has been visited before. "
                                     f"This may indicate invalid mesh topology or VEW string definition.")
                
                visited_nodes.add(curr_node)
                
                found_element = False
                for eid in eids:
                    # Safety check: validate element existence
                    if eid not in elements.index:
                        print(f"Warning: Element ID {eid} not found in mesh elements. Skipping.")
                        continue
                        
                    try:
                        elem_nodes = [int(e) for e in elements.loc[eid].to_list()[1:4]]
                    except Exception as e:
                        print(f"Warning: Error reading element {eid}: {e}. Skipping.")
                        continue
                        
                    if curr_node in elem_nodes:
                        curr_index = elem_nodes.index(curr_node)
                        if elem_nodes[self._map_elem_node_prev[curr_index]] == node2:
                            found_element = True
                            break
                
                if not found_element:
                    raise RuntimeError(f"Cannot find valid element containing current node {curr_node} "
                                     f"while traversing from node {node1} to node {node3}. "
                                     f"This may indicate disconnected mesh regions or invalid VEW string.")
                
                eids_right.append(eid)
                curr_node = elem_nodes[self._map_elem_node_next[curr_index]]

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
        # Validate for duplicate nodes between VEW strings
        print("Validating VEW strings for duplicate nodes...")
        all_processed_nodes = set()
        duplicate_conflicts = []
        
        for string_idx, vewstring in enumerate(vewstrings, 1):
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
                for string_idx, vewstring in enumerate(vewstrings, 1):
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


def add_vews_to_mesh(f14file: str, vewfile: str, output_f14: str = None) -> None:
    """Add VEW boundaries to an ADCIRC mesh based on VEW string definitions."""
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

    # Add VEW boundaries
    processing_time = time.time()
    adder = VEWBoundaryAdder(mesh)
    mesh_new = adder.add_vew_strings(vewstrings)
    mesh_new.description = "Generated by vew_boundary_adder"

    # Write output
    if output_f14 is None:
        output_f14 = str(Path(f14file).with_suffix('')) + '_vew.grd'
    
    print(f"Writing output mesh to {output_f14}...")
    write_time = time.time()
    mesh_new.write(output_f14, overwrite=True)
    
    total_time = time.time() - start_time


def get_parser():
    parser = argparse.ArgumentParser(add_help=False, description="Add VEW boundaries to an ADCIRC mesh")
    parser.add_argument("f14file", help="Input fort.14 file")
    parser.add_argument("vewfile", help="Input YAML file containing VEW string definitions")
    parser.add_argument("-o", "--output", help="Output fort.14 file")
    return parser

def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    add_vews_to_mesh(args.f14file, args.vewfile, args.output)


if __name__ == "__main__":
    main() 