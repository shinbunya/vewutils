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


class VEWBoundaryAdder:
    """Class for adding VEW boundaries to an ADCIRC mesh."""

    def __init__(self, mesh: AdcircMesh):
        """Initialize with an ADCIRC mesh."""
        self._mesh = mesh
        self._map_elem_node_prev = {0: 2, 1: 0, 2: 1}
        self._map_elem_node_next = {0: 1, 1: 2, 2: 0}

    def add_vew_string(self, vewstring: List[Dict]) -> AdcircMesh:
        """Add a single VEW string to the mesh.
        
        Args:
            vewstring: List of dictionaries containing node information for the VEW string
            
        Returns:
            Modified AdcircMesh object
        """
        # Create bank nodes and add them to the mesh
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
        # if it is the same as the last node. This ensures the beginning node is also duplicated.
        if vewstring[0]['node_id'] == vewstring[-1]['node_id']:
            vewstring.insert(0, vewstring[-2])
            vewstring.append(vewstring[2])

        # Add the bank nodes
        for nodedata in vewstring[1:-1]:
            node = nodedata['node_id']
            if node in map_node.keys():  # Skip if node already added (for closed strings)
                continue
            nn += 1
            map_node.update({node: nn})
            id_new.append(nn)
            x_new.append(self._mesh.x[node])
            y_new.append(self._mesh.y[node])
            z_new.append(nodedata['bank_elevation'])

        nodes_new = pd.concat([
            self._mesh.nodes,
            pd.DataFrame(index=id_new, data={'x': x_new, 'y': y_new, 'value_1': z_new})
        ])

        # Update elements
        elements = self._mesh.elements.elements
        elements_new = elements.copy()
        nodestring = [vewstring[i]['node_id'] for i in range(len(vewstring))]
        node_elements = self._mesh.node_elements.copy()

        # Detach elements along the vewstring
        for i in range(1, len(nodestring)-1):
            node1 = nodestring[i-1]
            node2 = nodestring[i]
            node3 = nodestring[i+1]

            # Find elements on the right side of the line segment
            eids = [int(node) + 1 for node in node_elements[node2-1]]

            # Find elements counterclockwise from node1 to node3
            curr_node = node1
            eids_right = []
            while curr_node != node3:
                for eid in eids:
                    elem_nodes = [int(e) for e in elements.loc[eid].to_list()[1:4]]
                    if curr_node in elem_nodes:
                        curr_index = elem_nodes.index(curr_node)
                        if elem_nodes[self._map_elem_node_prev[curr_index]] == node2:
                            break
                eids_right.append(eid)
                curr_node = elem_nodes[self._map_elem_node_next[curr_index]]

            for eid in eids_right:
                elements_new.loc[eid, elements.columns[elements.loc[eid] == node2]] = map_node[node2]

        # Update boundaries
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

        # Return updated mesh
        return AdcircMesh(nodes=nodes_new, elements=elements_new, boundaries=boundaries_new)

    def add_vew_strings(self, vewstrings: List[List[Dict]]) -> AdcircMesh:
        """Add multiple VEW strings to the mesh.
        
        Args:
            vewstrings: List of VEW string definitions
            
        Returns:
            Modified AdcircMesh object
        """
        mesh = self._mesh
        for vewstring in vewstrings:
            mesh = self.add_vew_string(vewstring)
            self._mesh = mesh  # Update internal mesh for next iteration
        return mesh


def add_vews_to_mesh(f14file: str, vewfile: str, output_f14: str = None) -> None:
    """Add VEW boundaries to an ADCIRC mesh based on VEW string definitions.
    
    Args:
        f14file: Path to input fort.14 file
        vewfile: Path to input YAML file containing VEW string definitions
        output_f14: Path to output fort.14 file. If None, will append '_vew' to input filename.
    """
    # Read mesh and VEW strings
    mesh = AdcircMesh.open(f14file)
    with open(vewfile, 'r') as file:
        vewdata = yaml.safe_load(file)
    vewstrings = vewdata['vewstrings']

    # Add VEW boundaries
    adder = VEWBoundaryAdder(mesh)
    mesh_new = adder.add_vew_strings(vewstrings)
    mesh_new.description = "Generated by vew_boundary_adder"

    # Write output
    if output_f14 is None:
        output_f14 = str(Path(f14file).with_suffix('')) + '_vew.grd'
    
    mesh_new.write(output_f14, overwrite=True)
    print(f"Wrote mesh with VEW boundaries to {output_f14}")


def main():
    parser = argparse.ArgumentParser(description="Add VEW boundaries to an ADCIRC mesh")
    parser.add_argument("f14file", help="Input fort.14 file")
    parser.add_argument("vewfile", help="Input YAML file containing VEW string definitions")
    parser.add_argument("-o", "--output", help="Output fort.14 file")
    args = parser.parse_args()

    add_vews_to_mesh(args.f14file, args.vewfile, args.output)


if __name__ == "__main__":
    main() 