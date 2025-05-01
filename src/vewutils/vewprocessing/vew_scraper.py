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

    def _strip_vewstring(self, mesh: AdcircMesh, ivewboundary: int) -> Tuple[AdcircMesh, List[Dict]]:
        """Strip a single VEW boundary from the mesh and return the modified mesh and VEW string."""
        boundaries = mesh.boundaries.to_dict().copy()
        vewboundaries = boundaries['64']
        
        node_elements = mesh.node_elements
        node_neighbors = mesh.node_neighbors

        vewboundary = vewboundaries[ivewboundary]

        bank_nodes = [int(vewboundary['node_id'][i][0]) for i in range(len(vewboundary['node_id']))]
        channel_nodes = [int(vewboundary['node_id'][i][1]) for i in range(len(vewboundary['node_id']))]

        # Making a new mesh without the vew boundary
        nodes = mesh.nodes.copy()
        nodes_new = mesh.nodes.copy()
        elements_new = mesh.elements.elements.copy()
        boundaries_new = boundaries.copy()

        # Create a mapping table between the old and new node numbers
        idx = 1
        map_node = {}
        for i in nodes_new.index:
            if i in bank_nodes:
                idx_bank = bank_nodes.index(i)
                map_node[i] = channel_nodes[idx_bank]
            else:
                map_node[i] = idx
                idx += 1

        # Update node ids in the element table
        elements_new[['node_1', 'node_2', 'node_3']] = elements_new[['node_1', 'node_2', 'node_3']].replace(map_node)

        # Drop the bank nodes in the node table
        nodes_new = nodes_new.drop(index=bank_nodes)
        nodes_new.index = nodes_new.index.map(map_node)

        # Drop the targeted vew boundary
        boundaries_new['64'].pop(ivewboundary)

        # Update the node numbers in the boundary table
        for ibtype in boundaries.keys():
            for i in range(len(boundaries[ibtype])):
                if ibtype is None:
                    dfbnd = pd.DataFrame(boundaries[ibtype][i]['node_id'])
                    dfbnd.replace(map_node, inplace=True)
                    boundaries_new[ibtype][i]['node_id'] = dfbnd.iloc[:,0].values
                elif ibtype.endswith('4'):
                    dfbnd = pd.DataFrame(boundaries[ibtype][i]['node_id'])
                    dfbnd = dfbnd.replace(map_node)
                    boundaries_new[ibtype][i]['node_id'] = list(zip(dfbnd.iloc[:,0].values, dfbnd.iloc[:,1].values))
                else:
                    dfbnd = pd.DataFrame(boundaries[ibtype][i]['node_id'])
                    dfbnd.replace(map_node, inplace=True)
                    boundaries_new[ibtype][i]['node_id'] = dfbnd.values.tolist()
                    
        new_mesh = AdcircMesh(nodes=nodes_new, elements=elements_new, boundaries=boundaries_new)

        triangles = mesh.elements.elements.iloc[:, 1:4]

        # Making a vewstring
        # Reverse node strings if it's not seeing the bank side on the right
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
            
        # Create a vewstring
        vewstring = []
        for i in range(len(bank_nodes)):
            vewstring.append({
                'node_id': int(map_node[channel_nodes[i]]),
                'x': float(nodes.loc[bank_nodes[i], 'x']),
                'y': float(nodes.loc[channel_nodes[i], 'y']),
                'bank_elevation': float(nodes.loc[bank_nodes[i], 'value_1']),
                'bank_mannings_n': 0.03
            })
               
        return new_mesh, vewstring

    def strip_vewstrings(self) -> Tuple[AdcircMesh, Dict]:
        """Strip all VEW boundaries from the mesh and return the modified mesh and VEW strings."""
        boundaries = self._mesh.boundaries.to_dict().copy()
        vewboundaries = boundaries['64']

        vewstrings = []
        mesh = self._mesh
        for ivewboundary in range(len(vewboundaries)):
            mesh, vewstring = self._strip_vewstring(mesh, 0)
            vewstrings.append(vewstring)
        
        vewstrings = {'vewstrings': vewstrings}

        return mesh, vewstrings


def main():
    """Main function to handle command line arguments and process the mesh."""
    parser = argparse.ArgumentParser(
        description="Scrape VEW boundaries from an ADCIRC mesh and save them to YAML format"
    )
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
    
    args = parser.parse_args()
    
    # Read the mesh file
    mesh = AdcircMesh.open(args.input_mesh)
    
    # Strip VEW boundaries
    scraper = VEWScraper(mesh)
    mesh_new, vewstrings = scraper.strip_vewstrings()
    mesh_new.description = args.description
    
    # Save the new mesh file
    mesh_new.write(args.output_mesh, overwrite=True)
    
    # Write the YAML output to a file
    with open(args.output_yaml, 'w') as f:
        yaml.dump(vewstrings, f, sort_keys=False)
    
    return 0


if __name__ == "__main__":
    main() 