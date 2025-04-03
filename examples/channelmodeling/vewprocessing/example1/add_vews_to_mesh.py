# %%
import os
import yaml
import numpy as np
import pandas as pd
from adcircpy import AdcircMesh
from adcircpy.mesh.base import Elements

# %%
# Paths to the input/output files
f14file = "/home/sbunya/GitHub/adcircutils/adcircutils/vewchannel/examples/example1/fort.14"
vewfile = "/home/sbunya/GitHub/adcircutils/adcircutils/vewchannel/examples/example1/vewstings.yaml"

# %%
# Load the YAML file
with open(vewfile, 'r') as file:
    vewdata = yaml.safe_load(file)
vewstrings = vewdata['vewstrings']

# %%
# Load the mesh file
mesh = AdcircMesh.open(f14file)

# %%
# Creating a new mesh with bank nodes
def add_vewstring(mesh, vewstring) -> AdcircMesh:
    # Dictionary to store mapping of node orders in an element
    map_elem_node_prev = {0: 2, 1: 0, 2: 1}
    map_elem_node_next = {0: 1, 1: 2, 2: 0}

    # Create bank nodes and add them to the mesh
    nn = mesh.nodes.shape[0]
    nodemap = {}
    id_new = []
    x_new = []
    y_new = []
    z_new = []

    if len(vewstring) < 3:
        raise ValueError("The length of vewstring {:d} is {:d}. It should be greater than 2.".format(istring, len(vewstring)))

    for nodedata in vewstring:
        node = nodedata['node']
        nn += 1
        nodemap.update({node: nn})
        id_new.append(nn)
        x_new.append(mesh.x[node])
        y_new.append(mesh.y[node])
        z_new.append(nodedata['bank_elevation'])

    nodes_new = pd.concat([mesh.nodes, pd.DataFrame(index=id_new, data={'x': x_new, 'y': y_new, 'value_1': z_new})])

    # Update elements
    elements = mesh.elements.elements
    elements_new = elements.copy()
    nodestring = [vewstring[i]['node'] for i in range(len(vewstring))]
    node_elements = mesh.node_elements
    
    # Add the second last node to the beginning if it is the same as the last node
    # This ensures the beginning node is also duplicated.
    if nodestring[0] == nodestring[-1]:
        nodestring.insert(0, nodestring[-2])

    # Detaching the elements along the vewstring
    for i in range(1, len(nodestring)-1):
        node1 = nodestring[i-1]
        node2 = nodestring[i]
        node3 = nodestring[i+1]

        # Find the elements that are on the right side of the line segment that is defined by node1, node2, and node3
        eids = node_elements[node2]

        # Starting from the element that contains node1, find the elements counterclose-wise until it reaches node3
        curr_node = node1
        eids_right = []
        while curr_node != node3:
            for eid in eids:
                elem_nodes = [int(e) for e in elements.loc[eid].to_list()[1:4]]
                if curr_node in elem_nodes:
                    curr_index = elem_nodes.index(curr_node)
                    if elem_nodes[map_elem_node_prev[curr_index]] == node2:
                        break
            eids_right.append(eid)
            curr_node = elem_nodes[map_elem_node_next[curr_index]]
        
        for eid in eids_right:
            elements_new.loc[eid, elements.columns[elements.loc[eid] == node2]] = nodemap[node2]

    # Add vew boundaries
    

    # Return the updated mesh   
    return AdcircMesh(nodes=nodes_new, elements=elements_new)

def add_vewstrings(mesh, vewstrings) -> AdcircMesh:
    for vewstring in vewstrings:
        mesh = add_vewstring(mesh, vewstring)
    return mesh

# Load the mesh file
mesh = AdcircMesh.open(f14file)
mesh_new = add_vewstrings(mesh, vewstrings)

# %%
mesh.make_plot(vmin=-4, vmax=1.5)
mesh_new.make_plot(vmin=-4, vmax=1.5)

# %%



