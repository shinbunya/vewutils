# %%
import os
import yaml
import numpy as np
import geopandas as gpd
from shapely.geometry import Point, LineString
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from adcircpy import AdcircMesh

# %%
# Paths to the input/output files
f14file = "/home/sbunya/GitHub/adcircutils/adcircutils/vewchannel/examples/example1/fort.14"
vewlocation_geofile = "/home/sbunya/GitHub/adcircutils/adcircutils/vewchannel/examples/example1/example1_vew_locations.shp"
vewfile = "/home/sbunya/GitHub/adcircutils/adcircutils/vewchannel/examples/example1/vewstings.yaml"

# Parameters
dist_max = 10.0  # Maximum distance for nearest neighbor search in meters

# %%
# Read the shapefile
vewlocation_gdf = gpd.read_file(vewlocation_geofile)
print(vewlocation_gdf)

# %%
# Load the mesh file
mesh = AdcircMesh.open(f14file)

# %%
# Create nodestrings along the lines in the geospatial file
x = mesh.nodes.x.values
y = mesh.nodes.y.values
tri = mesh.triangles - 1

# Create a KDTree for fast nearest neighbor search
tree = cKDTree(np.c_[x, y])

neighs = mesh.node_neighbors
neighs = [list(neighs[i]) for i in range(len(neighs))]

nodestrings = []

for iline, line in enumerate(vewlocation_gdf.geometry):
    slx, sly = line.coords[0][0], line.coords[0][1]
    elx, ely = line.coords[-1][0], line.coords[-1][1]
    
    # Query the KDTree for the nearest neighbor
    distance, nearest_node_index = tree.query([slx, sly])

    # Get the coordinates of the nearest node
    nearest_node_x = x[nearest_node_index]
    nearest_node_y = y[nearest_node_index]

    ni = nearest_node_index
    xi = nearest_node_x
    yi = nearest_node_y
    distance_to_endpoint = 1e10
    
    nodestring = [nearest_node_index]
    
    ipos = 0
    while ipos < len(line.coords) - 1:
        # Get the points in the line segment
        point0 = line.coords[ipos]
        point1 = line.coords[ipos + 1]
        xl0, yl0 = point0[0], point0[1]
        xl1, yl1 = point1[0], point1[1]

        line_segment = LineString([(xl0, yl0), (xl1, yl1)])
        
        neigh = neighs[ni]
        neigh = [n for n in neigh if n not in nodestring]

        xnei = [x[n] for n in neigh]
        ynei = [y[n] for n in neigh]
        
        # Compute the distances from the line to the points defined by xnei and ynei
        distances = [line_segment.distance(Point(xn, yn)) for xn, yn in zip(xnei, ynei)]
        min_distance = np.min(distances)
        
        if min_distance > dist_max:
            ipos += 1
            continue
        
        min_distance_index = np.argmin(distances)
        ni = neigh[min_distance_index]
        xi = x[ni]
        yi = y[ni]
        
        nodestring.append(int(ni))
        
        distance_to_endpoint = np.sqrt((elx - xi)**2 + (ely - yi)**2)
        
        if distance_to_endpoint < dist_max:
            break

    xi, yi = x[nodestring[0]], y[nodestring[0]]
    distance_to_endpoint = np.sqrt((elx - xi)**2 + (ely - yi)**2)
    if distance_to_endpoint < dist_max:
        nodestring.append(int(nodestring[0]))
        
    nodestrings.append(nodestring)

# %%
# Visualize the nodestrings
fig, ax = plt.subplots()
for nodestring in nodestrings:
    x_nodes = [x[n] for n in nodestring]
    y_nodes = [y[n] for n in nodestring]
    ax.plot(x_nodes, y_nodes, marker='o')
    print(x_nodes)
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')

# %%
# Assign values to each node in the nodestrings
nodestrings_with_values = []
for nodestring in nodestrings:
    nodestring_with_value = \
        [{'node': node+1, 
          'x': float(mesh.nodes.x.loc[node]), 'y': float(mesh.nodes.y.loc[node]),
          'bank_elevation': 1.5 if mesh.nodes.y.loc[node] > 500.0 else -4.0,
          'bank_mannings_n': 0.03,
          } for node in nodestring]
    nodestrings_with_values.append(nodestring_with_value)

# Prepare the data for YAML output
data = {'vewstrings': nodestrings_with_values}

# Print the YAML output
print(yaml.dump(data, sort_keys=False))

# Write the YAML output to a file
with open(vewfile, 'w') as f:
    yaml.dump(data, f, sort_keys=False)


