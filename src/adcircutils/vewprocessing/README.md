# ADCIRC VEW Processing Example

This example demonstrates how to process VEW (Vegetated Edge of Water) boundaries in ADCIRC meshes. The workflow includes:

1. Converting polylines to node strings
2. Adding VEW boundaries to a mesh
3. Extracting VEW boundaries from an existing mesh

## Setup

First, let's import the necessary modules and set up the environment:

```python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from adcircutils import polyline_converter, vew_processor
from adcircutils.mesh import read_mesh, write_mesh
from adcircutils.plotting import plot_mesh, plot_polylines
```

## Step 1: Load and Plot Original Mesh

Let's load and visualize the original mesh without VEWs:

```python
# Load the original mesh
mesh = read_mesh('input/mesh_without_vews.14')

# Plot the mesh
plt.figure(figsize=(10, 10))
plot_mesh(mesh)
plt.title('Original Mesh without VEWs')
plt.show()
```

## Step 2: Read and Plot VEW Polylines

Next, we'll read and visualize the VEW polylines from a shapefile:

```python
# Read VEW polylines from shapefile
vew_polylines = read_polylines('input/vew_polylines.shp')

# Plot the polylines
plt.figure(figsize=(10, 10))
plot_polylines(vew_polylines)
plt.title('VEW Polylines')
plt.show()
```

## Step 3: Convert Polylines to Node Strings

Now, we'll convert the polylines to node strings:

```python
# Create output directory if it doesn't exist
os.makedirs('output', exist_ok=True)

# Convert polylines to node strings
polyline_converter.main(
    mesh_file='input/mesh_without_vews.14',
    shapefile='input/vew_polylines.shp',
    output_file='output/vewstrings.yaml',
    node_search_radius=100,  # meters
    elevation_threshold=0.5  # meters
)
```

## Step 4: Add VEW Boundaries to Mesh

After generating the node strings, we can add the VEW boundaries to the mesh:

```python
# Add VEW boundaries to mesh
vew_processor.add_vew_boundaries(
    mesh=mesh,
    vew_strings_file='output/vewstrings.yaml',
    output_file='output/mesh_with_vews.14'
)
```

## Step 5: Extract VEW Boundaries from Mesh

Finally, we can extract the VEW boundaries from the modified mesh:

```python
# Extract VEW boundaries
vew_processor.extract_vew_boundaries(
    mesh_file='output/mesh_with_vews.14',
    output_file='output/extracted_vewstrings.yaml'
)
```

## Notes

- The node search radius and elevation threshold parameters can be adjusted based on your specific requirements
- The output files will be saved in the `output` directory
- Make sure all input files are in the correct format and location before running the scripts