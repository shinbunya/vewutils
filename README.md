![Header](image/adcirc_utilities_log.png)

# ADCIRC Utilities

This repository provides a collection of tools for working with ADCIRC meshes and related data. These utilities are designed to simplify and enhance workflows for ADCIRC modelers.

## Tools

### 1. `channelpaving`
Tools for adding or editing 1D/2D channels in an ADCIRC mesh. Features include:
- Extracting channel widths and depths along stream centerlines.
- Generating a channel mesh and embedding them into a background mesh.

### 2. `dem2adcdp`
Extract depths at nodes in an ADCIRC mesh from a provided Digital Elevation Model (DEM) file. Key functionalities:
- Compute mean, maximum, or minimum values in the DEM around a node.
- Specify minimum or maximum depths.
- Update nodes selectively using polygons.

### 3. `hydrologyboundary`
Generates the flow boundary condition file (`fort.20`) for flow boundaries defined in a mesh file (`fort.14`). Features include:
- Extracting discharge data from USGS stations or NOAA National Water Model historical datasets.

### 4. `mesh`
Tools for manipulating ADCIRC meshes. Features include:
- Merging multiple meshes with support for Vertical Element Wall (VEW)
- Subtracting one mesh from another while preserving boundaries
- Handling mesh boundaries and node renumbering
- Supporting different merging strategies (VEW boundary and merged nodes)

### 5. `nodalattribute`
Creates nodal attribute values. Key functionality:
- Calculates Manning's n nodal attribute values from landuse data
  - Map landuse categories to Manning's n values.
  - Assign calculated values to nodes in an ADCIRC mesh.
  - Support for various landuse data formats.

### 6. `plot`
Provides plotting tools for ADCIRC simulation results. Key functionalities:
- Plot hydrographs from simulation results and observations. Observed water levels are downloaded from either NOAA or USGS data repository.
- Plot error histograms at stations.

### 7. `postprocess`
Provides postprocessing tools for ADCIRC simulation results. Key functionality:
- Compute differences in maxele.63.nc files from different simulations.
- Compute disturbance values.

### 8. `vewprocessing`
Tools for processing Vertical Element Wall (VEW) in ADCIRC meshes. Features include:
- Converting VEW polylines to node strings in YAML format
- Adding VEWs to the mesh
- Scraping VEWs from the mesh
- Managing bank elevations and Manning's n values for VEW nodes

## Getting Started

To use these tools, clone the repository and install it by executing `install.sh`. Install `conda` before using the installation script. The installation script creates a new conda environment, the name of which is specified as a command line arugment. It install all dependencies, including a fork of adcircpy. It then installs `adcircutils` using `pip` command. Run the following commands for installation.

```bash
git clone https://github.com/shinbunya/adcircutils.git adcircutils
cd adcircutils
./install.sh adcircutils # <-- You can specify your conda environment name.
```

## Links to Examples

### [Channel Mesh Merging (examples/channelmerging/)](examples/channelmerging/example.ipynb)

This example demonstrates the process of merging multiple meshes (channel, land, and background) while handling Vertical Element Wall (VEW). The workflow includes:
- Combining channel and land meshes with VEWs
- Subtracting channel+land coverage from background mesh
- Merging all components into a final mesh
- Adjusting VEW channel elevations and barrier heights
- Transferring nodal attributes and updating Manning's n values

![Channel Merging Example](examples/channelmerging/image/image01.png)

### [Channel Mesh Paving (examples/channelpaving/)](examples/channelpaving/README.md)

This example demonstrates the channel paving process using tools in `adcircutils/channelpaving`. It provides a sample setup to showcase how to model and pave a channel in a background mesh for ADCIRC simulations. The example includes:
- Adding depth and width attributes to channel centerlines
- Creating a channel mesh and embedding it into a background mesh
- Visualizing the results using MATLAB Live Scripts

![Channel Paving Example](examples/channelpaving/image/image01.png)

### [VEW Processing (examples/vewprocessing/)](examples/vewprocessing/example.ipynb)

This example demonstrates the complete workflow for processing Vertical Element Wall (VEW) in ADCIRC meshes using the `adcircutils.vewprocessing` module. The process includes:
- Converting VEW polylines to node strings
- Adding VEWs to the mesh
- Scraping VEWs from the mesh
- Managing bank elevations and Manning's n values for VEW nodes (Application of given Manning's n values has not been implemented yet.)

![VEW Processing Example](examples/vewprocessing/image/image01.png)

## References

1. Bunya, S., et al. (2023). Techniques to embed channels in finite element shallow water equation models. Advances in Engineering Software, 103516. https://doi.org/10.1016/j.advengsoft.2023.103516

## License

This project is licensed under the [Apatche 2.0 License](https://opensource.org/license/apache-2-0).

## Contact

For questions or support, please open an issue in the repository or contact sbunya@unc.edu.
