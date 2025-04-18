.. image:: _static/adcirc_utilities_log.png
   :alt: ADCIRC Utilities Logo
   :align: center

Overview
========

ADCIRC Utils provides a collection of tools for working with ADCIRC meshes and related data. These utilities are designed to simplify and enhance workflows for ADCIRC modelers.

Tools
-----

channelpaving
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tools for adding or editing 1D/2D channels in an ADCIRC mesh. Features include:

* Extracting channel widths and depths along stream centerlines.
* Generating a channel mesh and embedding them into a background mesh.

dem2adcdp
~~~~~~~~~~

Extract depths at nodes in an ADCIRC mesh from a provided Digital Elevation Model (DEM) file. Key functionalities:

* Compute mean, maximum, or minimum values in the DEM around a node.
* Specify minimum or maximum depths.
* Update nodes selectively using polygons.

hydrologyboundary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generates the flow boundary condition file (``fort.20``) for flow boundaries defined in a mesh file (``fort.14``). Features include:

* Extracting discharge data from USGS stations or NOAA National Water Model historical datasets.

mesh
~~~~

Tools for manipulating ADCIRC meshes. Features include:

* Merging multiple meshes with support for Vertical Element Wall (VEW)
* Subtracting one mesh from another while preserving boundaries
* Handling mesh boundaries and node renumbering
* Supporting different merging strategies (VEW boundary and merged nodes)

nodalattribute
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Creates nodal attribute values. Key functionality:

* Calculates Manning's n nodal attribute values from landuse data
* Map landuse categories to Manning's n values.
* Assign calculated values to nodes in an ADCIRC mesh.
* Support for various landuse data formats.

plot
~~~~

Provides plotting tools for ADCIRC simulation results. Key functionalities:

* Plot hydrographs from simulation results and observations. Observed water levels are downloaded from either NOAA or USGS data repository.
* Plot error histograms at stations.

postprocess
~~~~~~~~~~~~~~~~~~~~~~~~~~

Provides postprocessing tools for ADCIRC simulation results. Key functionality:

* Compute differences in maxele.63.nc files from different simulations.
* Compute disturbance values.

vewprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tools for processing Vertical Element Wall (VEW) in ADCIRC meshes. Features include:

* Converting VEW polylines to node strings in YAML format
* Adding VEWs to the mesh
* Scraping VEWs from the mesh
* Managing bank elevations and Manning's n values for VEW nodes 