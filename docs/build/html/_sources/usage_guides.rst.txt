Usage Guides
============

This section provides detailed usage guides for all command-line tools in ADCIRC Utils.

Channelpaving
----

The channelpaving module provides tools for adding or editing 1D/2D channels in an ADCIRC mesh. The primary interface is through Python API and MATLAB scripts rather than command-line tools. See the examples directory for detailed examples of using this module.

DEM2ADCDP
----

Maps Digital Elevation Model (DEM) data onto ADCIRC mesh nodes.

.. code-block:: bash

    python -m adcircutils.dem2adcdp.dem2adcdp meshfile outmeshfile tiffile [options]

Arguments:

* ``meshfile``: Input mesh file
* ``outmeshfile``: Output mesh file with updated depths
* ``tiffile``: DEM file in GeoTIFF format

Options:

* ``--f13file FILE``: F13 file containing VEW information
* ``--cachefile FILE``: Cache file for storing intermediate results
* ``--assign_channelnode_depths``: Force to assign depths at channel nodes
* ``--target_add_all``: Add all nodes to target list
* ``--target_remove_all``: Remove all nodes from target list
* ``--target_add_by_polygons POLYGONS``: Add nodes within specified polygons (comma-separated list of polygon files)
* ``--target_remove_by_polygons POLYGONS``: Remove nodes within specified polygons
* ``--target_remove_by_polygons_outside POLYGONS``: Remove nodes outside specified polygons
* ``--target_add_channelnodes``: Add channel nodes to target list
* ``--target_remove_channelnodes``: Remove channel nodes from target list
* ``--target_add_banknodes``: Add bank nodes to target list
* ``--target_remove_banknodes``: Remove bank nodes from target list
* ``--target_add_boundarynodes``: Add boundary nodes to target list
* ``--channel_deeper_by FLOAT``: Difference in depth by which channel node is deepened compared to bank (default: 0.001)
* ``--channel_deeper_by_threshold FLOAT``: Threshold for applying channel_deeper_by option (default: 10000)
* ``--deepen``: Assign an extracted value only when it is deeper than the current depth
* ``--land_only``: Assign depth only when the resulting elevation is above 0 m
* ``--submerged_only``: Assign depth only when both the current and resulting elevations are below 0 m
* ``--min_depth FLOAT``: Minimum depth value to be assigned to mesh nodes (default: -100000.0)
* ``--min_depth_tapering_end FLOAT``: Depth at which tapering of the minimum depth application ends
* ``--max_depth FLOAT``: Maximum depth value to be assigned to mesh nodes (default: 100000.0)
* ``--method {mean,max,min}``: Extract method for elevation values (default: mean)
* ``--ignore_tiff``: Ignore values in tiff file and apply min_depth/max_depth only
* ``--ncores INT``: Number of CPU cores to use (default: 1)
* ``--chunk_size_poly INT``: Chunk size for polygon operations (default: 1000)
* ``--chunk_size_zonalstats INT``: Chunk size for zonal statistics operations (default: 1000)

Hydrologyboundary
----

The hydrologyboundary module provides tools for generating flow boundary condition files for ADCIRC. It supports retrieving flow data from both USGS gauges and the National Water Model (NWM).

Primary components:
* ``fluxboundaries``: Core module for handling flux boundary conditions
* ``nwm``: Utilities for working with National Water Model data
* ``usgs``: Utilities for working with USGS streamflow data

Example usage in Python:

.. code-block:: python

    from adcircutils.hydrologyboundary import fluxboundaries as fluxb
    from adcircutils.hydrologyboundary import usgs
    from adcircutils.hydrologyboundary import nwm
    from datetime import datetime, timezone, timedelta
    
    # Create a flux boundary handler
    flux_boundaries = fluxb.FluxBoundaries()
    
    # Generate flux boundary segments from an ADCIRC mesh
    adcmesh_file = "path/to/mesh.grd"
    flux_boundaries.generate_fluxboundary_segments_from_adcmesh(adcmesh_file)
    
    # Initialize boundary definitions
    flux_boundaries.init_boundary_defs()
    
    # Set boundary conditions using different methods:
    
    # 1. Set USGS gauge at specific location
    flux_boundaries.set_boundary_defs_at_lonlat(-77.5825, 35.4289, "usgs_siteid", '02091500')
    
    # 2. Set NWM feature at specific location
    flux_boundaries.set_boundary_defs_at_lonlat(-77.9974, 35.3370, "nwm_featureid", 11239411)
    
    # 3. Set constant flow rate (in m³/s)
    flux_boundaries.set_boundary_defs_at_lonlat(-77.5078, 35.4230, "const_m3s", 0.0)
    
    # 4. Automatically assign NWM features to remaining boundaries
    hydrofabricfile = "path/to/nwm_channels.shp"
    flux_boundaries.set_boundary_defs_to_nwm_at_NA(hydrofabricfile, max_dist=100.0)
    
    # 5. Set remaining undefined boundaries to zero flow
    flux_boundaries.set_boundary_defs_at_NA("const_m3s", 0.0)
    
    # Process the boundary definitions
    flux_boundaries.digest_boundary_defs()
    
    # Set the time range and generate discharges
    start = datetime(2018, 8, 23, 0, 0, tzinfo=timezone.utc)
    end = datetime(2018, 10, 9, 3, 0, tzinfo=timezone.utc)
    step = timedelta(seconds=900)  # 15-minute intervals
    
    # Retrieve discharge data from database
    flux_boundaries.set_discharge_from_database(start, end)
    
    # Write the fort.20 file
    flux_boundaries.write_fort20("output.fort20", start, end, step)

Mesh
----

Mesh Merger
~~~~~~~~~~

Merges multiple ADCIRC meshes together.

.. code-block:: bash

    python -m adcircutils.mesh.mesh_merger channelmesh landmesh [-o OUTPUT] [-d DESCRIPTION] [-b {merge,vew}] [--use-land-values]

Arguments:

* ``channelmesh``: Path to the channel mesh file
* ``landmesh``: Path to the land mesh file
* ``-o, --output``: Path to save the merged mesh (default: merged_mesh.grd)
* ``-d, --description``: Description for the merged mesh (default: merged)
* ``-b, --boundary-mode``: How to handle mesh boundaries: 'merge' to merge nodes, 'vew' to use VEW boundaries (default: merge)
* ``--use-land-values``: Use land mesh values at matching nodes (default: use channel mesh values)

Mesh Subtractor
~~~~~~~~~~~~~~

Subtracts one mesh from another while preserving boundaries.

.. code-block:: bash

    python -m adcircutils.mesh.mesh_subtractor mesh_a mesh_b [-o OUTPUT] [-d DESCRIPTION]

Arguments:

* ``mesh_a``: Path to first mesh file (mesh to subtract from)
* ``mesh_b``: Path to second mesh file (mesh defining subtraction boundary)
* ``-o, --output``: Output mesh file path (default: subtracted_mesh.grd)
* ``-d, --description``: Description for the subtracted mesh (default: subtracted)

Add Land Boundaries
~~~~~~~~~~~~~~~~~

Adds land boundaries to an ADCIRC mesh.

.. code-block:: bash

    python -m adcircutils.mesh.add_land_boundaries input_mesh [-o OUTPUT] [-d DESCRIPTION]

Arguments:

* ``input_mesh``: Path to the input mesh file
* ``-o, --output``: Path to save the output mesh (default: mesh_with_land_boundaries.grd)
* ``-d, --description``: Description for the output mesh (default: mesh with land boundaries)

Adjust VEW Channel Elevations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adjusts the elevation of VEW channels in an ADCIRC mesh.

.. code-block:: bash

    python -m adcircutils.mesh.adjust_vew_channel_elevations input_mesh [-o OUTPUT] [-t TOLERANCE]

Arguments:

* ``input_mesh``: Path to the input mesh file
* ``-o, --output``: Path to save the output mesh file (default: adjusted_mesh.14)
* ``-t, --tolerance``: Amount to lower channel node elevations below bank node elevations in meters (default: 0.001)

Adjust VEW Barrier Heights
~~~~~~~~~~~~~~~~~~~~~~~~

Adjusts the height of VEW barriers in an ADCIRC mesh.

.. code-block:: bash

    python -m adcircutils.mesh.adjust_vew_barrier_heights input_mesh [-o OUTPUT] [-t TOLERANCE]

Arguments:

* ``input_mesh``: Path to the input mesh file
* ``-o, --output``: Path to save the output mesh file (default: adjusted_mesh.14)
* ``-t, --tolerance``: Minimum amount that barrier heights should be above bank elevations in meters (default: 0.001)

Nodal Attribute
--------------

Manning's n Extractor
~~~~~~~~~~~~~~~~~~~

Extracts Manning's n values from landuse data.

.. code-block:: bash

    python -m adcircutils.nodalattribute.manningsn_extractor mesh landcover class_to_mn [-o OUTPUT] [-s SELECTED-NODES] [-f FORT13] [--format {fort13,csv}]

Arguments:

* ``mesh``: Path to the mesh file
* ``landcover``: Path to the land cover raster file
* ``class_to_mn``: Path to the class to Manning's n mapping file
* ``-o, --output``: Output file path (default: mannings_n.fort.13)
* ``-s, --selected-nodes``: Path to CSV file containing selected nodes
* ``-f, --fort13``: Path to input fort.13 file
* ``--format``: Output format: 'fort13' or 'csv' (default: fort13)

Attribute Transfer
~~~~~~~~~~~~~~~~

Transfers attributes from one mesh to another.

.. code-block:: bash

    python -m adcircutils.nodalattribute.attribute_transfer source_mesh source_attrs target_mesh [-o OUTPUT]

Arguments:

* ``source_mesh``: Path to source grid file (fort.14)
* ``source_attrs``: Path to source nodal attribute file (fort.13)
* ``target_mesh``: Path to target grid file (fort.14)
* ``-o, --output``: Output nodal attribute file path (default: transferred_attributes.fort.13)

Plot
----

Plot Hydrograph at Station
~~~~~~~~~~~~~~~~~~~~~~~~

Plots hydrographs from simulation results and observations.

.. code-block:: bash

    python -m adcircutils.plot.plot_hydrograph_at_station --station_owner STATION_OWNER --station_id STATION_ID [--station_lon STATION_LON] [--station_lat STATION_LAT] --station_datum STATION_DATUM --date_start DATE_START --date_end DATE_END --f63files F63FILES [F63FILES ...] --f63starts F63STARTS [F63STARTS ...] --f63labels F63LABELS [F63LABELS ...] [--plot_movingaverage] --outputfile OUTPUTFILE

Arguments:

* ``--station_owner``: Station owner: NOAA or USGS
* ``--station_id``: Station ID
* ``--station_lon``: Station longitude (optional)
* ``--station_lat``: Station latitude (optional)
* ``--station_datum``: Station datum: MSL or NAVD
* ``--date_start``: Start date (YYYY-MM-DD)
* ``--date_end``: End date (YYYY-MM-DD)
* ``--f63files``: List of fort.63.nc files
* ``--f63starts``: List of f63 start times (YYYY-MM-DD)
* ``--f63labels``: List of labels for f63 files
* ``--plot_movingaverage``: Plot moving average (flag)
* ``--outputfile``: Output figure file name

Utils
-----

Node Selector
~~~~~~~~~~~~

Selects nodes from an ADCIRC mesh based on specified criteria.

.. code-block:: bash

    python -m adcircutils.utils.node_selector mesh_file [-o OUTPUT] [-p POLYGON] [-m BOUNDARY-MESH] [-c CSV] [-t TOLERANCE] [-b {channel,bank,both}] [-op {union,intersection}]

Arguments:

* ``mesh_file``: Path to the ADCIRC mesh file
* ``-o, --output``: Output CSV file path (default: selected_nodes.csv)
* ``-p, --polygon``: Path to geospatial file containing polygons
* ``-m, --boundary-mesh``: Path to mesh file defining boundary
* ``-c, --csv``: Path to CSV file containing node IDs
* ``-t, --tolerance``: Tolerance in meters for boundary proximity (default: 1e-6)
* ``-b, --boundary-type``: Type of boundary nodes to include: 'channel', 'bank', or 'both' (default: both)
* ``-op, --operation``: Operation to combine multiple selections: 'union' or 'intersection' (default: union)

VEW Processing
--------------

Polyline Converter
~~~~~~~~~~~~~~~~~

Converts polylines to VEW strings in YAML format.

.. code-block:: bash

    python -m adcircutils.vewprocessing.polyline_converter meshfile polylinefile [-o OUTPUT] [-d DISTANCE] [-e ELEVATION] [-n MANNINGS]

Arguments:

* ``meshfile``: Path to the ADCIRC mesh file (fort.14)
* ``polylinefile``: Path to the polyline file (shapefile, geojson, etc.)
* ``-o, --output``: Path to save the VEW strings (default: vewstrings.yaml)
* ``-d, --distance``: Maximum distance for nearest neighbor search in meters (default: 10.0)
* ``-e, --elevation``: Bank elevation for VEW strings (default: 1.0)
* ``-n, --mannings``: Manning's n value for VEW strings (default: 0.02)

VEW Adder
~~~~~~~~~

Adds Vertical Element Walls to a mesh using node strings from a YAML file.

.. code-block:: bash

    python -m adcircutils.vewprocessing.vew_adder f14file vewfile -o OUTPUT

Arguments:

* ``f14file``: Input fort.14 file
* ``vewfile``: Input YAML file containing VEW string definitions
* ``-o, --output``: Output fort.14 file

VEW Scraper
~~~~~~~~~~~

Extracts VEW boundaries from an ADCIRC mesh and saves them to YAML format.

.. code-block:: bash

    python -m adcircutils.vewprocessing.vew_scraper input_mesh -o OUTPUT_MESH -y OUTPUT_YAML [-d DESCRIPTION]

Arguments:

* ``input_mesh``: Path to the input mesh file with VEW boundaries
* ``-o, --output-mesh``: Path to save the mesh without VEW boundaries
* ``-y, --output-yaml``: Path to save the extracted VEW strings in YAML format
* ``-d, --description``: Description for the output mesh (default: Generated by vew_scraper) 