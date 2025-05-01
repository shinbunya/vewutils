#!/bin/bash

# Create output directory if it doesn't exist
mkdir -p output

echo "===================================================================="
echo "Processing VEWs in ADCIRC mesh"
echo "===================================================================="

echo ""
echo "1. Converting VEW polylines to node strings of a mesh in YAML format"
echo "-----------------------------------------------------------------------"
python -m vewutils.vewprocessing.polyline_converter input/mesh_without_vews.14 input/vew_polylines.shp -o output/vewstrings.yaml -n 0.03 -e 2.0

echo ""
echo "2. Adding VEWs to the mesh"
echo "-----------------------------------------------------------------------"
python -m vewutils.vewprocessing.vew_adder input/mesh_without_vews.14 output/vewstrings.yaml -o output/mesh_with_vews.14

echo ""
echo "3. Scraping VEWs from the mesh"
echo "-----------------------------------------------------------------------"
python -m vewutils.vewprocessing.vew_scraper output/mesh_with_vews.14 -o output/mesh_vews_scraped.14 -y output/vewstrings_scraped.yaml

echo ""
echo "All done!"