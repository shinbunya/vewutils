#!/bin/bash

echo "===================================================================="
echo "Generating mesh with VEWs from channel, land, and background meshes"
echo "===================================================================="


echo ""
echo "1. Combining channel and land meshes with VEWs"
echo "-----------------------------------------------------------------------"
python -m adcircutils.mesh.mesh_merger input/channel.14 input/land.14 -b vew -o output/ch_la.14 -d 'merged: channel + land with VEWs' 

echo ""
echo "2. Combining channel and land meshes by merging overlapping nodes"
echo "-----------------------------------------------------------------------"
python -m adcircutils.mesh.mesh_merger input/channel.14 input/land.14 -b merge -o output/ch_la_mgd.14 -d 'merged: channel + land with merged nodes' 

echo ""
echo ""
echo "3. Subtracting channel + land coverage from background"
echo "-----------------------------------------------------------------------"
python -m adcircutils.mesh.mesh_subtractor input/background.14 output/ch_la_mgd.14 -o output/background_subtracted.14 -d 'subtracted: background - (channel + land)'

echo ""
echo ""
echo "4. Merging channel, land, and subtracted background"
echo "-----------------------------------------------------------------------"
python -m adcircutils.mesh.mesh_merger output/ch_la.14 output/background_subtracted.14 -b merge -o output/ch_la_bg.14 -d 'merged: background_subtracted + channel + land'

echo ""
echo ""
echo "5. Adding land boundaries to the merged mesh"
echo "-----------------------------------------------------------------------"
python -m adcircutils.mesh.add_land_boundaries output/ch_la_bg.14 -o output/ch_la_bg_lb.14 -d 'merged: background_subtracted + channel + land with land boundaries'

echo ""
echo ""
echo "6. Ensuring elevations of VEW channel nodes to be lower than bank nodes"
echo "-----------------------------------------------------------------------"
python -m adcircutils.mesh.adjust_vew_channel_elevations output/ch_la_bg_lb.14 -o output/ch_la_bg_lb_adjusted1.14

echo ""
echo ""
echo "7. Adjusting VEW barrier heights to be above the bank nodes"
echo "-----------------------------------------------------------------------"
python -m adcircutils.mesh.adjust_vew_barrier_heights output/ch_la_bg_lb_adjusted1.14 -o output/ch_la_bg_lb_adjusted2.14

echo ""
echo ""
echo "8. Copying nodal attributes in the background mesh to the new mesh."
echo "-----------------------------------------------------------------------"
python -m adcircutils.nodalattribute.attribute_transfer input/background.14 input/background.13 output/ch_la_bg_lb_adjusted2.14 -o output/ch_la_bg_lb_adjusted2.13

echo ""
echo ""
echo "9. Updating Manning's n values in the new mesh."
echo "-----------------------------------------------------------------------"
echo "9.1 Selecting channel mesh nodes from the new mesh."
echo "-----------------------------------------------------------------------"
python -m adcircutils.utils.node_selector output/ch_la_bg_lb_adjusted2.14 -o output/ch_la_bg_lb_adjusted2_channel_mesh_nodes.csv -m output/ch_la_mgd.14

echo ""
echo "9.2 Updating Manning's n values at the selected nodes."
echo "-----------------------------------------------------------------------"
python -m  adcircutils.nodalattribute.manningsn_extractor output/ch_la_bg_lb_adjusted2.14 input/ccap_landuse_sample.tif -s output/ch_la_bg_lb_adjusted2_channel_mesh_nodes.csv input/ccap_class_to_mn_openwater0.02.csv -f output/ch_la_bg_lb_adjusted2.13 -o output/ch_la_bg_lb_adjusted2_mn_updated.13 --format fort13

echo ""
echo "All done!"
