#!/usr/bin/env python3
import argparse
import sys

# Import get_parser and main for each command
from vewutils.mesh.mesh_merger import get_parser as mesh_merge_get_parser, main as mesh_merge_main
from vewutils.mesh.mesh_subtractor import get_parser as mesh_subtract_get_parser, main as mesh_subtract_main
from vewutils.mesh.boundary_exporter import get_parser as mesh_export_boundaries_get_parser, main as mesh_export_boundaries_main
from vewutils.mesh.export_submesh import get_parser as mesh_export_submesh_get_parser, main as mesh_export_submesh_main
from vewutils.mesh.add_land_boundaries import get_parser as mesh_add_land_boundaries_get_parser, main as mesh_add_land_boundaries_main
from vewutils.mesh.adjust_vew_barrier_heights import get_parser as mesh_adjust_vew_barrier_heights_get_parser, main as mesh_adjust_vew_barrier_heights_main
from vewutils.mesh.adjust_vew_channel_elevations import get_parser as mesh_adjust_vew_channel_elevations_get_parser, main as mesh_adjust_vew_channel_elevations_main
from vewutils.mesh.adjust_vew_node_coords import get_parser as mesh_adjust_vew_node_coords_get_parser, main as mesh_adjust_vew_node_coords_main
from vewutils.mesh.vew_boundary_repairer import get_parser as mesh_repair_vew_boundaries_get_parser, main as mesh_repair_vew_boundaries_main
from vewutils.mesh.flow_boundary_repairer import get_parser as mesh_repair_flow_boundaries_get_parser, main as mesh_repair_flow_boundaries_main
from vewutils.mesh.bandwidth_reducer import get_parser as mesh_reduce_bandwidth_get_parser, main as mesh_reduce_bandwidth_main
from vewutils.mesh.copy_depths import get_parser as mesh_copy_depths_get_parser, main as mesh_copy_depths_main
from vewutils.plot.plot_hydrograph_at_station import get_parser as plot_hydrograph_get_parser, main as plot_hydrograph_main
from vewutils.plot.plot_errorhistogram_at_station import get_parser as plot_errorhistogram_get_parser, main as plot_errorhistogram_main
from vewutils.plot.get_obswl import get_parser as plot_get_obswl_get_parser, main as plot_get_obswl_main
from vewutils.plot.plot_solution_at import get_parser as plot_solution_at_get_parser, main as plot_solution_at_main
from vewutils.plot.plot_solution_2d import get_parser as plot_solution_2d_get_parser, main as plot_solution_2d_main
from vewutils.plot.plot_max_ele_2d import get_parser as plot_max_ele_2d_get_parser, main as plot_max_ele_2d_main
from vewutils.plot.plot_max_one_to_one_at_stations import get_parser as plot_max_one_to_one_get_parser, main as plot_max_one_to_one_main
from vewutils.post.maxele_max import get_parser as post_maxele_max_get_parser, main as post_maxele_max_main
from vewutils.post.maxele_diff import get_parser as post_maxele_diff_get_parser, main as post_maxele_diff_main
from vewutils.post.maxele_add_disturbance import get_parser as post_maxele_add_disturbance_get_parser, main as post_maxele_add_disturbance_main
from vewutils.post.maxele_add_departure import get_parser as post_maxele_add_departure_get_parser, main as post_maxele_add_departure_main
from vewutils.post.maxele_attribution import get_parser as post_maxele_attribution_get_parser, main as post_maxele_attribution_main
from vewutils.post.create_da_error_surface_nc import get_parser as post_create_da_error_surface_get_parser, main as post_create_da_error_surface_main
from vewutils.post.convert_maxelelike_to_fort63like import get_parser as post_convert_maxelelike_to_fort63like_get_parser, main as post_convert_maxelelike_to_fort63like_main
from vewutils.post.images_to_mp4 import get_parser as post_images_to_mp4_get_parser, main as post_images_to_mp4_main
from vewutils.post.reduce_timesteps import get_parser as post_reduce_timesteps_get_parser, main as post_reduce_timesteps_main
from vewutils.post.rechunk_netcdf import get_parser as post_rechunk_netcdf_get_parser, main as post_rechunk_netcdf_main
from vewutils.post.concat_fort61 import get_parser as post_concat_fort61_get_parser, main as post_concat_fort61_main
from vewutils.post.replace_f61_with_f63_at_station import get_parser as post_replace_f61_get_parser, main as post_replace_f61_main
from vewutils.post.fort63_add_departure import get_parser as post_fort63_add_departure_get_parser, main as post_fort63_add_departure_main
from vewutils.dem2adcdp.dem2adcdp import get_parser as dem2adcdp_get_parser, main as dem2adcdp_main
from vewutils.nodalattribute.attribute_transfer import get_parser as nodalattribute_transfer_get_parser, main as nodalattribute_transfer_main
from vewutils.nodalattribute.manningsn_extractor import get_parser as nodalattribute_extract_manningsn_get_parser, main as nodalattribute_extract_manningsn_main
from vewutils.nodalattribute.manningsn_setter import get_parser as nodalattribute_set_manningsn_get_parser, main as nodalattribute_set_manningsn_main
from vewutils.nodalattribute.vew_manningsn_setter import get_parser as nodalattribute_set_vew_manningsn_get_parser, main as nodalattribute_set_vew_manningsn_main
from vewutils.vewprocessing.vew_adder import get_parser as vewprocessing_add_get_parser, main as vewprocessing_add_main
from vewutils.vewprocessing.polyline_converter import get_parser as vewprocessing_convert_polylines_get_parser, main as vewprocessing_convert_polylines_main
from vewutils.vewprocessing.vew_scraper import get_parser as vewprocessing_scrape_get_parser, main as vewprocessing_scrape_main
from vewutils.vewprocessing.connect_vewstrings import get_parser as vewprocessing_connect_vewstrings_get_parser, main as vewprocessing_connect_vewstrings_main
from vewutils.vewprocessing.yaml2geojson import get_parser as vewprocessing_yaml2geojson_get_parser, main as vewprocessing_yaml2geojson_main
from vewutils.utils.node_selector import get_parser as utils_select_nodes_get_parser, main as utils_select_nodes_main
from vewutils.channelpaving.NCFRISCrossSect2Depth import get_parser as channelpaving_add_depth_get_parser, main as channelpaving_add_depth_main
from vewutils.channelpaving.NHDArea2Width import get_parser as channelpaving_add_width_get_parser, main as channelpaving_add_width_main
from vewutils.channelpaving.channel_mesh_generator import get_parser as channelpaving_generate_mesh_get_parser, main as channelpaving_generate_mesh_main

def main():
    parser = argparse.ArgumentParser(
        description='vewutils: Unified CLI for VEW utility programs'
    )
    subparsers = parser.add_subparsers(dest='group', help='Command group')

    # Mesh group
    mesh_parser = subparsers.add_parser('mesh', help='Mesh manipulation commands')
    mesh_subparsers = mesh_parser.add_subparsers(dest='mesh_cmd')
    mesh_merge_parser = mesh_subparsers.add_parser(
        'merge',
        help='Merge ADCIRC meshes',
        parents=[mesh_merge_get_parser()],
        add_help=True
    )
    mesh_merge_parser.set_defaults(func=mesh_merge_main)
    mesh_subtract_parser = mesh_subparsers.add_parser(
        'subtract',
        help='Subtract meshes',
        parents=[mesh_subtract_get_parser()],
        add_help=True
    )
    mesh_subtract_parser.set_defaults(func=mesh_subtract_main)
    mesh_export_parser = mesh_subparsers.add_parser(
        'export-boundaries',
        help='Export mesh boundaries to GeoPackage',
        parents=[mesh_export_boundaries_get_parser()],
        add_help=True
    )
    mesh_export_parser.set_defaults(func=mesh_export_boundaries_main)
    mesh_export_submesh_parser = mesh_subparsers.add_parser(
        'export-submesh',
        help='Export a detached submesh containing a specified element',
        parents=[mesh_export_submesh_get_parser()],
        add_help=True
    )
    mesh_export_submesh_parser.set_defaults(func=mesh_export_submesh_main)
    mesh_add_land_parser = mesh_subparsers.add_parser(
        'add-land-boundaries',
        help='Add land boundaries to mesh',
        parents=[mesh_add_land_boundaries_get_parser()],
        add_help=True
    )
    mesh_add_land_parser.set_defaults(func=mesh_add_land_boundaries_main)
    mesh_adjust_barrier_parser = mesh_subparsers.add_parser(
        'adjust-vew-barrier-heights',
        help='Ensure VEW barrier heights are above bank elevations',
        parents=[mesh_adjust_vew_barrier_heights_get_parser()],
        add_help=True
    )
    mesh_adjust_barrier_parser.set_defaults(func=mesh_adjust_vew_barrier_heights_main)
    mesh_adjust_channel_parser = mesh_subparsers.add_parser(
        'adjust-vew-channel-elevations',
        help='Lower channel node elevations below bank nodes in VEW boundaries',
        parents=[mesh_adjust_vew_channel_elevations_get_parser()],
        add_help=True
    )
    mesh_adjust_channel_parser.set_defaults(func=mesh_adjust_vew_channel_elevations_main)
    mesh_adjust_coords_parser = mesh_subparsers.add_parser(
        'adjust-vew-node-coords',
        help='Match coordinates of VEW boundary node pairs',
        parents=[mesh_adjust_vew_node_coords_get_parser()],
        add_help=True
    )
    mesh_adjust_coords_parser.set_defaults(func=mesh_adjust_vew_node_coords_main)
    mesh_repair_vew_parser = mesh_subparsers.add_parser(
        'repair-vew-boundaries',
        help='Repair corrupted VEW boundaries where channel and bank nodes are identical',
        parents=[mesh_repair_vew_boundaries_get_parser()],
        add_help=True
    )
    mesh_repair_vew_parser.set_defaults(func=mesh_repair_vew_boundaries_main)
    mesh_repair_flow_parser = mesh_subparsers.add_parser(
        'repair-flow-boundaries',
        help='Repair corrupted flow boundaries where bank nodes are used instead of channel nodes',
        parents=[mesh_repair_flow_boundaries_get_parser()],
        add_help=True
    )
    mesh_repair_flow_parser.set_defaults(func=mesh_repair_flow_boundaries_main)
    mesh_reduce_bandwidth_parser = mesh_subparsers.add_parser(
        'reduce-bandwidth',
        help='Reduce mesh bandwidth using node renumbering algorithms',
        parents=[mesh_reduce_bandwidth_get_parser()],
        add_help=True
    )
    mesh_reduce_bandwidth_parser.set_defaults(func=mesh_reduce_bandwidth_main)
    mesh_copy_depths_parser = mesh_subparsers.add_parser(
        'copy-depths',
        help='Copy depths at selected nodes from source mesh to base mesh',
        parents=[mesh_copy_depths_get_parser()],
        add_help=True
    )
    mesh_copy_depths_parser.set_defaults(func=mesh_copy_depths_main)

    # Plot group
    plot_parser = subparsers.add_parser('plot', help='Plotting commands')
    plot_subparsers = plot_parser.add_subparsers(dest='plot_cmd')
    plot_hydro_parser = plot_subparsers.add_parser(
        'hydrograph',
        help='Plot hydrograph at a station',
        parents=[plot_hydrograph_get_parser()],
        add_help=True
    )
    plot_hydro_parser.set_defaults(func=plot_hydrograph_main)
    plot_errorhist_parser = plot_subparsers.add_parser(
        'errorhistogram',
        help='Plot error histogram at a station',
        parents=[plot_errorhistogram_get_parser()],
        add_help=True
    )
    plot_errorhist_parser.set_defaults(func=plot_errorhistogram_main)
    plot_get_obswl_parser = plot_subparsers.add_parser(
        'get-obswl',
        help='Retrieve observed water level data from various sources',
        parents=[plot_get_obswl_get_parser()],
        add_help=True
    )
    plot_get_obswl_parser.set_defaults(func=plot_get_obswl_main)
    plot_solution_at_parser = plot_subparsers.add_parser(
        'solution-at',
        help='Plot solution variables at a specific point',
        parents=[plot_solution_at_get_parser()],
        add_help=True
    )
    plot_solution_at_parser.set_defaults(func=plot_solution_at_main)
    plot_solution_2d_parser = plot_subparsers.add_parser(
        'solution-2d',
        help='Plot CG ADCIRC water levels and velocity fields from NetCDF files',
        parents=[plot_solution_2d_get_parser(add_help=False)],
        add_help=True
    )
    plot_solution_2d_parser.set_defaults(func=plot_solution_2d_main)
    plot_max_ele_2d_parser = plot_subparsers.add_parser(
        'maxele-2d',
        help='Plot CG ADCIRC maximum water level fields from maxele NetCDF files',
        parents=[plot_max_ele_2d_get_parser(add_help=False)],
        add_help=True
    )
    plot_max_ele_2d_parser.set_defaults(func=plot_max_ele_2d_main)
    plot_max_one_to_one_parser = plot_subparsers.add_parser(
        'maxele-scatter',
        help='Plot one-to-one maximum values at multiple stations',
        parents=[plot_max_one_to_one_get_parser()],
        add_help=True
    )
    plot_max_one_to_one_parser.set_defaults(func=plot_max_one_to_one_main)

    # Post group
    post_parser = subparsers.add_parser('post', help='Postprocessing commands')
    post_subparsers = post_parser.add_subparsers(dest='post_cmd')
    post_maxele_max_parser = post_subparsers.add_parser(
        'maxele-max',
        help='Calculate max across multiple maxele files',
        parents=[post_maxele_max_get_parser()],
        add_help=True
    )
    post_maxele_max_parser.set_defaults(func=post_maxele_max_main)
    post_maxele_diff_parser = post_subparsers.add_parser(
        'maxele-diff',
        help='Calculate diff between two maxele files',
        parents=[post_maxele_diff_get_parser()],
        add_help=True
    )
    post_maxele_diff_parser.set_defaults(func=post_maxele_diff_main)
    post_maxele_add_dist_parser = post_subparsers.add_parser(
        'maxele-add-disturbance',
        help='Add disturbance to maxele file',
        parents=[post_maxele_add_disturbance_get_parser()],
        add_help=True
    )
    post_maxele_add_dist_parser.set_defaults(func=post_maxele_add_disturbance_main)
    post_maxele_add_departure_parser = post_subparsers.add_parser(
        'maxele-add-departure',
        help='Add departure field to maxele file',
        parents=[post_maxele_add_departure_get_parser()],
        add_help=True
    )
    post_maxele_add_departure_parser.set_defaults(func=post_maxele_add_departure_main)
    post_maxele_attr_parser = post_subparsers.add_parser(
        'maxele-attribution',
        help='Add attribution to maxele file',
        parents=[post_maxele_attribution_get_parser()],
        add_help=True
    )
    post_maxele_attr_parser.set_defaults(func=post_maxele_attribution_main)
    post_create_da_error_surface_parser = post_subparsers.add_parser(
        'create-da-error-surface',
        help='Convert data assimilation error surface to maxele.63.nc format',
        parents=[post_create_da_error_surface_get_parser()],
        add_help=True
    )
    post_create_da_error_surface_parser.set_defaults(func=post_create_da_error_surface_main)
    post_convert_maxelelike_parser = post_subparsers.add_parser(
        'convert-maxelelike-to-fort63like',
        help='Convert maxele.63.nc-like files to fort.63.nc-like format',
        parents=[post_convert_maxelelike_to_fort63like_get_parser()],
        add_help=True
    )
    post_convert_maxelelike_parser.set_defaults(func=post_convert_maxelelike_to_fort63like_main)
    post_images_to_mp4_parser = post_subparsers.add_parser(
        'images-to-mp4',
        help='Convert image files to MP4 video using ffmpeg',
        parents=[post_images_to_mp4_get_parser()],
        add_help=True
    )
    post_images_to_mp4_parser.set_defaults(func=post_images_to_mp4_main)
    post_reduce_ts_parser = post_subparsers.add_parser(
        'reduce-timesteps',
        help='Reduce time steps in NetCDF',
        parents=[post_reduce_timesteps_get_parser()],
        add_help=True
    )
    post_reduce_ts_parser.set_defaults(func=post_reduce_timesteps_main)
    post_rechunk_netcdf_parser = post_subparsers.add_parser(
        'rechunk-netcdf',
        help='Rechunk NetCDF file to optimize for time-series access',
        parents=[post_rechunk_netcdf_get_parser()],
        add_help=True
    )
    post_rechunk_netcdf_parser.set_defaults(func=post_rechunk_netcdf_main)
    post_concat_fort61_parser = post_subparsers.add_parser(
        'concat-fort61',
        help='Concatenate time series data from multiple fort.61.nc files',
        parents=[post_concat_fort61_get_parser()],
        add_help=True
    )
    post_concat_fort61_parser.set_defaults(func=post_concat_fort61_main)
    post_replace_f61_parser = post_subparsers.add_parser(
        'replace-f61-station',
        help='Replace water level data for a station in fort.61.nc files with data from fort.63.nc',
        parents=[post_replace_f61_get_parser()],
        add_help=True
    )
    post_replace_f61_parser.set_defaults(func=post_replace_f61_main)
    post_fort63_add_departure_parser = post_subparsers.add_parser(
        'fort63-add-departure',
        help='Add departure field to fort.63.nc file',
        parents=[post_fort63_add_departure_get_parser()],
        add_help=True
    )
    post_fort63_add_departure_parser.set_defaults(func=post_fort63_add_departure_main)

    # DEM2ADCDP group (single command)
    dem2adcdp_parser = subparsers.add_parser(
        'dem2adcdp',
        help='Map DEM onto ADCIRC mesh nodes',
        parents=[dem2adcdp_get_parser()],
        add_help=True
    )
    dem2adcdp_parser.set_defaults(func=dem2adcdp_main)

    # Nodal Attribute group
    nodalattribute_parser = subparsers.add_parser('nodalattribute', help='Nodal attribute commands')
    nodalattribute_subparsers = nodalattribute_parser.add_subparsers(dest='nodalattribute_cmd')
    nodalattribute_transfer_parser = nodalattribute_subparsers.add_parser(
        'transfer',
        help='Transfer nodal attributes',
        parents=[nodalattribute_transfer_get_parser()],
        add_help=True
    )
    nodalattribute_transfer_parser.set_defaults(func=nodalattribute_transfer_main)
    nodalattribute_extract_mn_parser = nodalattribute_subparsers.add_parser(
        'extract-manningsn',
        help="Extract Manning's n values",
        parents=[nodalattribute_extract_manningsn_get_parser()],
        add_help=True
    )
    nodalattribute_extract_mn_parser.set_defaults(func=nodalattribute_extract_manningsn_main)
    nodalattribute_set_mn_parser = nodalattribute_subparsers.add_parser(
        'set-manningsn',
        help="Set Manning's n values at specific nodes",
        parents=[nodalattribute_set_manningsn_get_parser()],
        add_help=True
    )
    nodalattribute_set_mn_parser.set_defaults(func=nodalattribute_set_manningsn_main)
    nodalattribute_set_vew_mn_parser = nodalattribute_subparsers.add_parser(
        'set-vew-manningsn',
        help="Copy Manning's n values from channel nodes to bank nodes along VEW boundaries",
        parents=[nodalattribute_set_vew_manningsn_get_parser()],
        add_help=True
    )
    nodalattribute_set_vew_mn_parser.set_defaults(func=nodalattribute_set_vew_manningsn_main)

    # VEW Processing group
    vewprocessing_parser = subparsers.add_parser('vewprocessing', help='VEW processing commands')
    vewprocessing_subparsers = vewprocessing_parser.add_subparsers(dest='vewprocessing_cmd')
    vewprocessing_add_parser = vewprocessing_subparsers.add_parser(
        'add',
        help='Add VEW boundaries to mesh',
        parents=[vewprocessing_add_get_parser()],
        add_help=True
    )
    vewprocessing_add_parser.set_defaults(func=vewprocessing_add_main)
    vewprocessing_convert_parser = vewprocessing_subparsers.add_parser(
        'convert-polylines',
        help='Convert polylines to VEW strings',
        parents=[vewprocessing_convert_polylines_get_parser()],
        add_help=True
    )
    vewprocessing_convert_parser.set_defaults(func=vewprocessing_convert_polylines_main)
    vewprocessing_scrape_parser = vewprocessing_subparsers.add_parser(
        'scrape',
        help='Scrape VEW boundaries from mesh',
        parents=[vewprocessing_scrape_get_parser()],
        add_help=True
    )
    vewprocessing_scrape_parser.set_defaults(func=vewprocessing_scrape_main)
    vewprocessing_connect_parser = vewprocessing_subparsers.add_parser(
        'connect-vewstrings',
        help='Connect VEW strings from two YAML files',
        parents=[vewprocessing_connect_vewstrings_get_parser()],
        add_help=True
    )
    vewprocessing_connect_parser.set_defaults(func=vewprocessing_connect_vewstrings_main)
    vewprocessing_yaml2geojson_parser = vewprocessing_subparsers.add_parser(
        'yaml2geojson',
        help='Convert VEW string YAML files to GeoJSON format',
        parents=[vewprocessing_yaml2geojson_get_parser()],
        add_help=True
    )
    vewprocessing_yaml2geojson_parser.set_defaults(func=vewprocessing_yaml2geojson_main)

    # Utils group (optional, if implemented)
    utils_parser = subparsers.add_parser('utils', help='Utility commands')
    utils_subparsers = utils_parser.add_subparsers(dest='utils_cmd')
    utils_select_nodes_parser = utils_subparsers.add_parser(
        'select-nodes',
        help='Select nodes from mesh',
        parents=[utils_select_nodes_get_parser()],
        add_help=True
    )
    utils_select_nodes_parser.set_defaults(func=utils_select_nodes_main)

    # Channelpaving group
    channelpaving_parser = subparsers.add_parser('channelpaving', help='Channel paving utilities')
    channelpaving_subparsers = channelpaving_parser.add_subparsers(dest='channelpaving_cmd')

    channelpaving_add_depth_parser = channelpaving_subparsers.add_parser(
        'add-depth-to-polyline',
        help='Annotate polyline features with pt_depth values from NCFRIS cross-sections',
        parents=[channelpaving_add_depth_get_parser()],
        add_help=True
    )
    channelpaving_add_depth_parser.set_defaults(func=channelpaving_add_depth_main)

    channelpaving_add_width_parser = channelpaving_subparsers.add_parser(
        'add-width-to-polyline',
        help='Annotate polyline features with pt_width derived from NHDArea polygons',
        parents=[channelpaving_add_width_get_parser()],
        add_help=True
    )
    channelpaving_add_width_parser.set_defaults(func=channelpaving_add_width_main)

    channelpaving_generate_channel_mesh_parser = channelpaving_subparsers.add_parser(
        'generate-channel-mesh',
        help='Generate ADCIRC channel mesh from flowline data',
        parents=[channelpaving_generate_mesh_get_parser()],
        add_help=True
    )
    channelpaving_generate_channel_mesh_parser.set_defaults(func=channelpaving_generate_mesh_main)

    args = parser.parse_args()
    if hasattr(args, 'func'):
        return args.func(args)
    else:
        parser.print_help()
        return 1

if __name__ == '__main__':
    sys.exit(main()) 