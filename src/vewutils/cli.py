#!/usr/bin/env python3
import argparse
import sys

# Import get_parser and main for each command
from vewutils.mesh.mesh_merger import get_parser as mesh_merge_get_parser, main as mesh_merge_main
from vewutils.mesh.mesh_subtractor import get_parser as mesh_subtract_get_parser, main as mesh_subtract_main
from vewutils.mesh.boundary_exporter import get_parser as mesh_export_boundaries_get_parser, main as mesh_export_boundaries_main
from vewutils.mesh.add_land_boundaries import get_parser as mesh_add_land_boundaries_get_parser, main as mesh_add_land_boundaries_main
from vewutils.mesh.adjust_vew_barrier_heights import get_parser as mesh_adjust_vew_barrier_heights_get_parser, main as mesh_adjust_vew_barrier_heights_main
from vewutils.mesh.adjust_vew_channel_elevations import get_parser as mesh_adjust_vew_channel_elevations_get_parser, main as mesh_adjust_vew_channel_elevations_main
from vewutils.plot.plot_hydrograph_at_station import get_parser as plot_hydrograph_get_parser, main as plot_hydrograph_main
from vewutils.plot.plot_errorhistogram_at_station import get_parser as plot_errorhistogram_get_parser, main as plot_errorhistogram_main
from vewutils.post.maxele_max import get_parser as post_maxele_max_get_parser, main as post_maxele_max_main
from vewutils.post.maxele_diff import get_parser as post_maxele_diff_get_parser, main as post_maxele_diff_main
from vewutils.post.maxele_add_disturbance import get_parser as post_maxele_add_disturbance_get_parser, main as post_maxele_add_disturbance_main
from vewutils.post.maxele_attribution import get_parser as post_maxele_attribution_get_parser, main as post_maxele_attribution_main
from vewutils.post.reduce_timesteps import get_parser as post_reduce_timesteps_get_parser, main as post_reduce_timesteps_main
from vewutils.dem2adcdp.dem2adcdp import get_parser as dem2adcdp_get_parser, main as dem2adcdp_main
from vewutils.nodalattribute.attribute_transfer import get_parser as nodalattribute_transfer_get_parser, main as nodalattribute_transfer_main
from vewutils.nodalattribute.manningsn_extractor import get_parser as nodalattribute_extract_manningsn_get_parser, main as nodalattribute_extract_manningsn_main
from vewutils.vewprocessing.vew_adder import get_parser as vewprocessing_add_get_parser, main as vewprocessing_add_main
from vewutils.vewprocessing.polyline_converter import get_parser as vewprocessing_convert_polylines_get_parser, main as vewprocessing_convert_polylines_main
from vewutils.vewprocessing.vew_scraper import get_parser as vewprocessing_scrape_get_parser, main as vewprocessing_scrape_main
from vewutils.utils.node_selector import get_parser as utils_select_nodes_get_parser, main as utils_select_nodes_main

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

    # Post group
    post_parser = subparsers.add_parser('post', help='Postprocessing commands')
    post_subparsers = post_parser.add_subparsers(dest='post_cmd')
    post_maxele_max_parser = post_subparsers.add_parser(
        'maxele-max',
        help='Calculate max between two maxele files',
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
    post_maxele_attr_parser = post_subparsers.add_parser(
        'maxele-attribution',
        help='Add attribution to maxele file',
        parents=[post_maxele_attribution_get_parser()],
        add_help=True
    )
    post_maxele_attr_parser.set_defaults(func=post_maxele_attribution_main)
    post_reduce_ts_parser = post_subparsers.add_parser(
        'reduce-timesteps',
        help='Reduce time steps in NetCDF',
        parents=[post_reduce_timesteps_get_parser()],
        add_help=True
    )
    post_reduce_ts_parser.set_defaults(func=post_reduce_timesteps_main)

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

    args = parser.parse_args()
    if hasattr(args, 'func'):
        return args.func(args)
    else:
        parser.print_help()
        return 1

if __name__ == '__main__':
    sys.exit(main()) 