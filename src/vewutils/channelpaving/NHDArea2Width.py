import pandas as pd
import shapely
from shapely.geometry import LineString, Point, Polygon
import plotly.express as px
import plotly.graph_objects as go
import geopandas as gpd
import pyproj
import numpy as np
import matplotlib.pyplot as plt
import sys

import vewutils.channelpaving.utils as utils
import argparse


def NHDArea2Width(flowline_file, nhdarea_shpfiles, default_width, min_width, max_width, median_window, output_file, nhdplusids=None, plot=False):
    print("Starting NHD area to width processing...")
    print(f"Input files:")
    print(f"  - Flowlines: {flowline_file}")
    print(f"  - NHD area files: {len(nhdarea_shpfiles)} file(s)")
    for i, (path, layer) in enumerate(nhdarea_shpfiles):
        if layer:
            print(f"    {i+1}. {path} (layer: {layer})")
        else:
            print(f"    {i+1}. {path}")
    if nhdplusids is not None and len(nhdplusids) > 0:
        print(f"  - NHDPlusIDs: {nhdplusids}")
    else:
        print(f"  - NHDPlusIDs: Not specified (will use all areas)")
    print(f"  - Output: {output_file}")
    print(f"Parameters:")
    print(f"  - Default width: {default_width} m")
    print(f"  - Min width: {min_width} m")
    print(f"  - Max width: {max_width} m")
    print(f"  - Median window: {median_window}")
    
    # Flowlines
    print("Loading flowlines...")
    target_flowline = gpd.read_file(flowline_file)
    target_flowline.to_crs(pyproj.CRS.from_epsg(4326), inplace=True)
    flowlines = []

    for feature in target_flowline.geometry:
        if isinstance(feature, shapely.geometry.linestring.LineString):
            linestrings = [feature]
        elif isinstance(feature, shapely.geometry.multilinestring.MultiLineString):
            linestrings = feature.geoms
        else:
            continue

        flowlines.extend(linestrings)

    print(f"  Loaded {len(flowlines)} flowline segments from {len(target_flowline)} features")

    flowline = flowlines[0]
    if len(flowline.xy) == 2:
        lon = flowline.xy[0][0]
        lat = flowline.xy[1][0]
    else:
        lon = flowline.xy[0][0]
        lat = flowline.xy[0][1]

    print(f"Determining UTM CRS for location ({lon:.6f}, {lat:.6f})...")
    utm_crs_list = pyproj.database.query_utm_crs_info(
        datum_name="WGS 84",
        area_of_interest=pyproj.aoi.AreaOfInterest(
            west_lon_degree=lon,
            south_lat_degree=lat,
            east_lon_degree=lon,
            north_lat_degree=lat,
        ),
    )
    utm_crs = pyproj.CRS.from_epsg(utm_crs_list[0].code)
    print(f"  Using UTM CRS: {utm_crs}")

    # NHDArea polygons
    print("Loading NHD area polygons...")
    gdf_areas = []
    for i, nhdarea_shpfile in enumerate(nhdarea_shpfiles):
        print(f"  Loading file {i+1}/{len(nhdarea_shpfiles)}: {nhdarea_shpfile[0]}")
        if nhdarea_shpfile[1]:
            gdf_areas.append(gpd.read_file(nhdarea_shpfile[0], layer=nhdarea_shpfile[1]))
            print(f"    Loaded layer '{nhdarea_shpfile[1]}'")
        else:
            gdf_areas.append(gpd.read_file(nhdarea_shpfile[0]))
    
    print("  Combining all area polygons...")
    gdf_area = gpd.GeoDataFrame(pd.concat(gdf_areas, ignore_index=True))
    gdf_area.to_crs(pyproj.CRS.from_epsg(4326), inplace=True)
    print(f"  Total area polygons loaded: {len(gdf_area)}")
    
    print("Filtering area polygons by NHDPlusID...")
    if nhdplusids is None or len(nhdplusids) == 0:
        # Use all area polygons if no NHDPlusIDs specified
        target_area = gdf_area
        print(f"  No NHDPlusIDs specified, using all {len(target_area)} area polygons")
    else:
        if 'NHDPlusID' in gdf_area.columns:
            target_area = gdf_area[
                gdf_area['NHDPlusID'].isin(nhdplusids)
            ]
            print(f"  Using 'NHDPlusID' column")
        elif 'permanent_' in gdf_area.columns:
            target_area = gdf_area[
                gdf_area['permanent_'].isin(nhdplusids)
            ]
            print(f"  Using 'permanent_' column")
        else:
            raise Exception('No NHDPlusID or permanent_ identifier found in the area polygon file.')
        
        if len(target_area) == 0:
            print('  No area polygon found. Check the NHDPlusID.')
            return
        print(f"  Found {len(target_area)} matching area polygons")
    print("Processing polygon geometries...")
    area_polygons_lonlat = []
    for mpoly in target_area.geometry:
        if type(mpoly) == Polygon:
            area_polygons_lonlat.append(mpoly)
        else:
            for poly in mpoly.geoms:
                area_polygons_lonlat.append(poly)
    print(f"  Extracted {len(area_polygons_lonlat)} individual polygons")

    print("Converting polygons to UTM coordinates...")
    area_polygonss_utm = []
    for i, area_polygon_lonlat in enumerate(area_polygons_lonlat):
        if i % 10 == 0:
            print(f"    Converting polygon {i+1}/{len(area_polygons_lonlat)}")
        pts = utils.lonlat2utm(area_polygon_lonlat.exterior.coords.xy[0],area_polygon_lonlat.exterior.coords.xy[1],utm_crs)
        area_polygons_utm = [Polygon(list(zip(pts[0],pts[1]))).exterior]
        for area_polygon_in_lonlat in area_polygon_lonlat.interiors:
            pts = utils.lonlat2utm(area_polygon_in_lonlat.coords.xy[0],area_polygon_in_lonlat.coords.xy[1],utm_crs)
            area_polygons_utm.append(Polygon(list(zip(pts[0],pts[1]))).exterior)
        area_polygonss_utm.append(area_polygons_utm)
    print("  UTM conversion complete")

    # Find distances from center lines to water body polygon boundaries
    print("Processing flowlines to calculate widths...")
    intersection_points1 = []
    widths = []
    cpoints = []
    total_flowlines = len(target_flowline)
    print(f"  Processing {total_flowlines} flowlines...")

    print_freq = max(1, min(100, len(flowlines)/100))
    
    for jl, flowline in enumerate(flowlines):
        if jl%print_freq == 0:
            print(f'  Processing flowline {jl+1}/{total_flowlines} ({(jl+1)/total_flowlines*100:.1f}%)')

        widths_jl = []
        
        for il in range(len(flowline.xy[0])):
            lon0 = flowline.xy[0][il]
            lat0 = flowline.xy[1][il]

            if lon0 == None or lat0 == None:
                continue

            clon = lon0
            clat = lat0
            cpoint = Point(clon, clat)
            cpoint_utm = Point(utils.lonlat2utm(cpoint.x, cpoint.y, utm_crs))

            cpoints.append(cpoint)

            found = False
            for ip in range(len(area_polygons_lonlat)):
                area_polygon_lonlat = area_polygons_lonlat[ip]
                if cpoint.within(area_polygon_lonlat):
                    area_polygons_utm = area_polygonss_utm[ip]

                    p1s, p2s = shapely.ops.nearest_points(area_polygons_utm, cpoint_utm)
                    dists = [shapely.distance(p1i, cpoint_utm) for p1i in p1s]
                    min_index = np.argmin(dists)
                    p1_utm = p1s[min_index]

                    p1 = Point(utils.utm2lonlat(p1_utm.x, p1_utm.y, utm_crs))
                    intersection_points1.append(p1)

                    width = dists[min_index]*2.0
                    
                    width = max(min_width, min(width, max_width))

                    found = True
                    break
            
            if found:
                widths_jl.append(width)
            else:
                if jl % 100 == 0:  # Only print for every 100th flowline to avoid spam
                    print(f'    No intersection found at (jl, il) = ({jl}, {il}). Using default width.')
                widths_jl.append(default_width)
        
        width_median = np.median(widths_jl)
        target_flowline.at[jl,'width'] = width_median

        widths_mediannei = []
        for il in range(len(flowline.xy[0])):
            i0 = max(0, il - median_window)
            i1 = min(len(flowline.xy[0]), il + median_window)
            widths_mediannei.append(np.median(widths_jl[i0:i1]))

        target_flowline.at[jl,'pt_width'] = ','.join(["{:.2f}".format(val) for val in widths_mediannei])

        widths.extend(widths_mediannei)

    print("  Flowline processing complete")
    print(f"  Successfully processed {total_flowlines} flowlines")
    print(f"  Calculated widths for {len(widths)} points")

    # Save to file
    print(f"Saving results to {output_file}...")
    target_flowline.to_file(output_file)
    print("  File saved successfully")

    # Plot
    if plot:
        print("Generating interactive plot...")
        d1 = []
        # for poly in area_polygons_lonlat:
        #     dd1 = go.Scattermapbox(lon=list(poly.exterior.coords.xy[0]),
        #                         lat=list(poly.exterior.coords.xy[1]),
        #                         fill='toself',
        #                         mode='lines',
        #                         line=dict(color=None,width=0),
        #                         # fillcolor='rgba(28,163,236,0.5)'
        #                         fillcolor='orange'
        #                         )
        #     d1.append(dd1)

        d2 = []
        # for line in flowlines:
        #     dd2 = go.Scattermapbox(
        #                     lon=list(line.xy[0]),
        #                     lat=list(line.xy[1]),
        #                     mode='lines',
        #                     line=dict(color='navy',width=1),
        #                     showlegend=False,
        #                     )
        #     d2.append(dd2)

        d3 = go.Scattermapbox(lon=[p.x for p in intersection_points1],
                            lat=[p.y for p in intersection_points1],
                            mode='markers',
                            showlegend=False,
                            )

        d5 = go.Scattermapbox(lon=[p.x for p in cpoints],
                            lat=[p.y for p in cpoints],
                            mode='markers',
                            marker=dict(color=widths, colorscale='Jet',
                                        colorbar=dict(title='Width (m)')),
                            showlegend=False,
                            text=["{:.1f}".format(x) for x in widths],
                            hoverinfo='text',
                            )
        px = [p.x for p in cpoints]
        py = [p.y for p in cpoints]
        center_lon = np.mean([np.min(px), np.max(px)])
        center_lat = np.mean([np.min(py), np.max(py)])
        zoom = 8.5
        fig = go.Figure(data=(d1+d2+[d3,d5]))
        layout = go.Layout(
            mapbox=dict(
                center=dict(lat=center_lat, lon=center_lon),
                style='open-street-map',
                zoom=zoom,
            ),
            width=1000,
            height=800,
        )
        fig.update_layout(layout)
        fig.show()
        print("  Interactive plot displayed")
    
    print("NHD area to width processing completed successfully!")


def get_parser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--flowline-file', required=True, help='Input polyline(s) file to annotate')
    parser.add_argument('--nhdarea-shpfiles', required=True, nargs='+', help='One or more NHDArea files. Use layer syntax path::layer for GPKG')
    parser.add_argument('--nhdplusids', required=False, nargs='*', help='One or more NHDPlusID/permanent_ values to select areas. If not specified, all areas will be used.')
    parser.add_argument('--default-width', required=True, type=float, help='Default width (m) when no polygon found')
    parser.add_argument('--min-width', required=True, type=float, help='Minimum width clamp (m)')
    parser.add_argument('--max-width', required=True, type=float, help='Maximum width clamp (m)')
    parser.add_argument('--median-window', required=True, type=int, help='Window half-size for neighbor median smoothing')
    parser.add_argument('--output-file', required=True, help='Output file path with pt_width attribute')
    parser.add_argument('--plot', action='store_true', help='Show interactive map of results')
    return parser


def _parse_nhdarea_args(nhdarea_shpfiles_args):
    parsed = []
    for token in nhdarea_shpfiles_args:
        if '::' in token:
            path, layer = token.split('::', 1)
            parsed.append((path, layer))
        else:
            parsed.append((token, None))
    return parsed


def main(args):
    nhdarea_pairs = _parse_nhdarea_args(args.nhdarea_shpfiles)
    
    # Handle optional nhdplusids
    if args.nhdplusids is not None and len(args.nhdplusids) > 0:
        # Accept numeric ids or strings transparently
        def _coerce_id(x):
            try:
                return int(x)
            except Exception:
                return x
        nhdplusids = [_coerce_id(x) for x in args.nhdplusids]
    else:
        nhdplusids = None

    return NHDArea2Width(
        flowline_file=args.flowline_file,
        nhdarea_shpfiles=nhdarea_pairs,
        default_width=args.default_width,
        min_width=args.min_width,
        max_width=args.max_width,
        median_window=args.median_window,
        output_file=args.output_file,
        nhdplusids=nhdplusids,
        plot=args.plot,
    )
