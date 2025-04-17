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

import adcircutils.channelpaving.utils as utils


def NHDArea2Width(flowline_file, nhdarea_shpfiles, nhdplusids, default_width, min_width, max_width, median_window, output_file, plot=False):
    # Flowlines
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

    flowline = flowlines[0]
    if len(flowline.xy) == 2:
        lon = flowline.xy[0][0]
        lat = flowline.xy[1][0]
    else:
        lon = flowline.xy[0][0]
        lat = flowline.xy[0][1]

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

    # NHDArea polygons
    gdf_areas = []
    for nhdarea_shpfile in nhdarea_shpfiles:
        if nhdarea_shpfile[1]:
            gdf_areas.append(gpd.read_file(nhdarea_shpfile[0], layer=nhdarea_shpfile[1]))
        else:
            gdf_areas.append(gpd.read_file(nhdarea_shpfile[0]))
    gdf_area = gpd.GeoDataFrame(pd.concat(gdf_areas, ignore_index=True))
    gdf_area.to_crs(pyproj.CRS.from_epsg(4326), inplace=True)
    if 'NHDPlusID' in gdf_area.columns:
        target_area = gdf_area[
            gdf_area['NHDPlusID'].isin(nhdplusids)
        ]
    elif 'permanent_' in gdf_area.columns:
        target_area = gdf_area[
            gdf_area['permanent_'].isin(nhdplusids)
        ]
    if len(target_area) == 0:
        print('No area polygon found. Check the NHDPlusID.')
        return
    area_polygons_lonlat = []
    for mpoly in target_area.geometry:
        if type(mpoly) == Polygon:
            area_polygons_lonlat.append(mpoly)
        else:
            for poly in mpoly.geoms:
                area_polygons_lonlat.append(poly)

    area_polygonss_utm = []
    for area_polygon_lonlat in area_polygons_lonlat:
        pts = utils.lonlat2utm(area_polygon_lonlat.exterior.coords.xy[0],area_polygon_lonlat.exterior.coords.xy[1],utm_crs)
        area_polygons_utm = [Polygon(list(zip(pts[0],pts[1]))).exterior]
        for area_polygon_in_lonlat in area_polygon_lonlat.interiors:
            pts = utils.lonlat2utm(area_polygon_in_lonlat.coords.xy[0],area_polygon_in_lonlat.coords.xy[1],utm_crs)
            area_polygons_utm.append(Polygon(list(zip(pts[0],pts[1]))).exterior)
        area_polygonss_utm.append(area_polygons_utm)

    # Find distances from center lines to water body polygon boundaries
    intersection_points1 = []
    widths = []
    cpoints = []

    for jl, flowline in enumerate(flowlines):
        if jl%100 == 0:
            print('Now at {}/{}'.format(jl+1, len(target_flowline)))

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
                print('No intersection found at (jl, il) = ({}, {}). The default width will be used.'.format(jl, il))
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

    # Save to file
    target_flowline.to_file(output_file)

    # Plot
    if plot:
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
