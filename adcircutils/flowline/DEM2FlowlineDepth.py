import pandas as pd
import shapely
from shapely.geometry import LineString, Point, Polygon
import plotly.express as px
import plotly.graph_objects as go
import plotly.express as px
import geopandas as gpd
import pyproj
import numpy as np
import matplotlib.pyplot as plt
import sys
import rasterio
from pyproj import CRS
from pyproj import Transformer

def get_geotiff_minvalue_in_neighbor_at(fgeotiffs, nneipixel, lon, lat):
    for src in fgeotiffs:
        src_crs = src.crs
        bounds = src.bounds

        if src_crs != 'EPSG:4326':
            transformer = Transformer.from_crs(4326, src_crs, always_xy=True)
            xy = transformer.transform(lon, lat)
        else:
            xy = (lon, lat)

        # Get the transformation matrix to convert longitude and latitude to pixel coordinates
        col, row = ~src.transform * xy

        # Convert pixel coordinates to integer indices
        col_idx, row_idx = int(col), int(row)

        if col_idx < 0 or col_idx >= src.width or row_idx < 0 or row_idx >= src.height:
            # print('out of bounds ... lon: {}, lat: {}, xy: {}, col: {}, row: {}, bounds: {}'.format(lon, lat, xy, col, row, src.bounds))
            continue

        # Read the pixel value at the specified location
        r0 = max(0,row_idx-nneipixel)
        r1 = min(src.height, row_idx+nneipixel)
        c0 = max(0,col_idx-nneipixel)
        c1 = min(src.width, col_idx+nneipixel)
        value = src.read(1, window=((r0, r1), (c0, c1)))
        
        if np.any(np.isnan(value)):
            # print('None found for raster value at ({}, {}), ({}, {}). Continue to the next raster.'.format(lon, lat, col_idx, row_idx))
            continue

        min_value = np.min(value)

        if abs(min_value) < 1000.0:
            break

    return min_value

def DEM2FlowlineDepth(flowline_geo, geotiff_files, hydroflattened_depth, hydroflattened_depth_after_corrected, depth_min, nneipixel, min_window, output_geo):
    # read flowline data
    hydroflattened_depth_shift = hydroflattened_depth_after_corrected - hydroflattened_depth
    target_flowline = gpd.read_file(flowline_geo)
    flowlines = []

    for feature in target_flowline.geometry:
        if isinstance(feature, shapely.geometry.linestring.LineString):
            linestrings = [feature]
        elif isinstance(feature, shapely.geometry.multilinestring.MultiLineString):
            linestrings = feature.geoms
        else:
            continue

        flowlines.extend(linestrings)

    # open geotiff files
    fgeotiffs = [rasterio.open(geotiff_file) for geotiff_file in geotiff_files]

    # assign depths from DEMs
    for jl, flowline in enumerate(flowlines):
        if jl%1 == 0:
            print('Now at {}/{}'.format(jl+1, len(target_flowline)))

        depths = []
        for il in range(len(flowline.xy[0])):
            lon = flowline.xy[0][il]
            lat = flowline.xy[1][il]

            if lon == None or lat == None:
                raise Exception('Missing coordinates in flowline at (jl, il) = ({}, {})'.format(jl, il))

            depth = - get_geotiff_minvalue_in_neighbor_at(fgeotiffs, nneipixel, lon, lat)
            if depth <= hydroflattened_depth:
                depth = depth + hydroflattened_depth_shift
            depth = max(depth_min, depth)
            depths.append(depth)

        # taking minimum of neighbors
        depths_minnei = []
        for il in range(len(flowline.xy[0])):
            i0 = max(0, il - min_window)
            i1 = min(len(flowline.xy[0]), il + min_window)
            depths_minnei.append(np.min(depths[i0:i1]))

        target_flowline.at[jl,'pt_depth'] = ','.join(["{:.2f}".format(val) for val in depths_minnei])

    # close geotiff files
    for fgeotiff in fgeotiffs:
        fgeotiff.close()

    # Plot
    lons_all = []
    lats_all = []
    point_depths_all = []
    d2 = []
    for jl, flowline in enumerate(flowlines):
        point_depths = target_flowline.at[jl,'pt_depth'].split(',')
        point_depths_all.extend(point_depths)

        lons = []
        lats = []
        
        for il in range(len(flowline.xy[0])):
            lon = flowline.xy[0][il]
            lat = flowline.xy[1][il]
            lons.append(lon)
            lats.append(lat)

        lons_all.extend(lons)
        lats_all.extend(lats)

        dd2 = go.Scattermapbox(
                        lon=lons,
                        lat=lats,
                        mode='lines',
                        line=dict(color='navy',width=1),
                        showlegend=False)
        d2.append(dd2)

    point_depths_all_float = [float(text) for text in point_depths_all]

    d3 = go.Scattermapbox(
                        lon=lons_all,
                        lat=lats_all,
                        text=point_depths_all,
                        mode='markers',
                        marker=dict(color=[-val for val in point_depths_all_float],
                                    colorscale='Jet',
                                    colorbar=dict(title='Depth (m)')),
                        showlegend=False
                        )

    mapbox_access_token = 'pk.eyJ1Ijoic2J1bnlhIiwiYSI6ImNsZ2ZwMDVwMDAzbzIzbW55dWpvd2R6bmwifQ.Y4ejFHryUnGcB1a66mIuqg'
    zoom = 12
    mapbox_style = 'mapbox://styles/mapbox/streets-v11'
    fig = go.Figure(data=d2+[d3])
    # fig = go.Figure(data=[])
    layout = go.Layout(
        mapbox=dict(
            accesstoken=mapbox_access_token,
            center=dict(lat=np.mean(lats_all), lon=np.mean(lons_all)),
            style=mapbox_style,
            # style='open-street-map',
            zoom=zoom,
        ),
        width=1000,
        height=800,
    )
    fig.update_layout(layout)
    fig.show()

    # Save to file
    target_flowline.to_file(output_geo)
