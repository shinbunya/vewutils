import numpy as np
import geopandas as gpd
import fiona
from shapely.geometry import LineString
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.express as px

def NCFRISCrossSect2Depth(ncfirs_xsect_file, ncfris_hydramodel_file, flowlines_file, area_coverage_file, flowlines_output_file, plot=False):
    # Set the coordinate reference system used for the later process to UTM zone 17N. This is specific to the region of interest.
    crs_utm = 'epsg:32617'

    # Obtain cross-sections
    with fiona.open(ncfirs_xsect_file) as src:
        gdf_xsect = gpd.GeoDataFrame.from_features(src)

    # Obtain Hydra Model domains
    with fiona.open(ncfris_hydramodel_file) as src:
        gdf_hydra = gpd.GeoDataFrame.from_features(src)

    # Obtain flowlines
    gdf_flowlines = gpd.read_file(flowlines_file)

    # Obtain area coverage
    area_coverage = gpd.read_file(area_coverage_file).to_crs(crs_utm).unary_union

    # Obtain reduced versions of gdf_xsect and gdf_hydra
    if gdf_xsect.crs:
        gdf_xsect = gdf_xsect.to_crs(crs_utm)
    else:
        gdf_xsect = gdf_xsect.set_crs('epsg:2264').to_crs(crs_utm)
    if gdf_hydra.crs:
        gdf_hydra = gdf_hydra.to_crs(crs_utm)
    else:
        gdf_hydra = gdf_hydra.set_crs('epsg:2264').to_crs(crs_utm)
    gdf_flowlines = gdf_flowlines.to_crs(crs_utm)
    
    
    # Reduce gdf_xsect and gdf_hydra to only those within the convex hull of gdf_flowlines
    gdf_xsect_reduced = gdf_xsect[area_coverage.contains(gdf_xsect.geometry) | area_coverage.overlaps(gdf_xsect.geometry)]
    gdf_hydra_reduced = gdf_hydra[area_coverage.contains(gdf_hydra.geometry) | area_coverage.overlaps(gdf_hydra.geometry)]

    # Create pt_depth along flowlines from NCFRIS cross-sections
    ft2m = 0.3048
    gdf_flowlines['pt_depth'] = None
    for ifl in gdf_flowlines.index:
        if ifl%100 == 0:
            print('Now at {}/{}'.format(ifl+1, len(gdf_flowlines)))
            
        fl = gdf_flowlines.loc[ifl]

        # Create segments from the LineString fl
        segments = [LineString([fl.geometry.coords[i], fl.geometry.coords[i + 1]]) for i in range(len(fl.geometry.coords) - 1)]

        # Create a GeoDataFrame from the segments
        gdf_segments = gpd.GeoDataFrame(geometry=segments, crs=crs_utm)

        # Identify the indexes of the crosssections in gdf_xsect that cross segments in LineString fl
        pt_dists = np.zeros(len(gdf_segments)+1)
        intersect_dists = []
        intersect_bed_elevs = []
        for i, segment in enumerate(gdf_segments.geometry):
            pt_dists[i+1] = pt_dists[i] + segment.length
            contains = gdf_hydra_reduced.contains(segment.centroid)
            if np.any(contains):
                if np.count_nonzero(contains) > 1:
                    raise ValueError("More than one HYDRAID contains the segment")
                target_hydraid = gdf_hydra_reduced[contains].HYDRAID.values[0]
                gdf_xsect_target = gdf_xsect_reduced[gdf_xsect_reduced.HYDRAID == target_hydraid]
                intersects = gdf_xsect_target.intersects(segment)
                if np.any(intersects):
                    intersection_points = gdf_xsect_target[intersects].intersection(segment)
                    if len(intersection_points) == 0:
                        raise ValueError("No intersection point")
                    for point in intersection_points:
                        bed_elev = gdf_xsect_target.loc[gdf_xsect_target[intersects].index, 'BED_ELEV'].values[0]
                        if bed_elev == -8888:
                            continue
                        dist = pt_dists[i] + segment.project(point)
                        intersect_dists.append(dist)
                        intersect_bed_elevs.append(bed_elev*ft2m)

        if len(intersect_dists) == 0:
            pt_depth = ','.join(['-99999.0'] * len(pt_dists))
        else:
            interpolated_bed_elevs = np.interp(pt_dists, intersect_dists, intersect_bed_elevs)
            pt_depth = ','.join([f'{val:.2f}' for val in interpolated_bed_elevs])

        gdf_flowlines.loc[ifl, 'pt_depth'] = pt_depth
                
    # Save to file
    gdf_flowlines.to_file(flowlines_output_file)
    
    # Plot
    if plot:
        gdf_flowlines = gdf_flowlines.to_crs('epsg:4326')

        lons_all = []
        lats_all = []
        point_depths_all = []
        d2 = []
        for jl in range(len(gdf_flowlines)):
            if gdf_flowlines.loc[jl, 'pt_depth'] == None:
                point_depths = ['-99999.0'] * len(gdf_flowlines.loc[jl].geometry.xy[0])
            else:
                point_depths = gdf_flowlines.loc[jl, 'pt_depth'].split(',')
                
            point_depths = [val if val != '-99999.0' else None for val in point_depths]

            point_depths_all.extend(point_depths)

            lons = []
            lats = []
            
            for il in range(len(gdf_flowlines.loc[jl].geometry.xy[0])):
                lon = gdf_flowlines.loc[jl].geometry.xy[0][il]
                lat = gdf_flowlines.loc[jl].geometry.xy[1][il]
                lons.append(lon)
                lats.append(lat)

            lons_all.extend(lons)
            lats_all.extend(lats)

        point_depths_all_float = [float(text) if text != None else 0.0 for text in point_depths_all]

        d3 = go.Scattermapbox(
                            lon=lons_all,
                            lat=lats_all,
                            text=point_depths_all,
                            mode='markers',
                            marker=dict(color=point_depths_all_float,
                                        colorscale='Jet',
                                        colorbar=dict(title='Bed Elev. (m)')),
                            showlegend=False
                            )

        zoom = 8.5
        fig = go.Figure(data=d2+[d3])
        # fig = go.Figure(data=[])
        layout = go.Layout(
            mapbox=dict(
                center=dict(lat=np.mean([np.min(lats_all), np.max(lats_all)]), 
                            lon=np.mean([np.min(lons_all), np.max(lons_all)])),
                style='open-street-map',
                zoom=zoom,
            ),
            width=1000,
            height=800,
        )
        fig.update_layout(layout)
        fig.show()
        