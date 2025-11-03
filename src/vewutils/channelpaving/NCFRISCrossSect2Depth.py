import numpy as np
import geopandas as gpd
import fiona
from shapely.geometry import LineString
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.express as px
import argparse

def NCFRISCrossSect2Depth(ncfirs_xsect_file, ncfris_hydramodel_file, flowlines_file, area_coverage_file, flowlines_output_file, plot=False):
    print("Starting NCFRIS cross-section to depth processing...")
    print(f"Input files:")
    print(f"  - Cross-sections: {ncfirs_xsect_file}")
    print(f"  - Hydra model: {ncfris_hydramodel_file}")
    print(f"  - Flowlines: {flowlines_file}")
    if area_coverage_file:
        print(f"  - Area coverage: {area_coverage_file}")
    print(f"  - Output: {flowlines_output_file}")
    
    # Set the coordinate reference system used for the later process to UTM zone 17N. This is specific to the region of interest.
    crs_utm = 'epsg:32617'
    print(f"Using UTM CRS: {crs_utm}")

    # Obtain cross-sections
    print("Loading cross-sections...")
    with fiona.open(ncfirs_xsect_file) as src:
        gdf_xsect = gpd.GeoDataFrame.from_features(src)
    print(f"  Loaded {len(gdf_xsect)} cross-sections")

    # Obtain Hydra Model domains
    print("Loading hydra model domains...")
    with fiona.open(ncfris_hydramodel_file) as src:
        gdf_hydra = gpd.GeoDataFrame.from_features(src)
    print(f"  Loaded {len(gdf_hydra)} hydra model domains")

    # Obtain flowlines
    print("Loading flowlines...")
    gdf_flowlines = gpd.read_file(flowlines_file)
    print(f"  Loaded {len(gdf_flowlines)} flowlines")

    # Obtain area coverage
    if area_coverage_file:
        print("Loading area coverage...")
        area_coverage = gpd.read_file(area_coverage_file).to_crs(crs_utm).unary_union
        print("  Area coverage loaded")
    else:
        area_coverage = None
        print("  No area coverage file provided")

    # Obtain reduced versions of gdf_xsect and gdf_hydra
    print("Converting to UTM coordinate system...")
    if gdf_xsect.crs:
        gdf_xsect = gdf_xsect.to_crs(crs_utm)
    else:
        gdf_xsect = gdf_xsect.set_crs('epsg:2264').to_crs(crs_utm)
    if gdf_hydra.crs:
        gdf_hydra = gdf_hydra.to_crs(crs_utm)
    else:
        gdf_hydra = gdf_hydra.set_crs('epsg:2264').to_crs(crs_utm)
    gdf_flowlines = gdf_flowlines.to_crs(crs_utm)
    print("  Coordinate system conversion complete")
    
    # Reduce gdf_xsect and gdf_hydra to only those within the convex hull of gdf_flowlines
    print("Filtering data to area of interest...")
    if area_coverage:
        gdf_xsect_reduced = gdf_xsect[area_coverage.contains(gdf_xsect.geometry) | area_coverage.overlaps(gdf_xsect.geometry)]
        gdf_hydra_reduced = gdf_hydra[area_coverage.contains(gdf_hydra.geometry) | area_coverage.overlaps(gdf_hydra.geometry)]
        print(f"  Reduced to {len(gdf_xsect_reduced)} cross-sections and {len(gdf_hydra_reduced)} hydra domains in area of interest")
    else:
        gdf_xsect_reduced = gdf_xsect
        gdf_hydra_reduced = gdf_hydra
        print("  Using all data (no area filtering)")

    # Create pt_depth along flowlines from NCFRIS cross-sections
    print("Processing flowlines to add depth information...")
    ft2m = 0.3048
    gdf_flowlines['pt_depth'] = None
    total_flowlines = len(gdf_flowlines)
    print(f"  Processing {total_flowlines} flowlines...")
    
    for ifl in gdf_flowlines.index:
        if ifl%100 == 0:
            print(f'  Processing flowline {ifl+1}/{total_flowlines} ({(ifl+1)/total_flowlines*100:.1f}%)')
            
        fl = gdf_flowlines.loc[ifl]

        # Create segments from the LineString fl
        # Handle both LineString and MultiLineString geometries
        if fl.geometry.geom_type == 'LineString':
            # Simple LineString
            coords = list(fl.geometry.coords)
        else:
            # MultiLineString or other multi-part geometry - get coordinates from all parts
            coords = []
            for part in fl.geometry.geoms:
                coords.extend(list(part.coords))
        
        segments = [LineString([coords[i], coords[i + 1]]) for i in range(len(coords) - 1)]

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
                gdf_hydra_reduced_contains = gdf_hydra_reduced[contains]
                if np.count_nonzero(contains) > 1:
                    print(f"Warning: More than one HYDRAID contains the segment. All matched HYDRAIDs: {gdf_hydra_reduced_contains.HYDRAID.values}")
                gdf_hydra_reduced_contains_first = gdf_hydra_reduced_contains.iloc[0]
                target_hydraid = gdf_hydra_reduced_contains_first.HYDRAID
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
    
    print("  Flowline processing complete")
    print(f"  Successfully processed {total_flowlines} flowlines")
                
    # Save to file
    print(f"Saving results to {flowlines_output_file}...")
    gdf_flowlines.to_file(flowlines_output_file)
    print("  File saved successfully")
    
    # Plot
    if plot:
        print("Generating interactive plot...")
        gdf_flowlines = gdf_flowlines.to_crs('epsg:4326')

        lons_all = []
        lats_all = []
        point_depths_all = []
        d2 = []
        for jl in range(len(gdf_flowlines)):
            # Get coordinates for this geometry (handle both LineString and MultiLineString)
            if gdf_flowlines.loc[jl].geometry.geom_type == 'LineString':
                coords = list(gdf_flowlines.loc[jl].geometry.coords)
            else:
                # MultiLineString or other multi-part geometry
                coords = []
                for part in gdf_flowlines.loc[jl].geometry.geoms:
                    coords.extend(list(part.coords))
            
            if gdf_flowlines.loc[jl, 'pt_depth'] == None:
                point_depths = ['-99999.0'] * len(coords)
            else:
                point_depths = gdf_flowlines.loc[jl, 'pt_depth'].split(',')
                
            point_depths = [val if val != '-99999.0' else None for val in point_depths]

            point_depths_all.extend(point_depths)

            lons = []
            lats = []
            
            for il in range(len(coords)):
                lon = coords[il][0]  # x coordinate
                lat = coords[il][1]  # y coordinate
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
        print("  Interactive plot displayed")
    
    print("NCFRIS cross-section to depth processing completed successfully!")
        

def get_parser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--ncfirs-xsect-file', required=True, help='Path to NCFRIS cross-section file (e.g., GeoPackage/GeoJSON/Shapefile)')
    parser.add_argument('--ncfris-hydramodel-file', required=True, help='Path to NCFRIS hydra model domains file')
    parser.add_argument('--flowlines-file', required=True, help='Input flowlines file to annotate with pt_depth')
    parser.add_argument('--area-coverage-file', required=False, default=None, help='Optional polygon(s) defining area of interest')
    parser.add_argument('--flowlines-output-file', required=True, help='Output file for flowlines with pt_depth attribute')
    parser.add_argument('--plot', action='store_true', help='Show interactive map of results')
    return parser


def main(args):
    return NCFRISCrossSect2Depth(
        ncfirs_xsect_file=args.ncfirs_xsect_file,
        ncfris_hydramodel_file=args.ncfris_hydramodel_file,
        flowlines_file=args.flowlines_file,
        area_coverage_file=args.area_coverage_file,
        flowlines_output_file=args.flowlines_output_file,
        plot=args.plot,
    )
