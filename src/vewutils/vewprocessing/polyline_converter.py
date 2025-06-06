#!/usr/bin/env python3
"""
Module for converting polylines to VEW strings in ADCIRC meshes.
"""

import argparse
import yaml
import numpy as np
import geopandas as gpd
from shapely.geometry import Point, LineString
from scipy.spatial import cKDTree
from adcircpy import AdcircMesh
from typing import List, Dict, Any, Tuple
import pyproj
from pyproj import CRS, Transformer
import pandas as pd


def get_utm_crs(lon: float, lat: float) -> CRS:
    """
    Determine the UTM CRS based on a longitude/latitude point.
    
    Args:
        lon: Longitude of the point
        lat: Latitude of the point
        
    Returns:
        CRS object for the appropriate UTM zone
    """
    utm_zone = int((lon + 180) / 6) + 1
    hemisphere = 'north' if lat >= 0 else 'south'
    utm_crs = CRS.from_dict({
        'proj': 'utm',
        'zone': utm_zone,
        'hemisphere': hemisphere,
        'ellps': 'WGS84',
        'datum': 'WGS84',
        'units': 'm'
    })
    return utm_crs


def transform_mesh_coordinates(x: pd.Series, y: pd.Series, from_crs: CRS, to_crs: CRS) -> Tuple[pd.Series, pd.Series]:
    """
    Transform mesh coordinates from one CRS to another.
    
    Args:
        x: Series of x coordinates
        y: Series of y coordinates
        from_crs: Source CRS
        to_crs: Target CRS
        
    Returns:
        Tuple of transformed (x, y) coordinates as Series
    """
    transformer = Transformer.from_crs(from_crs, to_crs, always_xy=True)
    x_trans, y_trans = transformer.transform(x.values, y.values)
    return pd.Series(x_trans, index=x.index), pd.Series(y_trans, index=y.index)


class PolylineToVEWConverter:
    """Class for converting polylines to VEW strings."""
    
    def __init__(self, mesh: AdcircMesh, polylines_gdf: gpd.GeoDataFrame, mesh_crs: CRS = None, dist_max: float = 10.0):
        """
        Initialize the converter.
        
        Args:
            mesh: ADCIRC mesh object
            mesh_crs: CRS of the mesh coordinates. If None, assumes coordinates are already in meters
            polylines_gdf: GeoDataFrame containing polylines
            dist_max: Maximum distance for nearest neighbor search in meters
        """
        self._mesh = mesh
        self._dist_max = dist_max
        
        # Store original coordinates
        self._x_orig = mesh.nodes.x
        self._y_orig = mesh.nodes.y
        
        # Get the center of all polylines
        polylines_center = polylines_gdf.geometry.unary_union.centroid
        
        # If mesh_crs is provided, use it to transform the polylines and mesh to UTM
        if mesh_crs is not None:
            # Determine source CRS from the GeoDataFrame
            source_crs = polylines_gdf.crs
            if source_crs is None:
                raise ValueError("Input polylines must have a defined coordinate reference system (CRS)")
                
            # Get UTM CRS based on polylines center
            utm_crs = get_utm_crs(polylines_center.x, polylines_center.y)
        
            # Transform mesh coordinates from specified CRS to UTM
            self._x, self._y = transform_mesh_coordinates(
                self._x_orig, self._y_orig,
                mesh_crs, utm_crs
            )
        
            # Transform polylines to UTM
            self._polylines_gdf = polylines_gdf.to_crs(utm_crs)
        
        else:
            # If mesh_crs is None, assume coordinates are in meters and copy original coordinates
            self._x = self._x_orig.copy()
            self._y = self._y_orig.copy()
        
            self._polylines_gdf = polylines_gdf.copy()
        
        # Initialize KD-tree with transformed coordinates
        self._tree = cKDTree(np.c_[self._x.values, self._y.values])
        self._neighs = mesh.node_neighbors.copy()
        
    def _find_nearest_node(self, x: float, y: float) -> int:
        """Find the nearest mesh node to a point."""
        distance, nearest_node_index = self._tree.query([x, y])
        return nearest_node_index + 1  # Convert to 1-based indexing
        
    def _create_nodestring(self, line) -> List[int]:
        """Create a nodestring along a line."""
        slx, sly = line.coords[0][0], line.coords[0][1]
        elx, ely = line.coords[-1][0], line.coords[-1][1]
        
        # Find nearest node to start point
        ni = self._find_nearest_node(slx, sly)
        xi = self._x.iloc[ni-1]
        yi = self._y.iloc[ni-1]
        
        nodestring = [ni]
        ipos = 0
        
        while ipos < len(line.coords) - 1:
            # Get the points in the line segment
            point0 = line.coords[ipos]
            point1 = line.coords[ipos + 1]
            xl0, yl0 = point0[0], point0[1]
            xl1, yl1 = point1[0], point1[1]
            
            # Create a line segment between consecutive points
            line_segment = LineString([(xl0, yl0), (xl1, yl1)])
            
            neigh = self._neighs[ni]
            neigh = [n for n in neigh if n not in nodestring]
            
            if not neigh:  # If no unvisited neighbors, move to next segment
                ipos += 1
                continue
            
            xnei = [self._x.loc[n] for n in neigh]
            ynei = [self._y.loc[n] for n in neigh]
            
            # Compute distances from line segment to neighbor points
            distances = [line_segment.distance(Point(xn, yn)) for xn, yn in zip(xnei, ynei)]
            min_distance = np.min(distances) if distances else float('inf')
            
            if min_distance > self._dist_max:
                ipos += 1
                continue
            
            min_distance_index = np.argmin(distances)
            ni = neigh[min_distance_index]
            xi = self._x[ni]
            yi = self._y[ni]
            
            nodestring.append(int(ni))
            
            distance_to_endpoint = np.sqrt((elx - xi)**2 + (ely - yi)**2)
            
            if distance_to_endpoint < self._dist_max:
                break
                
        # Check if the nodestring should be closed
        xi, yi = self._x.loc[nodestring[0]], self._y.loc[nodestring[0]]
        distance_to_endpoint = np.sqrt((elx - xi)**2 + (ely - yi)**2)
        if distance_to_endpoint < self._dist_max:
            nodestring.append(int(nodestring[0]))
            
        return nodestring
        
    def _create_vewstring(self, nodestring: List[int], bank_elevation: float = 1.0, bank_mannings_n: float = 0.02) -> List[Dict[str, Any]]:
        """Create a VEW string from a nodestring."""
        vewstring = []
        for node_id in nodestring:
            # Create dictionary with ordered fields using original coordinates
            node = {
                'node_id': int(node_id),
                'x': float(self._x_orig[node_id]),
                'y': float(self._y_orig[node_id]),
                'bank_elevation': float(bank_elevation),
                'bank_mannings_n': float(bank_mannings_n)
            }
            vewstring.append(node)
        return vewstring
        
    def convert(self, bank_elevation: float = 1.0, bank_mannings_n: float = 0.02) -> List[List[Dict[str, Any]]]:
        """
        Convert polylines to VEW strings.
        
        Args:
            bank_elevation: Bank elevation for VEW strings
            bank_mannings_n: Manning's n value for VEW strings
            
        Returns:
            List of VEW strings
        """
        vewstrings = []
        for _, line in self._polylines_gdf.iterrows():
            nodestring = self._create_nodestring(line.geometry)
            if len(nodestring) > 1:  # Only keep strings with at least 2 nodes
                vewstring = self._create_vewstring(nodestring, bank_elevation, bank_mannings_n)
                vewstrings.append(vewstring)
        return vewstrings


def get_parser():
    parser = argparse.ArgumentParser(add_help=False, description="Convert polylines to VEW strings in ADCIRC meshes")
    parser.add_argument(
        "meshfile",
        help="Path to the ADCIRC mesh file (fort.14)"
    )
    parser.add_argument(
        '-c', '--mesh-crs',
        help='EPSG code or PROJ string for mesh coordinates (e.g., EPSG:4326 for WGS84). If not provided, assumes coordinates are in meters',
        default=None
    )
    parser.add_argument(
        "polylinefile",
        help="Path to the polyline file (shapefile, geojson, etc.)"
    )
    parser.add_argument(
        '-o', '--output',
        default='vewstrings.yaml',
        help='Path to save the VEW strings (default: vewstrings.yaml)'
    )
    parser.add_argument(
        '-d', '--distance',
        type=float,
        default=10.0,
        help='Maximum distance for nearest neighbor search in meters (default: 10.0)'
    )
    parser.add_argument(
        '-e', '--elevation',
        type=float,
        default=-99999.0,
        help='Bank elevation for VEW strings (default: -99999.0 to use mesh elevation)'
    )
    parser.add_argument(
        '-n', '--mannings-n',
        type=float,
        default=0.02,
        help="Manning's n value for VEW strings (default: 0.02)"
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    # Read the mesh and polyline files
    mesh = AdcircMesh.open(args.meshfile)
    polylines_gdf = gpd.read_file(args.polylinefile)
    
    # Parse mesh CRS if provided
    mesh_crs = None
    if args.mesh_crs is not None:
        try:
            mesh_crs = CRS.from_string(args.mesh_crs)
        except Exception as e:
            print(f"Error parsing mesh CRS: {e}")
            print("Please provide a valid EPSG code (e.g., EPSG:4326) or PROJ string")
            return 1
    
    # Create converter and convert polylines to VEW strings
    converter = PolylineToVEWConverter(mesh, polylines_gdf, mesh_crs, args.distance)
    vewstrings = converter.convert(args.elevation, args.mannings_n)
    
    # Save VEW strings to YAML file
    with open(args.output, 'w') as f:
        yaml.dump({'vewstrings': vewstrings}, f, default_flow_style=False, sort_keys=False)
    
    print(f"Successfully converted {len(vewstrings)} polyline(s) to VEW string(s)")
    print(f"VEW strings saved to: {args.output}")
    
    return 0


if __name__ == "__main__":
    main() 