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
from typing import List, Dict, Any


class PolylineToVEWConverter:
    """Class for converting polylines to VEW strings."""
    
    def __init__(self, mesh: AdcircMesh, polylines_gdf: gpd.GeoDataFrame, dist_max: float = 10.0):
        """
        Initialize the converter.
        
        Args:
            mesh: ADCIRC mesh object
            polylines_gdf: GeoDataFrame containing polylines
            dist_max: Maximum distance for nearest neighbor search in meters
        """
        self._mesh = mesh
        self._polylines_gdf = polylines_gdf
        self._dist_max = dist_max
        self._x = mesh.nodes.x
        self._y = mesh.nodes.y
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
            # Create dictionary with ordered fields
            node = {
                'node_id': int(node_id),
                'x': float(self._x[node_id]),
                'y': float(self._y[node_id]),
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


def main():
    """Main function to handle command line arguments and process the polylines."""
    parser = argparse.ArgumentParser(
        description="Convert polylines to VEW strings in ADCIRC meshes"
    )
    parser.add_argument(
        "meshfile",
        help="Path to the ADCIRC mesh file (fort.14)"
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
        default=1.0,
        help='Bank elevation for VEW strings (default: 1.0)'
    )
    parser.add_argument(
        '-n', '--mannings',
        type=float,
        default=0.02,
        help="Manning's n value for VEW strings (default: 0.02)"
    )
    
    args = parser.parse_args()
    
    # Read the mesh and polyline files
    mesh = AdcircMesh.open(args.meshfile)
    polylines_gdf = gpd.read_file(args.polylinefile)
    
    # Create converter and convert polylines to VEW strings
    converter = PolylineToVEWConverter(mesh, polylines_gdf, args.distance)
    vewstrings = converter.convert(args.elevation, args.mannings)
    
    # Save VEW strings to YAML file
    with open(args.output, 'w') as f:
        yaml.dump({'vewstrings': vewstrings}, f, default_flow_style=False, sort_keys=False)
    
    print(f"Successfully converted {len(vewstrings)} polylines to VEW strings")
    print(f"VEW strings saved to: {args.output}")
    
    return 0


if __name__ == "__main__":
    main() 