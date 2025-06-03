#!/usr/bin/env python3
"""
Module for exporting ADCIRC mesh boundaries to GeoPackage (.gpkg) format.
Each linestring will have its corresponding IBTYPE value as an attribute.
"""

import argparse
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import LineString
from typing import Dict, List, Tuple
from adcircpy import AdcircMesh


class BoundaryExporter:
    """Class for exporting ADCIRC mesh boundaries to GeoPackage format."""

    def __init__(self, mesh: AdcircMesh):
        """Initialize with an ADCIRC mesh."""
        self._mesh = mesh
        print("Initializing boundary exporter...")

    def _get_boundary_nodes_coordinates(self) -> Dict[int, Tuple[float, float]]:
        """
        Get coordinates for all nodes in the mesh.
        
        Returns:
            Dictionary mapping node IDs to (x, y) coordinates
        """
        print("  Extracting node coordinates from mesh...")
        
        # Get the nodes dataframe
        nodes_df = self._mesh.nodes
        
        # Check if nodes dataframe is empty
        if nodes_df.empty:
            print("  Warning: Mesh nodes dataframe is empty")
            return {}
            
        # Initialize a dictionary to store node coordinates
        node_coords = {}
        
        # Debug information
        print(f"  Node dataframe shape: {nodes_df.shape}")
        print(f"  Node dataframe columns: {nodes_df.columns.tolist()}")
        print(f"  Node dataframe index name: {nodes_df.index.name}")
        
        # Extract the first few rows to understand the data structure
        sample_rows = nodes_df.head(3)
        print(f"  Sample node data (first 3 rows):\n{sample_rows}")
        
        try:
            # The dataframe appears to have node IDs as the index and 'x', 'y' as columns
            # Extract coordinates directly using column names
            for node_id, row in nodes_df.iterrows():
                try:
                    # Convert to appropriate types
                    node_id_int = int(node_id)
                    x = float(row['x'])
                    y = float(row['y'])
                    node_coords[node_id_int] = (x, y)
                except (ValueError, TypeError, KeyError) as e:
                    if len(node_coords) < 5:  # Only print for the first few rows to avoid flooding the console
                        print(f"  Warning: Could not process node {node_id}: {e}")
            
            print(f"  Successfully extracted coordinates for {len(node_coords)} nodes")
            
            # Verify if we got all nodes
            if len(node_coords) < len(nodes_df) * 0.9:  # If we got less than 90% of the nodes
                print(f"  Warning: Only extracted {len(node_coords)} coordinates out of {len(nodes_df)} nodes")
                
                # Try alternative approach if the first method didn't work well
                if len(node_coords) < 1000:  # If we got very few nodes
                    print("  Trying alternative method to extract coordinates...")
                    node_coords.clear()
                    
                    # Try using direct attribute access
                    if hasattr(self._mesh, 'x') and hasattr(self._mesh, 'y'):
                        print("  Using mesh.x and mesh.y attributes...")
                        try:
                            x_values = self._mesh.x
                            y_values = self._mesh.y
                            
                            # Check if the mesh node indices are 0-based or 1-based
                            # ADCIRC typically uses 1-based indexing
                            base_index = 1
                            
                            # Get sample boundary node to check indexing
                            sample_boundary_node = None
                            for ibtype, boundaries in self._mesh.boundaries.to_dict().items():
                                if boundaries:
                                    for boundary in boundaries:
                                        node_ids = boundary['node_id']
                                        if node_ids and not isinstance(node_ids[0], tuple):
                                            sample_boundary_node = int(node_ids[0])
                                            break
                                    if sample_boundary_node:
                                        break
                            
                            # If we found a sample boundary node, check if it's within array range
                            if sample_boundary_node is not None:
                                if sample_boundary_node > len(x_values):
                                    print(f"  Warning: Sample boundary node {sample_boundary_node} is out of range for coordinate arrays (length {len(x_values)})")
                                    # Node ID might not be continuous, try using node dataframe
                                    for node_id in nodes_df.index:
                                        node_coords[int(node_id)] = (float(nodes_df.loc[node_id, 'x']), 
                                                                     float(nodes_df.loc[node_id, 'y']))
                                else:
                                    # Handle 0-based or 1-based indexing
                                    base_index = 1 if sample_boundary_node >= 1 else 0
                                    print(f"  Using {base_index}-based indexing for node IDs")
                                    
                                    # Map node IDs to coordinates
                                    for i in range(len(x_values)):
                                        node_id = i + base_index
                                        node_coords[node_id] = (float(x_values[i]), float(y_values[i]))
                            else:
                                # If no sample boundary node, just use 1-based indexing (typical for ADCIRC)
                                for i in range(len(x_values)):
                                    node_id = i + 1
                                    node_coords[node_id] = (float(x_values[i]), float(y_values[i]))
                            
                            print(f"  Alternative method: extracted {len(node_coords)} node coordinates")
                        except Exception as e:
                            print(f"  Error using direct attribute access: {e}")
                    
                    # Try a third method using row iteration and iloc
                    if len(node_coords) < 1000:
                        print("  Trying third method to extract coordinates...")
                        try:
                            # If nodes_df has an index that doesn't start at 1, 
                            # get the actual node IDs from the DataFrame
                            if 'id' in nodes_df.columns:
                                for _, row in nodes_df.iterrows():
                                    node_id = int(row['id'])
                                    node_coords[node_id] = (float(row['x']), float(row['y']))
                            else:
                                for i, (idx, row) in enumerate(nodes_df.iterrows()):
                                    # Try to use the index as the node ID
                                    node_id = int(idx) if isinstance(idx, (int, float, str)) else i + 1
                                    node_coords[node_id] = (float(row['x']), float(row['y']))
                            
                            print(f"  Third method: extracted {len(node_coords)} node coordinates")
                        except Exception as e:
                            print(f"  Error using third method: {e}")
                
            # Verify we have node coordinates by checking boundary nodes
            boundaries = self._mesh.boundaries.to_dict()
            boundary_node_ids = set()
            
            for ibtype, boundary_list in boundaries.items():
                for boundary in boundary_list:
                    node_ids = boundary['node_id']
                    if isinstance(node_ids[0], tuple):  # Weir boundary
                        # First string
                        first_nodes = [int(pair[0]) for pair in node_ids]
                        boundary_node_ids.update(first_nodes)
                        # Second string
                        second_nodes = [int(pair[1]) for pair in node_ids]
                        boundary_node_ids.update(second_nodes)
                    else:  # Other boundary types
                        try:
                            node_ids = [int(node_id) for node_id in node_ids]
                            boundary_node_ids.update(node_ids)
                        except (ValueError, TypeError) as e:
                            print(f"  Warning: Could not process node_ids: {e}")
            
            # Check if boundary nodes have coordinates
            boundary_nodes_with_coords = sum(1 for node_id in boundary_node_ids if node_id in node_coords)
            print(f"  Found {len(boundary_node_ids)} unique boundary nodes")
            print(f"  {boundary_nodes_with_coords} of {len(boundary_node_ids)} boundary nodes have coordinates")
            
            # Report if still missing many boundary nodes
            if boundary_nodes_with_coords < len(boundary_node_ids) * 0.5:
                print("  WARNING: Still missing coordinates for many boundary nodes!")
                # Get a sample of missing node IDs to report
                missing_sample = list(set(boundary_node_ids) - set(node_coords))[:10]
                print(f"  Sample of missing node IDs: {missing_sample}")
                
                # As a last resort, try to directly access the original coordinate arrays
                if hasattr(self._mesh, '_nodes_x') and hasattr(self._mesh, '_nodes_y'):
                    print("  Using _nodes_x and _nodes_y attributes as last resort...")
                    x_array = self._mesh._nodes_x
                    y_array = self._mesh._nodes_y
                    
                    for node_id in boundary_node_ids:
                        # ADCIRC may use 1-based indexing, so adjust
                        array_idx = node_id - 1
                        if 0 <= array_idx < len(x_array):
                            node_coords[node_id] = (float(x_array[array_idx]), float(y_array[array_idx]))
                    
                    boundary_nodes_with_coords = sum(1 for node_id in boundary_node_ids if node_id in node_coords)
                    print(f"  After last resort: {boundary_nodes_with_coords} of {len(boundary_node_ids)} boundary nodes have coordinates")
        
        except Exception as e:
            print(f"  Error extracting node coordinates: {e}")
        
        return node_coords

    def export_boundaries(self) -> gpd.GeoDataFrame:
        """
        Export all boundaries from the mesh to a GeoDataFrame.
        
        Returns:
            GeoDataFrame containing all boundary linestrings with IBTYPE attribute
        """
        print("\nExporting mesh boundaries to GeoPackage...")
        
        # Get node coordinates
        print("Getting node coordinates...")
        node_coords = self._get_boundary_nodes_coordinates()
        print(f"Retrieved coordinates for {len(node_coords)} nodes")
        
        if not node_coords:
            print("ERROR: No node coordinates retrieved! Cannot create boundary linestrings.")
            # Return empty GeoDataFrame
            gdf = gpd.GeoDataFrame(
                columns=['IBTYPE', 'segment_id', 'geometry'],
                geometry='geometry',
                crs="EPSG:4326"
            )
            return gdf
        
        # Create a list to store boundary linestrings
        boundary_data = []
        
        # Process each boundary type
        print("Processing boundaries...")
        boundaries = self._mesh.boundaries.to_dict()
        
        for ibtype, boundary_list in boundaries.items():
            # Handle None ibtype (use -1 as a placeholder)
            if ibtype is None:
                ibtype_int = -1
                print(f"Processing boundary type None (using {ibtype_int} as placeholder) with {len(boundary_list)} segments")
            else:
                try:
                    ibtype_int = int(ibtype)
                    print(f"Processing boundary type {ibtype_int} with {len(boundary_list)} segments")
                except (ValueError, TypeError):
                    # Handle case where ibtype can't be converted to int
                    ibtype_int = -99
                    print(f"Warning: Could not convert ibtype '{ibtype}' to integer, using {ibtype_int} as placeholder")
            
            for i, boundary in enumerate(boundary_list, 1):
                node_ids = boundary['node_id']
                
                if isinstance(node_ids[0], tuple):  # Weir boundary
                    # Process first string
                    first_nodes = [int(pair[0]) for pair in node_ids]
                    first_coords = [node_coords[node_id] for node_id in first_nodes if node_id in node_coords]
                    if len(first_coords) >= 2:
                        first_line = LineString(first_coords)
                        boundary_data.append({
                            'geometry': first_line,
                            'IBTYPE': ibtype_int,
                            'segment_id': f"{ibtype_int}_{i}_1"
                        })
                        print(f"  Created linestring for weir boundary {i} first string with {len(first_coords)} nodes")
                    else:
                        print(f"  Warning: Not enough valid coordinates for weir boundary {i} first string")
                    
                    # Process second string
                    second_nodes = [int(pair[1]) for pair in node_ids]
                    second_coords = [node_coords[node_id] for node_id in second_nodes if node_id in node_coords]
                    if len(second_coords) >= 2:
                        second_line = LineString(second_coords)
                        boundary_data.append({
                            'geometry': second_line,
                            'IBTYPE': ibtype_int,
                            'segment_id': f"{ibtype_int}_{i}_2"
                        })
                        print(f"  Created linestring for weir boundary {i} second string with {len(second_coords)} nodes")
                    else:
                        print(f"  Warning: Not enough valid coordinates for weir boundary {i} second string")
                else:  # Other boundary types
                    try:
                        node_ids = [int(node_id) for node_id in node_ids]
                        coords = [node_coords[node_id] for node_id in node_ids if node_id in node_coords]
                        
                        if len(coords) >= 2:
                            line = LineString(coords)
                            boundary_data.append({
                                'geometry': line,
                                'IBTYPE': ibtype_int,
                                'segment_id': f"{ibtype_int}_{i}"
                            })
                            print(f"  Created linestring for boundary {ibtype_int}_{i} with {len(coords)} nodes")
                        else:
                            missing_nodes = [nid for nid in node_ids if nid not in node_coords]
                            print(f"  Warning: Not enough valid coordinates for boundary {ibtype_int}_{i} (only {len(coords)} of {len(node_ids)})")
                            if missing_nodes and len(missing_nodes) < 10:
                                print(f"  Missing node IDs: {missing_nodes}")
                    except (ValueError, TypeError, KeyError) as e:
                        print(f"  Warning: Could not process segment {i} of boundary type {ibtype_int}: {e}")
        
        print(f"Created {len(boundary_data)} boundary linestrings")
        
        # Create GeoDataFrame
        print("Creating GeoDataFrame...")
        if not boundary_data:
            # Create empty GeoDataFrame with the correct schema
            gdf = gpd.GeoDataFrame(
                columns=['IBTYPE', 'segment_id', 'geometry'],
                geometry='geometry',
                crs="EPSG:4326"  # Default to WGS84, will be overridden if CRS info available
            )
        else:
            gdf = gpd.GeoDataFrame(boundary_data, geometry='geometry', crs="EPSG:4326")
        
        # Try to set CRS from mesh if available
        try:
            if hasattr(self._mesh, 'crs') and self._mesh.crs is not None:
                gdf.crs = self._mesh.crs
                print(f"Set CRS from mesh: {gdf.crs}")
            else:
                print("Using default CRS: EPSG:4326 (WGS84)")
        except Exception as e:
            print(f"Warning: Could not set CRS from mesh: {e}")
            print("Using default CRS: EPSG:4326 (WGS84)")
        
        return gdf

    def export_to_geopackage(self, output_path: str) -> None:
        """
        Export mesh boundaries to GeoPackage file.
        
        Args:
            output_path: Path to save the GeoPackage file
        """
        gdf = self.export_boundaries()
        
        # Save to GeoPackage
        print(f"\nSaving to GeoPackage: {output_path}")
        if gdf.empty:
            print("Warning: No boundary segments to export. GeoPackage will be empty.")
        
        gdf.to_file(output_path, driver="GPKG")
        print(f"Successfully exported {len(gdf)} boundary segments to: {output_path}")
        
        # Print summary
        if not gdf.empty:
            ibtype_summary = gdf.groupby('IBTYPE').size()
            print("\nBoundary type summary:")
            for ibtype, count in ibtype_summary.items():
                ibtype_label = ibtype
                if ibtype == -1:
                    ibtype_label = "None"
                elif ibtype == -99:
                    ibtype_label = "Unknown"
                print(f"  IBTYPE {ibtype_label}: {count} segments")
        else:
            print("No boundary segments to export")


def main():
    """Main function to handle command line arguments and process the mesh."""
    parser = argparse.ArgumentParser(
        description="Export ADCIRC mesh boundaries to GeoPackage (.gpkg) format"
    )
    parser.add_argument(
        "input_mesh",
        help="Path to the input mesh file"
    )
    parser.add_argument(
        '-o', '--output',
        default='mesh_boundaries.gpkg',
        help='Path to save the output GeoPackage (default: mesh_boundaries.gpkg)'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable additional debug output'
    )
    
    args = parser.parse_args()
    
    # Read the mesh file
    print("Reading input mesh...")
    mesh = AdcircMesh.open(args.input_mesh)
    print(f"Successfully read mesh from: {args.input_mesh}")
    
    # Print mesh information
    print("\nInput mesh info:")
    print(f"Number of nodes: {len(mesh.nodes)}")
    print(f"Number of elements: {len(mesh.elements.elements)}")
    
    if args.debug:
        # Print more debugging information about the mesh
        print("\nMesh attributes:", dir(mesh))
        if hasattr(mesh, 'boundaries'):
            boundary_types = mesh.boundaries.to_dict().keys()
            print(f"Boundary types in mesh: {list(boundary_types)}")
    
    # Export boundaries
    exporter = BoundaryExporter(mesh)
    exporter.export_to_geopackage(args.output)
    
    return 0


if __name__ == "__main__":
    main() 