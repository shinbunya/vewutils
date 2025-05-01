import numpy as np
import pandas as pd
from shapely.geometry import Polygon, Point, LineString
from shapely.ops import unary_union
from adcircpy import AdcircMesh
import geopandas as gpd
from typing import Union, List, Set, Optional
import argparse

class NodeSelector:
    def __init__(self, mesh: AdcircMesh):
        """Initialize the NodeSelector with an ADCIRC mesh.
        
        Args:
            mesh: An AdcircMesh object containing the mesh data
        """
        self.mesh = mesh
        self.nodes_df = mesh.nodes
        self.boundaries = mesh.boundaries.to_dict()
        self.max_node_id = len(self.nodes_df)
        
    def _validate_node_id(self, node_id: int, context: str = "") -> int:
        """Validate a single node ID and raise an informative error if invalid.
        
        Args:
            node_id: Node ID to validate
            context: Context string for error message
            
        Returns:
            Validated node ID
            
        Raises:
            ValueError: If node ID is invalid
        """
        if not isinstance(node_id, (int, np.integer)):
            raise ValueError(f"Invalid node ID type in {context}: {node_id} (type: {type(node_id)})")
        if node_id < 1:
            raise ValueError(f"Node ID must be positive in {context}: {node_id}")
        if node_id > self.max_node_id:
            raise ValueError(f"Node ID exceeds mesh size in {context}: {node_id} (max: {self.max_node_id})")
        return node_id
        
    def select_by_polygon(self, polygon_file: str, tolerance: float = 1e-6) -> Set[int]:
        """Select nodes that are inside or near a polygon from a geospatial file.
        
        Args:
            polygon_file: Path to the geospatial file containing polygons
            tolerance: Distance in meters to include nodes near polygon boundaries
            
        Returns:
            Set of node IDs (1-based) that are inside or near the polygon
        """
        print(f"Reading polygon file: {polygon_file}")
        # Read the polygon file
        gdf = gpd.read_file(polygon_file)
        if gdf.crs is None:
            raise ValueError("Polygon file must have a defined CRS")
            
        print("Converting mesh nodes to GeoDataFrame...")
        # Convert mesh nodes to GeoDataFrame with actual node IDs
        nodes_gdf = gpd.GeoDataFrame(
            self.nodes_df,
            geometry=[Point(x, y) for x, y in zip(self.nodes_df['x'], self.nodes_df['y'])],
            crs=gdf.crs
        )
        
        print("Creating buffered polygon...")
        # Buffer the polygon by tolerance
        buffered_polygon = gdf.geometry.unary_union.buffer(tolerance)
        
        print("Finding nodes inside or near the polygon...")
        # Find nodes inside or near the polygon
        selected_nodes = set()
        for node_id, row in nodes_gdf.iterrows():
            if row.geometry.within(buffered_polygon) or row.geometry.distance(buffered_polygon) <= tolerance:
                try:
                    # Node IDs in ADCIRC mesh are already 1-based
                    validated_id = self._validate_node_id(node_id, "polygon selection")
                    selected_nodes.add(validated_id)
                except ValueError as e:
                    raise ValueError(f"Invalid node ID in polygon selection: {str(e)}")
        
        print(f"Found {len(selected_nodes)} nodes inside or near the polygon")
        return selected_nodes
        
    def select_by_mesh(self, mesh_file: str, tolerance: float = 1e-6) -> Set[int]:
        """Select nodes that are inside or near the boundary of another mesh.
        
        Args:
            mesh_file: Path to the mesh file defining the boundary
            tolerance: Distance in meters to include nodes near mesh boundaries
            
        Returns:
            Set of node IDs (1-based) that are inside or near the mesh boundary
        """
        print(f"Loading boundary mesh: {mesh_file}")
        # Load the boundary mesh
        boundary_mesh = AdcircMesh.open(mesh_file)
        
        print("Generating boundary polygon...")
        # Get the boundary polygon of the mesh
        boundary_polygon = self._get_mesh_boundary_polygon(boundary_mesh)
        
        print("Creating buffered boundary...")
        # Buffer the boundary by tolerance
        buffered_boundary = boundary_polygon.buffer(tolerance)
        
        print("Finding nodes inside or near the boundary...")
        # Find nodes inside or near the boundary
        selected_nodes = set()
        for node_id, row in self.nodes_df.iterrows():
            point = Point(row['x'], row['y'])
            if point.within(buffered_boundary) or point.distance(buffered_boundary) <= tolerance:
                try:
                    # Node IDs in ADCIRC mesh are already 1-based
                    validated_id = self._validate_node_id(node_id, "mesh boundary selection")
                    selected_nodes.add(validated_id)
                except ValueError as e:
                    raise ValueError(f"Invalid node ID in mesh boundary selection: {str(e)}")
        
        print(f"Found {len(selected_nodes)} nodes inside or near the boundary")
        return selected_nodes
        
    def _get_mesh_boundary_polygon(self, mesh: AdcircMesh) -> Polygon:
        """Get the boundary polygon of a mesh by finding unshared edges.
        
        Args:
            mesh: An AdcircMesh object
            
        Returns:
            A shapely Polygon representing the mesh boundary
        """
        print("Finding unshared edges...")
        # Find unshared edges (boundary edges)
        edges = []
        elems = []
        triangles = mesh.elements.elements
        for j, element in triangles.iterrows():
            nodes = element[['node_1', 'node_2', 'node_3']].astype(int).values
            for i in range(3):
                edge = tuple(sorted([nodes[i], nodes[(i + 1) % 3]]))
                edges.append(edge)
                elems.append((j, i))

        # Count edge occurrences
        edge_counts = {}
        for edge in edges:
            edge_counts[edge] = edge_counts.get(edge, 0) + 1

        # Find unshared edges (boundary edges)
        unshared_edges = []
        for i, edge in enumerate(edges):
            if edge_counts[edge] == 1:
                unshared_edges.append(edge)
        print(f"Found {len(unshared_edges)} boundary edges")

        print("Creating boundary segments...")
        # Split node string into segments
        segments = []
        processed_edges = set()
        
        while len(processed_edges) < len(unshared_edges):
            # Find first unprocessed edge
            start_idx = next(i for i, edge in enumerate(unshared_edges) if i not in processed_edges)
            current_segment = list(unshared_edges[start_idx])
            processed_edges.add(start_idx)

            # Extend segment forward
            while True:
                found_next = False
                for i, edge in enumerate(unshared_edges):
                    if i in processed_edges:
                        continue
                    if edge[0] == current_segment[-1]:
                        current_segment.append(edge[1])
                        processed_edges.add(i)
                        found_next = True
                        break
                    elif edge[1] == current_segment[-1]:
                        current_segment.append(edge[0])
                        processed_edges.add(i)
                        found_next = True
                        break
                if not found_next:
                    break

            # Add segment if it forms a closed loop
            if len(current_segment) >= 3 and current_segment[0] == current_segment[-1]:
                segments.append(current_segment[:-1])  # Remove duplicate end node
        print(f"Found {len(segments)} boundary segments")

        print("Creating boundary polygons...")
        # Create polygons from each segment
        polygons = []
        for segment in segments:
            coords = [(mesh.x[node], mesh.y[node]) for node in segment]
            coords.append(coords[0])  # Close the polygon
            polygon = Polygon(coords)
            
            # Only add valid polygons
            if polygon.is_valid and not polygon.is_empty:
                polygons.append(polygon)
        
        if not polygons:
            raise ValueError("No valid polygons were created from the segments")
            
        print("Merging boundary polygons...")
        # Merge all polygons using unary_union
        try:
            boundary = unary_union(polygons)
            if not boundary.is_valid:
                raise ValueError("Union operation produced an invalid geometry")
        except Exception as e:
            raise ValueError(f"Failed to create valid boundary geometry: {str(e)}")
            
        return boundary
        
    def select_by_csv(self, csv_file: str) -> Set[int]:
        """Select nodes from a CSV file.
        
        Args:
            csv_file: Path to the CSV file containing node IDs
            
        Returns:
            Set of node IDs (1-based) from the CSV file
            
        Raises:
            ValueError: If CSV contains invalid node IDs
        """
        df = pd.read_csv(csv_file)
        if 'node_id' not in df.columns:
            raise ValueError("CSV file must contain a 'node_id' column")
            
        selected_nodes = set()
        for node_id in df['node_id'].values:
            try:
                validated_id = self._validate_node_id(node_id, f"CSV file {csv_file}")
                selected_nodes.add(validated_id)
            except ValueError as e:
                raise ValueError(f"Invalid node ID in CSV file {csv_file}: {str(e)}")
                
        return selected_nodes
        
    def filter_by_boundary_type(self, nodes: Set[int], boundary_type: str = 'both') -> Set[int]:
        """Filter nodes based on their position relative to VEW boundaries.
        
        Args:
            nodes: Set of node IDs to filter
            boundary_type: Type of boundary nodes to include ('channel', 'bank', or 'both')
            
        Returns:
            Filtered set of node IDs
            
        Raises:
            ValueError: If boundary contains invalid node IDs
        """
        if boundary_type not in ['channel', 'bank', 'both']:
            raise ValueError("boundary_type must be 'channel', 'bank', or 'both'")
            
        if boundary_type == 'both':
            return nodes
            
        # Get VEW boundary nodes
        vew_nodes = set()
        for boundary in self.boundaries.get('4', []):  # Type 4 is VEW boundary
            for node_pair in boundary['node_id']:
                try:
                    if boundary_type == 'channel':
                        node_id = self._validate_node_id(int(node_pair[0]), "VEW boundary channel node")
                    else:  # bank
                        node_id = self._validate_node_id(int(node_pair[1]), "VEW boundary bank node")
                    vew_nodes.add(node_id)
                except ValueError as e:
                    raise ValueError(f"Invalid node ID in VEW boundary: {str(e)}")
                    
        return nodes.intersection(vew_nodes)
        
    def combine_selections(self, nodes_list: List[Set[int]], operation: str = 'union') -> Set[int]:
        """Combine multiple sets of selected nodes using union or intersection.
        
        Args:
            nodes_list: List of sets of node IDs
            operation: Operation to perform ('union' or 'intersection')
            
        Returns:
            Combined set of node IDs
        """
        if not nodes_list:
            return set()
            
        if operation == 'union':
            return set().union(*nodes_list)
        elif operation == 'intersection':
            return set.intersection(*nodes_list)
        else:
            raise ValueError("operation must be 'union' or 'intersection'")

def main():
    """Command line interface for node selection."""
    parser = argparse.ArgumentParser(
        description='Select nodes from an ADCIRC mesh based on various conditions.'
    )
    parser.add_argument(
        'mesh_file',
        help='Path to the ADCIRC mesh file'
    )
    parser.add_argument(
        '-o', '--output',
        default='selected_nodes.csv',
        help='Output CSV file path (default: selected_nodes.csv)'
    )
    parser.add_argument(
        '-p', '--polygon',
        help='Path to geospatial file containing polygons'
    )
    parser.add_argument(
        '-m', '--boundary-mesh',
        help='Path to mesh file defining boundary'
    )
    parser.add_argument(
        '-c', '--csv',
        help='Path to CSV file containing node IDs'
    )
    parser.add_argument(
        '-t', '--tolerance',
        type=float,
        default=1e-6,
        help='Tolerance in meters for boundary proximity (default: 1e-6)'
    )
    parser.add_argument(
        '-b', '--boundary-type',
        choices=['channel', 'bank', 'both'],
        default='both',
        help='Type of boundary nodes to include (default: both)'
    )
    parser.add_argument(
        '-op', '--operation',
        choices=['union', 'intersection'],
        default='union',
        help='Operation to combine multiple selections (default: union)'
    )
    
    args = parser.parse_args()
    
    print(f"\nLoading mesh: {args.mesh_file}")
    # Load mesh
    mesh = AdcircMesh.open(args.mesh_file)
    selector = NodeSelector(mesh)
    print(f"Mesh loaded: {len(mesh.nodes)} nodes, {len(mesh.elements.elements)} elements")
    
    # Perform selections
    selected_nodes = []
    
    if args.polygon:
        print("\nSelecting nodes by polygon...")
        nodes = selector.select_by_polygon(args.polygon, args.tolerance)
        selected_nodes.append(nodes)
        
    if args.boundary_mesh:
        print("\nSelecting nodes by mesh boundary...")
        nodes = selector.select_by_mesh(args.boundary_mesh, args.tolerance)
        selected_nodes.append(nodes)
        
    if args.csv:
        print(f"\nReading nodes from CSV: {args.csv}")
        nodes = selector.select_by_csv(args.csv)
        print(f"Found {len(nodes)} nodes in CSV file")
        selected_nodes.append(nodes)
        
    # Combine selections
    if selected_nodes:
        print(f"\nCombining selections using {args.operation}...")
        final_nodes = selector.combine_selections(selected_nodes, args.operation)
        print(f"Combined selection contains {len(final_nodes)} nodes")
    else:
        final_nodes = set()
        
    # Filter by boundary type
    if args.boundary_type != 'both':
        print(f"\nFiltering nodes by boundary type: {args.boundary_type}")
        original_count = len(final_nodes)
        final_nodes = selector.filter_by_boundary_type(final_nodes, args.boundary_type)
        print(f"Filtered from {original_count} to {len(final_nodes)} nodes")
    
    # Write results
    print(f"\nWriting {len(final_nodes)} nodes to: {args.output}")
    df = pd.DataFrame({'node_id': sorted(final_nodes)})
    df.to_csv(args.output, index=False)
    print("Done!")

if __name__ == '__main__':
    main() 