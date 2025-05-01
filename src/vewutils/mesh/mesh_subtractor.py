import numpy as np
import pandas as pd
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
from adcircpy import AdcircMesh
from shapely.vectorized import contains

class MeshSubtractor:
    def __init__(self):
        pass

    def _find_edges(self, mesh):
        """Find edges and count how many elements each edge belongs to."""
        edges = []
        elems = []
        # Access triangular elements through elements.elements
        triangles = mesh.elements.elements
        for j, element in triangles.iterrows():
            # Elements table has columns: ['id', 'node_1', 'node_2', 'node_3']
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
        unshared_elems = []
        for i, edge in enumerate(edges):
            if edge_counts[edge] == 1:
                unshared_edges.append(edge)
                unshared_elems.append(elems[i])

        return unshared_edges, unshared_elems

    def _split_node_string(self, mesh, unshared_edges, unshared_elems):
        """Split node string into segments based on domain boundaries."""
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

        return segments

    def _generate_boundary_polygons(self, mesh, segments):
        """Generate a polygon from the mesh boundary segments.
        
        Args:
            mesh: AdcircMesh object containing node coordinates
            segments: List of node segments defining boundary rings
            
        Returns:
            MultiPolygon or Polygon: A shapely geometry representing the union of all boundary polygons
            
        Raises:
            ValueError: If no valid polygons could be created from the segments
        """
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
            
        # Merge all polygons using unary_union
        try:
            boundary = unary_union(polygons)
            if not boundary.is_valid:
                raise ValueError("Union operation produced an invalid geometry")
        except Exception as e:
            raise ValueError(f"Failed to create valid boundary geometry: {str(e)}")
            
        return boundary

    def subtract(self, mesh_a, mesh_b):
        """Remove elements from mesh_a that are inside mesh_b's domain boundary."""
        print("Starting mesh subtraction process...")
        
        # Find boundary edges of mesh_b
        print("Finding boundary edges of mesh B...")
        unshared_edges_b, unshared_elems_b = self._find_edges(mesh_b)
        print(f"Found {len(unshared_edges_b)} boundary edges")
        
        # Split into boundary segments
        print("Splitting boundary into segments...")
        segments_b = self._split_node_string(mesh_b, unshared_edges_b, unshared_elems_b)
        print(f"Found {len(segments_b)} boundary segments")
        
        # Generate boundary polygon for mesh_b
        print("Generating boundary polygon for mesh B...")
        boundary_b = self._generate_boundary_polygons(mesh_b, segments_b)
        print("Boundary polygon generated")
        
        # Calculate centroids for all elements in mesh_a
        print("Calculating element centroids...")
        triangles_a = mesh_a.elements.elements
        nodes_df = mesh_a.nodes
        
        # Convert to numpy for faster processing
        node_cols = ['node_1', 'node_2', 'node_3']
        element_nodes = triangles_a[node_cols].values - 1  # Convert to 0-based indexing
        
        # Calculate centroids for all elements
        x_coords = nodes_df['x'].values
        y_coords = nodes_df['y'].values
        
        # Get x,y coordinates for each node in each element
        element_x = x_coords[element_nodes]
        element_y = y_coords[element_nodes]
        
        # Calculate centroids (mean of triangle vertices)
        centroid_x = np.mean(element_x, axis=1)
        centroid_y = np.mean(element_y, axis=1)
        
        # Test containment of element centroids
        print("Testing containment of element centroids...")
        elements_inside = contains(boundary_b, centroid_x, centroid_y)
        elements_to_keep = ~elements_inside  # Keep elements whose centroids are outside
        keep_indices = np.where(elements_to_keep)[0]
        
        print(f"Found {len(keep_indices)} elements to keep out of {len(triangles_a)}")
        
        # Create new mesh with kept elements
        print("Creating new mesh with remaining elements...")
        elements_new = triangles_a.iloc[keep_indices].copy()
        
        # Find unique nodes in kept elements
        print("Renumbering nodes...")
        unique_nodes = np.unique(elements_new[node_cols].values)  # These are 1-based indices
        node_map = {old: new for new, old in enumerate(unique_nodes, start=1)}  # Keep 1-based
        
        # Create new node arrays using 1-based indices directly
        new_x = mesh_a.nodes.loc[unique_nodes, 'x']
        new_y = mesh_a.nodes.loc[unique_nodes, 'y']
        new_z = mesh_a.nodes.loc[unique_nodes, 'value_1']
        
        # Renumber elements (already 1-based)
        for col in node_cols:
            elements_new[col] = elements_new[col].map(node_map)
        
        # Reset element IDs to be 1-based
        elements_new.index = range(1, len(elements_new) + 1)
        
        # Create nodes DataFrame with 1-based index for ADCIRC convention
        nodes_new = pd.DataFrame(
            index=range(1, len(unique_nodes) + 1),  # 1-based indexing
            data={
                'x': new_x.values,
                'y': new_y.values,
                'value_1': new_z.values
            }
        )
        
        # Renumbering node numbers in boundaries
        boundaries_a = mesh_a.boundaries.to_dict()
        unique_nodes_set = set(unique_nodes)  # Convert to set for faster lookup
        
        # Check for boundary nodes that would be removed
        boundaries_new = {}
        for ibtype in boundaries_a.keys():
            boundaries_ibtype = boundaries_a[ibtype]
            boundaries_new_ibtype = []
            
            for i in range(len(boundaries_ibtype)):
                to_be_removed = False
                node_ids = boundaries_ibtype[i]['node_id']
                if ibtype is None: # open boundary
                    node_ids1 = node_ids
                elif ibtype.endswith('4'): # weir boundary
                    node_ids1 = [node_id[0] for node_id in node_ids]
                    node_ids2 = [node_id[1] for node_id in node_ids]
                else: # other boundaries
                    node_ids1 = node_ids
                    
                for nid in node_ids1:
                    if nid not in unique_nodes_set:
                        to_be_removed = True
                        break

                if not to_be_removed:
                    boundary_new = boundaries_ibtype[i].copy()
                    if ibtype is None: # open boundary
                        new_node_ids = [str(node_map[nid]) for nid in node_ids]
                    elif ibtype.endswith('4'): # weir boundary
                        new_node_ids = [(str(node_map[node_ids1[j]]), str(node_map[node_ids2[j]])) for j in range(len(node_ids))]
                    else: # other boundaries
                        new_node_ids = [str(node_map[nid]) for nid in node_ids]
                    boundary_new['node_id'] = new_node_ids
                    
                    boundaries_new_ibtype.append(boundary_new)
            
            if len(boundaries_new_ibtype) > 0:
                boundaries_new[ibtype] = boundaries_new_ibtype
                
        # Create new mesh
        print("Creating final mesh object...")
        new_mesh = AdcircMesh(nodes=nodes_new, elements=elements_new, boundaries=boundaries_new)
        
        print("Mesh subtraction complete!")
        print(f"Original mesh: {len(triangles_a)} elements, {len(nodes_df)} nodes")
        print(f"Result mesh: {len(elements_new)} elements, {len(nodes_new)} nodes")
        
        return new_mesh

def main():
    """Command line interface for mesh subtraction."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Remove elements from mesh A that are inside mesh B\'s domain boundary.'
    )
    parser.add_argument(
        'mesh_a', 
        help='Path to first mesh file (mesh to subtract from)'
    )
    parser.add_argument(
        'mesh_b',
        help='Path to second mesh file (mesh defining subtraction boundary)'
    )
    parser.add_argument(
        '-o', '--output',
        default='subtracted_mesh.grd',
        help='Output mesh file path (default: subtracted_mesh.grd)'
    )
    parser.add_argument(
        '-d', '--description',
        default='subtracted',
        help='Description for the subtracted mesh (default: subtracted)'
    )
    
    args = parser.parse_args()
    
    # Load meshes
    mesh_a = AdcircMesh.open(args.mesh_a, crs=None)
    mesh_b = AdcircMesh.open(args.mesh_b, crs=None)
    
    # Perform subtraction
    subtractor = MeshSubtractor()
    result_mesh = subtractor.subtract(mesh_a, mesh_b)

    # Write result
    result_mesh.description = args.description
    result_mesh.write(args.output, overwrite=True)
    print(f"Subtracted mesh written to {args.output}")

if __name__ == '__main__':
    main()
