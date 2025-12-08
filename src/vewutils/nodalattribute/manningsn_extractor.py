# %%
import os
import rasterio
from rasterio.mask import mask
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely import Point, Polygon
from adcircpy import AdcircMesh
from adcircpy.mesh.fort13 import NodalAttributes
import argparse
import math

class ManningsnExtractor:
    """Class for extracting average Manning's n values from land cover data for ADCIRC mesh nodes."""
    
    def __init__(self, mesh_file, landcover_file, class_to_mn_file, output_file, 
                 selected_nodes_file=None, fort13_file=None, output_format='fort13'):
        """Initialize the ManningsnExtractor.
        
        Args:
            mesh_file (str): Path to the mesh file
            landcover_file (str): Path to the land cover raster file
            class_to_mn_file (str): Path to the class to Manning's n mapping file
            output_file (str): Path to the output file
            selected_nodes_file (str, optional): Path to CSV file containing selected nodes
            fort13_file (str, optional): Path to input fort.13 file
            output_format (str): Output format ('fort13' or 'csv')
        """
        self.mesh_file = mesh_file
        self.landcover_file = landcover_file
        self.class_to_mn_file = class_to_mn_file
        self.output_file = output_file
        self.selected_nodes_file = selected_nodes_file
        self.fort13_file = fort13_file
        self.output_format = output_format
        
        # Load mesh and class mapping
        self.mesh = AdcircMesh.open(mesh_file)
        self.class_to_mn = pd.read_csv(class_to_mn_file)
        self.map_class = np.array(self.class_to_mn.classid)
        self.map_mn = np.array(self.class_to_mn.mn)
        
        # Get open water values
        self.class_openwater = self.class_to_mn[self.class_to_mn.classification == "Open Water"].classid.values[0]
        self.mn_openwater = self.class_to_mn[self.class_to_mn.classification == "Open Water"].mn.values[0]
        
        # Load selected nodes if provided
        self.selected_nodes = None
        if selected_nodes_file:
            self.selected_nodes = pd.read_csv(selected_nodes_file)['node_id'].values
            
            # Check if any nodes were specified
            if len(self.selected_nodes) == 0:
                raise ValueError(
                    f"Error: No nodes found in the selected nodes file: {selected_nodes_file}\n"
                    "Please ensure the file contains valid node IDs in a 'node_id' column."
                )
            
            # Validate node IDs
            max_node_id = self.mesh.nodes.shape[0]
            invalid_nodes = set(self.selected_nodes[self.selected_nodes < 1]) | set(self.selected_nodes[self.selected_nodes > max_node_id])
            if len(invalid_nodes) > 0:
                raise ValueError(
                    f"Error: The following node IDs in {selected_nodes_file} exceed the total number of nodes ({max_node_id}) in the mesh:\n"
                    f"{invalid_nodes}\n"
                    "Please ensure all node IDs are valid."
                )
            
        # Load fort.13 file if provided
        self.fort13 = None
        if fort13_file:
            self.fort13 = NodalAttributes(self.mesh)
            self.fort13.import_fort13(fort13_file)

    def _build_neighbors_table(self, nodes_to_process):
        """Build neighbors table for all nodes in nodes_to_process using element connectivity.
        
        Args:
            nodes_to_process (list): List of 0-based node indices to process
            
        Returns:
            dict: Dictionary mapping node_idx -> list of neighbor indices (1-based, ordered counter-clockwise)
        """
        # Convert nodes_to_process to 1-based for element matching
        nodes_1based_set = set(node_idx + 1 for node_idx in nodes_to_process)
        
        # Build node-to-elements table for relevant nodes
        node_to_elements = {node_idx: [] for node_idx in nodes_to_process}
        
        # Single pass through all elements to build node-to-element mapping
        triangles_array = self.mesh.triangles.values
        
        for elem_idx, element_nodes in enumerate(triangles_array):
            # element_nodes are already 1-based as stored in ADCIRC format
            
            # Check if any nodes in this element are in our nodes_to_process
            for node_id_1based in element_nodes:
                if node_id_1based in nodes_1based_set:
                    node_idx_0based = node_id_1based - 1
                    node_to_elements[node_idx_0based].append((elem_idx, element_nodes))
        
        # Build ordered neighbors table using element connectivity
        neighbors_table = {}
        
        for node_idx_0based in node_to_elements:
            if len(node_to_elements[node_idx_0based]) < 3:
                # Not enough elements to form a proper polygon
                neighbors = set()
                for elem_idx, element_nodes in node_to_elements[node_idx_0based]:
                    for neighbor_id_1based in element_nodes:
                        if neighbor_id_1based != (node_idx_0based + 1):
                            neighbors.add(neighbor_id_1based)
                neighbors_table[node_idx_0based] = list(neighbors)
                continue
            
            # Use element connectivity to order neighbors
            ordered_neighbors = self._order_neighbors_by_elements(
                node_idx_0based + 1,  # Convert to 1-based
                node_to_elements[node_idx_0based]
            )
            neighbors_table[node_idx_0based] = ordered_neighbors
        
        return neighbors_table

    def _order_neighbors_by_elements(self, center_node_1based, node_elements):
        """Order neighbor nodes by traversing around the center node via element connectivity.
        Takes advantage of counter-clockwise node ordering in elements.
        
        Args:
            center_node_1based (int): The central node ID (1-based)
            node_elements (list): List of (element_index, element_nodes) tuples containing the center node
            
        Returns:
            list: Ordered list of neighbor node IDs (1-based)
        """
        if len(node_elements) == 0:
            return []
        
        # Build element connectivity map: for each element, what are the two neighbors of center_node
        # and which neighbor connects to which adjacent element
        element_neighbors = {}
        
        for elem_idx, elem_nodes in node_elements:
            center_pos = list(elem_nodes).index(center_node_1based)
            # Since nodes are counter-clockwise, next neighbor is at (center_pos + 1) % 3
            # and previous neighbor is at (center_pos + 2) % 3
            next_neighbor = elem_nodes[(center_pos + 1) % 3]
            prev_neighbor = elem_nodes[(center_pos + 2) % 3]
            element_neighbors[elem_idx] = (prev_neighbor, next_neighbor)  # (incoming, outgoing)
        
        # Start with any element and traverse around the node
        start_elem_idx = node_elements[0][0]
        ordered_neighbors = []
        current_elem = start_elem_idx
        visited_elements = set()
        
        while current_elem is not None and current_elem not in visited_elements:
            visited_elements.add(current_elem)
            prev_neighbor, next_neighbor = element_neighbors[current_elem]
            
            # Add the outgoing neighbor from this element
            ordered_neighbors.append(next_neighbor)
            
            # Find the next element that has this neighbor as its incoming edge
            next_elem = None
            for elem_idx, (incoming, outgoing) in element_neighbors.items():
                if elem_idx not in visited_elements and incoming == next_neighbor:
                    next_elem = elem_idx
                    break
            
            current_elem = next_elem
        
        # Handle case where we couldn't complete the loop (boundary node)
        if len(visited_elements) < len(node_elements):
            # Add any remaining neighbors that weren't included
            all_neighbors = set()
            for elem_idx, (prev_neighbor, next_neighbor) in element_neighbors.items():
                all_neighbors.add(prev_neighbor)
                all_neighbors.add(next_neighbor)
            
            for neighbor in all_neighbors:
                if neighbor not in ordered_neighbors:
                    ordered_neighbors.append(neighbor)
        
        return ordered_neighbors

    def _create_nodal_polygons(self):
        """Create polygons representing the extent of each node's neighboring elements.
        
        Returns:
            gpd.GeoDataFrame: GeoDataFrame containing node IDs and their corresponding polygons
        """
        node_ids = []
        node_polys = []
        
        print('Creating nodal extent polygons...')
        
        if self.selected_nodes is not None:
            # Only process selected nodes
            nodes_to_process = [node_id - 1 for node_id in self.selected_nodes]  # Convert to 0-based index
            total_nodes = len(nodes_to_process)
        else:
            # Process all nodes
            nodes_to_process = list(range(self.mesh.nodes.shape[0]))
            total_nodes = self.mesh.nodes.shape[0]
        
        # Build neighbors table for all nodes to process in a single pass
        print('  Building neighbors table...')
        neighbors_table = self._build_neighbors_table(nodes_to_process)
        print('  Neighbors table built.')
        
        for i, node_idx in enumerate(nodes_to_process):
            if i % 1000 == 0:
                print(f'  Processing node {i+1} of {total_nodes}')
                
            # Get pre-computed neighbors for this node
            nei = list(neighbors_table[node_idx])
            nnei = len(nei)
            xs = list(self.mesh.x.loc[nei])
            ys = list(self.mesh.y.loc[nei])
            if len(xs) >= 3:
                xs.append(xs[0])
                ys.append(ys[0])
                node_ids.append(node_idx)
                node_polys.append(Polygon(zip(xs, ys)))
                
        print('  Done.')
        
        gdf = gpd.GeoDataFrame({'node_id': node_ids, 'geometry': node_polys}, crs='EPSG:4326')

        return gdf

    def _extract_mannings_values(self, dfpoly):
        """Extract average Manning's n values for each node.
        
        Args:
            dfpoly (gpd.GeoDataFrame): GeoDataFrame containing node polygons
            
        Returns:
            list: List of average Manning's n values for each node
        """
        nn = self.mesh.nodes.shape[0]
        mn_avgs = [self.mn_openwater] * nn
        
        # Open the tif file to get its CRS
        with rasterio.open(self.landcover_file) as src:
            raster_crs = src.crs
        
        # Convert the CRS of dfpoly to match the raster CRS
        dfpoly_rc = dfpoly.to_crs(raster_crs)
        
        # Open the tif file and process each polygon
        with rasterio.open(self.landcover_file) as src:
            print('Extracting average mannings n values...')
            for idx, row in dfpoly_rc.iterrows():
                if idx % 1000 == 0:
                    print(f'  Processing node {idx + 1} of {len(dfpoly_rc)} total nodes...')
                node_id = row['node_id']
                
                # Skip nodes not in selected_nodes if provided
                if self.selected_nodes is not None and (node_id + 1) not in self.selected_nodes:
                    continue
                    
                polygon = [row['geometry'].__geo_interface__]
                
                try:
                    out_image, out_transform = mask(src, polygon, crop=True)
                except:
                    continue
                
                lcs = out_image.flatten()
                lcs = lcs[lcs != 0]
                
                if len(lcs) == 0:
                    continue
                else:
                    mn_sum = 0.0
                    for lc in lcs:
                        mn_sum += self.map_mn[self.map_class == lc][0]
                    mn_avg = mn_sum / len(lcs)
                    mn_avgs[node_id] = mn_avg
                    
            print('  Done.')
            
        return mn_avgs

    def _write_output(self, mn_avgs):
        """Write the extracted Manning's n values to the output file.

        Args:
            mn_avgs (list): List of average Manning's n values for each node.
        """
        # Initialize final_mn_avgs with open water values
        final_mn_avgs = np.full(len(self.mesh.nodes), self.mn_openwater)
        
        # If fort.13 file is provided, get existing values
        if self.fort13 is not None:
            try:
                existing_mn = self.fort13.get_attribute('mannings_n_at_sea_floor')['values']
                if existing_mn.ndim > 1:
                    existing_mn = existing_mn.flatten()
                final_mn_avgs = existing_mn.copy()
            except (KeyError, AttributeError):
                pass  # If attribute doesn't exist, use open water values
        
        # Convert mn_avgs to numpy array for array indexing
        mn_avgs = np.array(mn_avgs)
        
        # Update values for selected nodes or all nodes
        if self.selected_nodes is not None:
            # Convert selected nodes to 0-based indexing
            selected_nodes_0based = np.array(self.selected_nodes) - 1
            # Update values - mn_avgs already contains values at correct indices
            final_mn_avgs[selected_nodes_0based] = mn_avgs[selected_nodes_0based]
        else:
            # If no nodes are selected, update all nodes
            final_mn_avgs = mn_avgs.copy()
            
        if self.output_format == 'fort13':
            if self.fort13 is not None:
                # Update the fort.13 file
                # Reshape the array to 2D if needed
                final_mn_avgs_2d = final_mn_avgs.reshape(-1, 1)
                
                # Reset defaults and non_default_indexes so they get recomputed from new values
                if 'mannings_n_at_sea_floor' in self.fort13._attributes:
                    self.fort13._attributes['mannings_n_at_sea_floor']['defaults'] = None
                    self.fort13._attributes['mannings_n_at_sea_floor'].pop('non_default_indexes', None)
                
                self.fort13.set_attribute('mannings_n_at_sea_floor', final_mn_avgs_2d)
                self.fort13.write(self.output_file, overwrite=True)
            else:
                # Create new fort.13 file
                fort13 = NodalAttributes(self.mesh)
                fort13.add_attribute('mannings_n_at_sea_floor', 'unitless')
                # Reshape the array to 2D if needed
                final_mn_avgs_2d = final_mn_avgs.reshape(-1, 1)
                fort13.set_attribute('mannings_n_at_sea_floor', final_mn_avgs_2d)
                fort13.write(self.output_file, overwrite=True)
        else:  # csv format
            # Create DataFrame for output
            dfmn = pd.DataFrame({'node_id': np.arange(len(final_mn_avgs)), 'mn_avg': final_mn_avgs})
            # Only include non-open water values
            dfmn = dfmn[dfmn.mn_avg != self.mn_openwater]
            
            with open(self.output_file, "w") as f:
                f.write("{:.3f}\n".format(self.mn_openwater))
                f.write("{:d}\n".format(dfmn.shape[0]))
                for i in dfmn.index:
                    f.write("{:d}    {:.3f}\n".format(dfmn.node_id[i] + 1, dfmn.mn_avg[i]))

    def extract(self):
        """Main method to execute the Manning's n extraction process."""
        # Create nodal polygons
        dfpoly = self._create_nodal_polygons()
        
        # Extract Manning's n values
        mn_avgs = self._extract_mannings_values(dfpoly)
        
        # Write output
        self._write_output(mn_avgs)

def get_parser():
    parser = argparse.ArgumentParser(add_help=False, description="Extract the average Manning's n values of the neighboring elements of each node.")
    parser.add_argument(
        'mesh',
        help='Path to the mesh file'
    )
    parser.add_argument(
        'landcover',
        help='Path to the land cover raster file'
    )
    parser.add_argument(
        'class_to_mn',
        help="Path to the class to Manning's n mapping file"
    )
    parser.add_argument(
        'output',
        help='Path to the output file'
    )
    parser.add_argument(
        '--selected-nodes',
        help='Path to CSV file containing selected nodes',
        default=None
    )
    parser.add_argument(
        '--fort13',
        help='Path to input fort.13 file',
        default=None
    )
    parser.add_argument(
        '--output-format',
        help="Output format ('fort13' or 'csv')",
        default='fort13'
    )
    return parser

def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    extractor = ManningsnExtractor(
        args.mesh,
        args.landcover,
        args.class_to_mn,
        args.output,
        selected_nodes_file=args.selected_nodes,
        fort13_file=args.fort13,
        output_format=args.output_format
    )
    extractor.extract()

if __name__ == '__main__':
    main()


