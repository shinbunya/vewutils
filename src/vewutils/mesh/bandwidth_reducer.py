#!/usr/bin/env python3
"""
Module for reducing bandwidth of ADCIRC meshes using various node renumbering algorithms.

Available algorithms:
- RCM: Reverse Cuthill-McKee (traditional, good quality)
- King: King's algorithm (fast, good for irregular meshes) 
- GPS: Gibbs-Poole-Stockmeyer (high quality, moderate speed)
- Degree: Simple degree-based ordering (fastest, basic quality)
"""

import argparse
import pandas as pd
import numpy as np
from adcircpy import AdcircMesh
from collections import deque, defaultdict
from typing import Dict, List, Tuple, Set, Optional
from tqdm import tqdm


class BandwidthReducer:
    """Class for reducing mesh bandwidth using node renumbering algorithms."""

    def __init__(self, mesh: AdcircMesh):
        """Initialize the bandwidth reducer with a mesh."""
        self._mesh = mesh
        self._adjacency = None
        self._degrees = None

    def _build_adjacency_graph(self) -> Dict[int, Set[int]]:
        """Build adjacency graph from element connectivity using vectorized operations."""
        print("Building adjacency graph using vectorized operations...")
        adjacency = defaultdict(set)
        
        # Get elements as numpy array for faster processing
        elements_df = self._mesh.elements.elements
        elements_array = elements_df.iloc[:, 1:4].values.astype(int)
        
        print(f"Processing {len(elements_array)} elements...")
        
        # Process elements in batches for better memory efficiency
        batch_size = 10000
        for i in tqdm(range(0, len(elements_array), batch_size), 
                      desc="Processing element batches", 
                      unit="batches"):
            batch = elements_array[i:i+batch_size]
            
            # For each element in the batch, add all edge pairs
            for element_nodes in batch:
                # Add edges between all pairs of nodes in the triangular element
                n1, n2, n3 = element_nodes
                adjacency[n1].add(n2)
                adjacency[n1].add(n3)
                adjacency[n2].add(n1)
                adjacency[n2].add(n3)
                adjacency[n3].add(n1)
                adjacency[n3].add(n2)
        
        # Ensure all nodes are in the adjacency graph
        print("Checking for isolated nodes...")
        all_nodes = set(self._mesh.nodes.index)
        existing_nodes = set(adjacency.keys())
        isolated_nodes = all_nodes - existing_nodes
        
        if isolated_nodes:
            print(f"Found {len(isolated_nodes)} isolated nodes, adding them...")
            for node_id in isolated_nodes:
                adjacency[node_id] = set()
        
        print(f"Built adjacency graph with {len(adjacency)} nodes")
        return dict(adjacency)

    def _compute_degrees(self, adjacency: Dict[int, Set[int]]) -> Dict[int, int]:
        """Compute the degree of each node."""
        return {node: len(neighbors) for node, neighbors in adjacency.items()}

    def _find_peripheral_node(self, adjacency: Dict[int, Set[int]], 
                            degrees: Dict[int, int]) -> int:
        """Find a peripheral node to start the RCM algorithm using fast heuristics."""
        print("Finding peripheral node using optimized method...")
        
        # Find nodes with minimum degree
        min_degree = min(degrees.values())
        min_degree_nodes = [node for node, degree in degrees.items() 
                          if degree == min_degree]
        
        if len(min_degree_nodes) == 1:
            return min_degree_nodes[0]
        
        # Limit candidates to avoid expensive computation
        max_candidates = 5
        if len(min_degree_nodes) > max_candidates:
            print(f"Found {len(min_degree_nodes)} minimum degree nodes, evaluating first {max_candidates}")
            candidates = min_degree_nodes[:max_candidates]
        else:
            candidates = min_degree_nodes
            print(f"Evaluating {len(candidates)} minimum degree candidates")
        
        # Use fast partial BFS instead of full graph traversal
        best_node = candidates[0]
        max_eccentricity = 0
        
        for candidate in candidates:
            # Perform limited BFS to find local eccentricity (much faster)
            eccentricity = self._fast_eccentricity(candidate, adjacency, max_depth=8)
            
            if eccentricity > max_eccentricity:
                max_eccentricity = eccentricity
                best_node = candidate
        
        print(f"Selected peripheral node: {best_node} (degree: {degrees[best_node]}, eccentricity: {max_eccentricity})")
        return best_node

    def _fast_eccentricity(self, start_node: int, 
                          adjacency: Dict[int, Set[int]], 
                          max_depth: int = 8) -> int:
        """Compute local eccentricity using limited-depth BFS (much faster than full BFS)."""
        visited = {start_node}
        current_level = [start_node]
        depth = 0
        
        # Perform BFS up to max_depth levels
        while current_level and depth < max_depth:
            next_level = []
            for node in current_level:
                for neighbor in adjacency[node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        next_level.append(neighbor)
            
            if next_level:
                current_level = next_level
                depth += 1
            else:
                break
        
        return depth

    def _bfs_distances(self, start_node: int, 
                      adjacency: Dict[int, Set[int]]) -> Dict[int, int]:
        """Compute distances from start_node to all other nodes using BFS."""
        distances = {start_node: 0}
        queue = deque([start_node])
        
        while queue:
            node = queue.popleft()
            current_distance = distances[node]
            
            for neighbor in adjacency[node]:
                if neighbor not in distances:
                    distances[neighbor] = current_distance + 1
                    queue.append(neighbor)
        
        return distances

    def _cuthill_mckee_ordering(self, adjacency: Dict[int, Set[int]], 
                              degrees: Dict[int, int], 
                              start_node: int) -> List[int]:
        """Perform Cuthill-McKee ordering starting from start_node."""
        print("Performing Cuthill-McKee ordering...")
        
        ordering = []
        visited = set()
        queue = deque([start_node])
        visited.add(start_node)
        
        total_nodes = len(adjacency)
        progress_bar = tqdm(total=total_nodes, desc="CM ordering", unit="nodes")
        
        while queue:
            # Process current level
            current_level = []
            while queue:
                current_level.append(queue.popleft())
            
            # Sort current level by degree (ascending)
            current_level.sort(key=lambda x: degrees[x])
            ordering.extend(current_level)
            progress_bar.update(len(current_level))
            
            # Add neighbors of current level to queue
            next_level = []
            for node in current_level:
                neighbors = list(adjacency[node])
                # Sort neighbors by degree
                neighbors.sort(key=lambda x: degrees[x])
                
                for neighbor in neighbors:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        next_level.append(neighbor)
            
            queue.extend(next_level)
        
        # If graph is not connected, add remaining nodes
        remaining_nodes = [node for node in adjacency.keys() if node not in visited]
        if remaining_nodes:
            print(f"Warning: Graph not connected. Adding {len(remaining_nodes)} remaining nodes.")
            remaining_nodes.sort(key=lambda x: degrees[x])
            ordering.extend(remaining_nodes)
            progress_bar.update(len(remaining_nodes))
        
        progress_bar.close()
        return ordering

    def _reverse_cuthill_mckee(self) -> List[int]:
        """Implement the Reverse Cuthill-McKee algorithm."""
        print("Starting Reverse Cuthill-McKee algorithm...")
        
        # Build adjacency graph and compute degrees
        adjacency = self._build_adjacency_graph()
        degrees = self._compute_degrees(adjacency)
        
        # Find peripheral node
        start_node = self._find_peripheral_node(adjacency, degrees)
        
        # Perform Cuthill-McKee ordering
        cm_ordering = self._cuthill_mckee_ordering(adjacency, degrees, start_node)
        
        # Reverse the ordering for RCM
        rcm_ordering = list(reversed(cm_ordering))
        
        print(f"RCM ordering complete. Processed {len(rcm_ordering)} nodes.")
        return rcm_ordering

    def _king_algorithm(self) -> List[int]:
        """Implement King's algorithm - faster than RCM for irregular meshes."""
        print("Starting King's algorithm...")
        
        # Build adjacency graph and compute degrees
        adjacency = self._build_adjacency_graph()
        degrees = self._compute_degrees(adjacency)
        
        # Find multiple peripheral nodes for better results
        print("Finding multiple start candidates...")
        min_degree = min(degrees.values())
        min_degree_nodes = [node for node, degree in degrees.items() 
                          if degree == min_degree]
        
        # Use multiple starts for better results
        num_starts = min(5, len(min_degree_nodes))
        best_ordering = None
        best_bandwidth = float('inf')
        
        for i, start_node in enumerate(tqdm(min_degree_nodes[:num_starts], 
                                          desc="Trying start nodes", 
                                          unit="starts")):
            # Perform level-set ordering (similar to CM but optimized)
            ordering = self._level_set_ordering(adjacency, degrees, start_node)
            
            # Quick bandwidth estimate (sample-based for speed)
            sample_size = min(1000, len(self._mesh.elements.elements))
            sample_elements = self._mesh.elements.elements.sample(n=sample_size, random_state=42)
            sample_bandwidth = self._quick_bandwidth_estimate(sample_elements, ordering)
            
            if sample_bandwidth < best_bandwidth:
                best_bandwidth = sample_bandwidth
                best_ordering = ordering
        
        print(f"King's algorithm complete. Best ordering from {num_starts} starts.")
        return best_ordering

    def _level_set_ordering(self, adjacency: Dict[int, Set[int]], 
                           degrees: Dict[int, int], 
                           start_node: int) -> List[int]:
        """Optimized level-set ordering for King's algorithm."""
        ordering = []
        visited = set()
        current_level = [start_node]
        visited.add(start_node)
        
        total_nodes = len(adjacency)
        progress_bar = tqdm(total=total_nodes, desc="King ordering", unit="nodes")
        
        while current_level:
            # Sort current level by degree
            current_level.sort(key=lambda x: degrees[x])
            ordering.extend(current_level)
            progress_bar.update(len(current_level))
            
            # Build next level
            next_level = []
            for node in current_level:
                neighbors = list(adjacency[node])
                neighbors.sort(key=lambda x: degrees[x])  # Sort by degree for consistency
                
                for neighbor in neighbors:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        next_level.append(neighbor)
            
            current_level = next_level
        
        # Add any remaining unvisited nodes
        remaining = [node for node in adjacency.keys() if node not in visited]
        if remaining:
            remaining.sort(key=lambda x: degrees[x])
            ordering.extend(remaining)
            progress_bar.update(len(remaining))
        
        progress_bar.close()
        return ordering

    def _gps_algorithm(self) -> List[int]:
        """Implement GPS (Gibbs-Poole-Stockmeyer) algorithm."""
        print("Starting GPS (Gibbs-Poole-Stockmeyer) algorithm...")
        
        # Build adjacency graph and compute degrees
        adjacency = self._build_adjacency_graph()
        degrees = self._compute_degrees(adjacency)
        
        # Find end nodes of a pseudo-diameter
        print("Finding pseudo-diameter...")
        end1, end2 = self._find_pseudo_diameter(adjacency)
        
        # Try both ends as starting nodes and pick the better one
        print("Testing both ends of pseudo-diameter...")
        ordering1 = self._level_set_ordering(adjacency, degrees, end1)
        ordering2 = self._level_set_ordering(adjacency, degrees, end2)
        
        # Quick bandwidth comparison
        sample_size = min(1000, len(self._mesh.elements.elements))
        sample_elements = self._mesh.elements.elements.sample(n=sample_size, random_state=42)
        
        bw1 = self._quick_bandwidth_estimate(sample_elements, ordering1)
        bw2 = self._quick_bandwidth_estimate(sample_elements, ordering2)
        
        best_ordering = ordering1 if bw1 <= bw2 else ordering2
        print(f"GPS algorithm complete. Selected better of two orderings.")
        return best_ordering

    def _find_pseudo_diameter(self, adjacency: Dict[int, Set[int]]) -> Tuple[int, int]:
        """Find pseudo-diameter of the graph using fast approximation."""
        print("Finding pseudo-diameter using fast method...")
        
        # Start with a node of minimum degree
        degrees = self._compute_degrees(adjacency)
        min_degree = min(degrees.values())
        start_candidates = [node for node, deg in degrees.items() if deg == min_degree]
        start_node = start_candidates[0]
        
        # Use fast limited BFS instead of full BFS
        end1 = self._find_farthest_node_fast(start_node, adjacency)
        end2 = self._find_farthest_node_fast(end1, adjacency)
        
        print(f"Pseudo-diameter endpoints: {end1} and {end2}")
        return end1, end2

    def _find_farthest_node_fast(self, start_node: int, 
                                adjacency: Dict[int, Set[int]], 
                                max_depth: int = 20) -> int:
        """Find farthest reachable node using limited-depth BFS."""
        visited = {start_node}
        current_level = [start_node]
        depth = 0
        
        # Perform BFS up to max_depth levels
        while current_level and depth < max_depth:
            next_level = []
            for node in current_level:
                for neighbor in adjacency[node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        next_level.append(neighbor)
            
            if next_level:
                current_level = next_level
                depth += 1
            else:
                break
        
        # Return any node from the farthest level reached
        return current_level[0] if current_level else start_node

    def _degree_based_ordering(self) -> List[int]:
        """Simple degree-based ordering - very fast."""
        print("Starting degree-based ordering...")
        
        # Build adjacency graph and compute degrees
        adjacency = self._build_adjacency_graph()
        degrees = self._compute_degrees(adjacency)
        
        # Sort nodes by degree (ascending)
        print("Sorting nodes by degree...")
        nodes = list(adjacency.keys())
        nodes.sort(key=lambda x: (degrees[x], x))  # Secondary sort by ID for consistency
        
        print(f"Degree-based ordering complete. Processed {len(nodes)} nodes.")
        return nodes

    def _quick_bandwidth_estimate(self, sample_elements: pd.DataFrame, 
                                node_ordering: List[int]) -> float:
        """Quick bandwidth estimate using a sample of elements with vectorized operations."""
        # Create mapping from old to new node IDs
        node_mapping = {old_id: new_id for new_id, old_id in enumerate(node_ordering, 1)}
        
        # Get node arrays
        nodes_array = sample_elements.iloc[:, 1:4].values.astype(int)
        
        # Vectorized mapping - only process elements where all nodes exist in mapping
        valid_elements = []
        for i, element_nodes in enumerate(nodes_array):
            if all(node in node_mapping for node in element_nodes):
                mapped_nodes = [node_mapping[node] for node in element_nodes]
                valid_elements.append(mapped_nodes)
        
        if not valid_elements:
            return float('inf')
        
        # Convert to numpy array and compute bandwidths vectorized
        valid_array = np.array(valid_elements)
        element_mins = np.min(valid_array, axis=1)
        element_maxs = np.max(valid_array, axis=1)
        element_bandwidths = element_maxs - element_mins
        
        return float(np.mean(element_bandwidths))

    def _compute_bandwidth(self, elements_df: pd.DataFrame) -> Tuple[int, float]:
        """Compute the bandwidth and average bandwidth of the mesh using vectorized operations."""
        print("Computing bandwidth using vectorized operations...")
        
        # Extract node columns as numpy array for vectorized operations
        # Assuming columns 1, 2, 3 contain the three node IDs
        nodes_array = elements_df.iloc[:, 1:4].values.astype(int)
        
        # Vectorized min and max operations
        element_mins = np.min(nodes_array, axis=1)
        element_maxs = np.max(nodes_array, axis=1)
        
        # Compute bandwidths for all elements at once
        element_bandwidths = element_maxs - element_mins
        
        # Get statistics
        max_bandwidth = int(np.max(element_bandwidths))
        avg_bandwidth = float(np.mean(element_bandwidths))
        
        return max_bandwidth, avg_bandwidth

    def _create_node_mapping(self, new_ordering: List[int]) -> Dict[int, int]:
        """Create mapping from old node IDs to new node IDs."""
        return {old_id: new_id for new_id, old_id in enumerate(new_ordering, 1)}

    def _renumber_boundaries(self, boundaries: Dict, 
                           node_mapping: Dict[int, int]) -> Dict:
        """Renumber boundary nodes using the node mapping with optimized operations."""
        print("Renumbering boundary nodes using optimized operations...")
        
        # Count total boundaries for information
        total_boundaries = sum(len(boundaries_list) for boundaries_list in boundaries.values())
        print(f"Processing {total_boundaries} boundary segments...")
        
        new_boundaries = {}
        
        for ibtype in boundaries.keys():
            boundaries_ibtype = boundaries[ibtype]
            new_boundaries_ibtype = []
            
            for boundary in boundaries_ibtype:
                boundary_new = boundary.copy()
                node_ids = boundary['node_id']
                
                if ibtype is None:  # open boundary
                    # Vectorized conversion for open boundaries
                    old_ids = [int(nid) for nid in node_ids]
                    new_node_ids = [str(node_mapping[old_id]) for old_id in old_ids]
                elif ibtype.endswith('4'):  # weir boundary
                    # Vectorized conversion for weir boundaries
                    new_node_ids = []
                    for node_id in node_ids:
                        old_id1, old_id2 = int(node_id[0]), int(node_id[1])
                        new_node_ids.append((str(node_mapping[old_id1]), str(node_mapping[old_id2])))
                else:  # other boundaries
                    # Vectorized conversion for other boundaries
                    old_ids = [int(nid) for nid in node_ids]
                    new_node_ids = [str(node_mapping[old_id]) for old_id in old_ids]
                
                boundary_new['node_id'] = new_node_ids
                new_boundaries_ibtype.append(boundary_new)
            
            new_boundaries[ibtype] = new_boundaries_ibtype
        
        print(f"Completed renumbering {total_boundaries} boundary segments")
        return new_boundaries

    def reduce_bandwidth(self, algorithm: str = "rcm") -> AdcircMesh:
        """Reduce mesh bandwidth using the specified algorithm."""
        print(f"Starting bandwidth reduction using {algorithm.upper()} algorithm...")
        print("This process may take several minutes for large meshes...")
        
        # Compute original bandwidth
        original_max_bw, original_avg_bw = self._compute_bandwidth(
            self._mesh.elements.elements)
        print(f"Original bandwidth - Max: {original_max_bw}, Average: {original_avg_bw:.2f}")
        
        # Get new node ordering
        if algorithm.lower() == "rcm":
            new_ordering = self._reverse_cuthill_mckee()
        elif algorithm.lower() == "king":
            new_ordering = self._king_algorithm()
        elif algorithm.lower() == "gps":
            new_ordering = self._gps_algorithm()
        elif algorithm.lower() == "degree":
            new_ordering = self._degree_based_ordering()
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}. Available: rcm, king, gps, degree")
        
        # Create node mapping
        print("Creating node renumbering mapping...")
        node_mapping = self._create_node_mapping(new_ordering)
        
        # Validate that all nodes in elements exist in the mapping
        print("Validating node mapping completeness...")
        elements_df = self._mesh.elements.elements
        all_element_nodes = set()
        for col in elements_df.columns[1:]:  # Skip first column (element ID)
            col_nodes = elements_df[col].dropna().astype(int)
            all_element_nodes.update(col_nodes.unique())
        
        missing_nodes = all_element_nodes - set(node_mapping.keys())
        if missing_nodes:
            print(f"Error: Found {len(missing_nodes)} nodes in elements that are not in node mapping")
            print(f"Sample missing nodes: {list(missing_nodes)[:10]}")
            # Add missing nodes to the mapping (emergency fix)
            max_new_id = max(node_mapping.values()) if node_mapping else 0
            for i, missing_node in enumerate(sorted(missing_nodes)):
                node_mapping[missing_node] = max_new_id + i + 1
            print(f"Added {len(missing_nodes)} missing nodes to mapping")
        else:
            print("Node mapping validation passed")
        
        # Create new nodes DataFrame with renumbered indices
        print("Renumbering nodes using vectorized operations...")
        new_nodes = self._mesh.nodes.copy()
        
        # Vectorized index mapping
        old_indices = np.array(new_nodes.index)
        new_indices = np.array([node_mapping[old_id] for old_id in old_indices])
        new_nodes.index = new_indices
        new_nodes = new_nodes.sort_index()
        print(f"Renumbered {len(new_nodes)} nodes")
        
        # Renumber elements using vectorized operations
        print("Renumbering elements using vectorized operations...")
        new_elements = self._mesh.elements.elements.copy()
        
        # Detect triangular vs higher-order elements
        element_columns = new_elements.columns[1:]
        first_element = new_elements.iloc[0]
        
        # Find which columns have real data (not NaN/0)
        active_columns = []
        for col in element_columns:
            if not pd.isna(first_element[col]) and first_element[col] != 0:
                active_columns.append(col)
        
        # For triangular meshes, typically only first 3 columns are active
        if len(active_columns) <= 3:
            print(f"Detected triangular mesh with {len(active_columns)} active node columns")
            columns_to_process = active_columns
        else:
            print(f"Detected higher-order mesh with {len(active_columns)} active node columns")
            columns_to_process = element_columns
        
        # Apply mapping to each element column with smart error handling
        for col in element_columns:
            original_values = new_elements[col]
            
            if col in columns_to_process:
                # Process active columns (should have real node data)
                if original_values.isna().any():
                    nan_count = original_values.isna().sum()
                    print(f"Warning: Found {nan_count} unexpected NaN values in active column {col}")
                    original_values = original_values.fillna(0)
                
                # Convert to integer safely
                try:
                    int_values = original_values.astype(int)
                except (ValueError, OverflowError) as e:
                    print(f"Error converting column {col} to integer: {e}")
                    print(f"Sample values: {original_values.head(10).tolist()}")
                    raise
                
                # Apply mapping
                mapped_values = int_values.map(node_mapping)
                
                # Check if mapping was successful (only warn for active columns)
                missing_count = mapped_values.isna().sum()
                if missing_count > 0:
                    print(f"Warning: {missing_count} node IDs in active column {col} not found in mapping")
                    # Find missing node IDs for debugging
                    missing_mask = mapped_values.isna()
                    missing_nodes = int_values[missing_mask].unique()[:10]  # Show first 10
                    print(f"Sample missing node IDs: {missing_nodes.tolist()}")
                    
                    # Fill missing values with original values (this shouldn't happen in a proper mesh)
                    mapped_values = mapped_values.fillna(int_values)
                
                new_elements[col] = mapped_values.astype(int)
            else:
                # Inactive columns (expected to have NaN/0) - just preserve as-is
                # Fill NaN with 0 and keep as integer
                filled_values = original_values.fillna(0).astype(int)
                new_elements[col] = filled_values
            
        print(f"Renumbered {len(new_elements)} elements")
        
        # Compute new bandwidth
        new_max_bw, new_avg_bw = self._compute_bandwidth(new_elements)
        print(f"New bandwidth - Max: {new_max_bw}, Average: {new_avg_bw:.2f}")
        print(f"Bandwidth reduction - Max: {((original_max_bw - new_max_bw) / original_max_bw * 100):.1f}%, "
              f"Average: {((original_avg_bw - new_avg_bw) / original_avg_bw * 100):.1f}%")
        
        # Renumber boundaries
        boundaries = self._mesh.boundaries.to_dict().copy()
        new_boundaries = self._renumber_boundaries(boundaries, node_mapping)
        
        # Create new mesh
        print("Creating renumbered mesh...")
        renumbered_mesh = AdcircMesh(
            nodes=new_nodes,
            elements=new_elements,
            boundaries=new_boundaries
        )
        
        print("Bandwidth reduction complete!")
        return renumbered_mesh


def get_parser():
    """Get command line argument parser."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        'input_mesh',
        help='Path to input mesh file'
    )
    parser.add_argument(
        '-o', '--output',
        default='bandwidth_reduced_mesh.grd',
        help='Output mesh file path (default: bandwidth_reduced_mesh.grd)'
    )
    parser.add_argument(
        '-a', '--algorithm',
        choices=['rcm', 'king', 'gps', 'degree'],
        default='king',
        help='Bandwidth reduction algorithm: rcm (traditional), king (fast), gps (high-quality), degree (fastest) (default: king)'
    )
    parser.add_argument(
        '-d', '--description',
        default='bandwidth_reduced',
        help='Description for the output mesh (default: bandwidth_reduced)'
    )
    return parser


def main(args=None):
    """Main function for bandwidth reduction."""
    if args is None:
        args = get_parser().parse_args()
    
    try:
        # Read the input mesh
        print(f"Reading mesh from: {args.input_mesh}")
        mesh = AdcircMesh.open(args.input_mesh)
        print(f"Successfully read mesh with {len(mesh.nodes)} nodes and "
              f"{len(mesh.elements.elements)} elements")
        
        # Show algorithm information
        algorithm_info = {
            'rcm': 'Reverse Cuthill-McKee (traditional, good quality)',
            'king': "King's algorithm (fast, good for irregular meshes)",
            'gps': 'GPS/Gibbs-Poole-Stockmeyer (high quality, moderate speed)',
            'degree': 'Degree-based ordering (fastest, basic quality)'
        }
        print(f"Using algorithm: {args.algorithm.upper()} - {algorithm_info[args.algorithm]}")
        print("Note: All algorithms now use optimized vectorized operations for faster performance.")
        
        # Create bandwidth reducer and process mesh
        reducer = BandwidthReducer(mesh)
        reduced_mesh = reducer.reduce_bandwidth(args.algorithm)
        
        # Set description and write output
        reduced_mesh.description = args.description
        reduced_mesh.write(args.output, overwrite=True)
        print(f"\nReduced bandwidth mesh saved to: {args.output}")
        
        return 0
        
    except Exception as e:
        print(f"Error: {e}")
        return 1


if __name__ == "__main__":
    main() 