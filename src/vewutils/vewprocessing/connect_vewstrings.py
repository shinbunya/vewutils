#!/usr/bin/env python3
"""
Module for connecting VEW strings from two YAML files by matching coordinates.
"""

import argparse
import yaml
import numpy as np
import copy
from typing import Dict, List, Tuple, Optional


def _calculate_distance(node1: Dict, node2: Dict) -> float:
    """Calculate Euclidean distance between two nodes."""
    dx = node1['x'] - node2['x']
    dy = node1['y'] - node2['y']
    return np.sqrt(dx * dx + dy * dy)


def _find_connections(vewstring_a: List[Dict], vewstring_b: List[Dict], 
                     threshold_dist: float) -> Tuple[Optional[Tuple[str, int]], Optional[Tuple[str, int]]]:
    """
    Find head and tail connections between two vewstrings separately.
    
    Returns:
        Tuple of (head_connection, tail_connection) where each can be:
        - None if no connection
        - Tuple of (connection_type, index_in_b) if connection found
        connection_type can be 'head-head', 'head-tail', 'tail-head', or 'tail-tail'.
    """
    if not vewstring_a or not vewstring_b:
        return (None, None)
    
    head_a = vewstring_a[0]
    tail_a = vewstring_a[-1]
    head_b = vewstring_b[0]
    tail_b = vewstring_b[-1]
    
    head_connection = None
    tail_connection = None
    
    # Check head connections
    if _calculate_distance(head_a, head_b) <= threshold_dist:
        head_connection = ('head-head', 0)
    elif _calculate_distance(head_a, tail_b) <= threshold_dist:
        head_connection = ('head-tail', -1)
    
    # Check tail connections
    if _calculate_distance(tail_a, head_b) <= threshold_dist:
        tail_connection = ('tail-head', 0)
    elif _calculate_distance(tail_a, tail_b) <= threshold_dist:
        tail_connection = ('tail-tail', -1)
    
    return (head_connection, tail_connection)


def _copy_node(node: Dict) -> Dict:
    """Create a deep copy of a node dictionary to avoid YAML anchors/aliases."""
    return copy.deepcopy(node)


def _merge_vewstrings_single(vewstring_a: List[Dict], vewstring_b: List[Dict], 
                             connection_type: str) -> List[Dict]:
    """
    Merge two vewstrings based on a single connection type.
    
    Args:
        vewstring_a: First vewstring
        vewstring_b: Second vewstring
        connection_type: Type of connection ('head-head', 'head-tail', 'tail-head', 'tail-tail')
    
    Returns:
        Merged vewstring with duplicate node removed.
    """
    # Create deep copies to avoid YAML anchors/aliases
    vewstring_a_copy = [_copy_node(node) for node in vewstring_a]
    vewstring_b_copy = [_copy_node(node) for node in vewstring_b]
    
    if connection_type == 'head-head':
        # Reverse B, then reversed(B)[1:] + A
        reversed_b = list(reversed(vewstring_b_copy))
        result = reversed_b[1:] + vewstring_a_copy
    
    elif connection_type == 'head-tail':
        # B[1:] + A
        result = vewstring_b_copy[1:] + vewstring_a_copy
    
    elif connection_type == 'tail-head':
        # A + B[1:]
        result = vewstring_a_copy + vewstring_b_copy[1:]
    
    elif connection_type == 'tail-tail':
        # A + reversed(B)[1:]
        reversed_b = list(reversed(vewstring_b_copy))
        result = vewstring_a_copy + reversed_b[1:]
    
    else:
        raise ValueError(f"Unknown connection type: {connection_type}")
    
    return result


def _merge_vewstrings_multi(vewstring_a: List[Dict], vewstring_b1: List[Dict], 
                            head_conn_type: str, vewstring_b2: List[Dict], 
                            tail_conn_type: str) -> List[Dict]:
    """
    Merge vewstring A with two different vewstrings B1 (at head) and B2 (at tail).
    
    Args:
        vewstring_a: Middle vewstring
        vewstring_b1: Vewstring connected at head of A
        head_conn_type: Connection type at head ('head-head' or 'head-tail')
        vewstring_b2: Vewstring connected at tail of A
        tail_conn_type: Connection type at tail ('tail-head' or 'tail-tail')
    
    Returns:
        Merged vewstring: B1 + A + B2 (with appropriate reversals)
    """
    # Create deep copies to avoid YAML anchors/aliases
    vewstring_a_copy = [_copy_node(node) for node in vewstring_a]
    vewstring_b1_copy = [_copy_node(node) for node in vewstring_b1]
    vewstring_b2_copy = [_copy_node(node) for node in vewstring_b2]
    
    # Build the head part (B1 connected to head of A)
    if head_conn_type == 'head-head':
        # Reverse B1, then reversed(B1)[1:] + A
        reversed_b1 = list(reversed(vewstring_b1_copy))
        head_part = reversed_b1[1:] + vewstring_a_copy
    elif head_conn_type == 'head-tail':
        # B1[1:] + A
        head_part = vewstring_b1_copy[1:] + vewstring_a_copy
    else:
        raise ValueError(f"Invalid head connection type: {head_conn_type}")
    
    # Build the tail part (B2 connected to tail of A, which is now tail of head_part)
    if tail_conn_type == 'tail-head':
        # head_part + B2[1:]
        return head_part + vewstring_b2_copy[1:]
    elif tail_conn_type == 'tail-tail':
        # head_part + reversed(B2)[1:]
        reversed_b2 = list(reversed(vewstring_b2_copy))
        return head_part + reversed_b2[1:]
    else:
        raise ValueError(f"Invalid tail connection type: {tail_conn_type}")


def _merge_vewstrings_loop(vewstring_a: List[Dict], vewstring_b: List[Dict], 
                           head_conn_type: str, tail_conn_type: str) -> List[Dict]:
    """
    Merge vewstring A with the same vewstring B at both head and tail (forming a loop).
    
    The result should include all unique nodes from A and B exactly once.
    Format: Some orientation of B + A's internal nodes (skip both A's head and tail since they're on B)
    
    Args:
        vewstring_a: Vewstring connecting to B at both ends
        vewstring_b: Larger vewstring that A connects to
        head_conn_type: Connection type at A's head ('head-head' or 'head-tail')
        tail_conn_type: Connection type at A's tail ('tail-head' or 'tail-tail')
    
    Returns:
        Merged vewstring with all nodes from A and B, avoiding duplicates
    """
    # Create deep copies to avoid YAML anchors/aliases
    vewstring_a_copy = [_copy_node(node) for node in vewstring_a]
    vewstring_b_copy = [_copy_node(node) for node in vewstring_b]
    
    # Strategy: Create a loop by inserting A's internal nodes between the matching endpoints of B
    # The result should be: B_start → A_internal → B_end → B_internal → (loops back)
    
    # Determine which end of B connects to A's head
    if head_conn_type == 'head-head' and tail_conn_type == 'tail-tail':
        # A's head connects to B's head, A's tail connects to B's tail
        # Path: B_head → A_internal → B_tail → reversed(B_internal) (to return to B_head)
        b_internal_reversed = list(reversed(vewstring_b_copy[1:-1]))
        result = [vewstring_b_copy[0]] + vewstring_a_copy[1:-1] + vewstring_b_copy[-1:] + b_internal_reversed
        # For a closed loop, we need to add the head node at the tail to complete the loop
        if len(result) > 1 and result[0]['node_id'] != result[-1]['node_id']:
            result.append(_copy_node(result[0]))
        
    elif head_conn_type == 'head-head' and tail_conn_type == 'tail-head':
        # A's head connects to B's head, A's tail connects to B's head (need to reverse A)
        # Path: B_head → reversed(A_internal) → B_head (same point, forms tight loop with B)
        reversed_a_internal = list(reversed(vewstring_a_copy[1:-1]))
        result = vewstring_b_copy + reversed_a_internal
        # For a closed loop, we need to add the head node at the tail to complete the loop
        if len(result) > 1 and result[0]['node_id'] != result[-1]['node_id']:
            result.append(_copy_node(result[0]))
        
    elif head_conn_type == 'head-tail' and tail_conn_type == 'tail-head':
        # A's head connects to B's tail, A's tail connects to B's head
        # Path: B_head → reversed(A_internal) → B_tail → reversed(B_internal) (to return to B_head)
        reversed_a_internal = list(reversed(vewstring_a_copy[1:-1]))
        b_internal_reversed = list(reversed(vewstring_b_copy[1:-1]))
        result = [vewstring_b_copy[0]] + reversed_a_internal + vewstring_b_copy[-1:] + b_internal_reversed
        # For a closed loop, we need to add the head node at the tail to complete the loop
        if len(result) > 1 and result[0]['node_id'] != result[-1]['node_id']:
            result.append(_copy_node(result[0]))
        
    elif head_conn_type == 'head-tail' and tail_conn_type == 'tail-tail':
        # A's head connects to B's tail, A's tail connects to B's tail (need to reverse A)
        # This means A forms a tight loop at B's tail
        reversed_a_internal = list(reversed(vewstring_a_copy[1:-1]))
        result = vewstring_b_copy + reversed_a_internal
        # For a closed loop, we need to add the head node at the tail to complete the loop
        if len(result) > 1 and result[0]['node_id'] != result[-1]['node_id']:
            result.append(_copy_node(result[0]))
        
    else:
        raise ValueError(f"Invalid connection types: head={head_conn_type}, tail={tail_conn_type}")
    
    return result


def connect_vewstrings(vewstrings_a: List[List[Dict]], vewstrings_b: List[List[Dict]], 
                      threshold_dist: float) -> List[List[Dict]]:
    """
    Connect vewstrings from two lists using a generic graph-based approach.
    
    Algorithm:
    1. Combine all vewstrings into a single list with source tracking
    2. Build a connection graph (which endpoints connect to which)
    3. Find connected components using graph traversal
    4. For each component, merge all strings in the correct order
    
    Args:
        vewstrings_a: List of vewstrings from first file
        vewstrings_b: List of vewstrings from second file
        threshold_dist: Maximum distance for coordinate matching
    
    Returns:
        List of connected and unconnected vewstrings.
    """
    # Combine all vewstrings into one list with metadata
    all_strings = []
    string_meta = []  # (source: 'a' or 'b', original_index: int)
    
    for i, vew in enumerate(vewstrings_a):
        all_strings.append(vew)
        string_meta.append(('a', i))
    for j, vew in enumerate(vewstrings_b):
        all_strings.append(vew)
        string_meta.append(('b', j))
    
    n_strings = len(all_strings)
    
    # Build connection graph: connections[i] = list of (j, connection_type, endpoint_i, endpoint_j)
    # endpoint_i and endpoint_j are 'head' or 'tail'
    connections = {}  # i -> list of (j, conn_type, end_i, end_j)
    
    for i in range(n_strings):
        connections[i] = []
        for j in range(i + 1, n_strings):
            hc, tc = _find_connections(all_strings[i], all_strings[j], threshold_dist)
            
            if hc is not None:
                # i's head connects to j
                conn_type = hc[0]
                if conn_type == 'head-head':
                    connections[i].append((j, conn_type, 'head', 'head'))
                else:  # head-tail
                    connections[i].append((j, conn_type, 'head', 'tail'))
            
            if tc is not None:
                # i's tail connects to j
                conn_type = tc[0]
                if conn_type == 'tail-head':
                    connections[i].append((j, conn_type, 'tail', 'head'))
                else:  # tail-tail
                    connections[i].append((j, conn_type, 'tail', 'tail'))
    
    # Validate: each string's endpoint should connect to at most one other string's endpoint
    for i in range(n_strings):
        head_conns = [c for c in connections[i] if c[2] == 'head']
        tail_conns = [c for c in connections[i] if c[2] == 'tail']
        
        if len(head_conns) > 1:
            other_strings = [c[0] for c in head_conns]
            raise ValueError(
                f"Head of string {i} ({string_meta[i]}) connects to multiple strings: {other_strings}"
            )
        if len(tail_conns) > 1:
            other_strings = [c[0] for c in tail_conns]
            raise ValueError(
                f"Tail of string {i} ({string_meta[i]}) connects to multiple strings: {other_strings}"
            )
    
    # Also check reverse connections (j -> i)
    reverse_connections = {}  # j -> list of (i, conn_type, end_j, end_i)
    for i in range(n_strings):
        for j, conn_type, end_i, end_j in connections[i]:
            reverse_connections.setdefault(j, []).append((i, conn_type, end_j, end_i))
    
    for j in range(n_strings):
        if j in reverse_connections:
            head_conns = [c for c in reverse_connections[j] if c[2] == 'head']
            tail_conns = [c for c in reverse_connections[j] if c[2] == 'tail']
            
            if len(head_conns) > 1:
                other_strings = [c[0] for c in head_conns]
                raise ValueError(
                    f"Head of string {j} ({string_meta[j]}) is connected by multiple strings: {other_strings}"
                )
            if len(tail_conns) > 1:
                other_strings = [c[0] for c in tail_conns]
                raise ValueError(
                    f"Tail of string {j} ({string_meta[j]}) is connected by multiple strings: {other_strings}"
                )

    # Find connected components using graph traversal
    visited = set()
    components = []  # List of component lists: each component is a list of string indices
    
    def dfs_component(start_idx: int, component: List[int]):
        """DFS to find all strings in a connected component."""
        if start_idx in visited:
            return
        visited.add(start_idx)
        component.append(start_idx)
        
        # Follow connections from this string
        for j, conn_type, end_i, end_j in connections.get(start_idx, []):
            if j not in visited:
                dfs_component(j, component)
        
        # Follow reverse connections (strings that connect TO this one)
        for i, conn_type, end_j, end_i in reverse_connections.get(start_idx, []):
            if i not in visited:
                dfs_component(i, component)
    
    # Find all connected components
    for i in range(n_strings):
        if i not in visited:
            component = []
            dfs_component(i, component)
            if component:
                components.append(component)
    
    # Merge each component into a single string
    result = []
    for comp_idx, component in enumerate(components):
        if len(component) == 1:
            # Single string, no merging needed
            i = component[0]
            result.append([_copy_node(node) for node in all_strings[i]])
            continue
        
        # Multiple strings: need to merge them in the correct order
        # Build a path through all strings in the component
        # Start with any string, then follow connections
        
        # Find a starting string (one with only one connection, or any)
        start_strings = [i for i in component if len(connections.get(i, [])) + len(reverse_connections.get(i, [])) <= 1]
        if not start_strings:
            start_strings = [component[0]]
        
        start = start_strings[0]
        merged_path = []  # List of (string_idx, reversed: bool)
        used = set()
        
        def build_path(current: int, current_reversed: bool):
            """Build a path through connected strings."""
            if current in used:
                return
            used.add(current)
            merged_path.append((current, current_reversed))
            
            # Find next string connected to current's tail (if not reversed) or head (if reversed)
            current_end = 'tail' if not current_reversed else 'head'
            
            # Look for connections from current's end
            for j, conn_type, end_i, end_j in connections.get(current, []):
                if end_i == current_end and j not in used:
                    # Determine if j needs to be reversed
                    j_reversed = (end_j == 'head')  # If connecting to j's head, reverse j
                    build_path(j, j_reversed)
                    return
            
            # Look for reverse connections (strings connecting TO current)
            for i, conn_type, end_j, end_i in reverse_connections.get(current, []):
                if end_i == current_end and i not in used:
                    # Determine if i needs to be reversed
                    i_reversed = (end_j == 'tail')  # If i's tail connects to current, don't reverse
                    build_path(i, i_reversed)
                    return
        
        # Try building path from start
        build_path(start, False)
        
        # If we didn't get all strings, try from the other end
        if len(used) < len(component):
            # Find unused strings and try connecting them
            unused = [i for i in component if i not in used]
            for u in unused:
                build_path(u, False)
        
        # If still not all, there might be branches - handle them
        if len(used) < len(component):
            # For now, add remaining as separate segments
            for i in component:
                if i not in used:
                    merged_path.append((i, False))
        
        # Now merge the path
        if not merged_path:
            continue
        
        # Start with first string
        first_idx, first_rev = merged_path[0]
        current_string = list(all_strings[first_idx])
        if first_rev:
            current_string = list(reversed(current_string))
        current_string = [_copy_node(node) for node in current_string]
        
        # Merge remaining strings
        for string_idx, reversed_flag in merged_path[1:]:
            next_string = list(all_strings[string_idx])
            if reversed_flag:
                next_string = list(reversed(next_string))
            next_string = [_copy_node(node) for node in next_string]
            
            # Find connection between current_string and next_string
            hc, tc = _find_connections(current_string, next_string, threshold_dist)
            
            if hc is not None and tc is not None:
                # Both ends connect - this forms a loop
                # Use the loop merge function to properly handle the topology
                merged = _merge_vewstrings_loop(current_string, next_string, hc[0], tc[0])
                current_string = merged
            elif hc is not None:
                # Head connection only
                merged = _merge_vewstrings_single(current_string, next_string, hc[0])
                current_string = merged
            elif tc is not None:
                # Tail connection only
                merged = _merge_vewstrings_single(current_string, next_string, tc[0])
                current_string = merged
            else:
                # No connection found - append as separate segment (shouldn't happen)
                current_string = current_string + next_string
        
        result.append(current_string)
    
    return result


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False, 
        description="Connect VEW strings from two YAML files by matching coordinates"
    )
    parser.add_argument(
        "vewstring_file_a",
        help="Path to the first VEW string YAML file"
    )
    parser.add_argument(
        "vewstring_file_b",
        help="Path to the second VEW string YAML file"
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Path to save the connected VEW strings in YAML format'
    )
    parser.add_argument(
        '--threshold-dist',
        type=float,
        default=1e-6,
        help='Maximum allowable distance between two nodes to decide if they are collocated (default: 1e-6)'
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    import yaml
    
    print("=== VEW String Connector ===")
    print(f"Input file A: {args.vewstring_file_a}")
    print(f"Input file B: {args.vewstring_file_b}")
    print(f"Output file: {args.output}")
    print(f"Threshold distance: {args.threshold_dist}")
    print()
    
    # Read YAML files
    print("Loading VEW string files...")
    with open(args.vewstring_file_a, 'r') as f:
        data_a = yaml.safe_load(f)
    with open(args.vewstring_file_b, 'r') as f:
        data_b = yaml.safe_load(f)
    
    # Validate structure
    if 'vewstrings' not in data_a:
        raise ValueError(f"File {args.vewstring_file_a} does not contain 'vewstrings' key")
    if 'vewstrings' not in data_b:
        raise ValueError(f"File {args.vewstring_file_b} does not contain 'vewstrings' key")
    
    vewstrings_a = data_a['vewstrings']
    vewstrings_b = data_b['vewstrings']
    
    if not isinstance(vewstrings_a, list):
        raise ValueError(f"File {args.vewstring_file_a}: 'vewstrings' must be a list")
    if not isinstance(vewstrings_b, list):
        raise ValueError(f"File {args.vewstring_file_b}: 'vewstrings' must be a list")
    
    print(f"✓ File A loaded: {len(vewstrings_a)} VEW strings")
    print(f"✓ File B loaded: {len(vewstrings_b)} VEW strings")
    print()
    
    # Connect vewstrings
    print("Connecting VEW strings...")
    connected_vewstrings = connect_vewstrings(vewstrings_a, vewstrings_b, args.threshold_dist)
    print()
    
    # Write output
    print("Saving connected VEW strings to YAML...")
    output_data = {'vewstrings': connected_vewstrings}
    with open(args.output, 'w') as f:
        yaml.dump(output_data, f, sort_keys=False, default_flow_style=False)
    print(f"✓ Connected VEW strings saved to: {args.output}")
    
    total_vew_nodes = sum(len(vew) for vew in connected_vewstrings)
    print(f"  Total VEW strings: {len(connected_vewstrings)}")
    print(f"  Total VEW nodes: {total_vew_nodes}")
    print()
    print("=== Processing Complete ===")
    
    return 0


if __name__ == "__main__":
    main()

