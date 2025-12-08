#!/usr/bin/env python3
"""
Module for converting VEW string YAML files to GeoJSON format.
"""

import argparse
import yaml
import json
from typing import Dict, List


def yaml_to_geojson(vewstrings: List[List[Dict]]) -> Dict:
    """
    Convert vewstrings from YAML format to GeoJSON format.
    
    Args:
        vewstrings: List of vewstrings, where each vewstring is a list of node dictionaries
    
    Returns:
        GeoJSON FeatureCollection dictionary
    """
    features = []
    
    for i, vewstring in enumerate(vewstrings):
        if not vewstring:
            continue
        
        # Extract coordinates: [[x1, y1], [x2, y2], ...]
        coordinates = []
        # Extract node attributes
        nodes = []
        
        for node in vewstring:
            # Validate required fields
            if 'x' not in node or 'y' not in node:
                raise ValueError(f"Node in vewstring {i} is missing 'x' or 'y' coordinates: {node}")
            if 'node_id' not in node:
                raise ValueError(f"Node in vewstring {i} is missing 'node_id': {node}")
            if 'bank_elevation' not in node:
                raise ValueError(f"Node in vewstring {i} is missing 'bank_elevation': {node}")
            if 'bank_mannings_n' not in node:
                raise ValueError(f"Node in vewstring {i} is missing 'bank_mannings_n': {node}")
            
            # Extract coordinates
            coordinates.append([float(node['x']), float(node['y'])])
            
            # Extract node attributes
            nodes.append({
                'node_id': int(node['node_id']),
                'bank_elevation': float(node['bank_elevation']),
                'bank_mannings_n': float(node['bank_mannings_n'])
            })
        
        # Create LineString feature
        feature = {
            'type': 'Feature',
            'geometry': {
                'type': 'LineString',
                'coordinates': coordinates
            },
            'properties': {
                'nodes': nodes
            }
        }
        
        features.append(feature)
    
    # Create FeatureCollection
    geojson = {
        'type': 'FeatureCollection',
        'features': features
    }
    
    return geojson


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Convert VEW string YAML files to GeoJSON format"
    )
    parser.add_argument(
        "input_yaml",
        help="Path to the input VEW string YAML file"
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Path to save the output GeoJSON file'
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    import yaml
    import json
    
    print("=== VEW String YAML to GeoJSON Converter ===")
    print(f"Input YAML file: {args.input_yaml}")
    print(f"Output GeoJSON file: {args.output}")
    print()
    
    # Read YAML file
    print("Loading VEW string YAML file...")
    with open(args.input_yaml, 'r') as f:
        data = yaml.safe_load(f)
    
    # Validate structure
    if 'vewstrings' not in data:
        raise ValueError(f"File {args.input_yaml} does not contain 'vewstrings' key")
    
    vewstrings = data['vewstrings']
    
    if not isinstance(vewstrings, list):
        raise ValueError(f"File {args.input_yaml}: 'vewstrings' must be a list")
    
    print(f"✓ YAML file loaded: {len(vewstrings)} VEW strings")
    print()
    
    # Convert to GeoJSON
    print("Converting to GeoJSON format...")
    geojson = yaml_to_geojson(vewstrings)
    print(f"✓ Converted {len(geojson['features'])} VEW strings to GeoJSON features")
    print()
    
    # Write output
    print("Saving GeoJSON file...")
    with open(args.output, 'w') as f:
        json.dump(geojson, f, indent=2, sort_keys=False)
    print(f"✓ GeoJSON file saved to: {args.output}")
    
    total_nodes = sum(len(feature['properties']['nodes']) for feature in geojson['features'])
    print(f"  Total features: {len(geojson['features'])}")
    print(f"  Total nodes: {total_nodes}")
    print()
    print("=== Conversion Complete ===")
    
    return 0


if __name__ == "__main__":
    main()

