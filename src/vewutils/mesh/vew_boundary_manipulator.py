#!/usr/bin/env python3
"""
Module for manipulating VEW boundaries in ADCIRC meshes.
"""

from typing import Dict, Tuple
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from adcircpy import AdcircMesh

class VEWBoundaryManipulator:
    """Class for manipulating VEW boundaries in ADCIRC meshes."""
    
    @staticmethod
    def find_matching_nodes(nodes_df1: pd.DataFrame, nodes_df2: pd.DataFrame, 
                          tolerance: float) -> Dict[int, int]:
        """
        Find matching nodes between two meshes within tolerance.
        
        Args:
            nodes_df1: First mesh nodes DataFrame
            nodes_df2: Second mesh nodes DataFrame
            tolerance: Distance tolerance for matching nodes
            
        Returns:
            Dictionary mapping first mesh node IDs to second mesh node IDs
        """
        print("Finding matching nodes...")
        coords1 = nodes_df1[['x', 'y']].values
        coords2 = nodes_df2[['x', 'y']].values
        matching_nodes = {}

        # Build KD-tree for coords2
        tree = cKDTree(coords2)
        
        # Query KD-tree for all points in coords1 at once
        distances, indices = tree.query(coords1, distance_upper_bound=tolerance)
        
        # Create matching_nodes dictionary for valid matches
        valid_matches = distances < tolerance
        for i, (valid, idx) in enumerate(zip(valid_matches, indices)):
            if valid:
                # Convert to 1-based indexing
                matching_nodes[i + 1] = idx + 1
        
        print(f"Found {len(matching_nodes)} matching nodes")
        return matching_nodes

    @staticmethod
    def get_vew_node_pairs(mesh: AdcircMesh) -> Tuple[Dict[int, int], Dict[int, int]]:
        """Extract bank and channel node pairs from VEW boundaries.
        
        Args:
            mesh: AdcircMesh object
            
        Returns:
            tuple: (bank_to_channel, channel_to_bank) dictionaries mapping node IDs
        """
        bank_to_channel = {}
        channel_to_bank = {}
        boundaries = mesh.boundaries.to_dict()
        
        if '64' in boundaries:  # Check for VEW boundaries (ibtype 64)
            for vew_boundary in boundaries['64']:
                for node_pair in vew_boundary['node_id']:
                    bank_node = int(node_pair[0])
                    channel_node = int(node_pair[1])
                    bank_to_channel[bank_node] = channel_node
                    channel_to_bank[channel_node] = bank_node
                    
        return bank_to_channel, channel_to_bank
    
    @staticmethod
    def lower_channel_elevations_above_banks(mesh: AdcircMesh, tolerance: float) -> AdcircMesh:
        """Lower channel node elevations only when they are above their corresponding bank nodes.
        
        Args:
            mesh: AdcircMesh object
            tolerance: Amount to lower channel node elevations below bank node elevations (in meters)
            
        Returns:
            Modified AdcircMesh object
        """
        bank_to_channel, channel_to_bank = VEWBoundaryManipulator.get_vew_node_pairs(mesh)
        
        # Count number of altered channel nodes
        altered_count = 0
        
        # Adjust channel node elevations
        for channel_node, bank_node in channel_to_bank.items():
            bank_elevation = mesh.nodes.loc[bank_node, 'value_1']
            channel_elevation = mesh.nodes.loc[channel_node, 'value_1']
            target_elevation = bank_elevation - tolerance
            # Only lower the channel elevation if it's greater than the target elevation
            if channel_elevation > target_elevation:
                mesh.nodes.loc[channel_node, 'value_1'] = target_elevation
                altered_count += 1
        
        print(f"Number of channel nodes altered: {altered_count}")
        return mesh
    
    @staticmethod
    def get_vew_boundaries(mesh: AdcircMesh) -> list:
        """Get all VEW boundaries from the mesh.
        
        Args:
            mesh: AdcircMesh object
            
        Returns:
            List of VEW boundary dictionaries
        """
        boundaries = mesh.boundaries.to_dict()
        return boundaries.get('64', [])
    
    @staticmethod
    def add_vew_boundary(mesh: AdcircMesh, 
                        node_pairs: list,
                        barrier_heights: list,
                        subcritical_coefficients: list = None,
                        supercritical_coefficients: list = None) -> AdcircMesh:
        """Add a new VEW boundary to the mesh.
        
        Args:
            mesh: AdcircMesh object
            node_pairs: List of (bank_node, channel_node) tuples
            barrier_heights: List of barrier heights for each node pair
            subcritical_coefficients: List of subcritical flow coefficients (default: 1.0)
            supercritical_coefficients: List of supercritical flow coefficients (default: 1.0)
            
        Returns:
            Modified AdcircMesh object
        """
        if subcritical_coefficients is None:
            subcritical_coefficients = [1.0] * len(node_pairs)
        if supercritical_coefficients is None:
            supercritical_coefficients = [1.0] * len(node_pairs)
            
        vew_boundary = {
            'node_id': [(str(bank), str(channel)) for bank, channel in node_pairs],
            'barrier_height': barrier_heights,
            'subcritical_flow_coefficient': subcritical_coefficients,
            'supercritical_flow_coefficient': supercritical_coefficients
        }
        
        boundaries = mesh.boundaries.to_dict()
        if '64' in boundaries:
            boundaries['64'].append(vew_boundary)
        else:
            boundaries['64'] = [vew_boundary]
            
        return AdcircMesh(nodes=mesh.nodes, elements=mesh.elements.elements, boundaries=boundaries)
    
    @staticmethod
    def ensure_barrier_heights_above_banks(mesh: AdcircMesh, tolerance: float = 0.001) -> AdcircMesh:
        """Set barrier heights in VEW boundaries to bank elevations plus tolerance.
        
        Args:
            mesh: AdcircMesh object
            tolerance: Amount to add to bank elevations for barrier heights (in meters, default: 0.001)
            
        Returns:
            Modified AdcircMesh object
        """
        bank_to_channel, _ = VEWBoundaryManipulator.get_vew_node_pairs(mesh)
        boundaries = mesh.boundaries.to_dict()
        
        # Count number of barrier heights set
        set_count = 0
        
        if '64' in boundaries:  # Check for VEW boundaries (ibtype 64)
            for vew_boundary in boundaries['64']:
                for i, node_pair in enumerate(vew_boundary['node_id']):
                    bank_node = int(node_pair[0])
                    bank_elevation = mesh.nodes.loc[bank_node, 'value_1']
                    target_height = bank_elevation + tolerance
                    
                    # Set barrier height to bank elevation plus tolerance
                    vew_boundary['barrier_height'][i] = target_height
                    set_count += 1
        
        print(f"Number of barrier heights set: {set_count}")
        return mesh 