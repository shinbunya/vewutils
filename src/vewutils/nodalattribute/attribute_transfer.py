import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from adcircpy import AdcircMesh
from adcircpy.mesh.fort13 import NodalAttributes

class AttributeTransfer:
    def __init__(self):
        pass

    def _get_vew_node_pairs(self, mesh):
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

    def _build_conversion_tables(self, source_mesh, target_mesh):
        """Build conversion tables between source and target nodes.
        
        Args:
            source_mesh: Source mesh object
            target_mesh: Target mesh object
            
        Returns:
            tuple: (source_to_target, target_to_source) dictionaries mapping node IDs
        """
        # Get VEW node pairs for both meshes
        source_bank_to_channel, source_channel_to_bank = self._get_vew_node_pairs(source_mesh)
        target_bank_to_channel, target_channel_to_bank = self._get_vew_node_pairs(target_mesh)
        
        # Build KDTree for source nodes
        source_coords = np.column_stack((source_mesh.nodes['x'], source_mesh.nodes['y']))
        source_tree = cKDTree(source_coords)
        
        # Build KDTree for target nodes
        target_coords = np.column_stack((target_mesh.nodes['x'], target_mesh.nodes['y']))
        target_tree = cKDTree(target_coords)
        
        # Initialize conversion tables
        source_to_target = {}
        target_to_source = {}
        
        # Build source_to_target mapping
        for source_idx in range(len(source_mesh.nodes)):
            source_node_id = source_idx + 1  # Convert to 1-based ID
            
            # Check if source node is a bank node
            if source_node_id in source_bank_to_channel:
                # Find nearest target bank node
                source_coord = source_coords[source_idx]
                _, nearest_target_idx = target_tree.query(source_coord)
                nearest_target_id = nearest_target_idx + 1  # Convert to 1-based ID
                
                # If nearest target node is a channel node, find its paired bank node
                if nearest_target_id in target_channel_to_bank:
                    nearest_target_id = target_channel_to_bank[nearest_target_id]
                
                source_to_target[source_node_id] = nearest_target_id
            
            # Check if source node is a channel node
            elif source_node_id in source_channel_to_bank:
                # Find nearest target channel node
                source_coord = source_coords[source_idx]
                _, nearest_target_idx = target_tree.query(source_coord)
                nearest_target_id = nearest_target_idx + 1  # Convert to 1-based ID
                
                # If nearest target node is a bank node, find its paired channel node
                if nearest_target_id in target_bank_to_channel:
                    nearest_target_id = target_bank_to_channel[nearest_target_id]
                
                source_to_target[source_node_id] = nearest_target_id
            
            else:
                # Regular node, find nearest target node
                source_coord = source_coords[source_idx]
                _, nearest_target_idx = target_tree.query(source_coord)
                nearest_target_id = nearest_target_idx + 1  # Convert to 1-based ID
                source_to_target[source_node_id] = nearest_target_id
        
        # Build target_to_source mapping
        for target_idx in range(len(target_mesh.nodes)):
            target_node_id = target_idx + 1  # Convert to 1-based ID
            
            # Check if target node is a bank node
            if target_node_id in target_bank_to_channel:
                # Find nearest source bank node
                target_coord = target_coords[target_idx]
                _, nearest_source_idx = source_tree.query(target_coord)
                nearest_source_id = nearest_source_idx + 1  # Convert to 1-based ID
                
                # If nearest source node is a channel node, find its paired bank node
                if nearest_source_id in source_channel_to_bank:
                    nearest_source_id = source_channel_to_bank[nearest_source_id]
                
                target_to_source[target_node_id] = nearest_source_id
            
            # Check if target node is a channel node
            elif target_node_id in target_channel_to_bank:
                # Find nearest source channel node
                target_coord = target_coords[target_idx]
                _, nearest_source_idx = source_tree.query(target_coord)
                nearest_source_id = nearest_source_idx + 1  # Convert to 1-based ID
                
                # If nearest source node is a bank node, find its paired channel node
                if nearest_source_id in source_bank_to_channel:
                    nearest_source_id = source_bank_to_channel[nearest_source_id]
                
                target_to_source[target_node_id] = nearest_source_id
            
            else:
                # Regular node, find nearest source node
                target_coord = target_coords[target_idx]
                _, nearest_source_idx = source_tree.query(target_coord)
                nearest_source_id = nearest_source_idx + 1  # Convert to 1-based ID
                target_to_source[target_node_id] = nearest_source_id
        
        return source_to_target, target_to_source

    def transfer_attributes(self, source_mesh_path, source_attr_path, target_mesh_path, output_path):
        """Transfer nodal attributes from source mesh to target mesh.
        
        Args:
            source_mesh_path (str): Path to source grid file (fort.14)
            source_attr_path (str): Path to source nodal attribute file (fort.13)
            target_mesh_path (str): Path to target grid file (fort.14)
            output_path (str): Path to write output nodal attribute file (fort.13)
        """
        print("Loading source mesh and attributes...")
        # Load source mesh and attributes
        source_mesh = AdcircMesh.open(source_mesh_path, crs=None)
        source_attrs = NodalAttributes(source_mesh)
        source_attrs.import_fort13(source_attr_path)
        
        print("Loading target mesh...")
        # Load target mesh
        target_mesh = AdcircMesh.open(target_mesh_path, crs=None)
        
        # Build conversion tables
        print("Building conversion tables...")
        source_to_target, target_to_source = self._build_conversion_tables(source_mesh, target_mesh)
        
        # Create target nodal attributes object
        target_attrs = NodalAttributes(target_mesh)
        
        # Transfer each attribute
        print("Transferring attributes...")
        for attr_name in source_attrs.get_attribute_names():
            print(f"Processing attribute: {attr_name}")
            source_attr = source_attrs.get_attribute(attr_name)
            source_values = source_attr['values']
            
            # Create target values array with same shape
            target_values = np.zeros((len(target_mesh.nodes), source_values.shape[1]))
            
            # Copy values using conversion tables
            for target_node_id in range(1, len(target_mesh.nodes) + 1):
                target_idx = target_node_id - 1  # Convert to 0-based index
                source_node_id = target_to_source[target_node_id]
                source_idx = source_node_id - 1  # Convert to 0-based index
                
                if attr_name == 'condensed_nodes':
                    # Convert node IDs in condensed_nodes attribute
                    values = source_values[source_idx].copy()
                    nonzero_count = np.sum(values > 0)
                    if nonzero_count > 0:
                        for i in range(len(values)):
                            if values[i] > 0:
                                old_id = int(values[i])
                                new_id = source_to_target[old_id]
                                values[i] = new_id
                    target_values[target_idx] = values
                else:
                    # Regular attribute, just copy the value
                    target_values[target_idx] = source_values[source_idx]
            
            # Add attribute to target
            target_attrs.add_attribute(attr_name, source_attr['units'])
            target_attrs.set_attribute(
                attr_name, 
                target_values
            )
        
        # Write output
        print(f"Writing output to {output_path}...")
        target_attrs.write(output_path, overwrite=True)
        print("Done!")

def main():
    """Command line interface for attribute transfer."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Transfer nodal attributes from one mesh to another.'
    )
    parser.add_argument(
        'source_mesh',
        help='Path to source grid file (fort.14)'
    )
    parser.add_argument(
        'source_attrs',
        help='Path to source nodal attribute file (fort.13)'
    )
    parser.add_argument(
        'target_mesh',
        help='Path to target grid file (fort.14)'
    )
    parser.add_argument(
        '-o', '--output',
        default='transferred_attributes.fort.13',
        help='Output nodal attribute file path (default: transferred_attributes.fort.13)'
    )
    
    args = parser.parse_args()
    
    transfer = AttributeTransfer()
    transfer.transfer_attributes(
        args.source_mesh,
        args.source_attrs,
        args.target_mesh,
        args.output
    )

if __name__ == '__main__':
    main() 