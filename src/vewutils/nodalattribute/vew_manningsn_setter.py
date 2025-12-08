#!/usr/bin/env python3
import numpy as np
import argparse
from adcircpy import AdcircMesh
from adcircpy.mesh.fort13 import NodalAttributes

class VEWManningsnSetter:
    """Class for copying Manning's n values from channel nodes to bank nodes along VEW boundaries."""
    
    def __init__(self, mesh_file, fort13_file, output_file):
        """Initialize the VEWManningsnSetter.
        
        Args:
            mesh_file (str): Path to the mesh file (fort.14)
            fort13_file (str): Path to the input nodal attribute file (fort.13)
            output_file (str): Path to the output fort.13 file
        """
        self.mesh_file = mesh_file
        self.fort13_file = fort13_file
        self.output_file = output_file
        
        # Load mesh
        self.mesh = AdcircMesh.open(mesh_file)
        
        # Load nodal attributes
        self.fort13 = NodalAttributes(self.mesh)
        self.fort13.import_fort13(fort13_file)
        
        # Get VEW node pairs
        self.channel_to_bank = self._get_vew_node_pairs()
    
    def _get_vew_node_pairs(self):
        """Extract channel to bank node pairs from VEW boundaries.
        
        Returns:
            dict: Dictionary mapping channel node IDs to bank node IDs
        """
        channel_to_bank = {}
        boundaries = self.mesh.boundaries.to_dict()
        
        if '64' in boundaries:  # Check for VEW boundaries (ibtype 64)
            for vew_boundary in boundaries['64']:
                for node_pair in vew_boundary['node_id']:
                    bank_node = int(node_pair[0])
                    channel_node = int(node_pair[1])
                    channel_to_bank[channel_node] = bank_node
        
        if len(channel_to_bank) == 0:
            print("Warning: No VEW boundaries found in mesh")
        else:
            print(f"Found {len(channel_to_bank)} VEW node pairs")
        
        return channel_to_bank
    
    def copy_mannings_values(self):
        """Copy Manning's n values from channel nodes to their corresponding bank nodes."""
        try:
            # Get existing Manning's n values
            existing_mn = self.fort13.get_attribute('mannings_n_at_sea_floor')['values']
            if existing_mn.ndim > 1:
                existing_mn = existing_mn.flatten()
            
            # Create a copy to modify
            new_mn = existing_mn.copy()
            
        except (KeyError, AttributeError):
            print("Warning: 'mannings_n_at_sea_floor' attribute not found in fort.13 file.")
            print("Creating new Manning's n attribute with default value 0.025")
            # Create new array with default Manning's n value
            new_mn = np.full(self.mesh.nodes.shape[0], 0.025)
        
        # Count number of bank nodes updated
        updated_count = 0
        
        # Copy Manning's n from channel nodes to bank nodes
        for channel_node, bank_node in self.channel_to_bank.items():
            # Convert to 0-based indexing
            channel_node_0based = channel_node - 1
            bank_node_0based = bank_node - 1
            
            # Copy Manning's n value from channel to bank
            channel_mn = new_mn[channel_node_0based]
            new_mn[bank_node_0based] = channel_mn
            updated_count += 1
        
        print(f"Copied Manning's n values from {updated_count} channel nodes to their corresponding bank nodes")
        
        # Update the fort.13 file
        # Reshape the array to 2D if needed
        new_mn_2d = new_mn.reshape(-1, 1)
        
        # Set the attribute
        self.fort13.set_attribute('mannings_n_at_sea_floor', new_mn_2d)
        
        # Write the output file
        self.fort13.write(self.output_file, overwrite=True)
        
        print(f"Updated fort.13 file written to: {self.output_file}")

def get_parser():
    """Get argument parser for the VEW Manning's n setter."""
    parser = argparse.ArgumentParser(
        add_help=False, 
        description="Copy Manning's n values from channel nodes to bank nodes along VEW boundaries"
    )
    parser.add_argument(
        'mesh',
        help='Path to the mesh file (fort.14)'
    )
    parser.add_argument(
        'fort13',
        help='Path to the input nodal attribute file (fort.13)'
    )
    parser.add_argument(
        'output',
        help='Path to the output fort.13 file'
    )
    return parser

def main(args=None):
    """Main function for the VEW Manning's n setter."""
    if args is None:
        args = get_parser().parse_args()
    
    setter = VEWManningsnSetter(
        args.mesh,
        args.fort13,
        args.output
    )
    
    setter.copy_mannings_values()

if __name__ == '__main__':
    main()

