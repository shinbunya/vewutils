#!/usr/bin/env python3
import numpy as np
import pandas as pd
import argparse
from adcircpy import AdcircMesh
from adcircpy.mesh.fort13 import NodalAttributes

class ManningsnSetter:
    """Class for setting Manning's n values at specific nodes in ADCIRC fort.13 files."""
    
    def __init__(self, mesh_file, fort13_file, output_file, mannings_value, selected_nodes_file):
        """Initialize the ManningsnSetter.
        
        Args:
            mesh_file (str): Path to the mesh file (fort.14)
            fort13_file (str): Path to the input nodal attribute file (fort.13)
            output_file (str): Path to the output fort.13 file
            mannings_value (float): Manning's n value to set
            selected_nodes_file (str): Path to text file containing selected node IDs
        """
        self.mesh_file = mesh_file
        self.fort13_file = fort13_file
        self.output_file = output_file
        self.mannings_value = float(mannings_value)
        self.selected_nodes_file = selected_nodes_file
        
        # Load mesh
        self.mesh = AdcircMesh.open(mesh_file)
        
        # Load nodal attributes
        self.fort13 = NodalAttributes(self.mesh)
        self.fort13.import_fort13(fort13_file)
        
        # Load selected nodes
        self.selected_nodes = self._load_selected_nodes()
        
    def _load_selected_nodes(self):
        """Load selected nodes from file.
        
        Returns:
            np.ndarray: Array of selected node IDs (1-based)
        """
        try:
            # Try reading as CSV first (if it has a header)
            df = pd.read_csv(self.selected_nodes_file)
            if 'node_id' in df.columns:
                selected_nodes = df['node_id'].values
            else:
                # If no 'node_id' column, assume first column contains node IDs
                selected_nodes = df.iloc[:, 0].values
        except:
            # If CSV reading fails, try reading as plain text file
            try:
                selected_nodes = np.loadtxt(self.selected_nodes_file, dtype=int)
            except:
                raise ValueError(f"Could not read selected nodes from {self.selected_nodes_file}. "
                               "File should be either a CSV with 'node_id' column or a text file "
                               "with one node ID per line.")
        
        # Ensure it's a 1D array
        if selected_nodes.ndim > 1:
            selected_nodes = selected_nodes.flatten()
            
        # Check if any nodes were specified
        if len(selected_nodes) == 0:
            raise ValueError(f"No nodes found in {self.selected_nodes_file}")
        
        # Validate node IDs
        max_node_id = self.mesh.nodes.shape[0]
        invalid_nodes = set(selected_nodes[selected_nodes < 1]) | set(selected_nodes[selected_nodes > max_node_id])
        if len(invalid_nodes) > 0:
            raise ValueError(
                f"The following node IDs in {self.selected_nodes_file} are invalid "
                f"(must be between 1 and {max_node_id}):\n{invalid_nodes}"
            )
            
        print(f"Loaded {len(selected_nodes)} selected nodes")
        return selected_nodes
    
    def set_mannings_values(self):
        """Set Manning's n values at the selected nodes."""
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
        
        # Convert selected nodes to 0-based indexing
        selected_nodes_0based = self.selected_nodes - 1
        
        # Set Manning's n values at selected nodes
        new_mn[selected_nodes_0based] = self.mannings_value
        
        print(f"Setting Manning's n = {self.mannings_value} at {len(self.selected_nodes)} nodes")
        
        # Update the fort.13 file
        # Reshape the array to 2D if needed
        new_mn_2d = new_mn.reshape(-1, 1)
        
        # Set the attribute
        self.fort13.set_attribute('mannings_n_at_sea_floor', new_mn_2d)
        
        # Write the output file
        self.fort13.write(self.output_file, overwrite=True)
        
        print(f"Updated fort.13 file written to: {self.output_file}")

def get_parser():
    """Get argument parser for the Manning's n setter."""
    parser = argparse.ArgumentParser(
        add_help=False, 
        description="Set Manning's n values at specific nodes in a fort.13 file"
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
    parser.add_argument(
        'mannings_value',
        type=float,
        help="Manning's n value to set at selected nodes"
    )
    parser.add_argument(
        '--selected-nodes',
        required=True,
        help='Path to text file containing selected node IDs (one per line or CSV with node_id column)'
    )
    return parser

def main(args=None):
    """Main function for the Manning's n setter."""
    if args is None:
        args = get_parser().parse_args()
    
    setter = ManningsnSetter(
        args.mesh,
        args.fort13,
        args.output,
        args.mannings_value,
        args.selected_nodes
    )
    
    setter.set_mannings_values()

if __name__ == '__main__':
    main() 