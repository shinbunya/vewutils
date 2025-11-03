#!/usr/bin/env python3
"""
Script to copy depths (value_1) at selected nodes from a source mesh to a base mesh,
then ensure VEW boundary barrier heights are above bank elevations before writing.
"""

import argparse
import numpy as np
import pandas as pd
from adcircpy import AdcircMesh
from vewutils.mesh.vew_boundary_manipulator import VEWBoundaryManipulator


class SelectedNodeDepthCopier:
    """Copy node depths from a source mesh to a base mesh for selected nodes."""

    def __init__(self, base_mesh_path: str, source_mesh_path: str, selected_nodes_path: str):
        self.base_mesh_path = base_mesh_path
        self.source_mesh_path = source_mesh_path
        self.selected_nodes_path = selected_nodes_path

        # Load meshes
        self.base_mesh = AdcircMesh.open(base_mesh_path)
        self.source_mesh = AdcircMesh.open(source_mesh_path)

        # Load selected nodes
        self.selected_nodes = self._load_selected_nodes()

        # Validate meshes have same or compatible node indexing for selected nodes
        max_base_node_id = self.base_mesh.nodes.shape[0]
        max_source_node_id = self.source_mesh.nodes.shape[0]
        invalid_in_base = set(self.selected_nodes[self.selected_nodes < 1]) | set(
            self.selected_nodes[self.selected_nodes > max_base_node_id]
        )
        invalid_in_source = set(self.selected_nodes[self.selected_nodes < 1]) | set(
            self.selected_nodes[self.selected_nodes > max_source_node_id]
        )
        if invalid_in_base:
            raise ValueError(
                f"Selected nodes contain IDs not present in base mesh (1..{max_base_node_id}): {invalid_in_base}"
            )
        if invalid_in_source:
            raise ValueError(
                f"Selected nodes contain IDs not present in source mesh (1..{max_source_node_id}): {invalid_in_source}"
            )

    def _load_selected_nodes(self) -> np.ndarray:
        """Load selected nodes from CSV or text file (1-based node IDs)."""
        try:
            df = pd.read_csv(self.selected_nodes_path)
            if 'node_id' in df.columns:
                selected_nodes = df['node_id'].values
            else:
                selected_nodes = df.iloc[:, 0].values
        except Exception:
            try:
                selected_nodes = np.loadtxt(self.selected_nodes_path, dtype=int)
            except Exception:
                raise ValueError(
                    f"Could not read selected nodes from {self.selected_nodes_path}. "
                    "Provide a CSV with 'node_id' column or a text file with one node ID per line."
                )

        if isinstance(selected_nodes, np.ndarray) and selected_nodes.ndim > 1:
            selected_nodes = selected_nodes.flatten()

        if len(selected_nodes) == 0:
            raise ValueError(f"No nodes found in {self.selected_nodes_path}")

        print(f"Loaded {len(selected_nodes)} selected nodes")
        return selected_nodes.astype(int)

    def copy_depths(self):
        """Copy value_1 from source mesh to base mesh at selected nodes."""
        # Convert to 0-based indices to use DataFrame iloc if desired, but here we use index labels (1-based)
        for node_id in self.selected_nodes:
            self.base_mesh.nodes.loc[int(node_id), 'value_1'] = self.source_mesh.nodes.loc[int(node_id), 'value_1']


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description=(
            "Copy depths (value_1) at selected nodes from source mesh to base mesh "
            "and ensure VEW barrier heights are above banks before writing."
        ),
    )
    parser.add_argument('base_mesh', help='Path to the base mesh file (fort.14)')
    parser.add_argument('source_mesh', help='Path to the source mesh file (fort.14)')
    parser.add_argument('--selected-nodes', required=True, help="Path to selected node IDs file (CSV with 'node_id' or text, 1 per line)")
    parser.add_argument('-o', '--output', default='updated_mesh.14', help='Output mesh file path (default: updated_mesh.14)')
    parser.add_argument('-t', '--tolerance', type=float, default=0.001, help='Barrier height tolerance above bank elevations (meters, default: 0.001)')
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()

    copier = SelectedNodeDepthCopier(args.base_mesh, args.source_mesh, args.selected_nodes)
    copier.copy_depths()

    # Ensure VEW barrier heights are above banks before writing, following adjust_vew_barrier_heights.py
    mesh_out = VEWBoundaryManipulator.ensure_barrier_heights_above_banks(copier.base_mesh, args.tolerance)
    mesh_out.write(args.output, overwrite=True)
    print(f"Updated mesh saved to: {args.output}")


if __name__ == '__main__':
    main()


