#!/usr/bin/env python3
import argparse
import os
import tempfile

import numpy as np
import pandas as pd
from adcircpy import AdcircMesh
from adcircpy.mesh.fort13 import NodalAttributes, parse_fort13


def strip_fortran_comments(line: str) -> str:
    """Remove Fortran-style trailing comments (``!``) from a fort.13 line."""
    return line.split('!', 1)[0].rstrip()


def fort13_content_without_comments(fort13_file: str) -> str:
    """Read a fort.13 file with Fortran ``!`` comments stripped from each line."""
    with open(fort13_file, 'r') as f:
        return '\n'.join(strip_fortran_comments(line) for line in f)


def _parse_fort13_from_cleaned_content(content: str):
    """Parse fort.13 content via a temporary file (adcircpy expects a path)."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.13', delete=False) as tmp:
        tmp.write(content)
        tmp_path = tmp.name
    try:
        return parse_fort13(tmp_path)
    finally:
        os.unlink(tmp_path)


def import_fort13_strip_comments(fort13: NodalAttributes, fort13_file: str) -> None:
    """Import a fort.13 file after stripping Fortran ``!`` comments."""
    content = fort13_content_without_comments(fort13_file)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.13', delete=False) as tmp:
        tmp.write(content)
        tmp_path = tmp.name
    try:
        fort13.import_fort13(tmp_path)
    finally:
        os.unlink(tmp_path)


def parse_fort13_strip_comments(fort13_file: str):
    """Parse a fort.13 file after stripping Fortran ``!`` comments."""
    return _parse_fort13_from_cleaned_content(fort13_content_without_comments(fort13_file))


def invalidate_adcircpy_nodal_attribute_cache(fort13, attribute_name):
    """Clear adcircpy-internal caches so fort.13 output reflects updated values.

    NodalAttributes.get_attribute() computes and stores ``defaults`` and
    ``non_default_indexes``; set_attribute() updates ``values`` but does not
    invalidate them, so write()/__str__() can omit nodes that newly differ
    from defaults (or keep stale rows).
    """
    attr = fort13._attributes.get(attribute_name)
    if attr is None:
        return
    attr.pop("defaults", None)
    attr.pop("non_default_indexes", None)


class NodalAttributeSetter:
    """Set a single nodal attribute value at selected nodes in a fort.13 file."""

    def __init__(
        self,
        mesh_file,
        fort13_file,
        output_file,
        attribute_name,
        value,
        selected_nodes_file,
    ):
        self.mesh_file = mesh_file
        self.fort13_file = fort13_file
        self.output_file = output_file
        self.attribute_name = attribute_name
        self.value = value
        self.selected_nodes_file = selected_nodes_file

        self.mesh = AdcircMesh.open(mesh_file)
        self.fort13 = NodalAttributes(self.mesh)
        self.fort13.import_fort13(fort13_file)
        self.selected_nodes = self._load_selected_nodes()

    def _load_selected_nodes(self):
        try:
            df = pd.read_csv(self.selected_nodes_file)
            if "node_id" in df.columns:
                selected_nodes = df["node_id"].values
            else:
                selected_nodes = df.iloc[:, 0].values
        except Exception:
            try:
                selected_nodes = np.loadtxt(self.selected_nodes_file, dtype=int)
            except Exception as exc:
                raise ValueError(
                    f"Could not read selected nodes from {self.selected_nodes_file}. "
                    "File should be either a CSV with 'node_id' column or a text file "
                    "with one node ID per line."
                ) from exc

        if selected_nodes.ndim > 1:
            selected_nodes = selected_nodes.flatten()

        if len(selected_nodes) == 0:
            raise ValueError(f"No nodes found in {self.selected_nodes_file}")

        max_node_id = self.mesh.nodes.shape[0]
        invalid = set(selected_nodes[selected_nodes < 1]) | set(
            selected_nodes[selected_nodes > max_node_id]
        )
        if invalid:
            raise ValueError(
                f"The following node IDs in {self.selected_nodes_file} are invalid "
                f"(must be between 1 and {max_node_id}):\n{invalid}"
            )

        print(f"Loaded {len(selected_nodes)} selected nodes")
        return selected_nodes

    def set_values(self):
        try:
            attr = self.fort13.get_attribute(self.attribute_name)
        except (KeyError, AttributeError) as exc:
            raise ValueError(
                f"Nodal attribute '{self.attribute_name}' was not found in {self.fort13_file}."
            ) from exc

        existing = np.asarray(attr["values"])
        if existing.ndim == 1:
            existing = existing.reshape(-1, 1)
        n_nodes, ncols = existing.shape

        if n_nodes != self.mesh.nodes.shape[0]:
            raise ValueError(
                f"Attribute '{self.attribute_name}' has {n_nodes} rows but mesh has "
                f"{self.mesh.nodes.shape[0]} nodes."
            )
        if ncols != 1:
            raise ValueError(
                f"Attribute '{self.attribute_name}' has {ncols} values per node; "
                "this tool currently supports attributes with exactly one value per node."
            )

        new_vals = existing.copy()
        scalar = np.asarray(self.value, dtype=new_vals.dtype)
        selected_0 = self.selected_nodes - 1
        new_vals[selected_0, 0] = scalar

        print(
            f"Setting {self.attribute_name} = {scalar} at {len(self.selected_nodes)} nodes"
        )

        self.fort13.set_attribute(self.attribute_name, new_vals)
        invalidate_adcircpy_nodal_attribute_cache(self.fort13, self.attribute_name)
        self.fort13.write(self.output_file, overwrite=True)
        print(f"Updated fort.13 file written to: {self.output_file}")


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description=(
            "Set a single nodal attribute value at selected nodes in a fort.13 file "
            "(attributes with one value per node only)."
        ),
    )
    parser.add_argument("mesh", help="Path to the mesh file (fort.14)")
    parser.add_argument("fort13", help="Path to the input nodal attribute file (fort.13)")
    parser.add_argument("output", help="Path to the output fort.13 file")
    parser.add_argument(
        "attribute_name",
        help="Nodal attribute name as in fort.13 (e.g. mannings_n_at_sea_floor)",
    )
    parser.add_argument(
        "value",
        type=float,
        help="Scalar value to assign at each selected node",
    )
    parser.add_argument(
        "--selected-nodes",
        required=True,
        help=(
            "Path to text file containing selected node IDs "
            "(one per line or CSV with node_id column)"
        ),
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()

    setter = NodalAttributeSetter(
        args.mesh,
        args.fort13,
        args.output,
        args.attribute_name,
        args.value,
        args.selected_nodes,
    )
    setter.set_values()


if __name__ == "__main__":
    main()
