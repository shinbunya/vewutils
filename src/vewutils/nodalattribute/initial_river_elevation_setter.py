#!/usr/bin/env python3
"""
Set the initial_river_elevation nodal attribute in an ADCIRC fort.13 file.
"""

import argparse
from typing import Dict, Optional, Set, Tuple

import numpy as np
from adcircpy import AdcircMesh
from adcircpy.mesh.fort13 import NodalAttributes

from vewutils.nodalattribute.nodalattribute_setter import (
    import_fort13_strip_comments,
    invalidate_adcircpy_nodal_attribute_cache,
    parse_fort13_strip_comments,
)

ATTRIBUTE_NAME = 'initial_river_elevation'
ATTRIBUTE_UNITS = 'm'


class InitialRiverElevationSetter:
    """Compute and write initial_river_elevation nodal attribute values.

    For each node i the initial river elevation ``zeta_ini_i`` is determined by:
        1) If node i already has a non-default value in the input fort.13 file,
           keep that value (the rules below are not applied).
        2) Otherwise, only for channel-side nodes of VEW boundaries (ibtype 64),
           ``zeta_ini_i = depth_over_bed - b_i`` where ``b_i`` is the nodal depth
           from the mesh (fort.14).
        3) Then, for those channel-side nodes, clamp from below:
           ``zeta_ini_i = max(zeta_ini_i, min_elev)``.

    Nodes that are neither overridden (rule 1) nor channel-side VEW nodes keep
    the attribute's default value.
    """

    def __init__(self,
                 mesh_file: str,
                 output_file: str,
                 depth_over_bed: float,
                 min_elev: float,
                 fort13_file: Optional[str] = None):
        self.mesh_file = mesh_file
        self.output_file = output_file
        self.depth_over_bed = float(depth_over_bed)
        self.min_elev = float(min_elev)
        self.fort13_file = fort13_file

        self.mesh = AdcircMesh.open(mesh_file)
        self.fort13 = NodalAttributes(self.mesh)
        if self.fort13_file:
            import_fort13_strip_comments(self.fort13, self.fort13_file)

    def _get_nodal_depths(self) -> np.ndarray:
        """Return nodal depths b_i (fort.14 depth, positive downward) in node order.

        adcircpy negates the fort.14 depth column on read (``Fort14.open``), so
        ``mesh.values`` holds bed elevation (positive up). The raw fort.14 nodal
        depth b_i is therefore the negative of that.
        """
        return -np.asarray(self.mesh.values, dtype=float).flatten()

    def _get_vew_channel_nodes(self) -> Set[int]:
        """Return the set of 1-based channel-side node IDs from VEW boundaries.

        VEW boundaries use ibtype 64 and store node pairs as
        ``[bank_node, channel_node]``; the channel-side node is the second entry.
        """
        channel_nodes: Set[int] = set()
        boundaries = self.mesh.boundaries.to_dict()
        for vew_boundary in boundaries.get('64', []):
            for node_pair in vew_boundary.get('node_id', []):
                if not isinstance(node_pair, (list, tuple)) or len(node_pair) < 2:
                    continue
                channel_nodes.add(int(node_pair[1]))
        return channel_nodes

    def _read_input_overrides(self) -> Tuple[Dict[int, float], Optional[float]]:
        """Read explicitly listed (non-default) initial_river_elevation values.

        The values are taken directly from the input fort.13 per-node section so
        that nodes carrying an explicit (non-default) value are preserved exactly.
        Returns a tuple ``(overrides, default)`` where ``overrides`` maps 1-based
        node IDs to values and ``default`` is the attribute's default value in the
        input file (``None`` when the attribute is absent or no input was given).
        """
        if not self.fort13_file:
            return {}, None

        parsed = parse_fort13_strip_comments(self.fort13_file)
        if ATTRIBUTE_NAME not in parsed:
            print(
                f"Note: '{ATTRIBUTE_NAME}' not found in {self.fort13_file}; "
                "channel-side nodes will be computed from the rules and all "
                "other nodes set to the default value.",
                flush=True,
            )
            return {}, None

        attribute = parsed[ATTRIBUTE_NAME]
        node_ids = attribute['node_id']
        values = np.asarray(attribute['values'], dtype=float)
        if values.ndim == 1:
            values = values.reshape((values.shape[0], 1))
        default = float(np.asarray(attribute['defaults'], dtype=float).flatten()[0])

        overrides: Dict[int, float] = {}
        for i, node_id in enumerate(node_ids):
            overrides[int(node_id)] = float(values[i, 0])

        print(f"Preserving {len(overrides)} non-default node(s) from {self.fort13_file}.", flush=True)
        return overrides, default

    def _compute_values(self) -> np.ndarray:
        """Compute the dense initial_river_elevation values for every node."""
        depths = self._get_nodal_depths()
        node_count = depths.shape[0]

        overrides, input_default = self._read_input_overrides()

        # Nodes that are neither overridden nor channel-side VEW nodes keep the
        # attribute's default value (the input file's default when available).
        default_value = input_default if input_default is not None else 0.0
        zeta = np.full(node_count, default_value, dtype=float)

        # Rules 2 and 3 apply only to channel-side nodes of VEW boundaries.
        channel_nodes = self._get_vew_channel_nodes()
        channel_idx = np.array(
            sorted(node_id - 1 for node_id in channel_nodes if 1 <= node_id <= node_count),
            dtype=int,
        )
        if channel_idx.size:
            # Rule 2: zeta_ini_i = depth_over_bed - b_i
            channel_zeta = self.depth_over_bed - depths[channel_idx]
            # Rule 3: zeta_ini_i = max(zeta_ini_i, min_elev)
            channel_zeta = np.maximum(channel_zeta, self.min_elev)
            zeta[channel_idx] = channel_zeta
        else:
            print(
                "Warning: no VEW channel-side nodes found; rules 2 and 3 were "
                "not applied to any node.",
                flush=True,
            )

        # Rule 1: keep non-default values already set in the input fort.13 file.
        out_of_range = 0
        for node_id, value in overrides.items():
            if node_id < 1 or node_id > node_count:
                out_of_range += 1
                continue
            zeta[node_id - 1] = value
        if out_of_range:
            print(
                f"Warning: {out_of_range} override node ID(s) outside mesh range were ignored.",
                flush=True,
            )

        print(f"Computed initial_river_elevation for {node_count} nodes.", flush=True)
        print(f"  VEW channel-side nodes updated = {int(channel_idx.size)}", flush=True)
        print(f"  default value (other nodes)    = {default_value}", flush=True)
        print(f"  depth_over_bed = {self.depth_over_bed}", flush=True)
        print(f"  min_elev       = {self.min_elev}", flush=True)
        print(f"  value range    = [{zeta.min():.6f}, {zeta.max():.6f}]", flush=True)
        return zeta.reshape((node_count, 1))

    def write_initial_river_elevation_fort13(self) -> None:
        """Compute and write the initial_river_elevation fort.13 file."""
        values = self._compute_values()

        if ATTRIBUTE_NAME in self.fort13.get_attribute_names():
            self.fort13.set_attribute(ATTRIBUTE_NAME, values)
            invalidate_adcircpy_nodal_attribute_cache(self.fort13, ATTRIBUTE_NAME)
        else:
            self.fort13.add_attribute(ATTRIBUTE_NAME, ATTRIBUTE_UNITS)
            self.fort13.set_attribute(ATTRIBUTE_NAME, values)

        self.fort13.write(self.output_file, overwrite=True)
        print(f"Wrote {ATTRIBUTE_NAME} fort.13 to: {self.output_file}", flush=True)


def get_parser():
    """Create argument parser for initial_river_elevation setting."""
    parser = argparse.ArgumentParser(
        add_help=False,
        description='Set the initial_river_elevation nodal attribute in a fort.13 file.'
    )
    parser.add_argument(
        'mesh',
        help='Path to the input mesh file (fort.14)'
    )
    parser.add_argument(
        '-i', '--input',
        dest='fort13',
        help=(
            'Optional input fort.13 file. When given, its initial_river_elevation '
            'attribute is updated (non-default nodes preserved) and all other '
            'attributes are carried over to the output. When omitted, the output '
            'contains only the initial_river_elevation attribute.'
        )
    )
    parser.add_argument(
        '-o', '--output',
        default='initial_river_elevation.fort.13',
        help='Path to the output fort.13 file (default: initial_river_elevation.fort.13)'
    )
    parser.add_argument(
        '--depth-over-bed',
        type=float,
        required=True,
        help='Water depth over the bed used in zeta_ini_i = depth_over_bed - b_i'
    )
    parser.add_argument(
        '--min-elev',
        type=float,
        required=True,
        help='Minimum allowed elevation used in zeta_ini_i = max(zeta_ini_i, min_elev)'
    )
    return parser


def main(args=None):
    """Run initial_river_elevation setting."""
    if args is None:
        args = get_parser().parse_args()

    setter = InitialRiverElevationSetter(
        mesh_file=args.mesh,
        output_file=args.output,
        depth_over_bed=args.depth_over_bed,
        min_elev=args.min_elev,
        fort13_file=args.fort13,
    )
    setter.write_initial_river_elevation_fort13()


if __name__ == '__main__':
    main()
