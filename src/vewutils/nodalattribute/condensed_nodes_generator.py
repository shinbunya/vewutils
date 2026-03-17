#!/usr/bin/env python3
"""
Create condensed_nodes nodal attributes from an ADCIRC mesh.
"""

import argparse
from typing import Dict, List, Optional, Set

import geopandas as gpd
import numpy as np
from pyproj import CRS, Geod, Transformer
from adcircpy import AdcircMesh
from adcircpy.mesh.fort13 import NodalAttributes


class CondensedNodesGenerator:
    """Generate condensed_nodes groups from VEW channel-node seeds and mesh connectivity."""

    def __init__(self,
                 mesh_file: str,
                 output_file: str,
                 distance_threshold: float,
                 epsg_from: int,
                 channel_polygon_file: Optional[str] = None,
                 max_iterations: int = 10):
        self.mesh_file = mesh_file
        self.output_file = output_file
        self.distance_threshold = distance_threshold
        self.epsg_from = epsg_from
        self.channel_polygon_file = channel_polygon_file
        self.max_iterations = max_iterations

        if self.distance_threshold <= 0.0:
            raise ValueError("distance_threshold must be positive")
        if self.epsg_from <= 0:
            raise ValueError("epsg_from must be positive")
        if self.max_iterations <= 0:
            raise ValueError("max_iterations must be positive")

        self.mesh = AdcircMesh.open(mesh_file)
        self.boundaries = self.mesh.boundaries.to_dict()
        self._mesh_crs = CRS.from_epsg(self.epsg_from)
        self._to_lonlat = Transformer.from_crs(self._mesh_crs, CRS.from_epsg(4326), always_xy=True)
        self._geod = Geod(ellps='WGS84')
        self._vew_boundary_nodes: Dict[int, Set[int]] = {}
        self._vew_bank_nodes: Set[int] = set()
        self._vew_channel_nodes: Set[int] = set()

    @staticmethod
    def _boundary_key_is_vew(boundary_key) -> bool:
        """Return True when the boundary key represents VEW (ibtype 64)."""
        return str(boundary_key) == '64'

    @staticmethod
    def _to_int(node_id) -> int:
        """Convert a mesh node identifier to int."""
        return int(node_id)

    def _load_vew_boundary_info(self) -> None:
        """Collect VEW boundary node membership information."""
        self._vew_boundary_nodes = {}
        self._vew_bank_nodes = set()
        self._vew_channel_nodes = set()

        for boundary_key, boundary_list in self.boundaries.items():
            if not self._boundary_key_is_vew(boundary_key):
                continue

            for boundary_idx, boundary in enumerate(boundary_list):
                members: Set[int] = set()
                for node_pair in boundary.get('node_id', []):
                    if not isinstance(node_pair, (list, tuple)) or len(node_pair) < 2:
                        continue
                    bank_node = self._to_int(node_pair[0])
                    channel_node = self._to_int(node_pair[1])
                    self._vew_bank_nodes.add(bank_node)
                    self._vew_channel_nodes.add(channel_node)
                    members.add(bank_node)
                    members.add(channel_node)
                if members:
                    self._vew_boundary_nodes[boundary_idx] = members

    def _load_polygon_candidate_nodes(self) -> Set[int]:
        """Load candidate nodes from channel polygons."""
        if not self.channel_polygon_file:
            return set(self.mesh.nodes.index.astype(int).tolist())

        gdf = gpd.read_file(self.channel_polygon_file)
        if gdf.empty:
            print(f"Warning: No features found in polygon file: {self.channel_polygon_file}", flush=True)
            return set()

        if gdf.crs is None:
            raise ValueError("channel polygon file must have a defined CRS")

        geom_series = gdf.geometry
        if hasattr(geom_series, 'union_all'):
            polygon_union = geom_series.union_all()
        else:
            polygon_union = geom_series.unary_union

        node_coords = self.mesh.nodes[['x', 'y']]
        points = gpd.GeoSeries(
            gpd.points_from_xy(node_coords['x'], node_coords['y']),
            index=node_coords.index,
            crs=f'EPSG:{self.epsg_from}'
        )
        if str(gdf.crs) != str(points.crs):
            points = points.to_crs(gdf.crs)
        inside_mask = points.intersects(polygon_union)
        return set(points.index[inside_mask].astype(int).tolist())

    def _get_lonlat_coords(self, node_ids: np.ndarray) -> np.ndarray:
        """Transform mesh node coordinates to lon/lat for geodesic distance checks."""
        node_xy = self.mesh.nodes.loc[node_ids, ['x', 'y']]
        lon, lat = self._to_lonlat.transform(
            node_xy['x'].to_numpy(dtype=float),
            node_xy['y'].to_numpy(dtype=float)
        )
        return np.column_stack([lon, lat])

    def _build_traversable_nodes(self) -> Set[int]:
        """Build the set of nodes eligible for connected-neighbor expansion."""
        traversable_nodes = self._load_polygon_candidate_nodes()
        traversable_nodes.difference_update(self._vew_bank_nodes)
        traversable_nodes.update(self._vew_channel_nodes)

        valid_nodes = set(self.mesh.nodes.index.astype(int).tolist())
        traversable_nodes.intersection_update(valid_nodes)
        return traversable_nodes

    def _build_vew_membership_lookup(self, node_ids: np.ndarray) -> List[Set[int]]:
        """Map nodes to VEW boundary strings that contain them."""
        node_to_index = {int(node_id): idx for idx, node_id in enumerate(node_ids.tolist())}
        membership: List[Set[int]] = [set() for _ in range(len(node_ids))]
        for boundary_idx, members in self._vew_boundary_nodes.items():
            for node_id in members:
                if node_id not in node_to_index:
                    continue
                membership[node_to_index[node_id]].add(boundary_idx)
        return membership

    def _build_groups(self, traversable_nodes: Set[int]) -> List[List[int]]:
        """Build VEW-seeded groups by iteratively expanding through connected close neighbors."""
        if not traversable_nodes:
            return []

        node_ids = np.array(sorted(traversable_nodes), dtype=int)
        node_to_index = {int(node_id): idx for idx, node_id in enumerate(node_ids.tolist())}
        print(f"Transforming traversable nodes from EPSG:{self.epsg_from} to lon/lat for geodesic distance checks...", flush=True)
        coords = self._get_lonlat_coords(node_ids)
        vew_membership = self._build_vew_membership_lookup(node_ids)
        node_neighbors = self.mesh.node_neighbors

        seed_node_ids = sorted(node_id for node_id in self._vew_channel_nodes if node_id in node_to_index)
        if not seed_node_ids:
            print("No VEW channel nodes found in the traversable node set.", flush=True)
            return []

        seed_indices = [node_to_index[node_id] for node_id in seed_node_ids]
        parent = np.arange(len(node_ids), dtype=int)
        rank = np.zeros(len(node_ids), dtype=int)
        component_vew_memberships = [set(groups) for groups in vew_membership]
        active = np.zeros(len(node_ids), dtype=bool)
        frontier = sorted(set(seed_indices))
        active[frontier] = True
        close_edges_examined = 0
        skipped_same_string_edges = 0

        def find(idx: int) -> int:
            while parent[idx] != idx:
                parent[idx] = parent[parent[idx]]
                idx = parent[idx]
            return idx

        def union_roots(root_a: int, root_b: int) -> int:
            if root_a == root_b:
                return root_a
            if rank[root_a] < rank[root_b]:
                parent[root_a] = root_b
                component_vew_memberships[root_b].update(component_vew_memberships[root_a])
                component_vew_memberships[root_a].clear()
                return root_b
            if rank[root_a] > rank[root_b]:
                parent[root_b] = root_a
                component_vew_memberships[root_a].update(component_vew_memberships[root_b])
                component_vew_memberships[root_b].clear()
                return root_a
            parent[root_b] = root_a
            component_vew_memberships[root_a].update(component_vew_memberships[root_b])
            component_vew_memberships[root_b].clear()
            rank[root_a] += 1
            return root_a

        print(f"Starting from {len(frontier)} VEW channel seed nodes.", flush=True)
        print(f"Traversable nodes available: {len(traversable_nodes)}", flush=True)

        for iteration in range(1, self.max_iterations + 1):
            if not frontier:
                print(f"No new nodes found at iteration {iteration - 1}; stopping expansion.", flush=True)
                break

            print(f"Iteration {iteration}: expanding {len(frontier)} frontier nodes...", flush=True)
            next_frontier: Set[int] = set()

            for current_idx in frontier:
                current_node_id = int(node_ids[current_idx])
                lon0 = coords[current_idx, 0]
                lat0 = coords[current_idx, 1]

                for neighbor_node_id in node_neighbors.get(current_node_id, []):
                    neighbor_node_id = int(neighbor_node_id)
                    if neighbor_node_id not in node_to_index:
                        continue

                    neighbor_idx = node_to_index[neighbor_node_id]
                    lon1 = coords[neighbor_idx, 0]
                    lat1 = coords[neighbor_idx, 1]
                    _, _, distance_m = self._geod.inv(lon0, lat0, lon1, lat1)
                    if distance_m > self.distance_threshold:
                        continue

                    close_edges_examined += 1
                    root_current = find(current_idx)
                    root_neighbor = find(neighbor_idx)
                    if root_current == root_neighbor:
                        continue

                    if component_vew_memberships[root_current].intersection(component_vew_memberships[root_neighbor]):
                        skipped_same_string_edges += 1
                        continue

                    union_roots(root_current, root_neighbor)
                    if not active[neighbor_idx]:
                        active[neighbor_idx] = True
                        next_frontier.add(neighbor_idx)

            frontier = sorted(next_frontier)
        else:
            if frontier:
                print(f"Reached max iterations ({self.max_iterations}).", flush=True)

        group_map: Dict[int, List[int]] = {}
        for idx, node_id in enumerate(node_ids.tolist()):
            if not active[idx]:
                continue
            root_idx = find(idx)
            group_map.setdefault(root_idx, []).append(int(node_id))

        groups = [sorted(group) for group in group_map.values() if len(group) >= 2]

        print(f"Activated nodes: {int(np.count_nonzero(active))}", flush=True)
        print(f"Connected close edges examined: {close_edges_examined}", flush=True)
        print(f"Edges skipped by VEW-string rule: {skipped_same_string_edges}", flush=True)
        print(f"Condensed-node groups created: {len(groups)}", flush=True)
        return groups

    @staticmethod
    def _groups_to_attribute_values(node_count: int, groups: List[List[int]]) -> np.ndarray:
        """Convert groups to dense condensed_nodes attribute values."""
        max_group_size = max((len(group) for group in groups), default=1)
        values = np.zeros((node_count, max_group_size - 1 if max_group_size > 1 else 1), dtype=int)

        for group in groups:
            representative = int(group[0])
            if representative < 1 or representative > node_count:
                continue
            if len(group) > 1:
                values[representative - 1, :len(group) - 1] = np.array(group[1:], dtype=int)

        return values

    def write_condensed_nodes_fort13(self) -> None:
        """Generate and write the condensed_nodes fort.13 file."""
        self._load_vew_boundary_info()
        traversable_nodes = self._build_traversable_nodes()
        groups = self._build_groups(traversable_nodes)
        values = self._groups_to_attribute_values(len(self.mesh.nodes), groups)

        fort13 = NodalAttributes(self.mesh)
        fort13.add_attribute('condensed_nodes', 'unitless')
        fort13.set_attribute('condensed_nodes', values)
        fort13.write(self.output_file, overwrite=True)

        print(f"VEW bank nodes excluded: {len(self._vew_bank_nodes)}", flush=True)
        print(f"VEW channel nodes used as seeds: {len(self._vew_channel_nodes)}", flush=True)
        print(f"Wrote condensed_nodes fort.13 to: {self.output_file}", flush=True)


def get_parser():
    """Create argument parser for condensed_nodes generation."""
    parser = argparse.ArgumentParser(
        add_help=False,
        description='Generate condensed_nodes fort.13 attributes from an ADCIRC mesh.'
    )
    parser.add_argument(
        'mesh',
        help='Path to the input mesh file (fort.14)'
    )
    parser.add_argument(
        '-o', '--output',
        default='condensed_nodes.fort.13',
        help='Path to the output fort.13 file (default: condensed_nodes.fort.13)'
    )
    parser.add_argument(
        '-t', '--distance-threshold',
        type=float,
        required=True,
        help='Maximum distance between connected nodes in the same condensed_nodes group'
    )
    parser.add_argument(
        '--epsg-from',
        type=int,
        required=True,
        help='EPSG code of the input mesh coordinates'
    )
    parser.add_argument(
        '-p', '--channel-polygons',
        help='Optional polygon file (GeoJSON/Shapefile) used to limit traversable nodes'
    )
    parser.add_argument(
        '--max-iterations',
        type=int,
        default=10,
        help='Maximum connected-neighbor expansion iterations (default: 10)'
    )
    return parser


def main(args=None):
    """Run condensed_nodes generation."""
    if args is None:
        args = get_parser().parse_args()

    generator = CondensedNodesGenerator(
        mesh_file=args.mesh,
        output_file=args.output,
        distance_threshold=args.distance_threshold,
        epsg_from=args.epsg_from,
        channel_polygon_file=args.channel_polygons,
        max_iterations=args.max_iterations
    )
    generator.write_condensed_nodes_fort13()


if __name__ == '__main__':
    main()
