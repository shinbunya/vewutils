#!/usr/bin/env python3
"""
Channel Mesh Generator - Python port (subset) of MATLAB msh_pave_channels.m

Scope: generate a channel mesh only from centerline polylines (GeoJSON/Shapefile)
with per-point or per-line widths and depths. No base-mesh subsetting, no gap-fill,
no merging back. Output is a standalone ADCIRC mesh of bank-to-bank triangular
strips along each provided centerline.
"""

import argparse
import json
import math
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.spatial import cKDTree
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
from adcircpy import AdcircMesh


@dataclass
class CenterlineSequence:
    """Represents a processed centerline sequence in lon/lat space."""
    lonlat: np.ndarray  # [N, 2] array of (lon, lat)
    depths: np.ndarray  # [N] array; interpretation depends on config
    widths: np.ndarray  # [N] array of widths [m]
    original_id: int = 0


@dataclass
class ChannelStrip:
    """Bank-to-bank strip geometry sampled at uniform spacing in lon/lat."""
    left_lonlat: np.ndarray   # [N, 2]
    right_lonlat: np.ndarray  # [N, 2]
    center_depths: np.ndarray # [N] centerline-aligned depths/elevations
    # Optional global IDs for node reuse (only endpoints typically populated)
    left_gid: Optional[np.ndarray] = None  # [N] or None
    right_gid: Optional[np.ndarray] = None # [N] or None
    # Optional channel mid nodes per cross-section (lists aligned with i)
    mid_lonlat: Optional[list] = None
    mid_depths: Optional[list] = None
    mid_gids: Optional[list] = None


class FlowlineReader:
    """Reads and processes flowline data from GeoJSON or Shapefile."""

    def __init__(self):
        self.centerlines = []

    def read_flowline_file(self, filename: str, zshift: float = 0.0,
                          zshift_taper_elev_range: Optional[Tuple[float, float]] = None,
                          effective_width_ratio: float = 1.0,
                          max_width: float = np.inf,
                          min_width: float = 0.0) -> List[CenterlineSequence]:
        """
        Read flowline data from file.

        Args:
            filename: Path to GeoJSON or Shapefile
            zshift: Vertical shift to apply to depths
            zshift_taper_elev_range: Range for tapering zshift [min_elev, max_elev]
            effective_width_ratio: Ratio to apply to widths
            max_width: Maximum allowed width
            min_width: Minimum allowed width

        Returns:
            List of CenterlineSequence objects
        """
        if filename.lower().endswith('.geojson'):
            return self._read_geojson(filename, zshift, zshift_taper_elev_range,
                                    effective_width_ratio, max_width, min_width)
        else:
            return self._read_shapefile(filename, zshift, zshift_taper_elev_range,
                                      effective_width_ratio, max_width, min_width)

    def _read_geojson(self, filename: str, zshift: float,
                     zshift_taper_elev_range: Optional[Tuple[float, float]],
                     effective_width_ratio: float, max_width: float,
                     min_width: float) -> List[CenterlineSequence]:
        """Read flowline data from GeoJSON file."""
        with open(filename, 'r') as f:
            data = json.load(f)

        centerlines = []
        features = data.get('features', [])

        for i, feature in enumerate(features):
            geometry = feature.get('geometry', {})
            properties = feature.get('properties', {})

            gtype = geometry.get('type')
            if gtype not in ('LineString', 'MultiLineString'):
                continue

            parts = []
            if gtype == 'LineString':
                coords = geometry.get('coordinates', []) or []
                if len(coords) >= 2:
                    parts.append(coords)
            else:
                for part in (geometry.get('coordinates', []) or []):
                    if part and len(part) >= 2:
                        parts.append(part)

            for coords in parts:
                x_coords = [float(c[0]) for c in coords]
                y_coords = [float(c[1]) for c in coords]

            depths = self._extract_depth_data(properties, len(x_coords), zshift,
                                            zshift_taper_elev_range)
            widths = self._extract_width_data(properties, len(x_coords),
                                            effective_width_ratio, max_width, min_width)

            centerline = CenterlineSequence(
                    lonlat=np.column_stack([x_coords, y_coords]),
                depths=depths,
                widths=widths,
                original_id=i
            )
            centerlines.append(centerline)

        return centerlines

    def _read_shapefile(self, filename: str, zshift: float,
                       zshift_taper_elev_range: Optional[Tuple[float, float]],
                       effective_width_ratio: float, max_width: float,
                       min_width: float) -> List[CenterlineSequence]:
        """Read flowline data from Shapefile."""
        gdf = gpd.read_file(filename)
        centerlines = []

        for i, row in gdf.iterrows():
            if row.geometry is None:
                continue

            g = row.geometry
            parts = []
            if g.geom_type == 'LineString':
                parts.append(list(g.coords))
            elif g.geom_type == 'MultiLineString':
                for part in g.geoms:
                    parts.append(list(part.coords))
            else:
                continue

            for coords in parts:
                if len(coords) < 2:
                    continue
                x_coords = [float(c[0]) for c in coords]
                y_coords = [float(c[1]) for c in coords]

            depths = self._extract_depth_data(row, len(x_coords), zshift,
                                            zshift_taper_elev_range)
            widths = self._extract_width_data(row, len(x_coords),
                                            effective_width_ratio, max_width, min_width)

            centerline = CenterlineSequence(
                    lonlat=np.column_stack([x_coords, y_coords]),
                depths=depths,
                widths=widths,
                    original_id=int(i) if np.isscalar(i) else 0
            )
            centerlines.append(centerline)

        return centerlines

    def _extract_depth_data(self, properties, n_points: int, zshift: float,
                           zshift_taper_elev_range: Optional[Tuple[float, float]]) -> np.ndarray:
        """Extract depth data from properties."""
        depths = np.full(n_points, np.nan, dtype=float)

        # Try single depth value first
        if 'depth' in properties:
            depth_val = properties['depth']
            if isinstance(depth_val, str):
                depth_val = float(depth_val) if depth_val else np.nan
            if not np.isnan(depth_val):
                depths.fill(depth_val)

        # Try point-wise depth data
        if 'pt_depth' in properties:
            pt_depth_str = properties['pt_depth']
            if pt_depth_str:
                depth_values = [float(x.strip()) for x in pt_depth_str.split(',')]
                if len(depth_values) == n_points:
                    depths = np.array(depth_values)

        # Apply z-shift and tapering
        if not np.all(np.isnan(depths)):
            valid_depths = depths[~np.isnan(depths)]
            if len(valid_depths) > 0:
                # Apply tapering if specified
                if zshift_taper_elev_range is not None:
                    min_elev, max_elev = zshift_taper_elev_range
                    for i, depth in enumerate(depths):
                        if not np.isnan(depth):
                            # Calculate taper rate based on elevation
                            zshift_rate = 1 - (max_elev + depth) / (max_elev - min_elev)
                            zshift_rate = np.clip(zshift_rate, 0.0, 1.0)
                            depths[i] = depth + zshift_rate * zshift
                else:
                    depths = depths + zshift

        return depths

    def _extract_width_data(self, properties, n_points: int,
                           effective_width_ratio: float, max_width: float,
                           min_width: float) -> np.ndarray:
        """Extract width data from properties."""
        widths = np.full(n_points, np.nan, dtype=float)

        # Try single width value first
        if 'width' in properties:
            width_val = properties['width']
            if isinstance(width_val, str):
                width_val = float(width_val) if width_val else np.nan
            if not np.isnan(width_val):
                widths.fill(width_val)

        # Try BtmWdth (bottom width)
        if 'BtmWdth' in properties:
            width_val = properties['BtmWdth']
            if isinstance(width_val, str):
                width_val = float(width_val) if width_val else np.nan
            if not np.isnan(width_val):
                widths.fill(width_val)

        # Try point-wise width data
        if 'pt_width' in properties:
            pt_width_str = properties['pt_width']
            if pt_width_str:
                width_values = [float(x.strip()) for x in pt_width_str.split(',')]
                if len(width_values) == n_points:
                    widths = np.array(width_values)

        # Apply scaling and limits
        if not np.all(np.isnan(widths)):
            widths = widths * effective_width_ratio
            widths = np.clip(widths, min_width, max_width)

        return widths


class CenterlineProcessor:
    """Processes centerline widths/depths (optional smoothing)."""

    def process_centerlines(self, centerlines: List[CenterlineSequence], smoothing_span: int = 0) -> List[CenterlineSequence]:
        processed: List[CenterlineSequence] = []
        for cl in centerlines:
            if smoothing_span > 1:
                cl = self._smooth_values(cl, smoothing_span)
            processed.append(cl)
        return processed

    def _smooth_values(self, centerline: CenterlineSequence, span: int) -> CenterlineSequence:
        if span <= 1:
            return centerline
        if not np.all(np.isnan(centerline.depths)):
            valid = ~np.isnan(centerline.depths)
            if np.sum(valid) > span:
                from scipy.signal import savgol_filter
                sm = np.full_like(centerline.depths, np.nan)
                sm[valid] = savgol_filter(centerline.depths[valid], span, 2)
                centerline.depths = sm
        if not np.all(np.isnan(centerline.widths)):
            valid = ~np.isnan(centerline.widths)
            if np.sum(valid) > span:
                from scipy.signal import savgol_filter
                sm = np.full_like(centerline.widths, np.nan)
                sm[valid] = savgol_filter(centerline.widths[valid], span, 2)
                centerline.widths = sm
        return centerline

    def merge_endpoints(self, centerlines: List[CenterlineSequence], radius_m: float,
                        lon0: float, lat0: float) -> List[CenterlineSequence]:
        """Merge centerline segments whose endpoints are within radius_m.

        This approximates the MATLAB behavior controlled by radius_to_merge_shppoints.
        """
        if not centerlines or radius_m <= 0:
            return centerlines

        def deg2cart(lon: np.ndarray, lat: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
            R = 6378206.4
            lonr = np.deg2rad(lon)
            latr = np.deg2rad(lat)
            lon0r = math.radians(lon0)
            lat0r = math.radians(lat0)
            x = R * (lonr - lon0r) * math.cos(lat0r)
            y = R * (latr - 0.0)
            return x, y

        def seq_endpoints_xy(seq: CenterlineSequence) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
            lon = seq.lonlat[:, 0]
            lat = seq.lonlat[:, 1]
            x, y = deg2cart(lon, lat)
            return x[0], y[0], x[-1], y[-1]

        def reverse_seq(seq: CenterlineSequence) -> CenterlineSequence:
            return CenterlineSequence(
                lonlat=seq.lonlat[::-1].copy(),
                depths=seq.depths[::-1].copy() if seq.depths is not None else seq.depths,
                widths=seq.widths[::-1].copy() if seq.widths is not None else seq.widths,
                original_id=seq.original_id,
            )

        def concat(a: CenterlineSequence, b: CenterlineSequence, drop_duplicate: bool = True) -> CenterlineSequence:
            if drop_duplicate and a.lonlat.shape[0] > 0 and b.lonlat.shape[0] > 0 and np.allclose(a.lonlat[-1], b.lonlat[0]):
                lonlat = np.vstack([a.lonlat, b.lonlat[1:]])
                depths = np.concatenate([a.depths, b.depths[1:]]) if a.depths is not None and b.depths is not None else a.depths
                widths = np.concatenate([a.widths, b.widths[1:]]) if a.widths is not None and b.widths is not None else a.widths
            else:
                lonlat = np.vstack([a.lonlat, b.lonlat])
                depths = np.concatenate([a.depths, b.depths]) if a.depths is not None and b.depths is not None else a.depths
                widths = np.concatenate([a.widths, b.widths]) if a.widths is not None and b.widths is not None else a.widths
            return CenterlineSequence(lonlat=lonlat, depths=depths, widths=widths, original_id=a.original_id)

        def dedup_close_points(seq: CenterlineSequence, tol_m: float) -> CenterlineSequence:
            if seq.lonlat.shape[0] <= 2:
                return seq
            lon = seq.lonlat[:, 0]
            lat = seq.lonlat[:, 1]
            x, y = deg2cart(lon, lat)
            keep = [0]
            for i in range(1, len(x)):
                if math.hypot(x[i] - x[keep[-1]], y[i] - y[keep[-1]]) > tol_m * 0.5:
                    keep.append(i)
            lonlat = seq.lonlat[keep]
            depths = seq.depths[keep] if seq.depths is not None and len(seq.depths) == len(seq.lonlat) else seq.depths
            widths = seq.widths[keep] if seq.widths is not None and len(seq.widths) == len(seq.lonlat) else seq.widths
            return CenterlineSequence(lonlat=lonlat, depths=depths, widths=widths, original_id=seq.original_id)

        # Initial de-dup of each sequence
        seqs = [dedup_close_points(s, radius_m) for s in centerlines]

        changed = True
        while changed:
            changed = False
            # Build endpoint arrays
            ex = []
            ey = []
            refs = []  # (idx, end_flag 0=start,1=end)
            for idx, s in enumerate(seqs):
                if s.lonlat.shape[0] < 2:
                    continue
                x0, y0, x1, y1 = seq_endpoints_xy(s)
                ex.extend([x0, x1])
                ey.extend([y0, y1])
                refs.extend([(idx, 0), (idx, 1)])
            if not ex:
                break
            pts = np.column_stack([ex, ey])
            tree = cKDTree(pts)
            parent = list(range(len(ex)))

            def find_root(a: int) -> int:
                while parent[a] != a:
                    parent[a] = parent[parent[a]]
                    a = parent[a]
                return a

            def union_root(a: int, b: int) -> None:
                ra, rb = find_root(a), find_root(b)
                if ra != rb:
                    parent[rb] = ra

            # Build connectivity clusters of endpoints within radius
            for idx_pt in range(len(ex)):
                nbrs = tree.query_ball_point(pts[idx_pt], radius_m)
                for nb in nbrs:
                    if nb != idx_pt:
                        union_root(idx_pt, nb)

            cluster_seq_members: Dict[int, set] = {}
            for idx_pt in range(len(ex)):
                root = find_root(idx_pt)
                seq_idx_pt, _ = refs[idx_pt]
                if seq_idx_pt == -1 or seqs[seq_idx_pt].lonlat.shape[0] == 0:
                    continue
                cluster_seq_members.setdefault(root, set()).add(seq_idx_pt)

            cluster_is_multibranch = {root: (len(seq_ids) >= 3)
                                      for root, seq_ids in cluster_seq_members.items()}

            used = set()
            for i in range(len(ex)):
                if i in used:
                    continue
                idx_i, end_i = refs[i]
                if idx_i == -1 or seqs[idx_i].lonlat.shape[0] == 0:
                    continue
                # Find neighbors within radius (k=5 fallback)
                neigh = tree.query_ball_point(pts[i], radius_m)
                # Remove self
                neigh = [j for j in neigh if j != i]
                if not neigh:
                    continue
                # Choose the closest
                dists = [math.hypot(pts[j, 0] - pts[i, 0], pts[j, 1] - pts[i, 1]) for j in neigh]
                j = neigh[int(np.argmin(dists))]
                if j in used:
                    continue
                idx_j, end_j = refs[j]
                if idx_j == -1 or idx_j == idx_i or seqs[idx_j].lonlat.shape[0] == 0:
                    continue
                root_i = find_root(i)
                root_j = find_root(j)
                if cluster_is_multibranch.get(root_i, False) or cluster_is_multibranch.get(root_j, False):
                    members = sorted(cluster_seq_members.get(root_i, set()) | cluster_seq_members.get(root_j, set()))
                    print(f"INFO: Skipping merge of sequences in multi-branch cluster: seq_i={idx_i}, seq_j={idx_j}, members={members}")
                    continue
                A = seqs[idx_i]
                B = seqs[idx_j]

                # Decide how to connect based on which ends are close
                if end_i == 1 and end_j == 0:
                    merged = concat(A, B)
                elif end_i == 0 and end_j == 1:
                    merged = concat(B, A)
                elif end_i == 1 and end_j == 1:
                    merged = concat(A, reverse_seq(B))
                elif end_i == 0 and end_j == 0:
                    merged = concat(reverse_seq(A), B)
                else:
                    continue
                # Replace A with merged, delete B
                merged = dedup_close_points(merged, radius_m)
                seqs[idx_i] = merged
                # Mark B as removed
                seqs[idx_j] = CenterlineSequence(lonlat=np.zeros((0, 2)), depths=np.array([]), widths=np.array([]), original_id=seqs[idx_j].original_id)
                used.add(i)
                used.add(j)
                changed = True
            # Compact list
            if changed:
                seqs = [s for s in seqs if s.lonlat.shape[0] >= 2]
        return seqs

    def split_at_near_endpoints(self, centerlines: List[CenterlineSequence], radius_m: float,
                                lon0: float, lat0: float) -> List[CenterlineSequence]:
        """Split sequences where an endpoint is within radius of an interior node of another sequence.

        Mirrors the MATLAB clean_sequences behavior where endpoints embedded in other sequences trigger splits.
        """
        if not centerlines or radius_m <= 0:
            return centerlines

        R = 6378206.4
        lon0r = math.radians(lon0)
        lat0r = math.radians(lat0)

        def deg2cart_array(lon: np.ndarray, lat: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
            lonr = np.deg2rad(lon)
            latr = np.deg2rad(lat)
            x = R * (lonr - lon0r) * math.cos(lat0r)
            y = R * (latr - 0.0)
            return x, y

        # Build global list of points and endpoints
        point_records = []  # (seq_idx, pt_idx, x, y)
        endpoints = []      # (seq_idx, pt_idx, x, y)
        for si, s in enumerate(centerlines):
            if s.lonlat.shape[0] == 0:
                continue
            lon = s.lonlat[:, 0]
            lat = s.lonlat[:, 1]
            x, y = deg2cart_array(lon, lat)
            for pi in range(len(x)):
                point_records.append((si, pi, x[pi], y[pi]))
            endpoints.append((si, 0, x[0], y[0]))
            endpoints.append((si, len(x) - 1, x[-1], y[-1]))

        if not point_records:
            return centerlines

        pts_xy = np.array([[pr[2], pr[3]] for pr in point_records])
        tree = cKDTree(pts_xy)

        # Map seq -> split indices
        split_map: Dict[int, List[int]] = {}
        for e in endpoints:
            esi, epi, ex, ey = e
            # Query neighbors within radius
            idxs = tree.query_ball_point([ex, ey], radius_m)
            for gi in idxs:
                si, pi, px, py = point_records[gi]
                if si == esi:
                    continue
                # interior point only
                npts = centerlines[si].lonlat.shape[0]
                if pi <= 0 or pi >= npts - 1:
                    continue
                split_map.setdefault(si, []).append(pi)

        if not split_map:
            return centerlines

        # Perform splits per sequence
        new_sequences: List[CenterlineSequence] = []
        for si, s in enumerate(centerlines):
            if si not in split_map:
                new_sequences.append(s)
                continue
            # Unique, sorted split indices
            cuts = sorted(set(split_map[si]))
            # Ensure cuts are strictly increasing and within bounds
            cuts = [c for c in cuts if 0 < c < s.lonlat.shape[0] - 1]
            if not cuts:
                new_sequences.append(s)
                continue
            last = 0
            for c in cuts:
                seg_lonlat = s.lonlat[last:c + 1]
                seg_depths = s.depths[last:c + 1] if s.depths is not None and len(s.depths) == len(s.lonlat) else s.depths
                seg_widths = s.widths[last:c + 1] if s.widths is not None and len(s.widths) == len(s.lonlat) else s.widths
                if seg_lonlat.shape[0] >= 2:
                    new_sequences.append(CenterlineSequence(lonlat=seg_lonlat, depths=seg_depths, widths=seg_widths, original_id=s.original_id))
                last = c
            # tail
            seg_lonlat = s.lonlat[last:]
            seg_depths = s.depths[last:] if s.depths is not None and len(s.depths) == len(s.lonlat) else s.depths
            seg_widths = s.widths[last:] if s.widths is not None and len(s.widths) == len(s.lonlat) else s.widths
            if seg_lonlat.shape[0] >= 2:
                new_sequences.append(CenterlineSequence(lonlat=seg_lonlat, depths=seg_depths, widths=seg_widths, original_id=s.original_id))

        return new_sequences
class EndpointCrosspointBuilder:
    """Compute crosspoint coordinates at multi-branch endpoints following msh_pave_channels.m logic."""

    EARTH_RADIUS_M = 6378206.4

    def __init__(self):
        self.gid_counter = 1  # Global ID counter like MATLAB gid_cnt

    @staticmethod
    def _deg2cart(lon: np.ndarray, lat: np.ndarray, lon0: float, lat0: float) -> Tuple[np.ndarray, np.ndarray]:
        lonr = np.deg2rad(lon)
        latr = np.deg2rad(lat)
        lon0r = math.radians(lon0)
        lat0r = math.radians(lat0)
        x = EndpointCrosspointBuilder.EARTH_RADIUS_M * (lonr - lon0r) * math.cos(lat0r)
        y = EndpointCrosspointBuilder.EARTH_RADIUS_M * (latr - 0.0)
        return x, y

    @staticmethod
    def _unit(vx: float, vy: float) -> Tuple[float, float]:
        n = math.hypot(vx, vy)
        if n == 0:
            return 1.0, 0.0
        return vx / n, vy / n

    @staticmethod
    def _angle(vec: Tuple[float, float]) -> float:
        return math.atan2(vec[1], vec[0])

    def compute_overrides(self, resampled: List[Dict], lon0: float, lat0: float, radius_merge: float) -> Tuple[Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]], Dict[int, List[Dict]], Dict[Tuple[int, int], Tuple[int, int]], Dict[Tuple[int, int], int]]:
        """
        Build overrides for endpoint bank coordinates using crosspoints.
        Follows MATLAB shp.endpoints.ids logic more closely.

        Returns:
            Tuple of (overrides_dict, junction_crosspoints_dict, crosspoint_gids_dict)
            - overrides_dict: (seq_idx, end_flag 0=start,1=end) -> (cp1_lonlat, cp2_lonlat)
            - junction_crosspoints_dict: group_id -> list of crosspoint coordinates for junction elements
            - crosspoint_gids_dict: (seq_idx, end_flag) -> (left_gid, right_gid) for tracking node IDs
        """
        # Build endpoint records like MATLAB shp.endpoints.ids
        endpoint_records = []
        for si, item in enumerate(resampled):
            lonlat = item['lonlat']
            widths = item['widths']
            depths = item['depths']
            if lonlat.shape[0] < 2:
                continue
            
            # Start endpoint: endpoint=lonlat[0], neighbor=lonlat[1]
            endpoint_records.append({
                'endpoint': tuple(lonlat[0]),
                'neighbor': tuple(lonlat[1]),
                'seq_idx': si,
                'end_flag': 0,
                'width': float(widths[0] if len(widths) else 0.0),
                'depth': float(depths[0] if len(depths) else 0.0)
            })
            
            # End endpoint: endpoint=lonlat[-1], neighbor=lonlat[-2]
            endpoint_records.append({
                'endpoint': tuple(lonlat[-1]),
                'neighbor': tuple(lonlat[-2]),
                'seq_idx': si,
                'end_flag': 1,
                'width': float(widths[-1] if len(widths) else 0.0),
                'depth': float(depths[-1] if len(depths) else 0.0)
            })

        if not endpoint_records:
            return {}, {}, {}, {}

        # Group endpoints by spatial proximity (like MATLAB unique endpoint IDs)
        endpoint_coords = []
        
        for i, rec in enumerate(endpoint_records):
            x, y = EndpointCrosspointBuilder._deg2cart(
                np.array([rec['endpoint'][0]]), np.array([rec['endpoint'][1]]), lon0, lat0
            )
            endpoint_coords.append([x[0], y[0]])
        
        endpoint_coords = np.array(endpoint_coords)
        
        tree = cKDTree(endpoint_coords)
        n = len(endpoint_records)
        parent = list(range(n))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        for i in range(n):
            nbrs = tree.query_ball_point(endpoint_coords[i], radius_merge)
            for j in nbrs:
                union(i, j)

        groups: Dict[int, List[int]] = {}
        for i in range(n):
            r = find(i)
            groups.setdefault(r, []).append(i)

        overrides: Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]] = {}
        junction_crosspoints: Dict[int, List[Dict]] = {}
        crosspoint_gids: Dict[Tuple[int, int], Tuple[int, int]] = {}  # (seq_idx, end_flag) -> (gid1, gid2)
        junction_endpoints: Dict[Tuple[int, int], int] = {}  # (seq_idx, end_flag) -> group_id
        
        # Process each connected endpoint group
        for group_id, idxs in groups.items():
            if len(idxs) < 2:
                continue
                
            # Get the mean endpoint location in Cartesian
            x0 = float(np.mean([endpoint_coords[i][0] for i in idxs]))
            y0 = float(np.mean([endpoint_coords[i][1] for i in idxs]))
            
            # Build direction vectors for each branch
            entries = []
            for i in idxs:
                rec = endpoint_records[i]
                # Direction from endpoint to neighbor
                ep_x, ep_y = EndpointCrosspointBuilder._deg2cart(
                    np.array([rec['endpoint'][0]]), np.array([rec['endpoint'][1]]), lon0, lat0
                )
                nb_x, nb_y = EndpointCrosspointBuilder._deg2cart(
                    np.array([rec['neighbor'][0]]), np.array([rec['neighbor'][1]]), lon0, lat0
                )
                ux, uy = EndpointCrosspointBuilder._unit(nb_x[0] - ep_x[0], nb_y[0] - ep_y[0])
                ang = EndpointCrosspointBuilder._angle((ux, uy))
                entries.append((i, ang, ux, uy))
            entries.sort(key=lambda t: t[1])

            # Compute crosspoints between adjacent branch pairs
            crossps = []
            for k in range(len(entries)):
                i_idx = entries[k][0]
                j_idx = entries[(k + 1) % len(entries)][0]
                
                rec1 = endpoint_records[i_idx]
                rec2 = endpoint_records[j_idx]
                
                ux1, uy1 = entries[k][2], entries[k][3]
                ux2, uy2 = entries[(k + 1) % len(entries)][2], entries[(k + 1) % len(entries)][3]
                
                w1 = max(0.0, rec1['width'])
                w2 = max(0.0, rec2['width'])
                d1 = rec1['depth']
                d2 = rec2['depth']

                dotv = ux1 * ux2 + uy1 * uy2
                if dotv < -0.90:
                    # Opposite directions case
                    nb1_x, nb1_y = EndpointCrosspointBuilder._deg2cart(
                        np.array([rec1['neighbor'][0]]), np.array([rec1['neighbor'][1]]), lon0, lat0
                    )
                    nb2_x, nb2_y = EndpointCrosspointBuilder._deg2cart(
                        np.array([rec2['neighbor'][0]]), np.array([rec2['neighbor'][1]]), lon0, lat0
                    )
                    x1n, y1n = nb1_x[0], nb1_y[0]
                    x2n, y2n = nb2_x[0], nb2_y[0]
                    tx = x2n - x1n
                    ty = y2n - y1n
                    L = math.hypot(tx, ty)
                    if L == 0:
                        nx, ny = -uy1, ux1
                    else:
                        nx, ny = ty / L, -tx / L
                    ln = 0.25 * (w1 + w2)
                    p1x = x0 + nx * ln
                    p1y = y0 + ny * ln
                else:
                    # General case: intersect offset lines
                    n1x, n1y = -uy1, ux1
                    n2x, n2y = -uy2, ux2
                    p11x = x0 + 0.5 * w1 * n1x
                    p11y = y0 + 0.5 * w1 * n1y
                    p21x = x0 - 0.5 * w2 * n2x
                    p21y = y0 - 0.5 * w2 * n2y
                    dett = -ux1 * uy2 + ux2 * uy1
                    if abs(dett) < 1e-12:
                        # Fallback: average
                        p1x = 0.5 * (p11x + p21x)
                        p1y = 0.5 * (p11y + p21y)
                    else:
                        px = -p11x + p21x
                        py = -p11y + p21y
                        t1 = (-uy2 * px + ux2 * py) / dett
                        p1x = p11x + t1 * ux1
                        p1y = p11y + t1 * uy1
                        
                dp = 0.5 * (d1 + d2)
                gid = self.gid_counter
                self.gid_counter += 1
                crossps.append({'x': p1x, 'y': p1y, 'depth': dp, 'gid': gid})
                
            # Store crosspoints for junction element generation (if 3+ branches)
            if len(idxs) >= 3:
                crosspoint_coords = []
                for cp in crossps:
                    lon_cp, lat_cp = EndpointCrosspointBuilder._cart2deg_single(cp['x'], cp['y'], lon0, lat0)
                    crosspoint_coords.append({
                        'lon': lon_cp, 
                        'lat': lat_cp, 
                        'depth': cp['depth'],
                        'gid': cp['gid']
                    })
                junction_crosspoints[group_id] = crosspoint_coords
                # Mark all member endpoints in this group as part of a junction
                for i_ep in idxs:
                    rec_ep = endpoint_records[i_ep]
                    junction_endpoints[(rec_ep['seq_idx'], rec_ep['end_flag'])] = group_id

            # Assign crosspoint pairs to each branch following MATLAB logic
            for idx_order, (i_idx, _, _, _) in enumerate(entries):
                rec = endpoint_records[i_idx]
                seq_idx = rec['seq_idx']
                end_flag = rec['end_flag']
                
                # MATLAB logic: crossp1 is current, crossp2 is previous (with wraparound)
                cp1_idx = idx_order
                cp2_idx = (idx_order - 1) % len(crossps) if idx_order > 0 else len(crossps) - 1
                
                cp1 = crossps[cp1_idx]
                cp2 = crossps[cp2_idx]
                
                lon1, lat1 = EndpointCrosspointBuilder._cart2deg_single(cp1['x'], cp1['y'], lon0, lat0)
                lon2, lat2 = EndpointCrosspointBuilder._cart2deg_single(cp2['x'], cp2['y'], lon0, lat0)
                
                # MATLAB bank assignment logic:
                # Start (end_flag=0): [crossps1, crossps2] = [left, right] -> (cp1, cp2)
                # End (end_flag=1): [crosspe2, crosspe1] = [right, left] -> (cp2, cp1)
                if end_flag == 0:  # start endpoint
                    overrides[(seq_idx, end_flag)] = (np.array([lon1, lat1]), np.array([lon2, lat2]))  # (left=cp1, right=cp2)
                    crosspoint_gids[(seq_idx, end_flag)] = (cp1['gid'], cp2['gid'])  # (left_gid, right_gid)
                else:  # end endpoint
                    overrides[(seq_idx, end_flag)] = (np.array([lon2, lat2]), np.array([lon1, lat1]))  # (left=cp2, right=cp1)
                    crosspoint_gids[(seq_idx, end_flag)] = (cp2['gid'], cp1['gid'])  # (left_gid, right_gid)

        return overrides, junction_crosspoints, crosspoint_gids, junction_endpoints

    @staticmethod
    def _cart2deg_single(x: float, y: float, lon0: float, lat0: float) -> Tuple[float, float]:
        lon0r = math.radians(lon0)
        lat0r = math.radians(lat0)
        lon = math.degrees(lon0r + x / (EndpointCrosspointBuilder.EARTH_RADIUS_M * math.cos(lat0r)))
        lat = math.degrees(y / EndpointCrosspointBuilder.EARTH_RADIUS_M)
        return lon, lat


class ChannelStripGenerator:
    """Generate bank-to-bank strips for each centerline with uniform spacing."""

    EARTH_RADIUS_M = 6378206.4

    @staticmethod
    def _deg2cart(lon: np.ndarray, lat: np.ndarray, lon0: float, lat0: float) -> Tuple[np.ndarray, np.ndarray]:
        lonr = np.deg2rad(lon)
        latr = np.deg2rad(lat)
        lon0r = math.radians(lon0)
        lat0r = math.radians(lat0)
        x = ChannelStripGenerator.EARTH_RADIUS_M * (lonr - lon0r) * math.cos(lat0r)
        y = ChannelStripGenerator.EARTH_RADIUS_M * (latr - 0.0)
        return x, y

    @staticmethod
    def _cart2deg(x: np.ndarray, y: np.ndarray, lon0: float, lat0: float) -> Tuple[np.ndarray, np.ndarray]:
        lon0r = math.radians(lon0)
        lat0r = math.radians(lat0)
        lon = np.rad2deg(lon0r + x / (ChannelStripGenerator.EARTH_RADIUS_M * math.cos(lat0r)))
        lat = np.rad2deg(y / ChannelStripGenerator.EARTH_RADIUS_M)
        return lon, lat

    @staticmethod
    def _interp_nan_safe(s_old: np.ndarray, values: np.ndarray, s_new: np.ndarray, fill: float) -> np.ndarray:
        out = np.full_like(s_new, fill, dtype=float)
        valid = ~np.isnan(values)
        if np.any(valid):
            out = np.interp(s_new, s_old[valid], values[valid], left=values[valid][0], right=values[valid][-1])
        return out

    def resample_centerline(self, lonlat: np.ndarray, widths: np.ndarray, depths: np.ndarray,
                             spacing_m: float, lon0: float, lat0: float,
                             keep_bend_node_degree: Optional[float] = None,
                             endpoint_angle_start_rad: Optional[float] = None,
                             endpoint_angle_end_rad: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        lon = lonlat[:, 0]
        lat = lonlat[:, 1]
        x, y = self._deg2cart(lon, lat, lon0, lat0)
        dx = np.diff(x)
        dy = np.diff(y)
        seg = np.hypot(dx, dy)
        if len(seg) == 0:
            return lonlat.copy(), widths.copy(), depths.copy()
        s = np.concatenate([[0.0], np.cumsum(seg)])
        L = s[-1]
        if L == 0.0:
            return lonlat.copy(), widths.copy(), depths.copy()
        nsteps = max(2, int(np.floor(L / max(spacing_m, 1e-6))) + 1)
        # Build resampling arc-lengths
        if keep_bend_node_degree is not None and keep_bend_node_degree > 0.0 and lonlat.shape[0] >= 3:
            # Identify anchor indices: endpoints + sharp bends (angle >= threshold)
            thr_rad = math.radians(float(keep_bend_node_degree))
            anchor_idx = [0]
            anchor_angles = [0.0]
            for i in range(1, lonlat.shape[0] - 1):
                v1x = x[i] - x[i - 1]
                v1y = y[i] - y[i - 1]
                v2x = x[i + 1] - x[i]
                v2y = y[i + 1] - y[i]
                n1 = math.hypot(v1x, v1y)
                n2 = math.hypot(v2x, v2y)
                if n1 == 0.0 or n2 == 0.0:
                    continue
                dot = v1x * v2x + v1y * v2y
                cs = max(-1.0, min(1.0, dot / (n1 * n2)))
                ang = math.acos(cs)
                if ang >= thr_rad:
                    anchor_idx.append(i)
                    anchor_angles.append(ang)
            anchor_idx.append(len(lonlat) - 1)
            anchor_angles.append(0.0)
            # Sort anchor_idx and anchor_angles together according to anchor_idx
            anchor_pairs = sorted(zip(anchor_idx, anchor_angles), key=lambda pair: pair[0])
            anchor_idx, anchor_angles = zip(*anchor_pairs)
            anchor_idx = list(anchor_idx)
            anchor_angles = list(anchor_angles)

            # Equispace between consecutive anchors
            s_new_list = []
            for k in range(len(anchor_idx) - 1):
                a = anchor_idx[k]
                b = anchor_idx[k + 1]
                if s[b] <= s[a]:
                    continue
                # Lk = s[b] - s[a]
                # nsteps_k = max(2, int(np.floor(Lk / max(spacing_m, 1e-6))) + 1)
                # seg_s = np.linspace(s[a], s[b], nsteps_k)
                # Use endpoint override angles for endpoints, curvature angles for interior anchors
                ang_a = endpoint_angle_start_rad if (a == 0 and endpoint_angle_start_rad is not None) else anchor_angles[k]
                ang_b = endpoint_angle_end_rad if (b == len(lonlat) - 1 and endpoint_angle_end_rad is not None) else anchor_angles[k + 1]
                wk_a = widths[a] if widths[a] is not None and not np.isnan(widths[a]) else 0.0
                wk_b = widths[b] if widths[b] is not None and not np.isnan(widths[b]) else 0.0
                shrink_a = 0.5 * math.tan(0.5 * ang_a) * wk_a
                shrink_b = 0.5 * math.tan(0.5 * ang_b) * wk_b
                Lk = s[b] - s[a] - shrink_a - shrink_b
                nsteps_k = max(2, int(np.floor(Lk / max(spacing_m, 1e-6))) + 1)
                seg_s = np.linspace(s[a] + shrink_a, s[b] - shrink_b, nsteps_k)
                seg_s[0] = s[a]
                seg_s[-1] = s[b]
                if k > 0 and len(seg_s) > 0:
                    seg_s = seg_s[1:]  # avoid duplicating the shared anchor
                s_new_list.append(seg_s)
            s_new = np.concatenate(s_new_list) if s_new_list else np.array([s[0], s[-1]])
        else:
            s_new = np.linspace(0.0, L, nsteps)
        x_new = np.interp(s_new, s, x)
        y_new = np.interp(s_new, s, y)
        lon_new, lat_new = self._cart2deg(x_new, y_new, lon0, lat0)
        widths_new = self._interp_nan_safe(s, widths, s_new, fill=
                                           np.nan if len(widths) == 0 else np.nan)
        depths_new = self._interp_nan_safe(s, depths, s_new, fill=
                                           np.nan if len(depths) == 0 else np.nan)

        lonlat_new = np.column_stack([lon_new, lat_new])
        return lonlat_new, widths_new, depths_new

    def make_strip(self, lonlat: np.ndarray, widths_m: np.ndarray, lon0: float, lat0: float,
                   default_width: float, angles: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
        n = lonlat.shape[0]
        lon = lonlat[:, 0]
        lat = lonlat[:, 1]
        x, y = self._deg2cart(lon, lat, lon0, lat0)
        left_x = np.zeros(n)
        left_y = np.zeros(n)
        right_x = np.zeros(n)
        right_y = np.zeros(n)
        for i in range(n):
            if i == 0:
                tx = x[1] - x[0]
                ty = y[1] - y[0]
            elif i == n - 1:
                tx = x[i] - x[i - 1]
                ty = y[i] - y[i - 1]
            else:
                tx = (x[i + 1] - x[i - 1]) * 0.5
                ty = (y[i + 1] - y[i - 1]) * 0.5
            norm = math.hypot(tx, ty)
            if norm == 0:
                px, py = 0.0, 1.0
            else:
                px = -ty / norm
                py = tx / norm
            # Use provided angle array if given; else fall back to curvature-based estimate
            if angles is not None and len(angles) == n:
                ang_i = float(angles[i])
            else:
                if i == 0 or i == n - 1:
                    ang_i = 0.0
                else:
                    v1x = x[i] - x[i - 1]
                    v1y = y[i] - y[i - 1]
                    v2x = x[i + 1] - x[i]
                    v2y = y[i + 1] - y[i]
                    n1 = math.hypot(v1x, v1y)
                    n2 = math.hypot(v2x, v2y)
                    if n1 == 0.0 or n2 == 0.0:
                        ang_i = 0.0
                    else:
                        dot = v1x * v2x + v1y * v2y
                        cs = max(-1.0, min(1.0, dot / (n1 * n2)))
                        ang_i = math.acos(cs)
            # Adjust half width by curvature: half_w_adjusted = half_w / cos(0.5 * ang_i)
            half_w = 0.5 * (widths_m[i] if not np.isnan(widths_m[i]) else default_width)
            denom = math.cos(0.5 * ang_i)
            if abs(denom) < 1e-6:
                denom = 1e-6
            half_w_adjusted = half_w / denom
            left_x[i] = x[i] + half_w_adjusted * px
            left_y[i] = y[i] + half_w_adjusted * py
            right_x[i] = x[i] - half_w_adjusted * px
            right_y[i] = y[i] - half_w_adjusted * py
        left_lon, left_lat = self._cart2deg(left_x, left_y, lon0, lat0)
        right_lon, right_lat = self._cart2deg(right_x, right_y, lon0, lat0)
        return np.column_stack([left_lon, left_lat]), np.column_stack([right_lon, right_lat])

    @staticmethod
    def add_channel_midnodes(left_lonlat: np.ndarray, right_lonlat: np.ndarray,
                              center_depths: np.ndarray, lon0: float, lat0: float,
                              spacing_factor: Optional[float] = None,
                              spacing_meters: Optional[float] = None,
                              start_is_junction: bool = False,
                              end_is_junction: bool = False) -> Tuple[list, list, list]:
        """Add mid nodes across wide channel cross-sections.

        Args:
            spacing_factor: Multiplier applied to longitudinal spacing (mutually exclusive with spacing_meters)
            spacing_meters: Absolute spacing in meters for mid-nodes (mutually exclusive with spacing_factor)

        Returns lists (mid_lonlat_per_i, mid_depths_per_i, mid_gids_per_i) aligned to cross-section index.
        GIDs are placeholders (None) and will be assigned later in the builder.
        """
        # Validate mutual exclusivity and set defaults
        if spacing_factor is not None and spacing_meters is not None:
            raise ValueError("Cannot specify both spacing_factor and spacing_meters")
        if spacing_factor is None and spacing_meters is None:
            spacing_factor = 1.0  # Default behavior
        n = left_lonlat.shape[0]
        if n != right_lonlat.shape[0] or n < 3:
            return [None] * n, [None] * n, [None] * n
        # Build approximate longitudinal spacing from arc-length of the center (midpoint of banks)
        mid = 0.5 * (left_lonlat + right_lonlat)
        R = ChannelStripGenerator.EARTH_RADIUS_M
        lon0r = math.radians(lon0)
        lat0r = math.radians(lat0)
        def ll2xy(arr):
            lonr = np.deg2rad(arr[:, 0])
            latr = np.deg2rad(arr[:, 1])
            x = R * (lonr - lon0r) * math.cos(lat0r)
            y = R * (latr - 0.0)
            return x, y
        mx, my = ll2xy(mid)
        seg = np.hypot(np.diff(mx), np.diff(my))
        s = np.concatenate([[0.0], np.cumsum(seg)])
        
        # Compute per cross-section spacing based on mode
        if spacing_meters is not None:
            # Use absolute spacing in meters
            s_step = np.full_like(s, max(1e-6, float(spacing_meters)))
        else:
            # Use spacing factor mode (scale by longitudinal step)
            s_diff = np.diff(s)
            s_step = np.zeros_like(s)
            if len(s_diff) > 0:
                s_step[1:] = s_diff
                s_step[0] = s_diff[0]
            s_step = np.maximum(1e-6, s_step) * max(0.0, float(spacing_factor))
        out_lonlat = [None] * n
        out_depths = [None] * n
        out_gids = [None] * n
        # Helper to generate mid nodes for a given cross-section index
        def build_midnodes_at(i: int) -> None:
            pl = left_lonlat[i]
            pr = right_lonlat[i]
            lx, ly = ll2xy(np.array([pl]))
            rx, ry = ll2xy(np.array([pr]))
            dx = rx[0] - lx[0]
            dy = ry[0] - ly[0]
            wid = math.hypot(dx, dy)
            step = max(1e-6, s_step[i])
            nn_add = int(math.floor(wid / step))
            if nn_add >= 1:
                # generate mid points (exclude endpoints)
                pts = []
                for j in range(1, nn_add + 1):
                    t = j / (nn_add + 1)
                    xm = lx[0] + (rx[0] - lx[0]) * t
                    ym = ly[0] + (ry[0] - ly[0]) * t
                    lonm = math.degrees(lon0r + xm / (R * math.cos(lat0r)))
                    latm = math.degrees(ym / R)
                    pts.append([lonm, latm])
                out_lonlat[i] = np.array(pts)
                out_depths[i] = np.full(len(pts), float(center_depths[i]))
                out_gids[i] = [None] * len(pts)

        # Iterate cross-sections including endpoints (subject to junction flags)
        for i in range(1, n - 1):
            build_midnodes_at(i)
        if not start_is_junction:
            build_midnodes_at(0)
        if not end_is_junction:
            build_midnodes_at(n - 1)
        return out_lonlat, out_depths, out_gids

    @staticmethod
    def compute_angles_for_resampled(lonlat_res: np.ndarray, lon0: float, lat0: float,
                                     start_is_junction: bool, end_is_junction: bool,
                                     start_override: Optional[Tuple[np.ndarray, np.ndarray]],
                                     end_override: Optional[Tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
        """Compute angle at each resampled point.

        - Interior points: curvature-based interior angle (between adjacent segments).
        - Non-junction endpoints: 0.
        - Junction endpoints: angle between centerline direction at endpoint and normal to the crosspoint chord.
        """
        n = lonlat_res.shape[0]
        angles = np.zeros(n, dtype=float)
        if n <= 1:
            return angles
        # Cartesian for direction calculations
        lon = lonlat_res[:, 0]
        lat = lonlat_res[:, 1]
        lonr = np.deg2rad(lon)
        latr = np.deg2rad(lat)
        lon0r = math.radians(lon0)
        lat0r = math.radians(lat0)
        R = ChannelStripGenerator.EARTH_RADIUS_M
        x = R * (lonr - lon0r) * math.cos(lat0r)
        y = R * (latr - 0.0)
        # Interior angles
        for i in range(1, n - 1):
            v1x = x[i] - x[i - 1]
            v1y = y[i] - y[i - 1]
            v2x = x[i + 1] - x[i]
            v2y = y[i + 1] - y[i]
            n1 = math.hypot(v1x, v1y)
            n2 = math.hypot(v2x, v2y)
            if n1 == 0.0 or n2 == 0.0:
                angles[i] = 0.0
            else:
                dot = v1x * v2x + v1y * v2y
                cs = max(-1.0, min(1.0, dot / (n1 * n2)))
                angles[i] = math.acos(cs)
        # Endpoint angles
        # Helper to compute angle between direction vector and normal to crosspoint chord
        def endpoint_angle(p_idx: int, cp_pair: Optional[Tuple[np.ndarray, np.ndarray]]) -> float:
            if cp_pair is None:
                return 0.0
            # centerline direction at endpoint
            if p_idx == 0 and n >= 2:
                vx, vy = x[1] - x[0], y[1] - y[0]
            elif p_idx == n - 1 and n >= 2:
                vx, vy = x[n - 1] - x[n - 2], y[n - 1] - y[n - 2]
            else:
                return 0.0
            nv = math.hypot(vx, vy)
            if nv == 0.0:
                return 0.0
            vx /= nv
            vy /= nv
            # crosspoint chord and its normal
            cp1, cp2 = cp_pair
            # convert cp lonlat to Cartesian
            c1x = R * (math.radians(cp1[0]) - lon0r) * math.cos(lat0r)
            c1y = R * (math.radians(cp1[1]) - 0.0)
            c2x = R * (math.radians(cp2[0]) - lon0r) * math.cos(lat0r)
            c2y = R * (math.radians(cp2[1]) - 0.0)
            cx, cy = c2x - c1x, c2y - c1y
            cn = math.hypot(cx, cy)
            if cn == 0.0:
                return 0.0
            # Two possible normals; pick the one giving the smaller angle
            nx1, ny1 = -cy / cn, cx / cn
            nx2, ny2 = -nx1, -ny1
            # angle between v and n = arccos( dot(v, n) )
            a1 = math.acos(max(-1.0, min(1.0, vx * nx1 + vy * ny1)))
            a2 = math.acos(max(-1.0, min(1.0, vx * nx2 + vy * ny2)))
            return min(a1, a2)
        # Start endpoint
        if start_is_junction:
            start_pair = start_override if start_override is not None else None
            angles[0] = endpoint_angle(0, start_pair)
        else:
            angles[0] = 0.0
        # End endpoint
        if end_is_junction:
            end_pair = end_override if end_override is not None else None
            angles[n - 1] = endpoint_angle(n - 1, end_pair)
        else:
            angles[n - 1] = 0.0
        return angles

    @staticmethod
    def detect_bend_nodes(lonlat: np.ndarray, lon0: float, lat0: float, threshold_deg: float) -> List[Dict]:
        """Detect interior nodes with interior angle <= threshold (degrees).

        Returns list of dicts with keys: pt_idx, lon, lat, angle_deg.
        """
        if lonlat.shape[0] < 3 or threshold_deg is None or threshold_deg <= 0.0:
            return []
        lon = lonlat[:, 0]
        lat = lonlat[:, 1]
        x, y = ChannelStripGenerator._deg2cart(lon, lat, lon0, lat0)
        thr_rad = math.radians(float(threshold_deg))
        out = []
        for i in range(1, lonlat.shape[0] - 1):
            v1x = x[i] - x[i - 1]
            v1y = y[i] - y[i - 1]
            v2x = x[i + 1] - x[i]
            v2y = y[i + 1] - y[i]
            n1 = math.hypot(v1x, v1y)
            n2 = math.hypot(v2x, v2y)
            if n1 == 0.0 or n2 == 0.0:
                continue
            dot = v1x * v2x + v1y * v2y
            cs = max(-1.0, min(1.0, dot / (n1 * n2)))
            ang = math.acos(cs)
            if ang >= thr_rad:
                out.append({'pt_idx': i, 'lon': float(lon[i]), 'lat': float(lat[i]), 'angle_deg': math.degrees(ang)})
        return out


class JunctionMeshGenerator:
    """Generate triangular elements at junctions using crosspoints."""
    
    @staticmethod
    def generate_junction_elements(junction_crosspoints: Dict[int, List[Dict]], 
                                 lon0: float, lat0: float) -> Tuple[List[pd.DataFrame], List[List[int]]]:
        """
        Generate triangular elements at junctions using the crosspoints.

        Args:
            junction_crosspoints: Dict mapping group_id -> list of crosspoint coordinates
            lon0, lat0: Reference coordinates for coordinate conversion

        Returns:
            Tuple of (junction_nodes_list, junction_elements_list)
        """
        junction_nodes_list = []
        junction_elements_list = []
        
        for group_id, crosspoints in junction_crosspoints.items():
            if len(crosspoints) < 3:
                continue
                
            # For a 3-way junction, create a single triangle using the 3 crosspoints
            if len(crosspoints) == 3:
                # Extract coordinates and depths
                coords = []
                depths = []
                gids = []
                for cp in crosspoints:
                    coords.append([cp['lon'], cp['lat']])
                    depths.append(cp['depth'])
                    gids.append(cp.get('gid'))
                
                coords = np.array(coords)
                
                # Create nodes DataFrame
                nodes_df = pd.DataFrame({
                    'x': coords[:, 0],
                    'y': coords[:, 1], 
                    'value_1': depths,
                    'gid': gids
                })
                
                junction_nodes_list.append(nodes_df)
                
                # Create single triangular element (indices 0, 1, 2)
                junction_elements_list.append([0, 1, 2])
                
            elif len(crosspoints) > 3:
                # For more complex junctions, use Delaunay triangulation
                try:
                    from scipy.spatial import Delaunay
                    
                    # Convert to Cartesian for triangulation
                    R = 6378206.4
                    lon0r = math.radians(lon0)
                    lat0r = math.radians(lat0)
                    
                    coords_cart = []
                    depths = []
                    
                    gids = []
                    for cp in crosspoints:
                        lon, lat = cp['lon'], cp['lat']
                        lonr = math.radians(lon)
                        latr = math.radians(lat)
                        x = R * (lonr - lon0r) * math.cos(lat0r)
                        y = R * (latr - 0.0)
                        coords_cart.append([x, y])
                        depths.append(cp['depth'])
                        gids.append(cp.get('gid'))
                    
                    coords_cart = np.array(coords_cart)
                    
                    # Create Delaunay triangulation
                    tri = Delaunay(coords_cart)
                    
                    # Convert back to lon/lat for nodes
                    coords_lonlat = []
                    for cp in crosspoints:
                        coords_lonlat.append([cp['lon'], cp['lat']])
                    
                    coords_lonlat = np.array(coords_lonlat)
                    
                    # Filter out tiny elements (area threshold)
                    valid_elements = []
                    for simplex in tri.simplices:
                        # Calculate triangle area in Cartesian space
                        p1, p2, p3 = coords_cart[simplex]
                        area = 0.5 * abs((p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]))
                        if area > 1e-2:  # Same threshold as MATLAB
                            valid_elements.append(simplex.tolist())
                    
                    if valid_elements:
                        # Create nodes DataFrame
                        nodes_df = pd.DataFrame({
                            'x': coords_lonlat[:, 0],
                            'y': coords_lonlat[:, 1], 
                            'value_1': depths,
                            'gid': gids
                        })
                        
                        junction_nodes_list.append(nodes_df)
                        junction_elements_list.extend(valid_elements)
                        
                except ImportError:
                    # Fallback: create a simple fan triangulation from first point
                    coords = []
                    depths = []
                    gids = []
                    for cp in crosspoints:
                        coords.append([cp['lon'], cp['lat']])
                        depths.append(cp['depth'])
                        gids.append(cp.get('gid'))
                    
                    coords = np.array(coords)
                    
                    # Create nodes DataFrame
                    nodes_df = pd.DataFrame({
                        'x': coords[:, 0],
                        'y': coords[:, 1], 
                        'value_1': depths,
                        'gid': gids
                    })
                    
                    junction_nodes_list.append(nodes_df)
                    
                    # Create fan triangulation from point 0
                    for i in range(1, len(crosspoints) - 1):
                        junction_elements_list.append([0, i, i + 1])
        
        return junction_nodes_list, junction_elements_list

    @staticmethod
    def generate_junction_elements_from_nodes(junction_nodes_list: List[pd.DataFrame],
                                              lon0: float,
                                              lat0: float,
                                              interior_spacing: Optional[float] = None
                                              ) -> Tuple[List[pd.DataFrame], List[List[List[int]]]]:
        """Generate constrained triangulations per junction ring.

        - Uses provided ring order (including any chord mid-nodes).
        - Optionally seeds interior Steiner points on a grid with spacing=interior_spacing (meters).
        - Delaunay triangulates in local Cartesian and filters triangles by polygon containment.

        Returns elements grouped per ring to preserve mapping: List[ring_elems].
        """
        def ll2xy(arr: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
            R = 6378206.4
            lon0r = math.radians(lon0)
            lat0r = math.radians(lat0)
            lonr = np.deg2rad(arr[:, 0])
            latr = np.deg2rad(arr[:, 1])
            x = R * (lonr - lon0r) * math.cos(lat0r)
            y = R * (latr - 0.0)
            return x, y

        def point_in_polygon(px: float, py: float, poly_x: np.ndarray, poly_y: np.ndarray) -> bool:
            inside = False
            j = len(poly_x) - 1
            for i in range(len(poly_x)):
                xi, yi = poly_x[i], poly_y[i]
                xj, yj = poly_x[j], poly_y[j]
                intersect = ((yi > py) != (yj > py)) and (
                    px < (xj - xi) * (py - yi) / (yj - yi + 1e-20) + xi
                )
                if intersect:
                    inside = not inside
                j = i
            return inside

        elements_grouped: List[List[List[int]]] = []
        updated_nodes_list: List[pd.DataFrame] = []

        for ring_df in junction_nodes_list:
            if ring_df is None or ring_df.empty or len(ring_df) < 3:
                elements_grouped.append([])
                updated_nodes_list.append(ring_df)
                continue

            # Ring in lon/lat
            coords_lonlat = ring_df[['x', 'y']].to_numpy(dtype=float)
            xv, yv = ll2xy(coords_lonlat)
            poly_x, poly_y = xv.copy(), yv.copy()

            # Determine spacing in meters for interior seeding
            # Fallback to mean edge length if not provided
            edges_len = np.hypot(np.diff(poly_x, append=poly_x[0]), np.diff(poly_y, append=poly_y[0]))
            mean_edge = float(np.mean(edges_len)) if len(edges_len) else 0.0
            spacing_m = float(interior_spacing) if interior_spacing and interior_spacing > 0 else max(mean_edge, 1e-6)

            # Seed interior points on a coarse grid if polygon is large
            minx, maxx = float(np.min(poly_x)), float(np.max(poly_x))
            miny, maxy = float(np.min(poly_y)), float(np.max(poly_y))
            nx = max(0, int(math.floor((maxx - minx) / max(spacing_m, 1e-6))) - 1)
            ny = max(0, int(math.floor((maxy - miny) / max(spacing_m, 1e-6))) - 1)
            steiner: List[Tuple[float, float]] = []
            if nx > 0 and ny > 0:
                for ix in range(1, nx + 1):
                    gx = minx + ix * (maxx - minx) / (nx + 1)
                    for iy in range(1, ny + 1):
                        gy = miny + iy * (maxy - miny) / (ny + 1)
                        if point_in_polygon(gx, gy, poly_x, poly_y):
                            steiner.append((gx, gy))

            # Build combined point set for triangulation
            pts_x = np.concatenate([poly_x, np.array([s[0] for s in steiner], dtype=float)]) if steiner else poly_x
            pts_y = np.concatenate([poly_y, np.array([s[1] for s in steiner], dtype=float)]) if steiner else poly_y

            # Triangulate
            try:
                from scipy.spatial import Delaunay
                tri = Delaunay(np.column_stack([pts_x, pts_y]))
                simplices = tri.simplices
            except Exception:
                simplices = np.empty((0, 3), dtype=int)

            # Keep triangles whose centroids are inside polygon
            kept: List[List[int]] = []
            if simplices.size > 0:
                cx = np.mean(pts_x[simplices], axis=1)
                cy = np.mean(pts_y[simplices], axis=1)
                mask = np.array([point_in_polygon(cx[i], cy[i], poly_x, poly_y) for i in range(len(cx))])
                for tri_idx in simplices[mask]:
                    kept.append(tri_idx.tolist())

            # Update nodes df to include interior points (Steiner)
            if steiner:
                # Convert back to lon/lat for steiner points
                R = 6378206.4
                lon0r = math.radians(lon0)
                lat0r = math.radians(lat0)
                def xy2ll(x, y):
                    return math.degrees(lon0r + x / (R * math.cos(lat0r))), math.degrees(y / R)
                new_lonlat = [xy2ll(s[0], s[1]) for s in steiner]
                if len(new_lonlat) > 0:
                    dp_max = float(ring_df['value_1'].max()) if 'value_1' in ring_df.columns else 0.0
                    add_df = pd.DataFrame({'x': [p[0] for p in new_lonlat],
                                           'y': [p[1] for p in new_lonlat],
                                           'value_1': [dp_max] * len(new_lonlat),
                                           'gid': [None] * len(new_lonlat)})
                    ring_df = pd.concat([ring_df, add_df], axis=0, ignore_index=True)

            updated_nodes_list.append(ring_df)
            elements_grouped.append(kept)

        return updated_nodes_list, elements_grouped


class ChannelMeshBuilder:
    """Build an AdcircMesh from a list of ChannelStrip objects."""

    @staticmethod
    def build_mesh(strips: List[ChannelStrip], center_depth_is_elevation: bool, 
                   junction_nodes: List[pd.DataFrame] = None, 
                   junction_elements: List[List[List[int]]] = None) -> AdcircMesh:
        nodes_list: List[pd.DataFrame] = []
        elements_list: List[List[int]] = []
        node_offset = 0
        # Global gid registry
        gid_to_index: Dict[int, int] = {}
        
        # Build channel strip nodes and elements first
        junction_node_coords = []  # Store junction node coordinates for deduplication
        if junction_nodes:
            for junc_nodes in junction_nodes:
                junction_node_coords.extend([(row['x'], row['y']) for _, row in junc_nodes.iterrows()])
        for strip in strips:
            n = strip.left_lonlat.shape[0]
            if n < 2:
                continue
            left = strip.left_lonlat
            right = strip.right_lonlat
            # If mid nodes present, append them to nodes and build per-cross-section index mapping
            mid_idx_per_i: List[List[int]] = [[] for _ in range(n)]
            if strip.mid_lonlat is not None:
                # Prepare flattened dataframe with cross-section index and GIDs
                mid_records = []
                for i in range(n):
                    pts = strip.mid_lonlat[i]
                    if pts is None:
                        continue
                    gids_i = strip.mid_gids[i] if strip.mid_gids and strip.mid_gids[i] is not None else [None] * len(pts)
                    for k in range(len(pts)):
                        mid_records.append({
                            'sec_i': i,
                            'x': float(pts[k][0]),
                            'y': float(pts[k][1]),
                            'value_1': float(strip.mid_depths[i][k]) if strip.mid_depths and strip.mid_depths[i] is not None else 0.0,
                            'gid': None if gids_i is None or k >= len(gids_i) else gids_i[k]
                        })
                if mid_records:
                    mid_df = pd.DataFrame(mid_records)
                    if not mid_df.empty:
                        # Assign provisional indices
                        start_idx = node_offset + 1
                        mid_df.index = range(start_idx, start_idx + len(mid_df))
                        # Reuse indices by gid where possible
                        if 'gid' in mid_df.columns:
                            for idx, row in list(mid_df.iterrows()):
                                gid = row['gid']
                                if gid is None or (isinstance(gid, float) and np.isnan(gid)):
                                    continue
                                if gid in gid_to_index:
                                    # Repoint to existing node index
                                    mid_df.rename(index={idx: gid_to_index[gid]}, inplace=True)
                                else:
                                    gid_to_index[gid] = idx
                        # Remove duplicate indices keeping first
                        mid_df = mid_df.loc[~mid_df.index.duplicated(keep='first')]
                        # Build per-cross-section index mapping from sec_i
                        for i in range(n):
                            idx_i = mid_df.index[mid_df['sec_i'] == i].tolist()
                            if idx_i:
                                mid_idx_per_i[i] = idx_i
                        # Append to nodes
                        nodes_list.append(mid_df[['x', 'y', 'value_1']])
                        node_offset = int(max(node_offset, int(max(mid_df.index))))
            # Prepare gid arrays if provided
            left_gid_arr = strip.left_gid if strip.left_gid is not None else np.array([None] * n, dtype=object)
            right_gid_arr = strip.right_gid if strip.right_gid is not None else np.array([None] * n, dtype=object)
            # Node ordering: left[0..n-1], right[0..n-1]
            nodes = np.vstack([left, right])
            gids = np.concatenate([left_gid_arr, right_gid_arr])
            # Elevations (value_1)
            depths = strip.center_depths
            if depths.shape[0] != n:
                depths = np.interp(np.arange(n), np.linspace(0, n - 1, len(strip.center_depths)), strip.center_depths)
            if center_depth_is_elevation:
                z_left = depths.copy()
                z_right = depths.copy()
            else:
                z_left = -depths.copy()
                z_right = -depths.copy()
            z = np.concatenate([z_left, z_right])
            # Build strip node dataframe with gid column
            df = pd.DataFrame({'x': nodes[:, 0], 'y': nodes[:, 1], 'value_1': z, 'gid': gids})
            # Assign indices deterministically, but we will attempt gid reuse
            df.index = range(node_offset + 1, node_offset + 1 + df.shape[0])
            # Reuse indices by gid where possible
            for rel_idx, (idx, row) in enumerate(df.iterrows()):
                gid = row['gid']
                if gid is None or (isinstance(gid, float) and np.isnan(gid)):
                    continue
                if gid in gid_to_index:
                    # Repoint this row's index to existing node index
                    df.rename(index={idx: gid_to_index[gid]}, inplace=True)
                else:
                    gid_to_index[gid] = idx
            # After potential index merges, ensure no duplicate indices in concatenation
            # Do NOT keep 'gid' in the final mesh nodes table
            nodes_list.append(df.loc[~df.index.duplicated(keep='first'), ['x', 'y', 'value_1']])
            # Elements: stitch with mid nodes where available; otherwise two bank-only triangles
            for i in range(n - 1):
                li = df.index[i]
                li1 = df.index[i + 1]
                ri = df.index[n + i]
                ri1 = df.index[n + i + 1]
                Ai = [li] + (mid_idx_per_i[i] if i < len(mid_idx_per_i) and mid_idx_per_i[i] else []) + [ri]
                Bi = [li1] + (mid_idx_per_i[i + 1] if (i + 1) < len(mid_idx_per_i) and mid_idx_per_i[i + 1] else []) + [ri1]
                if len(Ai) == 2 and len(Bi) == 2:
                    # No mid nodes on either cross-section; use two bank-only triangles
                    elements_list.append([li, li1, ri1])
                    elements_list.append([li, ri1, ri])
                else:
                    # Greedy triangulation between two polylines Ai and Bi
                    a = 0
                    b = 0
                    na = len(Ai)
                    nb = len(Bi)
                    while a < na - 1 or b < nb - 1:
                        if a == na - 1:
                            # Advance B
                            elements_list.append([Ai[a], Bi[b], Bi[b + 1]])
                            b += 1
                            continue
                        if b == nb - 1:
                            # Advance A
                            elements_list.append([Ai[a], Ai[a + 1], Bi[b]])
                            a += 1
                            continue
                        # Compare normalized next steps to interleave fairly
                        next_a = (a + 1) / max(1, (na - 1))
                        next_b = (b + 1) / max(1, (nb - 1))
                        if next_a <= next_b:
                            elements_list.append([Ai[a], Ai[a + 1], Bi[b]])
                            a += 1
                        else:
                            elements_list.append([Ai[a], Bi[b], Bi[b + 1]])
                            b += 1
            node_offset = max(node_offset, int(max(df.index)))
        if not nodes_list:
            empty_nodes = pd.DataFrame(columns=['x', 'y', 'value_1'])
            empty_elems = pd.DataFrame(columns=[0, 1, 2, 3])
            return AdcircMesh(nodes=empty_nodes, elements=empty_elems)
        
        nodes_df = pd.concat(nodes_list, axis=0)
        
        # Fix: Remove duplicate indices after concatenation (keep first occurrence)
        # This prevents AdcircMesh from renumbering duplicates which creates phantom nodes
        if not nodes_df.index.is_unique:
            num_before = len(nodes_df)
            nodes_df = nodes_df.loc[~nodes_df.index.duplicated(keep='first')]
            num_after = len(nodes_df)
            num_removed = num_before - num_after
            if num_removed > 0:
                print(f"Removed {num_removed} duplicate node indices from concatenated node list")
        
        # Ensure only x, y, value_1 columns exist for AdcircMesh
        nodes_df = nodes_df[['x', 'y', 'value_1']]
        elements_df = pd.DataFrame(elements_list, columns=[1, 2, 3])
        # Insert nNodes-per-element column (0) with value 3 for triangles
        elements_df.insert(0, 0, 3)
        elements_df.index = range(1, len(elements_list) + 1)
        
        # Create initial mesh
        mesh = AdcircMesh(nodes=nodes_df, elements=elements_df)
        
        # Add junction elements if provided, with node deduplication
        if junction_nodes and junction_elements:
            mesh = ChannelMeshBuilder._add_junction_elements_by_gid(
                mesh, junction_nodes, junction_elements, gid_to_index
            )
        # Ensure node indices are incremental and element connectivity updated
        mesh = ChannelMeshBuilder._renumber_nodes(mesh)

        return mesh

    @staticmethod
    def _add_junction_elements_by_gid(mesh: AdcircMesh, junction_nodes: List[pd.DataFrame], 
                                      junction_elements: List[List[List[int]]], gid_to_index: Dict[int, int]) -> AdcircMesh:
        """Add junction elements by reusing nodes via gid mapping (grouped per ring)."""
        existing_nodes = mesh.nodes.copy()
        existing_elements = mesh.elements.elements.copy()
        new_nodes = existing_nodes[['x', 'y', 'value_1']].copy()
        next_index = int(existing_nodes.index.max()) if len(existing_nodes) else 0

        # Map each ring's nodes to mesh indices
        rings_indices: List[List[int]] = []
        for junc_nodes in junction_nodes:
            if junc_nodes is None or junc_nodes.empty:
                rings_indices.append([])
                continue
            gids = junc_nodes['gid'].tolist() if 'gid' in junc_nodes.columns else [None] * len(junc_nodes)
            node_indices: List[int] = []
            for (_, row), gid in zip(junc_nodes.iterrows(), gids):
                if gid is not None and gid in gid_to_index:
                    node_indices.append(gid_to_index[gid])
                else:
                    next_index += 1
                    node_indices.append(next_index)
                    new_nodes.loc[next_index, ['x', 'y', 'value_1']] = [row['x'], row['y'], row['value_1']]
                    if gid is not None:
                        gid_to_index[gid] = next_index
            rings_indices.append(node_indices)

        # Build new elements per ring
        new_elements_list = []
        for ring_idx, elems in enumerate(junction_elements or []):
            if not elems:
                continue
            ring_map = rings_indices[ring_idx]
            for elem in elems:
                mapped = [ring_map[i] for i in elem]
                new_elements_list.append([3] + mapped)

        if new_elements_list:
            new_elements_df = pd.DataFrame(new_elements_list, columns=[0, 1, 2, 3])
            new_elements_df.index = range(len(existing_elements) + 1, len(existing_elements) + len(new_elements_list) + 1)
            updated_elements = pd.concat([existing_elements, new_elements_df])
        else:
            updated_elements = existing_elements

        return AdcircMesh(nodes=new_nodes, elements=updated_elements)

    @staticmethod
    def _renumber_nodes(mesh: AdcircMesh) -> AdcircMesh:
        """Renumber nodes to 1..N and remap element connectivity accordingly."""
        nodes_df = mesh.nodes.copy()
        elems_df = mesh.elements.elements.copy()
        old_indices = list(nodes_df.index)
        new_indices = list(range(1, len(old_indices) + 1))
        mapping = dict(zip(old_indices, new_indices))
        nodes_df.index = new_indices
        for col in [1, 2, 3]:
            if col in elems_df.columns:
                elems_df[col] = elems_df[col].map(mapping).astype(int)
        return AdcircMesh(nodes=nodes_df[['x', 'y', 'value_1']], elements=elems_df)

    @staticmethod
    def validate_mesh_nodes(mesh: AdcircMesh, tolerance: float = 1e-12) -> Dict:
        """
        Validate mesh for duplicated nodes and report findings.
        
        Args:
            mesh: AdcircMesh to validate
            tolerance: Coordinate tolerance for considering nodes as duplicates (meters)
            
        Returns:
            Dictionary with validation results including duplicate information
        """
        nodes_df = mesh.nodes.copy()
        
        if len(nodes_df) == 0:
            return {
                'has_duplicates': False,
                'duplicate_groups': [],
                'total_duplicates': 0,
                'unique_nodes': 0
            }
        
        # Convert to Cartesian coordinates for distance calculations
        # Use the first node as reference for local Cartesian projection
        ref_lon = float(nodes_df['x'].iloc[0])
        ref_lat = float(nodes_df['y'].iloc[0])
        
        R = 6378206.4  # Earth radius in meters
        lon0r = math.radians(ref_lon)
        lat0r = math.radians(ref_lat)
        
        def ll2xy(lon, lat):
            lonr = math.radians(lon)
            latr = math.radians(lat)
            x = R * (lonr - lon0r) * math.cos(lat0r)
            y = R * (latr - 0.0)
            return x, y
        
        # Convert all nodes to Cartesian
        coords_cart = []
        for _, row in nodes_df.iterrows():
            x, y = ll2xy(row['x'], row['y'])
            coords_cart.append([x, y])
        
        coords_cart = np.array(coords_cart)
        
        # Find duplicate groups using spatial clustering
        from scipy.spatial import cKDTree
        
        tree = cKDTree(coords_cart)
        duplicate_groups = []
        processed = set()
        
        for i in range(len(coords_cart)):
            if i in processed:
                continue
                
            # Find all nodes within tolerance of this node
            neighbors = tree.query_ball_point(coords_cart[i], tolerance)
            
            if len(neighbors) > 1:
                # This is a duplicate group
                group = []
                for j in neighbors:
                    if j not in processed:
                        node_idx = nodes_df.index[j]
                        group.append({
                            'node_index': int(node_idx),
                            'lon': float(nodes_df.loc[node_idx, 'x']),
                            'lat': float(nodes_df.loc[node_idx, 'y']),
                            'elevation': float(nodes_df.loc[node_idx, 'value_1']),
                            'cartesian_x': coords_cart[j][0],
                            'cartesian_y': coords_cart[j][1]
                        })
                        processed.add(j)
                
                if len(group) > 1:
                    duplicate_groups.append(group)
            else:
                processed.add(i)
        
        # Calculate statistics
        total_duplicates = sum(len(group) for group in duplicate_groups)
        unique_nodes = len(nodes_df) - total_duplicates + len(duplicate_groups)
        
        return {
            'has_duplicates': len(duplicate_groups) > 0,
            'duplicate_groups': duplicate_groups,
            'total_duplicates': total_duplicates,
            'unique_nodes': unique_nodes,
            'total_nodes': len(nodes_df),
            'tolerance_meters': tolerance
        }


class ChannelMeshGeneratorApp:
    """Main application class for channel mesh generation."""

    def __init__(self):
        self.flowline_reader = FlowlineReader()
        self.centerline_processor = CenterlineProcessor()
        self.strip_generator = ChannelStripGenerator()

    def generate_channel_mesh(self, flowline_file: str,
                            output_file: str,
                            channel_spacing: float = 100.0,
                            zshift: float = 0.0,
                            effective_width_ratio: float = 1.0,
                            max_width: float = np.inf,
                            min_width: float = 0.0,
                            smoothing_span: int = 0,
                            min_spacing: float = 10.0,
                            radius_to_merge_shppoints: float = 10.0,
                            split_radius_m: float = 2.0,
                            depth_is_elevation: bool = False,
                            keep_bend_node_degree: Optional[float] = None,
                            junction_midnode_spacing: Optional[float] = None,
                            channel_midnode_spacing_factor: Optional[float] = None,
                            channel_midnode_spacing: Optional[float] = None,
                            duplicate_tolerance: float = 1e-6) -> AdcircMesh:
        """
        Generate channel mesh from flowline file.

        Args:
            flowline_file: Path to flowline file (GeoJSON or Shapefile)
            output_file: Path for output mesh file
            channel_spacing: Base spacing for centerline nodes
            zshift: Vertical shift for depths
            effective_width_ratio: Width scaling factor
            max_width: Maximum channel width
            min_width: Minimum channel width
            smoothing_span: Smoothing span for depth/width values
            min_spacing: Minimum node spacing
            radius_to_merge_shppoints: Merge centerline segments whose endpoints are within this radius
            split_radius_m: Split when an endpoint lies within this radius of an interior point
            depth_is_elevation: Interpret depth values as bed elevation
            keep_bend_node_degree: Preserve nodes with interior angle >= this threshold
            junction_midnode_spacing: Spacing in meters for junction chord mid-nodes (defaults to channel midnode spacing if either channel option is specified)
            channel_midnode_spacing_factor: Multiplier for channel mid-node spacing (mutually exclusive with channel_midnode_spacing, default: 1.0)
            channel_midnode_spacing: Absolute spacing in meters for channel mid-nodes (mutually exclusive with channel_midnode_spacing_factor)
            duplicate_tolerance: Tolerance in meters for detecting duplicated nodes

        Returns:
            Generated AdcircMesh object
        """
        print(f"Reading flowline data from: {flowline_file}")
        centerlines = self.flowline_reader.read_flowline_file(
            flowline_file, zshift=zshift,
            effective_width_ratio=effective_width_ratio,
            max_width=max_width, min_width=min_width
        )

        print(f"Processing {len(centerlines)} centerline sequences...")        
        processed_centerlines = self.centerline_processor.process_centerlines(centerlines, smoothing_span=smoothing_span)

        # Global reference lat/lon for local Cartesian
        all_lon = np.concatenate([cl.lonlat[:, 0] for cl in processed_centerlines]) if processed_centerlines else np.array([0.0])
        all_lat = np.concatenate([cl.lonlat[:, 1] for cl in processed_centerlines]) if processed_centerlines else np.array([0.0])
        lon0 = float(np.mean(all_lon))
        lat0 = float(np.mean(all_lat))

        # Split where an endpoint lies near an interior point of another sequence
        if split_radius_m and split_radius_m > 0:
            before = len(processed_centerlines)
            processed_centerlines = self.centerline_processor.split_at_near_endpoints(
                processed_centerlines, split_radius_m, lon0, lat0
            )
            after = len(processed_centerlines)
            if after != before:
                print(f"Split sequences at near-endpoint intersections (radius={split_radius_m} m): {before} -> {after}")

        # Merge centerline segments at close endpoints (endpoint vs endpoint)
        if radius_to_merge_shppoints and radius_to_merge_shppoints > 0:
            print(f"Merging centerline segments within {radius_to_merge_shppoints} m...")
            processed_centerlines = self.centerline_processor.merge_endpoints(
                processed_centerlines, radius_to_merge_shppoints, lon0, lat0
            )

        # Prepass: compute junction endpoint angles on processed (unresampled) centerlines
        pre_items: List[Dict] = []
        for si, cl in enumerate(processed_centerlines):
            pre_items.append({'lonlat': cl.lonlat, 'widths': cl.widths, 'depths': cl.depths, 'seq_idx': si})
        pre_overrides, _, _, pre_junction_endpoints = EndpointCrosspointBuilder().compute_overrides(pre_items, lon0, lat0, radius_to_merge_shppoints)

        # Resample centerlines
        resampled: List[Dict] = []
        bend_reports = []
        for si, cl in enumerate(processed_centerlines):
            # Endpoint angles for spacing (radians), only for junction endpoints
            pre_start_pair = pre_overrides.get((si, 0)) if (si, 0) in pre_junction_endpoints else None
            pre_end_pair = pre_overrides.get((si, 1)) if (si, 1) in pre_junction_endpoints else None
            # Compute angles relative to crosspoint chord normals for spacing prepass
            def endpoint_angle_from_pair(lonlat_seq: np.ndarray, pair: Optional[Tuple[np.ndarray, np.ndarray]], at_start: bool) -> Optional[float]:
                if pair is None or lonlat_seq.shape[0] < 2:
                    return None
                # direction at endpoint
                p0 = lonlat_seq[0]
                p1 = lonlat_seq[1]
                pnm1 = lonlat_seq[-2]
                pn = lonlat_seq[-1]
                if at_start:
                    vx, vy = p1[0] - p0[0], p1[1] - p0[1]
                    ex, ey = p0[0], p0[1]
                else:
                    vx, vy = pn[0] - pnm1[0], pn[1] - pnm1[1]
                    ex, ey = pn[0], pn[1]
                # Convert to Cartesian
                R = ChannelStripGenerator.EARTH_RADIUS_M
                lon0r = math.radians(lon0)
                lat0r = math.radians(lat0)
                def ll2xy(pt):
                    return R * (math.radians(pt[0]) - lon0r) * math.cos(lat0r), R * (math.radians(pt[1]) - 0.0)
                vxn, vyn = ll2xy((ex + vx, ey + vy))
                vex, vey = ll2xy((ex, ey))
                vxm, vym = vxn - vex, vyn - vey
                nv = math.hypot(vxm, vym)
                if nv == 0.0:
                    return 0.0
                vxm /= nv
                vym /= nv
                cp1, cp2 = pair
                c1x, c1y = ll2xy(cp1)
                c2x, c2y = ll2xy(cp2)
                cx, cy = c2x - c1x, c2y - c1y
                cn = math.hypot(cx, cy)
                if cn == 0.0:
                    return 0.0
                nx1, ny1 = -cy / cn, cx / cn
                nx2, ny2 = -nx1, -ny1
                a1 = math.acos(max(-1.0, min(1.0, vxm * nx1 + vym * ny1)))
                a2 = math.acos(max(-1.0, min(1.0, vxm * nx2 + vym * ny2)))
                return min(a1, a2)

            ang_start = endpoint_angle_from_pair(cl.lonlat, pre_start_pair, at_start=True) if pre_start_pair is not None else None
            ang_end = endpoint_angle_from_pair(cl.lonlat, pre_end_pair, at_start=False) if pre_end_pair is not None else None

            lonlat_res, widths_res, depths_res = self.strip_generator.resample_centerline(
                cl.lonlat, cl.widths, cl.depths, channel_spacing, lon0, lat0,
                keep_bend_node_degree=keep_bend_node_degree,
                endpoint_angle_start_rad=ang_start,
                endpoint_angle_end_rad=ang_end
            )
            if np.any(np.isnan(widths_res)):
                widths_res = np.where(np.isnan(widths_res), max(min_width, 1.0), widths_res)
            if np.any(np.isnan(depths_res)):
                depths_res = np.where(np.isnan(depths_res), 0.0, depths_res)
            # Junction endpoint flags for this sequence
            start_is_junction = ((si, 0) in junction_endpoints) if 'junction_endpoints' in locals() else False
            end_is_junction = ((si, 1) in junction_endpoints) if 'junction_endpoints' in locals() else False
            resampled.append({'lonlat': lonlat_res, 'widths': widths_res, 'depths': depths_res, 'seq_idx': si,
                              'start_is_junction': start_is_junction, 'end_is_junction': end_is_junction})

            # Record sharp-bend nodes for reporting
            if keep_bend_node_degree is not None and keep_bend_node_degree > 0.0:
                bends = self.strip_generator.detect_bend_nodes(cl.lonlat, lon0, lat0, keep_bend_node_degree)
                for b in bends:
                    bend_reports.append({'seq_idx': si, **b})

        # Compute crosspoint overrides at merger endpoints
        cross_builder = EndpointCrosspointBuilder()
        overrides, junction_crosspoints, crosspoint_gids, junction_endpoints = cross_builder.compute_overrides(resampled, lon0, lat0, radius_to_merge_shppoints)
        
        # Optionally add mid nodes along junction chords before element generation
        # Default junction_midnode_spacing to channel midnode spacing if either channel option is specified
        if junction_midnode_spacing is not None:
            chord_spacing = junction_midnode_spacing
        elif channel_midnode_spacing is not None or channel_midnode_spacing_factor is not None:
            # Use channel midnode spacing when either channel option is specified
            if channel_midnode_spacing is not None:
                chord_spacing = channel_midnode_spacing
            else:
                # Use channel_spacing * factor when only factor is specified
                chord_spacing = channel_spacing * (channel_midnode_spacing_factor if channel_midnode_spacing_factor is not None else 1.0)
        else:
            # Fall back to channel_spacing (original behavior when no midnode options specified)
            chord_spacing = channel_spacing
        if junction_crosspoints and chord_spacing and chord_spacing > 0:
            # Augment junction_nodes/junction_elements by adding edge mid nodes
            aug_nodes = []
            group_ids_order: List[int] = []
            # start gid counter after existing crosspoint gids
            existing_gids: List[int] = []
            if 'crosspoint_gids' in locals() and crosspoint_gids:
                for pair in crosspoint_gids.values():
                    for g in pair:
                        if g is not None:
                            existing_gids.append(int(g))
            gid_counter_local = (max(existing_gids) + 1) if existing_gids else 1
            aug_elements = []
            total_mid_nodes_added = 0  # Track if any mid-nodes were actually added
            for group_id, cps in junction_crosspoints.items():
                group_ids_order.append(group_id)
                if len(cps) < 2:
                    continue
                # Build nodes df for this group including mid nodes along each chord
                group_nodes = []
                # depths per cp
                depths = [cp['depth'] for cp in cps]
                dp_max = max(depths) if depths else 0.0
                # Cartesian helpers
                R = ChannelStripGenerator.EARTH_RADIUS_M
                lon0r = math.radians(lon0)
                lat0r = math.radians(lat0)
                def ll2xy(lon, lat):
                    return R * (math.radians(lon) - lon0r) * math.cos(lat0r), R * (math.radians(lat) - 0.0)
                def xy2ll(x, y):
                    return math.degrees(lon0r + x / (R * math.cos(lat0r))), math.degrees(y / R)
                # Iterate chords
                for j in range(len(cps)):
                    j1 = j
                    j2 = (j + 1) % len(cps)
                    lon1, lat1 = cps[j1]['lon'], cps[j1]['lat']
                    lon2, lat2 = cps[j2]['lon'], cps[j2]['lat']
                    x1, y1 = ll2xy(lon1, lat1)
                    x2, y2 = ll2xy(lon2, lat2)
                    dx, dy = x2 - x1, y2 - y1
                    wid = math.hypot(dx, dy)
                    nn_add = int(math.floor(wid / max(chord_spacing, 1e-6)))
                    # Always push the chord start cp node on first chord
                    if j == 0:
                        group_nodes.append({'x': lon1, 'y': lat1, 'value_1': dp_max, 'gid': cps[j1].get('gid')})
                    if nn_add >= 1:
                        total_mid_nodes_added += nn_add  # Count mid-nodes added
                        for l in range(1, nn_add + 1):
                            xm = x1 + dx * (l / (nn_add + 1))
                            ym = y1 + dy * (l / (nn_add + 1))
                            lom, lam = xy2ll(xm, ym)
                            group_nodes.append({'x': lom, 'y': lam, 'value_1': dp_max, 'gid': gid_counter_local})
                            gid_counter_local += 1
                    # Push the chord end cp node
                    group_nodes.append({'x': lon2, 'y': lat2, 'value_1': dp_max, 'gid': cps[j2].get('gid')})
                # Deduplicate consecutive duplicates
                dedup = []
                for node in group_nodes:
                    if not dedup or (abs(dedup[-1]['x'] - node['x']) > 1e-12 or abs(dedup[-1]['y'] - node['y']) > 1e-12):
                        dedup.append(node)
                nodes_df = pd.DataFrame(dedup)
                aug_nodes.append(nodes_df)
                # Elements will be generated by JunctionMeshGenerator from this ring later; we pass nodes only
            # Only set junction_nodes if mid-nodes were actually added
            # If no mid-nodes added, fall back to using raw junction_crosspoints
            if aug_nodes and total_mid_nodes_added > 0:
                junction_nodes = aug_nodes
                junction_group_ids = group_ids_order
                print(f"Added {total_mid_nodes_added} chord mid-nodes to junction rings")

        # Generate junction elements (augmented rings if present)
        junction_generator = JunctionMeshGenerator()
        if 'junction_nodes' in locals() and junction_nodes:
            # Use ring nodes, perform constrained triangulation with optional interior points
            junction_nodes, junction_elements_grouped = junction_generator.generate_junction_elements_from_nodes(
                junction_nodes, lon0, lat0, interior_spacing=chord_spacing
            )
        else:
            j_nodes, j_elems = junction_generator.generate_junction_elements(junction_crosspoints, lon0, lat0)
            junction_nodes = j_nodes
            # Normalize to grouped format (one group per junction)
            junction_elements_grouped = [[elem] for elem in j_elems]
        
        total_junc_elems = sum(len(g) for g in junction_elements_grouped) if junction_nodes else 0
        if total_junc_elems:
            print(f"Generated {total_junc_elems} junction elements at {len(junction_nodes)} junctions")

        # Build strips, applying endpoint crosspoint overrides where present
        strips: List[ChannelStrip] = []
        for item in resampled:
            lonlat_res = item['lonlat']
            widths_res = item['widths']
            depths_res = item['depths']
            seq_idx = item['seq_idx']
            # Determine junction membership now that junction_endpoints is known
            start_is_junction = (seq_idx, 0) in junction_endpoints
            end_is_junction = (seq_idx, 1) in junction_endpoints
            # Prepare endpoint crosspoint pairs for angle calc at junction endpoints
            start_pair = overrides.get((seq_idx, 0)) if start_is_junction else None
            end_pair = overrides.get((seq_idx, 1)) if end_is_junction else None
            # Compute angles at resampled points including junction endpoints
            angles = self.strip_generator.compute_angles_for_resampled(
                lonlat_res, lon0, lat0,
                start_is_junction=start_is_junction,
                end_is_junction=end_is_junction,
                start_override=start_pair, end_override=end_pair
            )
            left, right = self.strip_generator.make_strip(lonlat_res, widths_res, lon0, lat0, default_width=max(min_width, 1.0), angles=angles)
            # Add channel mid nodes across wide sections
            mid_lonlat_list, mid_depths_list, mid_gids_list = self.strip_generator.add_channel_midnodes(
                left, right, depths_res, lon0, lat0,
                spacing_factor=channel_midnode_spacing_factor,
                spacing_meters=channel_midnode_spacing,
                start_is_junction=start_is_junction,
                end_is_junction=end_is_junction
            )

            # Apply overrides
            start_key = (seq_idx, 0)
            end_key = (seq_idx, 1)
            # Prepare gid arrays (None by default)
            left_gid = np.array([None] * left.shape[0], dtype=object)
            right_gid = np.array([None] * right.shape[0], dtype=object)

            # If junction chord mid-nodes were created for this endpoint, graft their GIDs into first/last span mid-lists
            def inject_chord_mid_gids(pair_key: Tuple[int, int], at_start: bool):
                nonlocal mid_gids_list, mid_lonlat_list, mid_depths_list
                if 'junction_group_ids' not in locals() or 'junction_nodes' not in locals():
                    return
                if pair_key not in junction_endpoints:
                    return
                grp_id = junction_endpoints[pair_key]
                if grp_id not in junction_group_ids:
                    return
                ring_idx = junction_group_ids.index(grp_id)
                ring_df = junction_nodes[ring_idx] if junction_nodes and ring_idx < len(junction_nodes) else None
                if ring_df is None or ring_df.empty or 'gid' not in ring_df.columns:
                    return
                gids_ring = ring_df['gid'].tolist()
                # identify this endpoint's two crosspoint gids (left,right)
                cp_pair_gids = crosspoint_gids.get(pair_key)
                if not cp_pair_gids:
                    return
                left_gid_cp, right_gid_cp = cp_pair_gids
                try:
                    i_left = gids_ring.index(left_gid_cp)
                    i_right = gids_ring.index(right_gid_cp)
                except ValueError:
                    return
                nring = len(gids_ring)

                def arc_indices(start_idx: int, end_idx: int, step: int) -> List[int]:
                    res = []
                    cur = (start_idx + step) % nring
                    while cur != end_idx:
                        res.append(cur)
                        cur = (cur + step) % nring
                    return res

                forward_idxs = arc_indices(i_left, i_right, +1)
                backward_idxs = arc_indices(i_left, i_right, -1)

                # choose the shorter non-empty arc (chord mid nodes should lie on that path)
                candidate_arcs = []
                if forward_idxs:
                    candidate_arcs.append((len(forward_idxs), forward_idxs, 'forward'))
                if backward_idxs:
                    candidate_arcs.append((len(backward_idxs), backward_idxs, 'backward'))
                if not candidate_arcs:
                    print(f"INFO chord seq={seq_idx} end={'start' if at_start else 'end'}: no candidate arcs (ring size={nring})")
                    return
                candidate_arcs.sort(key=lambda x: x[0])
                idxs = None
                mid_gids = None
                chosen_dir = None
                for arc_len, idxs_candidate, direction in candidate_arcs:
                    gids_candidate = [gids_ring[k] for k in idxs_candidate if gids_ring[k] is not None]
                    if gids_candidate:
                        idxs = idxs_candidate
                        mid_gids = gids_candidate
                        chosen_dir = direction
                        break
                if idxs is None or not mid_gids:
                    print(f"INFO chord seq={seq_idx} end={'start' if at_start else 'end'}: candidate arcs exist but mid_gids empty")
                    return
                # Decide which end of strip to attach
                if at_start:
                    idx_cross = 0
                else:
                    idx_cross = -1
                # Populate strip mid_gids for the adjacent cross-section index (0 or -1)
                if mid_gids_list is None:
                    # initialize local mid-node containers to allow gid injection
                    mid_gids_list = [None] * left.shape[0]
                if mid_lonlat_list is None:
                    mid_lonlat_list = [None] * left.shape[0]
                # If we have ring coordinates for these mid nodes, set them directly
                try:
                    mid_ll = ring_df.iloc[idxs][['x', 'y']].to_numpy(dtype=float)
                except Exception:
                    mid_ll = None
                try:
                    mid_dv = ring_df.iloc[idxs]['value_1'].to_numpy(dtype=float) if 'value_1' in ring_df.columns else None
                except Exception:
                    mid_dv = None
                if mid_ll is not None and len(mid_ll) == len(mid_gids):
                    mid_lonlat_list[idx_cross] = mid_ll.tolist()
                    if mid_depths_list is not None:
                        if mid_dv is not None and len(mid_dv) == len(mid_gids):
                            mid_depths_list[idx_cross] = mid_dv.tolist()
                        else:
                            use_depth = float(depths_res[0] if at_start else depths_res[-1]) if len(depths_res) else 0.0
                            mid_depths_list[idx_cross] = [use_depth] * len(mid_gids)
                    mid_gids_list[idx_cross] = list(mid_gids)
                else:
                    # Only set GIDs if we already have mid points and sizes match
                    mg = mid_gids_list[idx_cross] if mid_gids_list else None
                    if mg is None and mid_lonlat_list and mid_lonlat_list[idx_cross] is not None:
                        if len(mid_lonlat_list[idx_cross]) == len(mid_gids):
                            mid_gids_list[idx_cross] = list(mid_gids)

            if start_key in overrides:
                cp1, cp2 = overrides[start_key]
                gids = crosspoint_gids.get(start_key)
                left[0, :] = cp1
                right[0, :] = cp2
                if gids is not None:
                    left_gid[0] = gids[0]
                    right_gid[0] = gids[1]
                # Share chord mid-node GIDs with first span
                inject_chord_mid_gids(start_key, at_start=True)
            if end_key in overrides:
                cp1, cp2 = overrides[end_key]
                gids = crosspoint_gids.get(end_key)
                left[-1, :] = cp1
                right[-1, :] = cp2
                if gids is not None:
                    left_gid[-1] = gids[0]
                    right_gid[-1] = gids[1]
                # Share chord mid-node GIDs with last span
                inject_chord_mid_gids(end_key, at_start=False)

            strip_obj = ChannelStrip(left_lonlat=left, right_lonlat=right, center_depths=depths_res,
                                      left_gid=left_gid, right_gid=right_gid,
                                      mid_lonlat=mid_lonlat_list, mid_depths=mid_depths_list, mid_gids=mid_gids_list)
            setattr(strip_obj, 'seq_idx_debug', seq_idx)
            strips.append(strip_obj)

        print("Building mesh...")
        mesh = ChannelMeshBuilder.build_mesh(strips, center_depth_is_elevation=depth_is_elevation,
                                        junction_nodes=junction_nodes, junction_elements=junction_elements_grouped)
        # AdcircMesh.elements is an Elements wrapper; show count via underlying DataFrame
        elem_count = mesh.elements.elements.shape[0] if hasattr(mesh.elements, 'elements') else len(mesh.elements)
        print(f"Generated mesh with {len(mesh.nodes)} nodes and {elem_count} elements")

        # Validate mesh for duplicated nodes
        print("Validating mesh for duplicated nodes...")
        validation_results = ChannelMeshBuilder.validate_mesh_nodes(mesh, tolerance=duplicate_tolerance)
        
        if validation_results['has_duplicates']:
            print(f"WARNING: Found {len(validation_results['duplicate_groups'])} groups of duplicated nodes!")
            print(f"Total duplicate nodes: {validation_results['total_duplicates']}")
            print(f"Unique nodes: {validation_results['unique_nodes']}")
            print(f"Total nodes: {validation_results['total_nodes']}")
            print(f"Tolerance used: {validation_results['tolerance_meters']:.2e} meters")
            
            # Report details for each duplicate group
            for i, group in enumerate(validation_results['duplicate_groups']):
                print(f"\nDuplicate group {i+1} ({len(group)} nodes):")
                for j, node in enumerate(group):
                    print(f"  Node {j+1}: index={node['node_index']}, "
                          f"lon={node['lon']:.8f}, lat={node['lat']:.8f}, "
                          f"elevation={node['elevation']:.3f}")
            raise Exception("Duplicated nodes detected")
        else:
            print("No duplicated nodes detected.")

        if output_file:
            print(f"Writing mesh to: {output_file}")
            mesh.write(output_file, overwrite=True)

            # Also write bend report if requested
            if keep_bend_node_degree is not None and keep_bend_node_degree > 0.0 and bend_reports:
                import csv, os
                out_csv = os.path.splitext(output_file)[0] + "_bend_nodes.csv"
                with open(out_csv, 'w', newline='') as f:
                    w = csv.DictWriter(f, fieldnames=['seq_idx', 'pt_idx', 'lon', 'lat', 'angle_deg'])
                    w.writeheader()
                    for row in bend_reports:
                        w.writerow(row)
                print(f"Wrote bend-node report: {out_csv}")
            
            # Write duplicate nodes report if duplicates were found
            if validation_results['has_duplicates']:
                import csv, os
                out_csv = os.path.splitext(output_file)[0] + "_duplicate_nodes.csv"
                with open(out_csv, 'w', newline='') as f:
                    w = csv.DictWriter(f, fieldnames=['group_id', 'node_index', 'lon', 'lat', 'elevation', 'cartesian_x', 'cartesian_y'])
                    w.writeheader()
                    for group_id, group in enumerate(validation_results['duplicate_groups']):
                        for node in group:
                            w.writerow({
                                'group_id': group_id + 1,
                                'node_index': node['node_index'],
                                'lon': node['lon'],
                                'lat': node['lat'],
                                'elevation': node['elevation'],
                                'cartesian_x': node['cartesian_x'],
                                'cartesian_y': node['cartesian_y']
                            })
                print(f"Wrote duplicate-nodes report: {out_csv}")

        return mesh


def get_parser():
    """Create argument parser for CLI."""
    parser = argparse.ArgumentParser(
        description='Generate ADCIRC channel mesh from flowline data',
        prog='vewutils channelpaving generate-channel-mesh',
        add_help=False
    )

    parser.add_argument(
        'flowline_file',
        help='Path to flowline file (GeoJSON or Shapefile)'
    )

    parser.add_argument(
        'output_file',
        help='Path for output channel mesh file (.grd)'
    )

    parser.add_argument(
        '--channel-spacing', '-s',
        type=float,
        default=100.0,
        help='Base spacing for centerline nodes (default: 100.0)'
    )

    parser.add_argument(
        '--zshift', '-z',
        type=float,
        default=0.0,
        help='Vertical shift to apply to depths (default: 0.0)'
    )

    parser.add_argument(
        '--width-ratio', '-w',
        type=float,
        default=1.0,
        help='Width scaling ratio (default: 1.0)'
    )

    parser.add_argument(
        '--max-width',
        type=float,
        default=np.inf,
        help='Maximum channel width (default: unlimited)'
    )

    parser.add_argument(
        '--min-width',
        type=float,
        default=0.0,
        help='Minimum channel width (default: 0.0)'
    )

    parser.add_argument(
        '--smoothing-span',
        type=int,
        default=0,
        help='Smoothing span for depth/width values (default: 0, no smoothing)'
    )

    parser.add_argument(
        '--min-spacing',
        type=float,
        default=10.0,
        help='Minimum node spacing (default: 10.0)'
    )

    parser.add_argument(
        '--radius-to-merge-shppoints',
        type=float,
        default=10.0,
        help='Merge centerline segments whose endpoints are within this radius (meters)'
    )

    parser.add_argument(
        '--split-radius',
        type=float,
        default=2.0,
        help='Split when an endpoint lies within this radius (meters) of an interior point on another segment'
    )

    parser.add_argument(
        '--depth-is-elevation',
        action='store_true',
        help='Interpret provided depth values as bed elevation (value_1). If not set, depths are treated as positive depths and converted to negative elevations.'
    )

    parser.add_argument(
        '--keep-bend-node-degree',
        type=float,
        default=None,
        help='If set (degrees), preserve original centerline nodes whose interior angle >= this threshold during resampling.'
    )

    parser.add_argument(
        '--junction-midnode-spacing',
        type=float,
        default=None,
        help='Spacing in meters used to add mid nodes along junction chords. If not specified, uses --channel-midnode-spacing if set, otherwise uses --channel-spacing * --channel-midnode-spacing-factor (default 1.0), or falls back to --channel-spacing.'
    )

    # Mutually exclusive group for channel midnode spacing
    midnode_group = parser.add_mutually_exclusive_group()
    midnode_group.add_argument(
        '--channel-midnode-spacing-factor',
        type=float,
        default=None,
        help='Multiplier applied to longitudinal step when inserting mid nodes across channel transects (default: 1.0 if neither option is specified).'
    )
    midnode_group.add_argument(
        '--channel-midnode-spacing',
        type=float,
        default=None,
        help='Absolute spacing in meters for inserting mid nodes across channel transects (mutually exclusive with --channel-midnode-spacing-factor).'
    )

    parser.add_argument(
        '--duplicate-tolerance',
        type=float,
        default=1e-6,
        help='Tolerance in meters for detecting duplicated nodes (default: 1e-6)'
    )

    return parser


def main(args=None):
    """Main function for CLI."""
    if args is None:
        parser = get_parser()
        args = parser.parse_args()

    app = ChannelMeshGeneratorApp()

    # Handle defaults for mutually exclusive channel midnode spacing options
    channel_midnode_spacing_factor = args.channel_midnode_spacing_factor
    channel_midnode_spacing = args.channel_midnode_spacing
    if channel_midnode_spacing_factor is None and channel_midnode_spacing is None:
        # Default to spacing factor of 1.0 if neither is specified
        channel_midnode_spacing_factor = 1.0

    try:
        mesh = app.generate_channel_mesh(
            flowline_file=args.flowline_file,
            output_file=args.output_file,
            channel_spacing=args.channel_spacing,
            zshift=args.zshift,
            effective_width_ratio=args.width_ratio,
            max_width=args.max_width,
            min_width=args.min_width,
            smoothing_span=args.smoothing_span,
            min_spacing=args.min_spacing,
            radius_to_merge_shppoints=args.radius_to_merge_shppoints,
            split_radius_m=args.split_radius,
            depth_is_elevation=args.depth_is_elevation,
            keep_bend_node_degree=args.keep_bend_node_degree,
            junction_midnode_spacing=args.junction_midnode_spacing,
            channel_midnode_spacing_factor=channel_midnode_spacing_factor,
            channel_midnode_spacing=channel_midnode_spacing,
            duplicate_tolerance=args.duplicate_tolerance
        )

        print("Channel mesh generation completed successfully!")
        return 0

    except Exception as e:
        print(f"Error generating channel mesh: {e}")
        return 1


if __name__ == '__main__':
    exit(main())
