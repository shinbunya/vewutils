#!/usr/bin/env python3
"""
Convert VEW string YAML files to GeoJSON or Esri Shapefile.
Output format is chosen from the output file extension (.geojson, .json, or .shp).
"""

import argparse
import json
import os
from typing import Dict, List

import geopandas as gpd
import yaml


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

        coordinates = []
        nodes = []

        for node in vewstring:
            if "x" not in node or "y" not in node:
                raise ValueError(
                    f"Node in vewstring {i} is missing 'x' or 'y' coordinates: {node}"
                )
            if "node_id" not in node:
                raise ValueError(f"Node in vewstring {i} is missing 'node_id': {node}")
            if "bank_elevation" not in node:
                raise ValueError(
                    f"Node in vewstring {i} is missing 'bank_elevation': {node}"
                )
            if "bank_mannings_n" not in node:
                raise ValueError(
                    f"Node in vewstring {i} is missing 'bank_mannings_n': {node}"
                )

            coordinates.append([float(node["x"]), float(node["y"])])

            nodes.append(
                {
                    "node_id": int(node["node_id"]),
                    "bank_elevation": float(node["bank_elevation"]),
                    "bank_mannings_n": float(node["bank_mannings_n"]),
                }
            )

        # GeoJSON LineString and GEOS require at least two positions; duplicate a lone vertex.
        if len(coordinates) == 1:
            coordinates = [coordinates[0], coordinates[0]]

        feature = {
            "type": "Feature",
            "geometry": {"type": "LineString", "coordinates": coordinates},
            "properties": {"nodes": nodes},
        }

        features.append(feature)

    return {"type": "FeatureCollection", "features": features}


def output_format_from_path(output_path: str) -> str:
    """Return 'geojson' or 'shapefile' based on file extension."""
    ext = os.path.splitext(output_path)[1].lower()
    if ext in (".geojson", ".json"):
        return "geojson"
    if ext == ".shp":
        return "shapefile"
    raise ValueError(
        f"Unsupported output extension {ext!r}; use .geojson, .json, or .shp"
    )


def _geojson_to_shapefile_features(geojson: Dict) -> List[Dict]:
    """Build GeoJSON features with shapefile-safe properties (nodes as JSON string)."""
    out = []
    for f in geojson["features"]:
        nodes_json = json.dumps(
            f["properties"]["nodes"], separators=(",", ":")
        )
        out.append(
            {
                "type": "Feature",
                "geometry": f["geometry"],
                "properties": {"nodes_json": nodes_json},
            }
        )
    return out


def write_geo_output(geojson: Dict, output_path: str) -> None:
    """Write GeoJSON dict to a file; format from ``output_path`` extension."""
    fmt = output_format_from_path(output_path)
    if fmt == "geojson":
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(geojson, f, indent=2, sort_keys=False)
        return

    features = _geojson_to_shapefile_features(geojson)
    gdf = gpd.GeoDataFrame.from_features(features, crs=None)
    gdf.to_file(output_path, driver="ESRI Shapefile")


def get_parser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description=(
            "Convert VEW string YAML to GeoJSON or Shapefile; "
            "format is inferred from the output path (.geojson, .json, or .shp)"
        ),
    )
    parser.add_argument(
        "input_yaml", help="Path to the input VEW string YAML file"
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output path (.geojson / .json for GeoJSON, .shp for Esri Shapefile)",
    )
    return parser


def main(args=None):
    if args is None:
        args = get_parser().parse_args()

    fmt = output_format_from_path(args.output)
    fmt_label = "GeoJSON" if fmt == "geojson" else "Esri Shapefile"

    print("=== VEW String YAML to Geo / Shapefile Converter ===")
    print(f"Input YAML file: {args.input_yaml}")
    print(f"Output ({fmt_label}): {args.output}")
    print()

    print("Loading VEW string YAML file...")
    with open(args.input_yaml, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    if "vewstrings" not in data:
        raise ValueError(f"File {args.input_yaml} does not contain 'vewstrings' key")

    vewstrings = data["vewstrings"]

    if not isinstance(vewstrings, list):
        raise ValueError(f"File {args.input_yaml}: 'vewstrings' must be a list")

    print(f"✓ YAML file loaded: {len(vewstrings)} VEW strings")
    print()

    print("Converting...")
    geojson = yaml_to_geojson(vewstrings)
    print(f"✓ Converted {len(geojson['features'])} VEW strings to features")
    print()

    print(f"Saving {fmt_label}...")
    write_geo_output(geojson, args.output)
    print(f"✓ Saved to: {args.output}")

    total_nodes = sum(
        len(feature["properties"]["nodes"]) for feature in geojson["features"]
    )
    print(f"  Total features: {len(geojson['features'])}")
    print(f"  Total nodes: {total_nodes}")
    print()
    print("=== Conversion Complete ===")

    return 0


if __name__ == "__main__":
    main()
