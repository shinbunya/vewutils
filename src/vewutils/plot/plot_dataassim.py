#!/usr/bin/env python3
"""Plot nodal DWLC fields from an ADCIRC dataassim.dat file on a fort.14 mesh.

An analysis-cycle file looks like::

    # comment
    <dt seconds>
    0.0
    ## optional snapshot label
    <node>,<value>    or    <node> <value>
    ...
    ##
    <node>,<value>
    ...
    ##

A line that is only ``##`` starts the next record. A ``##`` line that also
contains text, such as ``## NCSCv3, 2026-10-12 12:00:00, 2026-10-14 12:00:00``,
is a label for the record that follows and is not a record break.

Each record may list every mesh node or only a subset. Listed nodes use either
a comma or whitespace between the 1-based node id and the value. Nodes that
are omitted, and nodes stored as the fill value -99999, are left blank.

The first record is the correction at t = 0, and each later record is one
``dt`` later.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from adcircpy import AdcircMesh
from cartopy.io import shapereader
from matplotlib.collections import LineCollection
from shapely.geometry import box

DEFAULT_ANALYSIS = (
    "/scratch/sbunya/floodwater_output/NCSCv3/gfs/NCSCv3_gfs/archive/"
    "20261001/hour_00/adcirc/analysis"
)
FILL_VALUE = -99999.0
# (xmin, xmax, ymin, ymax) in degrees.
# nc: Cape Fear through the Outer Banks.
# gulf-to-ec: Gulf of Mexico through the U.S. East Coast.
COVERAGES = {
    "nc": (-79.0, -74.8, 33.4, 36.7),
    "gulf-to-ec": (-98.0, -65.0, 18.0, 46.0),
}
DEFAULT_COVERAGE = "nc"


def parse_node_value(text):
    """Return (node_id, value) from a comma-separated or whitespace-separated line."""
    if "," in text:
        node_text, value_text = (part.strip() for part in text.split(",", 1))
    else:
        parts = text.split()
        if len(parts) != 2:
            raise ValueError(
                f"expected 'node,value' or 'node value', got {text!r}"
            )
        node_text, value_text = parts
    return int(node_text), float(value_text)


def is_record_separator(text):
    """A bare ``##`` line separates records. ``## text`` is a snapshot label."""
    return text.startswith("#") and text.strip("#").strip() == ""


def read_dataassim(path, nnodes):
    """Return (dt, times, records, labels) from a nodal correction file.

    Each record is a length-``nnodes`` array. Listed nodes are filled in;
    omitted nodes and the fill value stay NaN. Node ids in the file are 1-based.
    A record does not have to include every mesh node. ``labels`` holds the
    ``## text`` description that precedes each record, or ``""`` when absent.
    """
    records = []
    labels = []
    label = ""
    values = np.full(nnodes, np.nan, dtype=np.float64)
    with open(path, "r") as handle:
        comment = handle.readline()
        if not comment.startswith("#"):
            raise ValueError(f"{path}: expected a comment header, got {comment!r}")
        dt = float(handle.readline())
        time0 = float(handle.readline())
        for line_number, line in enumerate(handle, start=4):
            text = line.strip()
            if not text:
                continue
            if is_record_separator(text):
                if np.isfinite(values).any():
                    records.append(values)
                    labels.append(label)
                values = np.full(nnodes, np.nan, dtype=np.float64)
                label = ""
                continue
            if text.startswith("#"):
                if not np.isfinite(values).any():
                    label = text.lstrip("#").strip()
                continue
            try:
                node, value = parse_node_value(text)
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: {exc}") from exc
            if value <= FILL_VALUE:
                continue
            if node < 1 or node > nnodes:
                raise ValueError(f"{path}:{line_number}: node id {node} is outside 1..{nnodes}")
            values[node - 1] = value
    if np.isfinite(values).any():
        records.append(values)
        labels.append(label)
    if not records:
        raise ValueError(f"{path}: no finite nodal values found")
    times = [time0 + i * dt for i in range(len(records))]
    return dt, times, records, labels


def mesh_triangulation(mesh):
    """Longitude, latitude, and 0-based triangle connectivity from an AdcircMesh."""
    x = np.asarray(mesh.x, dtype=np.float64)
    y = np.asarray(mesh.y, dtype=np.float64)
    triangulation = mesh.triangulation
    if isinstance(triangulation, mtri.Triangulation):
        triangles = np.asarray(triangulation.triangles, dtype=np.int64)
    else:
        triangles = np.asarray(triangulation, dtype=np.int64)
        if triangles.min() >= 1:
            triangles -= 1
    if triangles.shape[1] != 3:
        raise ValueError(f"expected triangular elements, got shape {triangles.shape}")
    if len(x) != len(y):
        raise ValueError("mesh x and y lengths differ")
    return x, y, triangles


def _collect_lines(geom, segments):
    if geom.is_empty:
        return
    if geom.geom_type == "LineString":
        coords = np.asarray(geom.coords)
        if len(coords) >= 2:
            segments.append(coords[:, :2])
    elif hasattr(geom, "geoms"):
        for part in geom.geoms:
            _collect_lines(part, segments)


def coastline_segments(xmin, xmax, ymin, ymax, resolution="10m"):
    """Download a Natural Earth coastline and clip it to the map window.

    Cartopy retrieves the shapefile from the Natural Earth data service
    (https://naturalearth.s3.amazonaws.com) and caches it locally.
    """
    print(f"Downloading Natural Earth {resolution} coastline")
    path = shapereader.natural_earth(
        resolution=resolution, category="physical", name="coastline"
    )
    window = box(xmin, ymin, xmax, ymax)
    segments = []
    for geom in shapereader.Reader(path).geometries():
        if not geom.intersects(window):
            continue
        _collect_lines(geom.intersection(window), segments)
    print(f"Shoreline segments in view: {len(segments)}")
    return segments


def plot_record(
    x,
    y,
    triangles,
    values,
    output,
    title,
    vmin,
    vmax,
    levels,
    cmap,
    dpi,
    xmin,
    xmax,
    ymin,
    ymax,
    shorelines,
    cbar_ticks_increment,
):
    """Write one filled-contour figure. Triangles touching a fill node are masked."""
    finite = np.isfinite(values)
    triang = mtri.Triangulation(x, y, triangles)
    triang.set_mask(np.any(~finite[triangles], axis=1))

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_facecolor("0.92")
    levels_plot = np.linspace(vmin, vmax, levels)
    contour = ax.tricontourf(
        triang,
        values,
        levels=levels_plot,
        cmap=cmap,
        extend="both",
    )
    if shorelines:
        ax.add_collection(
            LineCollection(shorelines, colors="black", linewidths=0.6, zorder=5)
        )
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    mid_lat = 0.5 * (ymin + ymax)
    ax.set_aspect(1.0 / np.cos(np.deg2rad(mid_lat)), adjustable="box")
    colorbar = fig.colorbar(contour, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("DWLC (m)")
    if cbar_ticks_increment <= 0:
        raise ValueError(f"colorbar tick interval must be positive, got {cbar_ticks_increment}")
    start = np.ceil(vmin / cbar_ticks_increment - 1e-9) * cbar_ticks_increment
    stop = np.floor(vmax / cbar_ticks_increment + 1e-9) * cbar_ticks_increment
    ticks = np.round(np.arange(start, stop + cbar_ticks_increment * 0.5, cbar_ticks_increment), 10)
    colorbar.set_ticks(ticks)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    in_view = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax) & finite
    finite_values = values[in_view]
    print(
        f"Wrote {output}  "
        f"finite nodes in view={finite_values.size}  "
        f"min={finite_values.min():.4f}  max={finite_values.max():.4f}"
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Plot nodal dataassim.dat fields as color contours on a fort.14 mesh."
    )
    parser.add_argument(
        "--analysis-dir",
        default=DEFAULT_ANALYSIS,
        help="Directory containing fort.14 and dataassim.dat",
    )
    parser.add_argument("--fort14", help="fort.14 mesh (default: <analysis-dir>/fort.14)")
    parser.add_argument(
        "--dataassim", help="dataassim.dat file (default: <analysis-dir>/dataassim.dat)"
    )
    parser.add_argument(
        "--output-dir",
        help="Directory for PNG files (default: the dataassim.dat directory)",
    )
    parser.add_argument(
        "--output",
        help=(
            "Output PNG path. '{coverage}' and '{record}' are replaced for each figure. "
            "If '{record}' is omitted, '_recordN' is inserted before the extension. "
            "Default: <output-dir>/dataassim_{coverage}_record{record}.png"
        ),
    )
    parser.add_argument(
        "--flip-sign",
        action="store_true",
        help="Multiply the nodal values by -1 before plotting",
    )
    parser.add_argument("--vmin", type=float, default=-0.5, help="Color scale minimum (default: -0.5)")
    parser.add_argument("--vmax", type=float, default=0.5, help="Color scale maximum (default: 0.5)")
    parser.add_argument("--levels", type=int, default=21, help="Number of contour levels")
    parser.add_argument(
        "--cbar-ticks-increment",
        type=float,
        default=0.1,
        help="Interval between ticks on the contour colorbar (default: 0.1)",
    )
    parser.add_argument("--cmap", default="RdBu_r", help="Matplotlib colormap")
    parser.add_argument("--dpi", type=int, default=150, help="Output resolution")
    parser.add_argument(
        "--shoreline-resolution",
        choices=("10m", "50m", "110m"),
        default="10m",
        help="Natural Earth coastline resolution downloaded by Cartopy (default: 10m)",
    )
    parser.add_argument(
        "--coverage",
        choices=sorted(COVERAGES),
        default=DEFAULT_COVERAGE,
        help=(
            "Named map window. 'nc' is the North Carolina coast; "
            "'gulf-to-ec' spans the Gulf of Mexico to the U.S. East Coast. "
            f"Default: {DEFAULT_COVERAGE}."
        ),
    )
    parser.add_argument(
        "--xmin", type=float, help="Western longitude; overrides the west edge of --coverage"
    )
    parser.add_argument(
        "--xmax", type=float, help="Eastern longitude; overrides the east edge of --coverage"
    )
    parser.add_argument(
        "--ymin", type=float, help="Southern latitude; overrides the south edge of --coverage"
    )
    parser.add_argument(
        "--ymax", type=float, help="Northern latitude; overrides the north edge of --coverage"
    )
    return parser


def resolve_window(args):
    """Return (xmin, xmax, ymin, ymax), with explicit bounds overriding --coverage."""
    xmin, xmax, ymin, ymax = COVERAGES[args.coverage]
    if args.xmin is not None:
        xmin = args.xmin
    if args.xmax is not None:
        xmax = args.xmax
    if args.ymin is not None:
        ymin = args.ymin
    if args.ymax is not None:
        ymax = args.ymax
    if xmin >= xmax or ymin >= ymax:
        raise ValueError(
            f"map window must have increasing limits, got "
            f"lon [{xmin}, {xmax}], lat [{ymin}, {ymax}]"
        )
    return xmin, xmax, ymin, ymax


def resolve_output_path(output, output_dir, coverage, record):
    """Build the PNG path for one record.

    ``output`` may contain ``{coverage}`` and ``{record}``. Without those
    placeholders, the record index is appended so each snapshot gets its own file.
    """
    if output:
        template = output
    else:
        template = os.path.join(output_dir, "dataassim_{coverage}_record{record}.png")
    if "{record}" not in template:
        root, ext = os.path.splitext(template)
        if not ext:
            ext = ".png"
        template = f"{root}_record{{record}}{ext}"
    path = template.format(coverage=coverage, record=record)
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    return path


def main(argv=None):
    args = get_parser().parse_args(argv)
    fort14 = args.fort14 or os.path.join(args.analysis_dir, "fort.14")
    dataassim = args.dataassim or os.path.join(args.analysis_dir, "dataassim.dat")
    output_dir = args.output_dir or os.path.dirname(os.path.abspath(dataassim))
    if args.output is None:
        os.makedirs(output_dir, exist_ok=True)

    print(f"Reading mesh {fort14}")
    mesh = AdcircMesh.open(fort14)
    x, y, triangles = mesh_triangulation(mesh)
    nnodes = len(x)
    print(f"Mesh nodes={nnodes}  elements={len(triangles)}")

    xmin, xmax, ymin, ymax = resolve_window(args)
    print(f"Reading {dataassim}")
    dt, times, records, labels = read_dataassim(dataassim, nnodes)
    if args.flip_sign:
        records = [-record for record in records]
        print("Flipped the sign of the nodal values")
    in_view = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
    viewed = [
        record[in_view & np.isfinite(record)] for record in records
    ]
    viewed = [values for values in viewed if values.size]
    if not viewed:
        raise ValueError(
            f"no finite values inside lon [{xmin}, {xmax}], lat [{ymin}, {ymax}]"
        )
    vmin = args.vmin
    vmax = args.vmax
    if vmin >= vmax:
        raise ValueError(f"color limits must increase, got vmin={vmin}, vmax={vmax}")
    print(
        f"Records={len(records)}  dt={dt:g} s  coverage={args.coverage}  "
        f"window=[{xmin}, {xmax}] x [{ymin}, {ymax}]  "
        f"color scale=[{vmin:.4f}, {vmax:.4f}] m"
    )
    shorelines = coastline_segments(xmin, xmax, ymin, ymax, args.shoreline_resolution)

    for index, (time_seconds, values, label) in enumerate(zip(times, records, labels), start=1):
        hours = time_seconds / 3600.0
        output = resolve_output_path(args.output, output_dir, args.coverage, index)
        title = f"DWLC, record {index}, t = {time_seconds:g} s ({hours:g} h)"
        if label:
            title = f"{title}\n{label}"
        plot_record(
            x,
            y,
            triangles,
            values,
            output,
            title,
            vmin,
            vmax,
            args.levels,
            args.cmap,
            args.dpi,
            xmin,
            xmax,
            ymin,
            ymax,
            shorelines,
            args.cbar_ticks_increment,
        )


if __name__ == "__main__":
    main()
