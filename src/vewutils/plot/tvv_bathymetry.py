"""
Time-varying crest VEW (TVV, fort.142 ``VaryType=4``) wall-top depths.

ADCIRC moves the wall-top node of each ``VaryType=4`` VEW pair with its crest
table, ``DP_top(t) = delta - z_c(t)``, but output files such as fort.63.nc hold
only the initial ``depth``. These helpers rebuild the wall-top depth at a given
model time so plots can show the gate where it is.
"""
import os
import re

import numpy as np

_DEFAULT_SEARCH_RADIUS = 1.0e-6
_KEYVAL = re.compile(r"(\w+)\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s]+)")


def read_tvw_file(fort142_path):
    """
    Read the ``VaryType=4`` lines of a time-varying weir file (fort.142).

    Returns
    -------
    list of dict
        One dict per line with ``x``, ``y``, ``search_radius`` (None if not
        given) and ``schedule_file`` (resolved relative to the fort.142 directory).
        Lines of other VaryTypes are skipped.
    """
    base_dir = os.path.dirname(os.path.abspath(fort142_path))
    with open(fort142_path) as f:
        lines = [line.strip() for line in f if line.strip()]
    nlines = int(lines[0].split()[0])

    entries = []
    for line in lines[1:1 + nlines]:
        fields = {k.upper(): v.strip("'\"") for k, v in _KEYVAL.findall(line)}
        if int(float(fields.get("VARYTYPE", 0))) != 4:
            continue
        radius = fields.get("SEARCHRADIUS")
        entries.append(
            {
                "x": float(fields["X1"]),
                "y": float(fields["Y1"]),
                "search_radius": float(radius) if radius is not None else None,
                "schedule_file": os.path.join(base_dir, fields["SCHEDULEFILE"]),
            }
        )
    return entries


def read_crest_table(path):
    """
    Read a TVV crest table: a record count, then ``time  z_c`` per line.

    Returns
    -------
    times, crests : ndarray
    """
    with open(path) as f:
        lines = [line.split() for line in f if line.strip()]
    nrec = int(lines[0][0])
    data = np.array([[float(v) for v in rec[:2]] for rec in lines[1:1 + nrec]])
    times, crests = data[:, 0], data[:, 1]
    if np.any(np.diff(times) <= 0):
        raise ValueError(f"Crest table times must increase strictly: {path}")
    return times, crests


def crest_at(times, crests, t):
    """Crest elevation at time t: linear between records, held outside them."""
    return float(np.interp(t, times, crests))


def _locate_pair(entry, adc_x, adc_y, depth, default_radius):
    radius = entry["search_radius"]
    if radius is None:
        radius = default_radius
    dist = np.hypot(adc_x - entry["x"], adc_y - entry["y"])
    nodes = np.where(dist <= radius)[0]
    if len(nodes) != 2:
        raise ValueError(
            f"Expected the 2 nodes of one VEW pair within {radius:g} of "
            f"({entry['x']}, {entry['y']}), found {len(nodes)}: "
            f"{(nodes + 1).tolist()}. Adjust SearchRadius (or --tvv-search-radius)."
        )
    # The wall-top node is the shallower node of the pair.
    top, bed = nodes[np.argsort(depth[nodes])]
    return int(top), int(bed)


def tvv_depth(fort142_path, adc_x, adc_y, depth, model_time, default_radius=None):
    """
    Return a copy of the nodal depth with TVV wall-top nodes moved to their
    crest at model_time (seconds, the solution file's time axis).

    The wall-top depth follows ``DP_top(t) = DP_top(0) - (z_c(t) - z_c(0))``,
    which equals ADCIRC's ``delta - z_c(t)`` given its cold-start check that the
    crest table starts at the fort.14 barrier height. Crests below
    ``bed + delta`` are clamped as in ADCIRC.
    """
    if default_radius is None:
        default_radius = _DEFAULT_SEARCH_RADIUS
    new_depth = np.array(depth, dtype=float, copy=True)
    tables = {}
    for entry in read_tvw_file(fort142_path):
        sched = entry["schedule_file"]
        if sched not in tables:
            tables[sched] = read_crest_table(sched)
        times, crests = tables[sched]
        top, bed = _locate_pair(entry, adc_x, adc_y, depth, default_radius)
        z0 = crest_at(times, crests, 0.0)
        delta = z0 + depth[top]
        zc = max(crest_at(times, crests, model_time), -depth[bed] + delta)
        new_depth[top] = delta - zc
    return new_depth
