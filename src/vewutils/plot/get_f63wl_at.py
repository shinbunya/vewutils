# %%
import os
import pickle

import numpy as np
import xarray as xr
import pandas as pd
from netCDF4 import Dataset


# %%
def get_f63wl_at(f63file, stx, sty):
    """
    Extract water level time series at a given (stx, sty) from an ADCIRC f63
    NetCDF file by locating the containing triangle and doing barycentric
    interpolation of zeta at its three nodes for all times.
    Results are cached on disk for subsequent calls.
    """
    cachefile = f63file.replace(
        ".nc",
        "_{:9d}{:8d}_cache.pkl".format(int(stx * 1e6), int(sty * 1e6)),
    )

    # ------------------------------------------------------------------
    # 1. Cache short‑circuit
    # ------------------------------------------------------------------
    if os.path.exists(cachefile):
        with open(cachefile, "rb") as f:
            f63_time, f63_wl = pickle.load(f)
        return f63_time, f63_wl

    # ------------------------------------------------------------------
    # 2. Helper: point-in-triangle test (vectorized over many elements)
    # ------------------------------------------------------------------
    def is_point_in_triangle_multi(px, py, x1, y1, x2, y2, x3, y3):
        denominator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3))
        a = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / denominator
        b = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / denominator
        c = 1.0 - a - b
        return np.all(
            [
                np.all([0.0 <= a, a <= 1.0], axis=0),
                np.all([0.0 <= b, b <= 1.0], axis=0),
                np.all([0.0 <= c, c <= 1.0], axis=0),
            ],
            axis=0,
        )

    # ------------------------------------------------------------------
    # 3. Use xarray for mesh + time (cheap compared to zeta)
    # ------------------------------------------------------------------
    with xr.open_dataset(f63file) as ds:
        # Mesh
        adc_x = ds["x"].values
        adc_y = ds["y"].values
        adc_e = ds["element"].values - 1  # 0‑based

        # Find containing element
        isinside = is_point_in_triangle_multi(
            stx,
            sty,
            adc_x[adc_e[:, 0]],
            adc_y[adc_e[:, 0]],
            adc_x[adc_e[:, 1]],
            adc_y[adc_e[:, 1]],
            adc_x[adc_e[:, 2]],
            adc_y[adc_e[:, 2]],
        )

        # Time axis
        f63_time = pd.to_datetime(
            ds["time"].values.astype("datetime64[ms]")
        ).to_pydatetime()
        nt = len(f63_time)

    # ------------------------------------------------------------------
    # 4. Point not inside any triangle
    # ------------------------------------------------------------------
    if not np.any(isinside):
        print(f"Point ({stx}, {sty}) is not inside any triangle. [{f63file}]")
        f63_wl = np.full(nt, np.nan)
        with open(cachefile, "wb") as f:
            pickle.dump([f63_time, f63_wl], f)
        return f63_time, f63_wl

    # ------------------------------------------------------------------
    # 5. Element geometry + barycentric weights (only 3 nodes)
    # ------------------------------------------------------------------
    elemn = adc_e[isinside][0]
    elemx = adc_x[elemn]
    elemy = adc_y[elemn]

    x1, y1 = elemx[0], elemy[0]
    x2, y2 = elemx[1], elemy[1]
    x3, y3 = elemx[2], elemy[2]

    denominator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3))
    a = ((y2 - y3) * (stx - x3) + (x3 - x2) * (sty - y3)) / denominator
    b = ((y3 - y1) * (stx - x3) + (x1 - x3) * (sty - y3)) / denominator
    c = 1.0 - a - b

    # ------------------------------------------------------------------
    # 6. Heavy part: read zeta with netCDF4 directly (time, node)
    # ------------------------------------------------------------------
    with Dataset(f63file, mode="r") as nc:
        zeta_var = nc.variables["zeta"]  # dims typically ('time', 'node')
        
        # Detect chunking to choose optimal read strategy
        chunking = zeta_var.chunking()
        if chunking and len(chunking) == 2:
            time_chunk, node_chunk = chunking
            # If chunked by node (node_chunk is small, e.g., 1), fancy indexing is faster
            # If chunked by time (time_chunk is small, e.g., 1), individual reads are faster
            is_node_chunked = node_chunk <= 10 and time_chunk > 100
        else:
            # Default: assume bad chunking (chunked by time)
            is_node_chunked = False
        
        if is_node_chunked:
            # Optimized chunking: read all 3 nodes at once (faster for node-chunked files)
            elemv_all = zeta_var[:, elemn]
        else:
            # Bad chunking: read each node separately (faster for time-chunked files)
            # When file is chunked [1, 506054], reading nodes individually is much faster
            # than fancy indexing because it benefits from internal caching
            node1_data = zeta_var[:, elemn[0]]
            node2_data = zeta_var[:, elemn[1]]
            node3_data = zeta_var[:, elemn[2]]
            elemv_all = np.column_stack([node1_data, node2_data, node3_data])

    # ------------------------------------------------------------------
    # 7. Vectorized interpolation over time
    # ------------------------------------------------------------------
    f63_wl = a * elemv_all[:, 0] + b * elemv_all[:, 1] + c * elemv_all[:, 2]

    # ------------------------------------------------------------------
    # 8. Cache results
    # ------------------------------------------------------------------
    with open(cachefile, "wb") as f:
        pickle.dump([f63_time, f63_wl], f)

    return f63_time, f63_wl
