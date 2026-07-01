import os
import pickle

import numpy as np
import pandas as pd
import xarray as xr


def generate_transect_points(start_x, start_y, end_x, end_y, spacing):
    """
    Generate equispaced points along a straight-line transect.

    Returns
    -------
    distances : ndarray
        Distance along transect from the start point (m)
    x : ndarray
        X-coordinates of extraction points
    y : ndarray
        Y-coordinates of extraction points
    """
    dx = end_x - start_x
    dy = end_y - start_y
    length = np.hypot(dx, dy)

    if length == 0:
        return np.array([0.0]), np.array([start_x]), np.array([start_y])

    distances = np.arange(0.0, length + spacing * 0.5, spacing)
    if distances[-1] < length:
        distances = np.append(distances, length)

    t = distances / length
    x = start_x + t * dx
    y = start_y + t * dy
    return distances, x, y


def _is_point_in_triangle_multi(px, py, x1, y1, x2, y2, x3, y3):
    denominator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3))
    a = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / denominator
    b = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / denominator
    c = 1 - a - b
    return np.all(
        [
            np.all([0 <= a, a <= 1], axis=0),
            np.all([0 <= b, b <= 1], axis=0),
            np.all([0 <= c, c <= 1], axis=0),
        ],
        axis=0,
    )


def _interpolate_value(px, py, ex, ey, ev):
    x1, y1 = ex[0], ey[0]
    x2, y2 = ex[1], ey[1]
    x3, y3 = ex[2], ey[2]
    v1, v2, v3 = ev[0], ev[1], ev[2]
    denominator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3))
    a = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / denominator
    b = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / denominator
    c = 1 - a - b
    return a * v1 + b * v2 + c * v3


def _find_element(stx, sty, adc_x, adc_y, adc_e):
    isinside = _is_point_in_triangle_multi(
        stx,
        sty,
        adc_x[adc_e[:, 0]],
        adc_y[adc_e[:, 0]],
        adc_x[adc_e[:, 1]],
        adc_y[adc_e[:, 1]],
        adc_x[adc_e[:, 2]],
        adc_y[adc_e[:, 2]],
    )
    if not np.any(isinside):
        return None
    element_idx = int(np.where(isinside)[0][0])
    elemn = adc_e[element_idx]
    elemx = adc_x[elemn]
    elemy = adc_y[elemn]
    return element_idx, elemn, elemx, elemy


def resolve_timestep(time_values, timestep=None, datetime_str=None):
    """
    Resolve a time index from either a timestep index or a datetime string.
    """
    nt = len(time_values)
    if nt == 0:
        raise ValueError("Solution file contains no time steps")

    if timestep is not None:
        if timestep == -1:
            return nt - 1
        if timestep < 0 or timestep >= nt:
            raise ValueError(
                f"Time step {timestep} not found. Available time steps: 0 to {nt - 1}"
            )
        return timestep

    target = pd.to_datetime(datetime_str)
    times = pd.to_datetime(time_values.astype("datetime64[ms]"))
    return int(np.argmin(np.abs(times - target)))


def get_solution_along_transect(
    solution_file,
    start_x,
    start_y,
    end_x,
    end_y,
    spacing,
    var_names,
    timestep=None,
    datetime_str=None,
    use_cache=True,
):
    """
    Extract solution variables along a transect at a single time step.

    Parameters
    ----------
    solution_file : str
        Path to the NetCDF solution file
    start_x, start_y, end_x, end_y : float
        Transect endpoints
    spacing : float
        Spacing between extraction points in meters
    var_names : str or list of str
        Variable name(s) to extract
    timestep : int, optional
        0-based time step index (-1 for last time step)
    datetime_str : str, optional
        Target time in YYYY-MM-DDThh:mm:ss format
    use_cache : bool, optional
        If True, read from and write to on-disk cache files, by default True

    Returns
    -------
    dict
        Dictionary with distance, x, y, time, timestep, and requested variables
    """
    if isinstance(var_names, str):
        var_names = [var_names]

    if (timestep is None) == (datetime_str is None):
        raise ValueError("Exactly one of timestep or datetime_str must be provided")

    distances, px, py = generate_transect_points(
        start_x, start_y, end_x, end_y, spacing
    )

    cachefile = None
    if use_cache:
        var_str = "_".join(sorted(var_names))
        time_key = f"t{timestep}" if timestep is not None else datetime_str.replace(":", "")
        cachefile = solution_file.replace(
            ".nc",
            "_{:9d}{:8d}_{:9d}{:8d}_{:g}_{}_{}_transect_cache.pkl".format(
                int(start_x * 1e6),
                int(start_y * 1e6),
                int(end_x * 1e6),
                int(end_y * 1e6),
                spacing,
                var_str,
                time_key,
            ),
        )

        if os.path.exists(cachefile):
            with open(cachefile, "rb") as f:
                return pickle.load(f)

    with xr.open_dataset(solution_file) as ds:
        adc_x = ds.x.values
        adc_y = ds.y.values
        adc_e = ds.element.values - 1
        time_values = ds["time"].values
        timestep_to_use = resolve_timestep(time_values, timestep, datetime_str)
        selected_time = pd.to_datetime(
            time_values[timestep_to_use].astype("datetime64[ms]")
        ).to_pydatetime()

        result = {
            "distance": distances,
            "x": px,
            "y": py,
            "time": selected_time,
            "timestep": timestep_to_use,
        }

        for var_name in var_names:
            if var_name not in ds.variables:
                print(f"Warning: Variable '{var_name}' not found in {solution_file}")
                result[var_name] = np.full(len(distances), np.nan)
                continue

            var_shape = ds[var_name].shape
            on_elements = len(var_shape) == 3 and var_shape[2] == 1

            if len(var_shape) == 1:
                field = ds[var_name].values
            elif on_elements:
                field = ds[var_name][timestep_to_use, :, 0].values
            else:
                field = ds[var_name][timestep_to_use, :].values

            values = np.full(len(distances), np.nan)
            for i, (stx, sty) in enumerate(zip(px, py)):
                element = _find_element(stx, sty, adc_x, adc_y, adc_e)
                if element is None:
                    print(
                        f"Point ({stx}, {sty}) is not inside any triangle. [{solution_file}]"
                    )
                    continue

                element_idx, elemn, elemx, elemy = element
                if on_elements:
                    values[i] = field[element_idx]
                else:
                    elemv = field[elemn]
                    values[i] = _interpolate_value(stx, sty, elemx, elemy, elemv)

            result[var_name] = values

    if use_cache:
        with open(cachefile, "wb") as f:
            pickle.dump(result, f)

    return result
