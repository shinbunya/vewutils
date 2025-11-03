# %%
import os
import numpy as np
import xarray as xr
import matplotlib.tri as mtri
from scipy.spatial import KDTree
import pandas as pd
import pickle

# %%
def get_solution_at(solution_file, stx, sty, var_names):
    """
    Extract solution variables at a specific point from a NetCDF solution file.
    
    Parameters
    ----------
    solution_file : str
        Path to the NetCDF solution file
    stx : float
        X-coordinate of the point
    sty : float
        Y-coordinate of the point
    var_names : str or list of str
        Variable name(s) to extract (e.g., 'zeta', ['u-vel', 'v-vel'], etc.)
    
    Returns
    -------
    dict
        Dictionary with 'time' key containing datetime array and variable names as keys
        containing their respective time series arrays
    """
    # Handle single variable name as string
    if isinstance(var_names, str):
        var_names = [var_names]
    
    # Create cache filename based on coordinates and variable names
    var_str = '_'.join(sorted(var_names))
    cachefile = solution_file.replace('.nc', '_{:9d}{:8d}_{}_cache.pkl'.format(
        int(stx*1e6), int(sty*1e6), var_str))
    
    if os.path.exists(cachefile):
        with open(cachefile, 'rb') as f:
            cached_data = pickle.load(f)
        return cached_data

    # %%
    def is_point_in_triangle_multi(px, py, x1, y1, x2, y2, x3, y3):
        denominator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3))
        a = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / denominator
        b = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / denominator
        c = 1 - a - b
        return np.all([np.all([0 <= a, a <= 1], axis=0),np.all([0 <= b, b <= 1], axis=0),np.all([0 <= c, c <= 1], axis=0)], axis=0)

    is_point_in_triangle_multi_vec = np.vectorize(is_point_in_triangle_multi, signature='(n),(n),(n),(n),(n),(n),(n),(n)->(n)')

    # %%
    with xr.open_dataset(solution_file) as ds:
        adc_x = ds.x.values
        adc_y = ds.y.values
        adc_e = ds.element.values-1

    # %%
    isinside = is_point_in_triangle_multi(stx, sty, adc_x[adc_e[:,0]], adc_y[adc_e[:,0]], 
                                        adc_x[adc_e[:,1]], adc_y[adc_e[:,1]], 
                                        adc_x[adc_e[:,2]], adc_y[adc_e[:,2]])

    if not np.any(isinside):
        print(f"Point ({stx}, {sty}) is not inside any triangle. [{solution_file}]")
        with xr.open_dataset(solution_file) as ds:
            f63_time = pd.to_datetime(ds['time'].values.astype('datetime64[ms]')).to_pydatetime()
        
        # Create result dictionary with NaN values for all variables
        result = {'time': f63_time}
        for var_name in var_names:
            result[var_name] = np.array([np.nan]*len(f63_time))
        
        with open(cachefile, 'wb') as f:
            pickle.dump(result, f)
        return result
    
    # %%
    elemn = adc_e[isinside][0]
    elemx = adc_x[adc_e[isinside]][0]
    elemy = adc_y[adc_e[isinside]][0]

    # %%
    def interpolate_value(px, py, ex, ey, ev):
        x1 = ex[0]
        y1 = ey[0]
        x2 = ex[1]
        y2 = ey[1]
        x3 = ex[2]
        y3 = ey[2]
        v1 = ev[0]
        v2 = ev[1]
        v3 = ev[2]
        denominator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3))
        a = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / denominator
        b = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / denominator
        c = 1 - a - b
        return a * v1 + b * v2 + c * v3

    # %%
    with xr.open_dataset(solution_file) as ds:
        f63_time = pd.to_datetime(ds['time'].values.astype('datetime64[ms]')).to_pydatetime()
        nt = len(f63_time)
        
        # Initialize result dictionary
        result = {'time': f63_time}
        
        # Process each variable
        for var_name in var_names:
            if var_name not in ds.variables:
                print(f"Warning: Variable '{var_name}' not found in {solution_file}")
                result[var_name] = np.array([np.nan]*nt)
                continue
                
            var_data = np.array([np.nan]*nt)
            
            for t in range(nt):
                try:
                    # Get variable data for the element at time t
                    if var_name in ['u-vel', 'v-vel']:
                        # Handle velocity components which might be stored differently
                        elemv = ds[var_name][t][elemn].values
                    else:
                        # Handle other variables like zeta
                        elemv = ds[var_name][t][elemn].values
                    
                    var_data[t] = interpolate_value(stx, sty, elemx, elemy, elemv)
                except Exception as e:
                    print(f"Warning: Error processing {var_name} at time {t}: {e}")
                    var_data[t] = np.nan
            
            result[var_name] = var_data

    # %%
    with open(cachefile, 'wb') as f:
        pickle.dump(result, f)

    return result
