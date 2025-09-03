# %%
import os
import numpy as np
import xarray as xr
import matplotlib.tri as mtri
from scipy.spatial import KDTree
import pandas as pd
import pickle

# %%
def get_f64uv_at(f64file, stx, sty):
    cachefile = f64file.replace('.nc', '_{:9d}{:8d}_cache.pkl'.format(int(stx*1e6), int(sty*1e6)))
    if os.path.exists(cachefile):
        with open(cachefile, 'rb') as f:
            f64_time, f64_u, f64_v = pickle.load(f)
        return f64_time, f64_u, f64_v

    # %%
    def is_point_in_triangle_multi(px, py, x1, y1, x2, y2, x3, y3):
        denominator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3))
        a = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / denominator
        b = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / denominator
        c = 1 - a - b
        return np.all([np.all([0 <= a, a <= 1], axis=0),np.all([0 <= b, b <= 1], axis=0),np.all([0 <= c, c <= 1], axis=0)], axis=0)

    is_point_in_triangle_multi_vec = np.vectorize(is_point_in_triangle_multi, signature='(n),(n),(n),(n),(n),(n),(n),(n)->(n)')

    # %%
    with xr.open_dataset(f64file) as ds:
        adc_x = ds.x.values
        adc_y = ds.y.values
        adc_e = ds.element.values-1

    # %%
    isinside = is_point_in_triangle_multi(stx, sty, adc_x[adc_e[:,0]], adc_y[adc_e[:,0]], adc_x[adc_e[:,1]], adc_y[adc_e[:,1]], adc_x[adc_e[:,2]], adc_y[adc_e[:,2]])

    if not np.any(isinside):
        print(f"Point ({stx}, {sty}) is not inside any triangle. [{f64file}]")
        with xr.open_dataset(f64file) as ds:
            f64_time = pd.to_datetime(ds['time'].values.astype('datetime64[ms]')).to_pydatetime()
        f64_u = np.array([np.nan]*len(f64_time))
        f64_v = np.array([np.nan]*len(f64_time))
        with open(cachefile, 'wb') as f:
            pickle.dump([f64_time, f64_u, f64_v], f)
        return f64_time, f64_u, f64_v
    
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
    with xr.open_dataset(f64file) as ds:
        f64_time = pd.to_datetime(ds['time'].values.astype('datetime64[ms]')).to_pydatetime()
        nt = len(f64_time)
        f64_u = np.array([np.nan]*nt)
        f64_v = np.array([np.nan]*nt)
        for t in range(nt):
            elemv_u = ds['u-vel'][t][elemn].values
            elemv_v = ds['v-vel'][t][elemn].values
            f64_u[t] = interpolate_value(stx, sty, elemx, elemy, elemv_u)
            f64_v[t] = interpolate_value(stx, sty, elemx, elemy, elemv_v)

    # %%
    with open(cachefile, 'wb') as f:
        pickle.dump([f64_time, f64_u, f64_v], f)

    return f64_time, f64_u, f64_v 