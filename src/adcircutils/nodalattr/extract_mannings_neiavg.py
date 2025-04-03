# %%
import os
import rasterio
from rasterio.mask import mask
from rasterio.plot import show
import geowombat as gw
from pyproj import Transformer
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely import Point, Polygon
from adcircpy import AdcircMesh
import argparse

argparse = argparse.ArgumentParser(description='Extract the average Manning\'s n values of the neighboring elements of each node.')
argparse.add_argument('--landcoverfile', type=str, help='The path to the land cover raster file.')
argparse.add_argument('--meshfile', type=str, help='The path to the mesh file.')
argparse.add_argument('--class_to_mn_mapfile', type=str, help='The path to the class to mn map file.')
argparse.add_argument('--outputfile', type=str, help='The path to the output file.')
args = argparse.parse_args()

tiffile = args.landcoverfile
meshfile = args.meshfile
class_to_mn_file = args.class_to_mn_mapfile
outputfile = args.outputfile

# %%
mesh = AdcircMesh.open(meshfile)

# %%
class_to_mn = pd.read_csv(class_to_mn_file)
map_class = np.array(class_to_mn.classid)
map_mn = np.array(class_to_mn.mn)

# %%
class_to_mn

# %%
class_openwater = class_to_mn[class_to_mn.classification == "Open Water"].classid.values[0]
mn_openwater = class_to_mn[class_to_mn.classification == "Open Water"].mn.values[0]

# %%
nn = mesh.nodes.shape[0]
node_ids = []
node_polys = []
print('Creating nodal extent polygons...')
for i in range(nn):
    if i % 10000 == 0:
        print(i, '/', nn)
    nei = list(mesh.node_neighbors[i])
    nnei = len(nei)
    xs = list(mesh.x.iloc[nei])
    ys = list(mesh.y.iloc[nei])
    if len(xs) >= 3:
        xs.append(xs[0])
        ys.append(ys[0])
        node_ids.append(i)
        node_polys.append(Polygon(zip(xs, ys)))
print('done.')

# %%
dfpoly = gpd.GeoDataFrame({'node_id': node_ids, 'geometry': node_polys}, crs='EPSG:4326')

# %%
# Initialize a list to store the averaged mannings n values
mn_avgs = [mn_openwater]*nn

# Open the tif file to get its CRS
with rasterio.open(tiffile) as src:
    raster_crs = src.crs

# Convert the CRS of dfpoly to match the raster CRS
dfpoly_rc = dfpoly.to_crs(raster_crs)

# Open the tif file
with rasterio.open(tiffile) as src:
    # Loop through each polygon in dfpoly
    print('Extracting average mannings n values...')
    for idx, row in dfpoly_rc.iterrows():
        if idx % 1000 == 0:
            print(idx, '/', len(dfpoly_rc))
        node_id = row['node_id']
        polygon = [row['geometry'].__geo_interface__]
        
        # Mask the raster with the polygon
        try:
            out_image, out_transform = mask(src, polygon, crop=True)
        except:
            continue
        
        lcs = out_image.flatten()
        lcs = lcs[lcs != 0]

        if len(lcs) == 0:
            continue
        else:
            mn_sum = 0.0
            for lc in lcs:
                mn_sum += map_mn[map_class == lc][0]
            mn_avg = mn_sum / len(lcs)
            mn_avgs[node_id] = mn_avg
    print('done.')

# %%
node_id = np.arange(0, nn)
dfmn = pd.DataFrame({'node_id': node_id, 'mn_avg': mn_avgs})
dfmn.drop(dfmn[dfmn.mn_avg == mn_openwater].index, inplace=True)

# %%
with open(outputfile, "w") as f:
    f.write("{:.3f}\n".format(mn_openwater))
    f.write("{:d}\n".format(dfmn.shape[0]))
    for i in dfmn.index:
        f.write("{:d}    {:.3f}\n".format(dfmn.node_id[i] + 1, dfmn.mn_avg[i]))


