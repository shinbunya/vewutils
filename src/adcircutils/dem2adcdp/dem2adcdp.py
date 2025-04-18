# -*- coding: utf-8 -*-
"""
Assignment of depth to mesh nodes based on DEM data

@author: shinbunya
"""

import pandas as pd
from geopandas import GeoDataFrame
import matplotlib.tri as mtri
import numpy as np
from shapely import Polygon, Point
from shapely.ops import unary_union
import rasterio
from rasterstats import zonal_stats
import pickle
import re
from ipywidgets import IntProgress, Label, VBox
from IPython.display import display, clear_output
import time
from multiprocessing import Pool, Process, SimpleQueue
from functools import partial
import math
import geopandas as gpd

from adcircutils.dem2adcdp.adcmesh import Mesh
from adcircutils.dem2adcdp.adcmesh import F13

class DEM2DP:
    zonal_stats_node_cache = {}
    polygon_chunk_size = 100
    n2e = None
    mesh = None

    global wProgress, wLabel
    wProgress = None
    wLabel = None
    max_count = None
    time_start = None
    header = None
    
    def __init__(self):
        self.target_loaded = False

    def read_mesh(self, meshfile):
        self.mesh = Mesh()
        self.mesh.read(meshfile)

    def write_mesh(self, meshfile):
        self.mesh.write(meshfile)

    def read_f13(self, f13file):
        self.f13 = F13()
        self.f13.read(f13file)

    def gen_hash(lon, lat, area):
        return '{:+012.7f}{:+011.8f}{:+.0f}'.format(lon, lat, area)
    
    def init_progressbar(max_count, header):
        DEM2DP.max_count = max_count
        DEM2DP.time_start = time.time()
        DEM2DP.wProgress = IntProgress(min=0, max=max_count)
        DEM2DP.wLabel = Label()
        DEM2DP.header = header
        wVBox = VBox([DEM2DP.wProgress, DEM2DP.wLabel])

    def show_progress(i):
        percentage_done = max(1e-8, float(i)/float(DEM2DP.max_count)*100.0)
        elapsed_time = time.time() - DEM2DP.time_start
        elapsed_unit = 'sec'
        estimated_total_time = elapsed_time / (percentage_done/100.0)
        estimated_total_unit = 'sec'
        remaining_time = elapsed_time * (1/(percentage_done/100.0) - 1.0)
        remaining_unit = 'sec'
        if elapsed_time > 600:
            elapsed_time /= 60
            elapsed_unit = 'min'
        if estimated_total_time > 600:
            estimated_total_time /= 60
            estimated_total_unit = 'min'
        if remaining_time > 600:
            remaining_time /= 60
            remaining_unit = 'min'
        DEM2DP.wProgress.value = i
        msg = '{}: {:.2f}% done. elapsed time = {:.0f} {}. estimated remaining time = {:.0f} {}. estimated total time = {:.0f} {}.'.format(\
            DEM2DP.header, percentage_done, elapsed_time, elapsed_unit, remaining_time, remaining_unit, estimated_total_time, estimated_total_unit)
        DEM2DP.wLabel.value = msg 
        print(msg, flush=True)

    def update_progress():
        while True:
            DEM2DP.show_progress(DEM2DP.wProgress.value)
            time.sleep(1)

    def init_worker(shared_queue):
        global queue
        queue = shared_queue

    def DEM2DPElem(self, tiffile):
        with rasterio.open(tiffile) as src:
            src_crs = src.crs

        numEle = self.mesh.numEle
        elevs = np.zeros_like(self.mesh.coord['Longitude'])
        elevs[:] = np.nan
        areas = np.zeros_like(self.mesh.coord['Longitude'])
        areas[:] = np.nan
        polys = []
        for jelem in range(numEle):
            tri = self.mesh.tri[jelem,:]
            poly = Polygon(zip(self.mesh.coord.Longitude[tri],self.mesh.coord.Latitude[tri]))
            polys.append(poly)

        gdf = GeoDataFrame(index=list(range(len(polys))), crs='epsg:4326', geometry=polys)
        self.zonal_stats_elem = zonal_stats(gdf.to_crs(src_crs), tiffile)

        return

    def generate_nodal_poly(self, meshxy, bounds, nodes):
        global queue

        # determining chunks
        nnodes = len(nodes)
        if nnodes == 0: return []

        chunk_size = max(1, min(10000, nnodes))
        nodes_per_chunk = min(chunk_size, nnodes)
        chunks = [nodes[i*nodes_per_chunk:(i+1)*nodes_per_chunk] for i in range(math.floor(nnodes/nodes_per_chunk))]
        if max(chunks[-1]) != max(nodes):
            i = math.floor(nnodes/nodes_per_chunk)
            chunks.append(nodes[i*nodes_per_chunk:])

        polys = []
        for nodes_chunk in chunks:
            for index, i in enumerate(nodes_chunk):
                # check target
                if self.target_loaded and not self.target[i]:
                    polys.append(None)
                    continue

                # check bounds
                x = meshxy[i].x
                y = meshxy[i].y
                if x < bounds.left or x > bounds.right or y < bounds.bottom or y > bounds.top:
                    polys.append(None)
                    continue

                # generate poly
                pts = []
                elems = self.mesh.n2e[i]
                j = elems[0]
                jtri = list(self.mesh.tri[j,:])
                idx = jtri.index(i)
                i0 = jtri[(idx + 1) % 3]
                i1 = jtri[(idx + 2) % 3]
                pts = [i0, i1]
                elems = elems[1:]
                while True:
                    found = False
                    for jjj, jj in enumerate(elems):
                        jtri = list(self.mesh.tri[jj,:])
                        if i1 in jtri:
                            idx = jtri.index(i1)
                            i2 = jtri[(idx + 1) % 3]
                            if i2 != i and not i2 in pts:
                                found = True
                                break
                    if not found:
                        if i in pts:
                            break
                        else:
                            i1 = i
                            pts.append(i1)
                    else:
                        i1 = i2
                        pts.append(i1)
                        if i1 == i0:
                            break

                polys.append(Polygon(zip(self.mesh.coord.loc[pts, 'Longitude'],self.mesh.coord.loc[pts, 'Latitude'])))

            queue.put(len(nodes_chunk))

        return polys


    def obtain_zonal_stats_for_node(self, tiffile, gdfpolys, nodes):
        global queue

        nnodes = len(nodes)
        ret_zonal_stats_node = [None]*nnodes
        if nnodes == 0: return ret_zonal_stats_node

        # determining chunks
        chunk_size = max(1, min(1000, nnodes))
        nodes_per_chunk = min(chunk_size, nnodes)
        chunks = [nodes[i*nodes_per_chunk:(i+1)*nodes_per_chunk] for i in range(math.floor(nnodes/nodes_per_chunk))]
        if max(chunks[-1]) != max(nodes):
            i = math.floor(nnodes/nodes_per_chunk)
            chunks.append(nodes[i*nodes_per_chunk:])

        # generationg zonal stats
        cnt = 0
        for ichunk, nodes_chunk in enumerate(chunks):
            cnt_cached = 0
            target_nodes = []
            target_hashes = []
            target_indexes = []
            for index, i in enumerate(nodes_chunk):
                # check if poly is valid
                if gdfpolys.loc[i, 'geometry'] == None:
                    cnt += 1
                    continue

                # create hash
                lon = self.mesh.coord.loc[i,'Longitude']
                lat = self.mesh.coord.loc[i,'Latitude']
                area = gdfpolys.loc[i, 'geometry'].area
                hash = DEM2DP.gen_hash(lon, lat, area)

                # checking cache
                iscached = False
                if not self.zonal_stats_node_cache:
                    pass
                else:
                    if hash in self.zonal_stats_node_cache.keys():
                        iscached = True
                if not iscached:
                    target_nodes.append(i)
                    target_hashes.append(hash)
                    target_indexes.append(cnt)
                else:
                    ret_zonal_stats_node[cnt] = self.zonal_stats_node_cache[hash].to_frame().T
                    cnt_cached += 1
                cnt += 1

            # print("cached nodes: {:d}/{:d} at {:d}/{:d} ({:d},{:d})".format(cnt_cached, len(nodes_chunk), ichunk+1, len(chunks), len(nodes), nodes[0]), flush=True)

            if len(target_nodes) > 0:
                target_polys = gdfpolys.loc[target_nodes, 'geometry']
                zonal_stats_node = pd.DataFrame(zonal_stats(target_polys, tiffile, stats=['mean', 'max', 'min', 'count']))

                for ii in range(len(target_nodes)):
                    i = target_nodes[ii]
                    index = target_indexes[ii]
                    lon = self.mesh.coord['Longitude'][i]
                    lat = self.mesh.coord['Latitude'][i]
                    area = gdfpolys['geometry'].area[i]
                    hash = DEM2DP.gen_hash(lon, lat, area)

                    zonal_stats_node.loc[ii, 'lon'] = lon 
                    zonal_stats_node.loc[ii, 'lat'] = lat
                    zonal_stats_node.loc[ii, 'area'] = area
                    zonal_stats_node.loc[ii, 'hash'] = hash

                    ret_zonal_stats_node[index] = zonal_stats_node.loc[ii].to_frame().T

            # DEM2DP.show_progress(nodes_chunk[-1])
            queue.put(len(nodes_chunk))

        return ret_zonal_stats_node

    def add_all(self):
        print('adding all nodes', flush=True)
        self.target = [True] * self.mesh.numNod
        self.target_loaded = True
        print('- done', flush=True)

    def remove_all(self):
        print('removing all nodes', flush=True)
        self.target = [False] * self.mesh.numNod
        self.target_loaded = True
        print('- done', flush=True)

    def add_by_polygons(self, polygon_filename, refresh=False):
        print('adding target by polygon', flush=True)

        print('- reading polygon file: {}'.format(polygon_filename))
        # self.mesh.read_boundary_polygon(polygon_filename, msgout=False)
        # merged_poly = self.mesh.boundary_polygon.buffer(1e-6)

        print('- setting target flags', flush=True)
        gdf = gpd.read_file(polygon_filename)
        for poly in gdf['geometry']:
            merged_poly = unary_union(poly).buffer(1e-6)
            target = merged_poly.contains([Point(lon, lat) for lon, lat in zip(self.mesh.coord['Longitude'], self.mesh.coord['Latitude'])])

            if refresh or not hasattr(self, 'target'):
                self.target = target
            else:
                self.target = [a or b for a, b in zip(self.target, target)]

        self.target_loaded = True

        ntarget = np.count_nonzero(self.target)        
        print('- number of true values in mesh target: {} / {}'.format(ntarget, len(self.target)), flush=True)
        print('- done', flush=True)

    def remove_by_polygons(self, polygon_filename, refresh=False):
        print('removing target by polygon', flush=True)

        print('- reading polygon file: {}'.format(polygon_filename), flush=True)
        self.mesh.read_boundary_polygon(polygon_filename, msgout=False)
        merged_poly = self.mesh.boundary_polygon.buffer(-1e-6)

        print('- setting target flags', flush=True)
        target = merged_poly.contains([Point(lon, lat) for lon, lat in zip(self.mesh.coord['Longitude'], self.mesh.coord['Latitude'])])
        target = [not a for a in target]

        if refresh or not hasattr(self, 'target'):
            self.target = target
        else:
            self.target = [a and b for a, b in zip(self.target, target)]

        self.target_loaded = True

        ntarget = np.count_nonzero(self.target)
        print('- number of true values in mesh target: {} / {}'.format(ntarget, len(self.target)), flush=True)
        print('- done', flush=True)
        
    def remove_by_polygons_outside(self, polygon_filename, refresh=False):
        print('removing target outside of polygon', flush=True)

        print('- reading polygon file: {}'.format(polygon_filename), flush=True)
        self.mesh.read_boundary_polygon(polygon_filename, msgout=False)
        merged_poly = self.mesh.boundary_polygon.buffer(-1e-6)

        print('- setting target flags', flush=True)
        target = merged_poly.contains([Point(lon, lat) for lon, lat in zip(self.mesh.coord['Longitude'], self.mesh.coord['Latitude'])])

        if refresh or not hasattr(self, 'target'):
            self.target = target
        else:
            self.target = [a and b for a, b in zip(self.target, target)]

        self.target_loaded = True

        ntarget = np.count_nonzero(self.target)
        print('- number of true values in mesh target: {} / {}'.format(ntarget, len(self.target)), flush=True)
        print('- done', flush=True)
        
    def add_channelnodes(self):
        print('setting target flags for channel nodes', flush=True)
        if not hasattr(self, 'target'):
            target = [False] * self.mesh.numNod
        else:
            target = self.target

        for i in self.mesh.channelnodes:
            target[i] = True

        self.target = target        

        self.target_loaded = True

        ntarget = np.count_nonzero(self.target)        
        print('- number of true values in mesh target: {} / {}'.format(ntarget, len(self.target)), flush=True)
        print('- done', flush=True)

    def remove_channelnodes(self):
        print('removing target flags for channel nodes', flush=True)
        if not hasattr(self, 'target'):
            target = [True] * self.mesh.numNod
        else:
            target = self.target

        for i in self.mesh.channelnodes:
            target[i] = False

        self.target = target        

        self.target_loaded = True

        ntarget = np.count_nonzero(self.target)        
        print('- number of true values in mesh target: {} / {}'.format(ntarget, len(self.target)), flush=True)
        print('- done', flush=True)

    def add_banknodes(self):
        print('setting target flags for bank nodes', flush=True)
        if not hasattr(self, 'target'):
            target = [False] * self.mesh.numNod
        else:
            target = self.target

        for i in self.mesh.banknodes:
            target[i] = True

        self.target = target        

        self.target_loaded = True

        ntarget = np.count_nonzero(self.target)        
        print('- number of true values in mesh target: {} / {}'.format(ntarget, len(self.target)), flush=True)
        print('- done', flush=True)

    def remove_banknodes(self):
        print('removing target flags for bank nodes', flush=True)
        if not hasattr(self, 'target'):
            target = [True] * self.mesh.numNod
        else:
            target = self.target

        for i in self.mesh.banknodes:
            target[i] = False

        self.target = target        

        self.target_loaded = True

        ntarget = np.count_nonzero(self.target)        
        print('- number of true values in mesh target: {} / {}'.format(ntarget, len(self.target)), flush=True)
        print('- done', flush=True)

    def add_boundarynodes(self):
        print('setting target flags for boundary nodes', flush=True)
        if not hasattr(self, 'target'):
            target = [False] * self.mesh.numNod
        else:
            target = self.target

        for i in self.mesh.boundarynodes:
            target[i] = False

        self.target = target        

        self.target_loaded = True

        ntarget = np.count_nonzero(self.target)        
        print('- number of true values in mesh target: {} / {}'.format(ntarget, len(self.target)), flush=True)
        print('- done', flush=True)

    def extract_dem_vals_on_mesh_nodes(self, tiffile, ncores=8, chunk_size_poly=100000, chunk_size_zonalstats=1000):

        with rasterio.open(tiffile) as src:
            src_crs = src.crs
            bounds = src.bounds
        
        lons = self.mesh.coord['Longitude']
        lats = self.mesh.coord['Latitude']
        points = [Point(lon, lat) for lon, lat in zip(lons, lats)]
        meshxy = GeoDataFrame(index=list(range(len(lons))), crs='epsg:4326', geometry=points).to_crs(src_crs).geometry

        numEle = self.mesh.numEle
        numNod = self.mesh.numNod
        elevs = np.zeros_like(self.mesh.coord['Longitude'])
        elevs[:] = np.nan
        
        # determine target nodes
        if self.target_loaded:
            target_nodes = [i for i in range(numNod) if self.target[i]]
        else:
            target_nodes = list(range(numNod))

        # generate polygons
        print('generating nodal polygons', flush=True)
        max_count = len(target_nodes)
        chunk_size = chunk_size_poly
        nodes_per_task = min(chunk_size, math.floor(len(target_nodes)/ncores))
        chunks = [[target_nodes[j] for j in range(i*nodes_per_task,(i+1)*nodes_per_task)] for i in range(math.floor(len(target_nodes)/nodes_per_task))]
        if max(chunks[-1]) != len(target_nodes)-1:
            i = math.floor(len(target_nodes)/nodes_per_task)
            chunks.append([target_nodes[j] for j in range(i*nodes_per_task,len(target_nodes))])

        polys = [None] * numNod
        
        def show_processbar(max_count, description, queue):
            # DEM2DP.init_progressbar(max_count, description)
            sum = 0
            while True:
                item = queue.get()
                if item < 0:
                    break
                sum += item
                DEM2DP.show_progress(sum)

        DEM2DP.init_progressbar(max_count, "generating polygons")
        shared_queue = SimpleQueue()
        pb_process = Process(target=show_processbar, args=(max_count, "generating polygons", shared_queue))
        pb_process.start()

        with Pool(ncores, initializer=DEM2DP.init_worker, initargs=(shared_queue,)) as pool:
            cnt = 0
            for i, ret_polys in enumerate(pool.imap(partial(DEM2DP.generate_nodal_poly, self, meshxy, bounds), chunks)):
                cnt += len(chunks[i])
                # DEM2DP.show_progress(cnt)

                for index, poly in enumerate(ret_polys):
                    polys[chunks[i][index]] = poly
        
        shared_queue.put(-1)
        pb_process.join()

        # DEM2DP.show_progress(max_count)
        print('done', flush=True)

        print('converting polygons to gdf', flush=True)
        target_nodes = [i for i in range(numNod) if polys[i]]
        target_gdfpolys = GeoDataFrame(index=list(range(len(target_nodes))), crs='epsg:4326', geometry=[polys[i] for i in target_nodes]).to_crs(src_crs)
        gdfpolys = GeoDataFrame(index=list(range(numNod)), columns=target_gdfpolys.columns)
        gdfpolys.iloc[target_nodes] = target_gdfpolys.iloc[0:target_gdfpolys.shape[0]]
        print('done', flush=True)

        # generate zonal stats node by node
        print('generating zonal stats', flush=True)
        max_count = len(target_nodes)
        chunk_size = chunk_size_zonalstats
        nodes_per_task = min(chunk_size, math.floor(len(target_nodes)/ncores))
        print('target_nodes length: {}, chunk_size: {}, ncores: {}, nodes_per_task: {}'.format(len(target_nodes), chunk_size, ncores, nodes_per_task))
        chunks = [[target_nodes[j] for j in range(i*nodes_per_task,(i+1)*nodes_per_task)] for i in range(math.floor(len(target_nodes)/nodes_per_task))]
        if max(chunks[-1]) != len(target_nodes)-1:
            i = math.floor(len(target_nodes)/nodes_per_task)
            chunks.append([target_nodes[j] for j in range(i*nodes_per_task,len(target_nodes))])

        self.zonal_stats_node = pd.DataFrame()
        
        DEM2DP.init_progressbar(max_count, "generating zonal stats")
        shared_queue = SimpleQueue()
        pb_process = Process(target=show_processbar, args=(max_count, "generating zonal stats", shared_queue))
        pb_process.start()

        self.zonal_stats_node = pd.DataFrame(index=list(range(numNod)), columns=['mean', 'max', 'min', 'count', 'lon', 'lat', 'area', 'hash'])

        with Pool(ncores, initializer=DEM2DP.init_worker, initargs=(shared_queue,)) as pool:
            cnt = 0
            for i, zonal_stats_node_list in enumerate(pool.imap(partial(DEM2DP.obtain_zonal_stats_for_node, self, tiffile, gdfpolys), chunks)):
                cnt += len(chunks[i])
                # DEM2DP.show_progress(cnt)

                for j, zonal_stats_node in enumerate(zonal_stats_node_list):
                    if type(zonal_stats_node) != type(None):
                        self.zonal_stats_node.loc[chunks[i][j], 'mean'] = zonal_stats_node.iloc[0]['mean']
                        self.zonal_stats_node.loc[chunks[i][j], 'max'] = zonal_stats_node.iloc[0]['max']
                        self.zonal_stats_node.loc[chunks[i][j], 'min'] = zonal_stats_node.iloc[0]['min']
                        self.zonal_stats_node.loc[chunks[i][j], 'count'] = zonal_stats_node.iloc[0]['count']
                        self.zonal_stats_node.loc[chunks[i][j], 'lon'] = zonal_stats_node.iloc[0]['lon']
                        self.zonal_stats_node.loc[chunks[i][j], 'lat'] = zonal_stats_node.iloc[0]['lat']
                        self.zonal_stats_node.loc[chunks[i][j], 'area'] = zonal_stats_node.iloc[0]['area']
                        self.zonal_stats_node.loc[chunks[i][j], 'hash'] = zonal_stats_node.iloc[0]['hash']
        # DEM2DP.show_progress(max_count)

        shared_queue.put(-1)
        pb_process.join()

        self.zonal_stats_node.reset_index()
        print('done', flush=True)

        return

    def save_cache(self, filename):
        cache = {}
        for index, row in self.zonal_stats_node.iterrows():
            if not np.isnan(row['lon']):
                lon = row['lon']
                lat = row['lat']
                area = row['area']
                hash = DEM2DP.gen_hash(lon, lat, area)
                cache[hash] = row

        try:
            with open(filename, 'wb') as file:
                pickle.dump(cache, file)
        except FileNotFoundError:
            print('failed to write cache file: {}'.format(filename), flush=True)

    def load_cache(self, filename):
        try:
            with open(filename, 'rb') as file:
                self.zonal_stats_node_cache = pickle.load(file)
        except IOError:
            print("failed to open cache file: {}".format(filename), flush=True)

    def assign_mesh_depth(self, min_count=100, ignore_channelnodes=False, \
                          land_only=False, submerged_only=False, \
                          min_depth=-100000.0, min_depth_tapering_end=0.0, \
                          max_depth=100000.0, \
                          deepen=False, \
                          channel_deeper_by=0.0, channel_deeper_by_threshold=10000.0, \
                          method='mean', ignore_tiff=False):
        numNod = self.mesh.numNod

        valtype = method
        print('- assigning {} values'.format(valtype), flush=True)

        # Assign depth
        for i in range(numNod):
            if self.target_loaded and not self.target[i]:
                continue
            val = -self.mesh.coord.loc[i, 'Depth']
            if i in self.mesh.channelnodes:
                if not ignore_channelnodes:
                    count = self.zonal_stats_node.loc[i, 'count']
                    if count >= min_count:
                        if ignore_tiff:
                            val = -self.mesh.coord.loc[i, 'Depth']
                        else:
                            val = self.zonal_stats_node.loc[i, valtype]
                        if (not deepen or (self.mesh.coord.loc[i, 'Depth'] < min(max(-val, min_depth), max_depth))) and \
                           (not land_only or val > 0) and \
                           (not submerged_only or (min(max(-val, min_depth), max_depth) > 0 and self.mesh.coord.loc[i, 'Depth'] > 0)):
                            if min_depth_tapering_end < min_depth:
                                min_depth_tapered = max(min_depth_tapering_end, min(min_depth, min_depth - (min_depth - min_depth_tapering_end) * (min_depth - -val) / (min_depth - min_depth_tapering_end)))
                            else:
                                min_depth_tapered = min_depth
                            self.mesh.coord.loc[i, 'Depth'] = min(max(-val, min_depth_tapered), max_depth)
            else:
                if ignore_tiff:
                    val = -self.mesh.coord.loc[i, 'Depth']
                else:
                    val = self.zonal_stats_node.loc[i, valtype]
                count = self.zonal_stats_node.loc[i, 'count']
                if (ignore_tiff or count >= min_count) and \
                   (not deepen or self.mesh.coord.loc[i, 'Depth'] < min(max(-val, min_depth), max_depth)) and \
                   (not land_only or val > 0) and \
                   (not submerged_only or (min(max(-val, min_depth), max_depth) > 0 and self.mesh.coord.loc[i, 'Depth'] > 0)):
                    if min_depth_tapering_end < min_depth:
                        min_depth_tapered = max(min_depth_tapering_end, min(min_depth, min_depth - (min_depth - min_depth_tapering_end) * (min_depth - -val) / (min_depth - min_depth_tapering_end)))
                    else:
                        min_depth_tapered = min_depth
                    self.mesh.coord.loc[i, 'Depth'] = min(max(-val, min_depth_tapered), max_depth)

        # Update the weir node heights
        # - open boundary
        cnt = 0
        line = self.mesh.boundaries[cnt].strip()
        items = re.split(' +', line)
        nopb = int(items[0])
        cnt += 1
        
        line = self.mesh.boundaries[cnt].strip()
        items = re.split(' +', line)
        nopbt = int(items[0])
        cnt += 1
        
        for i in range(nopb):
            line = self.mesh.boundaries[cnt].strip()
            items = re.split(' +', line)
            nopbi = int(items[0])
            cnt += 1
            cnt += nopbi

        # - land boundary
        line = self.mesh.boundaries[cnt].strip()
        items = re.split(' +', line)
        nlb = int(items[0])
        cnt += 1
        
        line = self.mesh.boundaries[cnt].strip()
        items = re.split(' +', line)
        nlbt = int(items[0])
        cnt += 1

        cnt_bak = cnt
        landbnd = np.zeros((numNod), int)
        for i in range(nlb):
            line = self.mesh.boundaries[cnt].strip()
            items = re.split(' +', line)
            nlbi = int(items[0])
            ibtype = int(items[1])
            cnt += 1
            if ibtype == 20 or ibtype == 21:
                for j in range(nlbi):
                    line = self.mesh.boundaries[cnt].strip()
                    items = re.split(' +', line)
                    n1 = int(items[0])-1
                    landbnd[n1] = 1
                    cnt += 1
            else:
                cnt += nlbi

        cnt = cnt_bak
        for i in range(nlb):
            line = self.mesh.boundaries[cnt].strip()
            items = re.split(' +', line)
            nlbi = int(items[0])
            ibtype = int(items[1])
            cnt += 1
            if ibtype == 64:
                for j in range(nlbi):
                    line = self.mesh.boundaries[cnt].strip()
                    items = re.split(' +', line)
                    n1 = int(items[0])-1
                    n2 = int(items[1])-1
                    if (j == 0 or j == nlbi - 1) and landbnd[n1] == 1:
                        self.mesh.coord.loc[n1, 'Depth'] = self.mesh.coord.loc[n1, 'Depth'] - 10.0 # raise the bank elevation at the river head to avoid instability
                    dp1 = self.mesh.coord.loc[n1, 'Depth']
                    dp2 = self.mesh.coord.loc[n2, 'Depth']
                    if not self.target_loaded or self.target[n2]:
                        if (not land_only) and \
                           (not submerged_only or min(max(dp2, min_depth), max_depth) > 0):
                            if channel_deeper_by > 0 and dp2 - dp1 < channel_deeper_by_threshold:
                                dp2 = max(min(dp1+channel_deeper_by,channel_deeper_by), dp2)
                            if min_depth_tapering_end < min_depth:
                                min_depth_tapered = \
                                    max(min_depth_tapering_end,
                                        min(min_depth, 
                                            min_depth - (min_depth - min_depth_tapering_end) * \
                                                (min_depth - dp2) / (min_depth - min_depth_tapering_end)))
                            else:
                                min_depth_tapered = min_depth
                            dp2 = min(max(dp2, min_depth_tapered), max_depth)
                    self.mesh.coord.loc[n2, 'Depth'] = max(dp1+1e-4, dp2) # make sure river bottom is as low as or lower than bank
                    items[2] = '{:.10f}'.format(-min(dp1, dp2) + 1e-3)
                    self.mesh.boundaries[cnt] = ' '.join(items) + '\n'
                    cnt += 1
            else:
                cnt += nlbi

        # make sure condensed nodes have the same depth value
        if hasattr(self,"f13") and self.f13.loaded:
            for i in range(self.f13.condensed_nodes.shape[0]):
                nodes = self.f13.condensed_nodes[i,:]
                maxdepth = -99999.9
                for node in nodes:
                    if node >= 0:
                        maxdepth = max(maxdepth, self.mesh.coord.loc[node, 'Depth'])
                for node in nodes:
                    if node >= 0:
                        self.mesh.coord.loc[node, 'Depth'] = maxdepth
            
                
if __name__ == '__main__':
    import argparse
    import os

    parser = argparse.ArgumentParser(
        prog='dem2adcdp',
        description='Map DEM onto ADCIRC mesh nodes.',
    )
    parser.add_argument('meshfile', action='store', type=str)
    parser.add_argument('outmeshfile', action='store', type=str)
    parser.add_argument('tiffile', action='store', type=str)
    parser.add_argument('--f13file', action='store', required=False, type=str, default='')
    parser.add_argument('--cachefile', action='store', required=False, type=str, default='')
    parser.add_argument('--assign_channelnode_depths', action='store_true', required=False, help='Force to assign depths at channel nodes')
    parser.add_argument('--target_add_all', action='store_true', required=False, help='Add all nodes to target list')
    parser.add_argument('--target_remove_all', action='store_true', required=False, help='Remove all nodes from target list')
    parser.add_argument('--target_add_by_polygons', action='store', required=False, type=str, default=[])
    parser.add_argument('--target_remove_by_polygons', action='store', required=False, type=str, default=[])
    parser.add_argument('--target_remove_by_polygons_outside', action='store', required=False, type=str, default=[])
    parser.add_argument('--target_add_channelnodes', action='store_true', required=False, help='Add channel nodes to target list')
    parser.add_argument('--target_remove_channelnodes', action='store_true', required=False, help='Remove channel nodes to target list')
    parser.add_argument('--target_add_banknodes', action='store_true', required=False, help='Add bank nodes to target list')
    parser.add_argument('--target_remove_banknodes', action='store_true', required=False, help='Remove bank nodes to target list')
    parser.add_argument('--target_add_boundarynodes', action='store_true', required=False, help='Add boundary nodes to target list')
    parser.add_argument('--channel_deeper_by', action='store', required=False, type=float, default=1e-3, help='Difference in depth by which channel node is deepened as compared to bank')
    parser.add_argument('--channel_deeper_by_threshold', action='store', required=False, type=float, default=10000, help='Threshold on the diff between current channel depth and bank node for applying --channel_deeper_by option.')
    parser.add_argument('--deepen', action='store_true', required=False, help='Assign an extracted value only when it is deeper than the current depth')
    parser.add_argument('--land_only', action='store_true', required=False, help='Assign depth only when the resulting elevation is above 0 m')
    parser.add_argument('--submerged_only', action='store_true', required=False, help='Assign depth only when both the current and resulting elevations are below 0 m')
    parser.add_argument('--min_depth', action='store', required=False, type=float, default=-100000.0, help='Minimum depth value to be assigned to mesh nodes')
    parser.add_argument('--min_depth_tapering_end', action='store', required=False, type=float, default=-99999.0, help='Depth at which tapering of the minimum depth application ends')
    parser.add_argument('--max_depth', action='store', required=False, type=float, default=100000.0, help='Maximum depth value to be assigned to mesh nodes')
    parser.add_argument('--method', choices=['mean', 'max', 'min'], required=False, default='mean', help='Extract method for elevation values (mean, max, min)')
    parser.add_argument('--ignore_tiff', action='store_true', required=False, help='Ignore values in tiff file and apply min_depth/max_depth only')
    parser.add_argument('--ncores', action='store', required=False, type=int, default=1)
    parser.add_argument('--chunk_size_poly', action='store', required=False, type=int, default=1000)
    parser.add_argument('--chunk_size_zonalstats', action='store', required=False, type=int, default=1000)

    args = parser.parse_args()

    if args.min_depth_tapering_end == -99999.0:
        args.min_depth_tapering_end = args.min_depth

    d2d = DEM2DP()

    d2d.read_mesh(args.meshfile)

    if os.path.isfile(args.f13file):
        d2d.read_f13(args.f13file)

    if os.path.isfile(args.cachefile):
        d2d.load_cache(args.cachefile)

    if args.target_add_all:
        d2d.add_all()

    if args.target_remove_all:
        d2d.remove_all()

    if args.target_add_by_polygons:
        polygonfiles = args.target_add_by_polygons.strip().split(',')
        for polygonfile in polygonfiles:
            d2d.add_by_polygons(polygonfile.strip())

    if args.target_add_channelnodes:
        d2d.add_channelnodes()

    if args.target_add_banknodes:
        d2d.add_banknodes()

    if args.target_add_boundarynodes:
        d2d.target_add_boundarynodes()

    if args.target_remove_by_polygons:
        polygonfiles = args.target_remove_by_polygons.strip().split(',')
        for polygonfile in polygonfiles:
            d2d.remove_by_polygons(polygonfile.strip())

    if args.target_remove_by_polygons_outside:
        polygonfiles = args.target_remove_by_polygons_outside.strip().split(',')
        for polygonfile in polygonfiles:
            d2d.remove_by_polygons_outside(polygonfile.strip())

    if args.target_remove_channelnodes:
        d2d.remove_channelnodes()

    if args.target_remove_banknodes:
        d2d.remove_banknodes()

    d2d.extract_dem_vals_on_mesh_nodes( \
        args.tiffile, ncores=args.ncores, chunk_size_poly=args.chunk_size_poly,
        chunk_size_zonalstats=args.chunk_size_zonalstats)

    print('assigning mesh depth', flush=True)
    d2d.assign_mesh_depth( \
        ignore_channelnodes=(not args.assign_channelnode_depths),
        land_only=args.land_only,
        submerged_only=args.submerged_only,
        min_depth=args.min_depth,
        min_depth_tapering_end=args.min_depth_tapering_end,
        max_depth=args.max_depth,
        deepen=args.deepen,
        channel_deeper_by=args.channel_deeper_by,
        channel_deeper_by_threshold=args.channel_deeper_by_threshold,
        min_count=1,
        method=args.method, ignore_tiff=args.ignore_tiff)

    d2d.write_mesh(args.outmeshfile)

    if args.cachefile:
        print('saving cache', flush=True)
        d2d.save_cache(args.cachefile)

    print('all done', flush=True)

    exit(0)
