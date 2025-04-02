# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 08:50:00 2024
Module to read and write mesh and f13 files for the ADCIRC model

@author: shinbunya
"""

import sys
import argparse
import pandas as pd
import geopandas as gpd
import matplotlib.tri as mtri
import numpy as np
from shapely import Polygon
from shapely.geometry import MultiPolygon
from shapely.ops import unary_union
from rasterstats import zonal_stats
import re


class Mesh:
    """Class to store info in mesh file"""
    def __init__(self):
        self.loaded = False

    def read(self, meshFilename, msgout=True):
        if msgout:
            print("reading mesh file: {}".format(meshFilename), flush=True)

        # initialize variables 
        
        x = []
        y = []
        z = []
        nodNum = []
        eleNum = []
        
        with open(meshFilename) as gridFile:
    
            header = gridFile.readline()
            line = re.split(' +', gridFile.readline().strip())
            numEle = int(line[0])
            numNod = int(line[1])
            
            triangles = []

            # import coordinates of points and elevations
            
            for i in range(numNod):
                
                line = gridFile.readline().strip().split()
                nodNum.append(int(line[0]))
        
                x.append(float(line[1]))
                y.append(float(line[2]))
                z.append(float(line[3]))
        
            # import vertex order for each element
            # NOTE: the -1 in the triangles assembly is to make it 0 indexed
            
            for i in range(numEle):
                line = re.split(' +', gridFile.readline().strip())
                eleNum.append(int(line[0]))
                triangles.append([int(line[2])-1,int(line[3])-1,int(line[4])-1])

            boundaries = gridFile.readlines()

        # boundaries
        # - open boundary
        cnt = 0
        line = boundaries[cnt].strip()
        items = re.split(' +', line)
        nopb = int(items[0])
        cnt += 1
        
        line = boundaries[cnt].strip()
        items = re.split(' +', line)
        nopbt = int(items[0])
        cnt += 1
        
        for i in range(nopb):
            line = boundaries[cnt].strip()
            items = re.split(' +', line)
            nopbi = int(items[0])
            cnt += 1
            cnt += nopbi

        # - land boundary
        line = boundaries[cnt].strip()
        items = re.split(' +', line)
        nlb = int(items[0])
        cnt += 1
        
        line = boundaries[cnt].strip()
        items = re.split(' +', line)
        nlbt = int(items[0])
        cnt += 1

        boundarynodes = []
        banknodes = []
        channelnodes = []
        boundaries64 = []
        for i in range(nlb):
            line = boundaries[cnt].strip()
            items = re.split(' +', line)
            nlbi = int(items[0])
            ibtype = int(items[1])
            cnt += 1
            boundary64 = []
            for j in range(nlbi):
                line = boundaries[cnt].strip()
                items = re.split(' +', line)
                boundarynodes.append(int(items[0])-1)
                if ibtype == 64:
                    banknodes.append(int(items[0])-1)
                    channelnodes.append(int(items[1])-1)
                    boundary64.append(int(items[1])-1)
                cnt += 1
            if len(boundary64) > 0:
                boundaries64.append(boundary64)
                                
        # data conversion
        triang = mtri.Triangulation(x,y,triangles)
                
        # put xyz into dataframe for ease of use
        gridXYZ = pd.DataFrame({'Node Number' : nodNum,
                     'Latitude' : y,
                     'Longitude' : x,
                     'Depth' : z})
        
        # construct node2elem
        n2e = [[]] * numNod
        for j in range(numEle):
            tri = triang.triangles[j,:]
            n2e[tri[0]] = n2e[tri[0]] + [j]
            n2e[tri[1]] = n2e[tri[1]] + [j]
            n2e[tri[2]] = n2e[tri[2]] + [j]

        self.header = header
        self.coord = gridXYZ
        self.tri = triang.triangles
        self.numNod = numNod
        self.numEle = numEle
        self.boundarynodes = boundarynodes
        self.channelnodes = channelnodes
        self.banknodes = banknodes
        self.boundaries = boundaries
        self.boundaries64 = boundaries64
        self.n2e = n2e
        self.loaded = True
 
    def generate_boundary_polygon(self):
        print('-- extracting edges', flush=True)
        edges = []
        elems = []
        for j, triangle in enumerate(self.tri):
            for i in range(3):
                edge = tuple(sorted([triangle[i], triangle[(i + 1) % 3]]))
                edges.append(edge)
                elems.append((j,i))

        print('-- extracting boundary edges', flush=True)
        ndigits = len(str(self.numNod))            
        def hash(edge, ndigits):
            return "{i0:0{ndigits}}{i1:0{ndigits}}".format(i0=edge[0], i1=edge[1], ndigits=ndigits)
        dedges = {}
        for edge in edges:
            h = hash(edge, ndigits)
            dedges[h] = 0
        for edge in edges:
            h = hash(edge, ndigits)
            dedges[h] = dedges[h] + 1
        unshared_edges = []
        unshared_elems = []
        for i, edge in enumerate(edges):
            if dedges[hash(edge, ndigits)] == 1:
                unshared_edges.append(edge)
                unshared_elems.append(elems[i])

        print('-- constructing boundary polygons', flush=True)
        polys = []
        set_edges = set(list(range(len(unshared_edges))))
        set_edges_consumed = set()
        while len(set_edges) != len(set_edges_consumed):
            print('--- polygon {}'.format(len(polys)+1), flush=True)
            set_edges_notconsumed = set_edges.difference(set_edges_consumed)
            start_edge = list(set_edges_notconsumed)[0]
            edges_consumed = [start_edge]
            ielem = unshared_elems[start_edge][0]
            iedge = unshared_elems[start_edge][1]
            bndnodes = [self.tri[ielem][iedge], self.tri[ielem][(iedge+1)%3]]

            while bndnodes[0] != bndnodes[-1]:
                for i, edge in enumerate(unshared_edges):
                    ielem = unshared_elems[i][0]
                    iedge = unshared_elems[i][1]
                    if self.tri[ielem][iedge] == bndnodes[-1]:
                        inode = self.tri[ielem][(iedge+1)%3]
                        bndnodes.append(inode)
                        edges_consumed.append(i)
                        break
                
            set_edges_consumed.update(set(edges_consumed))

            if len(bndnodes) < 3:
                break
            
            poly = Polygon([(lon, lat) for lon, lat in zip(self.coord.loc[bndnodes,'Longitude'], self.coord.loc[bndnodes,'Latitude'])])
            polys.append(poly)

        print('-- merging polygons', flush=True)
        self.boundary_polygon = unary_union(polys)

    def generate_boundary64_polygons(self):
        print('-- constructing boundary polygons', flush=True)
        polys = []
        for boundary64 in self.boundaries64:
            poly = Polygon([(lon, lat) for lon, lat in zip(self.coord.loc[boundary64, 'Longitude'], self.coord.loc[boundary64, 'Latitude'])])
            polys.append(poly)
        self.boundary64_polygons = polys

    def write(self, mesh_filename, msgout=True):
        if msgout:
            print("writing mesh file: {}".format(mesh_filename), flush=True)

        with open(mesh_filename, 'w') as gridFile:
            gridFile.write(self.header)
            gridFile.write('{}  {}\n'.format(self.numEle, self.numNod))
            for i in range(self.numNod):
                gridFile.write('{:10d} {:15.10f} {:15.10f} {:15.10f}\n'.format(i+1, self.coord['Longitude'][i], self.coord['Latitude'][i], self.coord['Depth'][i]))
            for j in range(self.numEle):
                gridFile.write('{:d} 3 {:d} {:d} {:d}\n'.format(j+1, self.tri[j,0]+1, self.tri[j,1]+1, self.tri[j,2]+1))
            gridFile.writelines(''.join(self.boundaries).strip())

    def read_boundary_polygon(self, filename, msgout=True):
        if msgout:
            print("reading boundary polygon file: {}".format(filename), flush=True)

        gdf = gpd.read_file(filename)
        self.boundary_polygon = gdf.geometry[0]
        
    def write_boundary_polygon(self, filename, msgout=True):
        if msgout:
            print("writing boundary polygon file: {}".format(filename), flush=True)

        if not hasattr(self, 'boundary_polygon'):
            print('- boundary polygon is not available', flush=True)
            return
        
        gdf = gpd.GeoDataFrame(geometry=[self.boundary_polygon], crs='EPSG:4326')
        gdf.to_file(filename)
        
    def write_boundary64_polygons(self, filename, msgout=True):
        if msgout:
            print("writing boundary64 multipolygon file: {}".format(filename), flush=True)

        if not hasattr(self, 'boundary64_polygons'):
            print('- boundary64 polygons are not available', flush=True)
            return
        
        gdf = gpd.GeoDataFrame(geometry=[MultiPolygon(self.boundary64_polygons)], crs='EPSG:4326')
        gdf.to_file(filename)
        
class F13:
    """Class to store info in f13 file"""
    def __init__(self):
        self.loaded = False

    def read(self, f13Filename, msgout=True):
        if msgout:
            print("reading f13 file: {}".format(f13Filename), flush=True)

        with open(f13Filename) as file:
    
            header = file.readline()
            items = re.split(' +', file.readline().strip())
            numNod = int(items[0])

            found = False
            while True:
                line = file.readline().strip()
                if line == 'condensed_nodes':
                    found = True
                    break
                if not line:
                    break

            if not found:
                print("- condensed nodes not found", flush=True)
                return

            line = file.readline() # unit
            items = re.split(' +', file.readline().strip())
            numVars = int(items[0])
            defvals = [int(s) for s in re.split(' +', file.readline().strip())]
            
            while True:
                line = file.readline().strip()
                if line == 'condensed_nodes':
                    break
            
            items = re.split(' +', file.readline().strip())
            numPts = int(items[0])
            self.condensed_nodes = np.zeros((numPts,numVars + 1), int)
            
            for i in range(numPts):
                nodes = [int(s) - 1 for s in re.split(' +', file.readline().strip())]
                self.condensed_nodes[i,:] = nodes
        
            self.loaded = True


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog='adcmesh',
        description='''Tool to deal with ADCIRC mesh file.

    Typical usage: adcmesh generate boundary_polygon <input mesh file> <output file>
    ''')

    subparsers = parser.add_subparsers(dest='command')

    boundary_polygon_parser = subparsers.add_parser('generate')
    boundary_polygon_parser.add_argument('target', action='store', type=str)
    boundary_polygon_parser.add_argument('meshfile', action='store', type=str)
    boundary_polygon_parser.add_argument('output', action='store', type=str)

    args = parser.parse_args()

    if args.command == 'generate' and args.target == 'boundary_polygon':
        mesh = Mesh()
        mesh.read(args.meshfile)
        mesh.generate_boundary_polygon()
        mesh.write_boundary_polygon(args.output)

    if args.command == 'generate' and args.target == 'boundary_polygon_with_holes':
        mesh = Mesh()
        mesh.read(args.meshfile)
        mesh.generate_boundary_polygon_with_holes()
        mesh.write_boundary_polygon_with_holes(args.output)

    
