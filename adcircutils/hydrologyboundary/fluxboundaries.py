from datetime import datetime
import pytz

class FluxBoundaries:
    def __init__(self):
        pass

    def initialize(self, boundary_defs, fluxboundary_seg_file = "", adcmesh_file = ""):
        self.seg = FluxBoundarySegments()
        
        if fluxboundary_seg_file:
            self.seg.read_fluxboundary_segment_file(fluxboundary_seg_file)
        elif adcmesh_file:
            self.seg.generate_fluxboundary_segments_from_adcmesh(adcmesh_file)
        else:
            raise Exception('Either fluxboundary_seg_file or adcmesh_file must be specified')
        
        if len(self.seg.segment_lengths) != len(boundary_defs):
            print('# of boundaries in seg = {}'.format(len(self.seg.segment_lengths)))
            print('# of boundaries in boundary_defs = {}'.format(len(boundary_defs)))
            print("Numbers of boundaries do not match.")
            return

        self.flux_boundaries = []

        for i in range(len(boundary_defs)):
            bnddef = boundary_defs[i]
            discharge_timeseries_type = bnddef[0]
            seg_lens = self.seg.segment_lengths[i]
            
            if discharge_timeseries_type == 'const_m3s':
                discharge_val = bnddef[1]
                fluxbnd = FluxBoundary(seg_lens, discharge_timeseries_type=discharge_timeseries_type, discharge_val=discharge_val)
            elif discharge_timeseries_type == 'const_m2s':
                discharge_val = bnddef[1]
                fluxbnd = FluxBoundary(seg_lens, discharge_timeseries_type=discharge_timeseries_type, discharge_val=discharge_val)
            elif discharge_timeseries_type == 'usgs':
                discharge_file = bnddef[1]
                fluxbnd = FluxBoundary(seg_lens, discharge_timeseries_type=discharge_timeseries_type, discharge_file=discharge_file)
            elif discharge_timeseries_type == 'contrail':
                discharge_file = bnddef[1]
                fluxbnd = FluxBoundary(seg_lens, discharge_timeseries_type=discharge_timeseries_type, discharge_file=discharge_file)
            elif discharge_timeseries_type == 'usgs_siteid':
                usgs_siteid = bnddef[1]
                fluxbnd = FluxBoundary(seg_lens, discharge_timeseries_type=discharge_timeseries_type, usgs_siteid=usgs_siteid)
            elif discharge_timeseries_type == 'nwm_featureid':
                nwm_featureid = bnddef[1]
                fluxbnd = FluxBoundary(seg_lens, discharge_timeseries_type=discharge_timeseries_type, nwm_featureid=nwm_featureid)
            else:
                raise Exception('discharge_timeseries_type {} is not supported.'.format(discharge_timeseries_type))

            self.flux_boundaries.append(fluxbnd)

    def init_boundary_defs(self):
        if not hasattr(self, "seg"):
            raise Exception("generate_fluxboundary_segments_from_adcmesh should be called after seg attribute is generated.")
        self.boundary_defs = [['NA','NA'] for i in range(len(self.seg.segment_lengths))]

    def set_boundary_defs_at_lonlat(self, lon, lat, discharge_type, discharge_val):
        from geopy.distance import geodesic, lonlat
        import numpy as np

        if not hasattr(self, "seg"):
            raise Exception("generate_fluxboundary_segments_from_adcmesh should be called after seg attribute is generated.")
        if not hasattr(self, "boundary_defs"):
            raise Exception("init_boundary_defs should be called after seg attribute is generated.")

        nsegs = len(self.seg.segment_lengths)
        dists = np.full((nsegs),100000.0)
        for i in range(nsegs):
            seg_lonlat = self.seg.segment_lonlats[i]
            dists[i] = geodesic(lonlat(*seg_lonlat), (lat, lon)).km
        idx = np.argmin(dists)
        seg_lonlat = self.seg.segment_lonlats[idx]
        print('Segment for index {:d} found at ({:.4f}, {:.4f}) at a distance of {:.0f} m.'.format(idx,seg_lonlat[0],seg_lonlat[1],dists[idx]*1000.0))
        self.boundary_defs[idx][0] = discharge_type
        self.boundary_defs[idx][1] = discharge_val
        

    def set_boundary_defs_at_NA(self, discharge_type, discharge_val):
        if not hasattr(self, "seg"):
            raise Exception("generate_fluxboundary_segments_from_adcmesh should be called after seg attribute is generated.")
        if not hasattr(self, "boundary_defs"):
            raise Exception("init_boundary_defs should be called after seg attribute is generated.")

        cnt = 0
        for i in range(len(self.boundary_defs)):
            if self.boundary_defs[i][0] == 'NA':
                self.boundary_defs[i][0] = discharge_type
                self.boundary_defs[i][1] = discharge_val
                cnt = cnt + 1

        print('{:d} NA segments were set values.'.format(cnt))

                         
    def set_boundary_defs_to_nwm_at_NA(self, hydrofabricfile, max_dist = 40.0):
        import pandas as pd
        import numpy as np
        import geopandas as gpd
        from shapely.geometry import Point

        if not hasattr(self, "seg"):
            raise Exception("generate_fluxboundary_segments_from_adcmesh should be called after seg attribute is generated.")
        if not hasattr(self, "boundary_defs"):
            raise Exception("init_boundary_defs should be called after seg attribute is generated.")

        gdf_hydrofabric = gpd.read_file(hydrofabricfile).to_crs('32633')

        seglonlats = [Point(self.seg.segment_lonlats[i][0], \
                            self.seg.segment_lonlats[i][1]) \
                        for i in range(len(self.seg.segment_lonlats))]
        gdf_segs = gpd.GeoDataFrame(data=pd.DataFrame(seglonlats, columns=['geometry']), geometry='geometry', crs='4326').to_crs('32633')

        segx = [seg.x for seg in gdf_segs.geometry]
        segy = [seg.y for seg in gdf_segs.geometry]

        tgt = [np.all([ \
                    np.array(gdf_hydrofabric.loc[i,'geometry'].xy[0]) > min(segx) - 10000.0, \
                    np.array(gdf_hydrofabric.loc[i,'geometry'].xy[0]) < max(segx) + 10000.0, \
                    np.array(gdf_hydrofabric.loc[i,'geometry'].xy[1]) > min(segy) - 10000.0, \
                    np.array(gdf_hydrofabric.loc[i,'geometry'].xy[1]) < max(segy) + 10000.0] \
                ) for i in range(len(gdf_hydrofabric))]

        gdf_tgt = gdf_hydrofabric[tgt].reset_index()

        nvert = 0
        for i in range(len(gdf_tgt)):
            nvert += max(1, len(gdf_tgt.loc[i,'geometry'].xy[0])-1)

        nwmvert = pd.DataFrame(index=range(nvert),columns=['reach','vert','geometry'])
        cnt = 0
        for i in range(len(gdf_tgt)):
            n = max(1,len(gdf_tgt.loc[i,'geometry'].xy[0]) - 1)
            nwmvert.loc[cnt:(cnt+n-1),'reach'] = np.array([i]*n)
            nwmvert.loc[cnt:(cnt+n-1),'vert'] = np.array(range(n))
            nwmvert.loc[cnt:(cnt+n-1),'geometry'] = [Point(gdf_tgt.loc[i,'geometry'].xy[0][j], \
                                                        gdf_tgt.loc[i,'geometry'].xy[1][j]) for j in range(n)]
            cnt += n

        def calculate_distance(Point, gdf_hydrofabric):
            distances = []
            for linestring in gdf_hydrofabric.geometry:
                distance = Point.distance(linestring)
                distances.append(distance)
            return distances

        cnt = 0
        for i in range(len(gdf_segs)):
            if self.boundary_defs[i][0] == 'NA':
                distances = calculate_distance(gdf_segs.loc[i,'geometry'], gdf_tgt)
                imin = np.argmin(distances)
                if distances[imin] < max_dist:
                    dist = distances[imin]
                    self.boundary_defs[i][0] = 'nwm_featureid'
                    self.boundary_defs[i][1] = gdf_tgt.loc[imin,'featureID']
                    cnt = cnt + 1
                    print('NWM reach {} found for segment {:d} at ({:.8f},{:.8f}). Distance: {:.0f} m.'.format(gdf_tgt.loc[imin,'featureID'], i, self.seg.segment_lonlats[i][0], self.seg.segment_lonlats[i][1], dist))
                else:
                    print('No NWM reach found for segment {:d}. No boundary def is assigined'.format(i))

        print('{:d} NA segments were set as nwm_featureid.'.format(cnt))


    def digest_boundary_defs(self):
        self.flux_boundaries = []

        for i in range(len(self.boundary_defs)):
            bnddef = self.boundary_defs[i]
            discharge_timeseries_type = bnddef[0]
            seg_lens = self.seg.segment_lengths[i]
            
            if discharge_timeseries_type == 'const_m3s':
                discharge_val = bnddef[1]
                fluxbnd = FluxBoundary(seg_lens, discharge_timeseries_type=discharge_timeseries_type, discharge_val=discharge_val)
            elif discharge_timeseries_type == 'const_m2s':
                discharge_val = bnddef[1]
                fluxbnd = FluxBoundary(seg_lens, discharge_timeseries_type=discharge_timeseries_type, discharge_val=discharge_val)
            elif discharge_timeseries_type == 'usgs':
                discharge_file = bnddef[1]
                fluxbnd = FluxBoundary(seg_lens, discharge_timeseries_type=discharge_timeseries_type, discharge_file=discharge_file)
            elif discharge_timeseries_type == 'contrail':
                discharge_file = bnddef[1]
                fluxbnd = FluxBoundary(seg_lens, discharge_timeseries_type=discharge_timeseries_type, discharge_file=discharge_file)
            elif discharge_timeseries_type == 'usgs_siteid':
                usgs_siteid = bnddef[1]
                fluxbnd = FluxBoundary(seg_lens, discharge_timeseries_type=discharge_timeseries_type, usgs_siteid=usgs_siteid)
            elif discharge_timeseries_type == 'nwm_featureid':
                nwm_featureid = bnddef[1]
                fluxbnd = FluxBoundary(seg_lens, discharge_timeseries_type=discharge_timeseries_type, nwm_featureid=nwm_featureid)
            else:
                raise Exception('discharge_timeseries_type {} is not supported.'.format(discharge_timeseries_type))

            self.flux_boundaries.append(fluxbnd)

    def read_fluxboundary_segment_file(self, fluxboundary_seg_file):
        self.seg = FluxBoundarySegments()
        self.seg.read_fluxboundary_segment_file(fluxboundary_seg_file)

    def generate_fluxboundary_segments_from_adcmesh(self, adcmesh_file):
        self.seg = FluxBoundarySegments()
        self.seg.generate_fluxboundary_segments_from_adcmesh(adcmesh_file)

    def set_discharge_from_usgs_siteids(self, start, end):
        for fluxbnd in self.flux_boundaries:
            if fluxbnd.discharge_timeseries_type == 'usgs_siteid':
                fluxbnd.set_discharge_from_usgs_siteid(start, end)

    def set_discharge_from_nwm_featureids(self, start, end):
        from adcircutils.hydrologyboundary import nwm
        
        nwm_featureids = []
        for fluxbnd in self.flux_boundaries:
            if fluxbnd.discharge_timeseries_type == 'nwm_featureid':
                try:
                    nwm_featureid = fluxbnd.nwm_featureid
                except:
                    raise Exception('Exception raised during NWM feature ID collection')

                nwm_featureids.append(nwm_featureid)

        if nwm_featureids:
            dt = 3600
            # Select repository to use:
                    #1: "https://nwm-archive.s3.amazonaws.com/"   #NWM v1.2 1993-2020
                    #2: "https://noaa-nwm-retrospective-2-1-pds.s3.amazonaws.com/?prefix=model_output/"   #NWM v2.1 1993-2020
                    #3: "" # NWMv2.1 2020-2022
            repository_version = 4

            # Cache directory for the NWM output netCDF files
            cachedir = '/mnt/d/work/NWM_Coupling/nwm_res_cache/v3.0'

            # Get NWM outputs
            timqs_nwm, valqs_nwm = nwm.getDischarge(repository_version, cachedir, start, end, dt, nwm_featureids)

            # Set discharge values to each boundary
            cnt = 0
            for ibnd in range(len(self.flux_boundaries)):
                if self.flux_boundaries[ibnd].discharge_timeseries_type == 'nwm_featureid':
                    self.flux_boundaries[ibnd].set_discharge_values(timqs_nwm, valqs_nwm[cnt])
                    cnt = cnt + 1

    def set_discharge_from_database(self, start, end):
        # Set discharge from USGS site IDs
        self.set_discharge_from_usgs_siteids(start, end)
                
        # Set discharge from NWM feature IDs
        self.set_discharge_from_nwm_featureids(start, end)

    
    def write_fort20(self, fort20fname, start, end, step):        
        # Write fort.20
        with open(fort20fname, "w") as fp:
            fp.write('{:f}\n'.format(step.total_seconds()))

            dt = start
            while dt <= end:
                for fluxbnd in self.flux_boundaries:
                    discharge_m2s = fluxbnd.get_discharge_at_time_in_m2s(dt)
                    for i in range(len(discharge_m2s)):
                        fp.write('{:.6f}\n'.format(discharge_m2s[i]))
                dt = dt + step
                
    def plot_discharge_timeseries_m3s(self, start, end, step):
        for i in range(len(self.flux_boundaries)):
            self.flux_boundaries[i].plot_discharge_timeseries_m3s(start, end, step)

class FluxBoundarySegments:
    def __init__(self):
        pass

    def read_fluxboundary_segment_file(self, fname):
        self.nboundaries = 0
        self.segment_lengths = []
        with open(fname, "r") as fp:
            line = fp.readline()
            items = line.split()
            self.nboundaries = int(items[0])
            for i in range(self.nboundaries):
                line = fp.readline()
                pos = line.find("#")
                line = line[0:pos]
                line = line.strip()
                items = line.split()
                seglen = list(map(float, items))
                self.segment_lengths.append(seglen)

    def generate_fluxboundary_segments_from_adcmesh(self, adcmesh_file):
        # import matlab.engine
        # eng = matlab.engine.start_matlab()
        # seglens = eng.get_fluxboundary_segment_lengths_and_lonlats(adcmesh_file)
        # self.segment_lengths = []
        # self.segment_lonlats = []
        # for i in range(len(seglens)):
        #     sl = [seglens[i][0][j] for j in range(len(seglens[i][0]))]
        #     self.segment_lengths.append(sl[0:-2])
        #     self.segment_lonlats.append(sl[-2:])

        import numpy as np
        from adcircpy import AdcircMesh
        import geopy.distance

        mesh = AdcircMesh.open(adcmesh_file)
        
        self.segment_lengths = []
        self.segment_lonlats = []
        for stnodes in mesh.inflow_boundaries.node_id:
            nodes = [int(stnode) for stnode in stnodes]
            xs = mesh.x[nodes].values
            ys = mesh.y[nodes].values

            seglens = []
            for i in range(1,len(xs)):
                x1 = xs[i-1]
                x2 = xs[i]
                y1 = ys[i-1]
                y2 = ys[i]

                seglen = geopy.distance.geodesic((y1, x1), (y2, x2)).meters
                seglens.append(seglen)
                
            self.segment_lengths.append(seglens)
            self.segment_lonlats.append([np.mean(xs), np.mean(ys)])


class FluxBoundary:
    def __init__(self, segment_lengths, ibtype = 20, discharge_timeseries_type = 'const_m3s', discharge_val = 0.0, discharge_file = "",\
                usgs_siteid = "", nwm_featureid = ""):
        if len(segment_lengths) > 0:
            self.segment_lengths = segment_lengths
            self.tot_length = sum(segment_lengths)
            self.nseg_nodes = len(segment_lengths) + 1
            self.ibtype = ibtype
            self.discharge_timeseries_type = discharge_timeseries_type

            if discharge_timeseries_type == 'const_m3s':
                self.discharge_timeseries = float(discharge_val)
            elif discharge_timeseries_type == 'const_m2s':
                self.discharge_timeseries = list(map(float,discharge_val))
            elif discharge_timeseries_type == 'usgs':
                self.set_discharge_from_usgs_file(discharge_file)
            elif discharge_timeseries_type == 'contrail':
                self.set_discharge_from_contrail_file(discharge_file)
            elif discharge_timeseries_type == 'usgs_siteid':
                if usgs_siteid:
                    self.usgs_siteid = usgs_siteid
                else:
                    raise Exception('usgs_siteid is not specified.')
            elif discharge_timeseries_type == 'nwm_featureid':
                if nwm_featureid:
                    self.nwm_featureid = nwm_featureid
                else:
                    raise Exception('nwm_featureid is not specified.')
            else:
                raise Exception('discharge_timeseries_type {} is not supported.'.format(discharge_timeseries_type))
            
            self.segnode_lengths = [self.segment_lengths[0]]
            for i in range(len(segment_lengths)-1):
                l0 = self.segment_lengths[i]
                l1 = self.segment_lengths[i+1]
                self.segnode_lengths.append(0.5*(l0+l1))
            self.segnode_lengths.append(self.segment_lengths[-1])

        else:
            raise Exception("Empty boundary segment list")


    def set_discharge_from_usgs_file(self, usgs_filename):
        self.discharge_timeseries = []
        with open(usgs_filename, "r") as fp:
            line = fp.readline()
            while line.startswith("#"):
                line = fp.readline()
            if not line.startswith("agency_cd"):
                print('line = {:s}'.format(line))
                raise Exception("Format error (1)")
            line = fp.readline()
            if not line.startswith("5s"):
                print('line = {:s}'.format(line))
                raise Exception("Format error (2)")
            first = False
            line = fp.readline()
            while line:
                items = line.split()
                agency_cd = items[0]
                site_no = items[1]
                sdatetime = items[2] + " " + items[3]
                stz = items[4]
                discharge_cfs = float(items[5])

                dt = datetime.strptime(sdatetime,"%Y-%m-%d %H:%M")
                if stz == 'EDT':
                    stz = 'EST5EDT'
                tz = pytz.timezone(stz)

                dt_tz = tz.localize(dt)
                dt_utc = dt_tz.astimezone(pytz.utc)

                feet2m = 0.3048
                discharge_m3s = feet2m**3 * discharge_cfs 

                # print('{:.3f} # {:s}'.format(discharge_m3s,datetime.strftime(dt_utc,"%Y-%m-%d %H:%M %Z")))
                
                self.discharge_timeseries.append([dt_utc,discharge_m3s])

                line = fp.readline()

    def set_discharge_from_contrail_file(self, contrail_filename):
        self.discharge_timeseries = []
        with open(contrail_filename, "r") as fp:
            line = fp.readline() # header
            line = fp.readline()
            while line:
                items = line.split(',')
                sdatetime = items[0]
                stz = 'EDT'
                if len(items) == 5:
                    discharge_cfs = float(items[2])
                elif len(items) == 6:
                    discharge_cfs = float(items[2]+items[3])
                else:
                    raise Exception("Format error contrail file (1)")

                dt = datetime.strptime(sdatetime,"%Y-%m-%d %H:%M:%S")
                if stz == 'EDT':
                    stz = 'EST5EDT'
                tz = pytz.timezone(stz)

                dt_tz = tz.localize(dt)
                dt_utc = dt_tz.astimezone(pytz.utc)

                feet2m = 0.3048
                discharge_m3s = feet2m**3 * discharge_cfs 

                # print('{:.3f} # {:s}'.format(discharge_m3s,datetime.strftime(dt_utc,"%Y-%m-%d %H:%M %Z")))
                
                self.discharge_timeseries.append([dt_utc,discharge_m3s])

                line = fp.readline()
            
            self.discharge_timeseries.reverse()


    def set_discharge_from_usgs_siteid(self, start, end):
        from adcircutils.hydrologyboundary import usgs
        from datetime import timedelta

        # Cache directory for the NWM output netCDF files
        cachedir = '/mnt/d/work/NWM_Coupling/usgs_nwis_cache'

        timq, valq, timh, valh \
            = usgs.getQH(self.usgs_siteid, start, end+timedelta(hours=1), cachedir)

        self.discharge_timeseries = [[timq[0][i],valq[0][i]] for i in range(len(timq[0]))]


    def set_discharge_values(self, timq, valq):
        self.discharge_timeseries = [[timq[i],valq[i]] for i in range(len(timq))]


    def get_discharge_at_time_in_m2s(self, dt):
        if type(dt) == list:
            discharges = []
            for t in dt:
                discharge = self.get_discharge_at_time_in_m2s(t)
                discharges.append(discharge)
            return discharges

        if self.discharge_timeseries_type == 'const_m3s':
            if type(self.discharge_timeseries) == float:
                return [self.discharge_timeseries/self.tot_length] * self.nseg_nodes
            else:
                raise Exception("Type of discharge_timeseries must be float for 'const_m3s' type.") 

        elif self.discharge_timeseries_type == 'const_m2s':
            if type(self.discharge_timeseries) == list:
                return self.discharge_timeseries
            else:
                raise Exception("Type of discharge_timeseries must be list for 'const_2s' type.") 

        elif self.discharge_timeseries_type == 'usgs' or self.discharge_timeseries_type == 'contrail' or \
             self.discharge_timeseries_type == 'usgs_siteid' or self.discharge_timeseries_type == 'nwm_featureid':
            found = False
            
            try:
                dummy = self.count
            except:
                self.count = 0

            for i in range(len(self.discharge_timeseries)):
                idx = (i + self.count) % len(self.discharge_timeseries)
                if idx == len(self.discharge_timeseries):
                    continue
                dt0 = self.discharge_timeseries[idx][0]
                dt1 = self.discharge_timeseries[idx+1][0]
                if dt >= dt0 and dt <= dt1:
                    found = True
                    break
            if not found:
                print('datetime start     = {}'.format(self.discharge_timeseries[0][0]))
                print('datetime end       = {}'.format(self.discharge_timeseries[-1][0]))
                print('datetime specified = {}'.format(dt))
                raise Exception("Specified datetime is out of range.")
            
            self.count = idx
            i = idx

            delta = dt1 - dt0
            delta_at = dt - dt0
            r = delta_at.total_seconds() / delta.total_seconds()

            discharge0 = self.discharge_timeseries[i][1]
            discharge1 = self.discharge_timeseries[i+1][1]
            discharge = [((1-r)*discharge0 + r*discharge1)/self.tot_length] * self.nseg_nodes

            return discharge
        else:
            raise Exception('discharge_timeseries_type {} is not supported.'.format(self.discharge_timeseries_type))


    def get_discharge_at_time_in_m3s(self, dt):
        if type(dt) == list:
            discharges = []
            for t in dt:
                discharge = self.get_discharge_at_time_in_m3s(t)
                discharges.append(discharge)
            return discharges
                
        if self.discharge_timeseries_type == 'const_m3s':
            if type(self.discharge_timeseries) == float:
                return self.discharge_timeseries
            else:
                raise Exception("Type of discharge_timeseries must be float for 'const_m3s' type.") 

        elif self.discharge_timeseries_type == 'const_m2s':
            if type(self.discharge_timeseries) == list:
                dist = 0.0
                for i in range(len(self.discharge_timeseries)):
                    dist = dist + self.discharge_timeseries[i]*self.segnode_lengths[i]
                return dist
            else:
                raise Exception("Type of discharge_timeseries must be list for 'const_2s' type.") 

        elif self.discharge_timeseries_type == 'usgs' or self.discharge_timeseries_type == 'contrail' or \
             self.discharge_timeseries_type == 'usgs_siteid' or self.discharge_timeseries_type == 'nwm_featureid':
            found = False
            
            try:
                dummy = self.count
            except:
                self.count = 0

            for i in range(len(self.discharge_timeseries)):
                idx = (i + self.count) % len(self.discharge_timeseries)
                if idx == len(self.discharge_timeseries) - 1:
                    continue
                dt0 = self.discharge_timeseries[idx][0]
                dt1 = self.discharge_timeseries[idx+1][0]
                if dt >= dt0 and dt <= dt1:
                    found = True
                    break
            if not found:
                print('datetime start     = {}'.format(self.discharge_timeseries[0][0]))
                print('datetime end       = {}'.format(self.discharge_timeseries[-1][0]))
                print('datetime specified = {}'.format(dt))
                raise Exception("Specified datetime is out of range.")
            
            self.count = idx
            i = idx

            delta = dt1 - dt0
            delta_at = dt - dt0
            r = delta_at.total_seconds() / delta.total_seconds()

            discharge0 = self.discharge_timeseries[i][1]
            discharge1 = self.discharge_timeseries[i+1][1]
            discharge = (1-r)*discharge0 + r*discharge1

            return discharge
        else:
            raise Exception('discharge_timeseries_type {} is not supported.'.format(self.discharge_timeseries_type))

    def plot_discharge_timeseries_m3s(self, start, end, step):
        import matplotlib.pyplot as plt

        dts = []
        dt = start
        while dt <= end:
            dts.append(dt)
            dt = dt + step

        valq = self.get_discharge_at_time_in_m3s(dts)

        fig, ax = plt.subplots()
        ax.plot(dts,valq)
        ax.set_ylabel('Discharge (m$^3$/s)')
        ax.grid()
        ax.set_xlim([start,end])

