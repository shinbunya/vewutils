def searchNearestFeatureIDAtLonLat(hydrofabricfile,lon,lat):
    pass

def getDischarge(repository_version, init_date, end_date, dt, dest_feature_ids):
    import os
    import sys
    from pathlib import Path
    import numpy as np
    import wget
    from datetime import datetime
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config

    # init_date = datetime.strptime(init, '%Y%m%d %H:%M:%S')
    # end_date = datetime.strptime(end, '%Y%m%d %H:%M:%S')
    df = (end_date-init_date).total_seconds()
    ndt = int(np.round(df/dt) + 1)

    import datetime
    import netCDF4 as nc

    Q = np.zeros(shape=(len(dest_feature_ids),len(range(ndt))))
    t_name = [init_date + datetime.timedelta(seconds=i*dt) for i in range(ndt)]

    if repository_version == 1.2:
        # url='https://nwm-archive.s3.amazonaws.com/'   #NWM v1.2 retrospective
        # for i in range(ndt):
        #     #t_name[i] = init_date + datetime.timedelta(seconds=i*dt)
        #     yy = t_name[i].strftime('%Y')
        #     filename =  t_name[i].strftime('%Y%m%d%H%M')+'.CHRTOUT_DOMAIN1.comp'
        #     path = url + yy + '/' + filename 
        #     try:
        #         aux = wget.download(path)
        #     except:
        #         aux = wget.download(path)
        #         pass  
        #     ds = nc.Dataset(aux)
        #     feature_id = ds['feature_id'][:]
        #     lut = {v:i for i, v in enumerate(feature_id.data)}
        #     Q[:,i] = ds['streamflow'][ids]
        #     os.remove(aux)
        url='nwm-archive'   #NWM v1.2 retrospective
        print('Retrieving data from {:s}'.format(url))
        s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
        for i in range(ndt):
            #t_name[i] = init_date + datetime.timedelta(seconds=i*dt)
            yy = t_name[i].strftime('%Y')
            filename =  t_name[i].strftime('%Y%m%d%H%M')+'.CHRTOUT_DOMAIN1.comp'
            path = yy + '/' + filename
            cachedir = os.environ["NWM12_RES_CACHE_DIR"] 
            cache_filename = cachedir + '/' + path
            
            if Path(cache_filename).is_file():
                print('Found in local cache: {:s} ({:d}/{:d})'.format(path,(i+1),ndt))
            else:
                print('Downloading {:s} ({:d}/{:d})'.format(path,(i+1),ndt))
                Path(cachedir + '/' + yy).mkdir(parents=True, exist_ok=True)
                try:
                    s3.download_file(url,path,cache_filename)
                except:
                    s3.download_file(url,path,cache_filename)
                    pass
            ds = nc.Dataset(cache_filename)
            if i == 0:
                feature_id = ds['feature_id'][:]
                lut = {v:i for i, v in enumerate(feature_id.data)}
                ids = [lut[f] for f in dest_feature_ids]
            Q[:,i] = [ds['streamflow'][id] for id in ids]
            ds.close()
    elif repository_version == 2.1:
        url='noaa-nwm-retrospective-2-1-pds'   #NWM v2.1 retrospective
        print('Retrieving data from {:s}'.format(url))
        s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
        for i in range(ndt):
            #t_name[i] = init_date + datetime.timedelta(seconds=i*dt)
            yy = t_name[i].strftime('%Y')
            filename =  t_name[i].strftime('%Y%m%d%H%M')+'.CHRTOUT_DOMAIN1.comp'
            path = 'model_output/'+ yy + '/' + filename
            cachedir = os.environ["NWM21_RES_CACHE_DIR"] 
            cache_filename = cachedir + '/' + path
            
            if Path(cache_filename).is_file():
                print('Found in local cache: {:s} ({:d}/{:d})'.format(path,(i+1),ndt))
            else:
                print('Downloading {:s} ({:d}/{:d})'.format(path,(i+1),ndt))
                Path(cachedir + '/' + 'model_output/'+ yy).mkdir(parents=True, exist_ok=True)
                try:
                    s3.download_file(url,path,cache_filename)
                except:
                    s3.download_file(url,path,cache_filename)
                    pass
            ds = nc.Dataset(cache_filename)
            if i == 0:
                feature_id = ds['feature_id'][:]
                lut = {v:i for i, v in enumerate(feature_id.data)}
                ids = [lut[f] for f in dest_feature_ids]
            Q[:,i] = [ds['streamflow'][id] for id in ids]
            ds.close()
    elif repository_version == 3.0:
        url='noaa-nwm-retrospective-3-0-pds'   #NWM v3.0 retrospective
        print('Retrieving data from {:s}'.format(url))
        s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
        for i in range(ndt):
            #t_name[i] = init_date + datetime.timedelta(seconds=i*dt)
            yy = t_name[i].strftime('%Y')
            filename =  t_name[i].strftime('%Y%m%d%H%M')+'.CHRTOUT_DOMAIN1'
            path = 'CONUS/netcdf/CHRTOUT/'+ yy + '/' + filename 
            cachedir = os.environ["NWM30_RES_CACHE_DIR"] 
            cache_filename = cachedir + '/' + path
            
            if Path(cache_filename).is_file():
                print('Found in local cache: {:s} ({:d}/{:d})'.format(path,(i+1),ndt))
            else:
                print('Downloading {:s} ({:d}/{:d})'.format(path,(i+1),ndt))
                Path(cachedir + '/' + 'CONUS/netcdf/CHRTOUT/'+ yy).mkdir(parents=True, exist_ok=True)
                try:
                    s3.download_file(url,path,cache_filename)
                except:
                    s3.download_file(url,path,cache_filename)
                    pass
            ds = nc.Dataset(cache_filename)
            if i == 0:
                feature_id = ds['feature_id'][:]
                lut = {v:i for i, v in enumerate(feature_id.data)}
                ids = [lut[f] for f in dest_feature_ids]
            Q[:,i] = [ds['streamflow'][id] for id in ids]
            ds.close()
    else:
        sys.exit('Unrecognized repository option')

    return t_name, Q
            
