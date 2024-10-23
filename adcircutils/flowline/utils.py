import pyproj

def get_utm_crs(west_lon_degree,
                south_lat_degree,
                east_lon_degree,
                north_lat_degree):
    utm_crs_list = pyproj.database.query_utm_crs_info(
        datum_name="WGS 84",
        area_of_interest=pyproj.aoi.AreaOfInterest(
            west_lon_degree=west_lon_degree,
            south_lat_degree=south_lat_degree,
            east_lon_degree=east_lon_degree,
            north_lat_degree=north_lat_degree,
        ),
    )
    utm_crs = pyproj.CRS.from_epsg(utm_crs_list[0].code)
    return utm_crs    

def utmzone(lon):
    zone = int((lon + 180.0) / 6) + 1
    return zone

def lonlat2utm(lon,lat,utm_crs):
    transformer = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True)
    point_utm = transformer.transform(lon, lat)

    return point_utm

def utm2lonlat(x,y,utm_crs):
    transformer = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True)
    point_lonlat = transformer.transform(x, y)

    return point_lonlat