import geopandas as gpd
# import folium
from scipy.spatial import cKDTree
import numpy as np
import pandas as pd
from shapely.geometry import Point
from shapely.strtree import STRtree
import h5py
from datetime import datetime, timedelta
import netCDF4 as nc

# TODO make this script a function and each of the following inputs function parameters

# hiWayPts = gpd.read_file(r"Z:\GIS\Data_Highway\points_highway_4326.shp")
# hdf_filename = "Z:\SLaMM_Mar2022\SLaMM_2021.p06.hdf"
# output_nc = 'RAS_PointLocations.nc'
# stations_locations_txt = r"./system/postproc/ras/station_locations.txt"
# wse_gdf = gpd.read_file(r"C:\py\IIHR_RASpy\post_process_GIS\20220512_10\tempfiles\test_polygon_all_data.shp")


# Function to get the nearest point in B by location for each point in A.
def ckdnearest(gdA, gdB):
    nA = np.array(list(gdA.geometry.apply(lambda x: (x.x, x.y))))
    nB = np.array(list(gdB.geometry.apply(lambda x: (x.x, x.y))))
    btree = cKDTree(nB)
    dist, idx = btree.query(nA, k=1)
    gdB_nearest = gdB.iloc[idx].drop(columns="geometry").reset_index(drop=True)
    gdf = pd.concat(
        [
            gdA.reset_index(drop=True),
            gdB_nearest,
            pd.Series(dist, name='dist')
        ], 
        axis=1)

    return gdf


def wse_shp_to_nc(poly_wse_shp, stations_locations_txt, hdf_filename, output_nc):
    wse_gdf = gpd.read_file(poly_wse_shp)
    crs = wse_gdf.crs
    wse_gdf.to_crs(epsg=4326, inplace=True)

    station_locations = pd.read_csv(stations_locations_txt)

    hiway_gdf = gpd.GeoDataFrame(
        station_locations, geometry=gpd.points_from_xy(station_locations['x'], station_locations['y']))
    hiway_gdf.crs = 'epsg:4326'

    # Create the geometry x & y from the centroid of each feature in the wse shp file.
    wse_gdf["x"] = wse_gdf.centroid.x
    wse_gdf["y"] = wse_gdf.centroid.y

    # Create a point feature geoDataframe from the centroid geometry.
    wse_points = wse_gdf.copy()
    wse_points['geometry'] = wse_points['geometry'].centroid

    # Get a dataframe of the closest wse to each station location.
    nearest_gdf = ckdnearest(hiway_gdf, wse_points)

    # Use RAS output plan hdf to create a list of epoch times.
    hf = h5py.File(hdf_filename, 'r')
    hdf_startTime = hf['Plan Data']['Plan Information'].attrs['Simulation Start Time'].decode("utf-8")
    hdf_interval = hf['Plan Data']['Plan Information'].attrs['Base Output Interval'].decode("utf-8")
    hdf_timesteps = hf['Results']['Unsteady']['Output']['Output Blocks']['Base Output']['Unsteady Time Series']['Time']
    timesteps = hdf_timesteps.shape[0]

    startTime = datetime.strptime(hdf_startTime, '%d%b%Y %H:%M:%S')
    timesList = []

    if hdf_interval == '1HOUR':
        # Convert Interval from hours to days (1hr / 24hr / 1day)
        interval = 1/24
        t = datetime.strptime(hdf_startTime, '%d%b%Y %H:%M:%S')
        # append startTime to timeList then loop through timesteps to create full timesList
        timesList.append(t.timestamp())
        
        for step in range(timesteps-1):
            t = t + timedelta(hours = 1)
            timesList.append(t.timestamp())

    # Create the output netCDF file based on a TWI LFFS convetion for netCDF timeseries.
    ds = nc.Dataset(output_nc, "w", format="NETCDF4")
    time_dim = ds.createDimension("time", timesteps)
    station_dim = ds.createDimension("nstation", len(nearest_gdf))

    timevar = ds.createVariable("time", "i8", ("time"))
    timevar.units = "seconds since 1970-01-01 00:00:00"
    for i in range(len(timesList)):
        timevar[i] = timesList[i]

    lat = ds.createVariable("latitude", "f8", ("nstation"), zlib=True, complevel=2)
    lat.reference = "EPSG:4326"
    lon = ds.createVariable("longitude", "f8", ("nstation"), zlib=True, complevel=2)
    lon.reference = "EPSG:4326"

    datavar = ds.createVariable(
                    "values",
                    "f8",
                    ("nstation", "time"),
                    fill_value=-9999,
                    zlib=True,
                    complevel=2,
                )
    datavar.adcirc_type = "zeta"
    datavar.standard_name = "sea_surface_height_above_geoid"
    datavar.long_name = "water surface elevation above geoid"
    datavar.units = "ft"
    datavar.datum = "navd88 2009.55"
    
    point_id = ds.createVariable("point_id", "i4", ("nstation"), zlib=True, complevel=2)
    point_id[:] =nearest_gdf['point_id']

    tag = ds.createVariable("point_type", "i4", ("nstation"), zlib=True, complevel=2)
    tag.types = "0=gate, 1=levee, 2=roadway"
    tag[:]=nearest_gdf['type']

    geometry_list = np.array(list(nearest_gdf.geometry.apply(lambda x: (x.x, x.y))))
    for i in range(len(nearest_gdf)):
        lon[i] = geometry_list[i][0]
        lat[i] = geometry_list[i][1]

    ds.close()