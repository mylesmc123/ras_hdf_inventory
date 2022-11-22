import h5py, sys
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, LineString, Point
import folium

with h5py.File(r"Z:\OneDrive_1_11-9-2022\HydraulicsBuild_MB.p02.hdf",'r') as hf:
    
    df = pd.DataFrame()
    df['XS_Name'] = hf['/Results/Unsteady/Geometry Info/Cross Section Only'][:]
    
    gdf = gpd.GeoDataFrame()

    m = folium.Map([30.2, -90.3], zoom_start=8, tiles='cartodbdark_matter')
    
    coord_index_start = 0
    
    for i,v in enumerate(df['XS_Name']):
        # pick style based on even or odd index
        if (i % 2) == 0:
            style = {"color": "#228B22", "weight": 5, "opacity": 0.65}
        else:
            # style = {'fillColor': '#00FFFFFF', 'color': '#00FFFFFF'}
            style = {"color": "#ff7800", "weight": 5, "opacity": 0.65}
        parts = hf['/Geometry/Cross Sections/Polyline Parts'][i,1]
        # the index range of coordinates to pull based on the number of parts for each XS.
        coord_index_end = coord_index_start + parts
        coords_range = [coord_index_start, coord_index_end]
        bounds = [hf['/Geometry/Cross Sections/Polyline Points'][coord_index_start:coord_index_end].tolist()][0]
        # update start of range for next XS in loop.
        coord_index_start += parts
        line_geom = LineString(bounds)
        line = gpd.GeoDataFrame(index=[i], geometry=[line_geom])
        line['name'] = [v.decode("utf-8") ]
        # print(str(v))
        line['style'] = [style]
        line.crs = 'PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-96.0],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_origin",23.0],UNIT["Foot_US",0.3048006096012192]]'
        gdf = gdf.append(line)

g = folium.GeoJson(data=gdf).add_to(m)
folium.GeoJsonTooltip(fields=["name"]).add_to(g)
m.save(r'./out/hdf_XS_map.html')
gdf.to_file('./out/hdf_XS.shp')  