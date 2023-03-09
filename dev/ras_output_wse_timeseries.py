import os
import shutil

import numpy as np
import h5py
import shapefile
import time
import argparse
from utils import ras_output_extract_wse_to_shp
from utils import ras_output_wse_shp_to_nc

# TODO: create the following variables as Command Line Arguments.
# --outputdir RAS post processing directory. default = "./src/system/postproc/ras/"
# postp_dir = "./src/system/postproc/ras/"

# --forecastname output directory should be based on Forecast Name: "prjName_MetModel_DurationDays_ForecastDateTime_Advisory#"
# forecastName = 'Amite_GFS_5Day_06302022_0'

# --input p#.hdf file.
# hdf_filename = r"Z:\Amite\CalValModel\Amite_20200114.p36.hdf"

# --wkt coordinate system WKT of the RAS model. use Albers by default since its used for both Amite & SLAMM RAS models.
# coord_sys = 'PROJCS["NAD_1983_BLM_Zone_15N_ftUS",GEOGCS' \
#                 '["GCS_North_American_1983",DATUM["D_North_American_1983",'\
#                 'SPHEROID["GRS_1980",6378137.0,298.257222101]],' \
#                 'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],' \
#                 'PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",1640416.666666667],' \
#                 'PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-93.0],'\
#                 'PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],'\
#                 'UNIT["Foot_US",0.3048006096012192]]'

# coord_sys = 'PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS' \
#                 '["GCS_North_American_1983",DATUM["D_North_American_1983",' \
#                 'SPHEROID["GRS_1980",6378137.0,298.257222101]],' \
#                 'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],' \
#                 'PROJECTION["Albers"],PARAMETER["false_easting",0.0],' \
#                 'PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-96.0],' \
#                 'PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],' \
#                 'PARAMETER["latitude_of_origin",23.0],UNIT["Foot_US",0.3048006096012192]]'

# --stationlocations text file deafults to: ./system/postproc/ras/station_locations.txt
# stations_locations_txt = "./src/system/postproc/ras/station_locations.txt"

# --output The output nc file with timeseries for each point in the stationlocations list. Defaults to 'RAS_PointLocations.nc'
# output = 'RAS_WSE_Timeseries.nc'

# Command Line Arguments
p = argparse.ArgumentParser(description="Point extraction for RAS time series")
p.add_argument(
    "--file", help="Name of the ras hdf plan file to extract from (Ex: 'Z:/Amite/CalValModel/Amite_20200114.p36.hdf')", 
    required=True, type=str
)
p.add_argument(
    "--points",
    help="Optional. Name of point file and location for extracted locations. Defaults to: './src/system/postproc/ras/station_locations.txt'",
    required=False,
    default = "./src/system/postproc/ras/station_locations.txt",
    type=str,
)
# --outputdir RAS post processing directory. default = "./src/system/postproc/ras/"
p.add_argument(
    "--postprocessingdirectory", help="Optional. Name of post processing directory. Defaults to: './src/system/postproc/ras/'", 
    required=False, type=str,
    default="./src/system/postproc/ras/"
)
p.add_argument(
    "--output", help="Optional. Name of output filename to create in the post processing directory. Defaults to: RAS_WSE_Timeseries.nc", 
    required=False, type=str,
    default="RAS_WSE_Timeseries.nc"
)
p.add_argument(
    "--forecastname", help="files will be added to directory based on the Forecast Name (Ex: 'prjName_MetModel_DurationDays_ForecastDateTime_Advisory#' 'Amite_GFS_5Day_06302022_0')", 
    required=True, type=str
)
p.add_argument(
    "--wkt",
    help="Optional. The WKT to descripe the coordinate system of the RAS model to assume during geoprocessing. Defaults to: Albers_Equal_Area_Conic_USGS_version",
    required=False,
    default='PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS' \
                '["GCS_North_American_1983",DATUM["D_North_American_1983",' \
                'SPHEROID["GRS_1980",6378137.0,298.257222101]],' \
                'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],' \
                'PROJECTION["Albers"],PARAMETER["false_easting",0.0],' \
                'PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-96.0],' \
                'PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],' \
                'PARAMETER["latitude_of_origin",23.0],UNIT["Foot_US",0.3048006096012192]]',
    type=str,
)

args = p.parse_args()

#-------Locations of HDF file and export directory
# hdf_filename = r"Z:\SLaMM_Mar2022\SLaMM_2021.p06.hdf"
# hdf_filename = input
curr_date = time.strftime("%Y%m%d_%H")
home_dir = os.path.join(args.postprocessingdirectory,args.forecastname)
tempDir = os.path.join(home_dir,"tempfiles")
postp_area = os.path.join(args.postprocessingdirectory, "postp_area.shp")
output_nc = os.path.join(home_dir,args.output)

if not os.path.exists(home_dir):
    os.makedirs(home_dir)

ras_output_extract_wse_to_shp.tempDirSweep(tempDir)

#***************************Open HDF5 file and get 2DArea Names, Timestep
hf = h5py.File(args.file, 'r')
list_of_2DAreas = ras_output_extract_wse_to_shp.get2DAreaNames(hf)
timesteps = ras_output_extract_wse_to_shp.get_timesteps(hf)
all_data = np.empty((0, timesteps + 6), ) # Includes an extra wse column for maximum row value

#Initialize shapefile of all 2D flow cells
poly_wse_shp = os.path.join(tempDir, 'ras_wse.shp')
w = shapefile.Writer(poly_wse_shp)

# ***********Begin writing to polygon shapefile using rows of cell index pts
# w = shapefile.Writer('test_polygon_all_data')   <-------- Currently  initialize outside of 2DArea loop
# Start writing to Shapefile
# Writing field names and types
# "C": Characters, text.
# "N": Numbers, with or without decimals.
# "F": Floats(same as "N").
# "L": Logical, for boolean True / False values.
# "D": Dates.
# "M": Memo
w.field('Area2D', 'C')
w.field('Cell_Index', 'N')
w.field('Easting', 'N', decimal=2)
w.field('Northing', 'N', decimal=2)
w.field('min_elev', 'N', decimal=2)

#Creating Results fields, same number as timesteps
i = 0
while i < timesteps:
    w.field('wse_' + str(i), 'N', decimal=2)
    i += 1
# Shapefile is closed after 2DArea loop

#Add a wse for maximum water surface at end
w.field('wse_max', 'N', decimal=2)


#Loop through all 2D flow areas in HDF file and extract geometry pts and results
for curr_2DArea in list_of_2DAreas:
    print("Current 2D Area is: %s" % curr_2DArea)

    xy_pts = np.array(ras_output_extract_wse_to_shp.get2DArea_cellcenter_pts(curr_2DArea, hf))
    min_elev = np.array(ras_output_extract_wse_to_shp.get2DCells_min_elev(curr_2DArea, hf)).round(decimals=2)
    # transpose_min_elev = min_elev.T

    wse_data = np.array(ras_output_extract_wse_to_shp.get2DArea_wse_data(curr_2DArea, hf))
    transpose_wse = wse_data.T.round(decimals=2)

    #Find WSE values that are equal to cell min elev, set to NaN, all others set to 1
    repeats_cell_min_elev = np.tile(min_elev, (timesteps,1)).T
    cell_depths = transpose_wse - repeats_cell_min_elev
    cell_depths[cell_depths > 0] = 1
    cell_depths[cell_depths == 0] = 0

    #Remove zero depth values
    filtered_transpose_wse = cell_depths * transpose_wse
    filtered_transpose_wse[filtered_transpose_wse==0] = -9999
    filtered_transpose_wse.round(decimals=2)
    max_of_row = np.max(filtered_transpose_wse, axis=1)

    cell_index = np.arange(xy_pts.shape[0])
    curr_2DArea_index = [curr_2DArea.decode('UTF-8')]* (xy_pts.shape[0])

    #Adding columns to results array
    all_data_for_curr_2DArea = np.column_stack((curr_2DArea_index,cell_index, xy_pts, min_elev))
    all_data_for_curr_2DArea = np.concatenate((all_data_for_curr_2DArea, filtered_transpose_wse), axis=1)
    all_data_for_curr_2DArea = np.column_stack((all_data_for_curr_2DArea, max_of_row))

    #Save into the overall dataset
    all_data = np.append(all_data, all_data_for_curr_2DArea, axis=0)

    # Assemble 2D Cell Polygons
    cell_face_info = ras_output_extract_wse_to_shp.get_Cells_Face_Info(hf, curr_2DArea)
    cell_face_xy_pts = ras_output_extract_wse_to_shp.get_FacePoints_Coordinates(hf, curr_2DArea)
    cell_face_index_pts = ras_output_extract_wse_to_shp.get_Cells_FacePoints_Index(hf, curr_2DArea)

    #Assemble info about perimeter faces and facepoints
    cell_facept_is_perimeter = ras_output_extract_wse_to_shp.is_FacePoint_perimeter(hf, curr_2DArea)
    face_facept_index = ras_output_extract_wse_to_shp.get_faces_FacePoint_Index(hf, curr_2DArea)
    face_perimeter_info = ras_output_extract_wse_to_shp.get_faces_Perimeter_Info(hf, curr_2DArea)
    face_perimeter_values = ras_output_extract_wse_to_shp.get_faces_Perimeter_Values(hf, curr_2DArea)
    face_orientation_info = ras_output_extract_wse_to_shp.get_face_orientation_info(hf, curr_2DArea)
    face_orientation_values = ras_output_extract_wse_to_shp.get_face_orientation_values(hf, curr_2DArea)

    #Assemble current polygons
    cell_ids = []
    index_size = len(cell_face_index_pts[0])
    curr_2DArea_Polygon_xy_pts = []
    cell_id = 0
    cell_ids = []
    for row in cell_face_index_pts:
        #find if facepoints are perimeter
        perimeter_facepts = []
        
        for facept in row:
            
            if facept != -1:

                if cell_facept_is_perimeter[facept] == -1:
                    perimeter_facepts.append(facept)
        #print(perimeter_facepts)

        #Declare empty polygon list for 2D cell
        polygon = []
        i = 0
        while i < index_size:
            curr_facept = row[i]
            
            if curr_facept != -1:
                polygon.append(cell_face_xy_pts[curr_facept])

            if i < (index_size -1) :
                next_facept = row[i+1]

            if i == (index_size -1):
                next_facept = row[0]


            if curr_facept in perimeter_facepts:
                
                if next_facept in perimeter_facepts:
                    face_index=0
                    
                    for face in face_facept_index:
                        if curr_facept == face_facept_index[face_index][0]:
                            potential_face = face_index
                            
                            if next_facept == face_facept_index[potential_face][1]:
                                next_is_first = False
                                curr_face_index = face_index
                                # print("found face")
                                break

                        if next_facept == face_facept_index[face_index][0]:
                            potential_face = face_index
                            
                            if curr_facept == face_facept_index[potential_face][1]:
                                next_is_first = True
                                curr_face_index = face_index
                                # print("found face")
                                break

                        face_index +=1


                    perimeter_st_pt = face_perimeter_info[curr_face_index][0]
                    num_perimeter_pts = face_perimeter_info[curr_face_index][1]
                    perimeter_end_pt = perimeter_st_pt + num_perimeter_pts - 1
                    perimeter_pt_index = perimeter_st_pt

                    extra_perimeter_xy_pts = []

                    # print("...adding perimeter pts, for face %s" % curr_face_index )
                    while perimeter_pt_index <= perimeter_end_pt:
                        # polygon.append(face_perimeter_values[perimeter_pt_index])
                        extra_perimeter_xy_pts.append(face_perimeter_values[perimeter_pt_index])
                        perimeter_pt_index += 1

                    if next_is_first:
                        extra_perimeter_xy_pts = extra_perimeter_xy_pts[::-1]

                    polygon.extend(extra_perimeter_xy_pts)
                #polygon.append(cell_face_xy_pts[next_facept])

            i += 1

        #Append the first face pt coordinate
        polygon.append(cell_face_xy_pts[row[0]])

        #Append to the total 2D Area set if more than 2 points (there are some lateral weirs represented like this)
        if sum(1 for n in row if n != -1)>=3:
            curr_2DArea_Polygon_xy_pts.append(polygon)
            #Keep track of cell_ids that make it into the polygon set
            cell_ids.append(cell_id)

        cell_id += 1

    #--------------Saving polygons and records to shapefile-------------
    print ("writing %s polygons to shapefile..." %curr_2DArea.decode('UTF-8'))
    str_curr_2DArea = curr_2DArea.decode('UTF-8')
    
    for row_id, poly_row in enumerate(curr_2DArea_Polygon_xy_pts):
        
        if len(poly_row) > 2:
            w.poly([poly_row[::-1]]) #clockwise flip
            #w.record(INT=nr, LOWPREC=nr, MEDPREC=nr, HIGH)
            #w.record(Area2D=str_curr_2DArea,Cell_Index=cell_id, Easting=all_data_for_curr_2DArea[cell_id][])
            #w.record('Area2D', str_curr_2DArea)
            #w.record('Cell Index', cell_id)
            records = np.array(all_data_for_curr_2DArea[cell_ids[row_id]]).tolist()
            w.record(*records)

#Close 2DArea polygon shapefile
print("Closing shapefile with all 2D Area polygons.")
w.close()

print("Writing Projection file w/ hardcoded coordinate system.")
with open(os.path.join(tempDir,'ras_wse.prj'), 'w') as f:
    f.write(args.wkt)
    f.close()

ras_output_wse_shp_to_nc.wse_shp_to_nc(poly_wse_shp, args.points, args.file, output_nc)