import geopandas as gpd
from   shapely.geometry import Polygon
import shapefile # pyshed library
import shapely
import sys
import numpy as np
import warnings
import pandas as pd
import os

def check_easymore_input(
    temp_dir,
    output_dir,
    author_name,
    var_names,
    format_list,
    fill_value_list,
    remap_csv,
    var_names_remapped,
):
    if temp_dir != '':
        if temp_dir[-1] != '/':
            sys.exit('the provided temporary folder for EASYMORE should end with (/)')
        if not os.path.isdir(temp_dir):
            os.mkdir(temp_dir)
    if output_dir == '':
        sys.exit('the provided folder for EASYMORE remapped netCDF output is missing; please provide that')
    if output_dir != '':
        if output_dir[-1] != '/':
            sys.exit('the provided output folder for EASYMORE should end with (/)')
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
    if temp_dir == '':
        print("No temporary folder is provided for EASYMORE; this will result in EASYMORE saving the files in the same directory as python script")
    if author_name == '':
        print("no author name is provide and the author name is changed to (author name)!")
        author_name = "author name"
    if (len(var_names) != 1) and (len(format_list) == 1) and (len(fill_value_list) ==1):
        if (len(var_names) != len(fill_value_list)) and \
        (len(var_names) != len(format_list)) and \
        (len(format_list) == 1) and (len(fill_value_list) ==1):
            print('EASYMORE is given multiple varibales to be remapped but only on format and fill value'+\
                'EASYMORE repeat the format and fill value for all the variables in output files')
            format_list     = format_list     * len(var_names)
            fill_value_list = fill_value_list * len(var_names)
        else:
            sys.exit('number of varibales and fill values and formats do not match')
    if remap_csv != '':
        print('remap file is provided; EASYMORE will use this file and skip calculation of remapping')
    if len(var_names) != len(set(var_names)):
        sys.exit('the name of the variables you have provided from the source NetCDF file to be remapped are not unique')
    if var_names_remapped:
        if len(var_names_remapped) != len(set(var_names_remapped)):
            sys.exit('the name of the variables you have provided as the rename in the remapped file are not unique')
        if len(var_names_remapped) != len(var_names):
            sys.exit('the number of provided variables from the source file and names to be remapped are not the same length')
    if not var_names_remapped:
        var_names_remapped = var_names
    for i in np.arange(len(var_names)):
        print('EASYMORE will remap variable ',var_names[i],' from source file to variable ',var_names_remapped[i],' in remapped NeCDF file')
        
    return author_name, format_list, fill_value_list, var_names_remapped

def check_target_shp(shp, shp_ID, lat, lon, sort_ID):
        """
        @ author:                  Shervan Gharari
        @ Github:                  https://github.com/ShervanGharari/EASYMORE
        @ author's email id:       sh.gharari@gmail.com
        @ license:                 GNU-GPLv3
        this function check if target shapefile and add ID and centroid lat and lon is not provided
        Arguments
        ---------
        shp: geopandas dataframe, polygone, multipolygon, point, multipoint
        """
        # load the needed packages

        # sink/target shapefile check the projection
        if 'epsg:4326' not in str(shp.crs).lower():
            sys.exit('please project your shapefile to WGS84 (epsg:4326)')
        else: # check if the projection is WGS84 (or epsg:4326)
            print('EASYMORE detects that target shapefile is in WGS84 (epsg:4326)')
        # check if the ID, latitude, longitude are provided
        if shp_ID == '':
            print('EASYMORE detects that no field for ID is provided in sink/target shapefile')
            print('arbitarary values of ID are added in the field ID_t')
            shp['ID_t']  = np.arange(len(shp))+1
        else:
            print('EASYMORE detects that the field for ID is provided in sink/target shapefile')
            # check if the provided IDs are unique
            ID_values = np.array(shp[shp_ID])
            if len(ID_values) != len(np.unique(ID_values)):
                sys.exit('The provided IDs in shapefile are not unique; provide unique IDs or do not identify target_shp_ID')
            shp['ID_t'] = shp[shp_ID]
        if lat == '' or lon == '':
            print('EASYMORE detects that either of the fields for latitude or longitude is not provided in sink/target shapefile')
            # in WGS84
            print('calculating centroid of shapes in WGS84 projection;')
            print('for better appximation use the easymore equal area centroid function to preprocess target shapefile')
            df_point = pd.DataFrame()
            warnings.simplefilter('ignore') # silent the warning
            df_point ['lat'] = shp.centroid.y
            df_point ['lon'] = shp.centroid.x
            warnings.simplefilter('default') # back to normal
        if lat == '':
            print('EASYMORE detects that no field for latitude is provided in sink/target shapefile')
            print('latitude values are added in the field lat_t')
            shp['lat_t']  = df_point ['lat'] # centroid lat from target
        else:
            print('EASYMORE detects that the field latitude is provided in sink/target shapefile')
            shp['lat_t'] = shp[lat]
        if lon == '':
            print('EASYMORE detects that no field for longitude is provided in sink/target shapefile')
            print('longitude values are added in the field lon_t')
            shp['lon_t']  = df_point ['lon'] # centroid lon from target
        else:
            print('EASYMORE detects that the field longitude is provided in sink/target shapefile')
            shp['lon_t'] = shp[lon]
        # check other geometries and add buffer if needed
        detected_points = False
        detected_multipoints = False
        detected_lines = False
        for index, _ in shp.iterrows():
            polys = shp.geometry.iloc[index] # get the shape
            if polys.geom_type.lower() == "point"      or polys.geom_type.lower() == "points":
                detected_points = True
                polys = polys.buffer(10**-5).simplify(10**-5)
                shp.geometry.iloc[index] = polys
            if polys.geom_type.lower() == "multipoint" or polys.geom_type.lower() == "multipoints":
                detected_multipoints = True
                polys = polys.buffer(10**-5).simplify(10**-5)
                shp.geometry.iloc[index] = polys
            if polys.geom_type.lower() == "line"       or polys.geom_type.lower() == "lines":
                detected_lines = True
                polys = polys.buffer(10**-5).simplify(10**-5)
                shp.geometry.iloc[index] = polys
            if polys.geom_type.lower() == "polyline"   or polys.geom_type.lower() == "polylines":
                detected_lines = True
                polys = polys.buffer(10**-5).simplify(10**-5)
                shp.geometry.iloc[index] = polys
        # print messages
        if detected_points:
            print('EASYMORE detects point(s) as geometry of target shapefile and will apply small buffer to them')
        if detected_multipoints:
            print('EASYMORE detected multipoint as geometry of target shapefile and will considere it as multipolygone')
            print('hence EASYMORE will provide the average of all the point in each multipoint')
            print('if you mistakenly have given poitns as multipoints please correct the target shapefile')
        if detected_lines:
            print('EASYMORE detected line as geometry of target shapefile and will considere it as polygon (adding small buffer)')
        print('it seems everything is OK with the sink/target shapefile; added to EASYMORE object target_shp_gpd')
        if sort_ID:
            shp = shp.sort_values(by='ID_t')
            shp = shp.reset_index(drop=True)
        shp['order'] = np.arange(len(shp)) + 1 # order of the shapefile
        return shp

def check_source_nc(
    source_nc,
    var_names,
    var_lat,
    var_lon,
    var_time,
    tolerance,
 ):
    import glob
    import netCDF4 as nc4
    flag_do_not_match = False
    nc_names = glob.glob (source_nc)
    if not nc_names:
        sys.exit('EASYMORE detects no netCDF file; check the path to the soure netCDF files')
    else:
        ncid      = nc4.Dataset(nc_names[0])
        var_dim   = list(ncid.variables[var_names[0]].dimensions)
        lat_dim   = list(ncid.variables[var_lat].dimensions)
        lon_dim   = list(ncid.variables[var_lon].dimensions)
        lat_value = np.array(ncid.variables[var_lat])
        lon_value = np.array(ncid.variables[var_lon])
        # dimension check based on the first netcdf file
        if not (set(lat_dim) <= set(var_dim)):
            flag_do_not_match = True
        if not (set(lon_dim) <= set(var_dim)):
            flag_do_not_match = True
        if (len(lat_dim) == 2) and (len(lon_dim) == 2) and (len(var_dim) == 3): # case 2
            if not (set(lat_dim) == set(lon_dim)):
                flag_do_not_match = True
        if (len(lat_dim) == 1) and (len(lon_dim) == 1) and (len(var_dim) == 2): # case 3
            if not (set(lat_dim) == set(lon_dim)):
                flag_do_not_match = True
        # dimension check and consistancy for variable latitude
        for nc_name in nc_names:
            ncid = nc4.Dataset(nc_name)
            temp = list(ncid.variables[var_lat].dimensions)
            # fist check the length of the temp and lat_dim
            if len(temp) != len(lat_dim):
                flag_do_not_match = True
            else:
                for i in np.arange(len(temp)):
                    if temp[i] != lat_dim[i]:
                        flag_do_not_match = True
            temp = np.array(ncid.variables[var_lat])
            if np.sum(abs(lat_value-temp))>tolerance:
                flag_do_not_match = True
        # dimension check and consistancy for variable longitude
        for nc_name in nc_names:
            ncid = nc4.Dataset(nc_name)
            temp = list(ncid.variables[var_lon].dimensions)
            # fist check the length of the temp and lon_dim
            if len(temp) != len(lon_dim):
                flag_do_not_match = True
            else:
                for i in np.arange(len(temp)):
                    if temp[i] != lon_dim[i]:
                        flag_do_not_match = True
            temp = np.array(ncid.variables[var_lon])
            if np.sum(abs(lon_value-temp))>tolerance:
                flag_do_not_match = True
        # dimension check consistancy for variables to be remapped
        for var_name in var_names:
            # get the varibale information of lat, lon and dimensions of the varibale.
            for nc_name in nc_names:
                ncid = nc4.Dataset(nc_name)
                temp = list(ncid.variables[var_name].dimensions)
                # fist check the length of the temp and var_dim
                if len(temp) != len(var_dim):
                    flag_do_not_match = True
                else:
                    for i in np.arange(len(temp)):
                        if temp[i] != var_dim[i]:
                            flag_do_not_match = True
        # check varibale time and dimension time are the same name so time is coordinate
        for nc_name in nc_names:
            ncid = nc4.Dataset(nc_name)
            temp = ncid.variables[var_time].dimensions
            if len(temp) != 1:
                sys.exit('EASYMORE expects 1D time varibale, it seems time varibales has more than 1 dimension')
            if str(temp[0]) != var_time:
                sys.exit('EASYMORE expects time varibale and dimension to be different, they should be the same\
                for xarray to consider time dimension as coordinates')
    if flag_do_not_match:
        sys.exit('EASYMORE detects that all the provided netCDF files and varibale \
has different dimensions for the varibales or latitude and longitude')
    else:
        print('EASYMORE detects that the varibales from the netCDF files are identical\
in dimensions of the varibales and latitude and longitude')
        print('EASYMORE detects that all the varibales have dimensions of:')
        print(var_dim)
        print('EASYMORE detects that the longitude varibales has dimensions of:')
        print(lon_dim)
        print('EASYMORE detects that the latitude varibales has dimensions of:')
        print(lat_dim)

def check_source_nc_shp(source_nc, source_shp,
                        source_shp_lat, source_shp_lon,
                        lat, lon):
    """
    @ author:                  Shervan Gharari
    @ Github:                  https://github.com/ShervanGharari/EASYMORE
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 GNU-GPLv3
    This function checks the source netcdf file shapefile
    needs more development
    """
    # load the needed packages
    import geopandas as gpd
    from   shapely.geometry import Polygon
    import netCDF4 as nc4
    import glob
    #
    multi_source = False
    nc_names = glob.glob (source_nc)
    ncid = nc4.Dataset(nc_names[0])
    # sink/target shapefile is what we want the varibales to be remapped to
    shp = gpd.read_file(source_shp)
    if 'epsg:4326' not in str(shp.crs).lower():
        sys.exit('please project your source shapefile and varibales in source nc files to WGS84 (epsg:4326)')
    else: # check if the projection is WGS84 (or epsg:4326)
        print('EASYMORE detects that source shapefile is in WGS84 (epsg:4326)')
    # get the lat/lon from source shapfile and nc files
    lat_shp = np.array(shp[source_shp_lat]); lat_shp = lat_shp.astype(float)
    lon_shp = np.array(shp[source_shp_lon]); lon_shp = lon_shp.astype(float)
    #
    coord_shp         = pd.DataFrame()
    coord_shp ['lon'] = lon_shp
    coord_shp ['lat'] = lat_shp
    coord_shp         = coord_shp.sort_values(by='lon')
    coord_shp         = coord_shp.sort_values(by='lat')
    # check if all the lat/lon in shapefile are unique
    df_temp           = coord_shp.drop_duplicates()
    if len(df_temp)  != len(coord_shp):
        # print
        sys.exit('The latitude and longitude in source shapefile are not unique')
    # first check the model
    coord_nc          = pd.DataFrame()
    coord_nc  ['lon'] = lon
    coord_nc  ['lat'] = lat
    coord_nc          = coord_nc.sort_values(by='lon')
    coord_nc          = coord_nc.sort_values(by='lat')
    # check if all the lat/lon in shapefile are unique
    df_temp           = coord_nc.drop_duplicates()
    if len(df_temp)  != len(coord_nc):
        # print
        #sys.exit('The latitude and longitude in source NetCDF files are not unique')
        print('The latitude and longitude in source NetCDF files are not unique')
