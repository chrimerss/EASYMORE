"""
A parallel version of easymore
"""

import glob
import time
import netCDF4      as nc4
import numpy        as np
import pandas       as pd
import xarray       as xr
import geopandas as gpd
import dask
import sys
import os
import warnings
# warnings.filters()
from   datetime     import datetime
sys.path.append('/home/ZhiLi/US_routing/EASYMORE/easymore')
from check import check_target_shp, check_easymore_input, check_source_nc,\
                   check_source_nc_shp

def NetCDF_SHP_lat_lon(
    source_nc: str,
    var_names: list,
    var_ID: str,
    var_lon: str,
    var_lat: str,
    tolerance: float
) -> tuple:
    """
    @ author:                  Shervan Gharari
    @ Github:                  https://github.com/ShervanGharari/EASYMORE
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 GNU-GPLv3
    This function checks dimension of the source shapefile and checks the case of regular, rotated, and irregular
    also created the 2D array of lat and lon for creating the shapefile
    """
    from shapely.geometry import Polygon
    import shapefile # pyshed library
    #
    nc_names = glob.glob (source_nc)
    var_name = var_names[0]
    # open the nc file to read
    ncid = nc4.Dataset(nc_names[0])
    # deciding which case
    # case #1 regular latitude/longitude
    if (len(ncid.variables[var_lon].dimensions)==1) and\
    (len(ncid.variables[var_lon].dimensions)==1) and\
    (len(ncid.variables[var_names[0]].dimensions)==3):
        print('EASYMORE detects case 1 - regular lat/lon')
        case = 1
        # get the list of dimensions for the ncid sample varibale
        list_dim_name = list(ncid.variables[var_names[0]].dimensions)
        # get the location of lat dimensions
        location_of_lat = list_dim_name.index(list(ncid.variables[var_lat].dimensions)[0])
        locaiton_of_lon = list_dim_name.index(list(ncid.variables[var_lon].dimensions)[0])
        # det the dimensions of lat and lon
        len_of_lat = len(ncid.variables[var_lat][:])
        len_of_lon = len(ncid.variables[var_lon][:])
        if locaiton_of_lon > location_of_lat:
            lat = np.zeros([len_of_lat, len_of_lon])
            lon = np.zeros([len_of_lat, len_of_lon])
            for i in np.arange(len(ncid.variables[var_lon][:])):
                lat [:,i] = ncid.variables[var_lat][:]
            for i in np.arange(len(ncid.variables[var_lat][:])):
                lon [i,:] = ncid.variables[var_lon][:]
        else:
            lat = np.zeros([len_of_lon, len_of_lat])
            lon = np.zeros([len_of_lon, len_of_lat])
            for i in np.arange(len(ncid.variables[var_lon][:])):
                lat [i,:] = ncid.variables[var_lat][:]
            for i in np.arange(len(ncid.variables[var_lat][:])):
                lon [:,i] = ncid.variables[var_lon][:]
        # check if lat and lon are spaced equally
        lat_temp = np.array(ncid.variables[var_lat][:])
        lat_temp_diff = np.diff(lat_temp)
        lat_temp_diff_2 = np.diff(lat_temp_diff)
        max_lat_temp_diff_2 = max(abs(lat_temp_diff_2))
        print('max difference of lat values in source nc files are : ', max_lat_temp_diff_2)
        lon_temp = np.array(ncid.variables[var_lon][:])
        lon_temp_diff = np.diff(lon_temp)
        lon_temp_diff_2 = np.diff(lon_temp_diff)
        max_lon_temp_diff_2 = max(abs(lon_temp_diff_2))
        print('max difference of lon values in source nc files are : ', max_lon_temp_diff_2)
        # save lat, lon into the object
        lat      = np.array(lat).astype(float)
        lon      = np.array(lon).astype(float)
        ID= None
        # expanding just for the the creation of shapefile with first last rows and columns
        # if (max_lat_temp_diff_2<tolerance) and (max_lon_temp_diff_2<tolerance): # then lat lon are spaced equal
        # create expanded lat
        lat_expanded = np.zeros(np.array(lat.shape)+2)
        lat_expanded [1:-1,1:-1] = lat
        lat_expanded [:, 0]  = lat_expanded [:, 1] + (lat_expanded [:, 1] - lat_expanded [:, 2]) # populate left column
        lat_expanded [:,-1]  = lat_expanded [:,-2] + (lat_expanded [:,-2] - lat_expanded [:,-3]) # populate right column
        lat_expanded [0, :]  = lat_expanded [1, :] + (lat_expanded [1, :] - lat_expanded [2, :]) # populate top row
        lat_expanded [-1,:]  = lat_expanded [-2,:] + (lat_expanded [-2,:] - lat_expanded [-3,:]) # populate bottom row
        # create expanded lat
        lon_expanded = np.zeros(np.array(lon.shape)+2)
        lon_expanded [1:-1,1:-1] = lon
        lon_expanded [:, 0]  = lon_expanded [:, 1] + (lon_expanded [:, 1] - lon_expanded [:, 2]) # populate left column
        lon_expanded [:,-1]  = lon_expanded [:,-2] + (lon_expanded [:,-2] - lon_expanded [:,-3]) # populate right column
        lon_expanded [0, :]  = lon_expanded [1, :] + (lon_expanded [1, :] - lon_expanded [2, :]) # populate top row
        lon_expanded [-1,:]  = lon_expanded [-2,:] + (lon_expanded [-2,:] - lon_expanded [-3,:]) # populate bottom row
    # case #2 rotated lat/lon
    if (len(ncid.variables[var_lat].dimensions)==2) and (len(ncid.variables[var_lon].dimensions)==2):
        print('EASYMORE detects case 2 - rotated lat/lon')
        case = 2
        lat = ncid.variables[var_lat][:,:]
        lon = ncid.variables[var_lon][:,:]
        # creating/saving the shapefile
        lat = np.array(lat).astype(float)
        lon = np.array(lon).astype(float)
        ID= None
        lat_expanded= None
        lon_expanded= None
    # case #3 1-D lat/lon and 2 data for irregulat shapes
    if (len(ncid.variables[var_lat].dimensions)==1) and (len(ncid.variables[var_lon].dimensions)==1) and\
        (len(ncid.variables[var_names[0]].dimensions)==2):
        print('EASYMORE detects case 3 - irregular lat/lon; shapefile should be provided')
        case = 3
        lat = ncid.variables[var_lat][:]
        lon = ncid.variables[var_lon][:]
        #print(lat, lon)
        if var_ID  == '':
            print('EASYMORE detects that no varibale for ID of the source netCDF file; an arbitatiry ID will be provided')
            ID =  np.arange(len(lat))+1 # pass arbitarary values
        else:
            ID = ncid.variables[var_ID][:]
        # creating/saving the shapefile
        lat = np.array(lat).astype(float)
        lon = np.array(lon).astype(float)
        lon_expanded= None
        lat_expanded= None

    return (case, ID, lat, lon, lat_expanded, lon_expanded)


def lat_lon_SHP(lat, lon):
    """
    @ author:                  Shervan Gharari, Wouter Knoben
    @ Github:                  https://github.com/ShervanGharari/EASYMORE
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 GNU-GPLv3
    This function creates a shapefile for the source netcdf file
    Arguments
    ---------
    lat: the 2D matrix of lat_2D [n,m,]
    lon: the 2D matrix of lon_2D [n,m,]
    file_name: string, name of the file that the shapefile will be saved at
    """
    # check if lat/lon that are taken in has the same dimension
    from   shapely.geometry import Polygon
    from numba import jit
    # import shapefile # pyshed library
    # import shapely
    #
    lat_lon_shape = lat.shape
    # write the shapefile
    @jit(nopython=True)
    def bottleneck_geometry(lat: np.ndarray, lon: np.ndarray) -> tuple:
        '''
        replace two for-loop with C++
        '''
        lat_lon_shape= lon.shape
        geometry = []
        records= []
        m=0
        for i in range(1, lat_lon_shape[0] - 1):
            for j in range(1, lat_lon_shape[1] - 1):
                # checking is lat and lon is located inside the provided bo
                # empty the polygon variable
                # update records
                m += 1 # ID
                center_lat = lat[i,j] # lat value of data point in source .nc file
                center_lon = lon[i,j] # lon value of data point in source .nc file should be within [0,360]
                # Creating the lat of the shapefile
                Lat_Up       = (lat[i - 1, j] + lat[i, j]) / 2
                Lat_UpRright = (lat[i - 1, j] + lat[i - 1, j + 1] + lat[i, j + 1] + lat[i, j]) / 4
                Lat_Right    = (lat[i, j + 1] + lat[i, j]) / 2
                Lat_LowRight = (lat[i, j + 1] + lat[i + 1, j + 1] + lat[i + 1, j] + lat[i, j]) / 4
                Lat_Low      = (lat[i + 1, j] + lat[i, j]) / 2
                Lat_LowLeft  = (lat[i, j - 1] + lat[i + 1, j - 1] + lat[i + 1, j] + lat[i, j]) / 4
                Lat_Left     = (lat[i, j - 1] + lat[i, j]) / 2
                Lat_UpLeft   = (lat[i - 1, j - 1] + lat[i - 1, j] + lat[i, j - 1] + lat[i, j]) / 4
                # Creating the lon of the shapefile
                Lon_Up       = (lon[i - 1, j] + lon[i, j]) / 2
                Lon_UpRright = (lon[i - 1, j] + lon[i - 1, j + 1] + lon[i, j + 1] + lon[i, j]) / 4
                Lon_Right    = (lon[i, j + 1] + lon[i, j]) / 2
                Lon_LowRight = (lon[i, j + 1] + lon[i + 1, j + 1] + lon[i + 1, j] + lon[i, j]) / 4
                Lon_Low      = (lon[i + 1, j] + lon[i, j]) / 2
                Lon_LowLeft  = (lon[i, j - 1] + lon[i + 1, j - 1] + lon[i + 1, j] + lon[i, j]) / 4
                Lon_Left     = (lon[i, j - 1] + lon[i, j]) / 2
                Lon_UpLeft   = (lon[i - 1, j - 1] + lon[i - 1, j] + lon[i, j - 1] + lon[i, j]) / 4
                geometry.append([ (Lon_Up,        Lat_Up),\
                                (Lon_UpRright,  Lat_UpRright), \
                                (Lon_Right,     Lat_Right), \
                                (Lon_LowRight,  Lat_LowRight), \
                                (Lon_Low,       Lat_Low), \
                                (Lon_LowLeft,   Lat_LowLeft), \
                                (Lon_Left,      Lat_Left), \
                                (Lon_UpLeft,    Lat_UpLeft), \
                                (Lon_Up,        Lat_Up)])
                records.append([m, center_lat, center_lon])

        return (geometry, records)

    geometry, records= bottleneck_geometry(lat, lon)
    records= np.stack(records)
    gdf= gpd.GeoDataFrame(geometry=[Polygon(_geom) for _geom in geometry])
    gdf['ID_s']= records[:,0].astype(np.int32)
    gdf['lat_s']= records[:,1].astype(np.float32)
    gdf['lon_s']= records[:,2].astype(np.float32)

    return gdf

def add_lat_lon_source_SHP(
                            shp,
                            source_shp_lat,
                            source_shp_lon,
                            source_shp_ID):
    """
    @ author:                  Shervan Gharari
    @ Github:                  https://github.com/ShervanGharari/EASYMORE
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 GNU-GPLv3
    This function add lat, lon and ID from the source shapefile if provided
    Arguments
    ---------
    shp: geodataframe of source shapefile
    source_shp_lat: string, the name of the lat field in the source shapefile
    source_shp_lon: string, the name of the lon field in the source shapefile
    source_shp_ID: string, the name of the ID field in the source shapefile
    """
    import geopandas as gpd
    from   shapely.geometry import Polygon
    import shapefile # pyshed library
    import shapely
    shp['lat_s'] = shp [source_shp_lat].astype(float)
    shp['lon_s'] = shp [source_shp_lon].astype(float)
    if source_shp_ID != '':
        shp ['ID_s']  = shp [source_shp_ID]
    else:
        shp ['ID_s']  = np.arange(len(shp))+1
    return shp

def intersection_shp(
                    shp_1,
                    shp_2):
    from   shapely.geometry import Polygon
    from numba import jit
    # import shapefile # pyshed library
    # import shapely
    """
    @ author:                  Shervan Gharari
    @ Github:                  https://github.com/ShervanGharari/EASYMORE
    @ author's email id:       sh.gharari@gmail.com
    @license:                  GNU-GPLv3
    This fucntion intersect two shapefile. It keeps the fiels from the first and second shapefiles (identified by S_1_ and
    S_2_). It also creats other field including AS1 (area of the shape element from shapefile 1), IDS1 (an arbitary index
    for the shapefile 1), AS2 (area of the shape element from shapefile 1), IDS2 (an arbitary index for the shapefile 1),
    AINT (the area of teh intersected shapes), AP1 (the area of the intersected shape to the shapes from shapefile 1),
    AP2 (the area of teh intersected shape to the shapefes from shapefile 2), AP1N (the area normalized in the case AP1
    summation is not 1 for a given shape from shapefile 1, this will help to preseve mass if part of the shapefile are not
    intersected), AP2N (the area normalized in the case AP2 summation is not 1 for a given shape from shapefile 2, this
    will help to preseve mass if part of the shapefile are not intersected)
    Arguments
    ---------
    shp_1: geo data frame, shapefile 1
    shp_2: geo data frame, shapefile 2
    Returns
    -------
    result: a geo data frame that includes the intersected shapefile and area, percent and normalized percent of each shape
    elements in another one
    """
    @jit(nopython=True)
    def normalize(ids: np.ndarray, unique_ids: np.ndarray,
                        ap: np.ndarray) -> np.ndarray:
        ap_new= ap.copy()
        for i in unique_ids:
            idx= np.where(ids==i)
            ap_new[idx]= ap[idx]/ap[idx].sum()
        return ap_new

    # get the column name of shp_1
    column_names = shp_1.columns
    column_names = list(column_names)
    # removing the geometry from the column names
    column_names.remove('geometry')
    # renaming the column with S_1
    for i in range(len(column_names)):
        shp_1 = shp_1.rename(
            columns={column_names[i]: 'S_1_' + column_names[i]})
    # Caclulating the area for shp1
    shp_1['AS1']  = shp_1.area
    shp_1['IDS1'] = np.arange(shp_1.shape[0])+1
    # get the column name of shp_2
    column_names = shp_2.columns
    column_names = list(column_names)
    # removing the geometry from the colomn names
    column_names.remove('geometry')
    # renaming the column with S_2
    for i in range(len(column_names)):
        shp_2 = shp_2.rename(
            columns={column_names[i]: 'S_2_' + column_names[i]})
    # Caclulating the area for shp2
    shp_2['AS2'] = shp_2.area
    shp_2['IDS2'] = np.arange(shp_2.shape[0])+1
    # making intesection
    result = spatial_overlays(shp_1, shp_2, how='intersection')
    if result.empty:
        pass
    else:
        # Caclulating the area for shp2
        result['AINT'] = result['geometry'].area
        result['AP1'] = result['AINT']/result['AS1']
        result['AP2'] = result['AINT']/result['AS2']
        # taking the part of data frame as the numpy to incread the spead
        # finding the IDs from shapefile one
        ID_S1 = np.array (result['IDS1'])
        AP1 = np.array(result['AP1'])
        AP1N = AP1 # creating the nnormalized percent area
        ID_S1_unique = np.unique(ID_S1) #unique idea
        AP1N=normalize(ID_S1, ID_S1_unique, AP1)
        #replace computationally expensive part with C++
        # for i in ID_S1_unique:
        #     INDX = np.where(ID_S1==i) # getting the indeces
        #     AP1N[INDX] = AP1[INDX] / AP1[INDX].sum() # normalizing for that sum
        # taking the part of data frame as the numpy to incread the spead
        # finding the IDs from shapefile one
        ID_S2 = np.array (result['IDS2'])
        AP2 = np.array(result['AP2'])
        AP2N = AP2 # creating the nnormalized percent area
        ID_S2_unique = np.unique(ID_S2) #unique idea
        AP2N= normalize(ID_S2, ID_S2_unique, AP2)
        # replace computationally expensive part with C++
        # for i in ID_S2_unique:
        #     INDX = np.where(ID_S2==i) # getting the indeces
        #     AP2N[INDX] = AP2[INDX] / AP2[INDX].sum() # normalizing for that sum
        result ['AP1N'] = AP1N
        result ['AP2N'] = AP2N
    return result

def spatial_overlays(
                    df1,
                    df2,
                    how='intersection',
                    reproject=True):
    from   shapely.geometry import Polygon
    """
    Perform spatial overlay between two polygons.
    Currently only supports data GeoDataFrames with polygons.
    Implements several methods that are all effectively subsets of
    the union.
    author: Omer Ozak
    https://github.com/ozak
    https://github.com/geopandas/geopandas/pull/338
    license: GNU-GPLv3
    Parameters
    ----------
    df1: GeoDataFrame with MultiPolygon or Polygon geometry column
    df2: GeoDataFrame with MultiPolygon or Polygon geometry column
    how: string
        Method of spatial overlay: 'intersection', 'union',
        'identity', 'symmetric_difference' or 'difference'.
    use_sindex : boolean, default True
        Use the spatial index to speed up operation if available.
    Returns
    -------
    df: GeoDataFrame
        GeoDataFrame with new set of polygons and attributes
        resulting from the overlay
    """

    df1 = df1.copy()
    df2 = df2.copy()
    df1['geometry'] = df1.geometry.buffer(0)
    df2['geometry'] = df2.geometry.buffer(0)
    if df1.crs!=df2.crs and reproject:
        print('Data has different projections.')
        print('Converted data to projection of first GeoPandas DatFrame')
        df2.to_crs(crs=df1.crs, inplace=True)
    if how=='intersection':
        # Spatial Index to create intersections
        spatial_index = df2.sindex
        df1['bbox'] = df1.geometry.apply(lambda x: x.bounds)
        df1['sidx']=df1.bbox.apply(lambda x:list(spatial_index.intersection(x)))
        pairs = df1['sidx'].to_dict()
        nei = []
        for i,j in pairs.items():
            for k in j:
                nei.append([i,k])
        pairs = gpd.GeoDataFrame(nei, columns=['idx1','idx2'])
        pairs = pairs.merge(df1, left_on='idx1', right_index=True)
        pairs = pairs.merge(df2, left_on='idx2', right_index=True, suffixes=['_1','_2'])
        if pairs.empty:
            return pd.DataFrame()
        else:
            pairs['Intersection'] = pairs.apply(lambda x: (x['geometry_1'].intersection(x['geometry_2'])).buffer(0), axis=1)
            #pairs = gpd.GeoDataFrame(pairs, columns=pairs.columns, crs=df1.crs)
            pairs = gpd.GeoDataFrame(pairs, columns=pairs.columns)
            cols = pairs.columns.tolist()
            cols.remove('geometry_1')
            cols.remove('geometry_2')
            cols.remove('sidx')
            cols.remove('bbox')
            cols.remove('Intersection')
            dfinter = pairs[cols+['Intersection']].copy()
            dfinter.rename(columns={'Intersection':'geometry'}, inplace=True)
            #dfinter = gpd.GeoDataFrame(dfinter, columns=dfinter.columns, crs=pairs.crs)
            dfinter = gpd.GeoDataFrame(dfinter, columns=dfinter.columns, crs=df1.crs)
            dfinter = dfinter.loc[dfinter.geometry.is_empty==False]
            dfinter.drop(['idx1','idx2'], inplace=True, axis=1)
            return dfinter
    elif how=='difference':
        spatial_index = df2.sindex
        df1['bbox'] = df1.geometry.apply(lambda x: x.bounds)
        df1['sidx'] = df1.bbox.apply(lambda x:list(spatial_index.intersection(x)))
        df1['new_g'] = df1.apply(lambda x: reduce(lambda x, y: x.difference(y).buffer(0),
                                    [x.geometry]+list(df2.iloc[x.sidx].geometry)) , axis=1)
        df1.geometry = df1.new_g
        df1 = df1.loc[df1.geometry.is_empty==False].copy()
        df1.drop(['bbox', 'sidx', 'new_g'], axis=1, inplace=True)
        return df1
    elif how=='symmetric_difference':
        df1['idx1'] = df1.index.tolist()
        df2['idx2'] = df2.index.tolist()
        df1['idx2'] = np.nan
        df2['idx1'] = np.nan
        dfsym = df1.merge(df2, on=['idx1','idx2'], how='outer', suffixes=['_1','_2'])
        dfsym['geometry'] = dfsym.geometry_1
        dfsym.loc[dfsym.geometry_2.isnull()==False, 'geometry'] = dfsym.loc[dfsym.geometry_2.isnull()==False, 'geometry_2']
        dfsym.drop(['geometry_1', 'geometry_2'], axis=1, inplace=True)
        dfsym = gpd.GeoDataFrame(dfsym, columns=dfsym.columns, crs=df1.crs)
        spatial_index = dfsym.sindex
        dfsym['bbox'] = dfsym.geometry.apply(lambda x: x.bounds)
        dfsym['sidx'] = dfsym.bbox.apply(lambda x:list(spatial_index.intersection(x)))
        dfsym['idx'] = dfsym.index.values
        dfsym.apply(lambda x: x.sidx.remove(x.idx), axis=1)
        dfsym['new_g'] = dfsym.apply(lambda x: reduce(lambda x, y: x.difference(y).buffer(0),
                            [x.geometry]+list(dfsym.iloc[x.sidx].geometry)) , axis=1)
        dfsym.geometry = dfsym.new_g
        dfsym = dfsym.loc[dfsym.geometry.is_empty==False].copy()
        dfsym.drop(['bbox', 'sidx', 'idx', 'idx1','idx2', 'new_g'], axis=1, inplace=True)
        return dfsym
    elif how=='union':
        dfinter = spatial_overlays(df1, df2, how='intersection')
        dfsym = spatial_overlays(df1, df2, how='symmetric_difference')
        dfunion = dfinter.append(dfsym)
        dfunion.reset_index(inplace=True, drop=True)
        return dfunion
    elif how=='identity':
        dfunion = spatial_overlays(df1, df2, how='union')
        cols1 = df1.columns.tolist()
        cols2 = df2.columns.tolist()
        cols1.remove('geometry')
        cols2.remove('geometry')
        cols2 = set(cols2).intersection(set(cols1))
        cols1 = list(set(cols1).difference(set(cols2)))
        cols2 = [col+'_1' for col in cols2]
        dfunion = dfunion[(dfunion[cols1+cols2].isnull()==False).values]
        return dfunion

def create_row_col_df ( case,
                        lat_source,
                        lon_source,
                        lat_target_int,
                        lon_target_int):
    """
    @ author:                  Shervan Gharari
    @ Github:                  https://github.com/ShervanGharari/EASYMORE
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 GNU-GPLv3
    this fucntion gets the row and colomns of the source netcdf file and returns it
    ----------
    lat_source: numpy array of lat source
    lon_source: numpy array of lon source
    lat_target_int: numpy array of lat source
    lon_target_int: numpy array of lon source
    Returns
    -------
    rows: numpy array, rows from the source file based on the target lat/lon
    cols: numpy array, cols from the source file based on the target lat/lon
    """
    # from numba import jit
    # @jit(nopython=True)
    # def bottleneck(case: int, lat_target_int: np.ndarray, lon_target_int: np.ndarray,
    #                lat_source: np.ndarray, lon_source: np.ndarray) -> tuple:
    #     rows= []
    #     cols= []
    #     for i in np.arange(len(lat_target_int)):
    #         lat_lon_value_diff= np.abs(lat_target_int[i] - lat_source) + np.abs(lon_target_int - lon_source)
    #         if case == 1 or case == 2:
    #             row, col= np.where(lat_lon_value_diff == np.min(lat_lon_value_diff))
    #         else: #case==3
    #             row= np.where(lat_lon_value_diff==np.min(lat_lon_value_diff))[0]
    #             col= row
    #         # print(row,col)
    #         rows.append(row)
    #         cols.append(col)
    #     # out= np.stack([rows, cols], axis=-1)
    #     return (rows, cols)
    locs= [np.unravel_index(((lat_target_int[i]-lat_source)**2+(lon_target_int[i]-lon_source)**2).argmin(), lat_source.shape) for i in range(len(lat_target_int))]
    locs= np.stack(locs)

    # rows, cols= bottleneck(case, lat_target_int, lon_target_int, lat_source, lon_source)
    rows= locs[:,0]; cols= locs[:,1]
    # print(res)
    # rows=res[:,0];cols=res[:,1]
    return np.array(rows), np.array(cols)

def create_remap(case,
                int_df,
                lat_source,
                lon_source):
    """
    @ author:                  Shervan Gharari
    @ Github:                  https://github.com/ShervanGharari/EASYMORE
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 GNU-GPLv3
    this function add the corresponsing row and columns from the source NetCDF file
    Parameters
    ----------
    int_df: intersected data frame that includes the infromation for source and sink
    lat_source: numpy array of source lat
    lon_source: numpy array of source lon
    Returns
    -------
    int_df: dataframe, including the associated rows and cols and EASYMORE case
    """
    # the lat lon from the intersection/remap
    lat_source_int = np.array(int_df['lat_s'])
    lon_source_int = np.array(int_df['lon_s'])
    # call get row and col function
    rows, cols = create_row_col_df (case, lat_source, lon_source, lat_source_int, lon_source_int)
    # add rows and columns
    int_df['rows'] = rows
    int_df['cols'] = cols
    # pass the case to the remap_df
    int_df['easymore_case'] = case
    # save remap_df as csv for future use
    return int_df


def create_output_nc(
    nc_name,
    remap,
    var_time,
    case_name,
    remapped_dim_id,
    remapped_var_id,
    remapped_var_lat,
    remapped_var_lon,
    var_names,
    var_names_remapped,
    output_dir,
    format_list,
    fill_value_list,
    complevel,
    author_name,
    license,
    overwrite_existing_remap=True,
):
        """
        @ author:                  Shervan Gharari
        @ Github:                  https://github.com/ShervanGharari/EASYMORE
        @ author's email id:       sh.gharari@gmail.com
        @ license:                 GNU-GPLv3
        This funciton read different grids and sum them up based on the
        weight provided to aggregate them over a larger area
        """
        print('------REMAPPING------')
        case=int(remap.easymore_case.iloc[0])
        # creating the target_ID_lat_lon
        target_ID_lat_lon = pd.DataFrame()
        target_ID_lat_lon ['ID_t']  = remap ['ID_t']
        target_ID_lat_lon ['lat_t'] = remap ['lat_t']
        target_ID_lat_lon ['lon_t'] = remap ['lon_t']
        target_ID_lat_lon ['order_t'] = remap ['order_t']
        target_ID_lat_lon = target_ID_lat_lon.drop_duplicates()
        target_ID_lat_lon = target_ID_lat_lon.sort_values(by=['order_t'])
        target_ID_lat_lon = target_ID_lat_lon.reset_index(drop=True)
        # prepare the hru_id (here COMID), lat, lon
        hruID_var = np.array(target_ID_lat_lon['ID_t'])
        hruID_lat = np.array(target_ID_lat_lon['lat_t'])
        hruID_lon = np.array(target_ID_lat_lon['lon_t'])
        #
        rows = np.array(remap['rows']).astype(int)
        cols = np.array(remap['cols']).astype(int)
        number_of_target_elements = len(hruID_var)
        # check compression choice
        if isinstance(complevel, int) and (complevel >=1) and (complevel<=9):
            compflag = True
            complevel = complevel
            print('netcdf output file will be compressed at level', complevel)
        else:
            compflag = False
            complevel = 0
            print('netcdf output file will not be compressed.')
        # get the time unit and time var from source
        ncids = nc4.Dataset(nc_name)
        # Check data license, calendar and time units
        nc_att_list = ncids.ncattrs()
        nc_att_list = [each_att for each_att in nc_att_list]
        nc_att_list_lower = [each_att.lower() for each_att in nc_att_list]
        if license == '' and ('license' not in nc_att_list_lower):
            license == 'the original license of the source NetCDF file is not provided'
        if ('license' in nc_att_list_lower):
            if 'license' in nc_att_list:
                license == getattr(ncids, 'license')
            if 'License' in nc_att_list:
                license == getattr(ncids, 'License')
            if 'LICENSE' in nc_att_list:
                license == getattr(ncids, 'LICENSE')
            license == 'Original data license '+license
        if 'units' in ncids.variables[var_time].ncattrs():
            time_unit = ncids.variables[var_time].units
        else:
            sys.exit('units is not provided for the time varibale for source NetCDF of'+ nc_name)
        if 'calendar' in ncids.variables[var_time].ncattrs():
            time_cal = ncids.variables[var_time].calendar
        else:
            sys.exit('calendar is not provided for the time varibale for source NetCDF of'+ nc_name)
        time_var = ncids[var_time][:]
        length_of_time = len(time_var)
        target_date_times = nc4.num2date(time_var,units = time_unit,calendar = time_cal)
        target_name = output_dir + case_name + '_remapped_' + target_date_times[0].strftime("%Y-%m-%d-%H-%M-%S")+'.nc'
        if os.path.exists(target_name):
            if overwrite_existing_remap:
                print('Removing existing remapped .nc file.')
                os.remove(target_name)
            else:
                print('Remapped .nc file already exists. Going to next file.')
        for var in ncids.variables.values():
            if var.name == var_time:
                time_dtype =  str(var.dtype)
        time_dtype_code = 'f8' # initialize the time as float
        if 'float' in time_dtype.lower():
            time_dtype_code = 'f8'
        elif 'int' in time_dtype.lower():
            time_dtype_code = 'i4'
        # reporting
        print('Remapping '+nc_name+' to '+target_name)
        print('Started at date and time '+str(datetime.now()))
        with nc4.Dataset(target_name, "w", format="NETCDF4") as ncid: # creating the NetCDF file
            # define the dimensions
            dimid_N = ncid.createDimension(remapped_dim_id, len(hruID_var))  # limited dimensiton equal the number of hruID
            dimid_T = ncid.createDimension('time', None)   # unlimited dimensiton
            # Variable time
            time_varid = ncid.createVariable('time', time_dtype_code, ('time', ), zlib=compflag, complevel=complevel)
            # Attributes
            time_varid.long_name = var_time
            time_varid.units = time_unit  # e.g. 'days since 2000-01-01 00:00' should change accordingly
            time_varid.calendar = time_cal
            time_varid.standard_name = var_time
            time_varid.axis = 'T'
            time_varid[:] = time_var
            # Variables lat, lon, subbasin_ID
            # lat_varid = ncid.createVariable(remapped_var_lat, 'f8', (remapped_dim_id, ), zlib=compflag, complevel=complevel)
            # lon_varid = ncid.createVariable(remapped_var_lon, 'f8', (remapped_dim_id, ), zlib=compflag, complevel=complevel)
            hruId_varid = ncid.createVariable(remapped_var_id, 'f8', (remapped_dim_id, ), zlib=compflag, complevel=complevel)
            # Attributes
            # lat_varid.long_name = remapped_var_lat
            # lon_varid.long_name = remapped_var_lon
            hruId_varid.long_name = 'shape ID'
            # lat_varid.units = 'degrees_north'
            # lon_varid.units = 'degrees_east'
            hruId_varid.units = '1'
            # lat_varid.standard_name = remapped_var_lat
            # lon_varid.standard_name = remapped_var_lon
            # lat_varid[:] = hruID_lat
            # lon_varid[:] = hruID_lon
            hruId_varid[:] = hruID_var
            # general attributes for NetCDF file
            ncid.Conventions = 'CF-1.6'
            ncid.Author = 'The data were written by ' + author_name
            ncid.License = license
            ncid.History = 'Created ' + time.ctime(time.time())
            ncid.Source = 'Case: ' +case_name + '; remapped by script from library of Shervan Gharari (https://github.com/ShervanGharari/EASYMORE).'
            # write varibales
            for i in np.arange(len(var_names)):
                var_value  = _weighted_average( case,
                                                nc_name,
                                                var_time,
                                                length_of_time,
                                                number_of_target_elements,
                                                target_date_times,
                                                var_names[i],
                                                rows,
                                                cols,
                                                remap)
                # Variables writing
                varid = ncid.createVariable(var_names_remapped[i], format_list[i], ('time',remapped_dim_id ), fill_value = fill_value_list[i], zlib=compflag, complevel=complevel)
                varid [:] = var_value
                # Pass attributes
                if 'long_name' in ncids.variables[var_names[i]].ncattrs():
                    varid.long_name = ncids.variables[var_names[i]].long_name
                if 'units' in ncids.variables[var_names[i]].ncattrs():
                    varid.units = ncids.variables[var_names[i]].units
            # reporting
            print('Ended   at date and time '+str(datetime.now()))
            print('------')

def _weighted_average(case,
                    nc_name,
                    var_time,
                    length_of_time,
                    number_of_target_elements,
                    target_time,
                    varibale_name,
                    rows,
                    cols,
                    mapping_df):
    """
    @ author:                  Shervan Gharari
    @ Github:                  https://github.com/ShervanGharari/EASYMORE
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 GNU-GPLv3
    This function reads the data for a given time and calculates the weighted average
    Arguments
    ---------
    nc_name: string, name of the netCDF file
    target_time: string,
    varibale_name: string, name of varibale from source netcsf file to be remapped
    mapping_df: pandas dataframe, including the row and column of the source data and weight
    Returns
    -------
    weighted_value: a numpy array that has the remapped values from the nc file
    """
    # open dataset
    ds = xr.open_dataset(nc_name)
    # rename time varibale to time
    if var_time != 'time':
        ds = ds.rename({var_time:'time'})
    # prepared the numpy array for ouptut
    weighted_value = np.zeros([length_of_time, number_of_target_elements])
    m = 0 # counter
    for date in target_time: # loop over time
        ds_temp = ds.sel(time=date.strftime("%Y-%m-%d %H:%M:%S"),method="nearest")
        data = np.array(ds_temp[varibale_name])
        data = np.squeeze(data)
        # get values from the rows and cols and pass to np data array
        if case ==1 or case ==2:
            values = data [rows,cols]
        if case ==3:
            values = data [rows]
        values = np.array(values)
        # add values to data frame, weighted average and pass to data frame again
        mapping_df['values'] = values
        mapping_df['values_w'] = mapping_df['weight']*mapping_df['values']
        df_temp = mapping_df.groupby(['order_t'],as_index=False).agg({'values_w': 'sum'})
        df_temp = df_temp.sort_values(by=['order_t'])
        weighted_value [m,:] = np.array(df_temp['values_w'])
        m += 1
    return weighted_value

def process(geo,
            case,
            source_nc,
            source_shp,
            source_shp_ID,
            source_shp_lat,
            source_shp_lon,
            target_shp_gpd) -> pd.DataFrame:
    '''
    process partitioned longitudes and latitudes
    '''
    lat,lon=geo
    # print(lon.shape)
    if (case==1 or case==2) and (source_shp==''):
        source_shp_gpd= lat_lon_SHP(lat, lon)
    elif (case==1 or case==2) and (source_shp!=''):
        source_shp_gpd= gpd.read_file(source_shp).set_crs('epsg:4326')
        source_shp_gpd= add_lat_lon_source_SHP(source_shp_gpd, source_shp_lat,
                                               source_shp_lon, source_shp_ID)
    else: #case 3
        check_source_nc_shp(source_nc, source_shp, source_shp_lat, source_shp_lon,
                            lat, lon)
        source_shp_gpd= gpd.read_file(source_shp).set_crs('epsg:4326')
        source_shp_gpd=add_lat_lon_source_SHP(source_shp_gpd, source_shp_lat,
                                               source_shp_lon, source_shp_ID)
    shp_1= target_shp_gpd.set_crs('epsg:4326')
    shp_2= source_shp_gpd.set_crs('epsg:4326')

    #reproject to equal area
    if (str(shp_1.crs).lower() == str(shp_2.crs).lower()) and ('epsg:4326' in str(shp_1.crs).lower()):
        shp_1= shp_1.to_crs('epsg:6933')
        shp_2= shp_2.to_crs('epsg:6933')
    else:
        print('source projection: %s target projection: %s'%(str(shp_2.crs.lower()), str(shp_1.crs.lower())))
        sys.exit('The projection for source and target shapefile are not WGS84, please revise, assign')
    shp_int= intersection_shp(shp_1, shp_2)
    if shp_int.empty:
        return pd.DataFrame()
    else:
        shp_int = shp_int.sort_values(by=['S_1_ID_t'])
        shp_int = shp_int.to_crs ("epsg:4326")
        shp_int= shp_int.drop(columns=['geometry'])
        dict_rename = {
            'S_1_ID_t' : 'ID_t',
            'S_1_lat_t': 'lat_t',
            'S_1_lon_t': 'lon_t',
            'S_1_order': 'order_t',
            'S_2_ID_s' : 'ID_s',
            'S_2_lat_s': 'lat_s',
            'S_2_lon_s': 'lon_s',
            'AP1N'     : 'weight'
        }
        shp_int= shp_int.rename(columns=dict_rename)
        shp_int= pd.DataFrame(shp_int)
        int_df= shp_int.copy()
        lat_source= lat
        lon_source= lon
        int_df= create_remap(case, int_df, lat_source, lon_source)

        return int_df



def esmr(
        case_name                 =  'case_temp',   # name of the case
        target_shp                =  '',            # sink/target shapefile
        target_shp_ID             =  '',            # name of the column ID in the sink/target shapefile
        target_shp_lat            =  '',            # name of the column latitude in the sink/target shapefile
        target_shp_lon            =  '',            # name of the column longitude in the sink/target shapefile
        source_nc                 =  '',            # name of nc file to be remapped
        var_names                 =  [],            # list of varibale names to be remapped from the source NetCDF file
        var_lon                   =  '',            # name of varibale longitude in the source NetCDF file
        var_lat                   =  '',            # name of varibale latitude in the source NetCDF file
        var_time                  =  'time',        # name of varibale time in the source NetCDF file
        var_ID                    =  '',            # name of vriable ID in the source NetCDF file
        var_names_remapped        =  [],            # list of varibale names that will be replaced in the remapped file
        source_shp                =  '',            # name of source shapefile (essential for case-3)
        source_shp_lat            =  '',            # name of column latitude in the source shapefile
        source_shp_lon            =  '',            # name of column longitude in the source shapefile
        source_shp_ID             =  '',            # name of column ID in the source shapefile
        remapped_var_id           =  'ID',          # name of the ID variable in the new nc file; default 'ID'
        remapped_var_lat          =  'latitude',    # name of the latitude variable in the new nc file; default 'latitude'
        remapped_var_lon          =  'longitude',   # name of the longitude variable in the new nc file; default 'longitude'
        remapped_dim_id           =  'ID',          # name of the ID dimension in the new nc file; default 'ID'
        overwrite_existing_remap  =  True,          # Flag to automatically overwrite existing remapping files. If 'False', aborts the remapping procedure if a file is detected
        temp_dir                  =  './temp/',     # temp_dir
        output_dir                =  './output/',   # output directory
        format_list               =  ['f8'],        # float for the remapped values
        fill_value_list           =  ['-9999'],     # missing values set to -9999
        remap_csv                 =  '',            # name of the remapped file if provided
        author_name               =  '',            # name of the authour
        license                   =  '',            # data license
        tolerance                 =  1e-5,        # tolerance
        save_csv                  =  False,         # save csv
        save_shp                  =  False,         # save intermediate shape file
        ncpus                     =  8,             #number of cpus for dask
        sort_ID                   =  False,         # to sort the remapped based on the target shapfile ID; target_shp_ID should be given
        complevel                 =  4,             # netcdf compression level from 1 to 9. Any other value or object will mean no compression.
        version                   =  '0.0.3',       # version of the easymore
        
):
    
    from functools import partial
    import multiprocessing
    print('EASYMORE version '+version+ ' is initiated.')
    author_name, format_list, fill_value_list, var_names_remapped= check_easymore_input(
                                                                                    temp_dir,
                                                                                    output_dir,
                                                                                    author_name,
                                                                                    var_names,
                                                                                    format_list,
                                                                                    fill_value_list,
                                                                                    remap_csv,
                                                                                    var_names_remapped)
    # if var_names_remapped:
    if remap_csv == '' and not os.path.exists(temp_dir+case_name+'_remapping.csv'):
        target_shp_gpd = gpd.read_file(target_shp).set_crs('epsg:4326')
        target_shp_gpd = check_target_shp(target_shp_gpd, target_shp_ID, 
                                          target_shp_lat, target_shp_lon,
                                          sort_ID)
        #check source nc file
        check_source_nc(source_nc, var_names, var_lat, var_lon, var_time, tolerance)
        # find the case
        case, ID, lat, lon, lat_expanded, lon_expanded= NetCDF_SHP_lat_lon(source_nc,
                                                                           var_names,
                                                                           var_ID,
                                                                           var_lon,
                                                                           var_lat,
                                                                           tolerance)

        if lat_expanded is not None and lon_expanded is not None:
            lat= lat_expanded; lon= lon_expanded
        if case == 1 or case == 2:
            # split large lat/lon array into 2D subarray, so that we can apply parallel computing
            idx_lat= np.arange(lat.shape[0])
            idx_lon= np.arange(lon.shape[1])
            par_idx_lat= np.array_split(idx_lat, ncpus)
            par_idx_lon= np.array_split(idx_lon, ncpus)
            tasks= []
            tasks= [[lat[_idx_lat[0]:_idx_lat[-1], _idx_lon[0]:_idx_lon[-1]], lon[_idx_lat[0]:_idx_lat[-1], _idx_lon[0]:_idx_lon[-1]]] for _idx_lat in par_idx_lat for _idx_lon in par_idx_lon]
            # lat and lon are 2D array
            # loc= np.stack([lat,lon], axis=-1) # stack into (lat_dim, lon_dim, 2)
            # chunks_by_cores= (lon.shape[0]//ncpus, lon.shape[1]//ncpus, 2)
            # loc_da= da.from_array(loc, chunks=chunks_by_cores) #convert to dask array

            func= partial(process,  
                          case=case,
                          source_nc=source_nc, 
                          source_shp=source_shp,
                          source_shp_ID=source_shp_ID,
                          source_shp_lat=source_shp_lat,
                          source_shp_lon=source_shp_lon,
                          target_shp_gpd=target_shp_gpd)
            with multiprocessing.Pool(ncpus) as pool:
                res=pool.map(func, tasks)
            df= pd.concat(res)
            df.to_csv(temp_dir+case_name+'_remapping.csv')

    else:
        df= pd.read_csv(temp_dir+case_name+'_remapping.csv')
    # output nc file
    ncfiles= glob.glob(source_nc)
    func= partial(create_output_nc,
                remap= df,
                var_time=var_time,
                case_name= case_name,
                remapped_dim_id=remapped_dim_id,
                remapped_var_id=remapped_var_id,
                remapped_var_lat=remapped_var_lat,
                remapped_var_lon=remapped_var_lon,
                var_names= var_names,
                var_names_remapped=var_names_remapped,
                output_dir=output_dir,
                format_list=format_list,
                fill_value_list=fill_value_list,
                complevel=complevel,
                author_name=author_name,
                license=license,
                overwrite_existing_remap=overwrite_existing_remap
                )
    # start writing output
    with multiprocessing.Pool(ncpus) as pool:
        pool.map(func, ncfiles)

if __name__== '__main__':
    # test
    import time
    start= time.time()
    pfaf = 74
    esmr(
        'CREST_%s'%pfaf,
        '/home/ZhiLi/US_routing/cat_MERIT/cat_pfaf_%d_MERIT_Hydro_v07_Basins_v01_bugfix1.shp'%pfaf,
        'COMID',
        '',
        '',
        '/media/scratch/ZhiLi/CREST_output/nc_out/*.nc',
        ['runoff','subrunoff'],
        'lon',
        'lat',
        'time',
        temp_dir='/home/ZhiLi/US_routing/preprocess/temp/',
        output_dir='/media/scratch/ZhiLi/CREST_output/nc_out/',
        ncpus=8,
    )
    end= time.time()
    print('completed! processed time %.2f hours'%((end-start)/3600.))