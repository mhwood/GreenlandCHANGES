
from scipy.interpolate import interp2d
import xarray as xr
import numpy as np
import os
import gdal
import requests
from toolbox.time import YMD_to_DecYr

#######################################################################################
#These are the scripts for finding a list of GIMP files

def find_gimp_dem_files_in_domain(GD_object):

    def check_overlap(extents1, extents2):
        overlap = True
        if extents1[2] < extents2[0] or extents1[0] > extents2[2] or extents1[3] < extents2[1] or extents1[1] > \
                extents2[3]:
            overlap = False
        return (overlap)

    import changes.reference.gimp_domains as gid

    dem_files = []
    for row in range(6):
        for col in range(6):
            gimp_id = 'gimpdem'+str(row)+'_'+str(col)+'_v01.1'
            gimp_tile_extent = gid.gimp_id_to_extents[gimp_id]

            if check_overlap(GD_object.extents,gimp_tile_extent):
                dem_files.append(gimp_id)

    return(dem_files)


#######################################################################################
#These are the scripts for downloading the GIMP data

def download_GIMP_tile(dataFolder,dem_file):

    outputFile=os.path.join(dataFolder,dem_file+'.tif')
    url = 'https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0645.001/2003.02.20/'+dem_file+'.tif'
    resp = requests.get(url)
    open(outputFile, 'wb').write(resp.content)

def download_gimp_files(GD_object,dem_files):

    'https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0645.001/2003.02.20/gimpdem0_0_v01.1.tif'


    dataFolder = os.path.join(GD_object.data_folder,'Elevation','GIMP','Data')

    for dd in range(len(dem_files)):
        dem_file = dem_files[dd]
        if GD_object.print_sub_outputs:
            print('            Checking file '+dem_file)

        download_file=True
        if dem_file+'.tif' in os.listdir(dataFolder):
            download_file=False

        if download_file:
            if GD_object.print_sub_outputs:
                print('              Downloading file ' + dem_file)
            download_GIMP_tile(dataFolder, dem_file)
        else:
            if GD_object.print_sub_outputs:
                print('              File already downloaded')


#######################################################################################
#These are the scripts for creating the layer from the GIMP data

def create_gimp_dem_field(valid_points,regionX,regionY):

    #limit the points to the region domain
    dem_data = valid_points
    dem_data = dem_data[np.logical_and(dem_data[:, 0] <= np.max(regionX), dem_data[:, 0] >= np.min(regionX)), :]
    dem_data = dem_data[np.logical_and(dem_data[:, 1] <= np.max(regionY), dem_data[:, 1] >= np.min(regionY))]

    # create the glacier field
    dem_grid = np.zeros((len(regionY), len(regionX)))
    count_grid = np.zeros((len(regionY), len(regionX)))
    print('                Glacier field size: ' + str(np.shape(dem_grid)))

    #fill in the glacier field by averaging points based on their closest point
    for dd in range(np.shape(dem_data)[0]):
        x_index = np.argmin(np.abs(regionX-dem_data[dd,0]))
        y_index = np.argmin(np.abs(regionY-dem_data[dd,1]))
        dem_grid[y_index,x_index]+=dem_data[dd,2]
        count_grid[y_index,x_index]+=1

    dem_grid[count_grid>0] = dem_grid[count_grid>0]/count_grid[count_grid>0]
    dem_grid[count_grid == 0] = -99

    return(dem_grid)

def get_gimp_layer(GD_object,dem_file_names):

    dem_grid = -99 * np.ones((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))

    for dd in range(len(dem_file_names)):

        message = '            Adding in data for tile ' + dem_file_names[dd]
        if GD_object.print_sub_outputs:
            print(message)

        ds = gdal.Open(os.path.join(GD_object.data_folder,'Elevation','GIMP','Data',dem_file_names[dd]+'.tif'))
        dem = np.array(ds.GetRasterBand(1).ReadAsArray())
        trans = ds.GetGeoTransform()
        min_x = trans[0]
        max_x = trans[0] + trans[1] * np.shape(dem)[1]
        max_y = trans[3]
        min_y = trans[3] + trans[5] * np.shape(dem)[0]
        x=np.arange(min_x,max_x,trans[1])
        y=np.arange(max_y,min_y,trans[5])

        min_x_index = np.argmin(np.abs(np.min(GD_object.elevation_grid_x)-x))-2
        if min_x_index<0:
            min_x_index=0
        max_x_index = np.argmin(np.abs(np.max(GD_object.elevation_grid_x)-x))+2
        if max_x_index>len(x):
            max_x_index=len(x)
        max_y_index = np.argmin(np.abs(np.min(GD_object.elevation_grid_y)-y))+2
        if max_y_index>len(y):
            max_y_index=len(y)
        min_y_index = np.argmin(np.abs(np.max(GD_object.elevation_grid_y)-y))-2
        if min_y_index<0:
            min_y_index=0

        x = x[min_x_index:max_x_index]
        y = y[min_y_index:max_y_index]
        dem = dem[min_y_index:max_y_index,:]
        dem = dem[:, min_x_index:max_x_index]

        message = '              Creating an interpolation field'
        if GD_object.print_sub_outputs:
            print(message)
        set_int = interp2d(x,y,dem)

        message = '              Sampling the interpolation field'
        if GD_object.print_sub_outputs:
            print(message)
        for row in range(len(GD_object.elevation_grid_y)):
            if GD_object.elevation_grid_y[row] >= min_y and GD_object.elevation_grid_y[row] <= max_y:
                for col in range(len(GD_object.elevation_grid_x)):
                    if GD_object.elevation_grid_x[col]>=min_x and GD_object.elevation_grid_x[col]<=max_x:
                        interp_val = set_int(GD_object.elevation_grid_x[col],GD_object.elevation_grid_y[row])
                        if interp_val>-20:
                            dem_grid[row,col] = interp_val

    dem_date = '20050630'
    dem_decYr = YMD_to_DecYr(int(dem_date[:4]), int(dem_date[4:6]), int(dem_date[6:8]))

    return(dem_grid, dem_date, dem_decYr)


#######################################################################################
#These are the scripts for stacking the glacier fields into one

def save_gimp_layers(GD_object,dem_layer, dem_date, dem_decYr, dem_file_names):

    data_vars = {dem_date:(['y','x'],dem_layer)}

    swath = xr.Dataset(data_vars,coords={'y': GD_object.elevation_grid_y,'x': GD_object.elevation_grid_x})

    swath[dem_date].attrs['file_name'] = ','.join(dem_file_names)
    swath[dem_date].attrs['dec_yr'] = dem_decYr
    swath[dem_date].attrs['note'] = 'DEM date is an approximate center date in 2002-2009'

    output_file = os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Data',GD_object.region_name+' GIMP Elevation Grids.nc')
    message = '        Outputting data to ' + output_file
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    swath.to_netcdf(output_file)

#######################################################################################
#This is the routine to run the full compilation

def generate_gimp_dataset(GD_object):

    message = '    Running compilation for the GIMP data'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)

    # step 1: get a list of files for the region
    message = '        Finding a list of GIMP files which overlap the domain'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    dem_file_names = find_gimp_dem_files_in_domain(GD_object)

    # step 2: download the GIMP data
    message = '        Found ' + str(len(dem_file_names)) + ' files'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    download_gimp_files(GD_object, dem_file_names)

    if GD_object.create_elevation_stacks:

        message = '        Resampling the GIMP data onto to the regional domain'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

        # step 3: stack the gimp data into layers
        dem_layer, dem_date, dem_decYr = get_gimp_layer(GD_object,dem_file_names)

        # print('        Saving the GIMP data as an nc file')
        # step 4: save the gimp layer to an nc file
        save_gimp_layers(GD_object,dem_layer, dem_date, dem_decYr, dem_file_names)