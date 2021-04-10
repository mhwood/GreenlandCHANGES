
from datetime import datetime
import numpy as np
import requests
import os
from osgeo import gdal
from osgeo import osr
import xarray as xr
import matplotlib.pyplot as plt
from ....toolbox.time import YMD_to_DecYr
from ....toolbox.resample.interpolation import reproject_and_interpolate_onto_grid

def read_kms_dem(kms_dem_file_path):
    lon_0 = -75
    lon_spacing = 0.06
    n_cols = 1084
    lon = np.arange(lon_0, lon_0 + (n_cols) * lon_spacing, lon_spacing)

    lat_0 = 59.5
    lat_spacing = .02
    n_rows = 1226
    lat = np.arange(lat_0, lat_0 + (n_rows) * lat_spacing, lat_spacing)

    grid = np.fromfile(kms_dem_file_path, dtype='<f2')
    grid = np.reshape(grid, (n_rows, n_cols))
    grid = np.flipud(grid)
    grid[np.isnan(grid)] = 0

    grid_max = np.max(grid)
    reported_max = 3251.4

    grid = grid * reported_max / grid_max

    return(lon,lat,grid,lon_spacing,lat_spacing)

#######################################################################################
#These are the scripts for downloading the list of arctic dem files

def check_file_download(GD_object):

    # check if the data is there
    download_file = True
    if 'grnlnd_dem_wgs84.dat' in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','KMS')):
        download_file = False

    return(download_file)

def download_kms_dem(GD_object,download_url):
    output_file = os.path.join(GD_object.data_folder, 'Elevation','KMS','grnlnd_dem_wgs84.dat')

    with requests.get(download_url, stream=True) as r:
        r.raise_for_status()
        with open(output_file, 'wb') as f:                   #open local output file to write binary to, name returned filehandle f
            for chunk in r.iter_content(chunk_size=8192):   #for each 8192byte chunk in the incoming stream
                f.write(chunk)

    return(download_url)

def save_kms_dem_as_tif(GD_object):
    if 'grnlnd_dem_wgs84.tif' not in os.listdir(os.path.join(GD_object.data_folder, 'Elevation', 'KMS')):
        input_file = os.path.join(GD_object.data_folder, 'Elevation', 'KMS', 'grnlnd_dem_wgs84.dat')
        output_file = os.path.join(GD_object.data_folder, 'Elevation', 'KMS', 'grnlnd_dem_wgs84.tif')

        lon,lat,grid,lon_spacing,lat_spacing = read_kms_dem(input_file)

        geotransform = (np.min(lon), lon_spacing, 0, np.max(lat), 0, -lat_spacing)
        output_raster = gdal.GetDriverByName('GTiff').Create(output_file, len(lon), len(lat), 1,
                                                             gdal.GDT_Float32)
        output_raster.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        output_raster.SetProjection(srs.ExportToWkt())
        output_raster.GetRasterBand(1).WriteArray(np.flipud(grid))

def download_and_resave_kms_file(GD_object):

    download_file = check_file_download(GD_object)

    download_url = 'https://daacdata.apps.nsidc.org/DATASETS/nsidc0052_sar_mosaic_elevation_v1/grnlnd_dem_wgs84.dat'

    if download_file:
        if GD_object.print_sub_outputs:
            print('              Downloading KMS DEM')
        download_kms_dem(GD_object,download_url)
    else:
        if GD_object.print_sub_outputs:
            print('              File already obtained')

    if GD_object.save_tif_copy_of_kms_dem:
        save_kms_dem_as_tif(GD_object)

    return(download_url)

#######################################################################################
#These are the scripts for creating the layers in the KMS data

def get_kms_layer(GD_object):

    kms_file = os.path.join(GD_object.data_folder, 'Elevation', 'KMS', 'grnlnd_dem_wgs84.dat')
    lon,lat,grid,lon_spacing,lat_spacing = read_kms_dem(kms_file)

    dem_layer = reproject_and_interpolate_onto_grid(lon,lat,grid, 4326,
                                                         GD_object.elevation_grid_x,GD_object.elevation_grid_y,
                                                         3413, print_status_messages=True)
    plt.imshow(dem_layer)
    plt.show()

    start_datetime = datetime(1991,1,1)
    stop_datetime = datetime(1995, 12, 31)

    center_time = start_datetime + (stop_datetime - start_datetime) / 2
    start_time_str = str(start_datetime.year) + '{:02d}'.format(start_datetime.month) + '{:02d}'.format(start_datetime.day)
    stop_time_str = str(stop_datetime.year) + '{:02d}'.format(stop_datetime.month) + '{:02d}'.format(stop_datetime.day)
    center_time_str = str(center_time.year) + '{:02d}'.format(center_time.month) + '{:02d}'.format(center_time.day)

    print('Check the exact dates in Ekholm, 1996')

    dem_decYr = YMD_to_DecYr(int(center_time_str[:4]), int(center_time_str[4:6]), int(center_time_str[6:8]))

    print(dem_decYr)
    return(dem_layer, center_time_str, dem_decYr)

#######################################################################################
#These are the scripts to save the layers as a stack

def save_kms_layer(GD_object,dem_layer, dem_date, dem_decYr, url):

    data_vars = {}
    data_vars[dem_date]=(['y','x'],dem_layer)

    swath = xr.Dataset(data_vars,coords={'y': GD_object.elevation_grid_y,'x': GD_object.elevation_grid_x})

    swath[dem_date].attrs['file_name'] = url
    swath[dem_date].attrs['dec_yr'] = dem_decYr

    if len(GD_object.kms_output_file)>2:
        output_file = GD_object.kms_output_file
    else:
        output_file = os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Data',GD_object.region_name+' KMS Elevation Grids.nc')
    print('        Outputting data to ' + output_file)

    swath.to_netcdf(output_file)

#######################################################################################
#This is the main script to run the whole process

def generate_kms_dataset(GD_object):

    message = '    Running compilation for the KMS data'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)

    # step 2: download the data and down-sample it at the requested resolution
    message = '        Downloading DEM file (and creating nc if requested)'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    url = download_and_resave_kms_file(GD_object)

    if GD_object.create_elevation_stacks:

        if GD_object.overwrite_existing_elevation_data:
            continue_to_stack = True
        else:
            if len(GD_object.kms_output_file) > 2:
                output_file = GD_object.kms_output_file
            else:
                output_file = os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation', 'Data',
                                           GD_object.region_name + ' KMS Elevation Grid.nc')
            if os.path.isfile(output_file):
                continue_to_stack = False
            else:
                continue_to_stack = True

        if continue_to_stack:

            # step 3: stack the kms data into layers
            print('        Interpolating the KMS DEM onto to the regional domain')
            dem_layer, dem_dates, dem_decYr = get_kms_layer(GD_object)

            # step 4: save the kms layer to an nc file
            save_kms_layer(GD_object,dem_layer, dem_dates, dem_decYr, url)

