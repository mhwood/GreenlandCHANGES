
from datetime import datetime, timedelta
import numpy as np
from osgeo import ogr
import requests
import urllib
import os
from osgeo import gdal
import zipfile
from osgeo import osr
import xarray as xr
import shutil
from ....toolbox.time import YMD_to_DecYr
from ....toolbox.resample.interpolation import reproject_and_interpolate_onto_grid
from ....toolbox.reprojection import reproject_polygon
from ....toolbox.series import series_to_N_points
from ..reference.tandemx_domains import tandemx_id_to_extents

#######################################################################################
#These are the scripts for finding a list of TanDEMX files

def find_tandemx_dem_files_from_dict(GD_object):

    def bboxToWKT(bbox):
        output = "POLYGON (("
        for b in range(np.shape(bbox)[0]):
            output += str(bbox[b, 0]) + " " + str(bbox[b, 1]) + ", "
        output = output[:-2] + "))"
        return (output)

    bbox = np.array([[GD_object.extents[0],GD_object.extents[1]],
                     [GD_object.extents[2],GD_object.extents[1]],
                     [GD_object.extents[2],GD_object.extents[3]],
                     [GD_object.extents[0],GD_object.extents[3]],
                     [GD_object.extents[0],GD_object.extents[1]]])
    bbox = series_to_N_points(bbox,100)
    bboxWKT = bboxToWKT(bbox)
    poly1 = ogr.CreateGeometryFromWkt(bboxWKT)

    ####################################################################################################################
    # Loop through the metadata file to see which files pertain to the sample area

    message = '            Searching through the TanDEMX domains to find overlapping files'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    tandemx_file_names = list(tandemx_id_to_extents.keys())
    dem_files = []

    for s in range(len(tandemx_file_names)):
        tdx_extent = tandemx_id_to_extents[tandemx_file_names[s]]
        bbox = np.array([[tdx_extent[0], tdx_extent[1]],
                         [tdx_extent[2], tdx_extent[1]],
                         [tdx_extent[2], tdx_extent[3]],
                         [tdx_extent[0], tdx_extent[3]],
                         [tdx_extent[0], tdx_extent[1]]])
        bbox = reproject_polygon(bbox,4326,3413)
        bbox_tmp = np.copy(bbox)
        bbox_tmp[:,0] = bbox[:,1]
        bbox_tmp[:,1] = bbox[:,0]
        bboxWKTtdx = bboxToWKT(bbox_tmp)
        poly2 = ogr.CreateGeometryFromWkt(bboxWKTtdx)
        intersection = poly1.Intersection(poly2)
        intersection = intersection.ExportToWkt()
        if "EMPTY" not in intersection:
            fileID = tandemx_file_names[s]
            tile = fileID.split('_')[-1][:-1]+'0'
            lon_str = tile[3:]
            lat_str = tile[:3]
            url = 'https://download.geoservice.dlr.de/TDM90/files/'+lat_str+'/'+\
                  lon_str+'/'+fileID
            dem_files.append([tile, fileID, url])
    #         fileID = filePath.split('/')[-1][:-7]
    #         if '.edu' in filePath:
    #             try:
    #                 resp = urllib.request.urlopen(filePath)   #TanDEMX file list may have errors
    #             except urllib.error.HTTPError:                #This code filters out files which do not exist
    #                 pass
    #             else:
    #                 dem_files.append([tile, fileID, url])

    output = ''
    for strip in dem_files:
        output += ','.join(strip) + '\n'

    f = open(os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation','Metadata',GD_object.region_name + ' TanDEMX Files.csv'), 'w')
    f.write(output)
    f.close()
    return (dem_files)

def find_tandemx_dem_files_in_domain(GD_object):

    if GD_object.region_name + ' TanDEMX Files.csv' not in os.listdir(
            os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation','Metadata')):

        dem_files = find_tandemx_dem_files_from_dict(GD_object)

    else:
        message = '            Reading file list generated in a previous iteration of this process'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

        f=open(os.path.join(GD_object.project_folder, GD_object.region_name,
                            'Elevation','Metadata',GD_object.region_name + ' TanDEMX Files.csv'))
        lines=f.read()
        f.close()
        lines=lines.split('\n')
        dem_files=[]
        for line in lines:
            line=line.split(',')
            if len(line)>1:
                dem_files.append(line)

    return(dem_files)

#######################################################################################
#These are the scripts for downloading the list of arctic dem files

def check_file_download(GD_object,file_name,tile_name):
    download_file = True

    # check if the raw data is there
    raw_data_present = False
    if tile_name in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','TanDEMX','Data','90m_tiles')):
        if file_name + '.zip' in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','TanDEMX','Data','90m_tiles',tile_name)):
            raw_data_present = True
        if file_name + '_DEM.tif' in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','TanDEMX','Data','90m_tiles',tile_name)):
            raw_data_present = True

    # check if the resampled data is there
    resampled_data_present = False
    if tile_name in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','TanDEMX','Data','Regridded_'+str(int(GD_object.elevation_grid_posting))+'m_tiles')):
        if file_name + '_dem_regridded_' + str(int(GD_object.elevation_grid_posting)) + 'm.tif' in os.listdir(
                os.path.join(GD_object.data_folder, 'Elevation','TanDEMX','Data','Regridded_'+str(int(GD_object.elevation_grid_posting))+'m_tiles',tile_name)):
            resampled_data_present = True

    # if neither of these files is present, then download the data
    if raw_data_present or resampled_data_present:
        download_file = False

    return(download_file,resampled_data_present)

def download_tandemx_tile(GD_object, tile_name, file_name, download_url):

    if tile_name not in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','TanDEMX','Data','90m_tiles')):
        os.mkdir(os.path.join(GD_object.data_folder, 'Elevation','TanDEMX','Data','90m_tiles',tile_name))

    output_file = os.path.join(GD_object.data_folder, 'Elevation','TanDEMX','Data','90m_tiles',tile_name,file_name + '.zip')

    # resp = requests.get(download_url, auth=(GD_object.tandemx_username,GD_object.tandemx_password))
    # open(output_file, 'wb').write(resp.content)

    with requests.get(download_url,auth=(GD_object.tandemx_username,GD_object.tandemx_password), stream=True) as r:
        r.raise_for_status()
        with open(output_file, 'wb') as f:                   #open local output file to write binary to, name returned filehandle f
            for chunk in r.iter_content(chunk_size=8192):   #for each 8192byte chunk in the incoming stream
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                # if chunk:
                f.write(chunk)

def resample_tandemx_tile(GD_object, tile_name, file_name):


    tile_path = os.path.join(GD_object.data_folder, 'Elevation','TanDEMX','Data','90m_tiles',tile_name)
    #unxip the file if necessary
    if file_name+'_DEM.tif' not in os.listdir(tile_path):
        print('                Unzipping the data')
        zip = zipfile.ZipFile(os.path.join(tile_path,file_name+'.zip'))
        zip_files = zip.namelist()
        dir_name = zip_files[0].split('/')[0]

        with zipfile.ZipFile(os.path.join(tile_path, file_name + '.zip')) as zf:
            zf.extract(dir_name + '/DEM/' + file_name + '_DEM.tif', tile_path)
        os.rename(os.path.join(tile_path, dir_name, 'DEM', file_name + '_DEM.tif'),
                  os.path.join(tile_path, file_name + '_DEM.tif'))

        with zipfile.ZipFile(os.path.join(tile_path,file_name+'.zip')) as zf:
            zf.extract(dir_name+'/'+file_name+'.xml', tile_path)
        os.rename(os.path.join(tile_path,dir_name,file_name+'.xml'),
                  os.path.join(tile_path,file_name+'.xml'))

        shutil.rmtree(os.path.join(tile_path,dir_name))
        os.remove(os.path.join(tile_path,file_name+'.zip'))

    print('                Working on the down sample')
    print('                  Reading in the file')
    # read in the raster
    ds = gdal.Open(os.path.join(GD_object.data_folder, 'Elevation','TanDEMX','Data','90m_tiles',tile_name,file_name + '_DEM.tif'))
    dem = np.array(ds.GetRasterBand(1).ReadAsArray())
    transform = ds.GetGeoTransform()
    ds = None


    demX = np.arange(transform[0], transform[0] + transform[1] * np.shape(dem)[1], transform[1])
    demY = np.arange(transform[3], transform[3] + transform[5] * np.shape(dem)[0], transform[5])

    bbox = np.array([[np.min(demX), np.min(demY)],
                     [np.max(demX), np.min(demY)],
                     [np.max(demX), np.max(demY)],
                     [np.min(demX), np.max(demY)]])
    # bbox = series_to_N_points(bbox, 100)
    bbox = reproject_polygon(bbox,4326,3413)
    bbox_tmp = np.copy(bbox)
    bbox_tmp[:, 0] = bbox[:, 1]
    bbox_tmp[:, 1] = bbox[:, 0]
    bbox = bbox_tmp
    min_X = np.min(bbox[:, 0])
    # max_X = min_X+20000
    max_X = np.max(bbox[:, 0])
    min_Y = np.min(bbox[:, 1])
    # max_Y = min_Y + 20000
    max_Y = np.max(bbox[:, 1])

    rounding_level = int(np.ceil(np.log10(GD_object.elevation_grid_posting)))
    min_X = np.round(min_X, -rounding_level)
    min_Y = np.round(min_Y, -rounding_level)
    max_X = np.round(max_X, -rounding_level)
    max_Y = np.round(max_Y, -rounding_level)
    new_dem_X = np.arange(min_X, max_X + GD_object.elevation_grid_posting, GD_object.elevation_grid_posting)
    new_dem_Y = np.arange(min_Y, max_Y + GD_object.elevation_grid_posting, GD_object.elevation_grid_posting)

    resampled_grid = reproject_and_interpolate_onto_grid(demX, demY, dem, 4326,
                                        new_dem_X, new_dem_Y,
                                        3413, print_status_messages=True)


    #save the raster as a tif
    if tile_name not in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','TanDEMX',
                                                'Data','Regridded_'+str(int(GD_object.elevation_grid_posting))+'m_tiles')):
        os.mkdir(os.path.join(GD_object.data_folder, 'Elevation','TanDEMX',
                              'Data','Regridded_'+str(int(GD_object.elevation_grid_posting))+'m_tiles',tile_name))

    outputFile = os.path.join(GD_object.data_folder, 'Elevation','TanDEMX',
                              'Data','Regridded_'+str(int(GD_object.elevation_grid_posting))+'m_tiles',
                              tile_name,file_name+'_dem_regridded_'+str(int(GD_object.elevation_grid_posting))+'m.tif')
    geotransform = (np.min(new_dem_X), GD_object.elevation_grid_posting, 0, np.max(new_dem_Y), 0, -GD_object.elevation_grid_posting)
    output_raster = gdal.GetDriverByName('GTiff').Create(outputFile, len(new_dem_X), len(new_dem_Y), 1, gdal.GDT_Float32)
    output_raster.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3413)
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(np.flipud(resampled_grid))

def download_and_resample_tandemx_files(GD_object,dem_files):

    for dd in range(len(dem_files)):
        dem_file = dem_files[dd]
        tile_name = dem_file[0]
        file_name = dem_file[1]
        url = dem_file[2]+'.zip'
        if GD_object.print_sub_outputs:
            print('            Checking file '+file_name+' ('+str(tile_name)+', '+str(dd+1)+' of '+str(len(dem_files))+')')

        #check to see if the file needs to be downloaded
        if GD_object.overwrite_existing_elevation_data:
            download_file=True
            resampled_data_present=False
        else:
            download_file,resampled_data_present = check_file_download(GD_object,file_name,tile_name)

        #downloading...
        if download_file:
            if GD_object.print_sub_outputs:
                print('              Downloading file...')
            download_tandemx_tile(GD_object, tile_name, file_name, url)
        else:
            if GD_object.print_sub_outputs:
                print('              File already obtained')

        #resample the data if its not already done
        if not resampled_data_present:
            if GD_object.print_sub_outputs:
                print('              Downsampling file... ')
            resample_tandemx_tile(GD_object, tile_name, file_name)
            resampled_data_present=True
        else:
            if GD_object.print_sub_outputs:
                print('              File already resampled')

#######################################################################################
#These are the scripts for creating the layers in the TanDEMX data

def read_tandemx_date_span_from_xml(xml_path):

    f = open(xml_path)
    lines=f.read()
    f.close()
    lines = lines.split('\n')
    for line in lines:
        if 'startTime' in line:
            startLine = line
            startLine = startLine.split('<startTime>')[1]
            startLine = startLine.split('</startTime>')[0][:-1]
            startLine = startLine.replace('T',' ')
            start_datetime= datetime.strptime(startLine, '%Y-%m-%d %H:%M:%S.%f')
        if 'stopTime' in line:
            stopLine = line
            stopLine = stopLine.split('<stopTime>')[1]
            stopLine = stopLine.split('</stopTime>')[0][:-1]
            stopLine = stopLine.replace('T',' ')
            stop_datetime = datetime.strptime(stopLine, '%Y-%m-%d %H:%M:%S.%f')
    center_time = start_datetime + (stop_datetime-start_datetime)/2
    # time_span = (stop_datetime-start_datetime)

    start_time_str = str(start_datetime.year)+'{:02d}'.format(start_datetime.month)+'{:02d}'.format(start_datetime.day)
    stop_time_str = str(stop_datetime.year) + '{:02d}'.format(stop_datetime.month) + '{:02d}'.format(stop_datetime.day)
    center_time_str = str(center_time.year) + '{:02d}'.format(center_time.month) + '{:02d}'.format(center_time.day)

    return(center_time_str,start_time_str,stop_time_str)

def get_tandemx_dates_and_file_sets(GD_object,dem_files):

    all_dates = []
    for dem_file in dem_files:
        file_path=os.path.join(GD_object.data_folder,'Elevation','TanDEMX','Data',
                               'Regridded_50m_tiles',dem_file[0],
                               dem_file[1]+'_dem_regridded_50m.tif')
        xml_path = os.path.join(GD_object.data_folder, 'Elevation', 'TanDEMX', 'Data',
                                 '90m_tiles', dem_file[0],
                                 dem_file[1] + '.xml')
        center_time_str,start_time_str,stop_time_str = read_tandemx_date_span_from_xml(xml_path)
        all_dates.append(center_time_str)

    unique_dates = []
    for date in all_dates:
        if date not in unique_dates:
            unique_dates.append(date)

    file_sets = []
    for date in unique_dates:
        date_set = []
        for df in range(len(dem_files)):
            if all_dates[df] == date:
                date_set.append(dem_files[df])
        file_sets.append(date_set)
    return(unique_dates,file_sets)

def read_file_to_grid(GD_object,tile_name,file_name,dem_grid,count_grid):

    print_debug_messages = False

    regionYcopy=np.flipud(GD_object.elevation_grid_y)

    #read in the down sampled grid

    file_path = os.path.join(GD_object.data_folder, 'Elevation','TanDEMX','Data','Regridded_'+str(int(GD_object.elevation_grid_posting))+'m_tiles',tile_name,file_name+'_dem_regridded_'+str(int(GD_object.elevation_grid_posting))+'m.tif')
    ds = gdal.Open(file_path)
    dem = np.array(ds.GetRasterBand(1).ReadAsArray())
    transform = ds.GetGeoTransform()
    ds = None
    x=np.arange(transform[0],transform[0]+np.shape(dem)[1]*transform[1],transform[1])
    y = np.arange(transform[3], transform[3] + np.shape(dem)[0] * transform[5], transform[5])
    if print_debug_messages:
        print('                Array size: ' + str(np.shape(dem)))
        print('                Domain of array image: '+str((np.min(x),np.max(x),np.min(y),np.max(y))))

    #subset to the glacier domain
    xIndices = np.logical_and(x>=np.min(GD_object.elevation_grid_x),x<=np.max(GD_object.elevation_grid_x))
    yIndices = np.logical_and(y >= np.min(regionYcopy), y <= np.max(regionYcopy))
    dem=dem[yIndices,:]
    dem = dem[:, xIndices]
    x=x[xIndices]
    y=y[yIndices]

    X,Y = np.meshgrid(x,y)

    points = np.hstack([np.reshape(X,(np.size(X),1)),
                        np.reshape(Y, (np.size(Y),1))])
    values = np.reshape(dem,(np.size(dem),))
    points=points[values>-10,:]
    values=values[values>-10]
    if print_debug_messages:
        print('                Shape of points: ' + str(np.shape(points)))
        print('                Shape of values: ' + str(np.shape(values)))

    #put the points into the grid
    for pp in range(len(points)):
        x_index = np.argmin(np.abs(GD_object.elevation_grid_x-points[pp,0]))
        y_index = np.argmin(np.abs(GD_object.elevation_grid_y - points[pp, 1]))
        dem_grid[y_index,x_index]+=values[pp]
        count_grid[y_index,x_index]+=1

    return(dem_grid,count_grid)

def create_dem_grid(GD_object,file_set):

    dem_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))
    count_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))

    for i in range(len(file_set)):
        tile_name = file_set[i][0]
        file_name = file_set[i][1]
        print('                  Adding data from '+file_name+' ('+tile_name+')')
        dem_grid, count_grid = read_file_to_grid(GD_object,tile_name,file_name,dem_grid,count_grid)

    dem_grid[count_grid > 0] = dem_grid[count_grid > 0] / count_grid[count_grid > 0]
    dem_grid[count_grid == 0] = -99

    return(dem_grid)

def get_tandemx_layers(GD_object,dem_files):


    unique_dates, file_sets = get_tandemx_dates_and_file_sets(GD_object,dem_files)

    dem_layers=[]
    dem_dates=[]
    dem_decYrs = []
    dem_file_strings = []

    for dd in range(len(unique_dates)):
        print('          Adding data for date '+str(unique_dates[dd])+' ('+str(dd+1)+' of '+str(len(unique_dates))+')')
        dem_date = unique_dates[dd]
        file_set = file_sets[dd]
        dem_layer = create_dem_grid(GD_object,file_set)

        if np.any(dem_layer>-99):
            dem_layers.append(dem_layer)
            dem_dates.append(dem_date)
            dem_decYr = YMD_to_DecYr(int(dem_date[:4]), int(dem_date[4:6]), int(dem_date[6:8]))
            dem_decYrs.append(dem_decYr)
            dem_file_string=''
            for file_pair in file_set:
                dem_file_string+='/'.join(file_pair)+','
            dem_file_strings.append(dem_file_string[:-1])

    return(dem_layers, dem_dates, dem_decYrs,dem_file_strings)

#######################################################################################
#These are the scripts to save the layers as a stack

def save_tandemx_layers(GD_object,dem_layers, dem_dates, dem_decYrs, dem_file_names):

    data_vars = {}
    for dd in range(len(dem_dates)):
        data_vars[dem_dates[dd]]=(['y','x'],dem_layers[dd])

    swath = xr.Dataset(data_vars,coords={'y': GD_object.elevation_grid_y,'x': GD_object.elevation_grid_x})

    for dd in range(len(dem_dates)):
        swath[dem_dates[dd]].attrs['file_name'] = dem_file_names[dd]
        swath[dem_dates[dd]].attrs['dec_yr'] = dem_decYrs[dd]

    if len(GD_object.tandemx_output_file)>2:
        output_file = GD_object.tandemx_output_file
    else:
        output_file = os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Data',GD_object.region_name+' TanDEMX Elevation Grids.nc')
    print('        Outputting data to ' + output_file)

    swath.to_netcdf(output_file)

#######################################################################################
#This is the main script to run the whole process

def generate_TanDEMX_dataset(GD_object):

    message = '    Running compilation for the TanDEMX (90m DEM) data'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)

    # step 1: get a list of files for the region
    message = '        Finding a list of TanDEMX files which overlap the domain'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    dem_files = find_tandemx_dem_files_in_domain(GD_object)

    message = '            Found ' + str(len(dem_files)) + ' files'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    if isinstance(GD_object.max_number_of_tandemx_files,int):
        if len(dem_files)>GD_object.max_number_of_tandemx_files:
            dem_files = dem_files[:GD_object.max_number_of_tandemx_files]

    # step 2: download the data and down-sample it at the requested resolution
    message = '        Downloading files and down-sampling (if not already available)'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    download_and_resample_tandemx_files(GD_object,dem_files)

    if GD_object.create_elevation_stacks:

        if GD_object.overwrite_existing_elevation_data:
            continue_to_stack = True
        else:
            if len(GD_object.tandemx_output_file) > 2:
                output_file = GD_object.tandemx_output_file
            else:
                output_file = os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation', 'Data',
                                           GD_object.region_name + ' TanDEMX Elevation Grids.nc')
            if os.path.isfile(output_file):
                continue_to_stack = False
            else:
                continue_to_stack = True

        if continue_to_stack:

            # step 3: stack the tandemx data into layers
            print('        Subsetting the resampled TanDEMX data onto to the regional domain')
            print('        Resampling the TanDEMX data onto to the regional domain')
            dem_layers, dem_dates, dem_decYrs, dem_file_strings = get_tandemx_layers(GD_object,dem_files)

            # step 4: save the tandemx layer to an nc file
            save_tandemx_layers(GD_object,dem_layers, dem_dates, dem_decYrs, dem_file_strings)

