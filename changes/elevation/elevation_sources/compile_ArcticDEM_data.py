
from datetime import datetime
import numpy as np
from osgeo import ogr
import shapefile
import requests
import urllib
import os
import tarfile
from osgeo import gdal
from osgeo import osr
import xarray as xr
import shutil
from ....toolbox.time import YMD_to_DecYr

#######################################################################################
#These are the scripts for finding a list of ArcticDEM files

def find_arcticdem_dem_files_from_shapefile(GD_object):

    def stripOutlineToWKT(stripOutline):
        output = "POLYGON (("
        for point in stripOutline:
            output += str(point[0]) + " " + str(point[1]) + ", "
        output = output[:-2] + "))"
        return (output)

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
    bboxWKT = bboxToWKT(bbox)
    poly1 = ogr.CreateGeometryFromWkt(bboxWKT)

    shapefile_name = 'ArcticDEM_Strip_Index_Rel7'

    ####################################################################################################################
    # Download the metadata file if its not already done
    if shapefile_name not in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','Metadata')):
        message = '            Downloading the ArcticDEM metadata shapefile'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

        download_file = os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','Metadata',shapefile_name+'.zip')
        resp = requests.get('http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/ArcticDEM_Strip_Index_Rel7.zip')
        open(download_file, 'wb').write(resp.content)

        shutil.unpack_archive(download_file,
                              os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','Metadata',shapefile_name))

    ####################################################################################################################
    # Loop through the metadata file to see which files pertain to the sample area

    message = '            Searching through shapefile provided by PGC to find overlapping files'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)
        message = '                (This may take a minute or two)'
        print(message)

    r = shapefile.Reader(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','Metadata',shapefile_name,shapefile_name+'.shp'))
    shapes = r.shapes()
    records = r.records()
    dem_files = []

    for s in range(len(shapes)):
        stripOutline = shapes[s].points
        stripWKT = stripOutlineToWKT(stripOutline)   #convert stripOutline (object) to WKT formatted
        poly2 = ogr.CreateGeometryFromWkt(stripWKT)
        intersection = poly1.Intersection(poly2)
        intersection = intersection.ExportToWkt()
        if "EMPTY" not in intersection:
            filePath = records[s][9]
            tile = filePath.split('/')[10]
            fileID = filePath.split('/')[-1][:-7]
            if '.edu' in filePath:
                try:
                    resp = urllib.request.urlopen(filePath)   #ArcticDEM_Strip_Index_Rel7 has errors
                except urllib.error.HTTPError:                #This filters out files which do not exist
                    pass
                else:
                    dem_files.append([tile, fileID])

    output = ''
    for strip in dem_files:
        output += ','.join(strip) + '\n'

    f = open(os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation','Metadata',GD_object.region_name + ' ArcticDEM Files.csv'), 'w')
    f.write(output)
    f.close()
    return (dem_files)

def find_arcticdem_dem_files_in_domain(GD_object):

    if GD_object.region_name + ' ArcticDEM Files.csv' not in os.listdir(
            os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation','Metadata')):

        dem_files = find_arcticdem_dem_files_from_shapefile(GD_object)

    else:
        message = '            Reading file list generated in a previous iteration of this process'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

        f=open(os.path.join(GD_object.project_folder, GD_object.region_name,
                            'Elevation','Metadata',GD_object.region_name + ' ArcticDEM Files.csv'))
        lines=f.read()
        f.close()
        lines=lines.split('\n')
        dem_files=[]
        for line in lines:
            line=line.split(',')
            if len(line)>1:
                dem_files.append(line)

    output_dem_files = []

    for df in range(len(dem_files)):
        dem_date_str = dem_files[df][1].split('_')[2][:8]
        date_test = datetime(int(dem_date_str[:4]), int(dem_date_str[4:6]), int(dem_date_str[6:8]))
        if date_test >= GD_object.date_1 and date_test < GD_object.date_2:
            output_dem_files.append(dem_files[df])

    return(output_dem_files)

#######################################################################################
#These are the scripts for downloading the list of arctic dem files

def check_file_download(GD_object,file_name,tile_name):
    download_file = True

    # check if the raw data is there
    raw_data_present = False
    if tile_name in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles')):
        if file_name + '.tar.gz' in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles',tile_name)):
            raw_data_present = True
        if file_name + '_dem.tif' in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles',tile_name)):
            raw_data_present = True

    # check if the resampled data is there
    resampled_data_present = False
    if tile_name in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','Regridded_'+str(int(GD_object.elevation_grid_posting))+'m_tiles')):
        if file_name + '_dem_regridded_' + str(int(GD_object.elevation_grid_posting)) + 'm.tif' in os.listdir(
                os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','Regridded_'+str(int(GD_object.elevation_grid_posting))+'m_tiles',tile_name)):
            resampled_data_present = True

    # if neither of these files is present, then download the data
    if raw_data_present or resampled_data_present:
        download_file = False

    return(download_file,resampled_data_present)

def download_arcticdem_tile(GD_object, tile_name, file_name):

    if tile_name not in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles')):
        os.mkdir(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles',tile_name))

    download_url = 'http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/geocell/v3.0/2m/' + tile_name + '/' + file_name + '.tar.gz'
    output_file = os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles',tile_name,file_name + '.tar.gz')

    resp = requests.get(download_url)
    open(output_file, 'wb').write(resp.content)

def resample_arcticdem_tile_nearest_neighbor(GD_object, tile_name, file_name):

    #untar the file if necessary
    if file_name+'_dem.tif' not in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles',tile_name)):
        print('                Untaring the data')
        tar = tarfile.open(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles',tile_name,file_name + '.tar.gz'))
        tar.extract(file_name + '_dem.tif', path=os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles',tile_name))
        tar.close()

    print('                Working on the down sample')
    print('                  Reading in the file')
    #read in the raster
    ds = gdal.Open(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles',tile_name,file_name + '_dem.tif'))
    dem = np.array(ds.GetRasterBand(1).ReadAsArray())
    transform = ds.GetGeoTransform()
    ds = None


    demX = np.arange(transform[0], transform[0] + transform[1] * np.shape(dem)[1], transform[1])
    demY = np.arange(transform[3], transform[3] + transform[5] * np.shape(dem)[0], transform[5])

    #make a new grid for the raster
    rounding_level = int(np.ceil(np.log10(GD_object.elevation_grid_posting)))
    min_X = np.round(np.min(demX),-rounding_level)
    if min_X> np.min(demX):
        min_X-=10**rounding_level
    min_Y = np.round(np.min(demY), -rounding_level)
    if min_Y > np.min(demY):
        min_Y -= 10 ** rounding_level
    max_X = np.round(np.max(demX), -rounding_level)
    if max_X < np.max(demX):
        min_X += 10 ** rounding_level
    max_Y = np.round(np.max(demY), -rounding_level)
    if max_Y < np.max(demY):
        min_X += 10 ** rounding_level
    new_dem_X = np.arange(min_X, max_X + GD_object.elevation_grid_posting, GD_object.elevation_grid_posting)
    new_dem_Y = np.arange(min_Y, max_Y + GD_object.elevation_grid_posting, GD_object.elevation_grid_posting)

    print('                  Finding the nearest neighbors')
    # resample the raster onto the new grid
    #this is done via nearest neighbor - could be improved in future iterations
    newDem = -9999*np.ones((len(new_dem_Y), len(new_dem_X)))
    for i in range(len(new_dem_X)):
        demXindex = np.argmin(np.abs(new_dem_X[i] - demX))
        for j in range(len(new_dem_Y)):
            demYindex = np.argmin(np.abs(new_dem_Y[j] - demY))
            if dem[demYindex, demXindex]>-10:
                if ((new_dem_X[i]-demX[demXindex])**2 + (new_dem_Y[j]-demY[demYindex])**2)**0.5<GD_object.elevation_grid_posting:
                    newDem[j, i] = dem[demYindex, demXindex]

    #save the raster as a tif
    if tile_name not in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','Regridded_'+str(int(GD_object.elevation_grid_posting))+'m_tiles')):
        os.mkdir(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','Regridded_'+str(int(GD_object.elevation_grid_posting))+'m_tiles',tile_name))

    outputFile = os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','Regridded_'+str(int(GD_object.elevation_grid_posting))+'m_tiles',tile_name,file_name+'_dem_regridded_'+str(int(GD_object.elevation_grid_posting))+'m.tif')
    geotransform = (np.min(new_dem_X), GD_object.elevation_grid_posting, 0, np.max(new_dem_Y), 0, -GD_object.elevation_grid_posting)
    output_raster = gdal.GetDriverByName('GTiff').Create(outputFile, len(new_dem_X), len(new_dem_Y), 1, gdal.GDT_Float32)
    output_raster.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3413)
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(np.flipud(newDem))

def download_and_resample_arcticdem_files(GD_object,dem_files):

    for dd in range(len(dem_files)):
        dem_file = dem_files[dd]
        tile_name = dem_file[0]
        file_name = dem_file[1]
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
            download_arcticdem_tile(GD_object, tile_name, file_name)
        else:
            if GD_object.print_sub_outputs:
                print('              File already obtained')

        #resample the data if its not already done
        if GD_object.resample_high_resolution_arcticDEM_data:
            if not resampled_data_present:
                if GD_object.print_sub_outputs:
                    print('              Downsampling file... ')
                resample_arcticdem_tile_nearest_neighbor(GD_object, tile_name, file_name)
                resampled_data_present=True
            else:
                if GD_object.print_sub_outputs:
                    print('              File already resampled')

        if not GD_object.keep_high_resolution_arcticDEM_data:
            # remove the tar files and raw data after resampling because these files are large
            if resampled_data_present:
                if file_name + '.tar.gz' in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles',tile_name)):
                    os.remove(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles',tile_name,file_name + '.tar.gz'))
                if file_name + '_dem.tif' in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles',tile_name)):
                    os.remove(os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','2m_tiles',tile_name,file_name + '_dem.tif'))

#######################################################################################
#These are the scripts for creating the layers in the ArcticDEM data

def get_unique_dates_and_file_sets(dem_files):

    all_dates = []
    for dem_file in dem_files:
        file_name = dem_file[1]
        date = file_name.split('_')[2]
        all_dates.append(date)

    unique_dates = []
    for date in all_dates:
        if date not in unique_dates:
            unique_dates.append(date)

    file_sets = []
    for date in unique_dates:
        date_set = []
        for dem_file in dem_files:
            if dem_file[1].split('_')[2] == date:
                date_set.append(dem_file)
        file_sets.append(date_set)
    return(unique_dates,file_sets)

def read_file_to_grid(GD_object,tile_name,file_name,dem_grid,count_grid):

    print_debug_messages = False

    regionYcopy=np.flipud(GD_object.elevation_grid_y)

    #read in the down sampled grid

    file_path = os.path.join(GD_object.data_folder, 'Elevation','ArcticDEM','Regridded_'+str(int(GD_object.elevation_grid_posting))+'m_tiles',tile_name,file_name+'_dem_regridded_'+str(int(GD_object.elevation_grid_posting))+'m.tif')
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

def get_arcticdem_layers(GD_object,dem_files):

    unique_dates, file_sets = get_unique_dates_and_file_sets(dem_files)

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

def save_arcticdem_layers(GD_object,dem_layers, dem_dates, dem_decYrs, dem_file_names):

    data_vars = {}
    for dd in range(len(dem_dates)):
        data_vars[dem_dates[dd]]=(['y','x'],dem_layers[dd])

    swath = xr.Dataset(data_vars,coords={'y': GD_object.elevation_grid_y,'x': GD_object.elevation_grid_x})

    for dd in range(len(dem_dates)):
        swath[dem_dates[dd]].attrs['file_name'] = dem_file_names[dd]
        swath[dem_dates[dd]].attrs['dec_yr'] = dem_decYrs[dd]

    if len(GD_object.arcticdem_output_file)>2:
        output_file = GD_object.arcticdem_output_file
    else:
        output_file = os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Data',GD_object.region_name+' ArcticDEM Elevation Grids.nc')
    print('        Outputting data to ' + output_file)

    swath.to_netcdf(output_file)

#######################################################################################
#This is the main script to run the whole process

def generate_ArcticDEM_dataset(GD_object):

    message = '    Running compilation for the ArcticDEM (Worldview) data'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)

    # step 1: get a list of files for the region
    message = '        Finding a list of ArcticDEM files which overlap the domain'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    dem_files = find_arcticdem_dem_files_in_domain(GD_object)

    message = '            Found ' + str(len(dem_files)) + ' files'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    if isinstance(GD_object.max_number_of_arcticDEM_files,int):
        if len(dem_files)>GD_object.max_number_of_arcticDEM_files:
            dem_files = dem_files[:GD_object.max_number_of_arcticDEM_files]

    # step 2: download the data and down-sample it at the requested resolution
    message = '        Downloading files and down-sampling (if not already available)'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)
    download_and_resample_arcticdem_files(GD_object,dem_files)

    if GD_object.create_elevation_stacks:



        if GD_object.overwrite_existing_elevation_data:
            continue_to_stack = True
        else:
            if len(GD_object.arcticdem_output_file) > 2:
                output_file = GD_object.arcticdem_output_file
            else:
                output_file = os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation', 'Data',
                                           GD_object.region_name + ' ArcticDEM Elevation Grids.nc')
            if os.path.isfile(output_file):
                continue_to_stack = False
            else:
                continue_to_stack = True

        if continue_to_stack:

            # step 3: stack the arcticdem data into layers
            print('        Subsetting the resampled ArcticDEM data onto to the regional domain')
            print('        Resampling the ArcticDEM data onto to the regional domain')
            dem_layers, dem_dates, dem_decYrs, dem_file_strings = get_arcticdem_layers(GD_object,dem_files)

            # step 4: save the arcticdem layer to an nc file
            save_arcticdem_layers(GD_object,dem_layers, dem_dates, dem_decYrs, dem_file_strings)

