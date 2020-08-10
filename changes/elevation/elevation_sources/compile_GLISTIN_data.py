

import xarray as xr
import numpy as np
import os
import requests
from osgeo import ogr
from toolbox.reprojection import reproject_polygon
from toolbox.time import YMD_to_DecYr

#######################################################################################
# These files are to find the files in the domain

def find_glistin_dem_files_from_shapefile(GD_object):

    def stripOutlineToWKT(stripOutline): #converts strip outline to well known text ?? format
        output = "POLYGON (("
        for point in stripOutline:
            output += str(point[0]) + " " + str(point[1]) + ", "
        output = output[:-2] + "))"
        return (output)

    def bboxToWKT(bbox): #this takes the numpy array "bbox" and turns it into a wkt so ogr can make a polygon from the values
        output = "POLYGON (("
        for b in range(np.shape(bbox)[0]):
            output += str(bbox[b, 0]) + " " + str(bbox[b, 1]) + ", "
        output = output[:-2] + "))"
        #print('output of bboxToWKT:::', output)
        return (output)

    bbox = np.array([[np.min(GD_object.elevation_grid_x), np.min(GD_object.elevation_grid_y)],
                     [np.max(GD_object.elevation_grid_x), np.min(GD_object.elevation_grid_y)],
                     [np.max(GD_object.elevation_grid_x), np.max(GD_object.elevation_grid_y)],
                     [np.min(GD_object.elevation_grid_x), np.max(GD_object.elevation_grid_y)],
                     [np.min(GD_object.elevation_grid_x), np.min(GD_object.elevation_grid_y)]])
    bboxWKT = bboxToWKT(bbox)
    poly1 = ogr.CreateGeometryFromWkt(bboxWKT)

    import changes.reference.glistin_domains as gld
    glistin_ids = list(gld.glistin_id_to_polygon_outline.keys())

    output_ids = []
    for pr in range(len(glistin_ids)):
        glistin_id = glistin_ids[pr]
        glistin_id_polygon = np.array(gld.glistin_id_to_polygon_outline[glistin_id])
        glistinWKT = stripOutlineToWKT(glistin_id_polygon)
        poly2 = ogr.CreateGeometryFromWkt(glistinWKT)
        intersection = poly1.Intersection(poly2)
        intersection = intersection.ExportToWkt()
        if "EMPTY" not in intersection:
            output_ids.append(glistin_id)

    output = ''
    for strip in output_ids:
        output += ''.join(strip)+'.nc'+ '\n'

    f = open(os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Metadata',GD_object.region_name + ' GLISTIN Files.csv'), 'w')
    f.write(output)
    f.close()

    return (output_ids)


def find_glistin_dem_files_in_domain(GD_object):

    if GD_object.region_name + ' GLISTIN Files.csv' not in os.listdir(os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Metadata')):
        dem_files = find_glistin_dem_files_from_shapefile(GD_object)
    else:
        f=open(os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Metadata',GD_object.region_name + ' GLISTIN Files.csv'))
        lines=f.read()
        f.close()
        lines=lines.split('\n')
        dem_files = []
        for line in lines:
            if len(line)>1:
                dem_files.append(line)

    dem_dates= []
    for s in dem_files:
        date = s.split('_')[4][:8]
        dem_dates.append(date)


    return(dem_files, dem_dates)





#######################################################################################
#These are the scripts for downloading the GLISTIN data

def downloadStrip(GD_object,dataFolder,year,dem_file): #
    if year not in os.listdir(dataFolder):
        os.mkdir(os.path.join(dataFolder,year))
    url ='https://podaac-opendap.jpl.nasa.gov/opendap/allData/omg/L3/elevation/GLISTINA/'+year+'/'+dem_file
    if GD_object.print_sub_outputs:
        print('              URL: '+url) #
    outputFile=os.path.join(dataFolder,year,dem_file)

    with requests.get(url,auth=(GD_object.podaac_username,GD_object.podaac_password), stream=True) as r:
        r.raise_for_status()
        with open(outputFile, 'wb') as f:                   #open local output file to write binary to, name returned filehandle f
            for chunk in r.iter_content(chunk_size=8192):   #for each 8192byte chunk in the incoming stream
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                # if chunk:
                f.write(chunk)                              #write that 8192byte chunk to the filehandle f #

def download_glistin_files(GD_object,dem_files):

    if GD_object.podaac_username == '' or GD_object.podaac_password == '':
        print('Please identify your PO.DAAC username and password with the following commands in the init.py file:'+
                   '\n    GC.podaac_username = *********'+
                   '\n    GC.podaac_password = *********'+
                   '\n    Note: these credentials are NOT the same as your EARTHDATA credentials.'+
                   '\n    To find your username and password, login to EARTHDATA at https://urs.earthdata.nasa.gov/'+
                   '\n    using your EARTHDATA credentials. Then, navigate to https://podaac-tools.jpl.nasa.gov/drive/'+
                   '\n    to find your PODAAC credentials.')
        continue_to_stack = False

    else:

        dataFolder = os.path.join(GD_object.data_folder, 'Elevation', 'GLISTIN','Data','L3')
        for dd in range(len(dem_files)):
            dem_file = dem_files[dd]
            if GD_object.print_sub_outputs:
                print('            Checking file '+dem_file)
            year = dem_file.split('_')[4][:4] #
            download_file=True
            if year in os.listdir(dataFolder):
                if dem_file in os.listdir(os.path.join(dataFolder,year)):
                    download_file=False
            if download_file:
                if GD_object.print_sub_outputs:
                    print('              Downloading file ' + dem_file)
                downloadStrip(GD_object,dataFolder, year, dem_file)

        continue_to_stack = True
    return(continue_to_stack)

#######################################################################################
#These are the scripts for creating the layers in the GLISTIN data


def add_glistin_dem_field_to_grid(GD_object,dem_file,dem_date,dem_grid,count_grid):
    if GD_object.print_sub_outputs:
        print('                  Reading in the DEM')
    #read in the data
    subFolder=dem_date[:4]
    dem_file_path = os.path.join(GD_object.data_folder,'Elevation','GLISTIN','Data','L3',subFolder+'/'+dem_file)
    ds = xr.open_dataset(dem_file_path)
    x=np.array(ds['x'])
    y=np.array(ds['y'])
    dem_data=np.array(ds['elevation'])
    EPSG= ds['projection'].attrs['srid']
    EPSG=int(EPSG[-5:])

    if GD_object.print_sub_outputs:
        print('                  Limiting DEM to region domain')

    X,Y=np.meshgrid(x,y)

    points = np.hstack([np.reshape(X,(np.size(X),1)),
                        np.reshape(Y, (np.size(Y), 1))])
    dem_data = np.reshape(dem_data,(np.size(dem_data),1))

    points=points[~np.isnan(dem_data[:,0]),:]
    dem_data=dem_data[~np.isnan(dem_data[:,0])]
    points= reproject_polygon(points,EPSG,3413)

    good_x_indices = np.logical_and(points[:,0]>=np.min(GD_object.elevation_grid_x),points[:,0]<=np.max(GD_object.elevation_grid_x))
    dem_data = dem_data[good_x_indices]
    points = points[good_x_indices,:]
    good_y_indices = np.logical_and(points[:, 1] >= np.min(GD_object.elevation_grid_y), points[:, 1] <= np.max(GD_object.elevation_grid_y))
    dem_data = dem_data[good_y_indices]
    points = points[good_y_indices, :]

    if GD_object.print_sub_outputs:
        print('                  Adding points to domain via nearest neighbor interpolation')
    #fill in the glacier field on their closest points
    total_points_added=0
    for pp in range(np.shape(points)[0]):
        # if pp%10000==0:
        #     print('                      Point = '+str(pp)+' of '+str(np.shape(points)[0])+',  Total Points Added = '+str(total_points_added))
        x_index = np.argmin(np.abs(GD_object.elevation_grid_x-points[pp,0]))
        y_index = np.argmin(np.abs(GD_object.elevation_grid_y-points[pp,1]))
        if ((GD_object.elevation_grid_x[x_index]-points[pp,0])**2 + (GD_object.elevation_grid_y[y_index]-points[pp,1])**2)**0.5 < (50.0/2)*np.sqrt(2):
            dem_grid[y_index,x_index]+=dem_data[pp]
            count_grid[y_index,x_index]+=1
            total_points_added+=1

    return(dem_grid,count_grid)

def get_glistin_layers(GD_object,dem_file_names, dem_dates):

    dem_dates_unique = []
    for date in dem_dates:
        if date not in dem_dates_unique:
            dem_dates_unique.append(date)

    dem_layers=[]
    dem_decYrs = []
    for dd in range(len(dem_dates_unique)):
        dem_date_unique = dem_dates_unique[dd]
        if GD_object.print_sub_outputs:
            print('          Working on data for date ' + dem_date_unique + ' (' + str(dd + 1) + ' of ' + str(len(dem_dates_unique)) + ')')

        dem_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))
        count_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))
        for di in range(len(dem_dates)):
            if dem_dates[di]==dem_date_unique:
                dem_file = dem_file_names[di]
                if GD_object.print_sub_outputs:
                    print('              Adding in data from file '+dem_file+' ('+str(di+1)+' of '+str(len(dem_file_names))+')')
                dem_grid,count_grid= add_glistin_dem_field_to_grid(GD_object,dem_file,dem_date_unique,dem_grid,count_grid)

        dem_grid[count_grid > 0] = dem_grid[count_grid > 0] / count_grid[count_grid > 0]
        dem_grid[count_grid==0]=-99
        dem_layers.append(dem_grid)
        dem_decYr = YMD_to_DecYr(int(dem_date_unique[:4]), int(dem_date_unique[4:6]), int(dem_date_unique[6:8]))
        dem_decYrs.append(dem_decYr)

    return(dem_layers, dem_dates_unique, dem_decYrs)


#######################################################################################
#These are the scripts for stacking the glacier fields into one

def save_glistin_layers(GD_object, dem_layers, dem_dates_unique, dem_decYrs, dem_file_names):

    data_vars = {}
    for dd in range(len(dem_dates_unique)):
        data_vars[dem_dates_unique[dd]]=(['y','x'],dem_layers[dd])

    swath = xr.Dataset(data_vars,coords={'y': GD_object.elevation_grid_y,'x': GD_object.elevation_grid_x})

    for dd in range(len(dem_dates_unique)):
        swath[dem_dates_unique[dd]].attrs['file_name'] = dem_file_names[dd]
        swath[dem_dates_unique[dd]].attrs['dec_yr'] = dem_decYrs[dd]

    output_file = os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Data',GD_object.region_name+' GLISTIN Elevation Grids.nc')

    message = '        Outputting data to ' + output_file
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    swath.to_netcdf(output_file)


def generate_glistin_dataset(GD_object):

    message = '    Running compilation for the GLISTIN data'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)

    # step 1: get a list of files for the region
    message = '        Finding a list of GLISTIN files which overlap the domain'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    dem_file_names, dem_dates = find_glistin_dem_files_in_domain(GD_object)

    if isinstance(GD_object.max_number_of_glistin_files,int):
        if len(dem_file_names)>GD_object.max_number_of_glistin_files:
            dem_file_names = dem_file_names[:GD_object.max_number_of_glistin_files]
            dem_dates = dem_dates[:GD_object.max_number_of_glistin_files]

    message = '        Found ' + str(len(dem_file_names)) + ' files'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)


    #     #step 2: download files
    message = '        Downloading the GLISTIN files from PO.DAAC'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    continue_to_stack = download_glistin_files(GD_object,dem_file_names)

    if continue_to_stack and GD_object.create_elevation_stacks:

        message = '        Resampling the GLISTIN data onto to the regional domain'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

        # step 3: stack the glistin data into layers
        dem_layers, dem_dates_unique, dem_decYrs = get_glistin_layers(GD_object,dem_file_names, dem_dates)

        # step 4: save the glistin layer to an nc file
        save_glistin_layers(GD_object,dem_layers, dem_dates_unique, dem_decYrs, dem_file_names)