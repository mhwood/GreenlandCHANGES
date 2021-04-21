
import xarray as xr
import numpy as np
from pyproj import Proj, transform
import os
from osgeo import ogr, osr
import shapefile
import tarfile
import matplotlib.pyplot as plt
import requests
import shutil
import netCDF4 as nc4
from ....toolbox.series import series_to_N_points
from ....toolbox.reprojection import reproject_polygon
from ....toolbox.time import YMD_to_DecYr





#######################################################################################
#These are the scripts for finding a list of Sentinel3 files

def find_sentinel3_files_from_scihub(GD_object):

    search_url = 'https://scihub.copernicus.eu/dhus/search?q=producttype:SLC&rows=1&start=0&format=json'
    

    # def check_overlap(extents1, extents2):
    #     overlap = True
    #     if extents1[2] < extents2[0] or extents1[0] > extents2[2] or extents1[3] < extents2[1] or extents1[1] > \
    #             extents2[3]:
    #         overlap = False
    #     return (overlap)
    #
    # bbox = np.array([[GD_object.extents[0], GD_object.extents[1]],
    #                  [GD_object.extents[2], GD_object.extents[1]],
    #                  [GD_object.extents[2], GD_object.extents[3]],
    #                  [GD_object.extents[0], GD_object.extents[3]],
    #                  [GD_object.extents[0], GD_object.extents[1]]])
    #
    # bbox=series_to_N_points(bbox,200)
    # bbox=reproject_polygon(bbox,3413,4326)
    #
    # region_extents = [np.min(bbox[:,0]),np.min(bbox[:,1]),np.max(bbox[:,0]),np.max(bbox[:,1])]
    #
    # dem_files = []
    #
    # for file_set in file_names:
    #     swath_extents = extent_dictionary[file_set[1]]
    #     if check_overlap(region_extents,swath_extents):
    #         dem_files.append(file_set)
    #
    # output = ''
    # for strip in dem_files:
    #     output += ','.join(strip) + '\n'
    #
    # f = open(os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation','Metadata',GD_object.region_name+ ' Sentinel3 Files.csv'), 'w')
    # f.write(output[:-1])
    # f.close()

    dem_files = []

    return (dem_files)

def find_sentinel3_dem_files_in_domain(GD_object):
    if GD_object.region_name+' Sentinel3 Files.csv' not in os.listdir(os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation','Metadata')):
        dem_files = find_sentinel3_files_from_scihub(GD_object)
    else:
        f=open(os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation','Metadata',GD_object.region_name+ ' Sentinel3 Files.csv'))
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
# These are the scripts to download all of the sentinel3 files

def download_Sentinel3_file(GD_object,datemonth):

    download_file = True

    if 'L2S_'+datemonth+'nc.tgz' in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','Sentinel3','Data')):
        download_file = False

    if download_file:
        output_file = os.path.join(GD_object.data_folder, 'Elevation','Sentinel3','Data','L2S_'+datemonth+'nc.tgz')
        message = '         Downloading file: '+'L2S_'+datemonth+'nc.tgz'
        GD_object.output_summary += '\n' + message
        if GD_object.print_main_outputs:
            print(message)

        url = sentinel3_datemonth_to_url(datemonth)

        resp=requests.get(url)
        open(output_file, 'wb').write(resp.content)

        shutil.unpack_archive(output_file[:-4])
    else:
        message = '         Already downloaded data for ' + 'L2S_' + datemonth + 'nc.tgz'
        GD_object.output_summary += '\n' + message
        if GD_object.print_main_outputs:
            print(message)

def download_Sentinel3_files(GD_object):

    # make sure dir is set up
    if 'Sentinel3' not in os.listdir(os.path.join(GD_object.data_folder, 'Elevation')):
        os.mkdir(os.path.join(GD_object.data_folder, 'Elevation','Sentinel3'))
    if 'Data' not in os.listdir(os.path.join(GD_object.data_folder, 'Elevation','Sentinel3')):
        os.mkdir(os.path.join(GD_object.data_folder, 'Elevation','Sentinel3','Data'))

    message = '            (This may take a while but is only required on the first iteration)'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    datemonths = ['201507']

    for datemonth in datemonths:
        message = '      Checking for Sentinel-3 data in date-month '+datemonth
        GD_object.output_summary += '\n' + message
        if GD_object.print_main_outputs:
            print(message)
        download_Sentinel3_file(GD_object, datemonth)



    # for date,url in file_dictionary.items():
    #     fileID=url.split('/')[6][:12]           #file is downloaded as gzip
    #     fileIDtgz=url.split('/')[6][:12]+'.tgz'
    #     dataFolder=demFolder+'Sentinel3/Data/Swath Data/'
    #
    #     if fileID not in os.listdir(dataFolder):           #checking for unzipped directory
    #
    #         if fileIDtgz not in os.listdir(dataFolder):           #checking for zipped file
    #             download_Sentinel3_files(url,dataFolder,fileID,fileIDtgz)
    #
    #         else:
    #             print('         Zipped ',fileIDtgz,'is already in Swath Data')
    #             shutil.unpack_archive(dataFolder+fileIDtgz,dataFolder+fileID)
    #             print('         Unzipped',fileIDtgz)
    #
    #
    #     else:
    #         print('         Directory',fileID,'is already in Swath Data')

def read_Sentinel3_file_extents(GD_object):

    sentinel3_folder = os.path.join(GD_object.data_folder, 'Elevation', 'Sentinel3', 'Data')

    if 'Metadata' not in os.listdir(os.path.join(GD_object.data_folder, 'Elevation', 'Sentinel3')):
        os.mkdir(os.path.join(GD_object.data_folder, 'Elevation', 'Sentinel3', 'Metadata'))

    metadata_file = 'Sentinel3 Data Extents.csv'
    if metadata_file not in os.listdir(os.path.join(GD_object.data_folder, 'Elevation', 'Sentinel3', 'Metadata')):

        message = '    Creating a bounding box metadata file for the Sentinel-3 data\n        ' \
                  '(This may take a while but is only required on the first iteration)'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

        output = 'sub_dir,file_name,min_lon,min_lat,max_lon,max_lat'

        file_list = []
        extent_dictionary = {}

        for year in range(2012, 2017):
            year_files = []
            for month in range(1, 13):
                subFolder = 'L2S_' + str(year) + '{:02d}'.format(month) + 'nc'
                if subFolder in os.listdir(sentinel3_folder):
                    subFiles = os.listdir(os.path.join(sentinel3_folder, subFolder))
                    file_counter = 0
                    for fil in subFiles:
                        file_counter += 1
                        if fil[-3:] == '.nc':
                            ds = xr.open_dataset(os.path.join(sentinel3_folder, subFolder, fil))
                            lon = np.array(ds['lon'])
                            lat = np.array(ds['lat'])
                            swath_extents = [np.min(lon) - 360, np.min(lat), np.max(lon) - 360, np.max(lat)]
                            file_list.append([subFolder,fil])
                            output += '\n'+subFolder+','+fil+','+str(np.min(lon-360))+','\
                                      +str(np.min(lat))+','+str(np.max(lon) - 360)+','+str(np.max(lat))
                            extent_dictionary[fil] = swath_extents

        f = open(os.path.join(GD_object.data_folder, 'Elevation', 'Sentinel3', 'Metadata',metadata_file),'w')
        f.write(output)
        f.close()
    else:
        file_list = []
        extent_dictionary = {}

        f = open(os.path.join(GD_object.data_folder, 'Elevation', 'Sentinel3', 'Metadata', metadata_file))
        lines = f.read()
        f.close()
        lines = lines.split('\n')
        lines.pop(0)

        for line in lines:
            line = line.split(',')
            file_list.append([line[0],line[1]])
            extent_dictionary[line[1]] = [float(line[2]),float(line[3]),float(line[4]),float(line[5])]

    return(file_list,extent_dictionary)




#######################################################################################
#These are the scripts for adding the sentinel3 points to the domain

def get_unique_dates_and_file_sets(dem_files):

    all_dates = []
    for dem_file in dem_files:
        file_name = dem_file[1]
        date = file_name.split('_')[6][:8]
        all_dates.append(date)

    unique_dates = []
    for date in all_dates:
        if date not in unique_dates:
            unique_dates.append(date)

    file_sets = []
    for date in unique_dates:
        date_set = []
        for dem_file in dem_files:
            if dem_file[1].split('_')[6][:8] == date:
                date_set.append(dem_file)
        file_sets.append(date_set)
    return(unique_dates,file_sets)

def create_dem_grid(GD_object,file_set,region_extents):

    dem_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))
    count_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))

    for df in range(len(file_set)):
        file_name = file_set[df]

        message = '              Adding file ' + str(file_name[1]) + ' (' + str(df + 1) + ' of ' + str(len(file_set)) + ')'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

        #read in the file
        ds = xr.open_dataset(os.path.join(GD_object.data_folder, 'Elevation','Sentinel3','Data', file_name[0], file_name[1]))
        lon = np.array(ds['lon'])-360
        lat = np.array(ds['lat'])
        elev = np.array(ds['elev'])
        coh = np.array(ds['coh'])
        power = np.array(ds['powerdB'])

        #find points within 0.1 degrees of the box
        lon_indices = np.logical_and(lon>=region_extents[0]-0.1, lon<=region_extents[2]+0.1)
        lat_indices = np.logical_and(lat >= region_extents[1] - 0.1, lat <= region_extents[3] + 0.1)

        #find points with coherence greater than 0.8 and power greater than -150dB
        coh_indices = np.logical_and(coh >= 0.8, coh <= 1)
        power_indices = np.logical_and(power >= -150, power <100000)

        #find points which satisfy the domain, coherence, and power tests above
        new_latlon = np.logical_and(lon_indices == True,lat_indices==True)
        new_cohpower = np.logical_and(coh_indices==True,power_indices==True)

        #select the points which satisfy the domain, coherence, and power tests above
        elev = elev[np.logical_and(new_latlon,new_cohpower)]
        lon = lon[np.logical_and(new_latlon,new_cohpower)]
        lat = lat[np.logical_and(new_latlon,new_cohpower)]

    if len(elev)>0:
        #reproject these points
        points=np.hstack([np.reshape(lon, (np.size(lon),1)),
                          np.reshape(lat, (np.size(lat), 1))])
        points=reproject_polygon(points,4326,3413)

        # for the remaining points, add them to the closest region point if the distance is less than 25m
        for p in range(len(points)):
            x_index = np.argmin(np.abs(GD_object.elevation_grid_x-points[p,0]))
            y_index = np.argmin(np.abs(GD_object.elevation_grid_y - points[p, 1]))
            if np.abs(GD_object.elevation_grid_x[x_index]-points[p,0])<=25 and np.abs(GD_object.elevation_grid_y[y_index]-points[p,1])<=25:
                dem_grid[y_index,x_index]+=elev[p]
                count_grid[y_index,x_index]+=1


    if np.any(count_grid>0):
        #divide the grid by the count grid
        dem_grid[count_grid > 0] = dem_grid[count_grid > 0] / count_grid[count_grid > 0]
        dem_grid[count_grid == 0] = -99

    else:
        dem_grid=-99*np.ones_like(dem_grid)

    return(dem_grid)


def get_sentinel3_layers(GD_object,dem_file_names):

    unique_dates, file_sets = get_unique_dates_and_file_sets(dem_file_names)

    bbox = np.array([[GD_object.extents[0], GD_object.extents[1]],
                     [GD_object.extents[2], GD_object.extents[1]],
                     [GD_object.extents[2], GD_object.extents[3]],
                     [GD_object.extents[0], GD_object.extents[3]],
                     [GD_object.extents[0], GD_object.extents[1]]])

    bbox = series_to_N_points(bbox, 200)
    bbox = reproject_polygon(bbox, 3413, 4326)

    region_extents = [np.min(bbox[:, 0]), np.min(bbox[:, 1]), np.max(bbox[:, 0]), np.max(bbox[:, 1])]

    dem_layers=[]
    dem_dates=[]
    dem_decYrs = []
    dem_file_strings = []

    for dd in range(len(unique_dates)):
        message = '          Adding data for date ' + str(unique_dates[dd]) + ' (' + str(dd + 1) + ' of ' + str(len(unique_dates)) + ')'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)
        dem_date = unique_dates[dd]
        file_set = file_sets[dd]
        dem_layer = create_dem_grid(GD_object, file_set, region_extents)

        # C = plt.imshow(dem_layer)
        # plt.colorbar(C)
        # plt.title(dem_file_names[dd])
        # plt.show()

        if np.any(dem_layer>-99):
            message = '              Found '+str(np.sum(dem_layer>-99))+' valid points in these files'
            GD_object.output_summary += '\n' + message
            if GD_object.print_sub_outputs:
                print(message)
            dem_layers.append(dem_layer)
            dem_dates.append(dem_date)

            dem_decYr = YMD_to_DecYr(int(dem_date[:4]), int(dem_date[4:6]), int(dem_date[6:8]))
            dem_decYrs.append(dem_decYr)
            dem_file_string = ''
            for file_pair in file_set:
                dem_file_string += '/'.join(file_pair) + ','
            dem_file_strings.append(dem_file_string[:-1])
        else:
            message = '              No valid points found for these files'
            GD_object.output_summary += '\n' + message
            if GD_object.print_sub_outputs:
                print(message)

    return(dem_layers, dem_dates, dem_decYrs,dem_file_strings)



#######################################################################################
#These are the scripts for stacking the glacier fields into one

def save_sentinel3_layers_as_grids(GD_object,dem_layers, dem_dates, dem_decYrs, dem_file_names):

    data_vars = {}
    for dd in range(len(dem_dates)):
        data_vars[dem_dates[dd]]=(['y','x'],dem_layers[dd])

    swath = xr.Dataset(data_vars,coords={'y': GD_object.elevation_grid_y,'x': GD_object.elevation_grid_x})

    for dd in range(len(dem_dates)):
        swath[dem_dates[dd]].attrs['file_name'] = dem_file_names[dd]
        swath[dem_dates[dd]].attrs['dec_yr'] = dem_decYrs[dd]

    if len(GD_object.sentinel3_output_file)>2:
        output_file = GD_object.sentinel3_output_file
    else:
        output_file = os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Data',GD_object.region_name+' Sentinel3 Elevation Grids.nc')
    print('        Outputting data to ' + output_file)

    swath.to_netcdf(output_file)


def save_sentinel3_layers_as_points(GD_object,dem_layers, dem_dates, dem_decYrs, dem_file_names):

    ############################################################################
    # this first part converts the points to x,y,z,xi,yi columns grids
    message = '        Converting grids to points for efficient storage'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)


    X, Y = np.meshgrid(GD_object.elevation_grid_x, GD_object.elevation_grid_y)

    if len(dem_layers)>0:
        first_grid = dem_layers[0]
        xi = np.arange(np.shape(first_grid)[1])
        yi = np.arange(np.shape(first_grid)[0])
        XI, YI = np.meshgrid(xi, yi)

        point_stacks = []
        for grid in dem_layers:
            point_stack = np.hstack([np.reshape(X, (np.size(X), 1)),
                                     np.reshape(Y, (np.size(Y), 1)),
                                     np.reshape(grid, (np.size(grid), 1)),
                                     np.reshape(XI, (np.size(XI), 1)),
                                     np.reshape(YI, (np.size(YI), 1))])
            point_stack = point_stack[point_stack[:, 2] > -99, :]
            point_stacks.append(point_stack)

    ############################################################################
    # this part saves all the point columns grids to an nc file

    if len(GD_object.sentinel3_output_file) > 2:
        output_file = GD_object.sentinel3_output_file
    else:
        output_file = os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation', 'Data',
                                   GD_object.region_name + ' Sentinel3 Elevation Points.nc')

    ds = nc4.Dataset(output_file, "w", format="NETCDF4")

    ds.createDimension("len_x", len(GD_object.elevation_grid_x))
    xvar = ds.createVariable("x", "f4", ("len_x",))
    xvar[:] = GD_object.elevation_grid_x

    ds.createDimension("len_y", len(GD_object.elevation_grid_y))
    yvar = ds.createVariable("y", "f4", ("len_y",))
    yvar[:] = GD_object.elevation_grid_y

    for dd in range(len(dem_dates)):
        grp = ds.createGroup(dem_dates[dd])
        grp.createDimension("n_points", np.shape(point_stacks[dd])[0])

        xvar = grp.createVariable("x", "f4", ("n_points",))
        yvar = grp.createVariable("y", "f4", ("n_points",))
        pointsvar = grp.createVariable("points", "f4", ("n_points",))
        xivar = grp.createVariable("xi", "f4", ("n_points",))
        yivar = grp.createVariable("yi", "f4", ("n_points",))

        xvar[:] = point_stacks[dd][:, 0]
        yvar[:] = point_stacks[dd][:, 1]
        pointsvar[:] = point_stacks[dd][:, 2]
        xivar[:] = point_stacks[dd][:, 3]
        yivar[:] = point_stacks[dd][:, 4]

        grp.file_name = dem_file_names[dd]
        grp.dec_yr = dem_decYrs[dd]

    ds.close()

    message = '        Outputting data to ' + output_file
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)


def generate_Sentinel3_dataset(GD_object):

    # step 1: download all of the sentinel3 swaths
    message = '    Running compilation for the Sentinel-3 data'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)

    message = '        Finding files that overlap the domain'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)
    find_sentinel3_dem_files_in_domain(GD_object)
    # download_Sentinel3_files(GD_object)
    #
    # message = '        Creating a bounding box metadata file for the Sentinel-3 data (to speed up iterations)'
    # GD_object.output_summary += '\n' + message
    # if GD_object.print_main_outputs:
    #     print(message)
    # file_names, extent_dictionary = read_Sentinel3_file_extents(GD_object)
    #
    # # step 2: get a list of files for the region
    # print('        Finding a list of Sentinel3 files which overlap the domain')
    # dem_file_names = find_sentinel3_dem_files_in_domain(GD_object,file_names, extent_dictionary)
    # message = '              Found '+str(len(dem_file_names))+' overlapping files out of '+str(len(file_names))+' total files'
    # GD_object.output_summary += '\n' + message
    # if GD_object.print_sub_outputs:
    #     print(message)
    #
    # if GD_object.create_elevation_stacks:
    #
    #     if GD_object.overwrite_existing_elevation_data:
    #         continue_to_stack = True
    #     else:
    #         if len(GD_object.sentinel3_output_file) > 2:
    #             output_file = GD_object.sentinel3_output_file
    #         else:
    #             output_file = os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation', 'Data',
    #                                        GD_object.region_name + ' Sentinel3 Elevation Grids.nc')
    #         if os.path.isfile(output_file):
    #             continue_to_stack = False
    #         else:
    #             continue_to_stack = True
    #
    #     if continue_to_stack:
    #
    #         # step 3: match points in the sentinel3 swaths to points in the domain
    #         print('        Sampling the Sentinel-3 data onto to the regional domain')
    #         dem_layers, dem_dates, dem_decYrs, dem_file_strings = get_sentinel3_layers(GD_object,dem_file_names)
    #
    #         print('The length of the layers is '+str(len(dem_layers)))
    #
    #         # step 4: save the sentinel3 layer to an nc file
    #         if GD_object.save_sentinel3_points_as_grids:
    #             save_sentinel3_layers_as_grids(GD_object,dem_layers, dem_dates, dem_decYrs, dem_file_strings)
    #
    #         if GD_object.save_sentinel3_points_as_points:
    #             save_sentinel3_layers_as_points(GD_object,dem_layers, dem_dates, dem_decYrs, dem_file_strings)
    #     else:
    #         message = '        Sentinel-3 Elevation Stack has already been created'
    #         GD_object.output_summary += '\n' + message
    #         if GD_object.print_main_outputs:
    #             print(message)