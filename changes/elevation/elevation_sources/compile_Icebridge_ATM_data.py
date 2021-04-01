

import os
import csv
import xarray as xr
import numpy as np
import itertools
import requests
from requests.auth import HTTPBasicAuth
import time
import netCDF4 as nc4
from ....toolbox.time import YMD_to_DecYr
from ....toolbox.reprojection import reproject_polygon
from ....toolbox.series import series_to_N_points
import matplotlib.pyplot as plt

#######################################################################################
#These are the scripts for finding a list of icebridge_atm files and their associated links for downloading
def find_icebridge_atm_dem_files_from_nsidc_portal(GD_object):

    # this will make a lat/long string of the glacier domain
    def getPolygonString(regionX,regionY):

        x = regionX
        y = regionY

        bbox = np.array([[np.min(x), np.min(y)],
                         [np.max(x), np.min(y)],
                         [np.max(x), np.max(y)],
                         [np.min(x), np.max(y)],
                         [np.min(x), np.min(y)]])

        nSidePoints = 20
        bbox = series_to_N_points(bbox, 4 * nSidePoints + 4)

        print(bbox)

        bbox = reproject_polygon(bbox, 3413, 4326)

        print(bbox)

        polygon = ''
        for bb in range(np.shape(bbox)[0]):
            polygon += str(bbox[bb, 0]) + ',' + str(bbox[bb, 1]) + ','

        polygon = polygon[:-1]

        return (polygon)

    # this function is from the icebridge_atm data page
    def build_version_query_params(version):
        desired_pad_length = 3
        if len(version) > desired_pad_length:
            print('Version string too long: "{0}"'.format(version))
            quit()

        version = str(int(version))  # Strip off any leading zeros
        query_params = ''

        while len(version) <= desired_pad_length:
            padded_version = version.zfill(desired_pad_length)
            query_params += '&version={0}'.format(padded_version)
            desired_pad_length -= 1
        return query_params

    # this function is from the icebridge_atm data page
    def build_cmr_query_url(CMR_FILE_URL, short_name, version, time_start, time_end, polygon=None,
                            filename_filter=None):
        params = '&short_name={0}'.format(short_name)
        params += build_version_query_params(version)
        params += '&temporal[]={0},{1}'.format(time_start, time_end)

        if polygon:
            params += '&polygon={0}'.format(polygon)
        if filename_filter:
            params += '&producer_granule_id[]={0}&options[producer_granule_id][pattern]=true'.format(filename_filter)

        return CMR_FILE_URL + params

    # this function is from the icebridge_atm data page
    def cmr_filter_urls(search_results):
        """Select only the desired data files from CMR response."""
        if 'feed' not in search_results or 'entry' not in search_results['feed']:
            return []

        entries = [e['links']
                   for e in search_results['feed']['entry']
                   if 'links' in e]
        # Flatten "entries" to a simple list of links
        links = list(itertools.chain(*entries))

        urls = []
        unique_filenames = set()
        for link in links:
            if 'href' not in link:
                # Exclude links with nothing to download
                continue
            if 'inherited' in link and link['inherited'] is True:
                # Why are we excluding these links?
                continue
            if 'rel' in link and 'data#' not in link['rel']:
                # Exclude links which are not classified by CMR as "data" or "metadata"
                continue

            if 'title' in link and 'opendap' in link['title'].lower():
                # Exclude OPeNDAP links--they are responsible for many duplicates
                # This is a hack; when the metadata is updated to properly identify
                # non-datapool links, we should be able to do this in a non-hack way
                continue

            filename = link['href'].split('/')[-1]
            if filename in unique_filenames:
                # Exclude links with duplicate filenames (they would overwrite)
                continue
            unique_filenames.add(filename)

            urls.append(link['href'])

        return urls

    ####################################################################################################################
    # This part is for the ATM data
    short_name = 'ILATM2'
    version = '2'

    time_start = '1900-03-31T00:00:00Z'
    time_end = '2100-05-16T23:59:59Z'

    polygon = getPolygonString(GD_object.elevation_grid_x,GD_object.elevation_grid_y)
    filename_filter = '*'

    CMR_URL = 'https://cmr.earthdata.nasa.gov'
    URS_URL = 'https://urs.earthdata.nasa.gov'
    CMR_PAGE_SIZE = 2000
    CMR_FILE_URL = ('{0}/search/granules.json?provider=NSIDC_ECS'
                    '&sort_key[]=start_date&sort_key[]=producer_granule_id'
                    '&scroll=true&page_size={1}'.format(CMR_URL, CMR_PAGE_SIZE))


    cmr_query_url = build_cmr_query_url(CMR_FILE_URL, short_name=short_name, version=version,
                                        time_start=time_start, time_end=time_end,
                                        polygon=polygon, filename_filter=filename_filter)

    #Uses .netrc in home directory
    response = requests.get(cmr_query_url)
    search_page = response.json()
    urls = cmr_filter_urls(search_page)

    outputUrls = []
    outputDates = []
    outputFiles = []
    outputSensors = []

    for url in urls:
        if url[-4:] == '.csv':
            urlParts = url.split('/')
            fileID = urlParts[-1][:-4]
            outputFiles.append(fileID+'.csv')
            date = fileID.split('_')[1][:8]
            outputDates.append(date)
            outputUrls.append(url)
            outputSensors.append(short_name)



    textOutput = ''
    for ss in range(len(outputUrls)):
        textOutput += outputSensors[ss]+','+outputDates[ss]+','+outputFiles[ss]+','+outputUrls[ss] + '\n'

    f = open(os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation','Metadata',GD_object.region_name + ' IceBridge ATM Files.csv'),'w')
    f.write(textOutput[:-1])
    f.close()

    return (outputSensors,outputFiles,outputDates,outputUrls)

def find_icebridge_atm_dem_files_in_domain(GD_object):
    if GD_object.region_name + ' IceBridge ATM Files.csv' not in os.listdir(os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation','Metadata')):
        dem_file_sensors, dem_files, dem_file_dates, dem_file_links = find_icebridge_atm_dem_files_from_nsidc_portal(GD_object)
    else:
        f=open(os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation','Metadata',GD_object.region_name + ' IceBridge ATM Files.csv'))
        lines=f.read()
        f.close()
        lines=lines.split('\n')
        dem_files=[]
        dem_file_links=[]
        dem_file_dates=[]
        dem_file_sensors=[]
        for line in lines:
            line=line.split(',')
            if len(line)>1:
                dem_file_sensors.append(line[0])
                dem_file_dates.append(line[1])
                dem_files.append(line[2])
                dem_file_links.append(line[3])
    return(dem_files,dem_file_dates,dem_file_sensors,dem_file_links)


#######################################################################################
#These are the scripts for downloading the icebridge_atm data

def downloadStrip(stripFolder,sensor,date,file_name,url):

    outputFile = os.path.join(stripFolder,sensor,date,file_name)

    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(max_retries=20)
    session.mount('https://', adapter)
    session.mount('http://', adapter)
    resp = session.get(url)

    # resp = requests.get(url)
    open(outputFile, 'wb').write(resp.content)

def download_icebridge_atm_files(GD_object,dem_file_names, dem_file_dates, dem_file_sensors, dem_file_links):

    dataFolder = os.path.join(GD_object.data_folder,'Elevation','IceBridge','Data')

    for dd in range(len(dem_file_names)):
        dem_file_name = dem_file_names[dd]
        dem_date = dem_file_dates[dd]
        dem_file_sensor = dem_file_sensors[dd]
        if GD_object.print_sub_outputs:
            print('            Checking file '+dem_file_name)

        download_file=True

        if dem_file_sensor not in os.listdir(dataFolder):
            os.mkdir(os.path.join(dataFolder, dem_file_sensor))

        if dem_date not in os.listdir(os.path.join(dataFolder, dem_file_sensor)):
            os.mkdir(os.path.join(dataFolder, dem_file_sensor,dem_date))

        if dem_file_name in os.listdir(os.path.join(dataFolder, dem_file_sensor,dem_date)):
            download_file=False

        if download_file:
            if GD_object.print_sub_outputs:
                print('              Downloading file ' + dem_file_name+' ('+str(dd+1)+' of '+str(len(dem_file_names))+')')
            downloadStrip(dataFolder, dem_file_sensor,dem_date, dem_file_name, dem_file_links[dd])
        else:
            if GD_object.print_sub_outputs:
                print('              File already downloaded')


#######################################################################################
#These are the scripts for creating the layers in the icebridge_atm data


def read_valid_points_from_csv(filepath,yrmd):
    csv_file=open(filepath)
    read_csv_file=csv.reader(csv_file,delimiter=',')

    if int(yrmd) <= 20190000:
        for dd in range(10):   #before 2019, the csv files on nsidc have 10 header rows
            next(read_csv_file)
    else:
        for dd in range(11):   #after 2019, csv files have 11 header rows, one blank
            next(read_csv_file)

    lon_lat_elev=[]
    for row in read_csv_file:
        lat=float(row[1])
        lon=float(row[2])-360.0
        elev=float(row[3])
        lon_lat_elev.append([lon,lat,elev])

    lon_lat_elev=np.array(lon_lat_elev)

    return(lon_lat_elev)

def cumulate_icebridge_atm_dem_field(valid_points,regionX,regionY,dem_grid,count_grid):

    #limit the points to the region domain
    dem_data = valid_points
    dem_data = dem_data[np.logical_and(dem_data[:, 0] <= np.max(regionX), dem_data[:, 0] >= np.min(regionX)), :]
    dem_data = dem_data[np.logical_and(dem_data[:, 1] <= np.max(regionY), dem_data[:, 1] >= np.min(regionY))]


    #fill in the glacier field by averaging points based on their closest point
    for dd in range(np.shape(dem_data)[0]):
        x_index = np.argmin(np.abs(regionX-dem_data[dd,0]))
        y_index = np.argmin(np.abs(regionY-dem_data[dd,1]))
        dem_grid[y_index,x_index]+=dem_data[dd,2]
        count_grid[y_index,x_index]+=1

    return(dem_grid,count_grid)

def get_icebridge_atm_layers(GD_object,dem_file_names, dem_file_dates, dem_file_sensors):

    # get a list of unique dates and their corresponding files so that the measurements on the same day are in the same grid
    dem_dates = []
    dem_sensors = []
    dem_decYrs = []
    date_files = []
    filename_strings = []
    for dd in range(len(dem_file_names)):
        dem_file = dem_file_names[dd]
        dem_date = dem_file_dates[dd]
        if dem_date not in dem_dates:
            dem_dates.append(dem_date)
            dem_decYr = YMD_to_DecYr(int(dem_date[:4]), int(dem_date[4:6]), int(dem_date[6:8]))
            dem_decYrs.append(dem_decYr)
            date_files.append([dem_file])
            dem_file_sensor = dem_file_sensors[dd]
            dem_sensors.append(dem_file_sensor)
        else:
            date_files[dem_dates.index(dem_date)].append(dem_file)

    # loop through the dates, adding in points from each file on that date
    dem_layers = []
    for dd in range(len(dem_dates)):
        print('          Working on date ' + dem_dates[dd]+ ' (' + str(dd + 1) + ' of ' + str(len(dem_dates)) + ')')

        # create the glacier field
        dem_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))
        count_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))
        # print('                Glacier field size: ' + str(np.shape(dem_grid)))

        dem_sensor = dem_sensors[dd]
        dem_date = dem_dates[dd]
        filename_strings.append(', '.join(date_files[dd]))

        for ff in range(len(date_files[dd])):
            dem_file = date_files[dd][ff]

            print('              Adding file '+dem_file+' ('+str(ff+1)+' of '+str(len(date_files[dd]))+')')

            # #read in the points
            # print('            Reading in the file points')


            valid_points = read_valid_points_from_csv(os.path.join(GD_object.data_folder, 'Elevation', 'IceBridge', 'Data', dem_file_sensor, dem_date,dem_file),dem_date)

            # print(valid_points)

            # print('              Array size: ' + str(np.shape(valid_points)))
            # print('              Domain of array:  x = ' + str(np.min(valid_points[:, 0])) + ' to ' + str(np.max(valid_points[:, 0])) + ',  y = ' + str(np.min(valid_points[:, 1])) + ' to ' + str(np.max(valid_points[:, 1])))

            #reproject the points to polar stereo
            # print('            Reprojecting the points to polar stereo')
            valid_points = reproject_polygon(valid_points,4326,3413,x_column=1,y_column=0)

            # plt.plot(valid_points[:, 0], valid_points[:, 1], 'k.')
            # plt.show()

            # put the points on the region grid
            # print('            Sampling the points onto the regional domain')
            dem_grid,count_grid = cumulate_icebridge_atm_dem_field(valid_points, GD_object.elevation_grid_x,
                                                       GD_object.elevation_grid_y,dem_grid,count_grid)



        dem_grid[count_grid > 0] = dem_grid[count_grid > 0] / count_grid[count_grid > 0]
        dem_grid[count_grid == 0] = -99
        dem_layers.append(dem_grid)

    return(dem_layers, dem_dates, dem_decYrs, filename_strings)


#######################################################################################
#These are the scripts for stacking the glacier fields into one

def save_icebridge_atm_layers(GD_object,dem_layers, dem_dates, dem_decYrs, dem_file_names):

    data_vars = {}
    for dd in range(len(dem_dates)):
        data_vars[dem_dates[dd]]=(['y','x'],dem_layers[dd])

    swath = xr.Dataset(data_vars,coords={'y': GD_object.elevation_grid_x,'x': GD_object.elevation_grid_y})

    for dd in range(len(dem_dates)):
        swath[dem_dates[dd]].attrs['file_name'] = dem_file_names[dd]
        swath[dem_dates[dd]].attrs['dec_yr'] = dem_decYrs[dd]

    output_file = os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Data',GD_object.region_name+' IceBridge ATM Elevation Grids.nc')
    print('        Outputting data to ' + output_file)

    swath.to_netcdf(output_file)

def save_icebridge_atm_layers_as_points(GD_object, dem_layers, dem_dates, dem_decYrs, dem_file_name_strings):

    ############################################################################
    # this first part converts the points to x,y,z,xi,yi columns grids
    message = '        Converting grids to points for efficient storage'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    first_grid = dem_layers[0]
    xi = np.arange(np.shape(first_grid)[1])
    yi = np.arange(np.shape(first_grid)[0])
    XI, YI = np.meshgrid(xi, yi)
    X, Y = np.meshgrid(GD_object.elevation_grid_x, GD_object.elevation_grid_y)

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

    if len(GD_object.icebridge_atm_output_file)>2:
        output_file = GD_object.icebridge_atm_output_file
    else:
        output_file = os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation', 'Data',
                                   GD_object.region_name + ' IceBridge ATM Elevation Points.nc')

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

        grp.file_name = dem_file_name_strings[dd]
        grp.dec_yr = dem_decYrs[dd]

    ds.close()

    message = '        Outputting data to ' + output_file
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)




#######################################################################################
#This is the main routine to process the icebridge_atm data

def generate_Icebridge_ATM_dataset(GD_object):

    message = '    Running compilation for the Operation Icebridge ATM data'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)

    if GD_object.overwrite_existing_elevation_stacks:
        continue_with_generation = True
    else:
        if GD_object.region_name+' IceBridge ATM Elevation Points.nc' in os.listdir(os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Data')):
            continue_with_generation = False
            message = '        IceBridge ATM elevation compilation skipped for ' + GD_object.region_name + ' - already done'
            GD_object.output_summary += '\n' + message
            if GD_object.print_main_outputs:
                print(message)
        else:
            continue_with_generation = True


    if continue_with_generation:

        # step 1: get a list of files for the region
        message = '        Finding a list of IceBridge ATM files which overlap the domain'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)
        dem_file_names, dem_file_dates, dem_file_sensors, dem_file_links = find_icebridge_atm_dem_files_in_domain(GD_object)
        print('          Found '+str(len(dem_file_names))+' files')

        if isinstance(GD_object.max_number_of_icebridge_atm_files,int):
            if len(dem_file_names)>GD_object.max_number_of_icebridge_atm_files:
                dem_file_names = dem_file_names[:GD_object.max_number_of_icebridge_atm_files]
                dem_file_dates = dem_file_dates[:GD_object.max_number_of_icebridge_atm_files]
                dem_file_sensors = dem_file_sensors[:GD_object.max_number_of_icebridge_atm_files]
                dem_file_links = dem_file_links[:GD_object.max_number_of_icebridge_atm_files]

        # step 2: download the data
        message = '        Downloading files (if not already available)'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)
        download_icebridge_atm_files(GD_object,dem_file_names, dem_file_dates, dem_file_sensors, dem_file_links)

        if GD_object.create_elevation_stacks:
            message = '        Resampling the IceBridge ATM data onto to the regional domain'
            GD_object.output_summary += '\n' + message
            if GD_object.print_sub_outputs:
                print(message)

        # step 3: stack the icebridge_atm data into layers
        dem_layers, dem_dates, dem_decYrs, dem_filename_strings = get_icebridge_atm_layers(GD_object,dem_file_names, dem_file_dates, dem_file_sensors)

        # step 4: save the icebridge_atm layer to an nc file
        save_icebridge_atm_layers_as_points(GD_object,dem_layers, dem_dates, dem_decYrs, dem_filename_strings)