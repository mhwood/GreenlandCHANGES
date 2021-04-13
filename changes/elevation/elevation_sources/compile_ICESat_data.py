
import os
import h5py
from datetime import datetime
import xarray as xr
import netCDF4 as nc4
import numpy as np
import itertools
import requests
from ....toolbox.time import YMD_to_DecYr
from ....toolbox.reprojection import reproject_polygon
from ....toolbox.series import series_to_N_points

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


########################################################################################################################
#These are the scripts for finding a list of ICESat files and their associated links for downloading
#Note, these codes were adapted from the downloading scripts provided by NSIDC
#See https://nsidc.org/data/atl06 for more details

def find_icesat_dem_files_from_nsidc_portal(GD_object):

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

        bbox = reproject_polygon(bbox, 3413, 4326)

        polygon = ''
        for bb in range(np.shape(bbox)[0]):
            polygon += str(bbox[bb, 0]) + ',' + str(bbox[bb, 1]) + ','

        polygon = polygon[:-1]

        return (polygon)

    # this function is from the ICESat data page
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

    # this function is from the ICESat data page
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

    # this function is from the ICESat data page
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

    short_name = 'GLAH12'
    version = '034'
    time_start = '2003-02-20T00:00:00Z'
    time_end = '2009-10-11T23:59:59Z'

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
    outputFiles = []
    outputDates = []
    outputSets = []
    for url in urls:
        if url[-3:] == '.H5':
            urlParts = url.split('/')
            fileID = urlParts[-1][:-3]
            outputFiles.append(fileID+'.h5')
            date = url.split('/')[6].replace('.','')
            outputDates.append(date)
            # date = fileID.split('_')[1][:8]
            outputUrls.append(url)
            outputSets.append([date,fileID,url])

    textOutput = ''
    for outputSet in outputSets:
        textOutput += ','.join(outputSet) + '\n'

    f = open(os.path.join(GD_object.project_folder, GD_object.region_name,'Elevation','Metadata',GD_object.region_name + ' ICESat Files.csv'),'w')
    f.write(textOutput[:-1])
    f.close()

    return (outputDates,outputFiles,outputUrls)

def find_icesat_dem_files_in_domain(GD_object):
    if GD_object.region_name + ' ICESat Files.csv' not in os.listdir(os.path.join(GD_object.project_folder, GD_object.region_name,'Elevation','Metadata')):
        dem_dates, dem_files, dem_file_links = find_icesat_dem_files_from_nsidc_portal(GD_object)
    else:
        f=open(os.path.join(GD_object.project_folder, GD_object.region_name,'Elevation','Metadata',GD_object.region_name + ' ICESat Files.csv'))
        lines=f.read()
        f.close()
        lines=lines.split('\n')
        dem_dates = []
        dem_files =[]
        dem_file_links=[]
        for line in lines:
            line=line.split(',')
            if len(line)>1:
                dem_dates.append(line[0])
                dem_files.append(line[1]+'.h5')
                dem_file_links.append(line[2])

    output_dem_dates = []
    output_dem_files = []
    output_dem_file_links = []
    for df in range(len(dem_files)):
        dem_file=dem_files[df]
        dem_date_str = dem_dates[df]
        date_test = datetime(int(dem_date_str[:4]),int(dem_date_str[4:6]),int(dem_date_str[6:8]))
        if date_test >= GD_object.date_1 and date_test < GD_object.date_2:
            output_dem_dates.append(dem_date_str)
            output_dem_files.append(dem_file)
            output_dem_file_links.append(dem_file_links[df])


    return(output_dem_dates,output_dem_files,output_dem_file_links)


#######################################################################################
#These are the scripts for downloading the ICESat data

def downloadStrip(dataFolder,yearMonthDay,file_name,url):

    outputFile = os.path.join(dataFolder, yearMonthDay, file_name)

    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(max_retries=20)
    session.mount('https://', adapter)
    session.mount('http://', adapter)
    resp = session.get(url)

    # resp = requests.get(url)
    open(outputFile, 'wb').write(resp.content)

def download_icesat_files(GD_object,dem_dates,dem_files,dem_file_links):

    if 'ICESat' not in os.listdir(os.path.join(GD_object.data_folder,'Elevation')):
        os.mkdir(os.path.join(GD_object.data_folder,'Elevation','ICESat'))
    if 'Data' not in os.listdir(os.path.join(GD_object.data_folder,'Elevation','ICESat')):
        os.mkdir(os.path.join(GD_object.data_folder,'Elevation','ICESat','Data'))

    dataFolder = os.path.join(GD_object.data_folder,'Elevation','ICESat','Data')

    if GD_object.download_new_icesat_data:
        for dd in range(len(dem_files)):
            dem_file = dem_files[dd]
            if GD_object.print_sub_outputs:
                print('            Checking file '+dem_file)
            yearMonthDay = dem_dates[dd]
            if yearMonthDay not in os.listdir(dataFolder):
                os.mkdir(os.path.join(dataFolder, yearMonthDay))

            if dem_file not in os.listdir(os.path.join(dataFolder, yearMonthDay)):
                download_file=True
                if GD_object.print_sub_outputs:
                    print('              Downloading file ' + dem_file)
                    downloadStrip(dataFolder, yearMonthDay, dem_file, dem_file_links[dd])
            else:
                if GD_object.print_sub_outputs:
                    print('            ' + dem_file + ' file already downloaded')
                download_file = False



#######################################################################################
#These are the scripts for creating the layers in the ICESat data


def readXYZsFromh5file(GD_object,filepath):

    # see the page here for documentation:
    # https://nsidc.org/data/glas/data-dictionary-glah12

    # specifically, this part:
    # Group: Data_40HZ/Elevation_Surfaces/d_elev
    # Surface elevation with respect to the ellipsoid at the spot location determined by the ice-sheet
    # specific range after instrument corrections, atmospheric delays and tides have been applied.
    # The saturation elevation correction (d_satElevCorr) has not been applied and needs to be added
    # to this elevation. This can be over a one meter correction. If it is invalid then the elevation
    # should not be used. The saturation correction flag (sat_corr_flg) is an important flag to understand
    # the possible quality of the elevation data. The saturation index (i_satNdx) can be used for more
    # understanding of concerns on data quality from saturation effects. Also no correction for pulse
    # spreading from forward scatter has been applied.


    f = h5py.File(filepath, 'r')
    valid_points = np.array([[0,0],[0,0]])

    # #gtxx groups have the elevation data
    # #in each gtxx group, there is a land_ice_segments group
    # #this group has fields for longitude, latitude, h_li, and a boolean flag
    xyz_started=False

    elev_grp = f['Data_40HZ']['Elevation_Surfaces']
    elev = elev_grp.get('d_elev')[:]

    corr_grp = f['Data_40HZ']['Elevation_Corrections']
    corr = corr_grp.get('d_satElevCorr')[:]

    qual_grp = f['Data_40HZ']['Quality']
    flag = qual_grp.get('elev_use_flg')[:]

    loc_grp = f['Data_40HZ']['Geolocation']
    longitude = loc_grp.get('d_lon')[:]
    latitude = loc_grp.get('d_lat')[:]

    xyz = np.hstack([np.reshape(longitude,(len(longitude),1)),
                      np.reshape(latitude, (len(latitude), 1)),
                      np.reshape(elev, (len(elev), 1))])

    valid_points=xyz[np.logical_and(flag==0,np.abs(corr)<1e2),:]
    if GD_object.print_sub_outputs:
       print('              File has '+str(np.shape(valid_points)[0])+' valid points out of '+str(np.shape(xyz)[0])+' total')

    return(valid_points)

def create_icesat_dem_field(GD_object,valid_points):

    #limit the points to the region domain
    dem_data = valid_points
    dem_data = dem_data[np.logical_and(dem_data[:, 0] <= np.max(GD_object.elevation_grid_x), dem_data[:, 0] >= np.min(GD_object.elevation_grid_x)), :]
    dem_data = dem_data[np.logical_and(dem_data[:, 1] <= np.max(GD_object.elevation_grid_y), dem_data[:, 1] >= np.min(GD_object.elevation_grid_y))]

    # create the glacier field
    sum_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))
    count_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))
    # print('                Glacier field size: ' + str(np.shape(dem_grid)))

    #fill in the glacier field by averaging points based on their closest point
    for dd in range(np.shape(dem_data)[0]):
        x_index = np.argmin(np.abs(GD_object.elevation_grid_x-dem_data[dd,0]))
        y_index = np.argmin(np.abs(GD_object.elevation_grid_y-dem_data[dd,1]))
        sum_grid[y_index,x_index]+=dem_data[dd,2]
        count_grid[y_index,x_index]+=1

    return(sum_grid,count_grid)

def get_icesat_layers(GD_object,input_dem_dates,dem_file_names):

    #get a unique set of dates and the files that correspond to them
    unique_dem_dates = []
    dem_file_lists = []
    for df in range(len(dem_file_names)):
        dem_file = dem_file_names[df]
        dem_date = input_dem_dates[df]
        if dem_date not in unique_dem_dates:
            unique_dem_dates.append(dem_date)
            dem_file_lists.append([dem_file])
        else:
            dem_file_lists[unique_dem_dates.index(dem_date)].append(dem_file)

    dem_dates = []
    dem_layers = []
    dem_decYrs = []
    dem_filename_strings = []
    for dd in range(len(unique_dem_dates)):

        dem_date = unique_dem_dates[dd]
        dem_files = dem_file_lists[dd]
        dem_files_used = []

        if GD_object.print_sub_outputs:
            print('          Working on date ' + dem_date + ' (' + str(dd + 1) + ' of ' + str(len(unique_dem_dates)) + ')')

        sum_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))
        count_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))
        dem_decYr = YMD_to_DecYr(int(dem_date[:4]), int(dem_date[4:6]), int(dem_date[6:8]))

        # add data from each file to the dem grid
        for df in range(len(dem_files)):
            dem_file = dem_files[df]
            if GD_object.print_sub_outputs:
                print('            Working on file '+dem_file+' ('+str(df+1)+' of '+str(len(dem_files))+')')

            valid_points = readXYZsFromh5file(GD_object,os.path.join(GD_object.data_folder,'Elevation','ICESat','Data',dem_date,dem_file))

            #reproject the points to polar stereo
            if GD_object.print_sub_outputs:
                print('              Reprojecting the points to polar stereo')
            if len(valid_points)>3:
                valid_points = reproject_polygon(valid_points,4326,3413)

            #put the points on the region grid
            if GD_object.print_sub_outputs:
                print('              Sampling the points onto the regional domain')
            layer_sum_grid, layer_count_grid = create_icesat_dem_field(GD_object,valid_points)

            if np.sum(layer_count_grid)>0:
                sum_grid[layer_count_grid>0]+=layer_sum_grid[layer_count_grid>0]
                count_grid[layer_count_grid > 0] += layer_count_grid[layer_count_grid > 0]
                dem_files_used.append(dem_file)

        dem_grid = np.copy(sum_grid)
        dem_grid[count_grid > 0] = dem_grid[count_grid > 0] / count_grid[count_grid > 0]
        dem_grid[count_grid == 0] = -99

        if np.any(dem_grid > -99):
            if GD_object.print_sub_outputs:
                print('              Valid points found - layer added to main output')
            dem_decYrs.append(dem_decYr)
            dem_dates.append(dem_date)
            dem_layers.append(dem_grid)
            dem_filename_strings.append(','.join(dem_files_used))
        else:
            if GD_object.print_sub_outputs:
                print('              No valid points found within array')

    return(dem_layers, dem_dates, dem_decYrs,dem_filename_strings)


#######################################################################################
#These are the scripts for stacking the glacier fields into one

def save_icesat_layers_as_grids(GD_object, dem_layers, dem_dates, dem_decYrs, dem_file_names):

    data_vars = {}
    for dd in range(len(dem_dates)):
        data_vars[dem_dates[dd]]=(['y','x'],dem_layers[dd])

    swath = xr.Dataset(data_vars,coords={'y': GD_object.elevation_grid_y,'x': GD_object.elevation_grid_x})

    for dd in range(len(dem_dates)):
        swath[dem_dates[dd]].attrs['file_name'] = dem_file_names[dd]
        swath[dem_dates[dd]].attrs['dec_yr'] = dem_decYrs[dd]

    output_file = os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Data',GD_object.region_name+' ICESat Elevation Grids.nc')
    message = '        Outputting data to ' + output_file
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    swath.to_netcdf(output_file)

def save_icesat_layers_as_points(GD_object, dem_layers, dem_dates, dem_decYrs, dem_file_names):

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

    if len(GD_object.icesat_output_file)>2:
        output_file = GD_object.icesat_output_file
    else:
        output_file = os.path.join(GD_object.project_folder, GD_object.region_name, 'Elevation', 'Data',
                                   GD_object.region_name + ' ICESat Elevation Points.nc')

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



#######################################################################################
#This is the script to run the full ICESat compilation


def generate_ICESat_dataset(GD_object):

    message = '    Running compilation for the ICESat data'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)

    if not GD_object.save_icesat_points_as_grids and not GD_object.save_icesat_points_as_points:
        raise ValueError('    Need to specify either save_icesat_points_as_grids or save_icesat_points_as_points')

    if GD_object.overwrite_existing_elevation_stacks:
        continue_with_generation = True
    else:
        if GD_object.save_icesat_points_as_grids:
            if GD_object.region_name+' ICESat Elevation Grids.nc' in os.listdir(os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Data')):
                continue_with_generation = False
                message = '        ICESat elevation compilation skipped for ' + GD_object.region_name + ' - already done'
                GD_object.output_summary += '\n' + message
                if GD_object.print_main_outputs:
                    print(message)
            else:
                continue_with_generation = True
        else:
            if GD_object.region_name+' ICESat Elevation Points.nc' in os.listdir(os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Data')):
                continue_with_generation = False
                message = '        ICESat elevation compilation skipped for ' + GD_object.region_name + ' - already done'
                GD_object.output_summary += '\n' + message
                if GD_object.print_main_outputs:
                    print(message)
            else:
                continue_with_generation = True

    if continue_with_generation:

        # step 1: get a list of files for the region
        message = '        Finding a list of ICESat files which overlap the domain'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

        dem_dates, dem_file_names, dem_file_links = find_icesat_dem_files_in_domain(GD_object)

        message = '        Found '+str(len(dem_file_names))+' files'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

        if isinstance(GD_object.max_number_of_icesat_files,int):
            if len(dem_file_names)>GD_object.max_number_of_icesat_files:
                dem_dates = dem_dates[:GD_object.max_number_of_icesat_files]
                dem_file_names = dem_file_names[:GD_object.max_number_of_icesat_files]
                dem_file_links = dem_file_links[:GD_object.max_number_of_icesat_files]

        # step 2: download the data

        message = '        Downloading files (if not already available)'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)
        download_icesat_files(GD_object,dem_dates,dem_file_names,dem_file_links)

        if GD_object.create_elevation_stacks:
            message = '        Resampling the ICESat data onto to the regional domain'
            GD_object.output_summary += '\n' + message
            if GD_object.print_sub_outputs:
                print(message)

            # step 3: stack the icesat data into layers
            dem_layers, dem_dates, dem_decYrs, dem_filename_strings = get_icesat_layers(GD_object,dem_dates,dem_file_names)

            if len(dem_layers)>0:

                # step 4: save the icesat layer to an nc file
                if GD_object.save_icesat_points_as_grids:
                    save_icesat_layers_as_grids(GD_object, dem_layers, dem_dates, dem_decYrs, dem_filename_strings)

                if GD_object.save_icesat_points_as_points:
                    save_icesat_layers_as_points(GD_object, dem_layers, dem_dates, dem_decYrs, dem_filename_strings)