
import os
import h5py
import xarray as xr
import numpy as np
import itertools
import requests
from toolbox.time import YMD_to_DecYr
from toolbox.reprojection import reproject_polygon
from toolbox.series import series_to_N_points


########################################################################################################################
#These are the scripts for finding a list of ICESat2 files and their associated links for downloading
#Note, these codes were adapted from the downloading scripts provided by NSIDC
#See https://nsidc.org/data/atl06 for more details

def find_icesat2_dem_files_from_nsidc_portal(GD_object):

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

    # this function is from the ICESat2 data page
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

    # this function is from the ICESat2 data page
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

    # this function is from the ICESat2 data page
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

    short_name = 'ATL06'
    version = '003'

    time_start = '2018-01-01T00:00:00Z'
    time_end = '2021-01-01T00:00:00Z'

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
    outputSets = []
    for url in urls:
        if url[-3:] == '.h5':
            urlParts = url.split('/')
            fileID = urlParts[-1][:-3]
            outputFiles.append(fileID+'.h5')
            date = fileID.split('_')[1][:8]
            outputUrls.append(url)
            outputSets.append([date,fileID,url])

    textOutput = ''
    for outputSet in outputSets:
        textOutput += ','.join(outputSet) + '\n'

    f = open(os.path.join(GD_object.project_folder, GD_object.region_name,'Elevation','Metadata',GD_object.region_name + ' ICESat2 Files.csv'),'w')
    f.write(textOutput[:-1])
    f.close()

    return (outputFiles,outputUrls)

def find_icesat2_dem_files_in_domain(GD_object):
    if GD_object.region_name + ' ICESat2 Files.csv' not in os.listdir(os.path.join(GD_object.project_folder, GD_object.region_name,'Elevation','Metadata')):
        dem_files, dem_file_links = find_icesat2_dem_files_from_nsidc_portal(GD_object)
    else:
        f=open(os.path.join(GD_object.project_folder, GD_object.region_name,'Elevation','Metadata',GD_object.region_name + ' ICESat2 Files.csv'))
        lines=f.read()
        f.close()
        lines=lines.split('\n')
        dem_files=[]
        dem_file_links=[]
        for line in lines:
            line=line.split(',')
            if len(line)>1:
                dem_files.append(line[1]+'.h5')
                dem_file_links.append(line[2])
    return(dem_files,dem_file_links)


#######################################################################################
#These are the scripts for downloading the ICESat2 data

def downloadStrip(GD_object,stripFolder,yearMonth,stripID,url):

    if yearMonth not in os.listdir(stripFolder):
        os.mkdir(os.path.join(stripFolder,yearMonth))

    if stripID+'.h5' not in os.listdir(os.path.join(stripFolder,yearMonth)):
        outputFile=os.path.join(stripFolder,yearMonth,stripID+'.h5')
        resp = requests.get(url)
        open(outputFile, 'wb').write(resp.content)
    else:
        if GD_object.print_sub_outputs:
            print('            '+stripID+' file already downloaded and processed')

def download_icesat2_files(GD_object,dem_files,dem_file_links):

    dataFolder = os.path.join(GD_object.data_folder,'Elevation','ICESat2','Data')

    if GD_object.download_new_icesat2_data:
        for dd in range(len(dem_files)):
            dem_file = dem_files[dd]
            if GD_object.print_sub_outputs:
                print('            Checking file '+dem_file)
            yearMonth = dem_file.split('_')[1][:8]
            download_file=True

            if yearMonth in os.listdir(dataFolder):
                for fil in os.listdir(dataFolder+'/'+yearMonth):
                    if fil[:20]==dem_file[:20]: #this allows for older versions to be kept
                        download_file=False

            if download_file:
                if GD_object.print_sub_outputs:
                    print('              Downloading file ' + dem_file)
                downloadStrip(GD_object,dataFolder, yearMonth, dem_file[:-3], dem_file_links[dd])
            else:
                if GD_object.print_sub_outputs:
                    print('              File already downloaded')


#######################################################################################
#These are the scripts for creating the layers in the ICESat2 data


def readXYZsFromh5file(GD_object,filepath):
    #see page 53 of this file for a description:
    #https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL06_ATBD_r001.pdf

    f = h5py.File(filepath, 'r')
    valid_points = np.array([[0,0],[0,0]])


    # #gtxx groups have the elevation data
    # #in each gtxx group, there is a land_ice_segments group
    # #this group has fields for longitude, latitude, h_li, and a boolean flag
    xyz_started=False

    for beam in ['gt1r','gt1l','gt2r','gt2l','gt3r','gt3l']:
         grp = f[beam]['land_ice_segments']
         longitude = np.array(list(grp['longitude']))
         latitude = np.array(list(grp['latitude']))
         height = np.array(list(grp['h_li']))
         flag = np.array(list(grp['atl06_quality_summary']))
         xyz = np.hstack([np.reshape(longitude,(len(longitude),1)),
                          np.reshape(latitude, (len(latitude), 1)),
                          np.reshape(height, (len(height), 1))])

         xyz=xyz[flag==0,:]
         if GD_object.print_sub_outputs:
            print('              Beam '+beam+' has '+str(np.shape(xyz)[0])+' valid points')

         if not xyz_started:
             valid_points = np.copy(xyz)
             xyz_started=True
         else:
             valid_points=np.vstack([valid_points,xyz])

    return(valid_points)

def create_icesat2_dem_field(GD_object,valid_points):



    #limit the points to the region domain
    dem_data = valid_points
    dem_data = dem_data[np.logical_and(dem_data[:, 0] <= np.max(GD_object.elevation_grid_x), dem_data[:, 0] >= np.min(GD_object.elevation_grid_x)), :]
    dem_data = dem_data[np.logical_and(dem_data[:, 1] <= np.max(GD_object.elevation_grid_y), dem_data[:, 1] >= np.min(GD_object.elevation_grid_y))]

    # create the glacier field
    dem_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))
    count_grid = np.zeros((len(GD_object.elevation_grid_y), len(GD_object.elevation_grid_x)))
    # print('                Glacier field size: ' + str(np.shape(dem_grid)))

    #fill in the glacier field by averaging points based on their closest point
    for dd in range(np.shape(dem_data)[0]):
        x_index = np.argmin(np.abs(GD_object.elevation_grid_x-dem_data[dd,0]))
        y_index = np.argmin(np.abs(GD_object.elevation_grid_y-dem_data[dd,1]))
        dem_grid[y_index,x_index]+=dem_data[dd,2]
        count_grid[y_index,x_index]+=1

    dem_grid[count_grid>0] = dem_grid[count_grid>0]/count_grid[count_grid>0]
    dem_grid[count_grid == 0] = -99

    return(dem_grid)

def get_icesat2_layers(GD_object,dem_file_names):

    dem_layers=[]
    dem_dates=[]
    dem_decYrs = []
    for dd in range(len(dem_file_names)):
        dem_file = dem_file_names[dd]
        if GD_object.print_sub_outputs:
            print('          Working on file '+dem_file+' ('+str(dd+1)+' of '+str(len(dem_file_names))+')')
        dem_date = dem_file.split('_')[1][:8]
        dem_dates.append(dem_date)
        dem_decYr = YMD_to_DecYr(int(dem_date[:4]), int(dem_date[4:6]), int(dem_date[6:8]))
        dem_decYrs.append(dem_decYr)


        #read in the points
        if GD_object.print_sub_outputs:
            print('            Reading in the file points')
        valid_points = readXYZsFromh5file(GD_object,os.path.join(GD_object.data_folder,'Elevation','ICESat2','Data',dem_date,dem_file))
        # print('              Array size: ' + str(np.shape(valid_points)))
        # print('              Domain of array:  x = ' + str(np.min(valid_points[:, 0])) + ' to ' + str(np.max(valid_points[:, 0])) + ',  y = ' + str(np.min(valid_points[:, 1])) + ' to ' + str(np.max(valid_points[:, 1])))

        #reproject the points to polar stereo
        if GD_object.print_sub_outputs:
            print('            Reprojecting the points to polar stereo')
        valid_points = reproject_polygon(valid_points,4326,3413)

        #put the points on the region grid
        if GD_object.print_sub_outputs:
            print('            Sampling the points onto the regional domain')
        dem_layer = create_icesat2_dem_field(GD_object,valid_points)
        dem_layers.append(dem_layer)

    return(dem_layers, dem_dates, dem_decYrs)


#######################################################################################
#These are the scripts for stacking the glacier fields into one

def save_icesat2_layers(GD_object, dem_layers, dem_dates, dem_decYrs, dem_file_names):

    data_vars = {}
    for dd in range(len(dem_dates)):
        data_vars[dem_dates[dd]]=(['y','x'],dem_layers[dd])

    swath = xr.Dataset(data_vars,coords={'y': GD_object.elevation_grid_y,'x': GD_object.elevation_grid_x})

    for dd in range(len(dem_dates)):
        swath[dem_dates[dd]].attrs['file_name'] = dem_file_names[dd]
        swath[dem_dates[dd]].attrs['dec_yr'] = dem_decYrs[dd]

    output_file = os.path.join(GD_object.project_folder,GD_object.region_name,'Elevation','Data',GD_object.region_name+' ICESat2 Elevation Grids.nc')
    message = '        Outputting data to ' + output_file
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    swath.to_netcdf(output_file)




#######################################################################################
#This is the script to run the full ICESat-2 compilation


def generate_ICESat2_dataset(GD_object):

    message = '    Running compilation for the ICESat-2 data'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)

    # step 1: get a list of files for the region
    message = '        Finding a list of ICESat-2 files which overlap the domain'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    dem_file_names, dem_file_links = find_icesat2_dem_files_in_domain(GD_object)

    message = '        Found '+str(len(dem_file_names))+' files'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    if isinstance(GD_object.max_number_of_icesat2_files,int):
        if len(dem_file_names)>GD_object.max_number_of_icesat2_files:
            dem_file_names = dem_file_names[:GD_object.max_number_of_icesat2_files]
            dem_file_links = dem_file_links[:GD_object.max_number_of_icesat2_files]

    # step 2: download the data
    message = '        Downloading files (if not already available)'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)
    download_icesat2_files(GD_object,dem_file_names,dem_file_links)

    print('        Resampling the ICESat2 data onto to the regional domain')
    # step 3: stack the icesat2 data into layers
    dem_layers, dem_dates, dem_decYrs = get_icesat2_layers(GD_object,dem_file_names)

    # step 4: save the icesat2 layer to an nc file
    save_icesat2_layers(GD_object, dem_layers, dem_dates, dem_decYrs, dem_file_names)