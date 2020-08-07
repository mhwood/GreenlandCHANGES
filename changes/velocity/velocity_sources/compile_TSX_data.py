from toolbox.series import series_to_N_points
from toolbox.reprojection import reproject_polygon
import numpy as np
import itertools
import requests
import os
import netCDF4 as nc4
import gdal

def obtain_list_of_tsx_files_overlapping_domain(GD_object):

    if GD_object.glacier_name+' TSX Files.csv' in os.listdir(os.path.join(GD_object.project_folder,
                                                                                    GD_object.glacier_name, 'Velocity', 'Metadata')):
        f=open(os.path.join(GD_object.project_folder, GD_object.glacier_name, 'Velocity', 'Metadata',
                                      GD_object.glacier_name+' TSX Files.csv'))
        lines=f.read()
        f.close()
        lines=lines.split('\n')
        lines.pop(0)
        tsx_file_links = []
        tsx_file_names = []
        for line in lines:
            line=line.split(',')
            tsx_file_names.append(line[0])
            tsx_file_links.append(line[1])
    else:

        def create_polygon_string(GD_object):

            domain = np.array([[GD_object.extents[0],GD_object.extents[1]],
                               [GD_object.extents[2],GD_object.extents[1]],
                               [GD_object.extents[2],GD_object.extents[3]],
                               [GD_object.extents[0],GD_object.extents[3]],
                               [GD_object.extents[0],GD_object.extents[1]]])
            domain = series_to_N_points(domain,20)
            domain = reproject_polygon(domain,3413,4326)

            polygon = ''
            for pp in range(np.shape(domain)[0]):
                polygon += str(domain[pp, 0]) + ',' + str(domain[pp, 1]) + ','
            polygon = polygon[:-1]

            return (polygon)

        # this function is from NSIDC
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

        # this function is from NSIDC
        def build_cmr_query_url(CMR_FILE_URL, short_name, version, time_start, time_end, polygon=None,
                                filename_filter=None):
            params = '&short_name={0}'.format(short_name)
            params += build_version_query_params(version)
            params += '&temporal[]={0},{1}'.format(time_start, time_end)
            if polygon:
                params += '&polygon={0}'.format(polygon)
            if filename_filter:
                params += '&producer_granule_id[]={0}&options[producer_granule_id][pattern]=true'.format(
                    filename_filter)
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

        #this script below is adapted from the default downloading script provided by NSIDC
        short_name = 'NSIDC-0481'
        version = '1'
        time_start = '2015-01-01T00:00:00Z'
        time_end = '2025-01-01T23:59:59Z'
        polygon = create_polygon_string(GD_object)
        filename_filter = '*'

        CMR_URL = 'https://cmr.earthdata.nasa.gov'
        CMR_PAGE_SIZE = 2000
        CMR_FILE_URL = ('{0}/search/granules.json?provider=NSIDC_ECS'
                        '&sort_key[]=start_date&sort_key[]=producer_granule_id'
                        '&scroll=true&page_size={1}'.format(CMR_URL, CMR_PAGE_SIZE))

        cmr_query_url = build_cmr_query_url(CMR_FILE_URL, short_name=short_name, version=version,
                                            time_start=time_start, time_end=time_end,
                                            polygon=polygon, filename_filter=filename_filter)

        response = requests.get(cmr_query_url)
        search_page = response.json()
        urls = cmr_filter_urls(search_page)

        tsx_file_links = []
        tsx_file_names = []
        for url in urls:
            if 'vx_v1.2.tif' in url or 'vy_v1.2.tif' in url:
                tsx_file_links.append(url)
                url_parts = url.split('/')
                file_name = url_parts[-1]
                tsx_file_names.append(file_name)

        text_output = 'File_Name,URL'
        for ff in range(len(tsx_file_names)):
            text_output+='\n'+tsx_file_names[ff]+','+tsx_file_links[ff]

        f = open(os.path.join(GD_object.project_folder, GD_object.glacier_name, 'Velocity', 'Metadata',
                                        GD_object.glacier_name + ' TSX Files.csv'),'w')
        f.write(text_output)
        f.close()

    return(tsx_file_links, tsx_file_names)

def obtain_download_list(GD_object,tsx_file_names):
    golive_folder=os.path.join(GD_object.data_folder,'Velocity','TSX')
    download_list=[]
    for fil in tsx_file_names:
        location_id = fil.split('_')[1]
        if location_id not in os.listdir(golive_folder):
            download_list.append(fil)
        else:
            if fil not in os.listdir(os.path.join(golive_folder,location_id)):
                download_list.append(fil)
    n_existing_files = len(tsx_file_names)-len(download_list)
    return(download_list,n_existing_files)

def download_tsx_files(GD_object,tsx_file_links, tsx_file_names,download_list):

    for ff in range(len(tsx_file_names)):
        url = tsx_file_links[ff]
        file_name = tsx_file_names[ff]
        if file_name in download_list:
            if GD_object.print_sub_outputs:
                print('                Downloading '+file_name)
            location_id = file_name.split('_')[1]
            output_file = os.path.join(GD_object.data_folder,'Velocity','TSX',location_id,file_name)

            resp = requests.get(url)
            open(output_file, 'wb').write(resp.content)

def fileID_to_date_pair(fileID):
    def monthStringToInt(monthString):
        if monthString == 'Jan':
            month = 1
        if monthString == 'Feb':
            month = 2
        if monthString == 'Mar':
            month = 3
        if monthString == 'Apr':
            month = 4
        if monthString == 'May':
            month = 5
        if monthString == 'Jun':
            month = 6
        if monthString == 'Jul':
            month = 7
        if monthString == 'Aug':
            month = 8
        if monthString == 'Sep':
            month = 9
        if monthString == 'Oct':
            month = 10
        if monthString == 'Nov':
            month = 11
        if monthString == 'Dec':
            month = 12
        return (month)

    date1 = fileID.split('_')[2]
    date2 = fileID.split('_')[3]

    year1 = int('20' + date1[-2:])
    month1 = monthStringToInt(date1[2:-2])
    day1 = int(date1[:2])

    year2 = int('20' + date2[-2:])
    month2 = monthStringToInt(date2[2:-2])
    day2 = int(date2[:2])

    datePair = str(year1) + "{:02d}".format(int(month1)) + "{:02d}".format(int(day1)) + '-' + str(
        year2) + "{:02d}".format(int(month2)) + "{:02d}".format(int(day2))
    return(datePair)


def create_velocity_stack(GD_object,tsx_file_names):
    X, Y = np.meshgrid(GD_object.velocity_grid_x, GD_object.velocity_grid_y)

    #find list of file IDs which contain both x and y grids
    fileIDs=[]
    for file_name in tsx_file_names:
        fileID=file_name[:-12]
        if fileID+'_vy_v1.2.tif' in tsx_file_names and fileID+'_vx_v1.2.tif' in tsx_file_names and fileID not in fileIDs:
            fileIDs.append(fileID)

    vx_grids=[]
    vy_grids=[]
    v_grids=[]
    date_pairs=[]
    source_lists=[]

    for ff in range(len(fileIDs)):
        fileID=fileIDs[ff]
        location_id = fileID.split('_')[1]
        date_pair = fileID_to_date_pair(fileID)
        sources = fileID+'_vx_v1.2.tif,'+fileID + '_vy_v1.2.tif'
        message = '            Looking for velocity points in ' + fileID + ' (file '+str(ff+1)+' of '+str(len(fileIDs))+')'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

        vx_file = os.path.join(GD_object.data_folder, 'Velocity', 'TSX', location_id, fileID+'_vx_v1.2.tif')
        ds = gdal.Open(vx_file)
        vx_array = np.array(ds.GetRasterBand(1).ReadAsArray())
        transform = ds.GetGeoTransform()
        ds = None

        vy_file = os.path.join(GD_object.data_folder, 'Velocity', 'TSX', location_id, fileID + '_vy_v1.2.tif')
        ds = gdal.Open(vy_file)
        vy_array = np.array(ds.GetRasterBand(1).ReadAsArray())
        transform = ds.GetGeoTransform()
        ds = None

        x = np.arange(transform[0], transform[0] + np.shape(vx_array)[1] * transform[1], transform[1])
        y = np.arange(transform[3], transform[3] + np.shape(vx_array)[0] * transform[5], transform[5])

        xIndices = np.logical_and(x >= GD_object.extents[0], x <= GD_object.extents[2])
        yIndices = np.logical_and(y >= GD_object.extents[1], y <= GD_object.extents[3])

        x = x[xIndices]
        y = y[yIndices]

        vx_array = vx_array[yIndices, :]
        vx_array = vx_array[:, xIndices]
        vy_array = vy_array[yIndices, :]
        vy_array = vy_array[:, xIndices]

        vxSumGrid = np.zeros_like(X).astype(float)
        vySumGrid = np.zeros_like(X).astype(float)
        vSumGrid = np.zeros_like(X).astype(float)
        countGrid = np.zeros_like(X).astype(float)

        for xi in range(len(x)):
            xIndex = np.argmin(np.abs(GD_object.velocity_grid_x - x[xi]))
            for yi in range(len(y)):
                yIndex = np.argmin(np.abs(GD_object.velocity_grid_y - y[yi]))
                vx = vx_array[yi,xi]
                vy = vy_array[yi,xi]
                v = (vx**2 + vy**2)**0.5
                if v<1e6:
                    vxSumGrid[yIndex, xIndex] += vx
                    vySumGrid[yIndex, xIndex] += vy
                    vSumGrid[yIndex, xIndex] += v
                    countGrid[yIndex, xIndex] += 1

        VX = -99999.0 * np.ones_like(X)
        VY = -99999.0 * np.ones_like(X)
        V = -99999.0 * np.ones_like(X)
        VX[countGrid > 0] = vxSumGrid[countGrid > 0] / countGrid[countGrid > 0]
        VY[countGrid > 0] = vySumGrid[countGrid > 0] / countGrid[countGrid > 0]
        V[countGrid > 0] = vSumGrid[countGrid > 0] / countGrid[countGrid > 0]

        if np.any(V > -99999):
            vx_grids.append(VX)
            vy_grids.append(VY)
            v_grids.append(V)
            date_pairs.append(date_pair)
            source_lists.append(sources)
    return(vx_grids, vy_grids, v_grids, date_pairs, source_lists)

def output_data_stack(GD_object, vx_grids, vy_grids, v_grids, date_pairs, source_lists):
    output_file = os.path.join(GD_object.project_folder, GD_object.glacier_name, 'Velocity', 'Data',
                               GD_object.glacier_name + ' TSX Velocity Data.nc')

    data = nc4.Dataset(output_file, "w", format="NETCDF4")

    data.createDimension('y', len(GD_object.velocity_grid_y))
    data.createDimension('x', len(GD_object.velocity_grid_x))
    xvar = data.createVariable('x', 'f4', ("x",))
    yvar = data.createVariable('y', 'f4', ("y",))

    xvar[:] = GD_object.velocity_grid_x
    yvar[:] = GD_object.velocity_grid_y

    for dd in range(len(date_pairs)):
        grp = data.createGroup(date_pairs[dd])
        vx = grp.createVariable('VX', 'f4', ("y", "x"))
        vy = grp.createVariable('VY', 'f4', ("y", "x"))
        v = grp.createVariable('V', 'f4', ("y", "x"))

        vx[:, :] = vx_grids[dd]
        vy[:, :] = vy_grids[dd]
        v[:, :] = v_grids[dd]
        grp.source_files = source_lists[dd]

    data.close()

    message = '        Saved TSX file to ' + output_file
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)


def generate_TSX_dataset(GD_object,testing=False):

    #############################################################################################
    # This section is for downloading all of the files overlapping the domain

    message = '    Running compilation for the TSX data'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)

    tsx_file_links, tsx_file_names = obtain_list_of_tsx_files_overlapping_domain(GD_object)

    message = '        Found ' + str(len(tsx_file_names)) + ' TSX files overlapping the domain'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    if GD_object.download_new_tsx_data:

        if GD_object.overwrite_existing_velocity_data:
            download_list = tsx_file_names
            message = '            Overwriting old files by request - downloading ' + str(len(download_list)) + ' files'
            GD_object.output_summary += '\n' + message
            if GD_object.print_sub_outputs:
                print(message)
        else:
            download_list, n_existing_files = obtain_download_list(GD_object, tsx_file_names)
            message = '            Found ' + str(n_existing_files) + ' existing files - downloading ' + str(len(download_list)) + ' files'
            GD_object.output_summary += '\n' + message
            if GD_object.print_sub_outputs:
                print(message)

        download_tsx_files(GD_object,tsx_file_links, tsx_file_names,download_list)

    else:
        message = '        Skipping the download of new files by request'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

    if GD_object.overwrite_existing_golive_stack:
        stack_data = True
    else:
        if GD_object.glacier_name + ' TSX Velocity Data.nc' in os.listdir(os.path.join(GD_object.project_folder,
                                                                                         GD_object.glacier_name,
                                                                                         'Velocity', 'Data')):
            stack_data = False
        else:
            stack_data = True

    if stack_data:
        message = '        Stacking TSX data into a common grid'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)
        #############################################################################################
        # This section is for stacking the files into a common grid

        vx_grids, vy_grids, v_grids, date_pairs, source_lists = create_velocity_stack(GD_object, path_rows)

        #############################################################################################
        # This section is for stacking the files into a common grid
        output_data_stack(GD_object, vx_grids, vy_grids, v_grids, date_pairs, source_lists)
    else:
        message = '        TSX data stack has already been created'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)