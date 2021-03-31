from ....toolbox.series import series_to_N_points
from ....toolbox.reprojection import reproject_polygon
import numpy as np
import itertools
import requests
import os
import netCDF4 as nc4
import gdal
import matplotlib.pyplot as plt

def obtain_list_of_measures_optical_files_overlapping_domain(GD_object):

    if GD_object.region_name+' MEaSUREs Optical Files.csv' in os.listdir(os.path.join(GD_object.project_folder,
                                                                                    GD_object.region_name, 'Velocity', 'Metadata')):
        f=open(os.path.join(GD_object.project_folder, GD_object.region_name, 'Velocity', 'Metadata',
                                      GD_object.region_name+' MEaSUREs Optical Files.csv'))
        lines=f.read()
        f.close()
        lines=lines.split('\n')
        lines.pop(0)
        measures_optical_file_links = []
        measures_optical_file_names = []
        for line in lines:
            line=line.split(',')
            measures_optical_file_names.append(line[0])
            measures_optical_file_links.append(line[1])
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
        short_name = 'NSIDC-0646'
        version = '3'
        time_start = '2010-01-01T00:00:00Z'
        time_end = '2100-01-01T23:59:59Z'
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
        # print(urls)

        measures_optical_file_links = []
        measures_optical_file_names = []
        for url in urls:
            if 'vx_v0'+version+'.0.tif' in url or 'vy_v0'+version+'.0.tif' in url or 'ex_v0'+version+'.0.tif' in url or 'ey_v0'+version+'.0.tif' in url:
                measures_optical_file_links.append(url)
                url_parts = url.split('/')
                file_name = url_parts[-1]
                measures_optical_file_names.append(file_name)

        text_output = 'File_Name,URL'
        for ff in range(len(measures_optical_file_names)):
            text_output+='\n'+measures_optical_file_names[ff]+','+measures_optical_file_links[ff]

        f = open(os.path.join(GD_object.project_folder, GD_object.region_name, 'Velocity', 'Metadata',
                                        GD_object.region_name + ' MEaSUREs Optical Files.csv'),'w')
        f.write(text_output)
        f.close()

    return(measures_optical_file_links, measures_optical_file_names)

def obtain_download_list(GD_object,measures_optical_file_names):
    measures_folder=os.path.join(GD_object.data_folder,'Velocity','MEaSUREs','Optical','Data')
    download_list=[]
    for fil in measures_optical_file_names:
        location_id = fil.split('_')[1]
        if location_id not in os.listdir(measures_folder):
            download_list.append(fil)
        else:
            if fil not in os.listdir(os.path.join(measures_folder,location_id)):
                download_list.append(fil)
    n_existing_files = len(measures_optical_file_names)-len(download_list)
    return(download_list,n_existing_files)

def download_measures_optical_files(GD_object,measures_optical_file_links, measures_optical_file_names,download_list):

    for ff in range(len(measures_optical_file_names)):
        url = measures_optical_file_links[ff]
        file_name = measures_optical_file_names[ff]
        if file_name in download_list:
            if GD_object.print_sub_outputs:
                print('                Downloading '+file_name+' ('+str(download_list.index(file_name)+1)+' of '+str(len(download_list))+')')
            location_id = file_name.split('_')[1]
            if location_id not in os.listdir(os.path.join(GD_object.data_folder,'Velocity','MEaSUREs','Optical','Data')):
                os.mkdir(os.path.join(GD_object.data_folder,'Velocity','MEaSUREs','Optical','Data',location_id))
            output_file = os.path.join(GD_object.data_folder,'Velocity','MEaSUREs','Optical','Data',location_id,file_name)
            try:
                resp = requests.get(url)
                open(output_file, 'wb').write(resp.content)
            except:
                print('ISSUE WITH FILE '+file_name)

def create_velocity_stack(GD_object,measures_optical_file_names):
    measures_optical_folder = os.path.join(GD_object.data_folder, 'Velocity', 'MEaSUREs', 'Optical','Data')
    X, Y = np.meshgrid(GD_object.velocity_grid_x, GD_object.velocity_grid_y)

    #find list of file IDs which contain both x and y grids
    file_pairs=[]
    fileIDs = []
    for file_name in measures_optical_file_names:
        fileID=file_name[:19]
        location_id = file_name.split('_')[1]
        if fileID not in fileIDs:
            vx_found = False
            vy_found = False
            for file_name_check in os.listdir(os.path.join(measures_optical_folder, location_id)):
                if fileID + '_vx' in file_name_check:
                    vx_file = file_name_check
                    vx_found = True
                if fileID + '_vy' in file_name_check:
                    vy_file = file_name_check
                    vy_found = True
                if fileID + '_ex' in file_name_check:
                    ex_file = file_name_check
                    ex_found = True
                if fileID + '_ey' in file_name_check:
                    ey_file = file_name_check
                    ey_found = True
            if vx_found and vy_found and ex_found and ey_found:
                fileIDs.append(fileID)
                file_pairs.append([vx_file,vy_file,ex_file,ey_file])

    dates = []
    fileID_sets = []
    file_pairs_sets = []
    for fi in range(len(fileIDs)):
        date = fileIDs[fi].split('_')[2]
        if date not in dates:
            dates.append(date)
            fileID_sets.append([fileIDs[fi]])
            file_pairs_sets.append([file_pairs[fi]])
        else:
            fileID_sets[dates.index(date)].append(fileIDs[fi])
            file_pairs_sets[dates.index(date)].append(file_pairs[fi])

    vx_grids=[]
    vy_grids=[]
    ex_grids = []
    ey_grids = []
    v_grids = []
    e_grids = []
    source_lists=[]

    for dd in range(len(dates)):
        date=dates[dd]

        message = '            Working on date ' + dates[dd] + ' (date ' + str(dd + 1) + ' of ' + str(len(dates)) + ')'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

        fileIDs = fileID_sets[dd]
        file_pairs=file_pairs_sets[dd]
        sources=''

        vxSumGrid = np.zeros_like(X).astype(float)
        vySumGrid = np.zeros_like(X).astype(float)
        vSumGrid = np.zeros_like(X).astype(float)
        exSumGrid = np.zeros_like(X).astype(float)
        eySumGrid = np.zeros_like(X).astype(float)
        eSumGrid = np.zeros_like(X).astype(float)
        countGrid = np.zeros_like(X).astype(float)

        for ff in range(len(fileIDs)):
            fileID=fileIDs[ff]
            location_id = fileID.split('_')[1]
            vx_file_name = file_pairs[ff][0]
            vy_file_name = file_pairs[ff][1]
            ex_file_name = file_pairs[ff][2]
            ey_file_name = file_pairs[ff][3]
            sources += vx_file_name+','+vy_file_name
            message = '              Looking for velocity points in ' + fileID + ' (file '+str(ff+1)+' of '+str(len(fileIDs))+')'
            GD_object.output_summary += '\n' + message
            if GD_object.print_sub_outputs:
                print(message)

            vx_file = os.path.join(GD_object.data_folder, 'Velocity', 'MEaSUREs','Optical','Data', location_id, vx_file_name)
            ds = gdal.Open(vx_file)
            vx_array = np.array(ds.GetRasterBand(1).ReadAsArray())
            transform = ds.GetGeoTransform()
            ds = None

            vy_file = os.path.join(GD_object.data_folder, 'Velocity', 'MEaSUREs','Optical','Data', location_id, vy_file_name)
            ds = gdal.Open(vy_file)
            vy_array = np.array(ds.GetRasterBand(1).ReadAsArray())
            transform = ds.GetGeoTransform()
            ds = None

            ex_file = os.path.join(GD_object.data_folder, 'Velocity', 'MEaSUREs', 'Optical', 'Data', location_id,
                                   ex_file_name)
            ds = gdal.Open(ex_file)
            ex_array = np.array(ds.GetRasterBand(1).ReadAsArray())
            transform = ds.GetGeoTransform()
            ds = None

            ey_file = os.path.join(GD_object.data_folder, 'Velocity', 'MEaSUREs', 'Optical', 'Data', location_id,
                                   ey_file_name)
            ds = gdal.Open(ey_file)
            ey_array = np.array(ds.GetRasterBand(1).ReadAsArray())
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
            ex_array = ex_array[yIndices, :]
            ex_array = ex_array[:, xIndices]
            ey_array = ey_array[yIndices, :]
            ey_array = ey_array[:, xIndices]

            for xi in range(len(x)):
                xIndex = np.argmin(np.abs(GD_object.velocity_grid_x - x[xi]))
                for yi in range(len(y)):
                    yIndex = np.argmin(np.abs(GD_object.velocity_grid_y - y[yi]))
                    vx = vx_array[yi,xi]
                    vy = vy_array[yi,xi]
                    ex = ex_array[yi, xi]
                    ey = ey_array[yi, xi]
                    v = (vx ** 2 + vy ** 2) ** 0.5
                    e = (ex ** 2 + ey ** 2) ** 0.5
                    if v<1e6 and vx!=-99999:
                        vxSumGrid[yIndex, xIndex] += vx
                        vySumGrid[yIndex, xIndex] += vy
                        vSumGrid[yIndex, xIndex] += v
                        exSumGrid[yIndex, xIndex] += ex ** 2
                        eySumGrid[yIndex, xIndex] += ey ** 2
                        eSumGrid[yIndex, xIndex] += e ** 2
                        countGrid[yIndex, xIndex] += 1

        VX = -99999.0 * np.ones_like(X)
        VY = -99999.0 * np.ones_like(X)
        V = -99999.0 * np.ones_like(X)
        EX = -99999.0 * np.ones_like(X)
        EY = -99999.0 * np.ones_like(X)
        E = -99999.0 * np.ones_like(X)

        VX[countGrid > 0] = vxSumGrid[countGrid > 0] / countGrid[countGrid > 0]
        VY[countGrid > 0] = vySumGrid[countGrid > 0] / countGrid[countGrid > 0]
        V[countGrid > 0] = vSumGrid[countGrid > 0] / countGrid[countGrid > 0]
        EX[countGrid > 0] = (exSumGrid[countGrid > 0]) ** 0.5 / countGrid[countGrid > 0]
        EY[countGrid > 0] = (eySumGrid[countGrid > 0]) ** 0.5 / countGrid[countGrid > 0]
        E[countGrid > 0] = (eSumGrid[countGrid > 0]) ** 0.5 / countGrid[countGrid > 0]

        # plt.subplot(1,2,1)
        # V_copy = np.copy(V)
        # V_copy[V_copy == -99999] = 0
        # C = plt.contourf(V_copy)
        # plt.colorbar(C)
        # plt.subplot(1,2,2)
        # C2 = plt.contourf(countGrid)
        # plt.colorbar(C2)
        # plt.show()

        # if np.any(V>-9999):
        vx_grids.append(VX)
        vy_grids.append(VY)
        v_grids.append(V)
        ex_grids.append(EX)
        ey_grids.append(EY)
        e_grids.append(E)
        source_lists.append(sources)
        # else:
        #     dates.pop(dates.index(date))
    return(vx_grids, vy_grids, v_grids,ex_grids, ey_grids, e_grids, dates, source_lists)

def output_data_stack(GD_object, vx_grids, vy_grids, v_grids,ex_grids, ey_grids, e_grids, dates, source_lists):
    output_file = os.path.join(GD_object.project_folder, GD_object.region_name, 'Velocity', 'Data',
                               GD_object.region_name + ' MEaSUREs Optical Velocity Grids.nc')

    data = nc4.Dataset(output_file, "w", format="NETCDF4")

    data.createDimension('y', len(GD_object.velocity_grid_y))
    data.createDimension('x', len(GD_object.velocity_grid_x))
    xvar = data.createVariable('x', 'f4', ("x",))
    yvar = data.createVariable('y', 'f4', ("y",))

    xvar[:] = GD_object.velocity_grid_x
    yvar[:] = GD_object.velocity_grid_y

    for dd in range(len(dates)):
        grp = data.createGroup(dates[dd])
        vx = grp.createVariable('VX', 'f4', ("y", "x"))
        vy = grp.createVariable('VY', 'f4', ("y", "x"))
        v = grp.createVariable('V', 'f4', ("y", "x"))
        ex = grp.createVariable('EX', 'f4', ("y", "x"))
        ey = grp.createVariable('EY', 'f4', ("y", "x"))
        e = grp.createVariable('E', 'f4', ("y", "x"))

        vx[:, :] = vx_grids[dd]
        vy[:, :] = vy_grids[dd]
        v[:, :] = v_grids[dd]
        ex[:, :] = ex_grids[dd]
        ey[:, :] = ey_grids[dd]
        e[:, :] = e_grids[dd]
        grp.source_files = source_lists[dd]

    data.close()

    message = '        Saved MEaSUREs Optical file to ' + output_file
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)


def generate_MEaSUREs_Optical_dataset(GD_object,testing=False):

    if 'MEaSUREs' not in os.listdir(os.path.join(GD_object.data_folder,'Velocity')):
        os.mkdir(os.path.join(GD_object.data_folder,'Velocity','MEaSUREs'))

    if 'Optical' not in os.listdir(os.path.join(GD_object.data_folder,'Velocity','MEaSUREs')):
        os.mkdir(os.path.join(GD_object.data_folder,'Velocity','MEaSUREs','Optical'))

    if 'Data' not in os.listdir(os.path.join(GD_object.data_folder,'Velocity','MEaSUREs','Optical')):
        os.mkdir(os.path.join(GD_object.data_folder,'Velocity','MEaSUREs','Optical','Data'))

    if 'Metadata' not in os.listdir(os.path.join(GD_object.data_folder,'Velocity','MEaSUREs','Optical')):
        os.mkdir(os.path.join(GD_object.data_folder,'Velocity','MEaSUREs','Optical','Metadata'))

    #############################################################################################
    # This section is for downloading all of the files overlapping the domain

    message = '    Running compilation for the MEaSUREs Optical data'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)

    measures_optical_file_links, measures_optical_file_names = obtain_list_of_measures_optical_files_overlapping_domain(GD_object)

    message = '        Found ' + str(len(measures_optical_file_names)) + ' MEaSUREs Optical files overlapping the domain'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    if isinstance(GD_object.max_number_of_measures_optical_files, int):
        GD_object.max_number_of_measures_optical_files = 4 * GD_object.max_number_of_measures_optical_files
        if len(measures_optical_file_links) > GD_object.max_number_of_measures_optical_files:
            measures_optical_file_links = measures_optical_file_links[:GD_object.max_number_of_measures_optical_files]
            measures_optical_file_names = measures_optical_file_names[:GD_object.max_number_of_measures_optical_files]


    if GD_object.download_new_measures_optical_data:

        if GD_object.overwrite_existing_velocity_data:
            download_list = measures_optical_file_names
            message = '            Overwriting old files by request - downloading ' + str(len(download_list)) + ' files'
            GD_object.output_summary += '\n' + message
            if GD_object.print_sub_outputs:
                print(message)
        else:
            download_list, n_existing_files = obtain_download_list(GD_object, measures_optical_file_names)
            message = '            Found ' + str(n_existing_files) + ' existing files - downloading ' + str(len(download_list)) + ' files'
            GD_object.output_summary += '\n' + message
            if GD_object.print_sub_outputs:
                print(message)

        download_measures_optical_files(GD_object,measures_optical_file_links, measures_optical_file_names,download_list)

    else:
        message = '        Skipping the download of new files by request'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

    if GD_object.overwrite_existing_measures_optical_stack:
        stack_data = True
    else:
        if GD_object.region_name + ' MEaSUREs Optical Velocity Grids.nc' in os.listdir(os.path.join(GD_object.project_folder,
                                                                                         GD_object.region_name,
                                                                                         'Velocity', 'Data')):
            stack_data = False
        else:
            stack_data = True

    if stack_data and GD_object.create_velocity_stacks:
        message = '        Stacking MEaSUREs Optical data into a common grid'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)
        #############################################################################################
        # This section is for stacking the files into a common grid

        vx_grids, vy_grids, v_grids,ex_grids, ey_grids, e_grids, dates, source_lists = create_velocity_stack(GD_object, measures_optical_file_names)

        #############################################################################################
        # This section is for stacking the files into a common grid
        output_data_stack(GD_object, vx_grids, vy_grids, v_grids,ex_grids, ey_grids, e_grids, dates, source_lists)
    else:
        message = '        MEaSUREs Optical data stack has already been created'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)