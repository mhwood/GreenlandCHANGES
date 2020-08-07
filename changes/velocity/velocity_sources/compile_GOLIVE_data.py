
import os
import numpy as np
import ogr
import subprocess
import netCDF4 as nc4
from toolbox.reprojection import reproject_polygon

def obtain_list_of_overlapping_path_rows(GD_object):
    if GD_object.glacier_name + ' Landsat Path-Row List.csv' in os.listdir(os.path.join(GD_object.project_folder,
                                                                            GD_object.glacier_name, 'Velocity','Metadata')):
        f = open(os.path.join(GD_object.project_folder, GD_object.glacier_name, 'Velocity', 'Metadata',
                              GD_object.glacier_name + ' Landsat Path-Row List.csv'))
        lines = f.read()
        f.close()
        lines = lines.split('\n')
        lines.pop(0)
        path_rows = []
        for line in lines:
            path_rows.append(line)
    else:

        def numpyArrayToWKT(numpyArray):
            output = "POLYGON (("
            for b in range(np.shape(numpyArray)[0]):
                output += str(numpyArray[b, 0]) + " " + str(numpyArray[b, 1]) + ", "
            output = output[:-2] + "))"
            return (output)

        domainBoxArray = np.array([[GD_object.extents[0],GD_object.extents[1]],
                                   [GD_object.extents[2],GD_object.extents[1]],
                                   [GD_object.extents[2],GD_object.extents[3]],
                                   [GD_object.extents[0],GD_object.extents[3]],
                                   [GD_object.extents[0],GD_object.extents[1]]])
        domainWKT = numpyArrayToWKT(domainBoxArray)
        poly1 = ogr.CreateGeometryFromWkt(domainWKT)

        import changes.reference.landsat_path_row_domains as ld
        path_row_list = list(ld.path_row_string_to_polygon_outline.keys())

        path_rows = []
        for pr in range(len(path_row_list)):
            path_row = path_row_list[pr]
            path_row_polygon = np.array(ld.path_row_string_to_polygon_outline[path_row])
            landsatWKT = numpyArrayToWKT(path_row_polygon)
            poly2 = ogr.CreateGeometryFromWkt(landsatWKT)
            intersection = poly1.Intersection(poly2)
            intersection = intersection.ExportToWkt()
            if "EMPTY" not in intersection:
                path_rows.append(path_row)

        text_output = 'Path-Row'
        for ff in range(len(path_rows)):
            text_output += '\n' + path_rows[ff]

        f = open(os.path.join(GD_object.project_folder, GD_object.glacier_name, 'Velocity', 'Metadata',
                              GD_object.glacier_name + ' Landsat Path-Row List.csv'), 'w')
        f.write(text_output)
        f.close()

    return(path_rows)

def obtain_list_of_files_from_GOLIVE(path_row):
    urlHeader='ftp://dtn.rc.colorado.edu/work/nsidc0710/nsidc0710_landsat8_golive_ice_velocity_v1.1/'
    urlHeader+='p'+'{:03d}'.format(int(path_row[0]))+'_r'+'{:03d}'.format(int(path_row[1]))+'/'
    bashCommand='curl -b ~/.urs_cookies -c ~/.urs_cookies -l '+urlHeader
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    output=output.decode("utf-8")
    fileList=output.split('\n')
    files_list=[]
    links_list=[]
    for fil in fileList:
        if fil[-3:]=='.nc':
            files_list.append(fil)
            links_list.append(urlHeader+fil)

    return(files_list,links_list)

def obtain_list_of_tsx_files_overlapping_domain(GD_object,path_rows):
    if GD_object.glacier_name + ' GoLIVE Files.csv' in os.listdir(os.path.join(GD_object.project_folder,
                                                                            GD_object.glacier_name, 'Velocity','Metadata')):
        f = open(os.path.join(GD_object.project_folder, GD_object.glacier_name, 'Velocity', 'Metadata',
                              GD_object.glacier_name + ' GoLIVE Files.csv'))
        lines = f.read()
        f.close()
        lines = lines.split('\n')
        lines.pop(0)
        golive_file_links = []
        golive_file_names = []
        for line in lines:
            line = line.split(',')
            golive_file_names.append(line[0])
            golive_file_links.append(line[1])
    else:
        golive_file_links = []
        golive_file_names = []

        for pr in range(len(path_rows)):
            message = '            Gathering velocity files for path-row ' + str(path_rows[pr])
            GD_object.output_summary += '\n' + message
            if GD_object.print_sub_outputs:
                print(message)

            path_row = path_rows[pr].split('-')
            golive_files_for_path_row, golive_links_for_path_row = obtain_list_of_files_from_GOLIVE(path_row)

            message = '                Identified ' + str(len(golive_files_for_path_row)) + ' available files for '+ str(path_rows[pr])
            GD_object.output_summary += '\n' + message
            if GD_object.print_sub_outputs:
                print(message)

            golive_file_links+=golive_links_for_path_row
            golive_file_names += golive_files_for_path_row

        text_output = 'File_Name,URL'
        for ff in range(len(golive_file_names)):
            text_output += '\n' + golive_file_names[ff] + ',' + golive_file_links[ff]

        f = open(os.path.join(GD_object.project_folder, GD_object.glacier_name, 'Velocity', 'Metadata',
                              GD_object.glacier_name + ' GoLIVE Files.csv'), 'w')
        f.write(text_output)
        f.close()

    return(golive_file_links, golive_file_names)

def obtain_download_list(GD_object,golive_file_names):
    # if 'Path '+str(path_row[0])+' Row '+str(path_row[1]) not in os.listdir(os.join(self.)):
    #     os.mkdir(goLiveFolder+'/Path '+str(pathRow[0])+' Row '+str(pathRow[1]))
    golive_folder=os.path.join(GD_object.data_folder,'Velocity','GOLIVE')
    download_list=[]
    for fil in golive_file_names:
        path_row_str = fil[3:10]
        if path_row_str not in os.listdir(golive_folder):
            download_list.append(fil)
        else:
            if fil not in os.listdir(os.path.join(golive_folder,path_row_str)):
                download_list.append(fil)
    n_existing_files = len(golive_file_names)-len(download_list)
    return(download_list,n_existing_files)

def download_golive_files(GD_object,download_list,golive_file_links, golive_file_names):

    golive_folder = os.path.join(GD_object.data_folder, 'Velocity', 'GOLIVE')

    for ff in range(len(download_list)):
        pwd = os.getcwd()
        fil = download_list[ff]
        link = golive_file_links[golive_file_names.index(fil)]
        path_row_str = fil[3:10]
        subFolder = os.path.join(golive_folder,path_row_str)
        os.chdir(subFolder)
        if GD_object.print_sub_outputs:
            message = '                    Downloading ' + fil + ' (' + str(ff + 1) + ' of ' + str(len(download_list)) + ')'
            print(message)
        bashCommand='curl -b ~/.urs_cookies -c ~/.urs_cookies -L -n -O '+link
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
        if GD_object.print_sub_outputs:
            if 'RETR response: 550' in error.decode('utf-8'):
                print('                        Download failed (RETR response: 550)')
        os.chdir(pwd)

def find_date_pairs_in_path_rows(GD_object, path_rows):
    date_pairs=[]
    file_sets=[]
    path_row_sets=[]
    for path_row in path_rows:
        path = int(path_row.split('-')[0])
        row = int(path_row.split('-')[1])
        path_row_str = '{:03d}'.format(path)+'_'+'{:03d}'.format(row)
        if path_row_str in os.listdir(os.path.join(GD_object.data_folder,'Velocity','GOLIVE')):
            pathRowFolder=os.path.join(GD_object.data_folder,'Velocity','GOLIVE',path_row_str)
            checkFiles = os.listdir(pathRowFolder)
            for fil in checkFiles:
                if fil[-3:] == '.nc':
                    dateString=fil[15:32]
                    if dateString not in date_pairs:
                        date_pairs.append(dateString)
                        file_sets.append([fil])
                        path_row_sets.append([path_row_str])
                    else:
                        file_sets[date_pairs.index(dateString)].append(fil)
                        path_row_sets[date_pairs.index(dateString)].append(path_row_str)
    return(date_pairs,file_sets,path_row_sets)

def readFileToGlacierGrids(GD_object, path_row, fil):

    ncFile = os.path.join(GD_object.data_folder,'Velocity','GOLIVE',path_row,fil)
    data=nc4.Dataset(ncFile)
    times=str(data.variables['image_pair_times']).split('\n')
    startDate=''
    endDate=''
    for time in times:
        if 'start_date' in time:
            startDate=time.split(': ')[1][:10].replace('-','')
        if 'end_date' in time:
            endDate = time.split(': ')[1][:10].replace('-', '')
    outputString=startDate+'-'+endDate

    metadataLines = str(data.variables['input_image_details'])
    espg=metadataLines[metadataLines.index('"EPSG","326')+8:metadataLines.index('"EPSG","326')+13]

    X, Y = np.meshgrid(GD_object.velocity_grid_x,GD_object.velocity_grid_y)
    points=np.hstack([np.reshape(X,(np.size(X),1)),
                      np.reshape(Y, (np.size(Y), 1))])
    points=reproject_polygon(points,3413,int(espg))
    VX=-99999*np.ones((len(points),1))
    VY = -99999 * np.ones((len(points), 1))
    V = -99999 * np.ones((len(points), 1))

    x=data.variables['x'][:]
    y=data.variables['y'][:]
    vx=data.variables['vx_masked'][:,:]
    vy=data.variables['vy_masked'][:,:]
    vv=data.variables['vv_masked'][:,:]
    for pp in range(len(points)):
        xi = np.argmin(np.abs(x-points[pp,0]))
        yi = np.argmin(np.abs(y-points[pp,1]))
        dist=(x[xi]-points[pp,0])**2 + (y[yi]-points[pp,1])**2
        if dist<2*150**2:
            if vv[yi,xi]>0:
                VX[pp]=vx[yi,xi]
                VY[pp]=vy[yi,xi]
                V[pp]=vv[yi,xi]

    VX=np.reshape(VX,np.shape(X))
    VY=np.reshape(VY,np.shape(Y))
    V=np.reshape(V,np.shape(X))

    #this is to convert from m/d to m/yr
    VX=VX*365.25
    VY=VY*365.25
    V=V*365.25

    return(VX,VY,V,outputString)

def read_date_files_to_glacier_grid(GD_object, path_row_set,file_list):
    X, Y = np.meshgrid(GD_object.velocity_grid_x,GD_object.velocity_grid_y)
    vxSumGrid = np.zeros_like(X).astype(float)
    vySumGrid = np.zeros_like(X).astype(float)
    vSumGrid = np.zeros_like(X).astype(float)
    countGrid = np.zeros_like(X).astype(float)
    source_file_list=''
    for ff in range(len(file_list)):
        path_row=path_row_set[ff]
        fil=file_list[ff]
        source_file_list+=fil+','
        VX, VY, V, outputString = readFileToGlacierGrids(GD_object, path_row, fil)
        goodIndices=V>0
        vxSumGrid[goodIndices]+=VX[goodIndices]
        vySumGrid[goodIndices]+=VY[goodIndices]
        vSumGrid[goodIndices]+=V[goodIndices]
        countGrid[goodIndices]+=1
    VX = -99999.0*np.ones_like(X)
    VY = -99999.0 * np.ones_like(X)
    V = -99999.0 * np.ones_like(X)
    VX[countGrid>0]=vxSumGrid[countGrid>0]/countGrid[countGrid>0]
    VY[countGrid > 0] = vySumGrid[countGrid > 0] / countGrid[countGrid > 0]
    V[countGrid > 0] = vSumGrid[countGrid > 0] / countGrid[countGrid > 0]
    return(VX,VY,V,outputString,source_file_list[:-1])

def create_velocity_stack(GD_object, path_rows):

    # find date pairs in path-row set
    date_pairs, file_sets, path_row_sets = find_date_pairs_in_path_rows(GD_object, path_rows)

    # read the velocity field from each file to a stack for each date pair
    vx_grids=[]
    vy_grids=[]
    v_grids=[]
    date_strings=[]
    source_lists=[]
    for dd in range(len(date_pairs)):
        if GD_object.print_sub_outputs:
            print('            Adding data for date pair '+date_pairs[dd]+' ('+str(dd+1)+' of '+str(len(date_pairs))+')')
        file_list=file_sets[dd]
        path_row_set=path_row_sets[dd]
        VX,VY,V,output_date_string,source_files_string = read_date_files_to_glacier_grid(GD_object, path_row_set, file_list)
        vx_grids.append(VX)
        vy_grids.append(VY)
        v_grids.append(V)
        date_strings.append(output_date_string)
        source_lists.append(source_files_string)

    return(vx_grids, vy_grids, v_grids, date_strings, source_lists)

def output_data_stack(GD_object,vx_grids, vy_grids, v_grids, date_pairs, source_lists):
    output_file=os.path.join(GD_object.project_folder,'Data','Glaciers',GD_object.glacier_name,'Velocity',GD_object.glacier_name+' GOLIVE Velocity Data.nc')

    data=nc4.Dataset(output_file, "w", format="NETCDF4")

    data.createDimension('y',len(GD_object.velocity_grid_y))
    data.createDimension('x',len(GD_object.velocity_grid_x))
    xvar = data.createVariable('x','f4',("x",))
    yvar = data.createVariable('y','f4',("y",))

    xvar[:]=GD_object.velocity_grid_x
    yvar[:]=GD_object.velocity_grid_y

    for dd in range(len(date_pairs)):
        grp = data.createGroup(date_pairs[dd])
        vx = grp.createVariable('VX','f4',("y","x"))
        vy = grp.createVariable('VY','f4',("y","x"))
        v = grp.createVariable('V', 'f4', ("y", "x"))

        vx[:,:] = vx_grids[dd]
        vy[:,:] = vy_grids[dd]
        v[:,:] = v_grids[dd]
        grp.source_files = source_lists[dd]

    data.close()

    message = '        Saved GOLIVE file to ' + output_file
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)


def generate_GOLIVE_dataset(GD_object,testing=False):

    #############################################################################################
    # This section is for downloading all of the files overlapping the domain

    message = '    Running compilation for the GOLIVE data'
    GD_object.output_summary += '\n' + message
    if GD_object.print_main_outputs:
        print(message)

    path_rows = obtain_list_of_overlapping_path_rows(GD_object)

    message = '        Found '+str(len(path_rows))+' path-rows overlapping the domain: '+', '.join(path_rows)
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    golive_file_links, golive_file_names = obtain_list_of_tsx_files_overlapping_domain(GD_object,path_rows)

    message = '        Found ' + str(len(golive_file_names)) + ' total files in these path-rows'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    download_list, n_existing_files = obtain_download_list(GD_object, golive_file_names)
    message = '            Found ' + str(n_existing_files) + ' files already downloaded'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    if GD_object.download_new_golive_data:
        if GD_object.overwrite_existing_velocity_data:
            download_list = golive_file_names
            message = '            Overwriting old files by request'
            GD_object.output_summary += '\n' + message
            if GD_object.print_sub_outputs:
                print(message)

        message = '                Downloading ' + str(len(download_list)) + ' files'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

            download_golive_files(GD_object, download_list,golive_file_links, golive_file_names)
    else:
        message = '        Skipping the download of new files by request'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

    if GD_object.overwrite_existing_golive_stack:
        stack_data = True
    else:
        if GD_object.glacier_name + ' GoLIVE Velocity Data.nc' in os.listdir(os.path.join(GD_object.project_folder,
                                                                                         GD_object.glacier_name,
                                                                                         'Velocity', 'Data')):
            stack_data = False
        else:
            stack_data = True

    if stack_data:
        message = '        Stacking GoLIVE data into a common grid'
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
        message = '        GoLIVE data stack has already been created'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)