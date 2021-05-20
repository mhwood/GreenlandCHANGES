from ....toolbox.series import series_to_N_points
from ....toolbox.reprojection import reproject_polygon
import numpy as np
import itertools
import requests
import os
import netCDF4 as nc4
from osgeo import gdal


########################################################################################################################

def get_multiyear_fileNames():
    file_names = ['greenland_vel_mosaic250_vx_v1.tif','greenland_vel_mosaic250_vy_v1.tif',
                  'greenland_vel_mosaic250_ex_v1.tif','greenland_vel_mosaic250_ey_v1.tif']
    return(file_names)

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

    date1 = fileID.split('_')[4]
    date2 = fileID.split('_')[5]

    year1 = int('20' + date1[-2:])
    month1 = monthStringToInt(date1[2:-2])
    day1 = int(date1[:2])

    year2 = int('20' + date2[-2:])
    month2 = monthStringToInt(date2[2:-2])
    day2 = int(date2[:2])

    datePair = str(year1) + "{:02d}".format(int(month1)) + "{:02d}".format(int(day1)) + '-' + str(
        year2) + "{:02d}".format(int(month2)) + "{:02d}".format(int(day2))
    return(datePair)

def create_velocity_stack(GD_object,measures_multiyear_mosaic_file_names):
    measures_multiyear_mosaic_folder = os.path.join(GD_object.data_folder, 'Velocity', 'MEaSUREs', 'Multi-year_Mosaic')
    X, Y = np.meshgrid(GD_object.velocity_grid_x, GD_object.velocity_grid_y)

    message = '            Looking for velocity points'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    vxSumGrid = np.zeros_like(X).astype(float)
    vySumGrid = np.zeros_like(X).astype(float)
    vSumGrid = np.zeros_like(X).astype(float)
    exSumGrid = np.zeros_like(X).astype(float)
    eySumGrid = np.zeros_like(X).astype(float)
    eSumGrid = np.zeros_like(X).astype(float)
    countGrid = np.zeros_like(X).astype(float)

    vx_file_name = measures_multiyear_mosaic_file_names[0]
    vy_file_name = measures_multiyear_mosaic_file_names[1]
    ex_file_name = measures_multiyear_mosaic_file_names[2]
    ey_file_name = measures_multiyear_mosaic_file_names[3]

    vx_file = os.path.join(GD_object.data_folder, 'Velocity', 'MEaSUREs','Multi-year_Mosaic', vx_file_name)
    ds = gdal.Open(vx_file)
    vx_array = np.array(ds.GetRasterBand(1).ReadAsArray())
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

    vy_file = os.path.join(GD_object.data_folder, 'Velocity', 'MEaSUREs','Multi-year_Mosaic', vy_file_name)
    ds = gdal.Open(vy_file)
    vy_array = np.array(ds.GetRasterBand(1).ReadAsArray())
    ds = None

    vy_array = vy_array[yIndices, :]
    vy_array = vy_array[:, xIndices]

    ex_file = os.path.join(GD_object.data_folder, 'Velocity', 'MEaSUREs', 'Multi-year_Mosaic', ex_file_name)
    ds = gdal.Open(ex_file)
    ex_array = np.array(ds.GetRasterBand(1).ReadAsArray())
    ds = None

    ex_array = ex_array[yIndices, :]
    ex_array = ex_array[:, xIndices]

    ey_file = os.path.join(GD_object.data_folder, 'Velocity', 'MEaSUREs', 'Multi-year_Mosaic', ey_file_name)
    ds = gdal.Open(ey_file)
    ey_array = np.array(ds.GetRasterBand(1).ReadAsArray())
    ds = None

    ey_array = ey_array[yIndices, :]
    ey_array = ey_array[:, xIndices]

    for xi in range(len(x)):
        xIndex = np.argmin(np.abs(GD_object.velocity_grid_x - x[xi]))
        for yi in range(len(y)):
            yIndex = np.argmin(np.abs(GD_object.velocity_grid_y - y[yi]))
            vx = vx_array[yi, xi]
            vy = vy_array[yi, xi]
            ex = ex_array[yi, xi]
            ey = ey_array[yi, xi]
            v = (vx**2 + vy**2)**0.5
            e = (ex ** 2 + ey ** 2) ** 0.5
            if v<1e6 and vx!=-99999:
                vxSumGrid[yIndex, xIndex] += vx
                vySumGrid[yIndex, xIndex] += vy
                vSumGrid[yIndex, xIndex] += v
                exSumGrid[yIndex, xIndex] += ex**2
                eySumGrid[yIndex, xIndex] += ey**2
                eSumGrid[yIndex, xIndex] += e**2
                countGrid[yIndex, xIndex] += 1

    VX = -99999.0 * np.ones_like(X)
    VY = -99999.0 * np.ones_like(X)
    V = -99999.0 * np.ones_like(X)
    VX[countGrid > 0] = vxSumGrid[countGrid > 0] / countGrid[countGrid > 0]
    VY[countGrid > 0] = vySumGrid[countGrid > 0] / countGrid[countGrid > 0]
    V[countGrid > 0] = vSumGrid[countGrid > 0] / countGrid[countGrid > 0]
    EX = -99999.0 * np.ones_like(X)
    EY = -99999.0 * np.ones_like(X)
    E = -99999.0 * np.ones_like(X)
    EX[countGrid > 0] = (exSumGrid[countGrid > 0])**0.5 / countGrid[countGrid > 0]
    EY[countGrid > 0] = (eySumGrid[countGrid > 0])**0.5 / countGrid[countGrid > 0]
    E[countGrid > 0] = (eSumGrid[countGrid > 0])**0.5 / countGrid[countGrid > 0]
    print('              Found '+str(np.sum(countGrid))+' points')

    #     if np.any(V > -99999):
    #         vx_grids.append(VX)
    #         vy_grids.append(VY)
    #         v_grids.append(V)
    #         ex_grids.append(EX)
    #         ey_grids.append(EY)
    #         e_grids.append(E)
    #         source_lists.append(sources[:-1])
    #         output_date_pairs.append(date_pair)
    return(VX, VY, V, EX, EY, E)

########################################################################################################################

def output_data_stac_as_grids(GD_object, vx_grids, vy_grids, v_grids,ex_grids, ey_grids, e_grids, date_pairs, source_lists):

    if len(GD_object.measures_multiyear_mosaic_output_file)>2:
        output_file = GD_object.measures_multiyear_mosaic_output_file
    else:
        output_file = os.path.join(GD_object.project_folder, GD_object.region_name, 'Velocity', 'Data',
                                   GD_object.region_name + ' MEaSUREs Multi-year Mosaic Velocity Grids.nc')

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

    message = '        Saved MEaSUREs Multi-year Mosaic file to ' + output_file
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

def output_data_stack_as_points(GD_object, VX, VY, V, EX, EY, E):

    ############################################################################
    # this first part converts the points to x,y,z,xi,yi columns grids
    message = '        Converting grids to points for efficient storage'
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)

    first_grid = VX
    xi = np.arange(np.shape(first_grid)[1])
    yi = np.arange(np.shape(first_grid)[0])
    XI, YI = np.meshgrid(xi, yi)
    X, Y = np.meshgrid(GD_object.velocity_grid_x, GD_object.velocity_grid_y)

    point_stack = np.hstack([np.reshape(X, (np.size(X), 1)),
                             np.reshape(Y, (np.size(Y), 1)),
                             np.reshape(VX, (np.size(VX), 1)),
                             np.reshape(VY, (np.size(VY), 1)),
                             np.reshape(V, (np.size(V), 1)),
                             np.reshape(EX, (np.size(EX), 1)),
                             np.reshape(EY, (np.size(EY), 1)),
                             np.reshape(E, (np.size(E), 1)),
                             np.reshape(XI, (np.size(XI), 1)),
                             np.reshape(YI, (np.size(YI), 1))])
    point_stack = point_stack[point_stack[:, 4] > -99, :]

    ############################################################################
    # this part saves all the point columns grids to an nc file

    if len(GD_object.measures_multiyear_mosaic_output_file) > 2:
        output_file = GD_object.measures_multiyear_mosaic_output_file
    else:
        output_file = os.path.join(GD_object.project_folder, GD_object.region_name, 'Velocity', 'Data',
                                   GD_object.region_name + ' MEaSUREs Multi-year Mosaic Velocity Points.nc')

    ds = nc4.Dataset(output_file, "w", format="NETCDF4")

    ds.createDimension("len_x", len(GD_object.velocity_grid_x))
    xvar = ds.createVariable("x", "f4", ("len_x",))
    xvar[:] = GD_object.velocity_grid_x

    ds.createDimension("len_y", len(GD_object.velocity_grid_y))
    yvar = ds.createVariable("y", "f4", ("len_y",))
    yvar[:] = GD_object.velocity_grid_y

    grp = ds.createGroup('mosaic')
    grp.createDimension("n_points", np.shape(point_stack)[0])

    xvar = grp.createVariable("x", "f4", ("n_points",))
    yvar = grp.createVariable("y", "f4", ("n_points",))
    vx_pointsvar = grp.createVariable("vx_points", "f4", ("n_points",))
    vy_pointsvar = grp.createVariable("vy_points", "f4", ("n_points",))
    v_pointsvar = grp.createVariable("v_points", "f4", ("n_points",))
    ex_pointsvar = grp.createVariable("ex_points", "f4", ("n_points",))
    ey_pointsvar = grp.createVariable("ey_points", "f4", ("n_points",))
    e_pointsvar = grp.createVariable("e_points", "f4", ("n_points",))
    xivar = grp.createVariable("xi", "f4", ("n_points",))
    yivar = grp.createVariable("yi", "f4", ("n_points",))

    xvar[:] = point_stack[:, 0]
    yvar[:] = point_stack[:, 1]
    vx_pointsvar[:] = point_stack[:, 2]
    vy_pointsvar[:] = point_stack[:, 3]
    v_pointsvar[:] = point_stack[:, 4]
    ex_pointsvar[:] = point_stack[:, 5]
    ey_pointsvar[:] = point_stack[:, 6]
    e_pointsvar[:] = point_stack[:, 7]
    xivar[:] = point_stack[:, 8]
    yivar[:] = point_stack[:, 9]


    ds.close()

    message = '        Outputting data to ' + output_file
    GD_object.output_summary += '\n' + message
    if GD_object.print_sub_outputs:
        print(message)


########################################################################################################################

def generate_MEaSUREs_Multiyear_Mosaic_dataset(GD_object,testing=False):

    if GD_object.overwrite_existing_measures_multiyear_mosaic_stack:
        stack_data = True
    else:
        if len(GD_object.measures_multiyear_mosaic_output_file) > 2:
            output_file = GD_object.measures_multiyear_mosaic_output_file
        else:
            output_file = os.path.join(GD_object.project_folder,GD_object.region_name,
                                       'Velocity', 'Data', GD_object.region_name + ' MEaSUREs Multi-year Mosaic Velocity Grids.nc')
        if os.path.isfile(output_file):
            stack_data = False
        else:
            stack_data = True

    if stack_data and GD_object.create_velocity_stacks:
        message = '        Stacking MEaSUREs Multi-year Mosaic data into a common grid'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)

        measures_multiyear_mosaic_file_names = get_multiyear_fileNames()

        #############################################################################################
        # This section is for stacking the files into a common grid
        VX, VY, V, EX, EY, E = create_velocity_stack(GD_object, measures_multiyear_mosaic_file_names)

        #############################################################################################
        # This section is for stacking the files into a common grid
        output_data_stack_as_points(GD_object, VX, VY, V, EX, EY, E)
    else:
        message = '        MEaSUREs Multi-year Mosaic data stack has already been created'
        GD_object.output_summary += '\n' + message
        if GD_object.print_sub_outputs:
            print(message)