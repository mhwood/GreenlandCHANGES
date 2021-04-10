
# this function is used to interpolate a coarser resolution grid onto more dense grid

import numpy as np
import matplotlib.pyplot as plt
from ..series import series_to_N_points
from ..reprojection import reproject_polygon
from scipy.interpolate import interp2d
import matplotlib.path as mplPath

def reproject_and_interpolate_onto_grid(source_x,source_y,source_grid,source_epsg,
                                        dest_x,dest_y,dest_epsg,
                                        print_status_messages=False, interpolation_type = 'linear'):

    # plt.contourf(source_x,source_y,source_grid)
    # plt.show()

    if print_status_messages:
        print('                    Creating the interpolation object')
    set_int = interp2d(source_x,source_y,source_grid, kind=interpolation_type)

    X, Y = np.meshgrid(dest_x,dest_y)
    points = np.hstack([np.reshape(X, (np.size(X), 1)),
                        np.reshape(Y, (np.size(Y), 1))]).astype(float)

    if source_epsg!=dest_epsg:
        if print_status_messages:
            print('                    Reprojecting the points to ' + str(source_epsg)+' for interpolation')
        points = reproject_polygon(points, dest_epsg, source_epsg)

    if print_status_messages:
        print('                    Sampling the interpolation object on the grid')
    X = np.reshape(points[:, 0], np.shape(X))
    Y = np.reshape(points[:, 1], np.shape(X))
    grid = -99*np.ones_like(X)

    # spell this path out to make sure corners are preserved
    southern_line = np.array([[np.min(source_x),np.min(source_y)],
                            [np.max(source_x),np.min(source_y)]])
    southern_line = series_to_N_points(southern_line, 25)
    eastern_line = np.array([[np.max(source_x), np.min(source_y)],
                            [np.max(source_x), np.max(source_y)]])
    eastern_line = series_to_N_points(eastern_line, 25)
    northern_line = np.array([[np.max(source_x), np.max(source_y)],
                            [np.min(source_x), np.max(source_y)]])
    northern_line = series_to_N_points(northern_line, 25)
    western_line = np.array([[np.min(source_x), np.max(source_y)],
                            [np.min(source_x), np.min(source_y)]])
    western_line = series_to_N_points(western_line, 25)
    bbox = np.vstack([southern_line,
                     eastern_line,
                     northern_line,
                     western_line])

    # plt.plot(points[:,0],points[:,1],'k.')
    # plt.plot(bbox[:,0],bbox[:,1],'r-')
    # plt.show()

    source_path = mplPath.Path(bbox)

    for yi in range(np.shape(X)[0]):
        if print_status_messages:
            if yi == int(0.25*np.shape(X)[0]):
                print('                      25% complete...')
            if yi == int(0.5*np.shape(X)[0]):
                print('                      50% complete...')
            if yi == int(0.75*np.shape(X)[0]):
                print('                      75% complete...')
        for xi in range(np.shape(X)[1]):
            # make sure the point is inside the boundary
            # (interpolation seems to give "valid" values outside due to reprojection)
            # is_inside = source_path.contains_point((X[yi,xi],Y[yi,xi]))
            # int_val = set_int(X[yi, xi], Y[yi, xi])
            if source_path.contains_point((X[yi,xi],Y[yi,xi])):
                int_val = set_int(X[yi,xi],Y[yi,xi])
                if not np.isnan(int_val):
                    grid[yi,xi] = int_val

    return(grid)

