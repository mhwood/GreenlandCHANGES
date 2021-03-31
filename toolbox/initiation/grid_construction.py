
import numpy as np

def create_grid_from_extents(extents,posting,epsg):

    # the grid is created using the grid frame from bedmachine
    minBedMachineX = -652925
    maxBedMachineX = 879625
    minBedMachineY = -3384425
    maxBedMachineY = -632675

    x = np.arange(minBedMachineX, maxBedMachineX, posting)
    y = np.arange(minBedMachineY, maxBedMachineY, posting)

    x = x[x >= extents[0]]
    x = x[x <= extents[2]]

    y = y[y >= extents[1]]
    y = y[y <= extents[3]]
    return(x,y)

def create_grids_from_extents(GD_object):
    velocity_grid_x, velocity_grid_y = \
        create_grid_from_extents(GD_object.extents, GD_object.velocity_grid_posting, GD_object.velocity_grid_epsg)
    elevation_grid_x, elevation_grid_y = \
        create_grid_from_extents(GD_object.extents, GD_object.elevation_grid_posting, GD_object.elevation_grid_epsg)

    GD_object.velocity_grid_x = velocity_grid_x
    GD_object.velocity_grid_y = velocity_grid_y
    GD_object.elevation_grid_y = elevation_grid_y
    GD_object.elevation_grid_x = elevation_grid_x