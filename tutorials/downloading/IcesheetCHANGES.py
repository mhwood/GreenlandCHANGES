# this script defines a class used to hold all pertinent information to
# compile the elevation and velocity

import time
import datetime
import os

class IcesheetCHANGES:
    def __init__(self, project_folder, data_folder, region_name):

        #this initiates the global domain and all of its associated parameters
        self.project_folder = project_folder
        self.data_folder = data_folder
        self.region_name = region_name
        self.velocity_grid_x = []
        self.velocity_grid_y = []
        self.elevation_grid_y = []
        self.elevation_grid_x = []
        self.projection = 9999

class AntarcticCHANGES(IcesheetCHANGES):
    def __init__(self, project_folder, data_folder, region_name):
        super().__init__(project_folder, data_folder, region_name)
        self.projection = 3031

class GreenlandCHANGES(IcesheetCHANGES):
    def __init__(self, project_folder, data_folder, region_name):
        super().__init__(project_folder, data_folder, region_name)
        self.projection = 3413
