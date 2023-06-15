# this script defines a class used to hold all pertinent information to
# compile the elevation and velocity

import time
import datetime
import os

class IcesheetCHANGES:
    def __init__(self, project_folder, data_folder, collection_id, region_name, projection):

        #this initiates the global domain and all of its associated parameters
        self.project_folder = project_folder
        self.data_folder = data_folder
        self.collection_id = collection_id
        self.region_name = region_name
        self.projection = projection
        self.velocity_grid_x = []
        self.velocity_grid_y = []
        self.elevation_grid_y = []
        self.elevation_grid_x = []


    # Projections are defined by EPSG codes (https://epsg.io/)
    def set_projection(self, projection):
        self.projection = projection
    def get_projection(self):
        return self.projection

class AntarcticCHANGES(IcesheetCHANGES):
    def __init__(self, project_folder, data_folder, collection_id, region_name, projection = 3031):
        super().__init__(project_folder, data_folder, collection_id, region_name, projection)
        self.velocity_grid_x = []

class GreenlandCHANGES(IcesheetCHANGES):
    def __init__(self, project_folder, data_folder, collection_id, region_name, projection = 3413):
        super().__init__(project_folder, data_folder, collection_id, region_name, projection)
        self.velocity_grid_x = []


