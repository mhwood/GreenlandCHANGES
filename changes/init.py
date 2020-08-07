

##############################################################################################
# Specify the directories

project_folder = '/Users/mhwood/Documents/Research/Projects/CHANGES/Examples/'
data_folder='/Volumes/mhwood/Research/Data Repository/Greenland'

##############################################################################################
# Create the class object to initiate the routine
import changes.initiation.GreenlandCHANGES as gc
GC = gc.GreenlandCHANGES(project_folder,data_folder)

##############################################################################################
# Set custom parameters
GC.set_extents_by_glacier('Helheim')


##############################################################################################
# Run routines to compile velocity and elevation data
GC.execute_velocity_and_elevation_compilations()


