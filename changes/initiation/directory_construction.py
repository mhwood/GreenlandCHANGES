
import os

def construct_output_directory_structure(GC_object):

    ##########################################################################################################
    # these directories are constructed in the project_folder, where the final data will be stored

    #make the region folder
    if GC_object.region_name not in os.listdir(GC_object.project_folder):
        os.mkdir(os.path.join(GC_object.project_folder,GC_object.region_name))

    # make the velocity folder
    if GC_object.compile_velocity:
        if 'Velocity' not in os.listdir(os.path.join(GC_object.project_folder, GC_object.region_name)):
            os.mkdir(os.path.join(GC_object.project_folder, GC_object.region_name, 'Velocity'))

        if 'Data' not in os.listdir(os.path.join(GC_object.project_folder, GC_object.region_name, 'Velocity')):
            os.mkdir(os.path.join(GC_object.project_folder, GC_object.region_name, 'Velocity', 'Data'))

        if 'Metadata' not in os.listdir(os.path.join(GC_object.project_folder, GC_object.region_name, 'Velocity')):
            os.mkdir(os.path.join(GC_object.project_folder, GC_object.region_name, 'Velocity', 'Metadata'))

    if GC_object.compile_elevation:
        # make the elevation folder
        if 'Elevation' not in os.listdir(os.path.join(GC_object.project_folder, GC_object.region_name)):
            os.mkdir(os.path.join(GC_object.project_folder, GC_object.region_name,'Elevation'))

        if 'Data' not in os.listdir(os.path.join(GC_object.project_folder, GC_object.region_name,'Elevation')):
            os.mkdir(os.path.join(GC_object.project_folder, GC_object.region_name, 'Elevation','Data'))

        if 'Metadata' not in os.listdir(os.path.join(GC_object.project_folder, GC_object.region_name,'Elevation')):
            os.mkdir(os.path.join(GC_object.project_folder, GC_object.region_name, 'Elevation','Metadata'))



    ##########################################################################################################
    # these directories are constructed in the data_folder, where the downloaded data will be stored

    if GC_object.compile_elevation:
        # make the elevation folder
        if 'Elevation' not in os.listdir(os.path.join(GC_object.data_folder)):
            os.mkdir(os.path.join(GC_object.data_folder, 'Elevation'))

        # create subfolders for each source within the elevation folder
        if GC_object.compile_arcticDEM_data:
            if 'ArcticDEM' not in os.listdir(os.path.join(GC_object.data_folder, 'Elevation')):
                os.mkdir(os.path.join(GC_object.data_folder, 'Elevation','ArcticDEM'))

            if '2m_tiles' not in os.listdir(os.path.join(GC_object.data_folder, 'Elevation','ArcticDEM')):
                os.mkdir(os.path.join(GC_object.data_folder, 'Elevation','ArcticDEM','2m_tiles'))

            if 'Metadata' not in os.listdir(os.path.join(GC_object.data_folder, 'Elevation','ArcticDEM')):
                os.mkdir(os.path.join(GC_object.data_folder, 'Elevation','ArcticDEM','Metadata'))

            if 'Regridded_'+str(int(GC_object.elevation_grid_posting))+'m_tiles' not in os.listdir(os.path.join(GC_object.data_folder, 'Elevation','ArcticDEM')):
                os.mkdir(os.path.join(GC_object.data_folder, 'Elevation','ArcticDEM','Regridded_'+str(int(GC_object.elevation_grid_posting))+'m_tiles'))

    if GC_object.compile_icesat2_data:
        if 'ICESat2' not in os.listdir(os.path.join(GC_object.data_folder, 'Elevation')):
            os.mkdir(os.path.join(GC_object.data_folder, 'Elevation', 'ICESat2'))

        if 'Data' not in os.listdir(os.path.join(GC_object.data_folder, 'Elevation', 'ICESat2')):
            os.mkdir(os.path.join(GC_object.data_folder, 'Elevation', 'ICESat2','Data'))





