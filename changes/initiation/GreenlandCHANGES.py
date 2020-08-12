
# this script defines a class used to hold all pertinent information to compile the elevation and velocity
import time
import datetime
import os


class GreenlandCHANGES:

    def __init__(self, project_folder, data_folder):

        #this initiates the global domain and all of its associated parameters
        self.set_project_folder(project_folder)
        self.set_data_folder(data_folder)
        self.set_default_parameters()

    # ####################################################################################################################
    # # these are functions to set the metadata information

    def set_project_folder(self, project_folder):
        self.project_folder = project_folder

    def set_data_folder(self, data_folder):
        self.data_folder = data_folder

    # ####################################################################################################################
    # # these are functions to set parameters to prepare the retrieval and regridding scripts

    def create_directory_structure(self):
        import changes.initiation.directory_construction as dc
        dc.construct_output_directory_structure(self)

    def initiate_process_metadata_output(self):
        output_summary_header = 'Processes initiated on '+str(datetime.date.today())+'\n'
        self.output_summary += output_summary_header
        self.start_time = time.time()

    def set_default_parameters(self):

        self.print_main_outputs = True
        self.print_sub_outputs = True
        self.save_outputs_summary = True
        self.delete_downloaded_files = False

        self.region_initiated = False
        self.region_name = 'untitled_region'
        self.extents = []

        # this is metadata information pertaining to the process output
        self.output_summary = ''
        self.initiate_process_metadata_output()

        # this is metadata pertaining to the velocity compilation
        self.compile_velocity = True
        self.velocity_grid_posting = 300
        self.velocity_grid_epsg = 3413
        self.overwrite_existing_velocity_data = False
        self.create_velocity_stacks = True

        # GOLIVE
        self.compile_golive_data = True
        self.download_new_golive_data = True

        # TSX
        self.compile_tsx_data = True
        self.download_new_tsx_data = True

        # this is metadata pertaining to the elevation compilation
        self.compile_elevation = True
        self.elevation_grid_posting = 50
        self.elevation_grid_epsg = 3413
        self.overwrite_existing_elevation_data = False
        self.create_elevation_stacks = True

        # Arctic DEM
        self.compile_arcticDEM_data = True
        self.download_new_arcticDEM_data = True
        self.keep_high_resolution_arcticDEM_data = False
        self.resample_high_resolution_arcticDEM_data = True
        self.max_number_of_arcticDEM_files = 'all'

        # GIMP
        self.compile_gimp_data = True

        # GLISTIN
        self.compile_glistin_data = True
        self.download_new_glistin_data = True
        self.max_number_of_glistin_files = 'all'
        self.podaac_username = ''
        self.podaac_password = ''

        # ICESat-2
        self.compile_icesat2_data = True
        self.download_new_icesat2_data = True
        self.max_number_of_icesat2_files = 'all'

        self.compile_oib_data = True

    def set_custom_extents(self,extents):
        self.extents = extents

    def set_extents_by_glacier(self,glacier_name):
        self.region_name = glacier_name
        import changes.reference.glacier_domains as gd
        self.extents = gd.glacier_to_domain_extents[glacier_name]

    def deactivate_all_sources(self):
        self.compile_velocity = False
        self.compile_golive_data = False
        self.compile_tsx_data = False

        self.compile_elevation = False
        self.compile_arcticDEM_data = False
        self.compile_gimp_data = False
        self.compile_glistin_data = False
        self.compile_icesat2_data = False
        self.compile_oib_data = False

    # ####################################################################################################################
    # # these are functions to print out parameters for the user to view

    def print_initiation_parameters(self):
        self.check_region_intiated()
        print('Region Parameters:')
        print('    region_initiated: ', self.region_initiated)
        print('    region_name: ', self.region_name)
        print('    extents: ', self.extents)
        print(' ')
        print('Velocity Parameters:')
        print('    compile_velocity: ',self.compile_velocity)
        if self.compile_velocity:
            print('    velocity_grid_posting: ',self.velocity_grid_posting)
            print('    velocity_grid_epsg: ', self.velocity_grid_epsg)
            print('    create_velocity_stacks: ', self.create_velocity_stacks)
            print('    Velocity Sources:')
            print('        compile_golive_data:',self.compile_golive_data)
            print('        compile_tsx_data:', self.compile_tsx_data)
        print(' ')
        print('Elevation Parameters:')
        print('    compile_elevation: ', self.compile_elevation)
        if self.compile_elevation:
            print('    elevation_grid_posting: ', self.elevation_grid_posting)
            print('    elevation_grid_epsg: ', self.elevation_grid_epsg)
            print('    create_elevation_stacks: ', self.create_elevation_stacks)
            print('    Elevation Sources:')
            print('        compile_arcticDEM_data:', self.compile_arcticDEM_data)
            print('        compile_gimp_data:', self.compile_gimp_data)
            print('        compile_glistin_data:', self.compile_glistin_data)
            print('        compile_icesat2_data:', self.compile_icesat2_data)
            print('        compile_oib_data:', self.compile_oib_data)

    def print_arcticDEM_parameters(self):
        print('ArcticDEM Parameters:')
        print('    compile_arcticDEM_data: ', self.compile_arcticDEM_data)
        if self.compile_arcticDEM_data:
            print('    download_new_arcticDEM_data: ', self.download_new_arcticDEM_data)
            print('    keep_high_resolution_arcticDEM_data: ', self.keep_high_resolution_arcticDEM_data)
            print('    resample_high_resolution_arcticDEM_data: ', self.resample_high_resolution_arcticDEM_data)
            print('    max_number_of_arcticDEM_files: ', self.max_number_of_arcticDEM_files)

    def print_icesat2_parameters(self):
        print('ICESat-2 Parameters:')
        print('    compile_icesat2_data: ', self.compile_icesat2_data)
        if self.compile_icesat2_data:
            print('    download_new_icesat2_data: ', self.download_new_icesat2_data)
            print('    max_number_of_icesat2_files: ', self.max_number_of_icesat2_files)

    def print_glistin_parameters(self):
        print('GLISTIN-A Parameters:')
        print('    compile_glistin_data: ', self.compile_glistin_data)
        if self.compile_glistin_data:
            print('    download_new_glistin_data: ', self.download_new_glistin_data)
            print('    max_number_of_glistin_files: ', self.max_number_of_glistin_files)
            print('    podaac_username: ',self.podaac_username)
            print('    podaac_password: ', self.podaac_password)

    def print_gimp_parameters(self):
        print('GIMP Parameters:')
        print('    compile_gimp_data: ', self.compile_gimp_data)

    # ####################################################################################################################
    # # these are functions to set parameters to run the retrieval and regridding scripts

    def check_region_intiated(self):
        if len(self.extents)==4:
            extents_initiated = True
        else:
            extents_initiated = False
        if self.region_name == 'untitled_region':
            region_name_initiated = False
        else:
            region_name_initiated = True

        if extents_initiated and region_name_initiated:
            self.region_initiated = True


    def initiate_grids(self):
        import changes.initiation.grid_construction as gc
        gc.create_grids_from_extents(self)

    def write_process_metadata_output(self):
        self.end_time = time.time()
        total_time = self.end_time - self.start_time
        self.output_summary+='\n\nTotal time elapsed for compilation: '+str(datetime.timedelta(seconds=total_time))

        f=open(os.path.join(self.project_folder,self.region_name,self.region_name+' CHANGES Process Metadata.txt'),'w')
        f.write(self.output_summary)
        f.close()

    def execute_velocity_and_elevation_compilations(self):
        self.check_region_intiated()
        if not self.region_initiated:
            print('!! Alert !!')
            print('!! Please define a region name and valid extents before running the compilation !!')
            print('   There are two options:')
            print('   Option 1: Pre-defined glacier extent:')
            print('      GC.set_extents_by_glacier([glacier name here])')
            print('   Option 2: Custom region and extents:')
            print('      GC.region_name = [your region here]')
            print('      GC.extents = [min_x,min_y,max_x,max_y]')
        else:
            self.create_directory_structure()
            self.initiate_grids()

            if self.compile_velocity:
                import changes.velocity.velocity_compilation as vc
                vc.download_and_regrid_velocity_data(self)

            if self.compile_elevation:
                import changes.elevation.elevation_compilation as ec
                ec.download_and_regrid_elevation_data(self)

            self.write_process_metadata_output()

