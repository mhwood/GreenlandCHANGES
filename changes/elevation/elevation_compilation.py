

def download_and_regrid_elevation_data(GD_object):

    message = 'Creating elevation compilation for '+GD_object.region_name
    GD_object.output_summary += '\n'+message
    if GD_object.print_main_outputs:
        print(message)

    if GD_object.compile_arcticDEM_data:
        import changes.elevation.elevation_sources.compile_ArcticDEM_data as ArcticDEM
        ArcticDEM.generate_ArcticDEM_dataset(GD_object)