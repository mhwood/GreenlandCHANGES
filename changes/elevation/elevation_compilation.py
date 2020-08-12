

def download_and_regrid_elevation_data(GD_object):

    message = 'Creating elevation compilation for '+GD_object.region_name
    GD_object.output_summary += '\n'+message
    if GD_object.print_main_outputs:
        print(message)

    if GD_object.compile_arcticDEM_data:
        import changes.elevation.elevation_sources.compile_ArcticDEM_data as ArcticDEM
        ArcticDEM.generate_ArcticDEM_dataset(GD_object)

    if GD_object.compile_glistin_data:
        import changes.elevation.elevation_sources.compile_GLISTIN_data as GLISTIN
        GLISTIN.generate_glistin_dataset(GD_object)

    if GD_object.compile_gimp_data:
        import changes.elevation.elevation_sources.compile_GIMP_data as GIMP
        GIMP.generate_gimp_dataset(GD_object)

    if GD_object.compile_icesat2_data:
        import changes.elevation.elevation_sources.compile_ICESat2_data as ICESat2
        ICESat2.generate_ICESat2_dataset(GD_object)

