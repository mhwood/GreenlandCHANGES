


def download_and_regrid_velocity_data(GD_object):

    message = 'Creating velocity compilation for '+GD_object.region_name
    GD_object.output_summary += '\n'+message
    if GD_object.print_main_outputs:
        print(message)


    if GD_object.compile_golive_data:
        from .velocity_sources import compile_GOLIVE_data as GOLIVE
        GOLIVE.generate_GOLIVE_dataset(GD_object)

    if GD_object.compile_measures_insar_data:
        from .velocity_sources import compile_MEaSUREs_InSAR_data as MI
        MI.generate_MEaSUREs_InSAR_dataset(GD_object)

    if GD_object.compile_measures_optical_data:
        from .velocity_sources import compile_MEaSUREs_Optical_data as MO
        MO.generate_MEaSUREs_Optical_dataset(GD_object)

    if GD_object.compile_measures_quarterly_mosaic_data:
        from .velocity_sources import compile_MEaSUREs_Quarterly_Mosaic_data as MQ
        MQ.generate_MEaSUREs_Quarterly_Mosaic_dataset(GD_object)

    if GD_object.compile_measures_multiyear_mosaic_data:
        from .velocity_sources import compile_MEaSUREs_Multiyear_Mosaic_data as MQ
        MQ.generate_MEaSUREs_Multiyear_Mosaic_dataset(GD_object)