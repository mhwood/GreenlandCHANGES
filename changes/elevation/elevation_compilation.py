

def download_and_regrid_elevation_data(GD_object):

    message = 'Creating elevation compilation for '+GD_object.region_name
    GD_object.output_summary += '\n'+message
    if GD_object.print_main_outputs:
        print(message)

    if GD_object.compile_arcticDEM_data:
        from .elevation_sources import compile_ArcticDEM_data as ArcticDEM
        ArcticDEM.generate_ArcticDEM_dataset(GD_object)

    if GD_object.compile_glistin_data:
        from .elevation_sources import compile_GLISTIN_data as GLISTIN
        GLISTIN.generate_glistin_dataset(GD_object)

    if GD_object.compile_gimp_data:
        from .elevation_sources import compile_GIMP_data as GIMP
        GIMP.generate_gimp_dataset(GD_object)

    if GD_object.compile_icebridge_atm_data:
        from .elevation_sources import compile_Icebridge_ATM_data as IATM
        IATM.generate_Icebridge_ATM_dataset(GD_object)

    if GD_object.compile_icesat2_data:
        from .elevation_sources import compile_ICESat2_data as ICESat2
        ICESat2.generate_ICESat2_dataset(GD_object)

    if GD_object.compile_kms_data:
        from .elevation_sources import compile_KMS_data as KMS
        KMS.generate_kms_dataset(GD_object)

    if GD_object.compile_tandemx_data:
        from .elevation_sources import compile_TanDEMX_data as TanDEMX
        TanDEMX.generate_TanDEMX_dataset(GD_object)

