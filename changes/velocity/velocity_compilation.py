


def download_and_regrid_velocity_data(GD_object):

    message = 'Creating velocity compilation for '+GD_object.glacier_name
    GD_object.output_summary += '\n'+message
    if GD_object.print_main_outputs:
        print(message)

    if GD_object.run_preliminary_velocity_tests:
        a=1
        #this is where the scripts will be run first with the testing option turned on


    if GD_object.compile_golive_data:
        import changes.velocity.velocity_sources.compile_GOLIVE_data as GOLIVE
        GOLIVE.generate_GOLIVE_dataset(GD_object)

    if GD_object.compile_tsx_data:
        import changes.velocity.velocity_sources.compile_TSX_data as TSX
        TSX.generate_TSX_dataset(GD_object)