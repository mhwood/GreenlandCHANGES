
import numpy as np
import os
from toolbox.file_io import read_timeseries_nc
from toolbox.series import smooth_timeseries_seasonally
from scipy.interpolate import interp1d
from toolbox.file_io import write_timeseries_nc

def filter_by_seasonally_smoothing(timeseries,sources,n_sigma=2):
    differences = []
    smooth_timeseries = smooth_timeseries_seasonally(timeseries)

    set_int = interp1d(smooth_timeseries[:, 0], smooth_timeseries[:, 1])
    for vv in range(np.shape(timeseries)[0]):
        if timeseries[vv, 0] >= np.min(smooth_timeseries[:, 0]) and timeseries[vv, 0] <= np.max(smooth_timeseries[:, 0]):
            diff = np.abs(timeseries[vv, 1] - set_int(timeseries[vv, 0]))
            differences.append(diff)
    stdev = np.std(differences)

    filtered_timeseries = []
    filtered_sources = []
    for vv in range(np.shape(timeseries)[0]):
        if timeseries[vv, 0] >= np.min(smooth_timeseries[:, 0]) and timeseries[vv, 0] <= np.max(smooth_timeseries[:, 0]):
            diff = np.abs(timeseries[vv, 1] - set_int(timeseries[vv, 0]))

            if diff < n_sigma * stdev:
                new_line = []
                for col in range(np.shape(timeseries)[1]):
                    new_line.append(timeseries[vv,col])
                filtered_timeseries.append(new_line)
                filtered_sources.append(sources[vv])
    filtered_timeseries = np.array(filtered_timeseries)

    return(filtered_timeseries,filtered_sources)





def filter_timeseries(project_folder,glacier_name,variable,steps,output_type='nc'):
    file_path = os.path.join(project_folder, glacier_name, variable, 'Timeseries',
                 glacier_name + ' Median ' + variable + ' Timeseries.nc')

    if output_type=='nc':
        timeseries,attr_dict = read_timeseries_nc(file_path,variable)
        sources = attr_dict['sources'].split(',')


    for step in steps:
        if step == 'from_seasonal_smoothing':
            timeseries,sources = filter_by_seasonally_smoothing(timeseries,sources)
        else:
            TypeError('filtering step not recognized')

    file_path = os.path.join(project_folder, glacier_name, variable, 'Timeseries',
                             glacier_name + ' Median ' + variable + ' Timeseries - Filtered.nc')

    attr_dict['sources'] = sources
    attr_dict['filter_steps'] = ','.join(steps)

    write_timeseries_nc(file_path,variable,timeseries,attr_dict)



# This is for testing
# glacier_name = 'Jakobshavn'
# project_folder = '/Users/mhwood/Documents/Research/Projects/CHANGES/Examples/'
# variable = 'Velocity'
# steps = ['from_seasonal_smoothing','from_seasonal_smoothing']
# filter_timeseries(project_folder,glacier_name,variable,steps,output_type='nc')