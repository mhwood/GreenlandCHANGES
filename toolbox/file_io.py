
import numpy as np
import xarray as xr

def read_timeseries_csv(file_path):
    f = open(file_path)
    lines = f.read()
    f.close()
    lines = lines.split('\n')
    lines.pop(0)

    timeseries = []
    sources = []
    for line in lines:
        line = line.split(',')
        timeseries.append([float(line[0]), float(line[1])])
        sources.append(line[2])

    timeseries=np.array(timeseries)

    return(timeseries,sources)

def read_timeseries_nc(file_path,variable):

    ds = xr.open_dataset(file_path)

    time = np.array(ds['time'])
    value = np.array(ds[variable.lower()])
    timeseries = np.column_stack([time,value])

    attr_dict = {}
    for att in list(ds.attrs.keys()):
        attr_dict[att] = ds.attrs[att]

    return(timeseries,attr_dict)

def write_timeseries_nc(file_path,variable,timeseries,attr_dict):
    time = timeseries[:,0]
    values = timeseries[:,1]

    ds = xr.Dataset({variable.lower(): (['time'], values)}, coords={'time': time})

    ds['time'].attrs['units'] = 'decimal_years'

    if variable == 'Velocity':
        ds[variable.lower()].attrs['units'] = 'meters_per_year'
    if variable == 'Elevation':
        ds[variable.lower()].attrs['units'] = 'meters'

    for att in list(attr_dict.keys()):
        ds.attrs[att] = attr_dict[att]

    ds.to_netcdf(file_path)