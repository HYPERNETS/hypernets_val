import argparse
import os.path
from datetime import timedelta
from datetime import datetime as dt

import numpy as np
from netCDF4 import Dataset

parser = argparse.ArgumentParser(
    description="Obtaining information for running MDB_builder.")

parser.add_argument('-m', "--mode", help='Mode option', choices=["add_instrument_id"], required=True)
parser.add_argument('-i', "--input_path", help="Input path.")
parser.add_argument('-o', "--output", help="Output file.")
args = parser.parse_args()

def main():
    print(f'[INFO] Started MDB_utils!')
    if args.mode=='add_instrument_id':
        if not check_required_params(['input_path']):
            return
        if not os.path.isfile(args.input_path):
            print(f'[ERROR] {args.input_path} does not exist or is not a valid directory')
            return
        output_path = args.output if args.output else None
        add_instrument_id(args.input_path,output_path)

def add_instrument_id(input_path,output_path):
    rename_file = False
    if output_path is None:
        rename_file = True
        output_path =os.path.join(os.path.dirname(input_path),'Temp.nc')
    ##WRITTING OUPPUT DATASET
    dataset = Dataset(input_path)
    ncout = Dataset(output_path, 'w', format='NETCDF4')
    # copy global attributes all at once via dictionary
    ncout.setncatts(dataset.__dict__)
    # copy dimensions
    for name, dimension in dataset.dimensions.items():
        name_new = 'insitu_bands' if name=='insitu_original_bands' else name
        ncout.createDimension(name_new, (len(dimension) if not dimension.isunlimited() else None))
    ncout.createDimension('instrument_id',1)

    # copy variables
    for name, variable in dataset.variables.items():
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue
        dims_l = list(variable.dimensions)
        if 'insitu_original_bands' in dims_l:
            dims_l[dims_l.index('insitu_original_bands')] = 'insitu_bands'
        if name=='insitu_original_bands':
            dims_l = ['instrument_id','insitu_bands']
        print(f'[INFO] Variable: {name} Dimensions: {dims_l}')
        ncout.createVariable(name, variable.datatype, tuple(dims_l), fill_value=fill_value, zlib=True,complevel=6)
        ncout[name].setncatts(dataset[name].__dict__)
        ncout[name][:] = dataset[name][:]

    #new variable
    var_instrument_id = ncout.createVariable('institu_instrument_id','i4',('satellite_id','insitu_id'),fill_value=-999,zlib=True,complevel=6)
    var_instrument_id.description = 'Instrument id'
    var_instrument_id.flag_values = [0]
    var_instrument_id.flag_meanings = ['N/A']
    instrument_id = dataset.variables['insitu_time'][:]
    instrument_id[instrument_id.mask==False]=0
    var_instrument_id[:] = instrument_id[:]



    ncout.close()
    dataset.close()
    if rename_file:
        os.rename(output_path,input_path)

def check_required_params(param_list):
    b = True
    for param in param_list:
        if not args.__dict__[param]:
            print(f'[ERROR] {param} is required for mode {args.mode}')
            b = False
    return b

# %%
if __name__ == '__main__':
    main()
