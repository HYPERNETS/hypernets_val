import os
import shutil

from MDBFile import MDBFile
from netCDF4 import Dataset

class MDBWritter:

    def __init__(self,input_dataset,output_file):
        self.input_dataset = None
        if isinstance(input_dataset,MDBFile):
            self.input_dataset = input_dataset.nc
        elif isinstance(input_dataset,str):
            if os.path.isfile(input_dataset):
                try:
                    self.input_dataset = Dataset(input_dataset,'r')
                except:
                    print(f'[ERROR] {input_dataset} is not a valid NetCDF fiel')

        try:
            if os.path.exists(output_file):
                self.output_dataset = Dataset(output_file,'a',format='NETCDF4')
            else:
                self.output_dataset = Dataset(output_file, 'w', format='NETCDF4')
        except Exception as ex:
            print(f'[ERROR] Exception open/creating output dataset: {ex}')
            self.output_dataset = None



    def close(self):
        if self.input_dataset is not None:
            self.input_dataset.close()
        if self.output_dataset is not None:
            self.output_dataset.close()

    def copy_global_attributes(self):
        if self.output_dataset is None or self.input_dataset is None:
            return
        # copy global attributes all at once via dictionary
        self.output_dataset.setncatts(self.input_dataset.__dict__)

    def copy_dimensions(self):
        if self.output_dataset is None or self.input_dataset is None:
            return
        # copy dimensions
        for name, dimension in self.input_dataset.dimensions.items():
            self.output_dataset.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))


    def copy_variables(self,variables_keep,variables_remove):
        if self.output_dataset is None or self.input_dataset is None:
            return
        for name, variable in self.input_dataset.variables.items():
            if len(variables_keep) > 0:
                if not name in variables_keep:
                    continue
            if len(variables_remove) > 0:
                if name in variables_remove:
                    continue
            fill_value = None
            if '_FillValue' in list(variable.ncattrs()):
                fill_value = variable._FillValue

            # create variable
            self.output_dataset.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True, complevel=6)
            # copy variable attributes all at once via dictionary
            self.output_dataset[name].setncatts(self.input_dataset[name].__dict__)
            # copy data
            self.output_dataset[name][:] = self.input_dataset[name][:]

    def get_dims_from_shape(self,shape):
        if self.output_dataset is None:
            return
        dim_names = list(self.output_dataset.dimensions.keys())
        dim_sizes = [len(self.output_dataset.dimensions[d]) for d in self.output_dataset.dimensions]
        try:
            ldims = [dim_names[dim_sizes.index(s)] for s in shape]
            return tuple(ldims)
        except:
            return None

    def add_variable(self,var_name,array,dims,dtype,fill_value):
        if self.output_dataset is None:
            return
        if dtype is None:
            dtype = 'f4'
        if dims is None:
            dims = self.get_dims_from_shape(array.shape)
        if var_name in self.output_dataset.variables:
            print(f'[ERROR] Ouput variable {var_name} already exists')
            return
        try:
            # create variable
            self.output_dataset.createVariable(var_name, dtype, dims, fill_value=fill_value,zlib=True, complevel=6)
            self.output_dataset[var_name][:] = array[:]
        except:
            print(f'[ERROR] Variable {var_name} could not be created')
            return


def copy_nc(file_in,file_out):
    shutil.copy(file_in,file_out)










