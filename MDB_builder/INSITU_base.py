
class INSITUBASE:

    def __init__(self, mdb_options):
        self.mdb_options = mdb_options
        self.new_MDB = None


    def start_add_insitu_no_rrs(self ,extract_path ,ofile):
        self.new_MDB = self.copy_nc(extract_path, ofile)
        self.new_MDB.time_diff = self.mdb_options.insitu_options['time_window']

        if 'insitu_id' not in self.new_MDB.dimensions:
            n_insitu_id = self.mdb_options.insitu_options['n_insitu_id']
            self.new_MDB.createDimension('insitu_id', n_insitu_id)

        ##TIME VARIABLE
        if not 'insitu_time' in self.new_MDB.variables:
            insitu_time = self.new_MDB.createVariable('insitu_time', 'f8', ('satellite_id', 'insitu_id',), zlib=True,
                                                  complevel=6)
            insitu_time.units = "Seconds since 1970-1-1"
            insitu_time.description = 'In situ time in ISO 8601 format (UTC).'

        # TIME DIFFERENCE VARIABLE
        if not 'time_difference' in self.new_MDB.variables:
            time_difference = self.new_MDB.createVariable('time_difference', 'f4', ('satellite_id', 'insitu_id')
                                                          ,fill_value=-999, zlib=True, complevel=6)
            time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
            time_difference.units = "seconds"


    def add_insitu(self, extract_path, ofile):
        self.start_add_insitu(extract_path, ofile)

    def start_add_insitu(self, extract_path, ofile):

        self.new_MDB = self.copy_nc(extract_path, ofile)

        # time_window = 2  # del mdb_options
        # self.new_MDB.time_diff = f'{time_window * 60 * 60}'  # in seconds
        self.new_MDB.time_diff = self.mdb_options.insitu_options['time_window']

        n_insitu_id = self.mdb_options.insitu_options['n_insitu_id']
        n_insitu_bands = self.mdb_options.insitu_options['n_insitu_bands']
        self.new_MDB.createDimension('insitu_id', n_insitu_id)
        self.new_MDB.createDimension('insitu_original_bands',  n_insitu_bands)

        ##TIME VARIABLE
        insitu_time = self.new_MDB.createVariable('insitu_time', 'f8', ('satellite_id', 'insitu_id',), zlib=True,
                                                  complevel=6)
        insitu_time.units = "Seconds since 1970-1-1"
        insitu_time.description = 'In situ time in ISO 8601 format (UTC).'

        ##FILENAME VARIABLE->DEPRECATED
        # insitu_filename = self.new_MDB.createVariable('insitu_filename', 'S2', ('satellite_id', 'insitu_id'), zlib=True,
        #                                               complevel=6)
        # insitu_filename.description = 'In situ filename.'

        # ORIGINAL BANDS VARIABLE
        insitu_original_bands = self.new_MDB.createVariable('insitu_original_bands', 'f4', ('insitu_original_bands'),
                                                            fill_value=-999, zlib=True, complevel=6)
        insitu_original_bands.description = 'In situ bands in nm.'

        # RRS VARIABLE
        insitu_Rrs = self.new_MDB.createVariable('insitu_Rrs', 'f4', ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                            fill_value=-999, zlib=True, complevel=6)
        insitu_Rrs.description = 'In situ Rrs'

        # TIME DIFFERENCE VARIABLE
        time_difference = self.new_MDB.createVariable('time_difference', 'f4', ('satellite_id', 'insitu_id'),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
        time_difference.units = "seconds"

        # self.new_MDB.close()


    def add_shipborne_variables(self):

        if not 'insitu_spatial_index' in self.new_MDB.variables:
            index_spatial = self.new_MDB.createVariable('insitu_spatial_index', 'i2', ('satellite_id', 'insitu_id')
                                                        ,fill_value=-999 ,zlib=True, complevel=6)
            index_spatial.long_name = "Distance to the central pixel starting from zero"

        if not 'insitu_latitude' in self.new_MDB.variables:
            insitu_latitude = self.new_MDB.createVariable('insitu_latitude', 'f8', ('satellite_id', 'insitu_id',), zlib=True, complevel=6, fill_value=-999)
            insitu_latitude.long_name = "In situ latitude"

        if not 'insitu_longitude' in self.new_MDB.variables:
            insitu_longitude = self.new_MDB.createVariable('insitu_longitude', 'f8', ('satellite_id', 'insitu_id',)
                                                           ,zlib=True ,complevel=6, fill_value=-999)
            insitu_longitude.long_name = "In situ longitude"


    def add_general_variable(self ,name_var, dims, data_type ,ats):
        if not name_var.startswith('satellite_'):
            name_var = f'satellite_{name_var}'
        fill_value = -999
        if data_type == 'i1':
            fill_value = -1
        variable = self.new_MDB.createVariable(name_var, data_type, dims, zlib=True,
                                               complevel=6, fill_value=fill_value)
        if len(ats ) >0:
            for at in ats:
                if at=='_FillValue': continue
                variable.setncattr(at ,ats[at])

    def add_data_variable(self,file_in,variable_in,variable_out):
        if not variable_out in self.new_MDB.variables:
            print(f'[WARNING] {variable_out} is not set in the input file')
            return
        from netCDF4 import Dataset
        din = Dataset(file_in)
        if variable_in in din.variables:
            try:
                self.new_MDB.variables[variable_out][:] = din.variables[variable_in][:]
            except:
                print(f'[WARNING] Dimensions for var {variable_in} in {file_in} are not the same for variable out {variable_out}')
                return
        else:
            print(f'[WARNING] {variable_in} is not available in {file_in}')
            return
        din.close()


    def get_name_satellite_variable(self,name_var):
        return f'satellite_{name_var}' if not name_var.startswith('satellite_') else name_var


    def add_insitu_variable(self ,name_var ,data_type ,ats):
        if not name_var.startswith('insitu_'):
            name_var = f'insitu_{name_var}'
        fill_value = -999
        if data_type=='i1':
            fill_value = -1
        variable = self.new_MDB.createVariable(name_var ,data_type ,('satellite_id', 'insitu_id',) ,zlib=True
                                               ,complevel=6 ,fill_value=fill_value)
        if ats is None:
            return
        if len(ats )>0:
            for at in ats:
                variable.setncattr(at ,ats[at])

    def copy_nc_reduced(self ,ifile ,ofile):
        from netCDF4 import Dataset
        with Dataset(ifile) as src:
            dst = Dataset(ofile, 'w', format='NETCDF4')
            # copy global attributes all at once via dictionary
            dst.setncatts(src.__dict__)
            # copy dimensions, setting a longitude of 1 for insitu_id
            for name, dimension in src.dimensions.items():
                if name=='insitu_id':
                    dst.createDimension(name ,1)
                else:
                    dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

            # copy all file data except for the excluded
            for name, variable in src.variables.items():
                dst.createVariable(name, variable.datatype, variable.dimensions)
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)
                # copy data for all the variables which does not incluide insitu_id
                if 'insitu_id' not in variable.dimensions:
                    dst[name][:] = src[name][:]
        return dst


    def copy_nc(self, ifile, ofile):
        from netCDF4 import Dataset
        with Dataset(ifile) as src:
            dst = Dataset(ofile, 'w', format='NETCDF4')
            # copy global attributes all at once via dictionary
            dst.setncatts(src.__dict__)
            # copy dimensions
            for name, dimension in src.dimensions.items():
                dst.createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else None))
            # copy all file data except for the excluded
            for name, variable in src.variables.items():
                dst.createVariable(name, variable.datatype, variable.dimensions)
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)

                dst[name][:] = src[name][:]
        return dst

    def close(self):
        self.new_MDB.close()
        self.new_MDB = None
