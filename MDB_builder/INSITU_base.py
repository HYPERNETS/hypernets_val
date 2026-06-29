import numpy as np
from datetime import datetime as dt
from netCDF4 import Dataset
import pytz,os

class INSITUBASE:

    def __init__(self, mdb_options):
        self.mdb_options = mdb_options
        self.new_MDB = None


    def start_add_insitu_no_rrs(self ,extract_path ,ofile, variables_to_exclude):
        if variables_to_exclude is not None and len(variables_to_exclude)>0:
            self.new_MDB = self.copy_nc_excluding_variables(extract_path,ofile,variables_to_exclude)
        else:
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
        instrument_id_list = self.mdb_options.insitu_options['instrument_ids']
        if not instrument_id_list[0]=='N/A':
            instrument_id_list.insert(0,'N/A')
        n_instrument_id = len(instrument_id_list)
        self.new_MDB.createDimension('insitu_id', n_insitu_id)
        #self.new_MDB.createDimension('insitu_original_bands',  n_insitu_bands)
        self.new_MDB.createDimension('insitu_bands', n_insitu_bands)
        self.new_MDB.createDimension('instrument_id',n_instrument_id)


        ##TIME VARIABLE
        insitu_time = self.new_MDB.createVariable('insitu_time', 'f8', ('satellite_id', 'insitu_id'), zlib=True,
                                                  complevel=6,fill_value=-999.0)
        insitu_time.units = "Seconds since 1970-01-01 00:00:00"
        insitu_time.description = 'In situ time in ISO 8601 format (UTC).'

        #INSTRUMENT_ID VARIABLE
        instrument_id_var = self.new_MDB.createVariable('insitu_instrument_id','i2',('satellite_id', 'insitu_id'),fill_value=-999,zlib=True,complevel=6)
        instrument_id_var.description = 'Instrument id'
        instrument_id_var.flag_values = np.arange(n_instrument_id).tolist()
        instrument_id_var.flag_meanings = " ".join(instrument_id_list)



        # ORIGINAL BANDS VARIABLE
        insitu_original_bands = self.new_MDB.createVariable('insitu_original_bands', 'f4', ('instrument_id','insitu_bands'),
                                                            fill_value=-999, zlib=True, complevel=6)
        insitu_original_bands.description = 'In situ bands in nm.'


        # RRS VARIABLE
        insitu_Rrs = self.new_MDB.createVariable('insitu_Rrs', 'f4', ('satellite_id', 'insitu_bands', 'insitu_id'),
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

    def copy_nc_excluding_variables(self, ifile, ofile,variablesToExclude):
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
                if name in variablesToExclude:
                    continue
                dst.createVariable(name, variable.datatype, variable.dimensions)
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)

                dst[name][:] = src[name][:]
        return dst

    def close(self):
        self.new_MDB.close()
        self.new_MDB = None


class Mini_MDB_Builder():

    def __init__(self,mdb_options,verbose):
        self.verbose = verbose
        self.mdb_options = mdb_options
        self.new_MDB = None

    def start_mini_mdb(self, extract_path, ofile):
        if self.verbose:
            print(f'[INFO]->Starting mini MDB...')

        self.new_MDB = copy_nc(extract_path, ofile)

        # time_window = 2  # del mdb_options
        # self.new_MDB.time_diff = f'{time_window * 60 * 60}'  # in seconds
        self.new_MDB.time_diff = self.mdb_options['time_diff_match_up']
        if self.mdb_options['time_diff_temporal_variability'] is not None:
            self.new_MDB.time_diff_temporal_variability = self.mdb_options['time_diff_temporal_variability']

        n_insitu_id = self.mdb_options['ninsitu_max']
        n_insitu_bands = self.mdb_options['n_insitu_bands']
        instrument_id_list = self.mdb_options['instrument_ids']
        n_instrument_id = len(instrument_id_list)
        self.new_MDB.createDimension('insitu_id', n_insitu_id)
        if n_insitu_bands>0:
            self.new_MDB.createDimension('insitu_bands', n_insitu_bands)
        self.new_MDB.createDimension('instrument_id',n_instrument_id+1)##a first "row" is reserved for 'default'


        ##TIME VARIABLE
        insitu_time = self.new_MDB.createVariable('insitu_time', 'f8', ('satellite_id', 'insitu_id'), zlib=True,
                                                  complevel=6,fill_value=-999.0)
        insitu_time.units = "Seconds since 1970-01-01 00:00:00"
        insitu_time.description = 'In situ time in ISO 8601 format (UTC).'

        #INSTRUMENT_ID VARIABLE
        instrument_id_var = self.new_MDB.createVariable('insitu_instrument_id','i2',('satellite_id', 'insitu_id'),fill_value=-999,zlib=True,complevel=6)
        instrument_id_var.description = 'Instrument id'
        flag_values_array = np.arange(n_instrument_id)+1
        instrument_id_var.flag_values = flag_values_array.tolist()
        instrument_id_var.flag_meanings = " ".join(instrument_id_list)


        if n_insitu_bands>0:
            # ORIGINAL BANDS VARIABLE
            insitu_original_bands = self.new_MDB.createVariable('insitu_original_bands', 'f4', ('instrument_id','insitu_bands'),fill_value=-999, zlib=True, complevel=6)
            insitu_original_bands.description = 'In situ bands in nm.'

            # RRS VARIABLE
            insitu_Rrs = self.new_MDB.createVariable('insitu_Rrs', 'f4', ('satellite_id', 'insitu_bands', 'insitu_id'),fill_value=-999, zlib=True, complevel=6)
            insitu_Rrs.description = 'In situ Rrs'

        # TIME DIFFERENCE VARIABLE
        time_difference = self.new_MDB.createVariable('time_difference', 'f4', ('satellite_id', 'insitu_id'),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
        time_difference.units = "seconds"

    def add_rrs_uncentainty_variable(self):
        insitu_Rrs = self.new_MDB.createVariable('insitu_Rrs_unc', 'f4', ('satellite_id', 'insitu_bands', 'insitu_id'),
                                                 fill_value=-999, zlib=True, complevel=6)
        insitu_Rrs.description = 'In situ Rrs uncentainties'

    def add_non_spectral_variable(self,var_name,attrs):
        var = self.new_MDB.createVariable(var_name,'f4',('satellite_id', 'insitu_id'),fill_value=-999,zlib=True,complevel=6)
        if attrs is not None:
            for at in attrs:
                value = attrs[at]
                if value is not None:
                    var.setncattr(at,value)

    def add_spectral_variable(self,var_name,attrs):
        var = self.new_MDB.createVariable(var_name, 'f4', ('satellite_id', 'insitu_bands','insitu_id'), fill_value=-999, zlib=True,complevel=6)
        if attrs is not None:
            for at in attrs:
                value = attrs[at]
                if value is not None:
                    var.setncattr(at, value)



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

    def set_insitu_wavelengths(self,instrument_id,wl_array):
        if not 'insitu_original_bands' in self.new_MDB.variables:
            print(f'[ERROR] insitu_original_band is not set in the input MDB mini file')
            return False
        var_wl = self.new_MDB.variables['insitu_original_bands']
        if instrument_id>var_wl.shape[0]:
            print(f'[ERROR] Instrument id {instrument_id}>{var_wl.shape[0]}. Wavelength array could not be set')
            return False
        if len(wl_array)!=var_wl.shape[1]:
            print(f'[ERROR] Number of wavelengths retrieved from in situ file different from the expected ({len(wl_array)} != {var_wl.shape[1]})')
            return False
        var_wl[instrument_id,:] = wl_array[:]
        return True

    def set_insitu_wavelengths_all_instruments(self, wl_array):
        if not 'insitu_original_bands' in self.new_MDB.variables:
            print(f'[ERROR] insitu_original_bands is not set in the input MDB mini file')
            return False
        var_wl = self.new_MDB.variables['insitu_original_bands']

        if var_wl.shape[0]==wl_array.shape[0]+1:##we add and extra row for default

            wl_array = np.concat([np.ma.masked_all((1,wl_array.shape[1])),np.ma.array(wl_array)],axis=0)
            wl_array[0,:] = np.ma.masked ##make sure the row added is made by masked values


        if wl_array.shape!=var_wl.shape:
            print(f'[ERROR] shape of the insitu_original_bands variable ({var_wl.shape}) does not match the shape of the wavelengths array {wl_array.shape}')
            return False
        else:
            var_wl[:] = wl_array[:]
            return True

    def set_instrument_id(self,ninsitu_real,instrument_id):
        if not 'insitu_instrument_id' in self.new_MDB.variables:
            print(f'[ERROR] insitu_instrument_id is not set in the input  MDB mini file')
            return False
        var = self.new_MDB.variables['insitu_instrument_id']
        if ninsitu_real>var.shape[1]:
            print(f'[ERROR] {ninsitu_real} is greater than the ninsitu_max set to {var.shape[1]}')
            return False
        var[0,0:ninsitu_real] = instrument_id
        return True

    #Non-spectral variables include insitu_time,insitu_lat,insitu_lon,insitu_spatial_index,time_difference or other data variables
    def set_non_spectral_variables(self,name_var,array):
        if not name_var in self.new_MDB.variables:
            print(f'[ERROR] {name_var} is not set in the input  MDB mini file')
            return False
        ninsitu_real = array.shape[0]
        var = self.new_MDB.variables[name_var]
        if ninsitu_real>var.shape[1]:
            print(f'[ERROR] {ninsitu_real} is greater than the ninsitu_max set to {var.shape[1]}')
            return False
        var[0,0:ninsitu_real] = array[:]
        return True

    def set_insitu_rrs(self,array):
        return self.set_spectral_variables('insitu_Rrs',array)

    ##spectral variables include insitu_Rrs, insitu_Rrs_unc or other
    def set_spectral_variables(self,name_var,array):
        if not name_var in self.new_MDB.variables:
            print(f'[ERROR] {name_var} is not set in the input  MDB mini file')
            return False
        nwl = array.shape[0]
        ninsitu_real = array.shape[1]
        var = self.new_MDB.variables[name_var]
        if nwl != var.shape[1]:
            print(f'[ERROR] Number of wavelengths in the data array {nwl} is different from the number of bands in the {name_var} variable({var.shape[1]})')
            return False
        if ninsitu_real>var.shape[2]:
            print(f'[ERROR] {ninsitu_real} is greater than the n_insitu in the variable {name_var} ({var.shape[1]})')
            return False
        var[0,:,0:ninsitu_real] = array[:]
        return True

    def set_insitu_basic_variables_from_dict(self,arrays,is_fixed_site=False):
        if 'insitu_time' in arrays:
            if not self.set_non_spectral_variables('insitu_time',arrays['insitu_time']):
                print(f'[ERROR] Error setting variable insitu_time')
                return False
        if 'insitu_lat' in arrays and not is_fixed_site:
            if not self.set_non_spectral_variables('insitu_latitude',arrays['insitu_lat']):
                print(f'[ERROR] Error setting variable insitu_lat')
                return False
        if 'insitu_lon' in arrays and not is_fixed_site:
            if not self.set_non_spectral_variables('insitu_longitude',arrays['insitu_lon']):
                print(f'[ERROR] Error setting variable insitu_lon')
                return False
        if 'insitu_spatial_index' in arrays and not is_fixed_site:
            if not self.set_non_spectral_variables('insitu_spatial_index', arrays['insitu_spatial_index']):
                print(f'[ERROR] Error setting variable insitu_spatial_index')
                return False
        if 'time_diff' in arrays:
            if not self.set_non_spectral_variables('time_difference', arrays['time_diff']):
                print(f'[ERROR] Error setting variable time_diff')
                return False
        return True
    # def set_data_from_array(self, array, variable_out):
    #     if not variable_out in self.new_MDB.variables:
    #         print(f'[WARNING] {variable_out} is not set in the input file')
    #         return False
    #
    #     try:
    #         self.new_MDB.variables[variable_out][:] = array[:]
    #         return True
    #     except:
    #         print(f'[WARNING] Array shape {array.shape} does not cast with variable {variable_out}')
    #         return False

    def close_mini_mdb_file(self):
        self.new_MDB.close()
        self.new_MDB = None


def add_line_csv_with_mdbm_info(fw,file_nc,started):
    if not started:
        first_line = 'name;satellite_id;insitu_id;instrument_id;satellite_bands;insitu_bands;rows;columns'
        fw.write(first_line)
    dset  = Dataset(file_nc)
    nsat = len(dset.dimensions['satellite_id'])
    ninsitu = len(dset.dimensions['insitu_id'])
    ninstrument = len(dset.dimensions['instrument_id'])
    nwlsat = len(dset.dimensions['satellite_bands'])
    nwlinsitu = len(dset.dimensions['insitu_bands'])
    rows = len(dset.dimensions['rows'])
    cols = len(dset.dimensions['columns'])
    line = f'{os.path.basename(file_nc)};{nsat};{ninsitu};{ninstrument};{nwlsat};{nwlinsitu};{rows};{cols}'
    dset.close()
    fw.write('\n')
    fw.write(line)
    dims = np.array([nsat,ninsitu,ninstrument,nwlsat,nwlinsitu,rows,cols])

    return fw,dims

def get_mini_mdb_file_path(extract_path,info):
    time_min_diff = dt.fromtimestamp(info['time_min_diff']).astimezone(pytz.utc)
    ref_time = time_min_diff.strftime('%Y%m%dT%H%M%S')
    site = info['site'].replace(' ','_')
    name = f'MDBm_{info["satellite"]}{info["platform"]}_{info["sensor"]}_{site}_{ref_time}.nc'
    file_m = os.path.join(extract_path,name)
    return file_m

def copy_nc(ifile, ofile):

    with Dataset(ifile) as src:
        dst = Dataset(ofile, 'w', format='NETCDF4')
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            dst.createVariable(name, variable.datatype, variable.dimensions)
             # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)

            dst[name][:] = src[name][:]
    return dst

def get_insitu_object(insitu_type,insitu_options,verbose):
    insituBase = None
    if insitu_type=='MULTIPLE_CSV':
        try:
            from INSITU_multiplecsv import INSITU_MULTIPLE_CSV
        except:
            from MDB_builder.INSITU_multiplecsv import INSITU_MULTIPLE_CSV
        insituBase = INSITU_MULTIPLE_CSV(insitu_options,verbose)
    if insitu_type=='SINGLE_SEABASS':
        try:
            from INSITU_SeaBass import INSITU_SEABASS
        except:
            from MDB_builder.INSITU_SeaBass import INSITU_SEABASS
        insituBase = INSITU_SEABASS(insitu_options,verbose)
    if insitu_type=='HYPSTAR_L2':
        try:
            from INSITU_hypernets import HYPSTAR_L2
        except:
            from MDB_builder.INSITU_hypernets import HYPSTAR_L2
        insituBase = HYPSTAR_L2(insitu_options,verbose)
    if insitu_type=='SO_RAD':
        try:
            from INSITU_sorad import SO_RAD
        except:
            from MDB_builder.INSITU_sorad import SO_RAD
        insituBase = SO_RAD(insitu_options,verbose)
    if insitu_type=='SINGLE_CSV':
        try:
            from INSITU_singlecsv import INSITU_SINGLE_CSV
        except:
            from MDB_builder.INSITU_singlecsv import INSITU_SINGLE_CSV
        insituBase = INSITU_SINGLE_CSV(insitu_options,verbose)

    if insitu_type=='AERONET_OC':
        try:
            from INSITU_aeronet import INSITU_AERONET
        except:
            from MDB_builder.INSITU_aeronet import INSITU_AERONET
        insituBase = INSITU_AERONET(insitu_options,verbose)

    if insituBase is None:
        print(f'[ERROR] In situ class for in situ data type {insitu_type} is not available. Please check method get_insitu_object() in INSITU_base.py')
    return insituBase