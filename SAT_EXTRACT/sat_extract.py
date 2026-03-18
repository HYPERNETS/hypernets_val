
import numpy as np
from netCDF4 import Dataset
from datetime import datetime as dt
import numpy.ma as ma
import argparse,os,subprocess,pytz,__init__,sys,shutil
code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
import COMMON.common_functions as cfs
import COMMON.args_functions as arf
from OPTIONS.OptionsManager import OptionsManager
from datetime import timedelta
from multiprocessing import Pool
import MDB_builder.INSITU_base as ibase


class SatExtractOptions:

    def __init__(self,config_file,verbose):
        self.verbose = verbose
        self.omanager = None
        self.gmanager = None

        general_options_file = os.path.join(code_home,'OPTIONS','general_options.ini')
        self.gmanager = OptionsManager(general_options_file,None)
        self.omanager = OptionsManager(config_file, None)

        self.is_valid = self.gmanager.is_valid() and self.omanager.is_valid()

    def get_general_options(self,section):
        if not self.is_valid:
            return None
        soptions, required = self.gmanager.get_retrieve_options(section)

        if soptions is not None:
            options = self.omanager.get_options_as_dict(section, soptions,required)
            return options
        else:
            return None

    def check_insitu_options(self,insitu_type):
        insitu_options_file = os.path.join(code_home, 'OPTIONS', 'insitu_options.ini')
        om = OptionsManager(insitu_options_file, None)
        if not om.is_valid():
            return None
        insitu_type_list = om.get_option('GLOBAL','type_list',None,None,'strlist')
        if insitu_type_list is None:
            print(f'[ERROR] In situ type list (option GLOBAL/type_list) in file {insitu_options_file} is unavailable or not valid')
            return [None]*2

        if insitu_type is not None:
            if insitu_type not in insitu_type_list:
                print(f'[ERROR] In situ type {insitu_type} is not defined. Please choose among: {insitu_type_list}')
                return [None]*2
            else:
                sopt, required = self.get_retrive_options_insitu(insitu_type)
                if sopt is None:
                    print(f'[ERROR] Section {insitu_type} defining in situ type {insitu_type} is not available in {insitu_options_file}.')
                    return [None] * 2
                insitu_options = self.omanager.get_options_as_dict(insitu_type, sopt, required)
                if insitu_options is None:
                    print(f'[ERROR] In situ options for in situ type {insitu_type} could not be retrieved')
                    return [None] * 2
                return insitu_type,insitu_options


        insitu_options = None
        for opt in insitu_type_list:
            sopt,required = self.get_retrive_options_insitu(opt)
            if sopt is None:
                print(f'[ERROR] Section {opt} defining in situ type {opt} is not available in {insitu_options_file}.')
                return [None]*2
            insitu_options = self.omanager.get_options_as_dict(opt,sopt,required)
            if insitu_options is not None:
                insitu_type = opt
                break
        if insitu_type is None:
            print(f'[ERROR] In situ type is not defined in the configuration file {insitu_options_file}.')
            print(f'[ERROR] Please choose among: {insitu_type_list}')

        return insitu_type,insitu_options

    def check_sat_type(self,sat_type):
        sat_options_file = os.path.join(code_home, 'OPTIONS', 'satellite_options.ini')
        om = OptionsManager(sat_options_file, None)
        if not om.is_valid():
            return None
        sat_type_list = om.get_option('GLOBAL', 'type_list', None, None, 'strlist')
        if sat_type_list is None:
            print(
                f'[ERROR] Satellite type list (option GLOBAL/type_list) in file {sat_options_file} is unavailable or not valid')
            return None
        if sat_type is not None:
            if sat_type in sat_type_list:
                if self.verbose:
                    print(f'[INFO] Satellite data type: {sat_type}')
                return sat_type
            else:
                print(f'[ERROR] Satellite type {sat_type} is not defined. Please choose among: {sat_type_list}')
                return None
        for opt in sat_type_list:
            if self.omanager.get_options_list(opt) is not None:
                sat_type = opt
                break
        if sat_type is None:
            print(f'[ERROR] Satellite type is not defined in the configuration file {sat_options_file}.')
            print(f'[ERROR] Please choose among: {sat_type_list}')
        else:
            if self.verbose:
                print(f'[INFO] Satellite data type: {sat_type}')
        return sat_type

    def get_satellite_options(self,section):
        if not self.is_valid:
            return None
        soptions, required = self.get_retrive_options_sat(section)
        if soptions is not None:
            options = self.omanager.get_options_as_dict(section, soptions, required)
            return options
        else:
            print(f'[ERROR] Section for satellite type {section} is not available')
            return None

    def get_retrive_options_insitu(self,section):
        insitu_options_file = os.path.join(code_home, 'OPTIONS', 'insitu_options.ini')
        om = OptionsManager(insitu_options_file, None)
        if om.is_valid():
            return om.get_retrieve_options(section)
        else:
            return None

    def get_retrive_options_sat(self,section):
        insitu_options_file = os.path.join(code_home, 'OPTIONS', 'satellite_options.ini')
        om = OptionsManager(insitu_options_file, None)
        if om.is_valid():
            return om.get_retrieve_options(section)
        else:
            return None


    def prepare_config_file_for_sh_slurm(self,path_csv,file_copy):
        if path_csv is not None:
            self.omanager.add_value('MULTIPLE_CSV','path_csv',path_csv)
            self.omanager.add_value('MULTIPLE_CSV', 'col_date', 'date')
            self.omanager.add_value('MULTIPLE_CSV', 'col_lat', 'lat')
            self.omanager.add_value('MULTIPLE_CSV', 'col_lon', 'lon')
            self.omanager.add_value('MULTIPLE_CSV', 'col_sep', ';')
            self.omanager.add_value('MULTIPLE_CSV', 'format_date', '%%Y-%%m-%%d %%H:%%M:%%S.%%f')
            self.omanager.remove_option('MULTIPLE_CSV','col_time')
            self.omanager.remove_option('MULTIPLE_CSV', 'format_time')


        self.omanager.add_value('multiprocessing', 'use_slurm_sh', 'False')
        self.omanager.save_copy_as_file(file_copy)

    def get_output_path(self):
        options = self.get_general_options('file_path')
        return options['output_dir'] if options is not None else None

    def get_box_size(self):
        options = self.get_general_options('satellite_options')
        return options['extract_size'] if (options is not None and 'extract_size' in options) else None

    def get_overwrite(self):
        options = self.get_general_options('file_path')
        return options['overwrite'] if options is not None else None


    # def get_lat_long_var_names(self):
    #     options = self.get_general_options('satellite_options')
    #     if options is None:
    #         return [None]*2
    #     return options['lat_variable'],options['lon_variable']


    def get_input_path_info(self):
        options = self.get_general_options('file_path')
        if options is None:
            return None
        path_source = options['sat_source_dir']
        org = options['sat_source_dir_organization']
        wce = options['wce']
        if path_source is None:
            print(f'[ERROR] path_source')
            return None


        if org is not None:
            if org.lower() == 'none':
                org = None
            elif org == 'YYYYmmdd':
                org = '%Y/%m/%d'
            elif org == 'YYYYjjj':
                org = '%Y/%j'
            elif org == 'YYYYmm':
                org = '%Y/%m'
            else:
                org = org.replace('YYYY', '%Y')
                org = org.replace('mm', '%m')
                org = org.replace('dd', '%d')
                org = org.replace('jjj', '%j')

            if org is not None:
                try:
                    dt.now().strftime(org)
                except:
                    print(f'[ERROR] Satellite path date organization {org} is not valid')
                    return None

        input_path_info = {
            'path_source': path_source,
            'org': org,
            'wce': wce,
            'unzip_dir': options['unzip_dir']
        }
        return input_path_info


    def get_multiprocessing_options(self):
        options = self.get_general_options('multiprocessing')
        return options

    def get_satellite_global_atrib(self):
        attr_keys = ['satellite', 'platform', 'sensor', 'res', 'aco_processor', 'proc_version','cmems_region']
        options_cfile = self.get_general_options('satellite_options')
        at = {}
        for key in attr_keys:
            if key in options_cfile and options_cfile[key] is not None:
                at[key] = options_cfile[key]
            else:
                at[key] = ''

        return at



class SatExtract:
    def __init__(self, ofname):

        if ofname is None:
            self.EXTRACT = None
            self.FILE_CREATED = False
            pass

        self.FILE_CREATED = True
        try:
            self.EXTRACT = Dataset(ofname, 'w', format='NETCDF4')
        except:
            print('[ERROR] Error creating: ', ofname)
            self.FILE_CREATED = False

        self.geometry_variables = {
            'satellite_SZA': {
                'long_name': 'Sun Zenith Angle',
                'units': 'degrees'
            },
            'satellite_SAA': {
                'long_name': 'Sun Azimuth Angle',
                'units': 'degrees'
            },
            'satellite_OZA': {
                'long_name': 'Observation Zenith Angle',
                'units': 'degrees'
            },
            'satellite_OAA': {
                'long_name': 'Observation Azimuth Angle',
                'units': 'degrees'
            }
        }

    def set_variable_attributes(self,name_var,attrs):
        if not name_var in self.EXTRACT.variables:
            print(f'[WARNING] Variable {name_var} was not found. Attributed can not be set')
            return
        if '_FillValue' in attrs:
            del attrs['_FillValue']
        self.EXTRACT.variables[name_var].setncatts(attrs)

    def set_global_attributes(self, at):
        # Atributes
        self.EXTRACT.creation_time = dt.now().strftime("%Y-%m-%dT%H:%M:%SZ")
        satellite = at['satellite']
        platform = at['platform']
        sensor = at['sensor']
        res_str = at['res']
        cmems_region = at['cmems_region']
        self.EXTRACT.satellite = satellite
        self.EXTRACT.platform = platform
        self.EXTRACT.sensor = sensor
        if len(res_str)>0:
            self.EXTRACT.resolution = res_str
        if len(cmems_region)>0:
            self.EXTRACT.cmems_region = cmems_region
        self.EXTRACT.description = f'{satellite}{platform} {sensor.upper()} {res_str} L2 extract'
        self.EXTRACT.satellite_aco_processor = at['aco_processor']  # 'Atmospheric Correction processor: xxx'
        self.EXTRACT.satellite_proc_version = at['proc_version']  # proc_version_str


        self.EXTRACT.site = at['site']

    def create_dimensions_basic(self, size_box):
        self.EXTRACT.createDimension('satellite_id', None)
        self.EXTRACT.createDimension('rows', size_box)
        self.EXTRACT.createDimension('columns', size_box)

    def create_dimensions_basic_includinginsitu(self, size_box, n_insitubands, n_insituid):
        self.EXTRACT.createDimension('satellite_id', None)
        self.EXTRACT.createDimension('rows', size_box)
        self.EXTRACT.createDimension('columns', size_box)

        self.EXTRACT.createDimension('insitu_original_bands', n_insitubands)
        self.EXTRACT.createDimension('insitu_id', n_insituid)

    def create_dimensions(self, size_box, n_bands):
        # dimensions
        self.EXTRACT.createDimension('satellite_id', None)
        self.EXTRACT.createDimension('rows', size_box)
        self.EXTRACT.createDimension('columns', size_box)
        if n_bands > 0:
            self.EXTRACT.createDimension('satellite_bands', n_bands)

    def create_dimensions_incluidinginsitu(self, size_box, n_bands, n_insitubands, n_insituid):
        # dimensions
        self.EXTRACT.createDimension('satellite_id', None)
        self.EXTRACT.createDimension('rows', size_box)
        self.EXTRACT.createDimension('columns', size_box)
        self.EXTRACT.createDimension('satellite_bands', n_bands)

        self.EXTRACT.createDimension('insitu_original_bands', n_insitubands)
        self.EXTRACT.createDimension('insitu_id', n_insituid)

    def create_dimensions_includingbasicinsitu(self, size_box, n_bands, n_insituid):
        self.EXTRACT.createDimension('satellite_id', None)
        self.EXTRACT.createDimension('rows', size_box)
        self.EXTRACT.createDimension('columns', size_box)
        self.EXTRACT.createDimension('satellite_bands', n_bands)

        self.EXTRACT.createDimension('insitu_id', n_insituid)

    def create_dimension_insitu(self, n_insituid):
        self.EXTRACT.createDimension('insitu_id', n_insituid)

    def create_geometry_variables(self):
        for name_geom in self.geometry_variables:
            self.create_geometry_variable(name_geom)

    def create_geometry_variable(self, name_geom):
        if not name_geom in self.geometry_variables.keys():
            return None
        gvar = self.EXTRACT.createVariable(name_geom, 'f4', ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                           zlib=True, complevel=6)
        gvar.units = self.geometry_variables[name_geom]['units']
        gvar.long_name = self.geometry_variables[name_geom]['long_name']
        return gvar

    def create_satellite_time_variable(self, satellite_start_time):
        satellite_time = self.EXTRACT.createVariable('satellite_time', 'f8', ('satellite_id'), fill_value=-999,
                                                     zlib=True, complevel=6)

        satellite_time[0] = float(satellite_start_time.replace(tzinfo=pytz.utc).timestamp())
        satellite_time.units = "Seconds since 1970-1-1"

    def create_pdu_variable(self, pdu, sensor):

        satellite_PDU = self.EXTRACT.createVariable('satellite_PDU', 'S1', ('satellite_id'), zlib=True,
                                                    complevel=6)  # string

        satellite_PDU[0] = pdu
        satellite_PDU.long_name = f'{sensor} source PDU name'

    def create_lat_long_variables(self, lat, lon, window):

        start_idx_y = window[0]
        stop_idx_y = window[1]
        start_idx_x = window[2]
        stop_idx_x = window[3]
        nrows = (window[1] - window[0])
        ncols = (window[3] - window[2])

        # latitude
        satellite_latitude = self.EXTRACT.createVariable('satellite_latitude', 'f8',
                                                         ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                                         zlib=True, complevel=6)

        if lat.ndim == 1:
            for c in range(ncols):
                satellite_latitude[0, :, c] = [lat[start_idx_y:stop_idx_y]]
        else:
            satellite_latitude[0, :, :] = [lat[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]]

        satellite_latitude.short_name = 'latitude'

        # longitude
        satellite_longitude = self.EXTRACT.createVariable('satellite_longitude', 'f8',
                                                          ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                                          zlib=True, complevel=6)
        if lon.ndim == 1:
            for r in range(nrows):
                satellite_longitude[0, r, :] = [lon[start_idx_x:stop_idx_x]]
        else:
            satellite_longitude[0, :, :] = [lon[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]]
        satellite_longitude.short_name = 'longitude'

    def create_satellite_bands_variable(self, wavelengths):
        # Variable satellite_bands (wavelenghts for Rrs)
        satellite_bands = self.EXTRACT.createVariable('satellite_bands', 'f4', ('satellite_bands'), fill_value=-999,
                                                      zlib=True, complevel=6)
        satellite_bands[:] = wavelengths
        satellite_bands.units = 'nm'

    def create_rrs_variable(self, sensor):
        # Variable satellite_Rrs (NOT BRDF-corrected remote sensing reflectance)
        satellite_Rrs = self.EXTRACT.createVariable('satellite_Rrs', 'f4',
                                                    ('satellite_id', 'satellite_bands', 'rows', 'columns'),
                                                    fill_value=-999, zlib=True, complevel=6)
        satellite_Rrs.short_name = 'Satellite Rrs'
        satellite_Rrs.long_name = f'Above water Remote Sensing Reflectance for {sensor} acquisition'
        satellite_Rrs.units = "sr-1"
        return satellite_Rrs

    def create_rrs_unc_variable(self, sensor):
        # Variable satellite_Rrs (NOT BRDF-corrected remote sensing reflectance uncertainty)
        satellite_Rrs = self.EXTRACT.createVariable('satellite_Rrs_unc', 'f4',
                                                    ('satellite_id', 'satellite_bands', 'rows', 'columns'),
                                                    fill_value=-999, zlib=True, complevel=6)
        satellite_Rrs.short_name = 'Satellite Rrs uncertainty'
        satellite_Rrs.long_name = f'Uncertainty in above water remote sensing reflectance for {sensor} acquisition'
        satellite_Rrs.units = "sr-1"
        return satellite_Rrs

    def create_aot_variable(self, aot, window):
        start_idx_y = window[0]
        stop_idx_y = window[1]
        start_idx_x = window[3]
        stop_idx_x = window[4]
        # AOT
        satellite_AOT_0865p50_box = self.EXTRACT.createVariable('satellite_AOT_0865p50', 'f4',
                                                                ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                                                zlib=True, complevel=6)
        satellite_AOT_0865p50_box[0, :, :] = ma.array(aot[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])
        satellite_AOT_0865p50_box.description = 'Satellite Aerosol optical thickness'

    def create_flag_variable(self, var_name, var_array, description, flag_masks, flag_meanings, window,var_dtype='f4'):
        start_idx_y = window[0]
        stop_idx_y = window[1]
        start_idx_x = window[2]
        stop_idx_x = window[3]
        # Quality Flags
        satellite_flag = self.EXTRACT.createVariable(var_name, var_dtype, ('satellite_id', 'rows', 'columns'),
                                                     fill_value=-999, zlib=True, complevel=6)
        satellite_flag[0, :, :] = [ma.array(var_array[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])]
        satellite_flag.description = description
        satellite_flag.flag_masks = flag_masks
        satellite_flag.flag_meanings = flag_meanings

    def create_2D_variable_from_info_dict(self,var_info,var_name):
        if var_name is not None:
            if var_name in var_info.keys():
                var =  self.create_2D_variable_from_info_dict_impl(var_info[var_name],var_name)
                return True if var is not None else False
            else:
                print(f'[ERROR] Info for variable {var_name} is not available')
                return False
        else:
            for var_name in var_info:
                var = self.create_2D_variable_from_info_dict_impl(var_info[var_name],var_name)
                if var is None:
                    return False
            return True



    def create_2D_variable_from_info_dict_impl(self,var_info,var_name):
        try:

            var = self.EXTRACT.createVariable(var_name,var_info['data_type'],('satellite_id', 'rows', 'columns'),fill_value=var_info['fill_value'],zlib=True, complevel=6)
            attrs = var_info['attrs']
            if attrs is not None:
                var.setncatts(attrs)
            var_array = var_info['array']
            var[0, :, :] = var_array[:,:]
        except Exception as ex:
            print(f'[ERROR] Variable {var_name} could not be created. Exception: {ex}')
            return None
        return var

    def create_3D_variable_from_info_dict_impl(self,var_info,var_name):
        try:
            var = self.EXTRACT.createVariable(var_name,var_info['data_type'],('satellite_id','satellite_bands', 'rows', 'columns'),fill_value=var_info['fill_value'],zlib=True, complevel=6)
            attrs = var_info['attrs']
            if attrs is not None:
                var.setncatts(attrs)
            var_array = var_info['array']
            var[0, :, :,:] = var_array[:,:,:]
        except Exception as ex:
            print(f'[ERROR] Variable {var_name} could not be created. Exception: {ex}')
            return None
        return var




    def create_2D_variable_general(self, var_name, var_array, window):
        start_idx_y = window[0]
        stop_idx_y = window[1]
        start_idx_x = window[2]
        stop_idx_x = window[3]

        satellite_2d_band = self.EXTRACT.createVariable(var_name, 'f4', ('satellite_id', 'rows', 'columns'),
                                                        fill_value=-999.0, zlib=True, complevel=6)

        if len(var_array.shape) == 2:
            satellite_2d_band[0, :, :] = var_array[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
        elif len(var_array.shape) == 3:
            satellite_2d_band[0, :, :] = var_array[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]


        return satellite_2d_band

    def create_insitu_time_variable(self):
        insitu_time = self.EXTRACT.createVariable('insitu_time', 'f8', ('satellite_id', 'insitu_id',), zlib=True,
                                                  complevel=6)
        insitu_time.units = "Seconds since 1970-1-1"
        insitu_time.description = 'In situ time in ISO 8601 format (UTC).'
        return insitu_time

    def create_insitu_original_bands_variable(self):
        insitu_original_bands = self.EXTRACT.createVariable('insitu_original_bands', 'f4', ('insitu_original_bands'),
                                                            fill_value=-999, zlib=True, complevel=6)
        insitu_original_bands.description = 'In situ bands in nm'
        return insitu_original_bands

    def create_insitu_exact_wavelengths_variable(self):
        insitu_exact_wavelenghts = self.EXTRACT.createVariable('insitu_exact_wavelenghts', 'f4',
                                                               ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                                               fill_value=-999, zlib=True, complevel=6)
        insitu_exact_wavelenghts.description = 'In situ bands in nm'

        return insitu_exact_wavelenghts

    def create_insitu_rrs_variable(self):
        insitu_Rrs = self.EXTRACT.createVariable('insitu_Rrs', 'f4',
                                                 ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                                 fill_value=-999, zlib=True, complevel=6)
        insitu_Rrs.description = 'In situ Rrs'

        return insitu_Rrs

    def create_insitu_time_difference_variable(self):

        time_difference = self.EXTRACT.createVariable('time_difference', 'f4', ('satellite_id', 'insitu_id'),
                                                      fill_value=-999,
                                                      zlib=True, complevel=6)
        time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
        time_difference.units = "seconds"

        return time_difference

    def create_insitu_lat_long_variables(self):

        insitu_lat = self.EXTRACT.createVariable('insitu_latitude', 'f8', ('satellite_id', 'insitu_id'),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        insitu_lat.short_name = "latitude"
        insitu_lat.units = "degrees"

        insitu_lon = self.EXTRACT.createVariable('insitu_longitude', 'f8', ('satellite_id', 'insitu_id'),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        insitu_lon.short_name = "longitude"
        insitu_lon.units = "degrees"

        return insitu_lat, insitu_lon

    def create_insitu_variables_for_single_insitu_data(self):

        insitu_lat = self.EXTRACT.createVariable('insitu_latitude', 'f8', ('satellite_id',),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        insitu_lat.short_name = "latitude"
        insitu_lat.units = "degrees"

        insitu_lon = self.EXTRACT.createVariable('insitu_longitude', 'f8', ('satellite_id',),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        insitu_lon.short_name = "longitude"
        insitu_lon.units = "degrees"

        time_difference = self.EXTRACT.createVariable('time_difference', 'f4', ('satellite_id',),
                                                      fill_value=-999,
                                                      zlib=True, complevel=6)
        time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
        time_difference.units = "seconds"

        insitu_time = self.EXTRACT.createVariable('insitu_time', 'f8', ('satellite_id',), zlib=True,
                                                  complevel=6)
        insitu_time.units = "Seconds since 1970-1-1"
        insitu_time.description = 'In situ time in ISO 8601 format (UTC).'

        return insitu_lat, insitu_lon, insitu_time, time_difference

    def create_insitu_variable_for_single_insitu_data(self, name_var, units, desc):

        name = f'insitu_{name_var}'

        insitu_var = self.EXTRACT.createVariable(name, 'f8', ('satellite_id',),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        insitu_var.units = units
        insitu_var.short_name = name_var
        insitu_var.description = desc

        return insitu_var

    def create_insitu_flag_variable(self, name_var, flag_values, flag_meanings):

        name = f'insitu_flag_{name_var}'
        insitu_var = self.EXTRACT.createVariable(name, 'i4', ('satellite_id',), fill_value=-999, zlib=True, complevel=6)
        insitu_var.short_name = name_var
        insitu_var.flag_masks = flag_values
        insitu_var.flag_meanings = flag_meanings

        return insitu_var

    def create_insitu_flag_variable_version2(self, name_var, flag_meanings):

        flag_values = []
        flag_meanings_str = ' '.join(flag_meanings)
        for idx in range(len(flag_meanings)):
            flag_values.append(int(np.power(2, float(idx))))

        insitu_var = self.EXTRACT.createVariable(name_var, 'i4', ('satellite_id', 'insitu_id'), fill_value=-999,
                                                 zlib=True, complevel=6)
        insitu_var.short_name = name_var
        insitu_var.flag_masks = flag_values
        insitu_var.flag_meanings = flag_meanings_str

        return insitu_var

    def close_file(self):
        self.EXTRACT.close()


# class SatExtractOptions:
#     def  __init__(self,config_file):
#         self.options = None
#         if config_file is not None:
#             try:
#                 self.options = configparser.ConfigParser()
#                 self.options.read(config_file)
#             except:
#                 self.options = None



# def get_basic_options_from_file_config(args,options):
#     # size box
#     size_box = 25
#     if options.has_option('satellite_options', 'extract_size'):
#         size_box = int(options['satellite_options']['extract_size'])
#     # resolution
#     res = 'WFR'
#     if options.has_option('satellite_options', 'resolution'):
#         res = options['satellite_options']['resolution']
#     # path_out
#     path_out = None
#     if options.has_option('file_path', 'output_dir'):
#         path_out = options['file_path']['output_dir']
#     else:
#         if args.output:
#             path_out = args.output
#     if path_out is None:
#         print(f'[ERROR] compulsory option path_out was not defined in the config file or argument output')
#         return None
#     if not os.path.isdir(path_out):
#         path_out = create_dir(path_out)
#         if path_out is None:
#             print(f'[ERROR] path_out: {path_out} does not exist and could not be created')
#             return None
#     # makd_brdb
#     make_brdf = False
#     if options.has_option('satellite_options', 'brdf'):
#         if options['satellite_options']['brdf'].upper() == 'T' or options['satellite_options'][
#             'brdf'].upper() == 'TRUE':
#             make_brdf = True
#     # satellite_path_source
#     satellite_path_source = None
#     if options.has_option('file_path', 'sat_source_dir'):
#         satellite_path_source = options['file_path']['sat_source_dir']
#     else:
#         if args.path_to_sat:
#             satellite_path_source = args.path_to_sat
#     if satellite_path_source is None:
#         print(
#             '[ERROR] compulsory option satellite_path_source was not defined in the config file or argument path_to_sat')
#         return None
#     if not os.path.exists(satellite_path_source):
#         print(f'ERROR path: {satellite_path_source} does not exit')
#         return None
#     # tmp_path
#     tmp_path = None
#     if options.has_option('file_path', 'tmp_dir'):
#         tmp_path = options['file_path']['tmp_dir']
#         if os.path.exists(tmp_path) and os.path.isdir(tmp_path):
#             tmp_path_del = os.path.join(tmp_path, '*')
#             cmd = f'rm -r {tmp_path_del}'
#             prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
#             out, err = prog.communicate()
#             if err:
#                 print(err)
#         if not os.path.exists(tmp_path) or not os.path.isdir(tmp_path):
#             tmp_path = create_dir(tmp_path)
#
#     if tmp_path is None:
#         tmp_path = satellite_path_source
#         print(
#             f'[WARNING] tmp_path was not defined or does not exist. tmp path set to sat_path: {satellite_path_source}')
#     # org
#     org = None
#     if options.has_option('file_path', 'sat_source_dir_organization'):
#         org = options['file_path']['sat_source_dir_organization']
#
#     # Ohter bands
#     extra_bands = None
#     if options.has_option('satellite_options', 'extra_bands'):
#         str_val = options['satellite_options']['extra_bands']
#         extra_bands = [x.strip() for x in str_val.split(',')]
#
#     wce = None
#     if options.has_option('satellite_options','wce'):
#         wce = options['satellite_options']['wce'].strip()
#         wce = f'"{wce}' if not wce.startswith(f'"') else f'{wce}'
#         wce = f'{wce}"' if not wce.endswith(f'"') else f'{wce}'
#
#     basic_options = {
#         'satellite_path_source': satellite_path_source,
#         'path_out': path_out,
#         'tmp_path': tmp_path,
#         'size_box': size_box,
#         'make_brdf': make_brdf,
#         'org': org,
#         'resolution': res,
#         'extra_bands': extra_bands,
#         'wce': wce
#     }
#     return basic_options



def create_dir(path):
    if os.path.isdir(path):
        return path
    try:
        os.mkdir(path)
        return path
    except:
        return None

# def get_basic_options_from_arguments(args):
#     ##parameters with default values
#     size_box = 25
#     make_brdf = False
#     org = 'YYYY/jjj'
#
#     res = 'WFR'
#     if 'resolution' in args:
#         res='WRR' if args.resolution=='WRR' else 'WFR'
#
#     path_out = None
#     if args.output:
#         path_out = args.output
#     if path_out is None:
#         print(f'[ERROR] compulsory option path_out was not defined in the config file or argument output')
#         return None
#     if not os.path.isdir(path_out):
#         path_out = create_dir(path_out)
#         if path_out is None:
#             print(f'[ERROR] path_out: {path_out} does not exist and could not be created')
#             return None
#
#     satellite_path_source = None
#     if args.path_to_sat:
#         satellite_path_source = args.path_to_sat
#     if satellite_path_source is None:
#         print(
#             '[ERROR] compulsory option satellite_path_source was not defined in the config file or argument path_to_sat')
#         return None
#     if not os.path.exists(satellite_path_source):
#         print(f'ERROR path: {satellite_path_source} does not exit')
#         return None
#
#     tmp_path = satellite_path_source
#
#     wce = None
#     if args.wce_expression:
#         wce = args.wce_expression
#         wce = f'"{wce}' if not wce.startswith(f'"') else f'{wce}'
#         wce = f'{wce}"' if not wce.endswith(f'"') else f'{wce}'
#
#
#     basic_options = {
#         'satellite_path_source': satellite_path_source,
#         'path_out': path_out,
#         'tmp_path': tmp_path,
#         'size_box': size_box,
#         'make_brdf': make_brdf,
#         'org': org,
#         'resolution': res,
#         'wce': wce
#     }
#
#     wce = f'"PACE_OCI*.L2.OC_AOP.V2_0.NRT.nc"'  # wild card expression
#
#     return basic_options




def start_extract(ofname, extract_info,verbose):
    newEXTRACT = SatExtract(ofname)
    if not newEXTRACT.FILE_CREATED:
        print(f'[ERROR] File {ofname} could not be created')
        return None

    if verbose:
        print(f'[INFO] Starting file: {ofname}')

    window = extract_info['limits']

    newEXTRACT.set_global_attributes(extract_info['global_at'])

    newEXTRACT.create_dimensions(extract_info['size_box'], extract_info['n_bands'])
    newEXTRACT.create_lat_long_variables(extract_info['lat_array'], extract_info['lon_array'], window)
    newEXTRACT.create_satellite_time_variable(extract_info['satellite_time'][0])

    return newEXTRACT


def get_insitu_site(args,options, path_out):
    in_situ_lat = None
    in_situ_lon = None
    if args.sitename:
        site_name = args.sitename
        in_situ_lat, in_situ_lon = cfs.get_lat_lon_ins(site_name)  # in situ location based on the station name
        if in_situ_lat is None or in_situ_lon is None:
            print(f'[ERROR] {site_name} is not defined in the site list. Geographic coordinates are unknow')
            return None

    elif args.config_file:
        if  options.has_option('Time_and_sites_selection', 'site'):
            site_name = options['Time_and_sites_selection']['site'].strip()
            in_situ_lat, in_situ_lon = cfs.get_lat_lon_ins(site_name)  # in situ location based on the station name
            if in_situ_lat is None or in_situ_lon is None:
                print(f'[ERROR] {site_name} is not defined in the site list. Geographic coordinates are unknow')
                return None
        else:
            site_name = 'UNKNOWN'
            if options.has_option('Time_and_sites_selection', 'in_situ_lat') and options.has_option('Time_and_sites_selection', 'in_situ_lon'):
                try:
                    in_situ_lat = float(options['Time_and_sites_selection']['in_situ_lat'].strip())
                    in_situ_lon = float(options['Time_and_sites_selection']['in_situ_lon'].strip())
                except:
                    pass
            if in_situ_lat is None or in_situ_lon is None:
                print(f'[ERROR] In situ latitude/longitude or a valid site name must be defined in the configuration file')
                return None

    path_out_site = path_out
    if os.path.basename(path_out_site) != site_name:
        path_out_site = os.path.join(path_out, site_name)
    if not os.path.isdir(path_out_site):
        path_out_site = create_dir(path_out_site)
        if path_out_site is None:
            print(f'[ERROR] {path_out_site} does not exist and could not be created')
            return None
    in_situ_site = {
        'site_name': site_name,
        'latitude': in_situ_lat,
        'longitude': in_situ_lon,
        'path_out': path_out_site
    }
    return in_situ_site

# def get_start_end_date_from_args(args):
#     start_date = None
#     end_date = None
#     if args.startdate and args.enddate:
#         try:
#             start_date = dt.strptime(args.startdate, '%Y-%m-%d')
#             end_date = dt.strptime(args.enddate, '%Y-%m-%d')
#             if args.verbose:
#                 print(f'[INFO] Start date: {args.startdate} End date: {args.enddate}')
#         except:
#             print(f'[WARNING] --startdate (-sd) and/or --enddate (-ed) could not be parsed from {args.startdate} and/or {args.enddate}. Format should be: YYYY-mm-dd')
#             start_date = None
#             end_date = None
#     return start_date,end_date

def get_params_time(args,options):
    date_list = None
    if options is not None:
        datetime_start =dt.strptime('2000-01-01', '%Y-%m-%d')
        datetime_end = dt.today()
        if options.has_option('Time_and_sites_selection', 'time_start'):
            try:
                datetime_start =dt.strptime(options['Time_and_sites_selection']['time_start'], '%Y-%m-%d')
            except:
                print(f'WARNING: time_start format is not valid. Usind dafult value: 2000-01-01')

        if options.has_option('Time_and_sites_selection', 'time_end'):
            try:
                datetime_end =dt.strptime(options['Time_and_sites_selection']['time_stop'],
                                                 '%Y-%m-%d') + timedelta(seconds=59, minutes=59, hours=23)
            except:
                print(f'WARNING: time_end format is not valid. Usind dafult value: today')

        if options.has_option('Time_and_sites_selection', 'time_list_file'):
            time_list_file = options['Time_and_sites_selection']['time_list_file']
            if time_list_file:
                date_list, datetime_start, datetime_end = get_date_list_from_file(time_list_file, datetime_start,
                                                                                  datetime_end)

    else:
        if args.startdate:
            datetime_start =dt.strptime(args.startdate, '%Y-%m-%d')
        else:
            datetime_start =dt.strptime('2000-01-01', '%Y-%m-%d')
        if args.enddate:
            datetime_end =dt.strptime(args.enddate, '%Y-%m-%d') + timedelta(seconds=59, minutes=59, hours=23)
        else:
            datetime_end = dt.today()
        if args.date_list_file and os.path.isfile(args.date_list_file):
            date_list, datetime_start, datetime_end = get_date_list_from_file(args.date_list_file, datetime_start,
                                                                              datetime_end)

    if date_list is None:
        date_list = get_date_list_from_start_end_date(datetime_start, datetime_end)
    ndates = len(date_list)
    if args.verbose:
        print(f'[INFO] Start date: {datetime_start}')
        print(f'[INFO] End date: {datetime_end}')
        print(f'[INFO] # of dates: {ndates}')

    return datetime_start, datetime_end, date_list


def get_date_list_from_file(file_list, dt_start, dt_end):
    if not os.path.exists(file_list):
        return None
    date_list = []
    f1 = open(file_list, 'r')
    dt_start_real = None
    dt_end_real = None
    for line in f1:
        dateherestr = line.strip()
        try:
            datehere =dt.strptime(dateherestr, '%Y-%m-%d')
            if dt_start <= datehere <= dt_end:
                date_list.append(datehere)
                if dt_start_real is None:
                    dt_start_real = datehere
                    dt_end_real = datehere
                else:
                    if datehere < dt_start_real:
                        dt_start_real = datehere
                    if datehere > dt_end_real:
                        dt_end_real = datehere
        except:
            pass
    f1.close()
    if len(date_list) == 0:
        return None
    return date_list, dt_start_real, dt_end_real

def get_date_list_from_start_end_date(dt_start, dt_end):
    date_list = []
    dt = dt_start
    while dt <= dt_end:
        date_list.append(dt)
        dt = dt + timedelta(hours=24)
    return date_list

def get_list_products_day(path_source, date_here, wce , org):
    path_source_date = path_source
    orgs = {
        'YYYY':'%Y',
        'mm': '%m',
        'dd': '%d',
        'jjj': '%j'
    }
    if org is not None:
        for o in org.split('/'):
            o = o.strip()
            if o in orgs.keys():
                o=orgs[o]
            path_source_date = os.path.join(path_source_date,date_here.strftime(o))

    cmd = f'find {path_source_date} -name {wce}|sort|uniq'#>> {path_to_list}'
    prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE,stdout = subprocess.PIPE)
    out, err = prog.communicate()
    product_list = []
    for sout in out.decode().split('\n'):
        if len(sout.strip())>0:
            product_list.append(sout)

    return product_list


# def overwrite(options):
#     ow = False
#     if options.has_option('file_path', 'overwrite'):
#         ow = get_boolean_val(options['file_path']['overwrite'])
#     return ow

def get_boolean_val(sval):
    sval = sval.strip().upper()
    if sval == 'TRUE' or sval == 'T' or sval == '1':
        return True
    else:
        return False

def get_int_val(sval):
    sval = sval.strip()
    try:
        val = int(sval)
    except:
        val = None
    return val










def get_seabass_options(options, section):
    var_date = 'date'
    var_time = 'time'
    var_lat = 'lat'
    var_lon = 'lon'
    format_date = '%Y%m%d'
    format_time = '%H:%M:%S'

    if options.has_option(section, 'var_date'):
        var_date = options[section]['var_date'].strip()
    if options.has_option(section, 'var_lat'):
        var_lat = options[section]['var_lat'].strip()
    if options.has_option(section, 'var_lon'):
        var_lon = options[section]['var_lon'].strip()
    if options.has_option(section, 'format_date'):
        format_date = options[section]['format_date'].strip()

    if options.has_option(section, 'var_time'):
        var_time = options[section]['var_time'].strip()
    if options.has_option(section, 'format_time'):
        format_time = options[section]['format_time']

    return var_date, var_time, var_lat, var_lon, format_date, format_time




def get_geo_info(size_box,insitu_lat, insitu_lon, lat, lon):
    contain_flag = 0
    limits = None
    rc = None
    if cfs.contain_location(lat, lon, insitu_lat, insitu_lon) == 1:
        if lat.ndim == 1 and lon.ndim == 1:
            r = np.argmin(np.abs(lat.astype(np.float64) - insitu_lat))
            c = np.argmin(np.abs(lon.astype(np.float64) - insitu_lon))
        else:
            r, c = cfs.find_row_column_from_lat_lon(lat.astype(np.float64), lon.astype(np.float64), insitu_lat, insitu_lon)


        start_idx_x = (c - int(size_box / 2))  # lon
        stop_idx_x = (c + int(size_box / 2) + 1)  # lon
        start_idx_y = (r - int(size_box / 2))  # lat
        stop_idx_y = (r + int(size_box / 2) + 1)  # lat

        if lat.ndim == 1 and lon.ndim == 1:
            if start_idx_y >= 0 and (stop_idx_y + 1) < lat.shape[0] and start_idx_x >= 0 and (stop_idx_x + 1) < \
                    lon.shape[0]:
                contain_flag = 1
        else:
            if start_idx_y >= 0 and (stop_idx_y + 1) < lat.shape[0] and start_idx_x >= 0 and (stop_idx_x + 1) < \
                    lat.shape[1]:
                contain_flag = 1
        if contain_flag == 1:
            limits = [start_idx_y, stop_idx_y, start_idx_x, stop_idx_x]
            rc = [r, c]



    return limits, rc

def get_date_list_from_seabass(sb,var_date,var_time,format_date,format_time,start_date,end_date):
    if not var_date in sb.variables:
        return [None]*2


    date_list_orig = sb.data[var_date]
    time_list_orig = sb.data[var_time] if var_time in sb.variables else None
    format_date_time = f'{format_date}T{format_time}'
    date_array = []
    time_list = []
    for idx,x in enumerate(date_list_orig):

        try:
            date_here = dt.strptime(str(x), format_date)
        except:
            print(f'[ERROR] Error parsing dates in SeaBass file: {x} could not be parsed using {format_date} format.')
            print(
                f'[ERROR] Plase review SEABASS_SELECTION/format_date in the config. file. Expected format: {sb.variables[var_date][1]}')
            return [None]*2
        if start_date is not None and end_date is not None:
            if date_here<start_date or date_here>end_date:
                continue

        date_array.append(date_here.strftime('%Y-%m-%d'))
        if time_list_orig is None:
            time_list.append(dt.strptime(str(x), format_date))
        else:
            try:
                val_s = f'{str(x)}T{str(time_list_orig[idx])}'
                time_list.append(dt.strptime(val_s, format_date_time))
            except:
                print(
                    f'[ERROR] Error parsing dates in SeaBass file: {val_s} could not be parsed using {format_date_time} format.')
                print(
                    f'[ERROR] Plase review SEABASS_SELECTION/format_date in the config. file. Expected format: {sb.variables[var_date][1]}T{sb.variables[var_time][1]}')
                return [None] * 2

    if len(date_array)==0:
        print(f'[WARNING] No data was found for the given temporal range')
        return [None]*2
    date_array = np.array(date_array)
    time_list = np.array(time_list)

    return date_array,time_list

def get_date_list_from_dataframe(df,col_date,format_date,col_time,format_time):
    date_array_ts = df[col_date]
    time_array_ts = None
    if col_time is not None:
        format_datetime = f'{format_date}T{format_time}'
        time_array_ts = df[col_time]
    else:
        format_datetime = format_date
    only_date_array = []
    time_list = []
    for idx,xd in enumerate(date_array_ts):
        x = f'{xd}T{time_array_ts[idx]}' if time_array_ts is not None else xd
        try:
            time_here = dt.strptime(x, format_datetime)
            time_list.append(time_here)
            only_date_array.append(time_here.strftime('%Y-%m-%d'))
        except:
            print(f'[WARNING] {x} could not be parsed with format {format_datetime}')
            continue

    return only_date_array, time_list


def get_satellite_time_from_global_attributes(fproduct):
    dataset = Dataset(fproduct)
    start_date, end_date, sat_time = [None]*3
    if 'time_coverage_start' in dataset.ncattrs() and 'time_coverage_end' in dataset.ncattrs():
        try:
            start_date = dt.strptime(dataset.time_coverage_start, '%d-%b-%Y %H:%M:%S.%f').replace(tzinfo=pytz.utc)
            end_date = dt.strptime(dataset.time_coverage_end, '%d-%b-%Y %H:%M:%S.%f').replace(tzinfo=pytz.utc)
        except:
           pass
        if start_date is None and end_date is None:
            try:
                start_date = dt.strptime(dataset.time_coverage_start, '%Y%m%dT%H%M%SZ').replace(tzinfo=pytz.utc)
                end_date = dt.strptime(dataset.time_coverage_end, '%Y%m%dT%H%M%SZ').replace(tzinfo=pytz.utc)
            except:
                pass

        if start_date is not None and end_date is not None:
            sat_time = (start_date.timestamp() + end_date.timestamp()) / 2
            sat_time = dt.fromtimestamp(sat_time).astimezone(pytz.utc)

    dataset.close()
    return sat_time

def get_satellite_global_atrib_from_options(options):
    compulsory_keys = ['satellite', 'platform', 'sensor', 'res', 'aco_processor', 'proc_version']
    section = 'satellite_options'
    at = {}
    for key in compulsory_keys:
        if options.has_option(section, key):
            at[key] = options[section][key].strip()
        else:
            at[key] = ''

    return at

def add_insitu_global_atrib(at, site, latitude, longitude, other):
    at['site'] = site
    at['in_situ_lat'] = latitude
    at['in_situ_lon'] = longitude
    if other is not None:
        for key in other:
            at[key] = other[key]
    return at

def get_satellite_ref(global_at):
    if 'satellite' in global_at.keys():
        if 'platform' in global_at.keys():
            global_at['satelliteplatform'] = f'{global_at["satellite"]}{global_at["platform"]}'
        else:
            global_at['satelliteplatform'] = global_at['satellite']
    ref_keys = ['satelliteplatform', 'sensor', 'res', 'aco_processor', 'proc_version']
    ref = None
    for key in ref_keys:
        if key in global_at.keys() and len(global_at[key]) > 0:
            if ref is None:
                ref = global_at[key]
            else:
                ref = f'{ref}_{global_at[key]}'
    return ref

def get_path_date(path_base,org,date_here,createIfNotExist):
    if org is None:
        path_date = path_base
    else:
        try:
            path_date = os.path.join(path_base,date_here.strftime(org))
        except:
            print(f'[WARNING] Path date for date {date_here.strftime("%Y-%m-%d")} could not be parsed using {org} organization')
            return None

    if not os.path.isdir(path_date):
        if not createIfNotExist:
            print(f'[WARNING] Path date {path_date} does not exist or is not a valid directory')
            return None
        else:
            path_expected  = path_date
            path_date = path_base
            for time_ref in org.split('/'):
                path_date = os.path.join(path_date,date_here.strftime(time_ref))
                path_date = create_dir(path_date)
                if path_date is None:
                    print(f'[WARNING] Path date {path_expected} does not exist and could not be created. Please review permissions')
                    break
    return path_date

def concatenate_csv(file_list,file_out,remove_files):
    import pandas as pd
    objs = []
    for file in file_list:
        objs.append(pd.read_csv(file,sep=';'))
        if remove_files:
            try:
                os.remove(file)
            except:
                pass
    df = pd.concat(objs)
    df.to_csv(file_out,sep=';',index=False)

class SatExtractBase:
    def __init__(self,type_sat_extract,verbose):
        self.verbose  =verbose
        if type_sat_extract is None:
            type_sat_extract = 'OLCI'

        self.sat_extract_sensor = None
        if type_sat_extract == 'CNR':
            from sat_extract_CNR import SatExtractCNR
            self.sat_extract_sensor = SatExtractCNR(self.verbose)
        elif type_sat_extract == 'OCI':
            from sat_extract_OCI import SatExtractOCI
            self.sat_extract_sensor = SatExtractOCI(self.verbose)
        elif type_sat_extract == 'OLCI':
            from sat_extract_OLCI import SatExtractOLCI
            self.sat_extract_sensor = SatExtractOLCI(self.verbose)
        elif type_sat_extract == 'CMEMS':
            from sat_extract_CMEMS import SatExtractCMEMS
            self.sat_extract_sensor  = SatExtractCMEMS(self.verbose)
        else:
            print(f'[ERROR] Satellite extract sensor class for  {type_sat_extract} is not implemented: Please check SatExtractBase in sat_extract.py')

        self.start_date = None
        self.end_date = None


    def run(self,sat_extract_options,insituBase,ncores):
        ##common for all days for l3 products
        lat_array = None
        lon_array = None
        ##param list is used for ncores>0
        params_list = []
        insituBase.start_date = self.start_date
        insituBase.end_date = self.end_date

        if args.verbose:
            print(f'[INFO] Getting output path...')
        output_path = sat_extract_options.get_output_path()
        if output_path is None:
            return
        if args.verbose:
            print(f'[INFO] Output path: {output_path}')

        if args.verbose:
            print(f'[INFO] Preparing in situ data...')
        try:
            if not insituBase.prepare_data():
                return
        except Exception as ex:
            print(f'[ERROR] prepare_data() method is not available or not valid in class {type(insituBase)}')
            print(f'[ERROR] Exception: {ex}')
            return
        if args.verbose:
            print(f'[INFO] Preparing in situ data: Completed')

        if args.verbose:
            print(f'[INFO] Getting input path info...')
        input_path_info = sat_extract_options.get_input_path_info()
        if input_path_info is None:
            return
        if args.verbose:
            print(f'[INFO] Getting input path info: Completed')

        if args.verbose:
            print(f'[INFO] Getting extract options...')
        extract_options = self.sat_extract_sensor.get_extract_options(sat_extract_options)
        if extract_options is None:
            return
        if args.verbose:
            print(f'[INFO] Getting extract options: Completed')

        if args.verbose:
            print(f'[INFO] Getting overwrite...')
        overwrite = sat_extract_options.get_overwrite()
        if overwrite is None:
            return
        if args.verbose:
            print(f'[INFO] Overwrite: {overwrite}')

        for date_here in insituBase.date_list:
            if self.start_date is not None and self.end_date is not None:
                if date_here<self.start_date or date_here>self.end_date:
                    if self.verbose:
                        print(f'[INFO] Skipping date {date_here.strftime("%Y-%m-%d")} as it not in the range: {self.start_date.strftime("%Y-%m-%d")} to {self.end_date.strftime("%Y-%m-%d")}')
                    continue
            if args.verbose:
                print(f'[INFO] Working with date: {date_here.strftime("%Y-%m-%d")}')
            list_files,satellite_time = self.sat_extract_sensor.get_files_day(date_here, input_path_info,sat_extract_options)
            if list_files is None:
                print(f'[WARNING] Data files for {date_here.strftime("%Y-%m-%d")} could not be retrieved. Skipping...')
                continue
            if satellite_time is None:
                print(f'[WARNING] Satellite time for date could not be defined. Skipping...')
                continue

            #insitu_time,insitu_lat,insitu_lon,insitu_indices = insituBase.get_metadata_date(date_here)
            if insituBase.fixed_site:
                try:
                    insitu_time, insitu_lat, insitu_lon, insitu_indices = insituBase.get_metadata_date_basic(date_here)
                except:
                    print(f'[ERROR] Method get_metadata_date_basic is not available in class {type(insituBase)}')
            else:
                try:
                    insitu_time, insitu_lat, insitu_lon, insitu_indices = insituBase.get_metadata_date(date_here)
                except Exception as ex:
                    print(f'[ERROR] Method get_metadata_date is not available in class {type(insituBase)}. Exception: {ex}')

            path_extract_output_here = os.path.join(output_path,f'{insituBase.get_ref_date(date_here)}_extracts.csv')
            if os.path.exists(path_extract_output_here):
                os.remove(path_extract_output_here)
            extract_info = {
                    'global_at':sat_extract_options.get_satellite_global_atrib(),
                    'list_files': list_files,
                    'insitu_time': insitu_time,
                    'insitu_lat': insitu_lat,
                    'insitu_lon': insitu_lon,
                    'insitu_indices': insitu_indices,
                    'satellite_time': satellite_time,
                    'path_extract_output': path_extract_output_here,
                    'unzip_dir': input_path_info['unzip_dir']
            }
            if lat_array is None and lon_array is None and self.sat_extract_sensor.is_l3_product():
                if self.sat_extract_sensor.sat_type=='CMEMS':
                    pass
                else:
                    lat_array, lon_array = self.sat_extract_sensor.get_lat_lon_arrays(extract_options, list_files[0])
            if ncores == 0:
                self.sat_extract_sensor.run_extract_day(extract_options, extract_info, lat_array, lon_array,output_path, overwrite)
            else:
                if ncores > 0 or ncores == -1:
                    params_list.append([extract_options, extract_info, lat_array, lon_array, output_path, overwrite])


        if ncores > 0 or ncores == -1:
            if len(params_list)==0:
                print(f'[WARNING] No data for running multprocessing.')
                return
            if self.verbose:
                print(f'[INFO] Starting parallel processing. Number of dates: {len(params_list)}')
                print(f'[INFO] CPUs: {os.cpu_count()}')
                print(f'[INFO] Parallel processes: {ncores}')
            npool = os.cpu_count() if ncores < 0 else ncores
            poolhere = Pool(npool)
            poolhere.map(self.sat_extract_sensor.run_parallel_extract_day, params_list)



    def create_sh_slurm(self,sat_extract_options,mp_options,insituBase):
        import COMMON.sbatch_scripter as sbs
        insituBase.start_date = self.start_date
        insituBase.end_date = self.end_date
        output_path = sat_extract_options.get_output_path()
        if output_path is None:
            return
        if args.verbose:
            print(f'[INFO] Preparing in situ data...')
        if not insituBase.prepare_data():
            return
        temp_path = create_dir(os.path.join(output_path, f'Temp_{str(dt.now().timestamp()).replace(".", "_")}'))
        if temp_path is None:
            print(f'[ERROR] Temporary path for sbatch files could not be created. Please review permissions.')
            return
        if self.verbose:
            print(f'[INFO] Temporary path: {temp_path}')

        npool = os.cpu_count() if mp_options['ncores'] < 0 else 1 if mp_options['ncores'] == 0 else mp_options['ncores']
        nincluded = 0
        index_folder = 1
        sbatch_files = []
        sbatch_log_files = []
        dir_code = os.path.dirname(__init__.__file__)
        if insituBase.fixed_site:
            folder_csv = None
            line_py_base = f'python "{dir_code}/sat_extract.py" -sat {self.sat_extract_sensor.sat_type} -insitu {insituBase.insitu_type}'
        else:
            folder_csv = os.path.join(temp_path, f'CSV_{index_folder}')
            create_dir(folder_csv)
            line_py_base = f'python "{dir_code}/sat_extract.py" -sat {self.sat_extract_sensor.sat_type} -insitu MULTIPLE_CSV'
        insitu_date_list = insituBase.date_list

        slurm_start = None
        slurm_end = None

        for insitu_date in insitu_date_list:
            if self.start_date is not None and self.end_date is not None:
                if insitu_date<self.start_date or insitu_date>self.end_date:
                    if self.verbose:
                        print(f'[INFO] Skipping date {insitu_date.strftime("%Y-%m-%d")} as it not in the range: {self.start_date.strftime("%Y-%m-%d")} to {self.end_date.strftime("%Y-%m-%d")}')
                    continue
            if nincluded == npool:
                if self.verbose:
                    print(f'[NFO] Creating sbatch file...')
                file_config = os.path.join(temp_path, f'config_file_{index_folder}.ini')
                sat_extract_options.prepare_config_file_for_sh_slurm(folder_csv,file_config)
                file_sbatch = os.path.join(temp_path, f'sbatch_script_{index_folder}.slurm')
                scripter = sbs.SBATCH_SCRIPTER(file_sbatch)
                scripter.start_script(mp_options, True)
                line_py_here = f'{line_py_base} -sd {slurm_start.strftime("%Y-%m-%d")} -ed {slurm_end.strftime("%Y-%m-%d")} -c'
                scripter.add_line(f'{line_py_here} "{file_config}" -v')
                scripter.close_script()
                sbatch_files.append(file_sbatch)
                sbatch_log_files.append(os.path.join(temp_path, f'sbatch_script_log_{index_folder}.log'))
                # preparing next---
                nincluded = 0
                index_folder = index_folder + 1
                if not insituBase.fixed_site:
                    folder_csv = os.path.join(temp_path, f'CSV_{index_folder}')
                    create_dir(folder_csv)
                slurm_start = None
                slurm_end = None
                #-------------------------------


            if self.verbose:
                print(f'[INFO] Adding date {insitu_date.strftime("%Y%m%d")}')

            if slurm_start is None and slurm_end is None:
                slurm_start = insitu_date
                slurm_end = insitu_date
            else:
                slurm_end = insitu_date
            if not insituBase.fixed_site:
                name_csv = f'{insituBase.get_ref_date(insitu_date)}.csv'
                file_csv_out = os.path.join(folder_csv, name_csv)
                insituBase.create_csv_metadata_for_date(insitu_date,file_csv_out)
            nincluded = nincluded + 1

        # final sbatch file
        if self.verbose:
            print(f'[INFO] Creating final sbatch file...')
        file_config = os.path.join(temp_path, f'config_file_{index_folder}.ini')
        sat_extract_options.prepare_config_file_for_sh_slurm(folder_csv,file_config)
        file_sbatch = os.path.join(temp_path, f'sbatch_script_{index_folder}.slurm')
        scripter = sbs.SBATCH_SCRIPTER(file_sbatch)
        scripter.start_script(mp_options, True)
        line_py_here = f'{line_py_base} -sd {slurm_start.strftime("%Y-%m-%d")} -ed {slurm_end.strftime("%Y-%m-%d")} -c'
        scripter.add_line(f'{line_py_here} "{file_config}" -v')
        scripter.close_script()
        sbatch_files.append(file_sbatch)
        sbatch_log_files.append(os.path.join(temp_path, f'sbatch_script_log_{index_folder}.log'))

        # file out sh
        file_out_sh = os.path.join(temp_path, f'launch_multiple_sbatch.sh')
        sbs.prepare_sh_script_with_multiple_sbatch(file_out_sh, sbatch_files, sbatch_log_files,mp_options['slurm_sh_max_cores'])

        if self.verbose:
            print(f'[INFO] SH file: {file_out_sh} has been created. ')

        if mp_options['slurm_sh_launch']:
            if self.verbose:
                print(f'[INFO] Ready to launch SH file...')
            import subprocess
            cmd = f'sh {file_out_sh}'
            prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            out, err = prog.communicate()
            if err:
                print(f'[ERROR]Error lunching script: {err}')
            elif self.verbose:
                print(f'[INFO] CMD {file_out_sh} have been launched')

        else:
            if self.verbose:
                print(f'[INFO] SH script has not been launched. To do it, use: sh {file_out_sh}')
    # def create_sh_slurm_multiple_csv(self,options, output_path, mp_options):
    #     import COMMON.sbatch_scripter as sbs
    #
    #     temp_path = create_dir(
    #         os.path.join(output_path, f'Temp_{str(dt.now().timestamp()).replace(".", "_")}'))
    #     if temp_path is None:
    #         print(f'[ERROR] Temporary path for sbatch files could not be created. Please review permissions.')
    #         return
    #     if self.verbose:
    #         print(f'[INFO] Temporary path: {temp_path}')
    #     path_csv = options['MULTIPLE_CSV_SELECTION']['path_csv']
    #     if not os.path.isdir(path_csv):
    #         print(f'[ERROR] Path to csv files {path_csv} was not found or is not a valid directory')
    #         return
    #
    #     col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = get_csv_options(
    #         options, 'MULTIPLE_CSV_SELECTION')
    #
    #     npool = os.cpu_count() if mp_options['ncores'] < 0 else 1 if mp_options['ncores'] == 0 else mp_options['ncores']
    #     nincluded = 0
    #     index_folder = 1
    #     folder_csv = os.path.join(temp_path, f'CSV_{index_folder}')
    #     create_dir(folder_csv)
    #     sbatch_files = []
    #     sbatch_log_files = []
    #     dir_code = os.path.dirname(__init__.__file__)
    #     line_py_base = f'python {dir_code}/sat_extract_CNR.py -c'
    #
    #     for name in os.listdir(path_csv):
    #         if not name.endswith('csv'):
    #             continue
    #         try:
    #             file_csv = os.path.join(path_csv, name)
    #             pd.read_csv(file_csv, sep=col_sep)
    #         except:
    #             print(f'[ERROR] File {path_csv} is not a valid csv separated by {col_sep}')
    #             return
    #         if nincluded == npool:
    #             if self.verbose:
    #                 print(f'[NFO] Creating sbatch file...')
    #             file_config = os.path.join(temp_path, f'config_file_{index_folder}.ini')
    #             options.set('MULTIPLE_CSV_SELECTION', 'path_csv', folder_csv)
    #             options.set('multiprocessing', 'use_slurm_sh', 'False')
    #             with open(file_config, 'w') as configw:
    #                 options.write(configw)
    #             file_sbatch = os.path.join(temp_path, f'sbatch_script_{index_folder}.slurm')
    #             scripter = sbs.SBATCH_SCRIPTER(file_sbatch)
    #             scripter.start_script(mp_options, True)
    #             scripter.add_blank_lines(2)
    #             scripter.add_line(f'{line_py_base} {file_config} -v')
    #             scripter.close_script()
    #             sbatch_files.append(file_sbatch)
    #             sbatch_log_files.append(os.path.join(temp_path, f'sbatch_script_log_{index_folder}.log'))
    #             # preparing next---
    #             nincluded = 0
    #             index_folder = index_folder + 1
    #             folder_csv = os.path.join(temp_path, f'CSV_{index_folder}')
    #             create_dir(folder_csv)
    #             ##-------------------------------
    #
    #         if self.verbose:
    #             print(f'[INFO] Adding CSV file {name}')
    #         file_csv_copy = os.path.join(folder_csv, name)
    #         shutil.copy(file_csv, file_csv_copy)
    #         nincluded = nincluded + 1
    #
    #     # final sbatch file
    #     if self.verbose:
    #         print(f'[INFO] Creating final sbatch file...')
    #     file_config = os.path.join(temp_path, f'config_file_{index_folder}.ini')
    #     options.set('MULTIPLE_CSV_SELECTION', 'path_csv', folder_csv)
    #     print(options['MULTIPLE_CSV_SELECTION']['path_csv'])
    #     with open(file_config, 'w') as configw:
    #         options.write(configw)
    #     file_sbatch = os.path.join(temp_path, f'sbatch_script_{index_folder}.slurm')
    #     scripter = sbs.SBATCH_SCRIPTER(file_sbatch)
    #     scripter.start_script(mp_options, True)
    #     scripter.add_line(f'{line_py_base} {file_config} -v')
    #     scripter.close_script()
    #     sbatch_files.append(file_sbatch)
    #     sbatch_log_files.append(os.path.join(temp_path, f'sbatch_script_log_{index_folder}.log'))
    #
    #     # file out sh
    #     file_out_sh = os.path.join(temp_path, f'launch_multiple_sbatch.sh')
    #     sbs.prepare_sh_script_with_multiple_sbatch(file_out_sh, sbatch_files, sbatch_log_files,
    #                                                mp_options['slurm_sh_max_cores'])
    #     print(f'[INFO] SH file: {file_out_sh} has been created.')
    #     if mp_options['slurm_sh_launch']:
    #         import subprocess
    #         cmd = f'sh {file_out_sh}'
    #         prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    #         out, err = prog.communicate()
    #         if err:
    #             print(f'[ERROR]Error lunching script: {err}')

def get_utm_zone(latitude,longitude):
    #from pyproj import CRS
    #longitude = -8.85
    #latitude = 42.59
    # Determine the UTM zone number
    zone_number = int((longitude + 180) / 6) + 1
    ZONE_LETTERS = "CDEFGHJKLMNPQRSTUVWXX"
    zone_letter = ZONE_LETTERS[int(latitude + 80) >> 3]
    # Determine the hemisphere
    #hemisphere = 'north' if latitude >= 0 else 'south'
    #utm_crs = CRS.from_dict({'proj': 'utm', 'zone': zone_number, 'south': hemisphere == 'south'})
    return f'{zone_number}{zone_letter}'


def main():
    print('[INFO] Creating satellite extracts')
    if not args.config_file:
        return
    options = SatExtractOptions(args.config_file, args.verbose)
    if not options.is_valid:
        if not options.gmanager.is_valid():
            print('[ERROR] Problem processing global options')
        if not options.omanager.is_valid():
            print(f'[ERROR] Problem processing the configuration file {args.config_file}')
        return


    mp_options = options.get_multiprocessing_options()
    if mp_options is None:
        return
    if args.verbose:
        print(f'[INFO] Multiprocessing options:')
        for mp in mp_options:
            print(f'[INFO]   {mp} -> {mp_options[mp]}')


    sat_type = options.check_sat_type(args.sat_type)
    if sat_type is None:
        return
    satExtractBase = SatExtractBase(sat_type, args.verbose)
    if satExtractBase.sat_extract_sensor is None:
        return


    insitu_type, insitu_options = options.check_insitu_options(args.insitu_type)
    if insitu_type is None:
        return

    insituBase = ibase.get_insitu_object(insitu_type, insitu_options, args.verbose)
    if insituBase is None:
        return

    satExtractBase.start_date, satExtractBase.end_date = arf.get_start_end_date_from_args(args)


    if mp_options['use_slurm_sh']:
        satExtractBase.create_sh_slurm(options,mp_options,insituBase)
    else:
        satExtractBase.run(options,insituBase,mp_options['ncores'])



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Satellite extracts from multiple files (one file for variable) available in the CNR server.")
    parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
    parser.add_argument('-c', "--config_file", help="Config File.", required=True)
    parser.add_argument('-sd', "--start_date", help="The Start Date - format YYYY-MM-DD ")
    parser.add_argument('-ed', "--end_date", help="The End Date - format YYYY-MM-DD ")
    parser.add_argument('-insitu',"--insitu_type",help="In situ type")
    parser.add_argument('-sat', "--sat_type", help="Satellite type")
    parser.add_argument('-no_concat', "--no_concatenate", help="Use internaly for sbatch mode", action="store_true")
    parser.add_argument('-make_concat', "--make_concatenate",
                        help="Use internaly for sbatch mode to make final concatenation", action="store_true")
    parser.add_argument('-p', "--product_file", help="Image file.")
    args = parser.parse_args()

    main()
    # from sat_extract_OLCI import SatSourceOlci
    # dir_base = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/TARA_WORK'
    # file_olci = os.path.join(dir_base,'S3A_OL_2_WFR____20230519T101357_20230519T101657_20230520T221344_0179_099_065_1980_MAR_O_NT_003.SEN3')
    # osource = SatSourceOlci(file_olci,args.verbose)
    # # print(osource.timeliness,osource.collection)
    # # polygon = osource.get_polygon_from_manifest()
    # # print(polygon)
    # osource.get_geometry_array('OZA',[190,215,73,98])
