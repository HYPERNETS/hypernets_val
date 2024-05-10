import pytz

from MDB_reader.MDBFile import MDBFile
import os
import numpy as np
from datetime import datetime as dt
from netCDF4 import Dataset
import configparser
from MDB_reader.OptionsManager import OptionsManager
import math


class FlagBuilder:

    def __init__(self, path_mdb_file, config_file, options):
        self.path_mdb_file = path_mdb_file
        self.mfile = None
        self.VALID = False
        if path_mdb_file is not None and os.path.isfile(path_mdb_file):
            try:
                self.mfile = MDBFile(path_mdb_file)
                self.VALID = self.mfile.VALID
                if not self.VALID:
                    print(
                        f'[WARNING] {path_mdb_file} is not a MDB file, it will be used a regular NetCDF file, some functions could not work')
                    self.mfile = None
                    self.nc_file = path_mdb_file
                    self.VALID = True
            except:
                self.VALID = False

        if not self.VALID:
            print(f'[ERROR] {path_mdb_file} is not a valid  file')

        self.omanager = OptionsManager(config_file, options)
        self.flag_list = self.omanager.get_section_list(None)
        if self.flag_list is None:
            self.VALID = False
            print(f'[ERROR] Flag configuration could not be started')

        if self.VALID:
            print(f'[INFO] Configuration file with flags correctly started, including:')
            for f in self.flag_list:
                print(f'[INFO] {f}')
            print(f'[INFO]----------------------------------')

        # flag_insitu_spatial
        self.flag_options = {
            'type': {'type_param': 'str', 'list_values': ['spatial', 'temporal', 'ranges']},
            'typevirtual': {'type_param': 'str', 'list_values': ['spatial', 'temporal', 'ranges']},
            'use_pow2_vflags': {'type_param': 'boolean', 'default': False},
            'var_limit_to_valid': {'type_param': 'str', 'default': None},
            'limit_to_central_pixel': {'type_param': 'boolean', 'default': False},
            'flag_spatial_index': {'type_group': 'spatial', 'type_param': 'strlist'},
            'lat_variable': {'type_group': 'spatial', 'type_param': 'str', 'default': 'insitu_latitude'},
            'lon_variable': {'type_group': 'spatial', 'type_param': 'str', 'default': 'insitu_longitude'},
            'time_variable': {'type_group': 'temporal', 'type_param': 'str', 'default': 'insitu_time'},
            'time_ini': {'type_group': 'temporal', 'type_param': 'str', 'default': None},
            'time_fin': {'type_group': 'temporal', 'type_param': 'str', 'default': None},
            'time_flag_type': {'type_group': 'temporal', 'type_param': 'str', 'default': 'ranges',
                               'list_values': ['ranges', 'yearjday', 'yearmonth', 'jday', 'month', 'yearmonthday',
                                               'monthday', 'flag_satellite_yearmonthday']},
            'var_ranges': {'type_group': 'ranges', 'type_param': 'str'},
            'flag_ranges_indexm': {'type_group': 'ranges', 'type_param': 'strlist'}
        }

    def get_virtual_flag_list(self):
        return self.omanager.get_virtual_flag_list()

    def get_virtual_flags_options(self):
        vf_options = {}
        fl = self.get_virtual_flag_list()
        if fl is not None:
            for f in fl:
                vf_options[f] = self.get_options_dict(f)
        return vf_options

    def get_options_dict(self, flag_ref):
        options_dict = self.omanager.read_options_as_dict(flag_ref, self.flag_options)
        return options_dict

    def create_flag_array(self, flag_ref, create_copy):
        options_dict = self.omanager.read_options_as_dict(flag_ref, self.flag_options)
        # print(options_dict)
        type = options_dict['type']
        if type is None:
            type = options_dict['typevirtual']
        if type is None:
            return
        if type == 'spatial':
            array, dims, flag_names, flag_values = self.create_flag_array_spatial(options_dict)
            if create_copy:
                self.create_copy_with_flag_band(flag_ref, array, flag_names, flag_names, dims)

        if type == 'temporal':
            array, dims, flag_names, flag_values = self.create_flag_array_temporal(options_dict)
            if create_copy:
                self.create_copy_with_flag_band(flag_ref, array, flag_names, flag_names, dims)

        if type == 'ranges':
            array, dims, flag_names, flag_values = self.create_flag_array_ranges(options_dict)
            if create_copy:
                self.create_copy_with_flag_band(flag_ref, array, flag_names, flag_names, dims)

        return flag_values, flag_names, array

    def create_flag_array_ranges_v2(self, options_dict):

        default_value = 0
        flag_names = []
        flag_values = []
        fl_info = {}
        array = None
        nscans = 1
        dataset = Dataset(self.nc_file)
        for fl in options_dict:
            if fl.startswith('flag_ranges'):
                if options_dict[fl]['is_default']:
                    default_value = int(options_dict[fl]['flag_value'])
                value = options_dict[fl]['flag_value']
                if value not in flag_values:
                    flag_values.append(value)
                    flag_names.append(options_dict[fl]['flag_name'])
                value_s = str(options_dict[fl]['flag_value'])
                if value_s not in fl_info.keys():
                    fl_info[value_s] = [fl]
                else:
                    fl_info[value_s].append(fl)
                if array is None:
                    r_array = np.array(dataset.variables[options_dict[fl]['flag_var']][:])
                    array = r_array.copy().astype(np.int64)
                if array is not None and len(dataset.variables[options_dict[fl]['flag_var']].dimensions)==2:

                    r_array = np.array(dataset.variables[options_dict[fl]['flag_var']][:])
                    array = r_array.copy().astype(np.int64)
                    nscans = array.shape[1]

        array[:] = default_value

        for svalue in fl_info:
            value = int(svalue)
            nranges = len(fl_info[svalue])
            if nranges == 1:
                fl = fl_info[svalue][0]
                r_array = np.array(dataset.variables[options_dict[fl]['flag_var']][:])
                if nscans>1 and len(r_array.shape)==1:
                    r_array = self.replicate_scans(r_array,nscans)
                v_min = options_dict[fl]['min_range']
                v_max = options_dict[fl]['max_range']
                if v_min is not None and v_max is not None:
                    array[np.logical_and(r_array >= v_min, r_array <= v_max)] = value
                elif v_min is None and v_max is not None:
                    array[r_array <= v_max] = value
                elif v_min is not None and v_max is None:
                    array[r_array >= v_min] = value
            else:

                array_check = np.zeros(array.shape)
                for fl in fl_info[svalue]:
                    r_array = np.array(dataset.variables[options_dict[fl]['flag_var']][:])
                    if nscans > 1 and len(r_array.shape) == 1:
                        r_array = self.replicate_scans(r_array,nscans)
                    v_min = options_dict[fl]['min_range']
                    v_max = options_dict[fl]['max_range']

                    if v_min is not None and v_max is not None:
                        array_check[np.logical_and(r_array >= v_min, r_array <= v_max)] = array_check[np.logical_and(
                            r_array >= v_min, r_array <= v_max)] + 1
                    elif v_min is None and v_max is not None:
                        array_check[r_array <= v_max] = array_check[r_array <= v_max] + 1
                    elif v_min is not None and v_max is None:
                        array_check[r_array >= v_min] = array_check[r_array >= v_min] + 1
                condition = options_dict[fl]['flag_condition']

                if condition=='and':
                    array[array_check == nranges] = value
                elif condition=='or':
                    array[array_check>0] = value
        dataset.close()

        return array, flag_names, flag_values


    def replicate_scans(self,array,nscans):
        ndata = array.shape[0]
        new_array = np.zeros((ndata,nscans))
        for iscan in range(nscans):
            new_array[:,iscan] = array[:]
        return new_array
    def create_flag_array_ranges(self, options_dict):
        var_ranges = options_dict['var_ranges']
        r_array = self.mfile.get_full_array(var_ranges)
        array = r_array.copy().astype(np.int32)
        default_value = 0
        flag_names = []
        flag_values = []
        for fl in options_dict:
            if fl.startswith('flag_ranges'):
                if options_dict[fl]['is_default']:
                    default_value = int(options_dict[fl]['flag_value'])
                flag_values.append(options_dict[fl]['flag_value'])
                flag_names.append(options_dict[fl]['flag_name'])
        ##Default value
        array[:] = default_value
        ##Flag limites by lat-long values
        for fl in options_dict:
            if fl.startswith('flag_ranges') and not options_dict[fl]['is_default']:
                v_min = options_dict[fl]['min_range']
                v_max = options_dict[fl]['max_range']
                value = int(options_dict[fl]['flag_value'])
                array[np.logical_and(r_array >= v_min, r_array < v_max)] = value
        dims = self.mfile.get_dims_variable(var_ranges)

        if default_value != 0:
            fillValue = self.mfile.get_fill_value(var_ranges)
            if fillValue is not None:
                array[r_array == fillValue] = 0

        # for idx in range(1800):
        #     print(idx,':',r_array[0][idx],'-->',array[0][idx],' with default: ',default_value)

        return array, dims, flag_names, flag_values

    def create_flag_array_temporal(self, options_dict):

        var_time = options_dict['time_variable']

        flag_type = options_dict['time_flag_type']
        use_pow2_vflags = options_dict['use_pow2_vflags']
        time_array = self.mfile.get_full_array(var_time)
        array = time_array.copy().astype(np.int32)
        flag_names = []
        flag_values = []

        # default_value = 0

        orig_shape = array.shape
        if len(orig_shape) > 1:
            array1D = array.flatten()
        else:
            array1D = array

        ##limit to central pixels
        array1Dcentral = None
        if options_dict['limit_to_central_pixel']:
            array_central = self.mfile.get_full_array('insitu_spatial_index')
            if array_central.shape == orig_shape:
                if len(orig_shape) > 1:
                    array1Dcentral = array_central.flatten()
                else:
                    array1Dcentral = array_central

        ##limit to valid
        array1Dvalid = None
        var_limit_to_valid = options_dict['var_limit_to_valid']
        if var_limit_to_valid is not None:
            array_valid = self.mfile.get_full_array(var_limit_to_valid)
            if array_valid.shape == orig_shape:
                if len(orig_shape) > 1:
                    array1Dvalid = array_valid.flatten()
                else:
                    array1Dvalid = array_valid

        # print('°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°')
        # array1Dvalid[6990] = 1

        ##time ini + time_fin
        time_ini = None
        time_fin = None
        if options_dict['time_ini'] is not None and options_dict['time_fin'] is not None:
            time_ini = dt.strptime(options_dict['time_ini'], '%Y-%m-%d %H:%M').replace(tzinfo=pytz.utc)
            time_fin = dt.strptime(options_dict['time_fin'], '%Y-%m-%d %H:%M').replace(tzinfo=pytz.utc)
        fvalue = 0
        for idx in range(len(array1D)):
            val = array1D[idx]

            if val < 0:
                array1D[idx] = 0
                continue

            if array1Dvalid is not None:
                if array1Dvalid[idx] != 1:
                    continue

            if array1Dcentral is not None:
                if array1Dcentral[idx] != 0:
                    continue
            time_here = dt.utcfromtimestamp(val).replace(tzinfo=pytz.utc)

            if time_ini is not None and time_fin is not None:
                if time_here < time_ini or time_here > time_fin:
                    continue
            # if time_ini is not None and time_fin is not None:
            #     if time_ini < time_here < time_fin:
            #         continue

            if flag_type == 'yearjday':
                time_here_str = time_here.strftime('%Y%j')
            elif flag_type == 'yearmonth':
                time_here_str = time_here.strftime('%Y%m')
            elif flag_type == 'jday':
                time_here_str = time_here.strftime('%j')
            elif flag_type == 'month':
                time_here_str = time_here.strftime('%m')
            elif flag_type.endswith('yearmonthday'):
                time_here_str = time_here.strftime('%Y-%m-%d')
            elif flag_type == 'monthday':
                time_here_str = time_here.strftime('%m-%d')

            if flag_type.startswith('flag_satellite'):
                array_flag_satellite = self.mfile.get_full_array('flag_satellite')
                flag_satellite_here = array_flag_satellite[idx]
                if flag_satellite_here == 1:
                    time_here_str = f'S3A_{time_here_str}'
                if flag_satellite_here == 2:
                    time_here_str = f'S3B_{time_here_str}'

            if time_here_str not in flag_names:
                flag_names.append(time_here_str)
                vflag = self.get_flag_value(fvalue, use_pow2_vflags)
                flag_values.append(vflag)
                print(f'[INFO] Adding temporal flag: {time_here_str} with value: {vflag}----->{val}')
                fvalue = fvalue + 1
            else:
                index_flag = flag_names.index(time_here_str)
                vflag = flag_values[index_flag]

            array1D[idx] = vflag

        if len(orig_shape) > 1:
            array = array1D.reshape(orig_shape)
        else:
            array = array1D

        dims = self.mfile.get_dims_variable(var_time)
        # if default_value != 0:
        #     array[lat_array == -999.0] = 0

        return array, dims, flag_names, flag_values

    def get_flag_value(self, vorig, use_pow2_vflags):
        if use_pow2_vflags:
            fvalue = int(math.pow(2, vorig))
        else:
            fvalue = int(vorig + 1)
        return fvalue

    def create_flag_array_spatial(self, options_dict):
        var_latitude = options_dict['lat_variable']
        var_longitude = options_dict['lon_variable']
        lat_array = self.mfile.get_full_array(var_latitude)
        lon_array = self.mfile.get_full_array(var_longitude)
        array = lat_array.copy().astype(np.int32)
        default_value = 0
        flag_names = []
        flag_values = []
        for fl in options_dict:
            if fl.startswith('flag_spatial'):
                if options_dict[fl]['is_default']:
                    default_value = int(options_dict[fl]['flag_value'])
                flag_values.append(options_dict[fl]['flag_value'])
                flag_names.append(options_dict[fl]['flag_name'])
        ##Default value
        array[:] = default_value
        ##Flag limites by lat-long values
        for fl in options_dict:
            if fl.startswith('flag_spatial') and not options_dict[fl]['is_default']:
                lat_min = options_dict[fl]['lat_min']
                lat_max = options_dict[fl]['lat_max']
                lon_min = options_dict[fl]['lon_min']
                lon_max = options_dict[fl]['lon_max']
                value = int(options_dict[fl]['flag_value'])
                array[np.logical_and(np.logical_and(lat_array >= lat_min, lat_array <= lat_max),
                                     np.logical_and(lon_array >= lon_min, lon_array <= lon_max))] = value
        dims = self.mfile.get_dims_variable(var_longitude)
        if default_value != 0:
            array[lat_array == -999.0] = 0

        return array, dims, flag_names, flag_values

    def create_copy_with_flag_band(self, name_flag, array, flag_list, flag_values, dims):
        folder = os.path.dirname(self.path_mdb_file)
        time_now = dt.now().strftime('%Y%m%d%H%M%S')
        temp_file = os.path.join(folder, f'Temp_{time_now}.nc')
        ncout = Dataset(temp_file, 'w', format='NETCDF4')
        file_nc = Dataset(self.path_mdb_file, 'r')

        # copy global attributes all at once via dictionary
        ncout.setncatts(file_nc.__dict__)

        # copy dimensions
        for name, dimension in file_nc.dimensions.items():
            ncout.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))

        # copy variables
        for name, variable in file_nc.variables.items():
            fill_value = None
            if '_FillValue' in list(file_nc.ncattrs()):
                fill_value = variable._FillValue

            ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                                 complevel=6)
            # copy variable attributes all at once via dictionary
            ncout[name].setncatts(file_nc[name].__dict__)

            ncout[name][:] = file_nc[name][:]

        # flag_variable
        if name_flag in file_nc.variables:
            print(f'[INFO] Flag name {name_flag} is already in the original file')
        else:
            var = ncout.createVariable(name_flag, 'i4', dims, zlib=True, shuffle=True, complevel=6)

        flag_meanings = ' '.join(flag_list)
        var.flag_meanings = flag_meanings
        var.flag_values = flag_values
        var[:] = array[:]

        file_nc.close()
        ncout.close()
        # os.rename(temp_file,self.path_mdb_file)
