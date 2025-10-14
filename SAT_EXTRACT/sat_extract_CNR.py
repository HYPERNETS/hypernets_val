import argparse, os, sys, shutil,pytz,__init__,warnings
import pandas as pd
from netCDF4 import Dataset
from datetime import datetime as dt
import numpy.ma as ma
import numpy as np
from multiprocessing import Pool

# import user defined functions from other .py
import sat_extract as sextract

code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
import COMMON.common_functions as cfs
warnings.simplefilter("ignore",UserWarning)

class SatExtractCNR:
    def __init__(self,verbose):
        self.verbose = verbose
        self.sat_type = 'CNR'

    ##Method to retrieve file for a single date. Deprecated: get_cmems_multiple_product_day
    def get_files_day(self,datehere, input_path_info, sat_extract_options):


        extract_options = self.get_extract_options(sat_extract_options)

        path_source = input_path_info['path_source']
        org = input_path_info['org']
        path_day = sextract.get_path_date(path_source, org, datehere, False)
        if path_day is None:
            return [None]*2

        dataset_name_format_date = extract_options['dataset_name_format_date']
        dataset_var_list = extract_options['dataset_var_list']
        dataset_name_file = extract_options['dataset_name_file']
        #strdate = datehere.strftime(dataset_name_format_date)
        ncfiles = []

        for var in dataset_var_list:
            name = cfs.get_name_for_date(dataset_name_file,dataset_name_format_date,datehere)
            #name = dataset_name_file
            #name = name.replace('$DATE$', strdate)
            var = var.replace('.', '_')
            name = name.replace('$BAND$', var)
            fname = os.path.join(path_day, name)
            if os.path.exists(fname):
                ncfiles.append(fname)
            else:
                print(f'[WARNING] File: {fname} was not found')

        if len(ncfiles) == len(dataset_var_list):
            satellite_time = [self.get_satellite_time(sat_extract_options,ncfiles[0],datehere)]
            return ncfiles,satellite_time
        else:
            return [None]*2

    def is_l3_product(self):
        return True

    def get_lat_lon_arrays(self,extract_options,file_ref):
        nc_sat = Dataset(file_ref, 'r')
        var_lat = extract_options['lat_dim']
        var_lon = extract_options['lon_dim']
        lat = nc_sat.variables[var_lat][:]
        lon = nc_sat.variables[var_lon][:]
        nc_sat.close()

        return lat, lon

    def get_extract_options(self,sat_extract_options):

        options_out = sat_extract_options.get_satellite_options('CNR')
        if options_out is None:
            return None
        if options_out['dataset_var_list_out'] is None:
            options_out['dataset_var_list_out'] = options_out['dataset_var_list']

        n_bands = 0
        rrs_list = []
        is_reflectance = True
        for var in options_out['dataset_var_list']:
            try:
                value = float(var)
                rrs_list.append(value)
            except:
                is_reflectance = False
        if is_reflectance:
            n_bands = len(options_out['dataset_var_list'])

        options_out['n_bands'] = n_bands
        options_out['is_reflectance'] = is_reflectance
        options_out['size_box'] = sat_extract_options.get_box_size()
        options_out['rrs_list'] = rrs_list

        if options_out['otherfiles_var_list'] is not None:
            str_list = options_out['otherfiles_var_list']
            ofvl = {}
            for str_val in str_list:
                str_val_s = [x.strip() for x in str_val.split(':')]
                if len(str_val_s)!=2:
                    print(f'[ERROR] Option otherfiles_var_list in CNR section should be a comma-separated list of namefile:varlist tokens. {str_val} is not correct')
                    return None
                ofvl[str_val_s[0]] = [x.strip() for x in str_val_s[1].split(';')]
            options_out['otherfiles_var_list'] = ofvl

        return options_out

    def get_satellite_time(self,options, fproduct, daydate):
        sat_options = options.get_general_options('satellite_options')
        cmems_time = sat_options['satellite_time'] if (sat_options is not None and 'satellite_time' in sat_options) else '12:00'

        satellite_time = sextract.get_satellite_time_from_global_attributes(fproduct)
        if satellite_time is None:
            try:
                satellite_time = dt.strptime(f'{daydate.strftime("%Y%m%d")}T{cmems_time}',
                                             '%Y%m%dT%H:%M').replace(tzinfo=pytz.utc)
            except:
                print(f'{cmems_time} is not a valid satellite time option. Skipping')
                return None
        return satellite_time

    def run_parallel_extract_day(self,params):
        self.run_extract_day(params[0], params[1], params[2], params[3], params[4], params[5])

    def run_extract_day(self,extract_options, extract_info, lat_array, lon_array, output_path, overwrite):
        insitu_lat = extract_info['insitu_lat']
        insitu_lon = extract_info['insitu_lon']
        insitu_time = extract_info['insitu_time']
        insitu_indices = extract_info['insitu_indices']

        path_extract_output = extract_info['path_extract_output']
        few = open(path_extract_output, 'w')
        few.write('extract_file;insitu_index;insitu_time;insitu_lat;insitu_lon')

        ntimes = len(insitu_time)

        site_list = []
        format_datetime = '%Y-%m-%dT%H:%M:%S'

        for itime in range(ntimes):

            limits, rc = sextract.get_geo_info(extract_options['size_box'], insitu_lat[itime], insitu_lon[itime],lat_array, lon_array)
            if limits is None:
                print(f'[WARNING] In situ location out of the limits of the satellite product. Skipping...')
                continue

            global_at = extract_info['global_at'].copy()
            datehere = extract_info['satellite_time'][0]
            datehere_str = datehere.strftime('%Y%m%d')
            site = f'{sextract.get_satellite_ref(global_at)}_{datehere_str}_{rc[0]}_{rc[1]}'
            ofname = os.path.join(output_path, f'extract_{site}.nc')

            if os.path.exists(ofname) and not overwrite:
                few.write('\n')
                few.write(
                    f'extract_{site}.nc;{insitu_indices[0][itime]};{insitu_time[itime].strftime(format_datetime)};{insitu_lat[itime]};{insitu_lon[itime]}')
                print(
                    f'[WARNING] [{itime + 1}/{ntimes}] Satellite extract extract_{site}.nc already exists. {itime + 1}/{ntimes} Skiping...')
                continue

            if overwrite and site in site_list:
                few.write('\n')
                few.write(
                    f'extract_{site}.nc;{insitu_indices[0][itime]};{insitu_time[itime].strftime(format_datetime)};{insitu_lat[itime]};{insitu_lon[itime]}')
                print(
                    f'[WARNING] [{itime + 1}/{ntimes}] Satellite extract for {site}  has been already done.  Skiping...')
                continue

            if self.verbose:
                print(f'[INFO] [{itime + 1}/{ntimes}] Preparing extract for site {site}...')

            global_at_here = sextract.add_insitu_global_atrib(global_at, site, insitu_lat[itime], insitu_lon[itime],
                                                                 None)
            extract_info_here = {
                'global_at': global_at_here,
                'limits': limits,
                'size_box': extract_options['size_box'],
                'lat_array': lat_array,
                'lon_array': lon_array,
                'satellite_time': extract_info['satellite_time'],
                'n_bands': extract_options['n_bands'],
                'list_files': extract_info['list_files']
            }
            newExtract = sextract.start_extract(ofname, extract_info_here,self.verbose)
            if extract_options['is_reflectance']:
                newExtract = self.add_reflectance_multiple(newExtract, extract_info_here, extract_options['rrs_list'])
            else:
                newExtract = self.add_variable_multiple(newExtract, extract_info_here,extract_options['dataset_var_list'],extract_options['dataset_var_list_out'])

            ##variables in other files
            if extract_options['otherfiles_var_list'] is not None:
                ofvl = extract_options['otherfiles_var_list']
                for name_of in ofvl:
                    if  newExtract is None:
                        break
                    name_here = cfs.get_name_for_date(name_of, extract_options['dataset_name_format_date'],datehere)
                    file_here = os.path.join(os.path.dirname(extract_info_here['list_files'][0]), name_here)
                    if not os.path.exists(file_here):
                        print(f'[ERROR] Error getting variables from the extra file: {name_here}')
                        newExtract.close_file()
                        newExtract = None
                        continue
                    list_var = ofvl[name_of]
                    if len(list_var) == 1 and list_var[0] == '*':
                        list_var = self.get_list_nondim_variables(file_here)
                    list_var_end = []
                    for var in list_var:
                        var_n = var if var.startswith('satellite_') else f'satellite_{var}'
                        if var_n not in newExtract.EXTRACT.variables:
                            list_var_end.append(var)
                    self.create_other_files_variables(newExtract, extract_info_here, file_here,list_var_end)

            if newExtract is None:
                print(f'[ERROR] Error creating extract for site {site}')
                os.remove(ofname)
                few.write('\n')
                few.write(
                    f'noextract;{insitu_indices[0][itime]};{insitu_time[itime].strftime(format_datetime)};{insitu_lat[itime]};{insitu_lon[itime]}')
                continue
            else:
                print(f'[INFO] Extract {ofname} completed.')
                site_list.append(site)
                few.write('\n')
                few.write(
                    f'extract_{site}.nc;{insitu_indices[0][itime]};{insitu_time[itime].strftime(format_datetime)};{insitu_lat[itime]};{insitu_lon[itime]}')
        few.close()

    def create_other_files_variables(self, newEXTRACT,extract,file_nc,var_list):
        limits = extract['limits']
        start_idx_y = limits[0]
        stop_idx_y = limits[1]
        start_idx_x = limits[2]
        stop_idx_x = limits[3]
        for variable_in in var_list:
            variable_out = variable_in if variable_in.startswith('satellite_') else f'satellite_{variable_in}'
            nc_in = Dataset(file_nc, 'r')
            array = np.ma.squeeze(nc_in.variables[variable_in][:])[start_idx_y:stop_idx_y,start_idx_x:stop_idx_x]
            attrs = cfs.get_attrs_variable(nc_in.variables[variable_in])
            var_info = {
                'data_type': 'f4',
                'fill_value': -999.0,
                'array': array,
                'attrs': attrs
            }
            var = newEXTRACT.create_2D_variable_from_info_dict_impl(var_info, variable_out)
            nc_in.close()
            if var is None:
                print(f'[ERROR] Error creating variables from file {os.path.basename(file_nc)}')
                return False
        return True



    def get_list_nondim_variables(self,file_nc):
        var_list = []
        nc_sat = Dataset(file_nc, 'r')
        for name_var in nc_sat.variables:
            dims = nc_sat.variables[name_var].get_dims()
            if len(dims)>1:
                var_list.append(name_var)
        nc_sat.close()
        return var_list


    def add_reflectance_multiple(self,newEXTRACT, extract, wl_list):
        if not 'satellite_bands' in newEXTRACT.EXTRACT.dimensions:
            print(f'[ERROR] Dimension satellite bands is not defined')
            return
        if 'global_at' in extract:
            global_at = extract['global_at']
        if '1' in extract:
            global_at = extract['1']['global_at']
        list_files = extract['list_files']
        limits = extract['limits']
        start_idx_y = limits[0]
        stop_idx_y = limits[1]
        start_idx_x = limits[2]
        stop_idx_x = limits[3]

        nwl = len(list_files)

        if 'satellite_Rrs' in newEXTRACT.EXTRACT.variables:
            satellite_Rrs = newEXTRACT.EXTRACT.variables['satellite_Rrs']
        else:
            satellite_Rrs = newEXTRACT.create_rrs_variable(global_at['sensor'])
        wavelengths = []
        for iwl in range(nwl):
            wl = float(wl_list[iwl])
            wavelengths.append(float(wl))
            input_dataset = Dataset(list_files[iwl])
            for name, variable in input_dataset.variables.items():
                # wls = str(wl).replace('.', '_')
                wls = f'{wl:.2f}'
                wls = wls.replace('.', '_')
                if wls.endswith('_00'):
                    wls = wls[:-3]
                if wls.find('_') > 0 and wls.endswith('0'):
                    wls = wls[:-1]
                ifind = name.find(wls)
                if ifind >= 0:
                    if variable.ndim == 3:
                        bandarray = ma.array(variable[:, :, :])
                        satellite_Rrs[0, iwl, :, :] = bandarray[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
                    elif variable.ndim == 2:
                        bandarray = ma.array(variable[:, :])
                        satellite_Rrs[0, iwl, :, :] = bandarray[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
            input_dataset.close()

        if 'satellite_bands' in newEXTRACT.EXTRACT.variables:
            satellite_bands = newEXTRACT.EXTRACT.variables['satellite_bands']
            satellite_bands[:] = wavelengths
        else:
            newEXTRACT.create_satellite_bands_variable(wavelengths)

        return newEXTRACT

    def add_variable_multiple(self,newEXTRACT, extract, variable_list, variable_list_out):
        list_files = extract['list_files']
        limits = extract['limits']
        start_idx_y = limits[0]
        stop_idx_y = limits[1]
        start_idx_x = limits[2]
        stop_idx_x = limits[3]
        nvar = len(list_files)
        for idx in range(nvar):
            file_in = list_files[idx]
            variable_in = variable_list[idx]
            variable_out = f'satellite_{variable_list_out[idx]}'
            nc_in = Dataset(file_in, 'r')
            if variable_in in nc_in.variables:
                var_in = nc_in.variables[variable_in]
            elif variable_in.upper() in nc_in.variables:
                var_in = nc_in.variables[variable_in.upper()]
            else:
                print(f'[ERROR] Variable {variable_in} is not available in the file. Extract not created')
                newEXTRACT.close_file()
                return None
            var_array = ma.array(var_in[:])
            var_array = np.array(var_array.filled(-999.0))

            if variable_out not in newEXTRACT.EXTRACT.variables:
                variable = newEXTRACT.create_2D_variable_general(variable_out, var_array, limits)
            else:
                variable = newEXTRACT.EXTRACT.variables[variable_out]
                variable[0, :, :] = var_array[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]

            for at in var_in.ncattrs():
                if at == '_FillValue' or at == 'add_offset' or at == 'scale_factor':
                    continue
                variable.setncattr(at, var_in.getncattr(at))
            nc_in.close()

        return newEXTRACT
