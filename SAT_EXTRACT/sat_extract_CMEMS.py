import argparse,configparser,os, sys, pytz, subprocess, __init__,warnings
import pandas as pd
from netCDF4 import Dataset
from datetime import datetime as dt
from datetime import timedelta
import numpy.ma as ma
import numpy as np
import sat_extract as sextract
from SAT_EXTRACT.sat_extract import SatExtract

code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
import COMMON.common_functions as cfs
warnings.simplefilter("ignore",UserWarning)



class SatExtractCMEMS:

    def __init__(self,verbose):
        self.verbose = verbose
        self.sat_type = 'CMEMS'


    def is_l3_product(self):
        return True


    def get_extract_options(self,sat_extract_options):

        options_out = sat_extract_options.get_satellite_options(self.sat_type)
        if options_out is None:
            return None

        if options_out['dataset_name_file'].find('$UTM$')>0:
            options_out['is_utm'] = True

        n_bands = 0
        is_reflectance = False
        rrs_list = []
        if options_out['rrs_var_list'] is not None:
            for var in options_out['rrs_var_list']:
                try:
                    rrs_list.append(float(var))
                except:
                    print(f'[ERROR] rrs_var_list option should be a list of valid float values')
                    return None
            n_bands = len(options_out['rrs_var_list'])
            is_reflectance = True

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
                    print(f'[ERROR] Option otherfiles_var_list in CMEMS section should be a comma-separated list of namefile:varlist tokens. {str_val} is not correct')
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
                print(f'{cmems_time} is not a valid satellite time option. Skipping....')
                return None
        return satellite_time

    def get_files_day(self,datehere, input_path_info, sat_extract_options):

        extract_options = self.get_extract_options(sat_extract_options)

        path_source = input_path_info['path_source']
        org = input_path_info['org']
        path_day = sextract.get_path_date(path_source, org, datehere, False)
        if path_day is None:
            return [None]*2

        ncfiles, satellite_time = [None]*2

        name = cfs.get_name_for_date(extract_options['dataset_name_file'],extract_options['dataset_name_format_date'],datehere)


        if name.find('$UTM$')>0:
            preffix = name[0:name.index('$UTM$')]
            suffix = name[name.index('$UTM$')+5:]
            ncfiles = []
            satellite_time = []
            for name_day in os.listdir(path_day):
                if name_day.startswith(preffix) and name_day.endswith(suffix):
                    ncfiles.append(os.path.join(path_day,name_day))
                    stime = self.get_satellite_time(sat_extract_options,os.path.join(path_day,name_day),datehere)
                    satellite_time.append(stime)
                    if stime is None:
                        return [None]*2

        else:
            fname = os.path.join(path_day,name)
            if os.path.isfile(fname):
                ncfiles = [fname]
                satellite_time = [self.get_satellite_time(sat_extract_options,fname,datehere)]
                if satellite_time[0] is None:
                    return [None]*2
            else:
                print(f'[WARNING] {fname} is not a valid CMEMS file')


        return ncfiles, satellite_time

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
        list_files = extract_info['list_files']
        nfiles = len(list_files)

        if extract_options['is_utm']:
            dataset_name_file = cfs.get_name_for_date(extract_options['dataset_name_file'],extract_options['dataset_name_format_date'], extract_info['satellite_time'][0])
            iini = dataset_name_file.index('$UTM$')
            suffix  = dataset_name_file[iini+5:]
            sat_file_indices = [-1] * ntimes
            utm_files = ['']* nfiles
            for ifile in range(nfiles):
                name_file = os.path.basename(list_files[ifile])
                utm_zone = name_file[iini:name_file.index(suffix)]
                utm_files[ifile] = utm_zone
                for itime in range(ntimes):
                    utm_time = sextract.get_utm_zone(insitu_lat[itime],insitu_lon[itime])
                    if utm_time==utm_zone:
                        sat_file_indices[itime] = ifile
        else:
            sat_file_indices = [0] * ntimes
            utm_files = None



        sat_file_indices_used = np.unique(sat_file_indices)
        if np.max(sat_file_indices_used) == -1:
            print(
                f'[ERROR] No source files were found for the in situ data avaialable on {insitu_time[0].strftime("%Y-%m-%d")}')
            return
        site_list = []
        format_datetime = '%Y-%m-%dT%H:%M:%S'
        for idx in range(len(sat_file_indices_used)):
            ifile = sat_file_indices_used[idx]
            cmems_source = SatSourceCMEMS(list_files[ifile],extract_options)
            if not cmems_source.valid:
                continue
            if self.verbose:
                print(f'[INFO] Working with CMEMS source {list_files[ifile]}')
            indices_here = np.where(sat_file_indices==ifile)[0]
            ntimes_here = len(indices_here)
            insitu_lat_h = insitu_lat[indices_here]
            insitu_lon_h = insitu_lon[indices_here]
            insitu_time_h = insitu_time[indices_here]
            insitu_indices_h = insitu_indices[0][indices_here]
            valid_h = np.ones(ntimes_here)
            lat_array, lon_array = cmems_source.get_lat_lon_arrays()

            for itime in range(ntimes_here):
                global_at = extract_info['global_at'].copy()
                datehere_str = extract_info['satellite_time'][ifile].strftime('%Y%m%d')
                limits, rc = sextract.get_geo_info(extract_options['size_box'], insitu_lat[itime], insitu_lon[itime],lat_array, lon_array)
                if rc is None:
                    valid_h[itime]=0
                    continue

                site = f'{sextract.get_satellite_ref(global_at)}_{datehere_str}_{rc[0]}_{rc[1]}'
                ofname = os.path.join(output_path, f'extract_{site}.nc')

                if os.path.exists(ofname) and not overwrite:
                    few.write('\n')
                    few.write(
                        f'extract_{site}.nc;{insitu_indices_h[itime]};{insitu_time_h[itime].strftime(format_datetime)};{insitu_lat_h[itime]};{insitu_lon_h[itime]}')
                    print(
                        f'[WARNING] [{itime + 1}/{ntimes}] Satellite extract extract_{site}.nc already exists. {itime + 1}/{ntimes} Skiping...')
                    continue

                if overwrite and site in site_list:
                    few.write('\n')
                    few.write(
                        f'extract_{site}.nc;{insitu_indices_h[itime]};{insitu_time_h[itime].strftime(format_datetime)};{insitu_lat_h[itime]};{insitu_lon_h[itime]}')
                    print(
                        f'[WARNING] [{itime + 1}/{ntimes}] Satellite extract for {site}  has been already done.  Skiping...')
                    continue

                if self.verbose:
                    print(f'[INFO] Creating new extract: {ofname}')

                global_at_here = sextract.add_insitu_global_atrib(global_at, site, insitu_lat_h[itime],insitu_lon_h[itime], None)
                #window = limits
                extract_info_here = {
                    'global_at': global_at_here,
                    'limits': limits,
                    'size_box': extract_options['size_box'],
                    'lat_array': lat_array,
                    'lon_array': lon_array,
                    'satellite_time': [extract_info['satellite_time'][ifile]],
                    'n_bands': extract_options['n_bands']
                }

                ##start extract: define dimension and global attributes and creates satellite latitude, longitude and time variables
                newExtract = sextract.start_extract(ofname, extract_info_here, self.verbose)
                if newExtract is None:
                    print(f'[ERROR] satellite extract {ofname} could not be started')
                    continue
                ##wavelength and spectral variables (rrs,rrs_unc)
                if extract_options['n_bands']>0:
                    newExtract.create_satellite_bands_variable(extract_options['rrs_list'])
                    if not self.create_rrs_variables(newExtract,extract_info_here,extract_options, cmems_source):
                        print(f'[ERROR] Error creating Rrs variables. Deleting {ofname} extract and skipping...')
                        newExtract.close_file()
                        os.remove(ofname)
                        continue
                    if not self.create_spectral_variables(newExtract,extract_info_here,extract_options,cmems_source):
                        print(f'[ERROR] Error creating spectral variables. Deleting {ofname} extract and skipping...')
                        newExtract.close_file()
                        os.remove(ofname)
                        continue

                ##non spectral variables in the same file
                if extract_options['nospectral_var_list'] is not None:
                    if not self.create_non_spectral_variables(newExtract,extract_info_here,extract_options['nospectral_var_list'],cmems_source):
                        print(f'[ERROR] Error creating non- spectral variables. Deleting {ofname} extract and skipping...')
                        newExtract.close_file()
                        os.remove(ofname)
                        continue

                ##variables in other files
                if extract_options['otherfiles_var_list'] is not None:
                    ofvl = extract_options['otherfiles_var_list']
                    for name_of in ofvl:
                        name_here = cfs.get_name_for_date(name_of,extract_options['dataset_name_format_date'],extract_info['satellite_time'][ifile])
                        if extract_options['is_utm']:
                            name_here = name_here.replace('$UTM$',utm_files[ifile])
                        file_here = os.path.join(os.path.dirname(list_files[ifile]),name_here)
                        csource_here = SatSourceCMEMS(file_here,extract_options)
                        if not csource_here.valid:
                            print(
                                f'[ERROR] Error getting variables from the extra file: {name_here}')
                            newExtract.close_file()
                            os.remove(ofname)
                            continue
                        list_var = ofvl[name_of]
                        if len(list_var)==1 and list_var[0]=='*':
                            list_var = csource_here.get_list_nonspectral_variables()
                        list_var_end = []
                        for var in list_var:
                            var_n = var if var.startswith('satellite_') else f'satellite_{var}'
                            if var_n not in newExtract.EXTRACT.variables:
                                list_var_end.append(var)
                        self.create_non_spectral_variables(newExtract,extract_info_here,list_var_end,csource_here)

                newExtract.close_file()
                few.write('\n')
                few.write(
                    f'extract_{site}.nc;{insitu_indices_h[itime]};{insitu_time_h[itime].strftime(format_datetime)};{insitu_lat_h[itime]};{insitu_lon_h[itime]}')
        few.close()

    def create_spectral_variables(self,newExtract,extract_info,options,cmems_source):
        spectral_var_list = options['spectral_var_list']
        if spectral_var_list is None:
            print(f'[INFO]  No spectral variables are included. Skipping....')
            return True
        for var_format_name in spectral_var_list:
            var_name = var_format_name.replace('$BAND$','')
            if var_name.endswith('_'):
                var_name = var_name[:-1]
            var_name = var_name.replace('__','_')
            if not var_name.startswith('satellite_'):
                var_name = f'satellite_{var_name}'
            array,attrs = cmems_source.get_spectral_array(options['rrs_var_list'],var_format_name,extract_info['n_bands'], extract_info['limits'],False)
            if array is None:
                return False
            var_info = {
                'data_type':'f4',
                'fill_value': -999.0,
                'array': array,
                'attrs': attrs
            }
            var = newExtract.create_3D_variable_from_info_dict_impl(var_info,var_name)
            if var is None:
                return False
        return True


    def create_rrs_variables(self,newExtract,extract_info,options,cmems_source):

        array_rrs,attrs = cmems_source.get_spectral_array(options['rrs_var_list'],options['rrs_format'],extract_info['n_bands'],extract_info['limits'],True)
        if array_rrs is None:
            return False
        else:
            var_rrs= newExtract.create_rrs_variable(f'CMEMS_{extract_info["global_at"]["sensor"]}')
            var_rrs[0,:,:,:] = array_rrs[:,:,:]


        if len(options['rrs_unc_format'])==0 or options['rrs_unc_format'].lower()=='none':
            return True

        array_rrs_unc,attrs = cmems_source.get_spectral_array(options['rrs_var_list'], options['rrs_unc_format'],extract_info['n_bands'],extract_info['limits'],False)

        if array_rrs_unc is not None:
            var_rrs_unc = newExtract.create_rrs_unc_variable(f'CMEMS_{extract_info["global_at"]["sensor"]}')
            var_rrs_unc[0, :, :, :] = array_rrs_unc[:, :, :]


        return True

    def create_non_spectral_variables(self,newExtract,extract_info,var_list,cmems_source):

        for var_name in var_list:
            array,attrs = cmems_source.get_nonspectral_array(var_name,extract_info['limits'])
            if array is None:
                return False
            var_info = {
                'data_type':'f4',
                'fill_value': -999.0,
                'array': array,
                'attrs': attrs
            }
            if not var_name.startswith('satellite_'):
                var_name = f'satellite_{var_name}'
            var = newExtract.create_2D_variable_from_info_dict_impl(var_info,var_name)
            if var is None:
                return False


        return True

class SatSourceCMEMS:

    def __init__(self, path_source,extract_options):
        if os.path.isfile(path_source):
            self.path_source = path_source
        else:
            print(f'[ERROR] CMEMS source could not be  started with path: {path_source}')
            self.path_source = None

        self.dim_lat, self.dim_lon, self.dim_time = [None]*3
        self.valid = self.check_dataset(extract_options)

    def check_dataset(self,extract_options):
        if self.path_source is None:
            return False

        if extract_options is not None:
            dim_lat = extract_options['lat_dim']
            dim_lon = extract_options['lon_dim']
            dim_time = extract_options['time_dim']
        else:
            dim_lat = 'lat'
            dim_lon = 'lon'
            dim_time = 'time'

        try:
            dataset = Dataset(self.path_source)
        except:
            print(f'[ERROR]{self.path_source} is not a valid NetCDF file. It could not be opened')
            return False
        dv_available = True
        if not dim_lat in dataset.dimensions:
            print(f'[ERROR] Latitude dimension {dim_lat} is not available in the file {self.path_source}')
            dv_available = False
        if not dim_lat in dataset.variables:
            print(f'[ERROR] Latitude variable {dim_lat} is not available in the file {self.path_source}')
            dv_available = False
        if not dim_lon in dataset.dimensions:
            print(f'[ERROR] Longitude dimension {dim_lon} is not available in the file {self.path_source}')
            dv_available = False
        if not dim_lon in dataset.variables:
            print(f'[ERROR] Longitude variable {dim_lon} is not available in the file {self.path_source}')
            dv_available = False

        if len(dim_time)==0 or dim_time.lower()=='none':
            print(f'[WARNING] NetCDF file does not include a time dimension')
            dim_time = None
        else:
            notime = False
            if not dim_time in dataset.dimensions:
                print(f'[ERROR] Time dimension {dim_time} is not available in the file {self.path_source}')
                dv_available = False
                notime = True
            if not dim_time in dataset.variables:
                print(f'[ERROR] Time variable {dim_time} is not available in the file {self.path_source}')
                dv_available = False
                notime = True
            if notime:
                print(f'[WARNING] If time dimension and variable are not available in the source NetCDF files, please use time_dim: or time_dim: none in the [CMEMS] section in the configuration file')
        dataset.close()

        if dv_available:
            self.dim_lat,self.dim_lon,self.dim_time = dim_lat,dim_lon,dim_time


        return dv_available

    def get_lat_lon_arrays(self):
        if not self.valid:
            return [None]*2
        nc_sat = Dataset(self.path_source, 'r')
        lat = nc_sat.variables[self.dim_lat][:]
        lon = nc_sat.variables[self.dim_lon][:]
        nc_sat.close()
        return lat, lon




    def get_spectral_array(self,var_list,name_format,nbands,window,all_bands_required):
        if not self.valid:
            return [None]*2
        ##checking number of bands
        if 0 < nbands != len(var_list):
            print(f'[ERROR] Inconsistency between the number of bands {nbands} and number of variables ({len(nbands)}) defined in the configuration file')
            return [None]*2
        if name_format.find('$BAND$')==-1:
            print(f'[ERROR] Format for spectral variables should contain the tag $BAND$. {name_format} is not valid')
            return None
        ny = window[1]-window[0]
        nx = window[3]-window[2]
        array = ma.masked_all((nbands,ny,nx))
        attrs = None
        data_available = False
        nc_sat = Dataset(self.path_source, 'r')
        for iband in range(nbands):
            band = var_list[iband].replace('.','_')
            name_var = name_format.replace('$BAND$',band)
            if name_var not in nc_sat.variables:
                if all_bands_required:
                    print(f'[ERROR] Variable {name_var} is not available in the dataset, spectral array could not be retrieved. ')
                    nc_sat.close()
                    return [None]*2
                else:
                    print(f'[WARNING] Variable {name_var} is not available in the dataset, spectral array for band {iband}: {band} is masked')
            else:
                if self.dim_time is None:
                    array[iband, :, :] = nc_sat.variables[name_var][window[0]:window[1], window[2]:window[3]]
                else:
                    array[iband,:,:] = nc_sat.variables[name_var][0,window[0]:window[1],window[2]:window[3]]
                if attrs is None:
                    attrs = cfs.get_attrs_variable(nc_sat.variables[name_var])
                data_available  = True

        nc_sat.close()

        if not data_available:
            print(f'[WARNING] No spectral data found for {name_format}')
            return [None]*2

        return array,attrs

    def get_nonspectral_array(self,name_var,window):
        if not self.valid:
            return [None]*2

        nc_sat = Dataset(self.path_source, 'r')
        if name_var not in nc_sat.variables:
            print(f'[ERROR] Variable {name_var} is not available in the dataset, spectral array could not be retrieved. ')
            nc_sat.close()
            return [None]*2
        else:
            if self.dim_time is None:
                array = nc_sat.variables[name_var][window[0]:window[1], window[2]:window[3]]
            else:
                array = nc_sat.variables[name_var][0,window[0]:window[1],window[2]:window[3]]
            attrs = cfs.get_attrs_variable(nc_sat.variables[name_var])

        nc_sat.close()

        return array,attrs

    def get_list_nonspectral_variables(self):
        if not self.valid:
            return None
        var_list = []
        nc_sat = Dataset(self.path_source, 'r')
        for name_var in nc_sat.variables:
            dims = nc_sat.variables[name_var].get_dims()
            if len(dims)==3 and dims[0].name==self.dim_time and dims[1].name==self.dim_lat and dims[2].name==self.dim_lon:
                var_list.append(name_var)
        nc_sat.close()
        return var_list



