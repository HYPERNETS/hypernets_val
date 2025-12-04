import os, sys, argparse, __init__

import pytz
import numpy as np

code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)

import COMMON.common_functions as cfs
import sat_extract as sextract
from netCDF4 import Dataset
from datetime import datetime as dt
from datetime import timedelta
# import user defined functions from other .py


os.environ['QT_QPA_PLATFORM'] = 'offscreen'  # to avoid error "QXcbConnection: Could not connect to display"
path2ncrcat = '/opt/local/bin/ncrcat'

# parser = argparse.ArgumentParser(description="Create sat extract from PACE-OCI AOP files.")
# parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
# parser.add_argument('-sd', "--startdate", help="The Start Date - format YYYY-MM-DD ")
# parser.add_argument('-ed', "--enddate", help="The End Date - format YYYY-MM-DD ")
# parser.add_argument('-site', "--sitename", help="Site name.", choices=['VEIT', 'BEFR', 'BSBE','GAIT','LAIT','TRIT'])
# parser.add_argument('-date_file',"--date_list_file",help="Date list file")
# parser.add_argument('-c', "--config_file", help="Config File.")
# parser.add_argument('-ps', "--path_to_sat", help="Path to satellite sources.")
# parser.add_argument('-wce', "--wce_expression", help="Wild card expression.")
# parser.add_argument('-o', "--output", help="Path to output")
# # parser.add_argument('-nl', "--nolist",help="Do not create initial satellite lists, checking day by day allowing download",action="store_true")
# # parser.add_argument('-adownload', "--allow_download", help="Allow download", action="store_true")
#
# args = parser.parse_args()


class SatSourceOCI():

    def __init__(self, path_source):
        self.path_source = path_source
        self.nrow = -1
        self.ncol = -1
        self.nbands = -1
        self.wavelength = None
        self.dimensions = ['number_of_lines', 'pixels_per_line', 'bands_per_pixel',
                           'number_of_reflectance_location_values', 'number_of_cloud_phases', 'number_of_bands',
                           'number_of_reflective_bands', 'wavelength_3d']
        self.file_check = {
            'sensor_band_parameters': ['wavelength_3d'],
            'scan_line_attributes': ['msec', 'csol_z'],
            'geophysical_data': ['Rrs', 'Rrs_unc', 'l2_flags', 'avw', 'aot_865', 'angstrom'],
            'navigation_data': ['latitude', 'longitude', 'tilt'],
            'global_attributes': ['time_coverage_start', 'time_coverage_end']
        }
        self.start_time = None
        self.end_time = None
        self.VALID = self.check_file()
        if not self.VALID:
            print(f'[ERROR] {self.path_source} is not a valid PACE file')

        self.file_geo = self.get_file_geo()

    def get_file_geo(self):
        name_file = os.path.basename(self.path_source)
        name_file_geo = name_file.replace('OC_AOP','OC_GEO')
        file_geo = os.path.join(os.path.dirname(self.path_source),name_file_geo)
        return file_geo if os.path.exists(file_geo) else None

    def set_file_geo_from_options(self,work_date,name_format,name_format_date_file):
        if work_date is None:
            name_file = os.path.basename(self.path_source)
            work_date_str = name_file.split('.')[1]
            try:
                work_date = dt.strptime(work_date_str,'%Y%m%dT%H%M%S')
            except:
                return
            name_file_geo = name_format.replace('$DATE$', work_date.strftime(name_format_date_file))
        else:
            name_file_geo = name_format.replace('$DATE$',work_date.strftime(name_format_date_file))
        file_geo = os.path.join(os.path.dirname(self.path_source),name_file_geo)
        self.file_geo = file_geo if os.path.exists(file_geo) else None

    def get_lat_lon_oza_arrays(self):
        lat_array,lon_array =self.get_lat_long_arrays()
        oza_array = None
        if self.file_geo is None:
            return lat_array, lon_array, oza_array

        try:
            dataset = Dataset(self.file_geo)
        except:
            print(f'[WARNING] {self.file_geo} is not a valid NetCDF dataset. OC_GEO file can not be used')
            return lat_array, lon_array,oza_array
        if 'sensor_zenith' in dataset.variables:
            print(f'[INFO] Retrieving sensor zenith from OC_GEO file')
            oza_array = dataset.variables['sensor_zenith'][:]
        else:
            print(f'[WARNING] geolocation_data is not a valid group in NetCDF dataset {self.path_source}')

        dataset.close()

        return lat_array, lon_array, oza_array

    def check_file(self):
        if not os.path.isfile(self.path_source):
            return False

        try:
            dataset = Dataset(self.path_source)
            for dim in self.dimensions:
                if dim not in dataset.dimensions:
                    dataset.close()
                    print(f'[WARNING] Dimension {dim} is not included in the dataset')
                    return False
                if dim=='number_of_lines':
                    self.nrow = dataset.dimensions[dim].size
                if dim=='pixels_per_line':
                    self.ncol = dataset.dimensions[dim].size
                if dim=='wavelength_3d':
                    self.nbands = dataset.dimensions[dim].size
            for group in self.file_check:
                for value in self.file_check[group]:
                    if group == 'global_attributes':
                        if value not in dataset.ncattrs():
                            dataset.close()
                            print(f'[WARNING] Attribute {value} is not included in the global attribute list')
                            return False
                    else:
                        if group not in dataset.groups:
                            dataset.close()
                            print(f'[WARNING] Group {group} is not included in the dataset')
                            return False
                        if value not in dataset[group].variables:
                            dataset.close()
                            print(f'[WARNING] Variable {value} is not included in the group {group} of the dataset')
                            return False
            sbp = dataset['sensor_band_parameters']
            self.wavelength = sbp.variables['wavelength_3d'][:]
            self.start_time = dt.strptime(dataset.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
            self.end_time = dt.strptime(dataset.time_coverage_end, '%Y-%m-%dT%H:%M:%S.%fZ')
            print(f'[INFO] Starting OCI dataset: {os.path.basename(self.path_source)}')
            print(f'[INFO]  -> Number of lines (rows): {self.nrow}')
            print(f'[INFO]  -> Pixels per line (cols): {self.ncol}')
            print(f'[INFO]  -> Number of bands {self.nbands} ({self.wavelength[0]} nm - {self.wavelength[-1]} nm)')
            print(f'[INFO]  -> Overpass time: {self.start_time.strftime("%Y-%m-%dT%H:%M:%S")} - {self.end_time.strftime("%Y-%m-%dT%H:%M:%S")}')

            dataset.close()
        except:
            return False

        return True

    def get_sat_time(self, row):
        if not self.VALID:
            return None
        dataset = Dataset(self.path_source)
        ref_time = self.start_time
        if row >= 0:
            ref_time = ref_time.replace(hour=0, minute=0, second=0, microsecond=0, tzinfo=pytz.UTC)
            seconds_day = float(dataset['scan_line_attributes'].variables['msec'][row]) / 1000
            ref_time = ref_time + timedelta(seconds=seconds_day)
        else:
            ref_time = ref_time.replace(tzinfo=pytz.utc)
        dataset.close()
        return ref_time

    def check_location(self, in_situ_lat, in_situ_lon, size_box):

        lat, lon = self.get_lat_long_arrays()

        return self.check_location_impl(in_situ_lat, in_situ_lon,lat,lon,size_box)

    def check_location_impl(self, in_situ_lat, in_situ_lon, lat, lon, size_box):

        #lat, lon = self.get_lat_long_arrays()
        contain_flag = cfs.contain_location(lat, lon, in_situ_lat, in_situ_lon)
        r = -1
        c = -1
        if contain_flag == 1:
            r, c = cfs.find_row_column_from_lat_lon(lat, lon, in_situ_lat, in_situ_lon)

            start_idx_y, stop_idx_y, start_idx_x, stop_idx_x = self.get_window_limits(r, c, size_box)

            if start_idx_y >= 0 and (stop_idx_y + 1) < lat.shape[0] and start_idx_x >= 0 and (stop_idx_x + 1) < \
                    lat.shape[
                        1]:
                pass
            else:
                r = -1
                c = -1

        return r, c

    def get_lat_long_arrays(self):
        lat_array = None
        lon_array = None
        try:
            dataset = Dataset(self.path_source)
        except:
            print(f'[ERROR] {self.path_source} is not a valid NetCDF dataset')
            return lat_array, lon_array
        try:
            ng = dataset['navigation_data']
            lat_array = ng.variables['latitude'][:]
            lon_array = ng.variables['longitude'][:]
        except:
            print(f'[ERROR] navigation_Date is not a valid group in NetCDF dataset {self.path_source}')

        dataset.close()

        return lat_array, lon_array

    def get_window_limits(self, r, c, size_box):
        start_idx_y = (r - int(size_box / 2))
        stop_idx_y = (r + int(size_box / 2) + 1)
        start_idx_x = (c - int(size_box / 2))
        stop_idx_x = (c + int(size_box / 2) + 1)

        return start_idx_y, stop_idx_y, start_idx_x, stop_idx_x


    def get_2D_subarray(self,group,variable,window):
        if not self.VALID:
            return None
        #window = [start_idx_y, stop_idx_y, start_idx_x, stop_idx_x]
        try:
            dataset = Dataset(self.path_source)
            var_array = dataset[group].variables[variable][window[0]:window[1],window[2]:window[3]]
            attrs = dataset[group].variables[variable].__dict__
            dataset.close()
        except Exception as ex:
            print(f'[ERROR] {group}/{variable} could not be found in {self.path_source}. Skipping...')
            return None,None
        return var_array,attrs

    def get_3D_subarray(self,group,variable,window):
        if not self.VALID:
            return None
        dataset = Dataset(self.path_source)
        var_array = dataset[group].variables[variable][window[0]:window[1], window[2]:window[3],:]
        attrs = dataset[group].variables[variable].__dict__
        dataset.close()
        return var_array, attrs

    def get2D_subarray_from_scan_line_variable(self,group,variable,window,size_box):
        if not self.VALID:
            return None
        import numpy as np
        dataset = Dataset(self.path_source)
        var_array_1D = dataset[group].variables[variable][window[0]:window[1]]
        var_array = np.ma.repeat(var_array_1D.reshape(size_box,1),size_box,axis=1)
        attrs = dataset[group].variables[variable].__dict__
        dataset.close()
        return var_array, attrs

    def get_geo_variable(self,var_name,window):
        array = None
        attrs = None
        if self.file_geo is None:
            return array,attrs

        try:
            dataset = Dataset(self.file_geo)
        except:
            print(f'[WARNING] {self.file_geo} is not a valid NetCDF dataset. Geometry variables can not be retrieved')
            return array,attrs

        if var_name in dataset.variables:
            all_array = dataset.variables[var_name][:]
            if window is None:
                array = all_array
            else:
                array = all_array[window[0]:window[1], window[2]:window[3]]
            attrs = dataset.variables[var_name].__dict__
        else:
            print(f'[WARNING] variable {var_name} is not available in the dataset {self.path_source}')


        dataset.close()

        return array,attrs

    def get_tilt_value(self,window):
        if not self.VALID:
            return None
        row = 0 if window is None else window[0]
        dataset = Dataset(self.path_source)
        value = dataset['navigation_data'].variables['tilt'][row]
        dataset.close()
        return value

    def get_global_attrs(self,site_name):
        if not self.VALID:
            return None
        from datetime import datetime as dt
        dataset = Dataset(self.path_source)
        if site_name is None:
            site_name = 'UNKNOWN'
        at = {
            'creation_time': dt.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
            'satellite': dataset.platform,
            'platform': '',
            'sensor': dataset.instrument,
            'res':  dataset.spatialResolution,
            'aco_processor':dataset['processing_control'].software_name,
            'proc_version':dataset['processing_control'].software_version,
            'site':site_name
        }
        return at


class SatExtractOCI:

    def __init__(self, verbose):
        self.verbose = verbose
        self.sat_type = 'OCI'

    def get_extract_options(self,sat_extract_options):
        options_out = sat_extract_options.get_satellite_options('OCI')
        if options_out is None:
            return None
        options_out['size_box'] = sat_extract_options.get_box_size()
        try:
            omv = options_out['other_main_variables']
            omv_dict ={}
            for group_vars in omv.split(';'):
                gv = group_vars.split(':')
                omv_dict[gv[0].strip()] = [x.strip() for x in gv[1].split(',')]
            if len(omv_dict):
                options_out['other_main_variables'] = omv_dict
            else:
                options_out['other_main_variables'] = None

        except Exception as ex:
            print(f'[WARNING] Option other_main_variables could not be parsed correctly {ex} Variables will not be used')
            options_out['other_main_variables'] = None
        return options_out

    ##Method to retrieve file for a single date. Deprecated: get_cmems_multiple_product_day
    def get_files_day(self,datehere, input_path_info, sat_extract_options):

        extract_options = self.get_extract_options(sat_extract_options)

        path_source = input_path_info['path_source']
        org = input_path_info['org']
        path_day = sextract.get_path_date(path_source, org, datehere, False)
        if path_day is None:
            return [None]*2

        timeliness = extract_options['timeliness']
        if timeliness is None:
            print(f'[ERROR] {timeliness} options is not valid. Choose among NRT, DT or ANY')
            return [None]*2
        version_r = extract_options['version']

        main_file_tag = extract_options['main_file_tag']##Default for RRS: OC_AOP


        list_files = {}
        for name in os.listdir(path_day):
            if not name.startswith('PACE_OCI'):
                continue
            if not name.find(main_file_tag)>0:
                continue
            name_s = name.split('.')
            datetime =  name_s[1]
            version = name_s[4]
            if version!=version_r:
                continue
            tn = 'NRT' if name_s[-2] == 'NRT' else 'DT'
            if datetime not in list_files.keys():
                list_files[datetime]= {
                    'NRT': None,
                    'DT': None
                }

            list_files[datetime][tn] = name

        satellite_time = []
        ncfiles = []
        for l in list_files:
            satellite_time_here = dt.strptime(l,'%Y%m%dT%H%M%S')
            nrt_name = list_files[l]['NRT']
            dt_name = list_files[l]['DT']
            if nrt_name is not None and (timeliness=='NRT' or (timeliness=='ANY' and dt_name is None)):
                satellite_time.append(satellite_time_here)
                ncfiles.append(os.path.join(path_day,nrt_name))
            elif dt_name is not None and (timeliness=='DT' or timeliness=='ANY'):
                satellite_time.append(satellite_time_here)
                ncfiles.append(os.path.join(path_day,dt_name))


        #print(ncfiles,satellite_time)
        if len(ncfiles)>0:
            return ncfiles,satellite_time
        else:
            return [None]*2

    def is_l3_product(self):
        return False

    def run_parallel_extract_day(self,params):
        self.run_extract_day(params[0], params[1], params[2], params[3], params[4], params[5])

    ##lat_array and lon_array are kept as parameteres for compability with other extracts
    def run_extract_day(self,extract_options, extract_info, lat_array, lon_array, output_path, overwrite):
        insitu_lat = extract_info['insitu_lat']
        insitu_lon = extract_info['insitu_lon']
        insitu_time = extract_info['insitu_time']
        insitu_indices = extract_info['insitu_indices']

        path_extract_output = extract_info['path_extract_output']
        few = open(path_extract_output, 'w')
        few.write('extract_file;insitu_index;insitu_time;insitu_lat;insitu_lon')

        if self.verbose:
            print(f'[INFO] Output path extract list: {path_extract_output}')

        ntimes = len(insitu_time)
        list_files = extract_info['list_files']
        sat_file_indices = [-1]*ntimes
        geo_info_array = [None]*ntimes

        for ifile in range(len(list_files)):
            if self.verbose:
                print(f'[INFO] Checking file...')
            oci_source = SatSourceOCI(list_files[ifile])
            if oci_source.file_geo is None:
                oci_source.set_file_geo_from_options(None,extract_options['geo_file'],extract_options['geo_file_date_format'])
            if self.verbose:
                print(f'[INFO] OCI source was loaded')
                print(f'[INFO] Geo info file: {oci_source.file_geo}')
            lat_array, lon_array, oza_array = oci_source.get_lat_lon_oza_arrays()
            if self.verbose:
                print(f'[INFO] Latitude, longitude and oza arrays were obtained. Checking {ntimes} data points')
                print(f'[INFO] Latitude range: {np.min(lat_array)} to {np.max(lat_array)}')
                print(f'[INFO] Longitude range: {np.min(lon_array)} to {np.max(lon_array)}')


            for itime in range(ntimes):
                # print(f'[INFO] Checking point... {itime}')
                limits,rc = sextract.get_geo_info(extract_options['size_box'], insitu_lat[itime], insitu_lon[itime],lat_array, lon_array)


                oza = None if rc is None else abs(oza_array[rc[0],rc[1]])

                if ifile==0:
                    geo_info_array[itime] = {'rc':rc,'limits':limits,'oza':oza}
                    if rc is not None:
                        sat_file_indices[itime] = ifile
                else:
                    if rc is not None:
                        geo_info_prev = geo_info_array[itime]
                        update_geo_info = False
                        if geo_info_prev['rc'] is None:
                            update_geo_info = True
                        elif geo_info_prev['oza'] is not None and oza is not None and oza<geo_info_prev['oza']:
                            print(f'[WARNING] More than one image is available for point: {insitu_time[itime]},{insitu_lat[itime]},{insitu_lon[itime]}. Choosing {list_files[ifile]} as oza {oza} is lower than {geo_info_prev["oza"]} in file {list_files[sat_file_indices[itime]]}')
                            update_geo_info = True
                        if update_geo_info:
                            geo_info_array[itime] = {'rc': rc, 'limits': limits, 'oza': oza}
                            sat_file_indices[itime]= ifile


        if len(sat_file_indices)==1:
            sat_file_indices_used = sat_file_indices.copy()
        else:
            sat_file_indices_used = np.unique(sat_file_indices[sat_file_indices!=-1])

        if np.max(sat_file_indices_used)==-1:
            print(f'[ERROR] No source granules were found for the in situ data avaialable on {insitu_time[0].strftime("%Y-%m-%d")}')
            return
        geo_info_array = np.array(geo_info_array)
        site_list = []
        format_datetime = '%Y-%m-%dT%H:%M:%S'
        print('---------------------------------------------->',sat_file_indices_used)
        for idx in range(len(sat_file_indices_used)):
            ifile = sat_file_indices_used[idx]
            oci_source = SatSourceOCI(list_files[ifile])
            if oci_source.file_geo is None:
                oci_source.set_file_geo_from_options(None,extract_options['geo_file'],extract_options['geo_file_date_format'])
            if self.verbose:
                print(f'[INFO] Working with OCI source {list_files[ifile]}')
            print('tal',sat_file_indices,sat_file_indices_used,ifile)
            indices_here = np.where(sat_file_indices==ifile)[0]
            print('-->',indices_here)
            ntimes_here = len(indices_here)
            insitu_lat_h = insitu_lat[indices_here]
            insitu_lon_h = insitu_lon[indices_here]
            insitu_time_h = insitu_time[indices_here]
            insitu_indices_h = insitu_indices[0][indices_here]
            geo_info_array_here = geo_info_array[indices_here]
            lat_array, lon_array = oci_source.get_lat_long_arrays()

            if self.verbose:
                print(f'[INFO] Number of in situ spectra: {ntimes_here}')

            for itime in range(ntimes_here):
                global_at = extract_info['global_at'].copy()
                datehere_str = extract_info['satellite_time'][0].strftime('%Y%m%d')
                rc = geo_info_array_here[itime]['rc']
                site = f'{sextract.get_satellite_ref(global_at)}_{datehere_str}_{rc[0]}_{rc[1]}'
                ofname = os.path.join(output_path, f'extract_{site}.nc')

                if os.path.exists(ofname) and not overwrite:
                    few.write('\n')
                    few.write(
                        f'extract_{site}.nc;{insitu_indices_h[itime]};{insitu_time_h[itime].strftime(format_datetime)};{insitu_lat_h[itime]};{insitu_lon_h[itime]}')
                    print(
                        f'[INFO] [{itime + 1}/{ntimes}] Satellite extract extract_{site}.nc already exists. {itime + 1}/{ntimes} Skipping...')
                    continue

                if overwrite and site in site_list:
                    few.write('\n')
                    few.write(
                        f'extract_{site}.nc;{insitu_indices_h[itime]};{insitu_time_h[itime].strftime(format_datetime)};{insitu_lat_h[itime]};{insitu_lon_h[itime]}')
                    print(
                        f'[INFO] [{itime + 1}/{ntimes}] Satellite extract for {site}  has been already done.  SkiPping...')
                    continue

                if self.verbose:
                    print(f'[INFO] Creating new extract: {ofname}')

                global_at_here = sextract.add_insitu_global_atrib(global_at, site, insitu_lat_h[itime],
                                                                     insitu_lon_h[itime], None)
                window = geo_info_array_here[itime]['limits']
                extract_info_here = {
                    'global_at': global_at_here,
                    'limits': geo_info_array_here[itime]['limits'],
                    'size_box': extract_options['size_box'],
                    'lat_array': lat_array,
                    'lon_array': lon_array,
                    'satellite_time': [extract_info['satellite_time'][ifile]],
                    'n_bands': oci_source.nbands
                }
                ##start extract: define dimension and global attributes and creates satellite latitude, longitude and time variables
                newExtract = sextract.start_extract(ofname, extract_info_here,self.verbose)
                if newExtract is None:
                    print(f'[ERROR] satellite extract {ofname} could not be started')
                    continue
                ##wavelength
                newExtract.create_satellite_bands_variable(oci_source.wavelength)
                ##rrs variables
                if extract_options['main_file_tag']=='OC_AOP':
                    if self.verbose:
                        print(f'[INFO] Adding Rrs variables...')
                    rrs_t = oci_source.get_3D_subarray('geophysical_data', 'Rrs', window)
                    rrs_unc_t = oci_source.get_3D_subarray('geophysical_data', 'Rrs_unc', window)
                    brrs = self.create_rrs_oci_variables(newExtract,rrs_t, rrs_unc_t)
                    if not brrs:
                        print(f'[ERROR] rrs variables could not be added to satellite extract {ofname}')
                        newExtract.close_file()
                        os.remove(ofname)
                        continue
                ##flag variable
                if self.verbose:
                    print(f'[INFO] Adding flag variable...')
                array, attrs = oci_source.get_2D_subarray('geophysical_data', 'l2_flags', window)
                self.create_l2_oci_flag_variable(newExtract,array, attrs)
                ##geometry variables
                if self.verbose:
                    print(f'[INFO] Adding geometry variables...')
                self.create_geometry_variables(oci_source,newExtract,window)

                ##other variables
                other_main_variables = extract_options['other_main_variables']
                if other_main_variables is not None:
                    for group in other_main_variables:
                        for var_name in other_main_variables[group]:
                            array, attrs = oci_source.get_2D_subarray(group, var_name, window)
                            if array is not None:
                                self.create_2D_oci_variable(newExtract,var_name, array, attrs)

                #array, attrs = oci_source.get2D_subarray_from_scan_line_variable('scan_line_attributes', 'csol_z',window, size_box)
                #newExtract.create_2D_oci_variable('csol_z', array, attrs)
                #newExtract.create_tilt_variable(oci_source.get_tilt_value(window))
                newExtract.close_file()
                # 'n_bands': extract_options['n_bands'],
                # 'list_files': extract_info['list_files']
                few.write('\n')
                few.write(
                    f'extract_{site}.nc;{insitu_indices_h[itime]};{insitu_time_h[itime].strftime(format_datetime)};{insitu_lat_h[itime]};{insitu_lon_h[itime]}')




        few.close()

        # if len(list_files)==1:
        #     sat_file_indices = np.zeros(ntimes)
        # else:
        #     sat_file_indices = []##I have to deduce the file for each point

        # for ifile in range(len(list_files)):
        #     indices_ifile = np.where(sat_file_indices==ifile)
        #     if len(indices_ifile[0])==0:
        #         continue
        #     oci_source = SatSourceOCI(list_files[ifile])
        #     lat_array,lon_array, oza_array = oci_source.get_lat_lon_oza_arrays()
        #     insitu_time_here = insitu_time[indices_ifile]
        #     insitu_lat_here = insitu_lat[indices_ifile]
        #     insitu_lon_here = insitu_lon[indices_ifile]

    def create_geometry_variables(self,oci_source,newExtract,window):
        geom_variables = ['satellite_SZA','satellite_SAA','satellite_OZA','satellite_OAA']
        orig_variables = ['solar_zenith','solar_azimuth','sensor_zenith','sensor_azimuth']
        for var_geo,var_orig in zip(geom_variables,orig_variables):
            array, attrs = oci_source.get_geo_variable(var_orig,window)
            if array is not None:
                self.create_2D_oci_variable(newExtract,var_geo,array,attrs)

    def create_rrs_oci_variables(self,newExtract,rrs_t,rrs_unc_t):
        newExtract.create_rrs_variable('OCI')
        newExtract.create_rrs_unc_variable('OCI')
        b1 = self.set_rrs_data(newExtract,'satellite_Rrs',rrs_t[0],rrs_t[1])
        b2 = self.set_rrs_data(newExtract,'satellite_Rrs_unc', rrs_unc_t[0], rrs_unc_t[1])
        if b1 and b2:
            return True
        else:
            return False

    def set_rrs_data(self,newExtract,var_name,var_array,attrs):
        nbands = newExtract.EXTRACT.dimensions['satellite_bands'].size
        if nbands==var_array.shape[2]:
            for iband in range(nbands):
                newExtract.EXTRACT.variables[var_name][0,iband,:,:] = var_array[:,:,iband]

            for at in attrs:##attributes, long_name is already set
                if at == '_FillValue' or at == 'scale_factor' or at == 'add_offset' or at=='long_name':
                    continue
                val = attrs[at]
                if (at == 'valid_min' or at == 'valid_max') and 'scale_factor' in attrs and 'add_offset' in attrs:
                    val = (attrs[at] * attrs['scale_factor']) + attrs['add_offset']
                newExtract.EXTRACT.variables[var_name].setncattr(at, val)
            return True
        else:
            print(f'[ERROR] {var_name} data could not be set, number of bands do not coincide with the satellelit_bands dimensions')
            return False

    def create_l2_oci_flag_variable(self,newExtract,array,attrs):
        satellite_flag = newExtract.EXTRACT.createVariable('satellite_l2_flags','i4',('satellite_id', 'rows', 'columns'),zlib = True, complevel=6)

        satellite_flag[0, :, :] = array[:,:]
        satellite_flag.setncatts(attrs)

    def create_2D_oci_variable(self,newExtract,name,array,attrs):
        if not name.startswith('satellite_'): name = f'satellite_{name}'
        satellite_variable = newExtract.EXTRACT.createVariable(name,'f4',('satellite_id', 'rows', 'columns'),fill_value=-999,zlib = True, complevel=6)
        satellite_variable[0,:,:] = array[:,:]
        for at in attrs:
            if at=='_FillValue' or  at=='scale_factor' or at=='add_offset':
                continue
            val = attrs[at]
            if (at=='valid_min' or at=='valid_max') and 'scale_factor' in attrs and 'add_offset' in attrs:
                val = (attrs[at] * attrs['scale_factor']) + attrs['add_offset']


            satellite_variable.setncattr(at,val)
    # def create_tilt_variable(self,tilt_value):
    #     satellite_tilt = self.EXTRACT.createVariable('satellite_tilt', 'f4', ('satellite_id'), fill_value=-999,
    #                                                  zlib=True, complevel=6)
    #     # print('Satellite start time es: ',satellite_start_time)
    #     satellite_tilt[0] = float(tilt_value)
    #     satellite_tilt.long_name = 'Sensor tilt angle'
    #     satellite_tilt.units = 'degrees'
    #     satellite_tilt.valid_min = -25.0
    #     satellite_tilt.valid_max = 25.0




def launch_create_extract(in_situ_site, path_product, basic_options):
    # try:
    basic_options['in_situ_site'] = in_situ_site
    if 'variable_list' not in basic_options:
        basic_options['variable_list'] = []
    path_output = in_situ_site['path_out']
    if not os.path.exists(path_output):
        os.mkdir(path_output)
    create_extractv2(path_product, path_output, basic_options)
    if args.verbose:
        print('----------------------------------------------')
        # extract_path = create_extract(size_box, site, path_source, path_output, in_situ_lat, in_situ_lon, res_str,
        #                               make_brdf, None)
        # if not extract_path is None:
        #     print(f'file created: {extract_path}')
    # except Exception as e:
    #     if args.verbose:
    #         print(f'Exception launching extract: {e}')
    #         pass


def create_extractv2(path_product, path_output, options):
    size_box = options['size_box']
    # res_str = options['resolution']
    # make_brdf = options['make_brdf']
    insitu_info = options['in_situ_site']
    site_name = insitu_info['site_name']
    in_situ_lat = insitu_info['latitude']
    in_situ_lon = insitu_info['longitude']
    variable_list = options['variable_list']
    if args.verbose:
        print(f'[INFO] Creating extract v.2 for {site_name} from {path_product}')

    if args.verbose:
        print(f'[INFO] Checking in situ  location -> lat: {in_situ_lat}; lon: {in_situ_lon}')

    oci_source = SatSourceOCI(path_product)
    r, c = oci_source.check_location(in_situ_lat, in_situ_lon, size_box)
    if r == -1 and c == -1:
        print('[WARNING] File does NOT contains the in situ location! Skipping...')
        return

    if site_name == 'SHIPBORNE':
        site_name = f'{r}_{c}'
        options['site_name'] = site_name

    filename = 'extract_' + os.path.basename(path_product)[:-2].replace('.', '_') + f'{site_name}' + '.nc'
    ofname = os.path.join(path_output, filename)
    if os.path.exists(ofname):
        print(f'[WARNING] Sat. extract file: {ofname} already exist. Skipping...')
        return

    start_idx_y, stop_idx_y, start_idx_x, stop_idx_x = oci_source.get_window_limits(r, c, size_box)
    window = [start_idx_y, stop_idx_y, start_idx_x, stop_idx_x]
    if args.verbose:
        print(f'[INFO] Getting extraction window: {stop_idx_y}-{stop_idx_y}:{start_idx_x}-{stop_idx_x}')

    if args.verbose:
        print(f'[INFO] Starting extract file: {ofname}')

    newExtract = SatExtractOCI(ofname, variable_list)

    if args.verbose:
        print(f'[INFO] Number of satellite bands: {oci_source.nbands}')
    newExtract.create_dimensions(size_box, oci_source.nbands)


    newExtract.set_global_attributes(oci_source.get_global_attrs(site_name))
    lat, lon = oci_source.get_lat_long_arrays()
    newExtract.create_lat_long_variables(lat, lon, window)
    newExtract.create_satellite_time_variable(oci_source.get_sat_time(r))
    newExtract.create_satellite_bands_variable(oci_source.wavelength)
    ##rrs variables
    rrs_t = oci_source.get_3D_subarray('geophysical_data','Rrs',window)
    rrs_unc_t = oci_source.get_3D_subarray('geophysical_data','Rrs_unc',window)
    brrs  = newExtract.create_rrs_oci_variables(rrs_t,rrs_unc_t)
    ##flag variable
    array,attrs = oci_source.get_2D_subarray('geophysical_data','l2_flags',window)
    newExtract.create_l2_oci_flag_variable(array,attrs)
    ##other variables
    for group in newExtract.variable_list:
        for var_name in newExtract.variable_list[group]:
            array,attrs = oci_source.get_2D_subarray(group,var_name,window)
            newExtract.create_2D_oci_variable(var_name,array,attrs)
    array, attrs = oci_source.get2D_subarray_from_scan_line_variable('scan_line_attributes', 'csol_z', window, size_box)
    newExtract.create_2D_oci_variable('csol_z', array, attrs)
    newExtract.create_tilt_variable(oci_source.get_tilt_value(window))


    #newExtract.cre
    # newExtract.create_geometry_variables()
    # b = newExtract.set_rrs_data(oci_source, window)
    # newExtract.set_geometry_data(path_product, size_box, window)
    newExtract.close_file()

    if not brrs:
        os.remove(ofname)
        return None

    return ofname


# # %%
# def main():
#     print('Creating satellite extracts.')
#
#     # Getting basic_options from config_file or arguments
#     basic_options = None
#     options = None
#     if args.config_file:
#         if os.path.isfile(args.config_file):
#             options = sextract.config_reader(args.config_file)
#             basic_options = sextract.get_basic_options_from_file_config(args, options)
#     if basic_options is None:
#         basic_options = sextract.get_basic_options_from_arguments(args)
#     if basic_options is None:
#         return
#     if args.verbose:
#         for option in basic_options:
#             print(f'[INFO] {option}:{basic_options[option]}')
#
#     wce = basic_options['wce']
#     if wce is None:
#         wce = f'"PACE_OCI*.L2.OC_AOP.V3_0.nc"'
#         #wce = f'"PACE_OCI*.L2.OC_AOP.V2_0.NRT.nc"'  # wild card expression
#
#
#     in_situ_site = sextract.get_insitu_site(args, options, basic_options['path_out'])
#     if args.verbose:
#         print(
#             f'[INFO] In situ site: {in_situ_site["site_name"]} {in_situ_site["latitude"]},{in_situ_site["longitude"]}')
#         print(f'[INFO] Path out: {in_situ_site["path_out"]}')
#
#     # time options
#     datetime_start, datetime_end, date_list = sextract.get_params_time(args, options)
#
#     for date in date_list:
#         product_list = sextract.get_list_products_day(basic_options['satellite_path_source'], date, wce,
#                                                       basic_options['org'])
#         for product in product_list:
#             launch_create_extract(in_situ_site, product, basic_options)
#
#
# # %%
# if __name__ == '__main__':
#     main()
