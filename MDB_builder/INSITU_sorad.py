import os, sys, __init__
from datetime import datetime as dt
from datetime import time
from datetime import timedelta

import pytz
from netCDF4 import Dataset
import numpy as np
import pandas as pd

try:
    import SB_support as SBr
    import INSITU_base as ISb
except:
    import MDB_builder.SB_support as SBr
    import MDB_builder.INSITU_base as ISb

code_home = os.path.dirname(os.path.dirname(__init__.__file__))


#sys.path.append(code_home)


class SO_RAD:

    def __init__(self, insitu_options, verbose):
        self.insitu_type = 'SINGLE_SEABASS'
        self.fixed_site = False
        self.verbose = verbose
        self.insitu_options = insitu_options
        self.site = self.insitu_options['site']
        self.valid = True
        code_monda = os.path.join(os.path.dirname(code_home), 'MONDA', 'src')
        if not os.path.isdir(code_monda):
            print(f'[ERROR] MONDA code folder {code_monda} is not available. To download the code use:')
            print(f'git clone https://github.com/monocle-h2020/MONDA.git')
            self.valid = False
        else:
            try:
                sys.path.append(code_monda)
                from monda.sorad import access, qc
            except:
                print(f'[ERROR] monda.sorad class could not be loaded')
                self.valid = False

        #self.platform_id = self.insitu_options['platform_id']
        self.layers = {
            'fp': 'rsg:sorad_public_view_fp_full',
            '3c': 'rsg:sorad_public_view_3c_full'
        }
        self.start_date = None
        self.end_date = None
        self.file_list = {}
        self.date_list = []

        ##RRS & RRS UNC
        self.wl_array = None
        self.nwl = 0
        #self.rrs_band_list = None
        #self.rrs_unc_band_list = None
        self.rrs_array = None
        #self.rrs_unc_array = None
        self.non_spectral_vars = None

        self.qc_params = None

    def prepare_data(self):
        if not self.valid:
            return False
        if self.start_date is None and self.end_date is None:
            print(
                f'[ERROR] Method to retrieve start_data and end_date automatically has not been implemented. Please use -sd(--start_date) and -ed (--end_date) arguments')
            return False
        date_array, time_list, lat_array, lon_array, indices_array = self.get_metadata_arrays()
        if date_array is None:
            return False
        sorad_folder = self.get_sorad_folder()
        if sorad_folder is None:
            return False
        date_array_unique = np.unique(date_array)
        for date_ts in date_array_unique:
            date_ts_obj = dt.strptime(date_ts, '%Y-%m-%d')
            if self.start_date is not None and self.end_date is not None:
                if date_ts_obj < self.start_date or date_ts_obj > self.end_date:
                    continue
            if self.verbose:
                print(f'[INFO] Getting metadata for date: {date_ts}')

            indices = np.where(date_array == date_ts)
            df = pd.DataFrame(index=np.arange(len(indices[0])),
                              columns=['insitu_index', 'insitu_time', 'insitu_lat', 'insitu_lon'])
            df['insitu_index'] = indices_array[indices]
            df['insitu_time'] = [x.strftime('%Y-%m-%d %H:%M:%S') for x in time_list[indices]]
            df['insitu_lat'] = lat_array[indices]
            df['insitu_lon'] = lon_array[indices]
            file_csv = os.path.join(sorad_folder, f'Metadata_{date_ts}.csv')
            df.to_csv(file_csv, sep=';', index=None)
            self.file_list[date_ts] = file_csv
            self.date_list.append(date_ts_obj)

        if len(self.date_list) == 0:
            print(f'[ERROR] No valid dates were retrieved from CSV files')
            return False

        self.date_list.sort()

        return True

    def get_metadata_date(self, datehere):
        date_ts = datehere.strftime('%Y-%m-%d')
        file_csv = self.file_list[date_ts]
        df = pd.read_csv(file_csv, sep=';')

        insitu_time = np.ma.array([dt.strptime(x, '%Y-%m-%d %H:%M:%S') for x in df['insitu_time']])
        insitu_lat = np.ma.array(df['insitu_lat'])
        insitu_lon = np.ma.array(df['insitu_lon'])
        insitu_indices = np.ma.array(df['insitu_index'])

        sorted_indices = np.argsort(insitu_time)
        insitu_time = insitu_time[sorted_indices]
        insitu_lat = insitu_lat[sorted_indices]
        insitu_lon = insitu_lon[sorted_indices]
        insitu_indices = (insitu_indices[sorted_indices],)

        return insitu_time, insitu_lat, insitu_lon, insitu_indices

    def prepare_csv_metadata(self, output_file):
        if not self.valid:
            return False
        if self.start_date is None and self.end_date is None:
            print(
                f'[ERROR] Method to retrieve start_data and end_date automatically has not been implemented. Please use -sd(--start_date) and -ed (--end_date) arguments')
            return False

        date_array, time_list, lat_array, lon_array, indices_array = self.get_metadata_arrays()
        if date_array is None:
            return False
        fw = open(output_file, 'w')
        started = False
        date_array_unique = np.unique(date_array)
        for date_ts in date_array_unique:
            date_here = dt.strptime(date_ts, '%Y-%m-%d')
            if self.start_date is not None and self.end_date is not None:
                if date_here < self.start_date or date_here > self.end_date:
                    if self.verbose:
                        print(
                            f'[INFO][SO_RAD] Skipping date {date_here.strftime("%Y-%m-%d")} as it not in the range: {self.start_date.strftime("%Y-%m-%d")} to {self.end_date.strftime("%Y-%m-%d")}')
                    continue
            if self.verbose:
                print(f'[INFO][SO_RAD] Getting metadata for date: {date_ts}')
            indices = np.where(date_array == date_ts)
            insitu_lat = lat_array[indices]
            insitu_lon = lon_array[indices]
            lat_min = np.ma.min(insitu_lat) - 1
            lat_max = np.ma.max(insitu_lat) + 1
            lon_min = np.ma.min(insitu_lon) - 1
            lon_max = np.ma.max(insitu_lon) + 1
            line = f'{date_ts},{lat_min},{lat_max},{lon_min},{lon_max}'
            if started:
                fw.write('\n')
            if not started:
                started = True
            fw.write(line)

        fw.close()
        if not started:
            print(f'[WARNING][SO_RAD] No data were found for the given data range in the SO_RAD dataset.')
            os.remove(output_file)
            return False
        else:
            return True

    def get_metadata_arrays(self):
        platform_id = self.insitu_options['platform_id']
        rrsalgorithm = self.insitu_options['rrs_algorithm']
        if rrsalgorithm in self.layers:
            layer = self.layers[rrsalgorithm]
        else:
            print(
                f'[ERROR] Layer is not defined for the rrs algorithm {rrsalgorithm}. Please choose the algorithm among {self.layers.keys()}')
            return [None] * 5
        from monda.sorad import access
        work_date = self.start_date
        date_array = []
        time_list = []
        lat_array = []
        lon_array = []
        indices_array = []
        while work_date <= self.end_date:
            day_start = dt.combine(work_date, time.min, pytz.utc)
            day_end = dt.combine(work_date, time.max, pytz.utc)
            try:
                response = access.get_wfs(platform=platform_id, timewindow=(day_start, day_end), layer=layer, bbox=None)
            except Exception as ex:
                print(f'[ERROR] Error retrieving So-Rad data. Please check you internet connection.')
                return [None] * 5
            if not response['result']:
                print(f'[WARNING] No So-Rad data retrieved  for {work_date.strftime("%Y-%m-%d)")}')
                work_date = work_date + timedelta(hours=24)

            for i in range(len(response['result'])):
                there = response['result'][i]['time']
                date_array.append(there.strftime("%Y-%m-%d"))
                time_list.append(there)
                lat_array.append(response['result'][i]['lat'])
                lon_array.append(response['result'][i]['lon'])
                indices_array.append(i)

            work_date = work_date + timedelta(hours=24)

        date_array = np.array(date_array)
        time_list = np.array(time_list)
        lat_array = np.array(lat_array)
        lon_array = np.array(lon_array)
        indices_array = np.array(indices_array)

        return date_array, time_list, lat_array, lon_array, indices_array

    def get_sorad_folder(self):
        path_s = os.path.join(self.insitu_options['path_sorad'], self.insitu_options['platform_id'])
        try:
            os.makedirs(path_s, exist_ok=True)
        except:
            print(f'[ERROR] So-Rad in situ path does not exists and could not be created. Please review permissions')
        return path_s

    def get_ref_date(self, datehere):
        name_ref = f'SO_RAD_{self.insitu_options["platform_id"]}_{self.insitu_options["rrs_algorithm"]}_{datehere.strftime("%Y%m%d")}'
        return name_ref

    ##checking are done in prepare_data
    def check_rrs_and_data_variables(self):
        if self.insitu_options['non_spectral_variables'] is not None:
            self.non_spectral_vars = self.insitu_options['non_spectral_variables']

        if self.insitu_options['builder_minimum_qc'] is None:
            print(f'[ERROR] builder_minimum_qc should be a interger value between -1 and 3')
            return False
        q_0_params = self.insitu_options['q0_params']
        q_1_params = self.insitu_options['q1_params']
        q_2_params = self.insitu_options['q2_params']
        q_3_params = self.insitu_options['q3_params']
        if q_0_params is None or len(q_0_params)!=4:
            return False
        if q_1_params is None or len(q_1_params)!=3:
            return False
        if q_2_params is None or len(q_2_params)!=2:
            return False
        if q_3_params is None or len(q_3_params)!=4:
            return False

        self.qc_params = {
            'q_0': {
                'lower_azi_bound': q_0_params[0],
                'upper_azi_bound': q_0_params[1],
                'upper_tilt_bound': q_0_params[2],
                'upper_tilt_std_bound': q_0_params[3]
            },
            'q_1': {
                'glint_th': q_1_params[0],
                'min_ed_th': q_1_params[1],
                'ls_th': q_1_params[2]
            },
            'q_2': {
                'rho_upperbound': q_2_params[0],
                'rmsd_tol': q_2_params[1]
            },
            'q_3': {
                'sim_upper_th': q_3_params[0],
                'sim_lower_th': q_3_params[1],
                'rrs_upper_th': q_3_params[2],
                'rrs_lower_th': q_3_params[3]
            }
        }


        return True

    def create_mini_mdb_files(self, options, extract_dir, extracts, time_extracts, overwrite):
        #print(options)
        fcsv = os.path.join(extract_dir, f'MDBm_{time_extracts[0].strftime("%Y%m%d")}.csv')
        fw = open(fcsv, 'w')
        started = False
        dims = None
        date_ref = ''
        platform_id = self.insitu_options['platform_id']
        rrs_algorithm = self.insitu_options['rrs_algorithm']
        layer = self.layers[rrs_algorithm]

        arrays_response = None
        for t in time_extracts:
            data_ref_here = t.strftime('%Y%m%d')
            if data_ref_here != date_ref:
                response = self.get_response_date(t,platform_id,layer)
                if response is None:
                    return
                if self.wl_array is None:
                    if not self.get_wavelength_from_response(rrs_algorithm,response):
                        return
                arrays_response = self.get_arrays_from_response(rrs_algorithm,response)
                if arrays_response is None:
                    return
                date_ref = data_ref_here

            options['n_insitu_bands'] = self.nwl
            ref = t.strftime('%Y%m%dT%H%M%S')

            file_out = ISb.get_mini_mdb_file_path(extract_dir, extracts[ref])
            write_line = False
            if os.path.isfile(file_out) and not overwrite:
                print(f'[INFO] Mini MDB file {os.path.basename(file_out)} already exists. Skipping...')
                write_line = True
            else:
                if self.verbose:
                    print(f'[INFO] Creating mini MDB file {os.path.basename(file_out)}')
                created = self.create_mini_mdb_file_impl(file_out, options, extracts[ref],arrays_response)
                if created and os.path.isfile(file_out):
                    write_line = True
                else:
                    print(f'[WARNING] File MDBm {os.path.basename(file_out)} could not be created. Skipping...')
            if write_line:
                fw, dims = self.add_line_csv_with_MDBm_info(fw, file_out, started)
                started = True

        fw.close()

        return dims

    def create_mini_mdb_file_impl(self, file_out, options, extract_info,arrays_response):
        #print(extract_info.keys())
        ninsitu_real = self.check_ninsitu_real(extract_info)
        if ninsitu_real < 0:
            return False
        if ninsitu_real == 0:
            print(f'[WARNING] No in situ data found for {extract_info["file"]}')
            return False

        if ninsitu_real > options['ninsitu_max']:
            print(
                f'[WARNING] {ninsitu_real} is greater than the maximum number of in situ data points {options["ninsitu_max"]}. ')
            return False
        sat_time = dt.fromtimestamp(extract_info['satellite_time']).astimezone(pytz.utc)
        if self.verbose:
            print(f'[INFO] Number of in situ data points for the extract: {ninsitu_real}')
            print(f'[INFO] Getting response for date: {sat_time.strftime("%Y-%m-%d")}')

        indices = extract_info['insitu_indices']
        minimum_qc = self.insitu_options['builder_minimum_qc']
        if minimum_qc>=0:
            if 'quality_control' in arrays_response:
                quality_control = arrays_response['quality_control']
                check = np.zeros(quality_control.shape)
                check[indices] = 1
                check[np.logical_and(check==1,quality_control>=minimum_qc)]=2
                indices = np.where(check==2)[0]
                ninsitu_real = len(indices)
                if self.verbose and ninsitu_real>0:
                    print(f'[INFO] Number of in situ data points for the extract after checking the quality control {ninsitu_real}')

        if ninsitu_real == 0:
            print(f'[WARNING] No in situ data found for {extract_info["file"]}')
            return False

        builder = ISb.Mini_MDB_Builder(options, self.verbose)
        builder.start_mini_mdb(extract_info['file'], file_out)
        builder.add_shipborne_variables()
        # if self.insitu_options['rrs_unc_format'] is not None:
        #     builder.add_rrs_uncentainty_variable()

        building_error = False
        if len(options['instrument_ids']) == 1:
            builder.set_insitu_wavelengths(0, self.wl_array)
            builder.set_instrument_id(ninsitu_real, 1)
        #basic variables: insitu_lat,insitu_lon,insitu_time,time_diff,insitu_spatial_index

        if not builder.set_insitu_basic_variables_from_dict(extract_info):
            print(f'[ERROR] Error setting the basic variables. MDBm file could not be created')
            building_error = True



        if 'quality_control' in arrays_response:
            var_name_nc = 'insitu_quality_control'
            array = arrays_response['quality_control']
            array = array[indices]
            qc_attrs = {
                'long_name': 'So-Rad Quality Control Flag Band',
                'flag_values': [-1,0,1,2,3],
                'flag_meanings': 'NO_QC QC_0 QC_1 QC_2 QC_3'
            }
            builder.add_non_spectral_variable(var_name_nc,qc_attrs)
            if not builder.set_non_spectral_variables(var_name_nc, array):
                print(f'[ERROR] Error setting the {var_name_nc} variable. MDBm file could not be created')
                building_error = True

        rrs_here = np.ma.transpose(self.rrs_array[indices, :])
        if not builder.set_spectral_variables('insitu_Rrs', rrs_here):
            print(f'[ERROR] Error setting the in situ Rrs variable. MDBm file could not be created')
            building_error = True

        spectral_vars = {
            'ed': {
                'long_name':'Downwelling irradiance (Ed)'
            },
            'ls':{
                'long_name': 'Sky radiance(Ls) or downwelling radiance(Ld)'
            },
            'lt':{
                'long_name': 'Total radiance from the sea (Lt) or upwelling radiance(Lu)'
            }

        }
        for var_name in spectral_vars:
            if var_name in arrays_response:
                var_name_nc = f'insitu_{var_name}'
                array = arrays_response[var_name]
                array_here = np.ma.transpose(array[indices,:])
                builder.add_spectral_variable(var_name_nc,spectral_vars[var_name])
                if not builder.set_spectral_variables(var_name_nc,array_here):
                    print(f'[ERROR] Error setting the {var_name_nc} variable. MDBm file could not be created')
                    building_error = True



        if self.non_spectral_vars is not None:
            for var_name in self.non_spectral_vars:
                if not var_name in arrays_response:
                    print(f'[WARNING] {var_name} could not be retrieved. Check potential values in hypernets_val/OPTIONS/insitu_options, [SO_RAD] section')
                    continue
                print(f'[INFO] Adding variable:',var_name)

                var_name_nc = f'insitu_{var_name}' if not var_name.startswith('insitu_') else var_name
                builder.add_non_spectral_variable(var_name_nc, None)
                array = arrays_response[var_name]
                array = array[indices]
                if not builder.set_non_spectral_variables(var_name_nc, array):
                    print(f'[ERROR] Error setting the {var_name_nc} variable. MDBm file could not be created')
                    building_error = True

        builder.close_mini_mdb_file()
        if building_error:
            os.remove(file_out)
            return False
        return True

    def check_ninsitu_real(self, extract_info):
        if not 'insitu_indices' in extract_info or extract_info['insitu_indices'] is None:
            print(f'[ERROR] insitu_indices is required in the extract info')
            return -1
        ninsitu_real = len(extract_info['insitu_indices'])
        if not 'insitu_lat' in extract_info or extract_info['insitu_lat'] is None:
            print(f'[ERROR] insitu_lat is required in the extract info')
            return -1
        if not 'insitu_lon' in extract_info or extract_info['insitu_lon'] is None:
            print(f'[ERROR] insitu_lon is required in the extract info')
            return -1
        if not 'insitu_time' in extract_info or extract_info['insitu_time'] is None:
            print(f'[ERROR] insitu_time is required in the extract info')
            return -1
        if not 'time_diff' in extract_info or extract_info['time_diff'] is None:
            print(f'[ERROR] time_diff is required in the extract info')
            return -1

        if len(extract_info['insitu_lat']) < ninsitu_real:
            print(
                f'[ERROR] Discrepancy in the number of in situ data points between insitu_lat {len(extract_info["insitu_lat"])} and insitu_indices {ninsitu_real}')
            return -1
        if len(extract_info['insitu_lon']) < ninsitu_real:
            print(
                f'[ERROR] Discrepancy in the number of in situ data points between insitu_lat {len(extract_info["insitu_lon"])} and insitu_indices {ninsitu_real}')
            return -1
        if len(extract_info['insitu_time']) < ninsitu_real:
            print(
                f'[ERROR] Discrepancy in the number of in situ data points between insitu_lat {len(extract_info["insitu_time"])} and insitu_indices {ninsitu_real}')
            return -1
        if len(extract_info['time_diff']) < ninsitu_real:
            print(
                f'[ERROR] Discrepancy in the number of in situ data points between time_diff {len(extract_info["time_diff"])} and insitu_indices {ninsitu_real}')
            return -1

        if 'insitu_spatial_index' in extract_info and extract_info['insitu_spatial_index'] is not None and len(
                extract_info['insitu_spatial_index']) < ninsitu_real:
            print(
                f'[ERROR] Discrepancy in the number of in situ data points between insitu_spatial_index {len(extract_info["insitu_spatial_index"])}and insitu_indices {ninsitu_real}')
            return -1

        return ninsitu_real

    def add_line_csv_with_MDBm_info(self, fw, file_nc, started):
        if not started:
            first_line = 'name;satellite_id;insitu_id;instrument_id;satellite_bands;insitu_bands;rows;columns'
            fw.write(first_line)
        dset = Dataset(file_nc)
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
        dims = np.array([nsat, ninsitu, ninstrument, nwlsat, nwlinsitu, rows, cols])

        return fw, dims

    def get_response_date(self, work_date, platform_id, layer):
        if platform_id is None:
            platform_id = self.insitu_options['platform_id']
        if layer is None:
            rrsalgorithm = self.insitu_options['rrs_algorithm']
            layer = self.layers[rrsalgorithm]

        day_start = dt.combine(work_date, time.min, pytz.utc)
        day_end = dt.combine(work_date, time.max, pytz.utc)
        from monda.sorad import access
        try:
            response = access.get_wfs(platform=platform_id, timewindow=(day_start, day_end), layer=layer, bbox=None)
            return response
        except Exception as ex:
            print(f'[ERROR] Error getting So-Rad data response from PML geo-server. Please review your internet connection.')
            return None


    def get_arrays_from_response(self, rrs_algorithm, response):
        arrays = None
        if self.get_rrs_array(rrs_algorithm, response):
            from monda.sorad import access
            arrays = {
                'ed': access.get_l1spectra(response,'ed_',self.wl_array),
                'ls': access.get_l1spectra(response,'ls_',self.wl_array),
                'lt': access.get_l1spectra(response,'lt_',self.wl_array)
            }
            ##non-spectral variables for quality control.
            qc_vars = ['time','rel_view_az','tilt_avg','tilt_std','c3_rmsd','c3_rho_ds','c3_rho_dd','c3_rho_s']
            for var_name in qc_vars:
                arrays[var_name] = self.get_non_spectral_variable(var_name, response)

            arrays['quality_control']= self.get_quality_control_array(rrs_algorithm,arrays)

            ##arrays in the list but not included in the quality control
            if self.non_spectral_vars is not None:
                for var_name in self.non_spectral_vars:
                    if var_name not in arrays:
                        arrays[var_name] = self.get_non_spectral_variable(var_name, response)

        return arrays

    def get_wavelength_from_response(self,rrs_algorithm,response):
        try:
            if rrs_algorithm == '3c':
                self.wl_array = np.arange(response['result'][0]['c3_wl_grid'][0], response['result'][0]['c3_wl_grid'][1],response['result'][0]['c3_wl_grid'][2])
                self.nwl = self.wl_array.shape[0]
        except Exception as ex:
            print(f'[ERROR] Exception {ex} retrieving the wavelenght array')
            return False
        if self.verbose:
            print(f'[INFO] Wavelenght array was extracted')
        return True

    def get_rrs_array(self, rrs_algorithm, response):
        try:
            if rrs_algorithm == '3c':
                self.rrs_array = np.array([response['result'][i]['c3_rrs'][:] for i in range(len(response['result']))])
        except Exception as ex:
            print(f'[ERROR] Exception {ex} retrieving the Rrs array')
            return False
        if self.verbose:
            print(f'[INFO] Rrs data were extracted')
        return True

    def get_non_spectral_variable(self, var_name, response):
        try:
            array = np.array([response['result'][i][var_name] for i in range(len(response['result']))])
            return array
        except:
            print(f'[ERROR] Exception {ex} retrieving the non-spectral variable {var_name}')
            return None

    def get_quality_control_array(self,rrs_algorithm,arrays):
        from monda.sorad import qc

        time = arrays['time']
        ed = arrays['ed']
        lt = arrays['lt']
        ls = arrays['ls']
        rel_view_az  = arrays['rel_view_az']
        tilt_avg  = arrays['tilt_avg']
        tilt_std = arrays['tilt_std']
        rmsd_3c = arrays['c3_rmsd']
        rho_ds = arrays['c3_rho_ds']
        rho_dd = arrays['c3_rho_dd']
        rho_s = arrays['c3_rho_s']


        #print(self.qc_params)

        q_az = qc.rel_az_filter(rel_view_az, lower_azi_bound=self.qc_params['q_0']['lower_azi_bound'], upper_azi_bound=self.qc_params['q_0']['upper_azi_bound'])
        q_tilt, q_tilt_std = qc.tilt_filter(tilt_avg, tilt_std, upper_tilt_bound=self.qc_params['q_0']['upper_tilt_bound'], upper_tilt_std_bound=self.qc_params['q_0']['upper_tilt_std_bound'])
        q_0 = qc.combined_filter(q_az, qc.combined_filter(q_tilt, q_tilt_std))

        ### Step (i): Radiometric quality control filters (QC applied to radiance or irradiance spectra)
        q_lt_ed = qc.qc_lt_ed_filter(ed, lt, time, self.wl_array, threshold=self.qc_params['q_1']['glint_th'])  # Lt/Ed ratio (glint) filtering
        q_ed = qc.qc_ed_filter(ed, min_ed_threshold=self.qc_params['q_1']['min_ed_th'])  # filters on =Rd and Ls anomalies
        q_ls = qc.qc_ls_filter(ls, self.wl_array, threshold=self.qc_params['q_1']['ls_th'])
        q_1 = qc.combined_filter(qc.combined_filter(qc.combined_filter(q_lt_ed, q_ed), q_ls),q_0)  # combined step 0 and step 1 mask


        #### Step (ii): QC filters specfic to 3C optimisation output
        if rrs_algorithm == '3c':
            q_rho = qc.qc_3c_rho_filter(rho_ds, rho_dd, rho_s,upperbound=self.qc_params['q_2']['rho_upperbound'])  # removes data where rho terminates at optimization bounds
            q_1_resid = qc.qc_3cresidual(q_1, rmsd_3c,tol=self.qc_params['q_2']['rmsd_tol'])  # removes data above threshold standard-deivation multiple
            q_2 = qc.combined_filter(q_rho, q_1_resid)  # combined step 0, step 1 and step 2 mask
        elif rrsalgorithm == 'fp':
            q_2 = np.nan * np.ones(len(q_1))

        #### Step (iii): Additional QC filters that apply to Rrs spectrum
        q_ss = qc.qc_ss_nir_filter(self.wl_array, self.rrs_array, upperthreshold=self.qc_params['q_3']['sim_upper_th'], lowerthreshold=self.qc_params['q_3']['sim_lower_th'])  # similarity spectrum filter
        q_maxrange = qc.qc_rrs_maxrange(self.rrs_array, upperthreshold=self.qc_params['q_3']['rrs_upper_th'], lowerthreshold=self.qc_params['q_3']['rrs_lower_th'])  # filters on max and min Rrs
        q_min = qc.qc_rrs_min(self.rrs_array, self.wl_array)

        # Optional filters
        # q_coastal = qc_coastalwater_rrsfilter(rrs, wl) #  filter based on expected shape of rrs - example from Warren 2019 used.
        # Users can input their own spectra here (will depend on water type)
        # q_var = qc_radiometric_variability(ed, lt, ls, time, wl, windowlength = 60, var_threshold =1.1, var_metric = 'zscore_max')

        if rrs_algorithm == '3c':  # recommended rrs qc mask for 3C method (combines steps (0), (i), (ii) and (iii))
            q_3 = qc.combined_filter(q_2, qc.combined_filter(q_min, (qc.combined_filter(q_ss, q_maxrange))))
        elif rrs_algorithm == 'fp':  # recommended rrs qc mask for fp method (combines steps (i) and (iii) QC)
            q_3 = qc.combined_filter(q_1, qc.combined_filter(q_min, (qc.combined_filter(q_ss, q_maxrange))))


        quality_control_array = q_0-1
        quality_control_array[q_1 == 1] = 1
        quality_control_array[q_2 == 1] = 2
        quality_control_array[q_3 == 1] = 3

        if self.verbose:
            n_q_0 = np.count_nonzero(quality_control_array >= 0)
            n_q_1 = np.count_nonzero(quality_control_array >= 1)
            n_q_2 = np.count_nonzero(quality_control_array >= 2)
            n_q_3 = np.count_nonzero(quality_control_array >= 3)
            print(f'[INFO] Quality control summary check. All: {quality_control_array.shape[0]} q_0: {n_q_0} q_1: {n_q_1} q_2: {n_q_2} q_3: {n_q_3}')

        return  quality_control_array
