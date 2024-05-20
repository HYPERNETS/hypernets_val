from INSITU_base import INSITUBASE
import subprocess, os
from datetime import datetime as dt
from datetime import timedelta


class INSITU_HYPERNETS_DAY(INSITUBASE):

    def __init__(self, mdb_options, rsync_user, verbose):
        self.mdb_options = mdb_options
        self.verbose = verbose
        self.new_MDB = None

        # Default: rsync_user: hypstar
        if mdb_options is not None:
            rsync_user = mdb_options.insitu_options['rsync_user']
        if rsync_user is None:
            rsync_user = 'hypstar'
        self.url_base = f'{rsync_user}@enhydra.naturalsciences.be'
        self.base_folder = '/waterhypernet/hypstar/processed_v2/'
        self.ssh_base = 'ssh -X -Y -p 9022'
        self.ls_base = 'ls processed_v2/'
        self.find_ref = 'HYPERNETS_W_SITE_L2A_REF*'

        self.rsync_base = f'rsync -a -e \'ssh -p 9022\' {self.url_base}:{self.base_folder}'
        self.rsync_url = f'rsync -a -e \'ssh -p 9022\' {self.url_base}'

        self.CHECK_SSH = self.check_ssh()

        self.insitu_extract_variables = {
            'insitu_quality_flag': {
                'name_orig': 'quality_flag',
                'type': 'u4',
                'standard_name': 'quality_flag',
                'long_name': 'A variable with the standard name of quality_flag contains an indication of assessed '
                             'quality information of another data variable. The linkage between the data variable and '
                             'the variable or variables with the standard_name of quality_flag is achieved using the '
                             'ancillary_variables attribute.',
                'flag_meanings': 'saturation nonlinearity bad_pointing placeholder1 lon_default lat_default outliers '
                                 'L0_thresholds L0_discontinuity dark_outliers vza_irradiance clear_sky_irradiance '
                                 'angles_missing lu_eq_missing fresnel_angle_missing fresnel_default '
                                 'temp_variability_ed temp_variability_lu min_nbred min_nbrlu min_nbrlsky '
                                 'def_wind_flag simil_fail',
                'flag_mask': '1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, '
                             '131072, 262144, 524288, 1048576, 2097152, 4194304 '
            },
            'insitu_viewing_azimuth_angle': {
                'name_orig': 'viewing_azimuth_angle',
                'type': 'f4',
                'standard_name': 'sensor_azimuth_angle',
                'long_name': 'sensor_azimuth_angle is the horizontal angle between the line of sight from the '
                             'observation point to the sensor and a reference direction at the observation point, '
                             'which is often due north. The angle is measured clockwise positive, starting from the '
                             'reference direction. A comment attribute should be added to a data variable with this '
                             'standard name to specify the reference direction.',
                'units': 'degrees',
                'reference': 'True North',
                'preferred_symbol': 'vaa'
            },
            'insitu_viewing_zenith_angle': {
                'name_orig': 'viewing_zenith_angle',
                'type': 'f4',
                'standard_name': 'sensor_zenith_angle',
                'long_name': 'sensor_zenith_angle is the angle between the line of sight to the sensor and the local '
                             'zenith at the observation target. This angle is measured starting from directly '
                             'overhead and its range is from zero (directly overhead the observation target) to 180 '
                             'degrees (directly below the observation target). Local zenith is a line perpendicular '
                             'to the Earth\'s surface at a given location. \'Observation target\' means a location on '
                             'the Earth defined by the sensor performing the observations.',
                'units': 'degrees',
                'preferred_symbol': 'vza'
            },
            'insitu_solar_azimuth_angle': {
                'name_orig': 'solar_azimuth_angle',
                'type': 'f4',
                'standard_name': 'solar_azimuth_angle',
                'long_name': 'Solar azimuth angle is the horizontal angle between the line of sight to the sun and a '
                             'reference direction which is often due north. The angle is measured clockwise.',
                'units': 'degrees',
                'reference': 'True North',
                'preferred_symbol': 'saa'
            },
            'insitu_solar_zenith_angle': {
                'name_orig': 'solar_zenith_angle',
                'type': 'f4',
                'standard_name': 'solar_zenith_angle',
                'long_name': 'Solar zenith angle is the the angle between the line of sight to the sun and the local '
                             'vertical.',
                'units': 'degrees',
                'preferred_symbol': 'sza'
            },
            'insitu_site_flag': {
                'name_orig': None,
                'type': 'i1',
                'standard_name': 'insitu_site_flag',
                'flag_meanings': 'INVALID',
                'flag_values': '1'
            }
        }
        self.insitu_spectral_variables = {
            'insitu_Rrs_nosc': {
                'name_orig': 'reflectance_nosc',
                'type': 'f4',
                'preferred_symbol': 'rhow_nosc',
                'standard_name': 'water_leaving_reflectance_nosc',
                'long_name': 'Reflectance of the water column at the surface without correction for the NIR similarity spectrum (see Ruddick et al., 2006)',
                'units': '-'

            }
        }

    # def add_insitu(self, extract_path, ofile):
    #     self.start_add_insitu(extract_path, ofile)

    def download_sequence_metadata(self,site,sequence_folder,output_folder):
        url_base_raw = f'hypstar@enhydra.naturalsciences.be'
        base_folder = '/home/hypstar/'
        rsync_url = f'rsync -a -e \'ssh -p 9022\' {url_base_raw}:{base_folder}{site}/DATA/{sequence_folder}/metadata.txt'
        cmd = f'{rsync_url} {output_folder}'
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(err)

    def create_mdb_insitu_extract(self, extract_path, ofile):
        self.start_add_insitu(extract_path, ofile)
        self.add_new_variables()

    def add_new_variables(self):
        for var_name in self.insitu_extract_variables:
            type = self.insitu_extract_variables[var_name]['type']
            fill_value_here = None
            if type == 'u4':
                fill_value_here = -999
            var = self.new_MDB.createVariable(var_name, type, ('satellite_id', 'insitu_id'), zlib=True, complevel=6,
                                              fill_value=fill_value_here)
            for at in self.insitu_extract_variables[var_name]:
                if at == 'type' or at == 'name_orig':
                    continue
                var.setncattr(at, self.insitu_extract_variables[var_name][at])
            if type == 'u4':
                var[:] = -999
            if self.insitu_extract_variables[var_name]['name_orig'] is None:
                var[:] = 0
        for var_name in self.insitu_spectral_variables:
            type = self.insitu_spectral_variables[var_name]['type']
            var = self.new_MDB.createVariable(var_name, type, ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                              zlib=True, complevel=6)
            for at in self.insitu_spectral_variables[var_name]:
                if at == 'type' or at == 'name_orig':
                    continue
                var.setncattr(at, self.insitu_spectral_variables[var_name][at])
            var[:] = -999
        # self.new_MDB.close()
        if self.verbose:
            print('[INFO] Added new variables')

    def set_data(self, inputpath, insitu_idx, sat_time, extract_info):
        from netCDF4 import Dataset
        import numpy as np
        import numpy.ma as ma
        nc_ins = Dataset(inputpath)
        file_name = inputpath.split('/')[-1]
        insitu_time_f = float(nc_ins.variables['acquisition_time'][0])
        if np.isnan(insitu_time_f):
            try:
                insitu_time = dt.strptime(file_name.split('_')[5], '%Y%m%dT%H%M')
            except:
                insitu_time = None
        else:
            insitu_time = dt.utcfromtimestamp(insitu_time_f)

        if insitu_time is None:
            print(f'[ERROR] In situ time was not defined for in situ file: {file_name}')
            nc_ins.close()
            return False

        time_diff = float(abs((sat_time - insitu_time).total_seconds()))
        time_max = self.mdb_options.insitu_options['time_window']  # time max in seconds
        time_diffh = time_diff / 3600
        time_maxh = time_max / 3600
        if time_diff > time_max:
            print(
                f'[WARNING] Time difference {time_diffh:.2f} hours is greater than max. time window ({time_maxh:.2f} hours). Skipping...')
            nc_ins.close()
            return False

        if self.verbose:
            print(f'[INFO] In situ path: {file_name}')
            print(f'[INFO] Sat. Time: {sat_time} Ins. Time: {insitu_time} Time diff.: {time_diffh:.2f} hours')

        # print(inputpath,insitu_time,sat_time,time_diff/3600)
        self.new_MDB.variables['insitu_time'][0, insitu_idx] = insitu_time_f
        self.new_MDB.variables['time_difference'][0, insitu_idx] = time_diff
        # self.new_MDB.variables['insitu_filename'][0, insitu_idx] = file_name DEPRECATED
        wini = 0
        wfin = 1600
        iini = 0
        ifin = 1600
        if nc_ins.variables['wavelength'].shape[0] < 1600:  # == 1537:
            iref = 1600 - nc_ins.variables['wavelength'].shape[0]
            wini = iref
            wfin = 1600
            iini = 0
            ifin = 1600 - iref
            # print(wini,wfin,iini,ifin)

        if insitu_idx == 0:
            self.new_MDB.variables['insitu_original_bands'][wini:wfin] = [nc_ins.variables['wavelength'][iini:ifin]]
            if wini > 0:
                wlref = self.new_MDB.variables['insitu_original_bands'][wini]
                for iw in range(wini, -1, -1):
                    self.new_MDB.variables['insitu_original_bands'][iw] = wlref
                    wlref = wlref - 0.50

        insitu_rhow_vec = [x for x, in nc_ins.variables['reflectance'][:]]
        insitu_RrsArray = ma.array(insitu_rhow_vec).transpose() / np.pi
        self.new_MDB.variables['insitu_Rrs'][0, wini:wfin, insitu_idx] = [insitu_RrsArray[iini:ifin]]

        for var_name in self.insitu_extract_variables:
            var_ins = self.insitu_extract_variables[var_name]['name_orig']
            if var_ins is not None:
                self.new_MDB.variables[var_name][0, insitu_idx] = [nc_ins.variables[var_ins][0]]

        for var_name in self.insitu_spectral_variables:
            var_ins = self.insitu_spectral_variables[var_name]['name_orig']
            var_array = ma.array(nc_ins.variables[var_ins][iini:ifin])  # [x for x in nc_ins.variables[var_ins][:]]
            # print('--->',var_array.shape)
            if var_name.find('Rrs') > 0:
                var_array = var_array / np.pi
            else:
                var_array = var_array
            # print('--->', var_array.shape)
            self.new_MDB.variables[var_name][0, wini:wfin, insitu_idx] = [var_array]
        nc_ins.close()

        if insitu_idx == 0:
            self.check_attributes(extract_info)

        return True

    def check_attributes(self, extract_info):
        if 'insitu_site_name' not in self.new_MDB.ncattrs():
            self.new_MDB.insitu_site_name = extract_info['insitu_site_name']
        if 'insitu_lat' not in self.new_MDB.ncattrs():
            self.new_MDB.insitu_site_name = extract_info['insitu_lat']
        if 'insitu_lon' not in self.new_MDB.ncattrs():
            self.new_MDB.insitu_site_name = extract_info['insitu_lon']
        if 'sensor' not in self.new_MDB.ncattrs():
            self.new_MDB.sensor = extract_info['sensor']
        if 'satellite' not in self.new_MDB.ncattrs():
            self.new_MDB.satellite = extract_info['satellite']
        if 'platform' not in self.new_MDB.ncattrs():
            self.new_MDB.platform = extract_info['platform']
        if 'resolution' not in self.new_MDB.ncattrs():
            self.new_MDB.resolution = extract_info['resolution']
        self.new_MDB.satellite_aco_processor = extract_info['ac']

    def close_mdb(self):
        self.new_MDB.close()

    def get_insitu_files(self, sat_time):
        site = self.mdb_options.param_insitu['station_name']
        level = self.mdb_options.insitu_options['level']

        pathbase = self.mdb_options.insitu_path_source
        if pathbase.split('/')[-1] != site:
            pathbase = os.path.join(pathbase, site)
        year_str = sat_time.strftime('%Y')
        month_str = sat_time.strftime('%m')
        day_str = sat_time.strftime('%d')
        path_day = os.path.join(pathbase, year_str, month_str, day_str)
        if not os.path.exists(path_day):
            return None

        list_files = []
        for name in os.listdir(path_day):
            if name.endswith('.nc') and name.find(level) > 0:
                list_files.append(os.path.join(path_day, name))
        if len(list_files) == 0:
            return None

        return list_files

    def create_path_day(self, time):
        site = self.mdb_options.param_insitu['station_name']
        pathbase = self.mdb_options.insitu_path_source
        if pathbase.split('/')[-1] == site:
            path_site = pathbase
        else:
            path_site = os.path.join(pathbase, site)
            if not os.path.exists(path_site):
                os.mkdir(path_site)
        year_str = time.strftime('%Y')
        path_year = os.path.join(path_site, year_str)
        if not os.path.exists(path_year):
            os.mkdir(path_year)
        month_str = time.strftime('%m')
        path_month = os.path.join(path_year, month_str)
        if not os.path.exists(path_month):
            os.mkdir(path_month)
        day_str = time.strftime('%d')
        path_day = os.path.join(path_month, day_str)
        if not os.path.exists(path_day):
            os.mkdir(path_day)
        return path_day

    def get_sequences_day_ssh(self,sitename,date_here):
        year_str = date_here.strftime('%Y')
        month_str = date_here.strftime('%m')
        day_str = date_here.strftime('%d')
        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{year_str}/{month_str}/{day_str}'
        sequence_list = self.get_list_sequence_folders(cmd)
        return sequence_list

    def get_files_day_ssh(self, sat_time, dotransfer):
        if self.verbose:
            print(f'[INFO] =====================================================================')
            print(f'[INFO] Getting Hypstar in situ files for date: {sat_time} via SSH...')
        sitename = self.mdb_options.param_insitu['station_name']
        level = self.mdb_options.insitu_options['level']
        year_str = sat_time.strftime('%Y')
        month_str = sat_time.strftime('%m')
        day_str = sat_time.strftime('%d')
        sat_time_min = sat_time - timedelta(hours=3)
        sat_time_max = sat_time + timedelta(hours=3)

        list_files_d = {}

        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{year_str}/{month_str}/{day_str}'
        sequence_list = self.get_list_sequence_folders(cmd)

        if len(sequence_list) > 0:
            for sequence in sequence_list:
                if self.verbose:
                    print(f'[INFO] Date: {year_str}-{month_str}-{day_str} Checking sequence folder: {sequence}')
                insitu_time_str = sequence[3:]
                insitu_time = dt.strptime(insitu_time_str, '%Y%m%dT%H%M%S')
                if sat_time_min <= insitu_time <= sat_time_max:
                    cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{year_str}/{month_str}/{day_str}/{sequence}/*.nc'
                    list_files = self.get_list_files(cmd)
                    for file in list_files:
                        if file.find(level) > 0:
                            list_files_d[insitu_time_str] = file

        if len(sequence_list) == 0:
            cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{year_str}/{month_str}/{day_str}/*.nc'
            list_files = self.get_list_files(cmd)
            if len(list_files) > 0:
                for file in list_files:
                    if file.find(level) > 0:
                        if self.verbose:
                            print(f'[INFO] Date: {year_str}-{month_str}-{day_str} Checking file: {file}')
                        insitu_time_str_basic = file.split('/')[-1].split('_')[5]
                        insitu_time = dt.strptime(insitu_time_str_basic, '%Y%m%dT%H%M').replace(second=0, microsecond=0)
                        insitu_time_str = insitu_time.strftime('%Y%m%dT%H%M%S')
                        if sat_time_min <= insitu_time <= sat_time_max:
                            list_files_d[insitu_time_str] = file

        if self.verbose:
            print(f'[INFO] {len(list_files_d)} files were found for {sat_time} via SSH')
            print(f'[INFO] =====================================================================')

        if len(list_files_d) > 0 and dotransfer:
            self.transfer_files_ssh(list_files_d)

        return list_files_d

    def transfer_files_ssh(self, list_files_d):
        if self.verbose:
            print(f'[INFO] =====================================================================')
            print(f'[INFO] Starting transfer of {len(list_files_d)} files via SSH...')

        for insitu_time_str in list_files_d:
            insitu_time = dt.strptime(insitu_time_str, '%Y%m%dT%H%M%S')
            path_day = self.create_path_day(insitu_time)
            input_path = list_files_d[insitu_time_str]
            name = input_path.split('/')[-1]
            output_path = os.path.join(path_day, name)
            if self.verbose:
                print(f'[INFO] Transfering file: {input_path} to {output_path}')
            cmd = f'{self.rsync_base}{input_path} {output_path}'
            prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            out, err = prog.communicate()
            if err:
                print(err)
        if self.verbose:
            print(f'[INFO] Transfering files completed')
            print(f'[INFO] =====================================================================')

    def transfer_files_to_output_folder_via_ssh(self, lista_files, output_folder):
        if self.verbose:
            print(f'[INFO] =====================================================================')
            print(f'[INFO] Starting transfer of {len(lista_files)} files via SSH...')
        for file in lista_files:
            if self.verbose:
                print(f'[INFO] Transfering file: {file} to {output_folder}')
            cmd = f'{self.rsync_url}:{file} {output_folder}'
            prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            out, err = prog.communicate()
            if err:
                print(err)
        if self.verbose:
            print(f'[INFO] Transfering files completed')
            print(f'[INFO] =====================================================================')

    def check_ssh(self):
        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}'
        try:
            subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT, timeout=10)
            return True
        except:
            print(f'[ERROR] Access to {self.url_base} via ssh is not allowed')
            return False

    def get_list_dates(self, sitename, start_date_ref, end_date_ref):
        list_dates = []
        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}'
        list_year = self.get_list_folder_dates(cmd)
        if len(list_year) > 0:
            for y in list_year:
                cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{y}'
                list_month = self.get_list_folder_dates(cmd)
                if len(list_month) > 0:
                    for m in list_month:
                        if self.verbose:
                            print(f'[INFO] Checking dates via SSH. Year: {y} Month: {m}')
                        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{y}/{m}'
                        list_days = self.get_list_folder_dates(cmd)
                        if len(list_days) > 0:
                            for d in list_days:
                                datehere_str = f'{y}{m}{d}'
                                datehere = dt.strptime(datehere_str, '%Y%m%d')
                                add_date = True
                                if start_date_ref is not None and datehere < start_date_ref:
                                    add_date = False
                                if end_date_ref is not None and datehere > end_date_ref:
                                    add_date = False
                                if add_date:
                                    list_dates.append(datehere)
        return list_dates

    def save_list_dates_to_file(self, fout, sitename, start_date_ref, end_date_ref, sat_extract_dir):
        list_dates = self.get_list_dates(sitename, start_date_ref, end_date_ref)
        if len(list_dates) == 0:
            print(f'[WARNING] Data were not found for site: {sitename}')
            return
        dates_extracts = []
        if sat_extract_dir is not None:
            for name in os.listdir(sat_extract_dir):
                if not name.endswith('.nc'):
                    continue
                fname = os.path.join(sat_extract_dir, name)
                from netCDF4 import Dataset
                dataset = Dataset(fname)
                if 'satellite_time' in dataset.variables:
                    sat_time = float(dataset.variables['satellite_time'][0])
                    sat_time_str = dt.utcfromtimestamp(sat_time).strftime('%Y-%m-%d')
                    dates_extracts.append(sat_time_str)
                dataset.close()
        f1 = open(fout, 'w')
        for date in list_dates:
            date_str = date.strftime('%Y-%m-%d')
            if date_str in dates_extracts:
                print(f'[INFO] Sat. extracts are already available for the date: {date_str}')
                continue
            f1.write(date_str)
            f1.write('\n')
        f1.close()

    def get_start_and_end_dates(self, sitename, start_date_ref, end_date_ref):
        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}'
        list_year = self.get_list_folder_dates(cmd)
        start_date = None
        end_date = None
        year_ini = start_date_ref.year
        month_ini = start_date_ref.month
        year_end = end_date_ref.year
        month_end = end_date_ref.month
        if len(list_year) > 0:
            for y in list_year:
                if int(y) < year_ini or int(y) > year_end:
                    continue
                cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{y}'
                list_month = self.get_list_folder_dates(cmd)
                if len(list_month) > 0:
                    for m in list_month:
                        if int(m) < month_ini or int(m) > month_end:
                            continue
                        if self.verbose:
                            print(f'[INFO] Checking dates via SSH. Year: {y} Month: {m}')
                        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{y}/{m}'
                        list_days = self.get_list_folder_dates(cmd)
                        if len(list_days) > 0:
                            for d in list_days:
                                datehere_str = f'{y}{m}{d}'
                                datehere = dt.strptime(datehere_str, '%Y%m%d')
                                if start_date is None:
                                    start_date = datehere
                                    end_date = datehere
                                else:
                                    if datehere < start_date:
                                        start_date = datehere
                                    if datehere > end_date:
                                        end_date = datehere

        return start_date, end_date

    def get_list_folder_dates(self, cmd):
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = prog.communicate()
        list = out.decode('utf-8').split('\n')
        listd = []
        for l in list:
            try:
                int(l)
                listd.append(l)
            except:
                pass
        return listd

    def get_list_sequence_folders(self, cmd):
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = prog.communicate()
        list = out.decode('utf-8').split('\n')
        listd = []
        for l in list:
            try:
                if l.startswith('SEQ'):
                    listd.append(l)
            except:
                pass
        return listd

    def get_list_files(self, cmd):
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = prog.communicate()
        list = out.decode('utf-8').split('\n')
        listd = []
        for l in list:
            if l == '':
                continue
            try:
                listd.append(l)
            except:
                pass
        return listd

    def get_files_download(self, date_here, site):
        folder_date = os.path.join(self.base_folder, site, date_here.strftime('%Y'), date_here.strftime('%m'),
                                   date_here.strftime('%d'))
        self.find_ref = self.find_ref.replace('SITE', site)

        cmd = f'{self.ssh_base} {self.url_base} find {folder_date} -name {self.find_ref}'

        list_files = self.get_list_files(cmd)

        if len(list_files) == 0:
            return None

        return list_files




