import os
from datetime import datetime as dt
from INSITU_hypernets import INSITU_HYPERNETS_DAY


class MDBBuilderOptions:
    def __init__(self, options, verbose):
        self.options = options
        self.verbose = verbose

        self.VALID = True

        # path_out
        self.path_out = None
        # satellite_path_source
        self.list_mdbfiles_pathout = None
        self.satellite_path_source = None
        # insitu_type
        self.insitu_type = None
        # insitu_sensors
        self.insitu_sensors = {
            'HYPERNETS': 'HYPSTAR'
        }
        # insitu path source
        self.insitu_path_source = None

        # dates
        self.start_date = dt.now().replace(month=1, day=1, hour=0, minute=0, second=0, microsecond=0)
        self.end_date = dt.now().replace(hour=23, minute=59, second=59)

        # param insitu
        self.param_insitu = None
        # param sat extract
        self.param_sat = None
        # in situ options
        self.insitu_options = None

    def set_compulsory_options(self, cfs):
        if self.VALID:
            self.get_path_out()
        if self.VALID:
            self.get_satellite_path_source()
        if self.VALID:
            self.get_insitu_type()
        if self.VALID:
            self.get_param_insitu_site(cfs)
        if self.VALID:
            self.get_insitu_path_source()
        return self.VALID

    def get_path_out(self):
        if self.options.has_option('file_path', 'output_dir'):
            self.path_out = self.options['file_path']['output_dir']
        else:
            print(f'[ERROR] Path out [file_path]/[output_dir] is not defined in the configuration file')
            self.VALID = False
            return
        if not os.path.isdir(self.path_out):
            try:
                os.mkdir(self.path_out)
            except OSError:
                print(f'[ERROR] Path out {self.path_out} doest not exist and could not be created')
                self.VALID = False
                return

        if self.verbose:
            print(f'[INFO] Path to output: {self.path_out}')

        self.list_mdbfiles_pathout = []
        for name in os.listdir(self.path_out):
            self.list_mdbfiles_pathout.append(name)

    def get_satellite_path_source(self):
        if self.options.has_option('file_path', 'sat_extract_dir'):
            self.satellite_path_source = self.options['file_path']['sat_extract_dir']
        else:
            print(
                f'[ERROR] Satellite extract path [file_path]/[sat_extract_dir] is not defined in the configuration file')
            self.VALID = False
            return
        if not os.path.isdir(self.satellite_path_source):
            print(f'[ERROR] Satellite path source {self.satellite_path_source} is not a valid directory')
            self.VALID = False
            return
        if self.verbose:
            print(f'[INFO] Path to satellite sources: {self.satellite_path_source}')

    def get_insitu_type(self):
        if self.options.has_option('Time_and_sites_selection', 'insitu_type'):
            self.insitu_type = self.options['Time_and_sites_selection']['insitu_type']
        else:
            print(
                f'[ERROR] insitu type [Time_and_sites_selection]/[insitu_type] is not defined in the configuration file')
            self.VALID = False

    def get_insitu_path_source(self):
        if self.options['file_path']['ins_source_dir']:
            self.insitu_path_source = self.options['file_path']['ins_source_dir']
        else:
            print(
                f'[ERROR] In situ path source [file_path]/[ins_source_dir] is not defined in the configuration file')
            self.VALID = False
            return
        if not os.path.isdir(self.insitu_path_source):
            try:
                os.mkdir(self.insitu_path_source)
            except OSError:
                print(f'[ERROR] Insitu path source {self.insitu_path_source} doest not exist and could not be created')
                self.VALID = False
                return

    def get_param_insitu_site(self, cfs):
        if self.options.has_option('Time_and_sites_selection', 'site'):
            station_name = self.options['Time_and_sites_selection']['site']
        else:
            print(f'[ERROR]In situ site is not defined in configuration file: Time_and_sites_selection/site')
            self.VALID = False
            return
        in_situ_lat = None
        in_situ_lon = None
        if self.options.has_option('Time_and_sites_selection', 'site_lat') and self.options.has_option(
                'Time_and_sites_selection', 'site_long'):
            try:
                in_situ_lat = float(self.options['Time_and_sites_selection']['site_lat'])
                in_situ_lon = float(self.options['Time_and_sites_selection']['site_lat'])
            except:
                pass
        if in_situ_lat is None and in_situ_lon is None:
            in_situ_lat, in_situ_lon = cfs.get_lat_lon_ins(station_name)  # in situ location based on the station name
        if self.verbose:
            print(f'[INFO] Station name: {station_name} with lat: {in_situ_lat}, long: {in_situ_lon}')

        self.param_insitu = {
            'station_name': station_name,
            'insitu_lat': in_situ_lat,
            'insitu_lon': in_situ_lon
        }

    def get_param_sat_extracts(self):
        sat_sensor = 'SENSOR'
        sat_satellite = 'SAT'
        sat_platform = 'PLATFORM'
        sat_res = 'RESOLUTION'
        atm_corr = 'STANDARD'
        prefix = None
        if self.options.has_option('satellite_options', 'sensor'):
            sat_sensor = self.options['satellite_options']['sensor']
        if self.options.has_option('satellite_options', 'satellite'):
            sat_satellite = self.options['satellite_options']['satellite']
        if self.options.has_option('satellite_options', 'platform'):
            sat_platform = self.options['satellite_options']['platform']
        if self.options.has_option('satellite_options', 'resolution'):
            sat_res = self.options['satellite_options']['resolution']
        if self.options.has_option('satellite_options', 'ac'):
            atm_corr = self.options['satellite_options']['ac']
        if self.options.has_option('satellite_options', 'prefix'):
            prefix = self.options['satellite_options']['prefix']
        else:
            if (sat_satellite.upper() == 'S3' or sat_satellite.upper() == 'SENTINEL-3') and \
                    (sat_platform.upper().startswith('A') or sat_platform.upper().startswith('B')):
                prefix = f'{sat_satellite}{sat_platform}'  # wild card expression
            elif atm_corr == 'OLCI-L3':
                prefix = f'CMEMS2_O'

        if self.verbose:
            print(
                f'[INFO] Satellite extracts options----------------------------------------------------------------START')
            print(f'[INFO] Satellite: {sat_satellite.upper()}')
            print(f'[INFO] Satellite sensor: {sat_sensor.upper()}')
            print(f'[INFO] Satellite platform: {sat_platform.upper()}')
            print(f'[INFO] Satellite resolution: {sat_res.upper()}')
            print(f'[INFO] Satellite atmospheric correction: {atm_corr.upper()}')
            print(f'[INFO] Prefix: {prefix}')
            print(
                f'[INFO] Satellite extracts options------------------------------------------------------------------END')

        self.param_sat = {
            'satellite': sat_satellite.upper(),
            'sensor': sat_sensor.upper(),
            'platform': sat_platform.upper(),
            'resolution': sat_res.upper(),
            'ac': atm_corr.upper(),
            'prefix': prefix
        }

    def get_insitu_options(self):
        section = 'insitu_options'
        self.insitu_options = {
            'level': self.get_value_param(section,'level','L2A','str'),
            'apply_rsync': self.get_value_param(section,'apply_rsync',True,'boolean'),
            'insitu_site_flags': self.get_value_param(section,'insitu_site_flag_flags','INVALID','str'),
            'bad_spectra_file_list': self.get_value_param(section,'insitu_bad_spectra_file_list',None,'file'),
            'bad_spectra_prefix': self.get_value_param(section,'insitu_bad_spectra_prefix',None,'str'),
            'bad_spectra_format_time': self.get_value_param(section,'insitu_bad_spectra_format_time','%Y%m%dT%H%M%S','str')
        }

    def get_value(self, section, key):
        value = None
        if self.options.has_option(section, key):
            value = self.options[section][key]
        return value

    def get_value_param(self, section, key, default, type):
        value = self.get_value(section, key)
        if value is None:
            return default
        if type == 'str':
            return value
        if type == 'file':
            if not os.path.exists(value.strip()):
                return default
            else:
                return value.strip()
        if type == 'int':
            return int(value)
        if type == 'boolean':
            if value == '1' or value.upper() == 'TRUE':
                return True
            elif value == '0' or value.upper() == 'FALSE':
                return False
            else:
                return True
        if type == 'rrslist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                vals = vals.strip().replace('.', '_')
                list.append(f'RRS{vals}')
            return list
        if type == 'strlist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                list.append(vals.strip())
            return list
        if type == 'floatlist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                vals = vals.strip()
                list.append(float(vals))
            return list

    def get_mdb_extract_path(self, path_extract, ins_sensor):
        filename = f'MDB_{ins_sensor}_{path_extract}'
        outputpath = os.path.join(self.path_out, 'MDB_EXTRACTS')
        if not os.path.exists(outputpath):
            os.mkdir(outputpath)
        filepath = os.path.join(outputpath, filename)
        return filepath

    def get_sat_extracts_info(self):
        from datetime import datetime as dt
        sat_extract_info = {}
        sat_extract_info['time_min'] = self.start_date.strftime('%Y%m%dT%H%M%S')
        sat_extract_info['time_max'] = self.end_date.strftime('%Y%m%dT%H%M%S')
        sat_extract_info['time_creation'] = dt.now().strftime('%Y%m%dT%H%M%S')
        sat_extract_info['insitu_site_name'] = self.param_insitu['station_name']
        sat_extract_info['insitu_lat'] = self.param_insitu['insitu_lat']
        sat_extract_info['insitu_lon'] = self.param_insitu['insitu_lon']
        sat_extract_info['sensor'] = self.param_sat['sensor']
        sat_extract_info['satellite'] = self.param_sat['satellite']
        sat_extract_info['platform'] = self.param_sat['platform']
        sat_extract_info['resolution'] = self.param_sat['resolution']
        sat_extract_info['ac'] = self.param_sat['ac']
        sensor_str = 'UNKNONW'
        if self.insitu_type in self.insitu_sensors:
            sensor_str = self.insitu_sensors[self.insitu_type]
        sat_extract_info['ins_sensor'] = sensor_str

        return sat_extract_info

    def get_mdb_path(self):
        sat_extract_info = self.get_sat_extracts_info()
        datetime_min = sat_extract_info['time_min']
        datetime_max = sat_extract_info['time_max']
        sensor_str = sat_extract_info['sensor']
        station_name = sat_extract_info['insitu_site_name']
        resolution = sat_extract_info['resolution']
        platform = sat_extract_info['satellite'] + sat_extract_info['platform']
        ins_sensor = sat_extract_info['ins_sensor']
        # datetime_creation = sat_extract_info['time_creation']
        filename = f'MDB_{platform}_{sensor_str}_{resolution}_{datetime_min}_{datetime_max}_{ins_sensor}_{station_name}.nc'
        filepath = os.path.join(self.path_out, filename)
        return filepath

    def get_dates(self):
        if self.options['Time_and_sites_selection']['time_start']:
            self.start_date = dt.strptime(self.options['Time_and_sites_selection']['time_start'], '%Y-%m-%d').replace(
                hour=0, minute=0, second=0, microsecond=0)
        else:
            if 'sensor' in self.param_sat:
                if self.param_sat['sensor'] == 'OLCI':
                    self.start_date = dt(2016, 4, 15).replace(hour=0, minute=0, second=0, microsecond=0)
                if self.param_sat['sensor'] == 'MULTI':
                    self.start_date = dt(1997, 9, 1).replace(hour=0, minute=0, second=0, microsecond=0)
        if self.options['Time_and_sites_selection']['time_stop']:
            self.end_date = dt.strptime(self.options['Time_and_sites_selection']['time_stop'], '%Y-%m-%d').replace(
                hour=23, minute=59, second=59)

        ##checking dates with available in situ dates
        if self.insitu_type is not None and self.param_insitu is not None and 'station_name' in self.param_insitu:
            check_rsync = False
            if self.insitu_options is not None:
                if 'apply_rsync' in self.insitu_options:
                    check_rsync = self.insitu_options['apply_rsync']
            if not check_rsync:
                return

            station_name = self.param_insitu['station_name']
            if self.insitu_type == 'HYPERNETS':

                hday = INSITU_HYPERNETS_DAY(None, self.verbose)
                if hday.CHECK_SSH:
                    sday, eday = hday.get_start_and_end_dates(station_name, self.start_date, self.end_date)
                    if sday is not None and eday is not None:
                        if sday > self.start_date:
                            self.start_date = sday.replace(hour=0, minute=0, second=0, microsecond=0)
                        if eday < self.end_date:
                            self.end_date = eday.replace(hour=23, minute=59, second=59)
