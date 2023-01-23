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
        # insitu path source
        self.insitu_path_source = None

        # dates
        self.start_date = dt.now().replace(month=1, day=1, hour=0, minute=0, second=0, microsecond=0)
        self.end_date = dt.now().replace(hour=23, minute=59, second=59)

        # param insitu
        self.param_insitu = None
        # param sat extract
        self.param_sat = None

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
        if self.options.has_option('Time_and_sites_selection', 'sites'):
            station_name = self.options['Time_and_sites_selection']['sites']
        else:
            print(f'[ERROR]In situ site is not defined in configuration file: Time_and_sites_selection/sites')
            self.VALID = False
            return
        in_situ_lat, in_situ_lon = cfs.get_lat_lon_ins(station_name)  # in situ location based on the station name
        if self.verbose:
            print(f'[INFO] Station name: {station_name} with lat: {in_situ_lat}, long: {in_situ_lon}')

        self.param_insitu = {
            'station_name': station_name,
            'in_situ_lat': in_situ_lat,
            'in_situ_lon': in_situ_lon
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
            station_name = self.param_insitu['station_name']
            if self.insitu_type == 'HYPERNETS':
                hday = INSITU_HYPERNETS_DAY(None)
                if hday.CHECK_SSH:
                    sday, eday = hday.get_start_and_end_dates(station_name)
                    if sday is not None and eday is not None:
                        if sday > self.start_date:
                            self.start_date = sday.replace(hour=0, minute=0, second=0, microsecond=0)
                        if eday < self.end_date:
                            self.end_date = eday.replace(hour=23, minute=59, second=59)