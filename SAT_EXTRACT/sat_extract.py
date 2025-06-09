import numpy as np
from netCDF4 import Dataset
from datetime import datetime as dt
import numpy.ma as ma
import configparser,os,subprocess,pytz,__init__,sys
code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
import COMMON.common_functions as cfs
from datetime import timedelta

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
        self.EXTRACT.satellite = satellite
        self.EXTRACT.platform = platform
        self.EXTRACT.sensor = sensor
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

    def create_flag_variable(self, var_name, var_array, description, flag_masks, flag_meanings, window):
        start_idx_y = window[0]
        stop_idx_y = window[1]
        start_idx_x = window[2]
        stop_idx_x = window[3]
        # Quality Flags
        satellite_flag = self.EXTRACT.createVariable(var_name, 'f4', ('satellite_id', 'rows', 'columns'),
                                                     fill_value=-999, zlib=True, complevel=6)
        satellite_flag[0, :, :] = [ma.array(var_array[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])]
        satellite_flag.description = description
        satellite_flag.flag_masks = flag_masks
        satellite_flag.flag_meanings = flag_meanings

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


class SatExtractOptions:
    def  __init__(self,config_file):
        self.options = None
        if config_file is not None:
            try:
                self.options = configparser.ConfigParser()
                self.options.read(config_file)
            except:
                self.options = None

def config_reader(FILEconfig):
    options = configparser.ConfigParser()
    options.read(FILEconfig)
    return options

def get_basic_options_from_file_config(args,options):
    # size box
    size_box = 25
    if options.has_option('satellite_options', 'extract_size'):
        size_box = int(options['satellite_options']['extract_size'])
    # resolution
    res = 'WFR'
    if options.has_option('satellite_options', 'resolution'):
        res = options['satellite_options']['resolution']
    # path_out
    path_out = None
    if options.has_option('file_path', 'output_dir'):
        path_out = options['file_path']['output_dir']
    else:
        if args.output:
            path_out = args.output
    if path_out is None:
        print(f'[ERROR] compulsory option path_out was not defined in the config file or argument output')
        return None
    if not os.path.isdir(path_out):
        path_out = create_dir(path_out)
        if path_out is None:
            print(f'[ERROR] path_out: {path_out} does not exist and could not be created')
            return None
    # makd_brdb
    make_brdf = False
    if options.has_option('satellite_options', 'brdf'):
        if options['satellite_options']['brdf'].upper() == 'T' or options['satellite_options'][
            'brdf'].upper() == 'TRUE':
            make_brdf = True
    # satellite_path_source
    satellite_path_source = None
    if options.has_option('file_path', 'sat_source_dir'):
        satellite_path_source = options['file_path']['sat_source_dir']
    else:
        if args.path_to_sat:
            satellite_path_source = args.path_to_sat
    if satellite_path_source is None:
        print(
            '[ERROR] compulsory option satellite_path_source was not defined in the config file or argument path_to_sat')
        return None
    if not os.path.exists(satellite_path_source):
        print(f'ERROR path: {satellite_path_source} does not exit')
        return None
    # tmp_path
    tmp_path = None
    if options.has_option('file_path', 'tmp_dir'):
        tmp_path = options['file_path']['tmp_dir']
        if os.path.exists(tmp_path) and os.path.isdir(tmp_path):
            tmp_path_del = os.path.join(tmp_path, '*')
            cmd = f'rm -r {tmp_path_del}'
            prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            out, err = prog.communicate()
            if err:
                print(err)
        if not os.path.exists(tmp_path) or not os.path.isdir(tmp_path):
            tmp_path = create_dir(tmp_path)

    if tmp_path is None:
        tmp_path = satellite_path_source
        print(
            f'[WARNING] tmp_path was not defined or does not exist. tmp path set to sat_path: {satellite_path_source}')
    # org
    org = None
    if options.has_option('file_path', 'sat_source_dir_organization'):
        org = options['file_path']['sat_source_dir_organization']

    # Ohter bands
    extra_bands = None
    if options.has_option('satellite_options', 'extra_bands'):
        str_val = options['satellite_options']['extra_bands']
        extra_bands = [x.strip() for x in str_val.split(',')]

    wce = None
    if options.has_option('satellite_options','wce'):
        wce = options['satellite_options']['wce'].strip()
        wce = f'"{wce}' if not wce.startswith(f'"') else f'{wce}'
        wce = f'{wce}"' if not wce.endswith(f'"') else f'{wce}'

    basic_options = {
        'satellite_path_source': satellite_path_source,
        'path_out': path_out,
        'tmp_path': tmp_path,
        'size_box': size_box,
        'make_brdf': make_brdf,
        'org': org,
        'resolution': res,
        'extra_bands': extra_bands,
        'wce': wce
    }
    return basic_options



def create_dir(path):
    try:
        os.mkdir(path)
        return path
    except:
        return None

def get_basic_options_from_arguments(args):
    ##parameters with default values
    size_box = 25
    make_brdf = False
    org = 'YYYY/jjj'

    res = 'WFR'
    if 'resolution' in args:
        res='WRR' if args.resolution=='WRR' else 'WFR'

    path_out = None
    if args.output:
        path_out = args.output
    if path_out is None:
        print(f'[ERROR] compulsory option path_out was not defined in the config file or argument output')
        return None
    if not os.path.isdir(path_out):
        path_out = create_dir(path_out)
        if path_out is None:
            print(f'[ERROR] path_out: {path_out} does not exist and could not be created')
            return None

    satellite_path_source = None
    if args.path_to_sat:
        satellite_path_source = args.path_to_sat
    if satellite_path_source is None:
        print(
            '[ERROR] compulsory option satellite_path_source was not defined in the config file or argument path_to_sat')
        return None
    if not os.path.exists(satellite_path_source):
        print(f'ERROR path: {satellite_path_source} does not exit')
        return None

    tmp_path = satellite_path_source

    wce = None
    if args.wce_expression:
        wce = args.wce_expression
        wce = f'"{wce}' if not wce.startswith(f'"') else f'{wce}'
        wce = f'{wce}"' if not wce.endswith(f'"') else f'{wce}'


    basic_options = {
        'satellite_path_source': satellite_path_source,
        'path_out': path_out,
        'tmp_path': tmp_path,
        'size_box': size_box,
        'make_brdf': make_brdf,
        'org': org,
        'resolution': res,
        'wce': wce
    }

    wce = f'"PACE_OCI*.L2.OC_AOP.V2_0.NRT.nc"'  # wild card expression

    return basic_options


def get_multiprocessing_options(options):
    mp_options ={
        'ncores': {
            'type':'int',
            'default': 0
        },
        'use_sbatch':{
            'type': 'boolean',
            'default': False
        },
        'sbatch_max_cores':{
            'type': 'int',
            'default': 8
        },
        'sbatch_launch':{
            'type': 'boolean',
            'default': True
        },
        'sbatch_partition':{
            'type': 'str',
            'default': 'octac_rep'
        },
        'sbatch_email':{
            'type': 'str',
            'default': None
        },
        'sbatch_email_type':{
            'type': 'str',
            'default': 'BEGIN,END,FAIL'
        },
        'conda_source':{
            'type':'str',
            'default': 'load_miniconda3.source'
        },
        'conda_env':{
            'type':'str',
            'default': 'op_proc_202211v2'
        }
    }
    options_out = {}

    if options.has_section('multiprocessing'):
        for key in mp_options:
            value = mp_options[key]['default']
            if options.has_option('multiprocessing',key):
                value = options['multiprocessing'][key].strip()
                type_v = mp_options[key]['type']
                if type_v=='int':
                    value = get_int_val(value)
                    if value is None:
                        print(f'[ERROR] Option {key} in multiprocessing section {value} is not a a valid integer')
                        return None
                if type_v=='boolean':
                    value = get_boolean_val(options['multiprocessing'][key])

            options_out[key] = value

    return options_out




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

def get_start_end_date_from_args(args):
    start_date = None
    end_date = None
    if args.startdate and args.enddate:
        try:
            start_date = dt.strptime(args.startdate, '%Y-%m-%d')
            end_date = dt.strptime(args.enddate, '%Y-%m-%d')
        except:
            print(f'[WARNING] --startdate (-sd) and/or --enddate (-ed) could not be parsed from {args.startdate} and/or {args.enddate}. Format should be: YYYY-mm-dd')
            start_date = None
            end_date = None
    return start_date,end_date

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


def overwrite(options):
    ow = False
    if options.has_option('file_path', 'overwrite'):
        ow = get_boolean_val(options['file_path']['overwrite'])
    return ow

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

def get_output_path(options,verbose):
    if options.has_option('file_path', 'output_dir'):
        output_path = options['file_path']['output_dir']
        if not os.path.isdir(output_path):
            if verbose:
                print(f'[WARNING] Trying to create the output path: {output_path}')
            output_path = create_dir(output_path)
        if output_path is None:
            print(f'[ERROR] Output {output_path} is not available and could not be created')
            return None
        if verbose:
            print(f'[INFO] Output path:{output_path}')
        return output_path
    else:
        print('[ERROR] section:file_path, option: output_dir is not included in the configuration file')
        return None

def get_input_path_info(options):
    path_source = options['file_path']['sat_source_dir'].strip()
    if not os.path.isdir(path_source):
        print(f'[ERROR] Source path {path_source} is not available or is not a valid directory')
        return None
    org = None
    if options.has_option('file_path', 'sat_source_dir_organization') and options['file_path'][
        'sat_source_dir_organization']:
        org = options['file_path']['sat_source_dir_organization']

    if org is not None:
        if org.lower()=='none':
            org = None
        elif org=='YYYYmmdd':
            org = '%Y/%m/%d'
        elif org == 'YYYYjjj':
            org = '%Y/%j'
        elif org == 'YYYYmm':
            org = '%Y/%m'
        else:
            org = org.replace('YYYY','%Y')
            org = org.replace('mm', '%m')
            org = org.replace('dd', '%d')
            org = org.replace('jjj', '%j')
        try:
           
            dt.now().strftime(org)
        except:
            print(f'[ERROR] Satellite path date organization {org} is not valid')
            return None


    wce = '*'
    if options.has_option('file_path', 'wce') and options['file_path']['wce']:
        wce = options['file_path']['wce']

    input_path_info = {
        'path_source': path_source,
        'org': org,
        'wce': wce
    }
    return input_path_info

def get_cnr_extract_options(options):
    section = 'CNR_OPTIONS'
    n_bands = 0

    dataset_name_file = None
    dataset_name_format_date = '%Y%j'
    if options.has_option(section,'dataset_name_file'):
        dataset_name_file = options[section]['dataset_name_file']
    if options.has_option(section,'dataset_name_format_date'):
        dataset_name_format_date = options[section]['dataset_name_format_date']

    dataset_var_list = None
    dataset_var_list_out = None
    if options.has_option(section, 'dataset_var_list'):
        s = options[section]['dataset_var_list']
        dataset_var_list = [x.strip() for x in s.split(',')]
        dataset_var_list_out = dataset_var_list
        if options.has_option(section, 'dataset_var_list_out'):
            s = options[section]['dataset_var_list_out']
            dataset_var_list_out = [x.strip() for x in s.split(',')]

    rrs_list = []
    is_reflectance = True
    for var in dataset_var_list:
        try:
            value = float(var)
            rrs_list.append(value)
        except:
            is_reflectance = False
    if is_reflectance:
            n_bands = len(dataset_var_list)

    options_out = {
        'is_reflectance': is_reflectance,
        'n_bands': n_bands,
        'dataset_name_file': dataset_name_file,
        'dataset_name_format_date': dataset_name_format_date,
        'dataset_var_list': dataset_var_list,
        'dataset_var_list_out': dataset_var_list_out,
        'rrs_list': rrs_list,
        'size_box': get_box_size(options)
    }

    return options_out

def get_csv_options(options, section):
    #if section == 'CSV_SELECTION':
    col_date = 'date'
    col_lat = 'lat'
    col_lon = 'lon'
    col_sep = ';'
    format_date = '%Y-%m-%dT%H:%M'
    if section == 'MULTIPLE_CSV_SELECTION':
        col_date = 'timestamp'
        col_lat = 'lat'
        col_lon = 'lon'
        format_date = '%Y-%m-%d %H:%M:%S.%f'
        col_sep = ','

    if options.has_option(section, 'col_date'):
        col_date = options[section]['col_date']
    if options.has_option(section, 'col_lat'):
        col_lat = options[section]['col_lat']
    if options.has_option(section, 'col_lon'):
        col_lon = options[section]['col_lon']
    if options.has_option(section, 'format_date'):
        format_date = options[section]['format_date']
    if options.has_option(section, 'col_sep'):
        col_sep = options[section]['col_sep']

    col_time = None
    format_time = '%H:%M:%S'
    if options.has_option(section, 'col_time'):
        col_time = options[section]['col_time']
    if options.has_option(section, 'format_time'):
        format_time = options[section]['format_time']

    return col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep

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


def get_box_size(options):
    if options.has_option('satellite_options', 'extract_size'):
        try:
            size_box = int(options['satellite_options']['extract_size'])
        except:
            size_box = 25
    else:
        size_box = 25
    return size_box

def get_geo_info(options, file_nc, insitu_lat, insitu_lon, lat, lon):
    if lat is None or lon is None:
        lat, lon = get_lat_lon_arrays(options, file_nc)

    contain_flag = 0
    limits = None
    rc = None
    if cfs.contain_location(lat, lon, insitu_lat, insitu_lon) == 1:
        if lat.ndim == 1 and lon.ndim == 1:
            r = np.argmin(np.abs(lat - insitu_lat))
            c = np.argmin(np.abs(lon - insitu_lon))
        else:
            r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
        if 'size_box' in options:
            size_box = options['size_box']
        else:
            size_box = get_box_size(options)
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
    only_date_array_unique = np.unique(only_date_array).tolist()
    return only_date_array_unique, time_list

def get_lat_lon_arrays(options, file_nc):
    nc_sat = Dataset(file_nc, 'r')
    var_lat, var_lon = get_lat_long_var_names(options)
    lat = nc_sat.variables[var_lat][:]
    lon = nc_sat.variables[var_lon][:]
    #lat, lon = get_lat_long_arrays(nc_sat, var_lat, var_lon)
    nc_sat.close()
    return lat, lon

def get_lat_long_var_names(options):
    var_lat = 'lat'
    var_lon = 'lon'
    if options.has_option('satellite_options', 'lat_variable'):
        var_lat = options['satellite_options']['lat_variable']
    if options.has_option('satellite_options', 'lon_variable'):
        var_lon = options['satellite_options']['lon_variable']
    return var_lat, var_lon

def get_satellite_time_l3(options,fproduct,daydate):
    cmems_time = '11:00'
    if options.has_option('satellite_options', 'satellite_time'):
        cmems_time = options['satellite_options']['satellite_time'].strip()
    satellite_time = get_satellite_time_from_global_attributes(fproduct)
    if satellite_time is None:
        try:
            satellite_time = dt.strptime(f'{daydate.strftime("%Y%m%d")}T{cmems_time}',
                                         '%Y%m%dT%H:%M').replace(tzinfo=pytz.utc)
        except:
            print(f'{cmems_time} is not a valid satellite time option. Skipping')
            return None
    return satellite_time

def get_satellite_time_from_global_attributes(fproduct):
    dataset = Dataset(fproduct)
    if 'time_coverage_start' in dataset.ncattrs() and 'time_coverage_end' in dataset.ncattrs():
        try:
            start_date = dt.strptime(dataset.time_coverage_start, '%d-%b-%Y %H:%M:%S.%f').replace(tzinfo=pytz.utc)
            end_date = dt.strptime(dataset.time_coverage_end, '%d-%b-%Y %H:%M:%S.%f').replace(tzinfo=pytz.utc)
            sat_time = (start_date.timestamp() + end_date.timestamp()) / 2
            sat_time = dt.fromtimestamp(sat_time).astimezone(pytz.utc)
            dataset.close()
            return sat_time

        except:
            dataset.close()
            return None
    dataset.close()
    return None

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
    ref_keys = ['satellite', 'platform', 'sensor', 'res', 'aco_processor', 'proc_version']
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
            print(f'[ERROR] Path date for date {date_here.strftime("%Y-%m-%d")} could not be parsed using {org} organization')
            return None

    if not os.path.isdir(path_date):
        if not createIfNotExist:
            print(f'[ERROR] Path date {path_date} does not exist or is not a valid directory')
            return None
        else:
            path_date = create_dir(path_date)
            if path_date is None:
                print(f'[ERROR] Path date {path_date} does not exist and could not be created. Please review permissions')
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

# def get_lat_long_arrays(nc_sat, var_lat, var_lon):
#     vlat = nc_sat.variables[var_lat]
#     vlon = nc_sat.variables[var_lon]
#     if vlat.ndim == 2 and vlon.ndim == 2:
#         lat = nc_sat.variables[var_lat][:, :]
#         lon = nc_sat.variables[var_lon][:, :]
#     if vlat.ndim == 1 and vlon.ndim == 1:
#         lat = nc_sat.variables[var_lat][:]
#         lon = nc_sat.variables[var_lon][:]
#     return lat, lon