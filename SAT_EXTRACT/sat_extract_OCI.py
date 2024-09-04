import os, sys, argparse, __init__

import pytz

code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
from SAT_EXTRACT.sat_extract import SatExtract
import COMMON.common_functions as cfs
import sat_extract as sextract
from netCDF4 import Dataset
from datetime import datetime as dt
from datetime import timedelta

os.environ['QT_QPA_PLATFORM'] = 'offscreen'  # to avoid error "QXcbConnection: Could not connect to display"
path2ncrcat = '/opt/local/bin/ncrcat'

parser = argparse.ArgumentParser(description="Create sat extract from PACE-OCI AOP files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-sd', "--startdate", help="The Start Date - format YYYY-MM-DD ")
parser.add_argument('-ed', "--enddate", help="The End Date - format YYYY-MM-DD ")
parser.add_argument('-site', "--sitename", help="Site name.", choices=['VEIT', 'BEFR', 'BSBE'])
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-ps', "--path_to_sat", help="Path to satellite sources.")
parser.add_argument('-o', "--output", help="Path to output")
# parser.add_argument('-nl', "--nolist",help="Do not create initial satellite lists, checking day by day allowing download",action="store_true")
# parser.add_argument('-adownload', "--allow_download", help="Allow download", action="store_true")

args = parser.parse_args()


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
        dataset = Dataset(self.path_source)
        var_array = dataset[group].variables[variable][window[0]:window[1],window[2]:window[3]]
        attrs = dataset[group].variables[variable].__dict__
        dataset.close()
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


class SatExtractOCI(SatExtract):

    def __init__(self, ofname, variable_list):
        if ofname is not None:
            SatExtract.__init__(self, ofname)
        if variable_list is not None and len(variable_list)>0:
            self.variable_list = variable_list
        else:
            self.variable_list = {
                'geophysical_data': ['angstrom','aot_865','avw']
            }
        self.sensor = 'oci'

    def create_rrs_oci_variables(self,rrs_t,rrs_unc_t):
        self.create_rrs_variable('OCI')
        self.create_rrs_unc_variable('OCI')
        b1 = self.set_rrs_data('satellite_Rrs',rrs_t[0],rrs_t[1])
        b2 = self.set_rrs_data('satellite_Rrs_unc', rrs_unc_t[0], rrs_unc_t[1])
        if b1 and b2:
            return True
        else:
            return False

    def set_rrs_data(self,var_name,var_array,attrs):
        nbands = self.EXTRACT.dimensions['satellite_bands'].size
        if nbands==var_array.shape[2]:
            for iband in range(nbands):
                self.EXTRACT.variables[var_name][0,iband,:,:] = var_array[:,:,iband]

            for at in attrs:##attributes, long_name is already set
                if at == '_FillValue' or at == 'scale_factor' or at == 'add_offset' or at=='long_name':
                    continue
                val = attrs[at]
                if (at == 'valid_min' or at == 'valid_max') and 'scale_factor' in attrs and 'add_offset' in attrs:
                    val = (attrs[at] * attrs['scale_factor']) + attrs['add_offset']
                self.EXTRACT.variables[var_name].setncattr(at, val)
            return True
        else:
            print(f'[ERROR] {var_name} data could not be set, number of bands do not coincide with the satellelit_bands dimensions')
            return False

    def create_l2_oci_flag_variable(self,array,attrs):
        satellite_flag = self.EXTRACT.createVariable('satellite_l2_flags','i4',('satellite_id', 'rows', 'columns'),zlib = True, complevel=6)

        satellite_flag[0, :, :] = array[:,:]
        satellite_flag.setncatts(attrs)

    def create_2D_oci_variable(self,name,array,attrs):
        if not name.startswith('satellite_'): name = f'satellite_{name}'
        satellite_variable = self.EXTRACT.createVariable(name,'f4',('satellite_id', 'rows', 'columns'),fill_value=-999,zlib = True, complevel=6)
        satellite_variable[0,:,:] = array[:,:]
        for at in attrs:
            if at=='_FillValue' or  at=='scale_factor' or at=='add_offset':
                continue
            val = attrs[at]
            if (at=='valid_min' or at=='valid_max') and 'scale_factor' in attrs and 'add_offset' in attrs:
                val = (attrs[at] * attrs['scale_factor']) + attrs['add_offset']


            satellite_variable.setncattr(at,val)

    def create_tilt_variable(self,tilt_value):
        satellite_tilt = self.EXTRACT.createVariable('satellite_tilt', 'f4', ('satellite_id'), fill_value=-999,
                                                     zlib=True, complevel=6)
        # print('Satellite start time es: ',satellite_start_time)
        satellite_tilt[0] = float(tilt_value)
        satellite_tilt.long_name = 'Sensor tilt angle'
        satellite_tilt.units = 'degrees'
        satellite_tilt.valid_min = -25.0
        satellite_tilt.valid_max = 25.0

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


# %%
def main():
    print('Creating satellite extracts.')

    # Getting basic_options from config_file or arguments
    basic_options = None
    options = None
    if args.config_file:
        if os.path.isfile(args.config_file):
            options = sextract.config_reader(args.config_file)
            basic_options = sextract.get_basic_options_from_file_config(args, options)
    if basic_options is None:
        basic_options = sextract.get_basic_options_from_arguments(args)
    if basic_options is None:
        return
    if args.verbose:
        for option in basic_options:
            print(f'[INFO] {option}:{basic_options[option]}')

    # size_box = basic_options['size_box']
    # res = basic_options['resolution']
    # path_out = basic_options['path_out']
    # make_brdf = basic_options['make_brdf']
    # satellite_path_source = basic_options['satellite_path_source']
    # tmp_path = basic_options['tmp_path']
    # org = basic_options['org']
    wce = f'"PACE_OCI*.L2.OC_AOP.V2_0.NRT.nc"'  # wild card expression

    in_situ_site = sextract.get_insitu_site(args, options, basic_options['path_out'])
    if args.verbose:
        print(
            f'[INFO] In situ site: {in_situ_site["site_name"]} {in_situ_site["latitude"]},{in_situ_site["longitude"]}')
        print(f'[INFO] Path out: {in_situ_site["path_out"]}')

    # time options
    datetime_start, datetime_end, date_list = sextract.get_params_time(args, options)

    for date in date_list:
        product_list = sextract.get_list_products_day(basic_options['satellite_path_source'], date, wce,
                                                      basic_options['org'])
        for product in product_list:
            launch_create_extract(in_situ_site, product, basic_options)


# %%
if __name__ == '__main__':
    main()
