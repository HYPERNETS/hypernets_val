import argparse
import configparser
import os
import sys
from netCDF4 import Dataset
from datetime import datetime as dt
from datetime import timedelta
import numpy.ma as ma

# import user defined functions from other .py
code_home = os.path.abspath('../')
sys.path.append(code_home)
from SAT_EXTRACT.sites_mng import SitesMng
from SAT_EXTRACT.sat_extract import SatExtract
import COMMON.common_functions as cfs

parser = argparse.ArgumentParser(description="Create Match-up DataBase files (MDB) files.")
parser.add_argument("-d", "--debug", help="Debugging mode.", action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-p', "--product_file", help="Image file")

args = parser.parse_args()


def launch_create_extract(filepath, options):
    # Retrieving sites
    smng = SitesMng()
    path_output = get_output_path(options)
    site_file, site_list, region_list = get_site_options(options)
    if not site_file is None:
        sites = smng.get_sites_from_file(site_file, site_list, region_list, path_output)
    else:
        sites = smng.get_sites_from_list(site_list, path_output)

    if len(sites) == 0:
        print('WARNING: Site is not defined...')

    # Start datase
    nc_sat = Dataset(filepath, 'r')

    # Retriving lat and long arrays
    lat, lon = get_lat_long_arrays(options, nc_sat)

    # Retrieving global atribbutes
    global_at = get_global_atrib(nc_sat)

    # Working for each site, checking if there is in the image
    for site in sites:
        insitu_lat = sites[site]['latitude']
        insitu_lon = sites[site]['longitude']
        contain_flag = 0
        if cfs.contain_location(lat, lon, insitu_lat, insitu_lon) == 1:
            r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
            if r >= 0 and r + 1 < lat.shape[0] and c >= 0 and c + 1 < lat.shape[1]:
                contain_flag = 1
        if contain_flag == 1:
            filename = filepath.split('/')[-1].replace('.', '_') + '_extract_' + site + '.nc'
            global_at['pdu'] = [filepath.split('/')[-1]]

            if not os.path.exists(sites[site]['path_out']):
                os.mkdir(sites[site]['path_out'])
            ofname = os.path.join(sites[site]['path_out'], filename)
            global_at['station_name'] = site
            global_at['in_situ_lat'] = insitu_lat
            global_at['in_situ_lon'] = insitu_lon
            create_extract(ofname, options, nc_sat, global_at, lat, lon, r, c)
        else:
            print('Station out of the image')


def create_extract(ofname, options, nc_sat, global_at, lat, long, r, c):
    size_box = get_box_size(options)
    start_idx_x = (c - int(size_box / 2))
    stop_idx_x = (c + int(size_box / 2) + 1)
    start_idx_y = (r - int(size_box / 2))
    stop_idx_y = (r + int(size_box / 2) + 1)
    window = [start_idx_y, stop_idx_y, start_idx_x, stop_idx_x]

    newEXTRACT = SatExtract(ofname)
    if not newEXTRACT.FILE_CREATED:
        print('File not created')
        return

    newEXTRACT.set_global_attributes(global_at)
    reflectance_bands, reflectance_bands_group, n_bands = get_reflectance_bands_info(options)
    newEXTRACT.create_dimensions(size_box, n_bands)

    newEXTRACT.create_lat_long_variables(lat, long, window)

    # Sat time start:  ,+9-2021-12-24T18:23:00.471Z
    newEXTRACT.create_satellite_time_variable(dt.strptime(nc_sat.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ'))

    # Rrs and wavelenghts
    satellite_Rrs = newEXTRACT.create_rrs_variable(global_at['sensor'])
    rrs_group = nc_sat.groups[reflectance_bands_group]
    rbands = list(reflectance_bands.keys())
    wavelenghts = []
    for index in range(len(rbands)):
        rband = rbands[index]
        bandarray = ma.array(rrs_group.variables[rband][:, :])
        satellite_Rrs[0, index, :, :] = bandarray[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
        wl = reflectance_bands[rband]['wavelenght']
        wavelenghts.append(wl)
    newEXTRACT.create_satellite_bands_variable(wavelenghts)

    # flags
    flag_band_name, flag_band_group = get_flags_info(options)
    flag_group = nc_sat.groups[flag_band_group]
    flag_band = flag_group.variables[flag_band_name]
    newEXTRACT.create_flag_variable(f'satellite_{flag_band_name}', flag_band, flag_band.long_name, flag_band.flag_masks,
                                    flag_band.flag_meanings, window)


def config_reader(FILEconfig):
    """
    Reads and checks configuration file for the validation
    Args:
        FILEconfig(str): configuration file path

    Return:
        options object
    """
    options = configparser.ConfigParser()
    options.read(FILEconfig)
    return options


def get_output_path(options):
    if 'file_path' in options and 'output_dir' in options['file_path']:
        return options['file_path']['output_dir']
    else:
        return None


def get_box_size(options):
    if 'satellite_options' in options and 'extract_size' in options['satellite_options']:
        size_box = int(options['satellite_options']['extract_size'])
    else:
        size_box = 25
    return size_box


def get_site_options(options):
    site_file = None
    site_list = None
    region_list = None
    if options['Time_and_sites_selection']['sites_file']:
        site_file = options['Time_and_sites_selection']['sites_file']
    if options['Time_and_sites_selection']['sites']:
        site_list = options['Time_and_sites_selection']['sites'].split(',')
        site_list = [s.strip() for s in site_list]
    if not options['Time_and_sites_selection']['sites'] and options['Time_and_sites_selection']['sites_region']:
        region_list = options['Time_and_sites_selection']['sites_region'].split(',')
        region_list = [r.strip() for r in region_list]
    return site_file, site_list, region_list


def get_lat_long_arrays(options, nc_sat):
    geo_group = nc_sat.groups['navigation_data']
    lat = geo_group.variables['latitude'][:, :]
    lon = geo_group.variables['longitude'][:, :]
    return lat, lon


def get_global_atrib(nc_sat):
    at = {'sensor': nc_sat.instrument}
    if nc_sat.platform == 'JPSS-1':
        at['satellite'] = 'NOAA-20/JPSS-1'
        at['platform'] = 'JPSS1'
    elif nc_sat.platform == 'Suomi-NPP':
        at['satellite'] = nc_sat.platform
        at['platform'] = 'SNPP'
    else:
        at['platform'] = nc_sat.platform
    at['res'] = ''
    proc_group = nc_sat.groups['processing_control']
    at['aco_processor'] = f'{proc_group.software_name} {proc_group.software_version}'
    at['proc_version'] = nc_sat.processing_version

    at['station_name'] = ''
    at['in_situ_lat'] = -999
    at['in_situ_lon'] = -999
    at['pdu'] = ''

    return at


def get_reflectance_bands_info(options):
    reflectance_bands = {
        'Rrs_411': {'wavelenght': 411},
        'Rrs_445': {'wavelenght': 445},
        'Rrs_489': {'wavelenght': 489},
        'Rrs_556': {'wavelenght': 556},
        'Rrs_667': {'wavelenght': 667},
        'Rrs_746': {'wavelenght': 746},
        'Rrs_868': {'wavelenght': 868},
        'Rrs_1238': {'wavelenght': 1238},
        'Rrs_1604': {'wavelenght': 1604},
        'Rrs_2258': {'wavelenght': 2258},
    }
    reflectance_bands_group = 'geophysical_data'
    nbands = len(reflectance_bands)
    return reflectance_bands, reflectance_bands_group, nbands


def get_flags_info(options):
    flag_band_name = 'l2_flags'
    flag_band_group = 'geophysical_data'
    return flag_band_name, flag_band_group


def get_find_product_info(options):
    path_source = options['file_path']['sat_source_dir']
    org = None
    if options.has_option('file_path', 'sat_source_dir_organization') and options['file_path'][
        'sat_source_dir_organization']:
        org = options['file_path']['sat_source_dir_organization']
    wce = '*'
    if options.has_option('file_path', 'wce') and options['file_path']['wce']:
        wce = options['file_path']['wce']
    time_start = dt.strptime(options['Time_and_sites_selection']['time_start'], '%Y-%m-%d')
    time_stop = dt.strptime(options['Time_and_sites_selection']['time_stop'], '%Y-%m-%d') + timedelta(hours=24)

    sensor = options['satellite_options']['sensor']
    platform = options['satellite_options']['platform']

    return path_source, org, wce, time_start, time_stop, sensor, platform


def check_product(filepath, sensor, platform, time_start, time_stop):
    nc_sat = Dataset(filepath, 'r')
    checkIns = False
    checkPlat = False
    checkTime = False
    if 'instrument' in nc_sat.ncattrs() and nc_sat.instrument == sensor:
        checkIns = True
    if 'platform' in nc_sat.ncattrs() and nc_sat.platform == platform:
        checkPlat = True
    if 'time_coverage_start' in nc_sat.ncattrs():
        time_sat = dt.strptime(nc_sat.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
        checkTime = True

    check_res = True
    if not checkIns:
        print(f'Attribute instrument is not available in the dataset, or not equal to: {sensor}')
        check_res = False
    if not checkPlat:
        fprint('Attribute platform is not available in the dataset, or not equal to: {platform}')
        check_res = False
    if not checkTime:
        fprint('Attribute time_coverage_start is not available in the dataset, or not equal to: {platform}')
        check_res = False
    else:
        if time_sat < time_start or time_sat > time_stop:
            fprint('Product is out of the temporal coverage')
            check_res = False

    return check_res


def create_list_products(path_source, org, wce):
    print('product list...')


def main():
    print('Creating satellite extracts.')

    if not args.config_file:
        return

    options = config_reader(args.config_file)

    # work only with the specified product file
    if args.product_file:
        filepath = args.product_file
        launch_create_extract(filepath, options)
    else:
        path_source, org, wce, time_start, time_stop, sensor, platform = get_find_product_info(options)
        print('aqui...')


if __name__ == '__main__':
    main()
