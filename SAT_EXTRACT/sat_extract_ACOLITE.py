import argparse
import configparser
import os
import sys
from netCDF4 import Dataset
from datetime import datetime as dt
from datetime import timedelta
import numpy.ma as ma
import numpy as np
import subprocess

# import user defined functions from other .py
code_home = os.path.abspath('../')
sys.path.append(code_home)
from SAT_EXTRACT.sat_extract import SatExtract
import COMMON.common_functions as cfs

parser = argparse.ArgumentParser(description="Create Match-up DataBase files (MDB) files.")
parser.add_argument("-d", "--debug", help="Debugging mode.", action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-p', "--product_file", help="Image file.")

args = parser.parse_args()


def launch_create_extract(filepath, options):
    if os.path.isdir((filepath)):
        pathbase = filepath
        for fname in os.listdir(pathbase):
            if fname.find('L2R.nc') > 0:
                filepath = os.path.join(pathbase, fname)

    ncreated = 0

    # Retrieving sites
    path_output = get_output_path(options)
    if path_output is None:
        print(f'ERROR: {path_output} is not valid')
        return ncreated

    site_file, site_list, region_list = get_site_options(options)
    if site_file is not None:
        sites = cfs.get_sites_from_file(site_file, site_list, region_list, path_output)
    else:
        sites = cfs.get_sites_from_list(site_list, path_output)

    if len(sites) == 0:
        print('[ERROR] No sites are defined')
        return ncreated

    # Start dataset
    if args.verbose:
        if args.verbose:
            print('[INFO] Starting dataset...')
    nc_sat = Dataset(filepath, 'r')

    # Retriving lat and long arrays
    if args.verbose:
        print('[INFO] Retrieving lat/long data...')
    lat, lon = get_lat_long_arrays(nc_sat)

    # Retrieving global atribbutes
    if args.verbose:
        print('[INFO] Retrieving global attributes...')
    global_at = get_global_atrib(filepath)

    # Working for each site, checking if there is in the image
    for site in sites:
        if args.verbose:
            print(f'[INFO]Working for site: {site}')
        insitu_lat = sites[site]['latitude']
        insitu_lon = sites[site]['longitude']
        contain_flag = 0
        if cfs.contain_location(lat, lon, insitu_lat, insitu_lon) == 1:
            r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
            size_box = get_box_size(options)
            start_idx_x = (c - int(size_box / 2))
            stop_idx_x = (c + int(size_box / 2) + 1)
            start_idx_y = (r - int(size_box / 2))
            stop_idx_y = (r + int(size_box / 2) + 1)
            if start_idx_y >= 0 and (stop_idx_y + 1) < lat.shape[0] and start_idx_x >= 0 and (stop_idx_x + 1) < \
                    lat.shape[1]:
                contain_flag = 1
        if contain_flag == 1:
            filename = filepath.split('/')[-1].replace('.', '_') + '_extract_' + site + '.nc'
            pdu = filepath.split('/')[-1]
            if not os.path.exists(sites[site]['path_out']):
                os.mkdir(sites[site]['path_out'])
            ofname = os.path.join(sites[site]['path_out'], filename)
            global_at['station_name'] = site
            global_at['in_situ_lat'] = insitu_lat
            global_at['in_situ_lon'] = insitu_lon
            res = create_extract(ofname, pdu, options, nc_sat, global_at, lat, lon, r, c)
            if res:
                ncreated = ncreated + 1
                print(f'[INFO] Extract file created: {ofname}')
        else:
            if args.verbose:
                print(f'[WARNING] Site {site} out of the image')

    return ncreated


def create_extract(ofname, pdu, options, nc_sat, global_at, lat, long, r, c):
    size_box = get_box_size(options)
    start_idx_x = (c - int(size_box / 2))
    stop_idx_x = (c + int(size_box / 2) + 1)
    start_idx_y = (r - int(size_box / 2))
    stop_idx_y = (r + int(size_box / 2) + 1)
    window = [start_idx_y, stop_idx_y, start_idx_x, stop_idx_x]

    reflectance_bands, n_bands = get_reflectance_bands_info(nc_sat)
    if n_bands == 0:
        print('[ERROR] reflectance bands are not defined')
        return False
    # flag_band_name = 'c2rcc_flags'
    # if args.atm_correction == 'FUB':
    #     flag_band_name = 'quality_flags'

    newEXTRACT = SatExtract(ofname)
    if not newEXTRACT.FILE_CREATED:
        print(f'[ERROR] File {ofname} could not be created')
        return False

    if args.verbose:
        print(f'[INFO]Creating file: {ofname}')
    newEXTRACT.set_global_attributes(global_at)
    newEXTRACT.create_dimensions(size_box, n_bands)

    newEXTRACT.create_lat_long_variables(lat, long, window)

    # Sat time start:  ,+9-2021-12-24T18:23:00.471Z
    # if 'isodate' in nc_sat.ncattrs():
    isodate = dt.strptime(nc_sat.isodate, '%Y-%m-%dT%H:%M:%S.%f+00:00')
    print(isodate)
    newEXTRACT.create_satellite_time_variable(isodate)

    # else:
    #     sat_time = get_sat_time_from_fname(pdu)
    #     if sat_time is not None:
    #         newEXTRACT.create_satellite_time_variable(sat_time)
    #     else:
    #         print(f'[ERROR] Satellite time is not defined...')
    #         newEXTRACT.close_file()
    #         return False

    # pdu variable
    newEXTRACT.create_pdu_variable(pdu, global_at['sensor'])

    # Rrs and wavelenghts
    satellite_Rrs = newEXTRACT.create_rrs_variable(global_at['sensor'])
    rbands = list(reflectance_bands.keys())
    wavelenghts = []
    for index in range(len(rbands)):
        rband = rbands[index]
        bandarray = ma.array(nc_sat.variables[rband][:, :])
        satellite_Rrs[0, index, :, :] = bandarray[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x] / np.pi
        wl = reflectance_bands[rband]['wavelenght']
        wavelenghts.append(wl)
    newEXTRACT.create_satellite_bands_variable(wavelenghts)

    # flags
    # flag_band = nc_sat.variables[flag_band_name]
    #
    # newEXTRACT.create_flag_variable(f'satellite_{flag_band_name}', flag_band, flag_band.long_name, flag_band.flag_masks,
    #                                 flag_band.flag_meanings, window)

    newEXTRACT.close_file()

    return True


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
    if options.has_option('file_path', 'output_dir'):
        output_path = options['file_path']['output_dir']
        if not os.path.exists(output_path):
            if args.verbose:
                print(f'[WARNING] Creating new output path: {output_path}')
                try:
                    os.mkdir(output_path)
                except OSError:
                    print(f'[ERROR] Folder: {output_path} could not be creaded')
                    return None
        if args.verbose:
            print(f'[INFO] Output path:{output_path}')
        return output_path
    else:
        print('[ERROR] section:file_path, option: output_dir is not included in the configuration file')
        return None


def get_box_size(options):
    if options.has_option('satellite_options', 'extract_size'):
        try:
            size_box = int(options['satellite_options']['extract_size'])
        except:
            size_box = 25
    else:
        size_box = 25
    return size_box


def get_site_options(options):
    site_file = None
    site_list = None
    region_list = None
    if not options.has_option('Time_and_sites_selection', 'sites'):
        print(f'[ERROR] section Time_and_sites_selection, option sites not found in configuration file')
    if options.has_option('Time_and_sites_selection', 'sites_file') and \
            options['Time_and_sites_selection']['sites_file']:
        site_file = options['Time_and_sites_selection']['sites_file']
    if options['Time_and_sites_selection']['sites']:
        site_list = options['Time_and_sites_selection']['sites'].split(',')
        site_list = [s.strip() for s in site_list]
    if not options['Time_and_sites_selection']['sites'] and \
            options.has_option('Time_and_sites_selection', 'sites_region') and \
            options['Time_and_sites_selection']['sites_region']:
        region_list = options['Time_and_sites_selection']['sites_region'].split(',')
        region_list = [r.strip() for r in region_list]
    return site_file, site_list, region_list


def get_lat_long_arrays(nc_sat):
    lat = nc_sat.variables['lat'][:, :]
    lon = nc_sat.variables['lon'][:, :]
    return lat, lon


def get_sat_time_from_fname(fname):
    val_list = fname.split('_')
    sat_time = None
    for v in val_list:
        try:
            sat_time = dt.strptime(v, '%Y%m%dT%H%M%S')
            break
        except ValueError:
            continue
    return sat_time


def get_global_atrib(file_path):
    if file_path.find('S3A') > 0:
        at = {'sensor': 'OLCI', 'satellite': 'S3', 'platform': 'A'}
    elif file_path.find('S3B') > 0:
        at = {'sensor': 'OLCI', 'satellite': 'S3', 'platform': 'B'}
    else:
        at = {'sensor': 'unknown', 'satellite': 'unknown', 'platform': 'unknown'}

    at['res'] = ''
    at['aco_processor'] = 'ACOLITE'
    at['proc_version'] = '20220222.0'

    at['station_name'] = ''
    at['in_situ_lat'] = -999
    at['in_situ_lon'] = -999

    return at


def get_reflectance_bands_info(nc_sat):
    search_pattern = 'rhos_'

    reflectance_bands = {}
    nbands = 0
    imin = 100000
    imax = 0
    for var in nc_sat.variables:
        if var.startswith(search_pattern):
            ival = int(var.split('_')[1].strip())
            if ival < imin:
                imin = ival
            if ival > imax:
                imax = ival

    for ival in range(imin, imax + 1):
        var_name = search_pattern + str(ival)
        if var_name in nc_sat.variables:
            wl_band = nc_sat.variables[var_name].wave_nm
            reflectance_bands[var_name] = {'wavelenght': wl_band}
            nbands = nbands + 1

    if nbands == 0:
        reflectance_bands = None

    return reflectance_bands, nbands


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

    return path_source, org, wce, time_start, time_stop


def check_product(filepath, time_start, time_stop):

    if os.path.isdir((filepath)):
        pathbase = filepath
        for fname in os.listdir(pathbase):
            if fname.find('L2R.nc') > 0:
                filepath = os.path.join(pathbase, fname)
    else:
        return False

    if not filepath.endswith('.nc'):
        return False

    try:
        nc_sat = Dataset(filepath, 'r')
    except OSError:
        if args.verbose:
            print(f'[WARNING] {filepath} is not a valid dataset. Skipping...')
        return False

    checkTime = False
    if 'isodate' in nc_sat.ncattrs():
        try:
            time_sat = dt.strptime(nc_sat.isodate, '%Y-%m-%dT%H:%M:%S.%f+00:00')
            checkTime = True
        except ValueError:
            checkTime = False

    check_res = True
    if not checkTime and check_res:
        if args.verbose:
            print(f'[WARNING] Attribute start_date is not available in the dataset. Skipping...')
        check_res = False
    if check_res:
        if time_sat < time_start or time_sat > time_stop:
            if args.verbose:
                print('[WARNING] Product is out of the temporal coverage. Skipping...')
            check_res = False

    return check_res


def get_path_list_products(options):
    path_out = get_output_path(options)
    if path_out is None:
        return None
    # path_to_list = f'{path_out}/file_{platform}_list.txt'
    path_to_list = f'{path_out}/file_list.txt'
    return path_to_list


def get_list_products(path_to_list, path_source, org, wce, time_start, time_stop):
    if args.verbose:
        print('[INFO]Creating list of products')

    if os.path.exists(path_to_list):
        if args.verbose:
            print('[INFO]Deleting previous file list...')
        cmd = f'rm {path_to_list}'
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]{err}')

    if org is None:
        if wce == '*':
            cmd = f'find {path_source} |sort|uniq> {path_to_list}'
        else:
            cmd = f'find {path_source} -name {wce}|sort|uniq> {path_to_list}'
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]{err}')

    product_path_list = []
    with open(path_to_list, 'r') as file:
        for cnt, line in enumerate(file):
            path_to_sat_source = line[:-1]
            if not os.path.isdir(path_to_sat_source):
                continue
            if args.verbose:
                print('-----------------')
                print(f'Checking: {path_to_sat_source}')
            if check_product(path_to_sat_source, time_start, time_stop):
                product_path_list.append(path_to_sat_source)

    return product_path_list


def main():
    print('[INFO]Creating satellite extracts')

    if not args.config_file:
        return

    options = config_reader(args.config_file)

    # work only with the specified product file
    if args.product_file:
        print('------------------------------')
        filepath = args.product_file
        if args.verbose:
            print(f'[INFO]Extracting sat extract for product: {filepath}')
        ncreated = launch_create_extract(filepath, options)
        print('------------------------------')
        print(f'[INFO] COMPLETED. {ncreated} sat extract files were created')
    else:
        print('------------------------------')
        path_source, org, wce, time_start, time_stop = get_find_product_info(options)
        path_to_list = get_path_list_products(options)
        if path_to_list is not None:
            product_path_list = get_list_products(path_to_list, path_source, org, wce, time_start, time_stop)
            if len(product_path_list) == 0:
                if args.verbose:
                    print('-----------------')
                print(f'[WARNING] No valid datasets found on:  {path_source}')
                print('------------------------------')
                print('[INFO] COMPLETED. 0 sat extract files were created')
                return
            ncreated = 0
            for filepath in product_path_list:
                if args.verbose:
                    print('------------------------------')
                    print(f'[INFO]Extracting sat extract for product: {filepath}')
                nhere = launch_create_extract(filepath, options)
                ncreated = ncreated + nhere
            print('------------------------------')
            print(f'COMPLETED. {ncreated} sat extract files were created')


if __name__ == '__main__':
    main()
