import argparse
import configparser
import os
import sys
import subprocess
from netCDF4 import Dataset
from datetime import datetime as dt
from datetime import timedelta
import numpy as np

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
            print(f'[INFO] Starting dataset: {filepath}')
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
            if os.path.exists(ofname):
                ncreated = ncreated + 1
                print(f'[INFO] Extract file already created: {ofname}. Skipping...')
                nc_sat.close()
                continue
            global_at['station_name'] = site
            global_at['in_situ_lat'] = insitu_lat
            global_at['in_situ_lon'] = insitu_lon
            res = create_extract(ofname, pdu, options, nc_sat, global_at, lat, lon, r, c, None, None)
            if res:
                ncreated = ncreated + 1
                print(f'[INFO] Extract file created: {ofname}')
        else:
            if args.verbose:
                print(f'[WARNING] Site {site} out of the image')

    return ncreated


def create_extract(ofname, pdu, options, nc_sat, global_at, lat, long, r, c, skie_file, irows):
    size_box = get_box_size(options)
    start_idx_x = (c - int(size_box / 2))
    stop_idx_x = (c + int(size_box / 2) + 1)
    start_idx_y = (r - int(size_box / 2))
    stop_idx_y = (r + int(size_box / 2) + 1)
    window = [start_idx_y, stop_idx_y, start_idx_x, stop_idx_x]

    newEXTRACT = SatExtract(ofname)
    if not newEXTRACT.FILE_CREATED:
        print(f'[ERROR] File {ofname} could not be created')
        return False

    if args.verbose:
        print(f'[INFO] Creating file: {ofname}')
    newEXTRACT.set_global_attributes(global_at)
    if skie_file is not None:
        newEXTRACT.create_dimensions_basic_includinginsitu(size_box,skie_file.get_n_bands(), 30)
        # newEXTRACT.create_dimensions_incluidinginsitu(size_box, n_bands, skie_file.get_n_bands(), 30)
    else:
        newEXTRACT.create_dimensions_basic(size_box)
        # newEXTRACT.create_dimensions(size_box, n_bands)

    newEXTRACT.create_lat_long_variables(lat, long, window)

    # Sat time start:  ,+9-2021-12-24T18:23:00.471Z
    if 'start_date' in nc_sat.ncattrs():
        newEXTRACT.create_satellite_time_variable(dt.strptime(nc_sat.start_date, '%d-%b-%Y %H:%M:%S.%f'))
    else:
        sat_time = get_sat_time_from_fname(pdu)
        if sat_time is not None:
            newEXTRACT.create_satellite_time_variable(sat_time)
        else:
            print(f'[ERROR] Satellite time is not defined...')
            newEXTRACT.close_file()
            return False

    # pdu variable
    newEXTRACT.create_pdu_variable(pdu, global_at['sensor'])

    # 2D variables
    for var in nc_sat.variables:
        ilat = var.find('lat')
        ilon = var.find('lon')
        if ilat >= 0 or ilon >= 0:
            continue
        band = nc_sat.variables[var]
        if 'flag_meanings' in band.ncattrs() and 'flag_masks' in band.ncattrs():
            if 'long_name' in band.ncattrs():
                long_name = band.long_name
            else:
                long_name = f'Flags for {var}'
            newEXTRACT.create_flag_variable(f'satellite_{var}', band, long_name, band.flag_masks,
                                            band.flag_meanings, window)
        else:
            newEXTRACT.create_2D_variable_general(var, band, window)


    if skie_file is not None:
        insitu_origbands_var = newEXTRACT.create_insitu_original_bands_variable()
        insitu_origbands_var[:] = skie_file.wavelengths
        insitutime_var = newEXTRACT.create_insitu_time_variable()
        insitu_exactwl_var = newEXTRACT.create_insitu_exact_wavelengths_variable()
        insitu_timedif_var = newEXTRACT.create_insitu_time_difference_variable()
        insitu_rrs_var = newEXTRACT.create_insitu_rrs_variable()
        insitu_lat, insitu_lon = newEXTRACT.create_insitu_lat_long_variables()
        index_var = 0
        for irow in irows:
            insitu_time_here = skie_file.get_time_at_subdb(irow, 'DT')
            insitutime_var[0, index_var] = insitu_time_here.timestamp()
            insitu_exactwl_var[0, :, index_var] = skie_file.wavelengths
            insitu_rrs_var[0, :, index_var] = np.array(skie_file.get_spectra_at_subdb(irow))

            sat_time_here = dt.strptime(nc_sat.start_date, '%d-%b-%Y %H:%M:%S.%f')
            time_diff = float(abs((insitu_time_here - sat_time_here).total_seconds()))
            insitu_timedif_var[0, index_var] = time_diff
            latp, lonp = skie_file.get_lat_lon_at_subdb(irow)
            insitu_lat[0, index_var] = latp
            insitu_lon[0, index_var] = lonp
            index_var = index_var + 1

    newEXTRACT.close_file()

    return True


def launch_create_extract_syke(filepath, skie_file, options):

    ncreated = 0

    path_output = get_output_path(options)
    if path_output is None:
        print(f'ERROR: {path_output} is not valid')
        return ncreated

    #Start dataset
    if args.verbose:
        if args.verbose:
            print('[INFO] Starting dataset...')
    nc_sat = Dataset(filepath, 'r')

    # Retriving lat and long arrays
    if args.verbose:
        print('[INFO] Retrieving lat/long data...')
    lat, lon = get_lat_long_arrays(nc_sat)

    sat_time = get_sat_time_from_fname(filepath)
    nminutes = 120
    if options.has_option('satellite_options', 'max_time_diff'):
        nminutes = int(options['satellite_options']['max_time_diff'])

    if args.verbose:
        print(f'[INFO] Maximum time diferrence set to: {nminutes}')


    # Retrieving global atribbutes
    if args.verbose:
        print('[INFO] Retrieving global attributes...')
    global_at = get_global_atrib(filepath)

    max_time_diff = nminutes * 60
    check_date = skie_file.get_subdf(sat_time, sat_time, max_time_diff)
    if not check_date:
        print(f'[WARNING] Spectral data for date: {sat_time} are not avaiable')
        return ncreated
    if args.verbose:
        print(f'[INFO] Extracting data for sat time: {sat_time}')

    extracts = {}
    for irow in range(skie_file.get_n_df_sub()):
        insitu_lat, insitu_lon = skie_file.get_lat_lon_at_subdb(irow)
        size_box = get_box_size(options)
        contain_flag = check_location(insitu_lat, insitu_lon, lat, lon, size_box)
        if not contain_flag:
            continue
        r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
        keyrc = f'{r}_{c}'
        if keyrc not in extracts.keys():
            if args.verbose:
                print(f'[INFO]   --> In situ spectra added for row: {r} and column: {c}')
            extracts[keyrc] = {
                'r': r,
                'c': c,
                'irows': [irow]
            }
        else:
            extracts[keyrc]['irows'].append(irow)

    for extract in extracts:
        r = int(extracts[extract]['r'])
        c = int(extracts[extract]['c'])
        if args.verbose:
            print(f'[INFO] Creating extract for row: {r} and column {c}')
        filename = filepath.split('/')[-1].replace('.', '_') + '_extract_skie_' + extract + '.nc'
        pdu = filepath.split('/')[-1]
        ofname = os.path.join(path_output, filename)
        global_at['station_name'] = 'skie'
        global_at['in_situ_lat'] = f'in situ points at pixel row={r},col={c}'
        global_at['in_situ_lon'] = f'in situ points at pixel row={r},col={c}'

        b = create_extract(ofname, pdu, options, nc_sat, global_at, lat, lon, r, c, skie_file,
                            extracts[extract]['irows'])

        if b:
            ncreated = ncreated + 1

    nc_sat.close()

    return ncreated




def check_location(insitu_lat, insitu_lon, lat, lon, size_box):
    contain_flag = False
    if cfs.contain_location(lat, lon, insitu_lat, insitu_lon) == 1:
        r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
        start_idx_x = (c - int(size_box / 2))
        stop_idx_x = (c + int(size_box / 2) + 1)
        start_idx_y = (r - int(size_box / 2))
        stop_idx_y = (r + int(size_box / 2) + 1)
        if start_idx_y >= 0 and (stop_idx_y + 1) < lat.shape[0] and start_idx_x >= 0 and (stop_idx_x + 1) < \
                lat.shape[1]:
            contain_flag = True
    return contain_flag
def get_global_atrib(file_path):
    if file_path.find('S3A') > 0:
        at = {'sensor': 'OLCI', 'satellite': 'S3', 'platform': 'A'}
    elif file_path.find('S3B') > 0:
        at = {'sensor': 'OLCI', 'satellite': 'S3', 'platform': 'B'}
    else:
        at = {'sensor': 'unknown', 'satellite': 'unknown', 'platform': 'unknown'}

    at['res'] = ''
    at['aco_processor'] = 'IDEPIX'
    at['proc_version'] = ''

    at['station_name'] = ''
    at['in_situ_lat'] = -999
    at['in_situ_lon'] = -999

    return at





def get_box_size(options):
    if options.has_option('satellite_options', 'extract_size'):
        try:
            size_box = int(options['satellite_options']['extract_size'])
        except:
            size_box = 25
    else:
        size_box = 25
    return size_box


def get_lat_long_arrays(nc_sat):
    lat = nc_sat.variables['lat'][:, :]
    lon = nc_sat.variables['lon'][:, :]
    return lat, lon


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

def config_reader(FILEconfig):
    options = configparser.ConfigParser()
    options.read(FILEconfig)
    return options


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
            if os.path.isdir(path_to_sat_source):
                continue
            if args.verbose:
                print('-----------------')
                print(f'Checking: {path_to_sat_source}')
            if check_product(path_to_sat_source, time_start, time_stop):
                product_path_list.append(path_to_sat_source)

    return product_path_list


def check_product(filepath, time_start, time_stop):
    try:
        nc_sat = Dataset(filepath, 'r')
    except OSError:
        if args.verbose:
            print(f'[WARNING] {filepath} is not a valid dataset. Skipping...')
        return False

    checkTime = False
    if 'start_date' in nc_sat.ncattrs():
        time_sat = dt.strptime(nc_sat.start_date, '%d-%b-%Y %H:%M:%S.%f')
        checkTime = True
    else:
        time_sat = get_sat_time_from_fname(filepath)
        if time_sat is not None:
            checkTime = True

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


def main():
    print('[INFO]Creating satellite extracts')

    if not args.config_file:
        return

    options = config_reader(args.config_file)

    # SYKE, NO IMPLEMENTED
    if options.has_option('file_path', 'path_skie') and options.has_option('file_path', 'path_skie_code'):
        path_skie = options['file_path']['path_skie']
        path_skie_code = options['file_path']['path_skie_code']
        if not os.path.exists(path_skie):
            print(f'[ERROR] Skie file: {path_skie} is not avialable')
            return
        if not os.path.exists(path_skie_code) or not os.path.isdir(path_skie_code):
            print(f'[ERROR] Skie code folder: {path_skie_code} is not avialable')
            return
        if args.verbose:
            print(f'[INFO] Path skie: {path_skie}')
            print(f'[INFO] Path skie: {path_skie_code}')
        sys.path.append(path_skie_code)
        try:
            from skie_csv import SKIE_CSV
        except ModuleNotFoundError:
            print(f'[ERROR] Skie module is not found')
            return
        skie_file = SKIE_CSV(path_skie)
        skie_file.start_list_dates()
        skie_file.extract_wl_colnames()

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
                    pname = filepath.split('/')[-1]
                    print(f'[INFO]Extracting sat extract for product: {pname}')
                nhere = launch_create_extract_syke(filepath, skie_file, options)
                ncreated = ncreated + nhere
            print('------------------------------')
            print(f'COMPLETED. {ncreated} sat extract files were created')
        return

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
