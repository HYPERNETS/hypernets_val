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
import shutil

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


def launch_create_extract_syke(filepath, syke_file, options):
    ncreated = 0

    path_output = get_output_path(options)
    if path_output is None:
        print(f'ERROR: {path_output} is not valid')
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

    sat_time = dt.strptime(nc_sat.start_time, '%Y-%m-%d %H:%M:%S')
    nminutes = 120
    if options.has_option('satellite_options', 'max_time_diff'):
        nminutes = int(options['satellite_options']['max_time_diff'])

    # Retrieving global atribbutes
    if args.verbose:
        print('[INFO] Retrieving global attributes...')
    global_at = get_global_atrib(nc_sat)

    max_time_diff = nminutes * 60
    check_date = syke_file.get_subdf(sat_time, sat_time, max_time_diff)
    if not check_date:
        print(f'[WARNING] Spectral data for date: {sat_time} are not avaiable')
        return ncreated
    if args.verbose:
        print(f'[INFO] Extracting data for sat time: {sat_time}')

    extracts = {}
    for irow in range(syke_file.get_n_df_sub()):
        insitu_lat, insitu_lon = syke_file.get_lat_lon_at_subdb(irow)
        size_box = get_box_size(options)
        contain_flag = check_location(insitu_lat, insitu_lon, lat, lon, size_box)
        if not contain_flag:
            continue
        r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
        # insitu_time = skie_file.get_time_at_subdb(irow, None)
        # distd, timed, speed = skie_file.get_dist_timedif_speed_at_subdb(irow)
        # index = skie_file.get_index_at_subdb(irow)
        # print(irow, index, insitu_time, insitu_lat, insitu_lon, r, c, distd, timed, speed)
        keyrc = f'{r}_{c}'
        if keyrc not in extracts.keys():
            extracts[keyrc] = {
                'r': r,
                'c': c,
                'irows': [irow]
            }
        else:
            extracts[keyrc]['irows'].append(irow)

    output_dir_prev = None
    if options.has_option('file_path', 'output_dir_prev'):
        output_dir_prev = options['file_path']['output_dir_prev']

    for extract in extracts:
        r = int(extracts[extract]['r'])
        c = int(extracts[extract]['c'])
        filename = filepath.split('/')[-1].replace('.', '_') + '_extract_skie_' + extract + '.nc'
        pdu = filepath.split('/')[-1]
        ofname = os.path.join(path_output, filename)

        if os.path.exists(ofname):
            if args.verbose:
                print(f'[INFO] {ofname} already exists. Skipping..')
            continue

        if output_dir_prev is not None:
            ofname_prev = os.path.join(output_dir_prev, filename)
            if os.path.exists(ofname_prev):
                if args.verbose:
                    print(f'[INFO] Copying {ofname} from previous directory')
                shutil.copy(ofname_prev, ofname)
                continue

        global_at['station_name'] = 'syke'
        global_at['in_situ_lat'] = f'in situ points at pixel row={r},col={c}'
        global_at['in_situ_lon'] = f'in situ points at pixel row={r},col={c}'

        b = create_extract(ofname, pdu, options, nc_sat, global_at, lat, lon, r, c, syke_file,
                           extracts[extract]['irows'])
        if b:
            ncreated = ncreated + 1

    nc_sat.close()

    return ncreated


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
            print('[INFO] Starting dataset...')
    nc_sat = Dataset(filepath, 'r')

    # Retriving lat and long arrays
    if args.verbose:
        print('[INFO] Retrieving lat/long data...')
    lat, lon = get_lat_long_arrays(nc_sat)

    # Retrieving global atribbutes
    if args.verbose:
        print('[INFO] Retrieving global attributes...')
    global_at = get_global_atrib(nc_sat)

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
            res = create_extract(ofname, pdu, options, nc_sat, global_at, lat, lon, r, c, None, None)
            if res:
                ncreated = ncreated + 1
                print(f'[INFO] Extract file created: {ofname}')
        else:
            if args.verbose:
                print(f'[WARNING] Site {site} out of the image')

    return ncreated


def create_extract_insitu(ofname, pdu, options, nc_sat, global_at, r, c, insitu_time, insitu_info, name_variables):
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
    flag_band_name = 'bitmask'

    newEXTRACT = SatExtract(ofname)
    if not newEXTRACT.FILE_CREATED:
        print(f'[ERROR] File {ofname} could not be created')
        return False

    if args.verbose:
        print(f'[INFO]Creating file: {ofname}')
    newEXTRACT.set_global_attributes(global_at)

    newEXTRACT.create_dimensions(size_box, n_bands)

    lat, long = get_lat_long_arrays(nc_sat)
    newEXTRACT.create_lat_long_variables(lat, long, window)

    # Sat time start:  ,+9-2021-12-24T18:23:00.471Z
    newEXTRACT.create_satellite_time_variable(dt.strptime(nc_sat.start_time, '%Y-%m-%d %H:%M:%S'))

    # pdu variable
    newEXTRACT.create_pdu_variable(pdu, global_at['sensor'])

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
    flag_band = nc_sat.variables[flag_band_name]
    desc = flag_band.description.split(',')
    flag_masks = ''
    flag_list = []
    flag_meanings = ''
    flag_started = False
    for d in desc:
        dval = d.split(':')
        flag_list.append(int(dval[1]))
        if not flag_started:
            flag_masks = dval[1]
            flag_meanings = dval[0]
            flag_started = True
        else:
            flag_masks = flag_masks + ' , ' + dval[1]
            flag_meanings = flag_meanings + ' ' + dval[0]
    newEXTRACT.create_flag_variable(f'satellite_{flag_band_name}', flag_band, 'Polymer Bitmask', flag_list,
                                    flag_meanings, window)

    ##add in situ variables
    formatdt = '%Y-%m-%d %H:%M:%S'
    sat_time = dt.strptime(nc_sat.start_time, formatdt)
    timediff = abs((insitu_time - sat_time).total_seconds())
    insitulat_var, insitulon_var, insitutime_var, time_difference_var = newEXTRACT.create_insitu_variables_for_single_insitu_data()
    insitulat_var[0] = global_at['in_situ_lat']
    insitulon_var[0] = global_at['in_situ_lon']
    insitutime_var[0] = insitu_time.timestamp()
    time_difference_var[0] = timediff
    for name_var in name_variables:
        insitu_var = newEXTRACT.create_insitu_variable_for_single_insitu_data(name_var, 'unknown', 'unknown')
        insitu_var[0] = insitu_info[name_var]

    newEXTRACT.close_file()

    return True


def create_extract(ofname, pdu, options, nc_sat, global_at, lat, long, r, c, syke_file, irows):
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
    flag_band_name = 'bitmask'

    newEXTRACT = SatExtract(ofname)
    if not newEXTRACT.FILE_CREATED:
        print(f'[ERROR] File {ofname} could not be created')
        return False

    if args.verbose:
        print(f'[INFO]Creating file: {ofname}')
    newEXTRACT.set_global_attributes(global_at)
    if syke_file is not None:
        newEXTRACT.create_dimensions_incluidinginsitu(size_box, n_bands, syke_file.get_n_bands(), 30)
    else:
        newEXTRACT.create_dimensions(size_box, n_bands)

    newEXTRACT.create_lat_long_variables(lat, long, window)

    # Sat time start:  ,+9-2021-12-24T18:23:00.471Z
    newEXTRACT.create_satellite_time_variable(dt.strptime(nc_sat.start_time, '%Y-%m-%d %H:%M:%S'))

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
    flag_band = nc_sat.variables[flag_band_name]
    desc = flag_band.description.split(',')
    flag_masks = ''
    flag_list = []
    flag_meanings = ''
    flag_started = False
    for d in desc:
        dval = d.split(':')
        flag_list.append(int(dval[1]))
        if not flag_started:
            flag_masks = dval[1]
            flag_meanings = dval[0]
            flag_started = True
        else:
            flag_masks = flag_masks + ' , ' + dval[1]
            flag_meanings = flag_meanings + ' ' + dval[0]

    newEXTRACT.create_flag_variable(f'satellite_{flag_band_name}', flag_band, 'Polymer Bitmask', flag_list,
                                    flag_meanings, window)

    if syke_file is not None:
        insitu_origbands_var = newEXTRACT.create_insitu_original_bands_variable()
        insitu_origbands_var[:] = syke_file.wavelengths
        insitutime_var = newEXTRACT.create_insitu_time_variable()
        insitu_exactwl_var = newEXTRACT.create_insitu_exact_wavelengths_variable()
        insitu_timedif_var = newEXTRACT.create_insitu_time_difference_variable()
        insitu_rrs_var = newEXTRACT.create_insitu_rrs_variable()
        insitu_lat, insitu_lon = newEXTRACT.create_insitu_lat_long_variables()
        index_var = 0
        for irow in irows:
            insitu_time_here = syke_file.get_time_at_subdb(irow, 'DT')
            insitutime_var[0, index_var] = insitu_time_here.timestamp()
            insitu_exactwl_var[0, :, index_var] = syke_file.wavelengths
            insitu_rrs_var[0, :, index_var] = np.array(syke_file.get_spectra_at_subdb(irow))
            sat_time_here = dt.strptime(nc_sat.start_time, '%Y-%m-%d %H:%M:%S')
            time_diff = float(abs((insitu_time_here - sat_time_here).total_seconds()))
            insitu_timedif_var[0, index_var] = time_diff
            latp, lonp = syke_file.get_lat_lon_at_subdb(irow)
            insitu_lat[0, index_var] = latp
            insitu_lon[0, index_var] = lonp
            index_var = index_var + 1

    newEXTRACT.close_file()

    return True


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
    lat = nc_sat.variables['latitude'][:, :]
    lon = nc_sat.variables['longitude'][:, :]
    return lat, lon


def get_global_atrib(nc_sat):
    at = {'sensor': nc_sat.sensor}
    if nc_sat.sensor == 'OLCI':
        at['satellite'] = 'S3'
        at['platform'] = ''
        path_origin = nc_sat.l2_filename
        if path_origin.find('S3A') > 0:
            at['platform'] = 'A'
        elif path_origin.find('S3B') > 0:
            at['platform'] = 'B'

    at['res'] = ''
    at['aco_processor'] = 'Polymer'
    at['proc_version'] = ''
    polymer_path = nc_sat.dir_base
    ipol = polymer_path.find('polymer-v')
    if ipol > 0:
        at['proc_version'] = polymer_path[ipol + 10:len(polymer_path)]

    at['station_name'] = ''
    at['in_situ_lat'] = -999
    at['in_situ_lon'] = -999

    return at


def get_reflectance_bands_info(nc_sat):
    reflectance_bands = None
    nbands = 0

    if nc_sat.bands_rw:
        reflectance_bands = {}
        bands_here = nc_sat.bands_rw[1:len(nc_sat.bands_rw) - 1].split(',')
        wl_list = nc_sat.central_wavelength[1:len(nc_sat.central_wavelength) - 1].split(',')
        wl_dict = {}
        for wl in wl_list:
            wls = wl.split(':')
            wl_dict[wls[0].strip()] = float(wls[1].strip())
        for band in bands_here:
            name_band = 'Rw' + band.strip()
            wl_band = wl_dict[band.strip()]
            reflectance_bands[name_band] = {'wavelenght': wl_band}
            nbands = nbands + 1

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
    try:
        nc_sat = Dataset(filepath, 'r')
    except OSError:
        if args.verbose:
            print(f'[WARNING] {filepath} is not a valid dataset. Skipping...')
        return False

    checkTime = False

    if 'start_time' in nc_sat.ncattrs():
        time_sat = dt.strptime(nc_sat.start_time, '%Y-%m-%d %H:%M:%S')
        checkTime = True

    check_res = True
    if not checkTime and check_res:
        if args.verbose:
            print(f'[WARNING] Attribute time_coverage_start is not available in the dataset. Skipping...')
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
            if os.path.isdir(path_to_sat_source):
                continue

            STATUS_STR = 'NO PASSED'
            if check_product(path_to_sat_source, time_start, time_stop):
                product_path_list.append(path_to_sat_source)
                STATUS_STR = 'OK'

            if args.verbose:
                # print('-----------------')
                pname = path_to_sat_source.split('/')[-1]
                print(f'[INFO ]Check: {pname}. STATUS: {STATUS_STR}')

    return product_path_list


def run_syke_option(options):
    if args.verbose:
        print(f'[INFO] STARTING SYKE MODE')
    path_syke = options['file_path']['path_syke']
    path_syke_code = options['file_path']['path_syke_code']
    if not os.path.exists(path_syke):
        print(f'[ERROR] Syke file: {path_syke} is not avialable')
        return
    if not os.path.exists(path_syke_code) or not os.path.isdir(path_syke_code):
        print(f'[ERROR] Syke code folder: {path_syke_code} is not avialable')
        return
    if args.verbose:
        print(f'[INFO] Path syke: {path_syke}')
        print(f'[INFO] Path syke: {path_syke_code}')
    sys.path.append(path_syke_code)
    try:
        from skie_csv import SKIE_CSV
    except ModuleNotFoundError:
        print(f'[ERROR] Syke module is not found')
        return
    syke_file = SKIE_CSV(path_syke)
    syke_file.start_list_dates()
    syke_file.extract_wl_colnames()

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
            nhere = launch_create_extract_syke(filepath, syke_file, options)
            ncreated = ncreated + nhere
        print('------------------------------')
        print(f'COMPLETED. {ncreated} sat extract files were created')


def run_insitu_option(options):
    path_insitu = options['file_path']['path_insitu']
    path_insitu_code = options['file_path']['path_insitu_code']
    if not os.path.exists(path_insitu):
        print(f'[ERROR] In situ file: {path_insitu} is not avialable')
        return
    if not os.path.exists(path_insitu_code) or not os.path.isdir(path_insitu_code):
        print(f'[ERROR] In situ code folder: {path_insitu_code} is not avialable')
        return
    if args.verbose:
        print(f'[INFO] Path in situ: {path_insitu}')
        print(f'[INFO] Path code: {path_insitu_code}')
    sys.path.append(path_insitu_code)
    try:
        from insitu.nc_insitu_file import NCInsituFile
        from insitu.csv_insitu_file import CSVInsituFile
    except ModuleNotFoundError:
        print(f'[ERROR] trims3images_frominsitunc module is not found')
        return

    if path_insitu.endswith('nc'):
        ifobj = NCInsituFile(path_insitu)
    if path_insitu.endswith('csv'):
        ifobj = CSVInsituFile(path_insitu)
    path_source, org, wce, time_start, time_stop = get_find_product_info(options)
    ifobj.get_valid_dates(None, time_start, time_stop)
    valid_dates = ifobj.valid_dates
    name_variables = ifobj.variables

    path_to_list = get_path_list_products(options)
    if path_to_list is None:
        return

    product_path_list = get_list_products(path_to_list, path_source, org, wce, time_start, time_stop)
    if len(product_path_list) == 0:
        print(f'[WARNING] No valid datasets found on:  {path_source}')
        print('[INFO] COMPLETED. 0 sat extract files were created')
        return
    path_output = get_output_path(options)
    if path_output is None:
        print(f'ERROR: {path_output} is not valid')
        return

    nminutes = 120
    if options.has_option('satellite_options', 'max_time_diff'):
        nminutes = int(options['satellite_options']['max_time_diff'])
    if args.verbose:
        print(f'[INFO] Maximum time difference betwee sat and in situ times: {nminutes} minutes')

    formatdt = '%Y-%m-%d %H:%M:%S'
    maxtimediff = nminutes * 60
    size_box = get_box_size(options)
    ncreated = 0
    for filepath in product_path_list:
        if args.verbose:
            print('------------------------------')
            pname = filepath.split('/')[-1]
            print(f'[INFO]Checking sat extracts for product: {pname}')
        nc_sat = Dataset(filepath, 'r')
        lat, lon = get_lat_long_arrays(nc_sat)
        sat_time = dt.strptime(nc_sat.start_time, formatdt)
        print(sat_time)
        sat_date_str = sat_time.strftime('%Y-%m-%d')

        if sat_date_str in valid_dates.keys():
            for h in valid_dates[sat_date_str]:
                insitu_time_str = f'{sat_date_str} {h}'
                insitu_time = dt.strptime(insitu_time_str, formatdt)
                insitu_lat = valid_dates[sat_date_str][h]['lat']
                insitu_lon = valid_dates[sat_date_str][h]['lon']
                timediff = abs((insitu_time - sat_time).total_seconds())
                contain_flag = check_location(insitu_lat, insitu_lon, lat, lon, size_box)
                if timediff < maxtimediff and contain_flag:
                    if args.verbose:
                        print(
                            f'[INFO] Preparing extract. In situ latitude: {insitu_lat} Longitude: {insitu_lon} Time: {insitu_time_str}')
                    r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
                    extract = f'{r}_{c}'
                    filename = filepath.split('/')[-1].replace('.', '_') + '_extract_insitu_' + extract + '.nc'
                    pdu = filepath.split('/')[-1]
                    ofname = os.path.join(path_output, filename)
                    global_at = get_global_atrib(nc_sat)
                    global_at['station_name'] = 'in situ dataset'
                    global_at['in_situ_lat'] = insitu_lat
                    global_at['in_situ_lon'] = insitu_lon
                    insitu_info = valid_dates[sat_date_str][h]
                    created = create_extract_insitu(ofname, pdu, options, nc_sat, global_at, r, c, insitu_time,
                                                    insitu_info, name_variables)
                    if created:
                        ncreated = ncreated + 1
        nc_sat.close()
        # nhere = 1  # nhere = launch_create_extract_syke(filepath, syke_file, options)
        # ncreated = ncreated + nhere
    print('------------------------------')
    print(f'COMPLETED. {ncreated} sat extract files were created')


def main():
    print('[INFO]Creating satellite extracts')

    if not args.config_file:
        return

    options = config_reader(args.config_file)

    if options.has_option('file_path', 'path_syke') and options.has_option('file_path', 'path_syke_code'):
        run_syke_option(options)
        return

    if options.has_option('file_path', 'path_insitu') and options.has_option('file_path', 'path_insitu_code'):
        run_insitu_option(options)
        return

    # work only with the specified product file
    if args.product_file:
        print('------------------------------')
        filepath = args.product_file
        if args.verbose:
            print(f'[INFO]Extracting sat extract for product: {filepath}')
        ncreated = launch_create_extract(filepath, options, None, None)
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
