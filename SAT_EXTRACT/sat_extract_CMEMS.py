import argparse
import configparser
import os
import sys

import pandas as pd
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
parser.add_argument('-ac', "--atm_correction", help="Atmospheric correction algorithm (Default: C2RCC)",
                    choices=["C2RCC", "FUB"])

args = parser.parse_args()


# no implemented yet
def launch_create_extract_skie(filepath, skie_file, options):
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
    var_lat, var_lon = get_lat_long_var_names(options)
    lat, lon = get_lat_long_arrays(nc_sat, var_lat, var_lon)
    if args.atm_correction == 'FUB':
        sat_time = get_sat_time_from_fname(filepath)
    else:
        sat_time = dt.strptime(nc_sat.start_date, '%Y-%m-%d')
    nminutes = 120
    if options.has_option('satellite_options', 'max_time_diff'):
        nminutes = int(options['satellite_options']['max_time_diff'])

    # Retrieving global atribbutes
    if args.verbose:
        print('[INFO] Retrieving global attributes...')
    global_at = get_global_atrib(nc_sat, options)

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
        contain_flag, r, c = check_location(insitu_lat, insitu_lon, lat, lon, size_box)
        if not contain_flag:
            continue
        # r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
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

    for extract in extracts:
        r = int(extracts[extract]['r'])
        c = int(extracts[extract]['c'])
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


def launch_create_extract_station(filepath, options, insitu_lat, insitu_lon):
    created = False
    path_output = get_output_path(options)
    if path_output is None:
        print(f'ERROR: {path_output} is not valid')
        return created

    # Start dataset
    if args.verbose:
        print(f'[INFO] Starting dataset. File: {filepath}')
    nc_sat = Dataset(filepath, 'r')

    # Retriving lat and long arrays
    if args.verbose:
        print('[INFO] Retrieving lat/long data...')
    var_lat, var_lon = get_lat_long_var_names(options)
    lat, lon = get_lat_long_arrays(nc_sat, var_lat, var_lon)

    # Retrieving global atribbutes
    if args.verbose:
        print('[INFO] Retrieving global attributes...')
    global_at = get_global_atrib(nc_sat, options)

    contain_flag = 0
    if cfs.contain_location(lat, lon, insitu_lat, insitu_lon) == 1:
        if lat.ndim == 1 and lon.ndim == 1:
            r = np.argmin(np.abs(lat - insitu_lat))
            c = np.argmin(np.abs(lon - insitu_lon))
        else:
            r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
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
        #filename = filepath.split('/')[-1].replace('.', '_') + '_extract.nc'
        filename = filepath.split('/')[-1][:-3]
        filename = f'{filename}_{r}_{c}_extract.nc'
        pdu = filepath.split('/')[-1]

        ofname = os.path.join(path_output, filename)
        # global_at['station_name'] =
        global_at['in_situ_lat'] = insitu_lat
        global_at['in_situ_lon'] = insitu_lon
        res = create_extract(ofname, pdu, options, nc_sat, global_at, lat, lon, r, c, None, None)
        if res:
            created = True
            print(f'[INFO] Extract file created: {ofname}')
    else:
        if args.verbose:
            print(f'[WARNING] Site out of the image')

    return created


def launch_create_extract(filepath, options):
    ncreated = 0

    path_output = get_output_path(options)
    if path_output is None:
        print(f'ERROR: {path_output} is not valid')
        return ncreated

    # Retrieving sites
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
    var_lat, var_lon = get_lat_long_var_names(options)
    lat, lon = get_lat_long_arrays(nc_sat, var_lat, var_lon)

    # Retrieving global atribbutes
    if args.verbose:
        print('[INFO] Retrieving global attributes...')
    global_at = get_global_atrib(nc_sat, options)

    # Working for each site, checking if there is in the image
    for site in sites:
        if args.verbose:
            print(f'[INFO]Working for site: {site}')
        insitu_lat = sites[site]['latitude']
        insitu_lon = sites[site]['longitude']
        contain_flag = 0
        if cfs.contain_location(lat, lon, insitu_lat, insitu_lon) == 1:
            if lat.ndim == 1 and lon.ndim == 1:
                r = np.argmin(np.abs(lat - insitu_lat))
                c = np.argmin(np.abs(lon - insitu_lon))
            else:
                r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
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


# def check_contain_flag():


def create_extract(ofname, pdu, options, nc_sat, global_at, lat, long, r, c, skie_file, irows):
    size_box = get_box_size(options)
    start_idx_x = (c - int(size_box / 2))
    stop_idx_x = (c + int(size_box / 2) + 1)
    start_idx_y = (r - int(size_box / 2))
    stop_idx_y = (r + int(size_box / 2) + 1)

    window = [start_idx_y, stop_idx_y, start_idx_x, stop_idx_x]
    print('=====================')
    print(size_box)
    print(window)
    print('lat variables',lat.shape)
    print('lon varaible:',long.shape)

    search_pattern = 'rrs_'
    wl_atrib = None
    if options.has_option('satellite_options', 'rrs_prefix'):
        search_pattern = options['satellite_options']['rrs_prefix']
        if not search_pattern.endswith('_'):
            search_pattern = f'{search_pattern}_'
    if options.has_option('satellite_options', 'wl_atrib'):
        wl_atrib = options['satellite_options']['wl_atrib']
    reflectance_bands, n_bands = get_reflectance_bands_info(nc_sat, search_pattern, wl_atrib)

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
        print(f'[INFO]    Creating file: {ofname}')

    newEXTRACT.set_global_attributes(global_at)

    if skie_file is not None:
        newEXTRACT.create_dimensions_incluidinginsitu(size_box, n_bands, skie_file.get_n_bands(), 30)
    else:
        newEXTRACT.create_dimensions(size_box, n_bands)

    newEXTRACT.create_lat_long_variables(lat, long, window)

    # Sat time start:  ,+9-2021-12-24T18:23:00.471Z
    if 'start_date' in nc_sat.ncattrs():
        sat_time = dt.strptime(nc_sat.start_date, '%Y-%m-%d')
        sat_time = sat_time.replace(hour=0, minute=0, second=0, microsecond=0)
        newEXTRACT.create_satellite_time_variable(sat_time)
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

    # Rrs and wavelenghts
    satellite_Rrs = newEXTRACT.create_rrs_variable(global_at['sensor'])
    rbands = list(reflectance_bands.keys())
    wavelenghts = []
    for index in range(len(rbands)):
        rband = rbands[index]
        rbandvar = nc_sat.variables[rband]
        if rbandvar.ndim == 3:
            bandarray = ma.array(rbandvar[:, :, :])
            satellite_Rrs[0, index, :, :] = bandarray[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
        elif rbandvar.ndim == 2:
            bandarray = ma.array(rbandvar[:, :])
            satellite_Rrs[0, index, :, :] = bandarray[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
        wl = reflectance_bands[rband]['wavelenght']
        wavelenghts.append(wl)
    newEXTRACT.create_satellite_bands_variable(wavelenghts)

    # flags
    # flag_band = nc_sat.variables[flag_band_name]
    #
    # newEXTRACT.create_flag_variable(f'satellite_{flag_band_name}', flag_band, flag_band.long_name, flag_band.flag_masks,
    #                                 flag_band.flag_meanings, window)

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
            if args.atm_correction == 'FUB':
                sat_time_here = get_sat_time_from_fname(pdu)
            else:
                sat_time_here = dt.strptime(nc_sat.start_date, '%Y-%m-%d')
            time_diff = float(abs((insitu_time_here - sat_time_here).total_seconds()))
            insitu_timedif_var[0, index_var] = time_diff
            latp, lonp = skie_file.get_lat_lon_at_subdb(irow)
            insitu_lat[0, index_var] = latp
            insitu_lon[0, index_var] = lonp
            index_var = index_var + 1

    newEXTRACT.close_file()

    return True


def create_extract_multiple(ofname, pdu, options, nc_files, global_at, lat, long, r, c):
    size_box = get_box_size(options)
    start_idx_x = (c - int(size_box / 2))
    stop_idx_x = (c + int(size_box / 2) + 1)
    start_idx_y = (r - int(size_box / 2))
    stop_idx_y = (r + int(size_box / 2) + 1)

    window = [start_idx_y, stop_idx_y, start_idx_x, stop_idx_x]

    n_bands = len(nc_files)

    if n_bands == 0:
        print('[ERROR] reflectance bands are not defined')
        return False

    newEXTRACT = SatExtract(ofname)
    if not newEXTRACT.FILE_CREATED:
        print(f'[ERROR] File {ofname} could not be created')
        return False

    if args.verbose:
        print(f'[INFO]    Creating file: {ofname}')

    newEXTRACT.set_global_attributes(global_at)

    newEXTRACT.create_dimensions(size_box, n_bands)

    newEXTRACT.create_lat_long_variables(lat, long, window)

    # Sat time start:  ,+9-2021-12-24T18:23:00.471Z
    nc_sat = Dataset(nc_files[0])
    if 'start_date' in nc_sat.ncattrs():
        sat_time = dt.strptime(nc_sat.start_date, '%Y-%m-%d')
        sat_time = sat_time.replace(hour=0, minute=0, second=0, microsecond=0)
        newEXTRACT.create_satellite_time_variable(sat_time)
    else:
        sat_time = get_sat_time_from_fname(pdu)
        if sat_time is not None:
            newEXTRACT.create_satellite_time_variable(sat_time)
        else:
            print(f'[ERROR] Satellite time is not defined...')
            newEXTRACT.close_file()
            return False
    nc_sat.close()

    # pdu variable
    newEXTRACT.create_pdu_variable(pdu, global_at['sensor'])

    # Rrs and wavelenghts
    satellite_Rrs = newEXTRACT.create_rrs_variable(global_at['sensor'])
    wavelenghts = []
    for index in range(len(nc_files)):
        f = nc_files[index]
        nc_sat = Dataset(f)
        name_file = f.split('/')[-1]
        rband = name_file.split('-')[1].upper()
        wls = rband[3:].replace('_', '.')
        wl = float(wls)
        wavelenghts.append(wl)
        rbandvar = nc_sat.variables[rband]
        if rbandvar.ndim == 3:
            bandarray = ma.array(rbandvar[:, :, :])
            satellite_Rrs[0, index, :, :] = bandarray[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
        elif rbandvar.ndim == 2:
            bandarray = ma.array(rbandvar[:, :])
            satellite_Rrs[0, index, :, :] = bandarray[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
        nc_sat.close()

    newEXTRACT.create_satellite_bands_variable(wavelenghts)

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


def check_location(insitu_lat, insitu_lon, lat, lon, size_box):
    contain_flag = False
    r = -1
    c = -1
    if cfs.contain_location(lat, lon, insitu_lat, insitu_lon) == 1:
        if lat.ndim == 1 and lon.ndim == 1:
            r = np.argmin(np.abs(lat - insitu_lat))
            c = np.argmin(np.abs(lon - insitu_lon))
        else:
            r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
        start_idx_x = (c - int(size_box / 2))
        stop_idx_x = (c + int(size_box / 2) + 1)
        start_idx_y = (r - int(size_box / 2))
        stop_idx_y = (r + int(size_box / 2) + 1)

        # print('Check location line 336: ', r, c, start_idx_y,stop_idx_y,start_idx_x,stop_idx_x)

        if lat.ndim == 1 and lon.ndim == 1:
            if start_idx_y >= 0 and (stop_idx_y + 1) < lat.shape[0] and start_idx_x >= 0 and (stop_idx_x + 1) < \
                    lon.shape[0]:
                contain_flag = True
        else:
            if start_idx_y >= 0 and (stop_idx_y + 1) < lat.shape[0] and start_idx_x >= 0 and (stop_idx_x + 1) < \
                    lat.shape[1]:
                contain_flag = True

    return contain_flag, r, c


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


def get_sat_time(nc_sat, pdu):
    sat_time = None
    if 'start_date' in nc_sat.ncattrs():
        sat_time = dt.strptime(nc_sat.start_date, '%Y-%m-%d')
        sat_time = sat_time.replace(hour=11, minute=0, second=0)
    else:
        sat_time = get_sat_time_from_fname(pdu)
    return sat_time


def get_site_options(options):
    site_file = None
    site_list = None
    region_list = None
    if not options.has_option('Time_and_sites_selection', 'sites'):
        print(f'[ERROR] section Time_and_sites_selection, option sites not found in configuration file')
    if options.has_option('Time_and_sites_selection', 'sites_file') and options['Time_and_sites_selection'][
        'sites_file']:
        site_file = options['Time_and_sites_selection']['sites_file']
    if options['Time_and_sites_selection']['sites']:
        site_list = options['Time_and_sites_selection']['sites'].split(',')
        site_list = [s.strip() for s in site_list]
    if not options['Time_and_sites_selection']['sites'] and options.has_option('Time_and_sites_selection',
                                                                               'sites_region') and \
            options['Time_and_sites_selection']['sites_region']:
        region_list = options['Time_and_sites_selection']['sites_region'].split(',')
        region_list = [r.strip() for r in region_list]
    return site_file, site_list, region_list


def get_lat_long_arrays(nc_sat, var_lat, var_lon):
    vlat = nc_sat.variables[var_lat]
    vlon = nc_sat.variables[var_lon]
    if vlat.ndim == 2 and vlon.ndim == 2:
        lat = nc_sat.variables[var_lat][:, :]
        lon = nc_sat.variables[var_lon][:, :]
    if vlat.ndim == 1 and vlon.ndim == 1:
        lat = nc_sat.variables[var_lat][:]
        lon = nc_sat.variables[var_lon][:]
    return lat, lon


def get_lat_long_var_names(options):
    var_lat = 'lat'
    var_lon = 'lon'
    if options.has_option('satellite_options', 'lat_variable'):
        var_lat = options['satellite_options']['lat_variable']
    if options.has_option('satellite_options', 'lon_variable'):
        var_lon = options['satellite_options']['lon_variable']
    return var_lat, var_lon


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


def get_global_atrib(nc_sat, options):
    at = {}
    equiv = {
        'project': 'aco_processor',
        'product_version': 'proc_version',
        'grid_resolution': 'res'
    }
    if options.has_option('satellite_options', 'global_atrib'):
        list_str = options['satellite_options']['global_atrib']
        list_attributes = list_str.split(',')
        for atrib in list_attributes:
            atrib_name = atrib.strip()
            if atrib_name in equiv.keys():
                atrib_name = equiv[atrib_name]
            if atrib_name in nc_sat.ncattrs():
                at[atrib_name] = nc_sat.getncattr(atrib_name)
            else:
                print(f'[WARNING] Attribute {atrib_name} is not available in dataset. Skyping...')

    compulsory_keys = ['satellite', 'platform', 'sensor', 'res', 'aco_processor', 'proc_version']
    for key in compulsory_keys:
        if key not in at.keys():
            at[key] = ''

    at['station_name'] = ''
    at['in_situ_lat'] = -999
    at['in_situ_lon'] = -999

    return at


def get_reflectance_bands_info(nc_sat, search_pattern, wl_atrib):
    # search_pattern = 'rrs_'
    # if args.atm_correction == 'FUB':
    #     search_pattern = 'reflec_'

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

    if imax == 0 and imin == 100000:
        search_pattern = search_pattern[:-1]
        lw = len(search_pattern)
        for var in nc_sat.variables:
            if var.startswith(search_pattern):
                ival = var[lw:].strip()
                ival = ival.replace('_', '.')
                if wl_atrib is not None:
                    wl_band = nc_sat.variables[var].getncattr(wl_atrib)
                else:
                    wl_band = float(ival)
                reflectance_bands[var] = {'wavelenght': wl_band}
                nbands = nbands + 1
    else:
        for ival in range(imin, imax + 1):
            var_name = search_pattern + str(ival)
            if var_name in nc_sat.variables:
                if wl_atrib is not None:
                    wl_band = nc_sat.variables[var_name].getncattr(wl_atrib)
                else:
                    wl_band = ival
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
    # print('temp')
    return path_source, org, wce, time_start, time_stop


def check_product(filepath, time_start, time_stop):
    try:
        nc_sat = Dataset(filepath, 'r')
    except OSError:
        if args.verbose:
            print(f'[WARNING] {filepath} is not a valid dataset. Skipping...')
        return False

    checkTime = False
    if args.atm_correction == 'FUB':
        time_sat = get_sat_time_from_fname(filepath)
        if time_sat is not None:
            checkTime = True
    else:
        if 'start_date' in nc_sat.ncattrs():
            time_sat = dt.strptime(nc_sat.start_date, '%Y-%m-%d')
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
    flag_variables = ifobj.flag_variables
    flag_info = ifobj.flag_info

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
        print(f'[INFO] Maximum time difference between sat and in situ times: {nminutes} minutes')

    maxtimediff = nminutes * 60
    size_box = get_box_size(options)
    ncreated = 0
    formatdt = '%Y-%m-%d %H:%M:%S'
    for filepath in product_path_list:
        if args.verbose:
            print('------------------------------')
            pname = filepath.split('/')[-1]
            print(f'[INFO]Checking sat extracts for product: {pname}')
        nc_sat = Dataset(filepath, 'r')
        var_lat, var_lon = get_lat_long_var_names(options)
        lat, lon = get_lat_long_arrays(nc_sat, var_lat, var_lon)

        # sat_time = dt.strptime(nc_sat.start_time, formatdt)
        pdu = filepath.split('/')[-1]
        sat_time = get_sat_time(nc_sat, pdu)
        sat_date_str = sat_time.strftime('%Y-%m-%d')

        if sat_date_str in valid_dates.keys():
            for h in valid_dates[sat_date_str]:
                insitu_time_str = f'{sat_date_str} {h}'
                insitu_time = dt.strptime(insitu_time_str, formatdt)
                insitu_lat = valid_dates[sat_date_str][h]['lat']
                insitu_lon = valid_dates[sat_date_str][h]['lon']
                timediff = abs((insitu_time - sat_time).total_seconds())
                contain_flag, r, c = check_location(insitu_lat, insitu_lon, lat, lon, size_box)
                if timediff < maxtimediff and contain_flag:
                    if args.verbose:
                        print(
                            f'[INFO] Preparing extract. In situ latitude: {insitu_lat} Longitude: {insitu_lon} Time: {insitu_time_str}')
                    itime = insitu_time.strftime('%H%M%S')
                    extract = f'{itime}_{r}_{c}'
                    filename = filepath.split('/')[-1].replace('.', '_') + '_extract_insitu_' + extract + '.nc'
                    ofname = os.path.join(path_output, filename)
                    global_at = get_global_atrib(nc_sat, options)
                    global_at['station_name'] = 'in situ dataset'
                    global_at['in_situ_lat'] = insitu_lat
                    global_at['in_situ_lon'] = insitu_lon
                    insitu_info = valid_dates[sat_date_str][h]

                    created = create_extract_insitu(ofname, pdu, options, nc_sat, global_at, r, c, insitu_time,
                                                    insitu_info, name_variables, flag_variables, flag_info)
                    if created:
                        ncreated = ncreated + 1
        nc_sat.close()
        # nhere = 1  # nhere = launch_create_extract_syke(filepath, syke_file, options)
        # ncreated = ncreated + nhere
    print('------------------------------')
    print(f'COMPLETED. {ncreated} sat extract files were created')


def create_extract_insitu(ofname, pdu, options, nc_sat, global_at, r, c, insitu_time, insitu_info, name_variables,
                          flag_variables, flag_info):
    size_box = get_box_size(options)
    start_idx_x = (c - int(size_box / 2))
    stop_idx_x = (c + int(size_box / 2) + 1)
    start_idx_y = (r - int(size_box / 2))
    stop_idx_y = (r + int(size_box / 2) + 1)
    window = [start_idx_y, stop_idx_y, start_idx_x, stop_idx_x]

    search_pattern = 'rrs_'
    wl_atrib = None
    if options.has_option('satellite_options', 'rrs_prefix'):
        search_pattern = options['satellite_options']['rrs_prefix']
        if not search_pattern.endswith('_'):
            search_pattern = f'{search_pattern}_'
    if options.has_option('satellite_options', 'wl_atrib'):
        wl_atrib = options['satellite_options']['wl_atrib']
    reflectance_bands, n_bands = get_reflectance_bands_info(nc_sat, search_pattern, wl_atrib)
    if n_bands == 0:
        print('[ERROR] reflectance bands are not defined')
        return False

    newEXTRACT = SatExtract(ofname)
    if not newEXTRACT.FILE_CREATED:
        print(f'[ERROR] File {ofname} could not be created')
        return False

    if args.verbose:
        print(f'[INFO]Creating file: {ofname}')
    newEXTRACT.set_global_attributes(global_at)

    newEXTRACT.create_dimensions(size_box, n_bands)

    var_lat, var_lon = get_lat_long_var_names(options)
    lat, long = get_lat_long_arrays(nc_sat, var_lat, var_lon)
    newEXTRACT.create_lat_long_variables(lat, long, window)

    # Sat time start:  ,+9-2021-12-24T18:23:00.471Z
    sat_time = get_sat_time(nc_sat, pdu)
    newEXTRACT.create_satellite_time_variable(sat_time)

    # pdu variable
    newEXTRACT.create_pdu_variable(pdu, global_at['sensor'])

    satellite_Rrs = newEXTRACT.create_rrs_variable(global_at['sensor'])

    rbands = list(reflectance_bands.keys())
    wavelenghts = []
    for index in range(len(rbands)):
        rband = rbands[index]
        bandarray = ma.array(nc_sat.variables[rband][:, :])
        satellite_Rrs[0, index, :, :] = bandarray[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
        wl = reflectance_bands[rband]['wavelenght']
        wavelenghts.append(wl)
    newEXTRACT.create_satellite_bands_variable(wavelenghts)

    ##add in situ variables
    timediff = abs((insitu_time - sat_time).total_seconds())
    insitulat_var, insitulon_var, insitutime_var, time_difference_var = newEXTRACT.create_insitu_variables_for_single_insitu_data()
    insitulat_var[0] = global_at['in_situ_lat']
    insitulon_var[0] = global_at['in_situ_lon']
    insitutime_var[0] = insitu_time.timestamp()
    time_difference_var[0] = timediff
    for name_var in name_variables:
        insitu_var = newEXTRACT.create_insitu_variable_for_single_insitu_data(name_var, 'unknown', 'unknown')
        insitu_var[0] = insitu_info[name_var]
    for name_var in flag_variables:
        insitu_flag_var = newEXTRACT.create_insitu_flag_variable(name_var, flag_info[name_var]['flag_masks'],
                                                                 flag_info[name_var]['flag_meanings_str'])
        insitu_flag_var[0] = insitu_info[name_var]
    newEXTRACT.close_file()

    return True


def run_cmems_option(options):
    if args.verbose:
        print('[INFO] Started CMEMS option...')

    ncreated = 0
    path_output = get_output_path(options)
    if path_output is None:
        print(f'ERROR: {path_output} is not valid')
        return ncreated

    # Retrieving sites
    site_file, site_list, region_list = get_site_options(options)
    if site_file is not None:
        sites = cfs.get_sites_from_file(site_file, site_list, region_list, path_output)
    else:
        sites = cfs.get_sites_from_list(site_list, path_output)

    if len(sites) == 0:
        print('[ERROR] No sites are defined')
        return ncreated

    # path_code_eistools = '/home/Luis.Gonzalezvilas/eistools'
    path_code_eistools = '/store/COP2-OC-TAC/CODE/eistools'
    sys.path.append(path_code_eistools)
    import product_info
    import reformatCMEMS_202207_class
    reformat = reformatCMEMS_202207_class.ReformatCMEMS()

    product_name = options['file_path']['cmems_product']
    dataset_name = options['file_path']['cmems_dataset']
    pinfo = product_info.ProductInfo()
    pinfo.MODE = 'REFORMAT'
    pinfo.set_dataset_info(product_name, dataset_name)
    flist = options['Time_and_sites_selection']['time_list']
    ff = open(flist, 'r')
    for line in ff:
        strdate = line.strip()
        try:
            date = dt.strptime(strdate, '%Y-%m-%d')
            if args.verbose:
                print('----------------------------------')
                print(f'[INFO] Date: {strdate}')
            ##checking if output files already exist
            filesExist = True
            expected_file_nc = pinfo.get_file_path_orig_name(None, date)
            for site in sites:
                path_output_site = os.path.join(path_output, site)
                filename = expected_file_nc.split('/')[-1].replace('.', '_') + '_extract_' + site + '.nc'
                ofname = os.path.join(path_output_site, filename)
                if not os.path.exists(ofname):
                    filesExist = False
            if filesExist:
                if args.verbose:
                    print(f'[INFO] Files for date: {strdate} already exist. Skipping...')
                continue

            # if not os.path.exists(expected_file_nc):
            #     reformat.make_reformat_daily_dataset(pinfo, date, date, args.verbose)
            # filenc = pinfo.get_file_path_orig(None, date)
            # if filenc is None:
            #     print(f'[WARNING] Refformatted file {expected_file_nc} could not be created. Skypping...')
            #     continue
            # if args.verbose:
            #     print(f'[INFO] Reformatted file {filenc}')
            # nhere = create_extract_cmems(filenc, options, sites, path_output)

            path_nc = os.path.dirname(expected_file_nc)
            if not os.path.exists(path_nc):
                print(f'[WARNING] Path nc files {path_nc} does not exist. Skypping...')
                continue
            nhere = create_extract_cmems_multiple(path_nc, date, options, sites, path_output)

            ncreated = ncreated + nhere
            # if args.verbose:
            #     print(f'[INFO] Removing file {filenc}')
            # os.remove(filenc)

        except:
            print('ERROR FILE')
            pass
    ff.close()
    if args.verbose:
        print('*********************************************************************************')
        print(f'[INFO] Extraction completed. Sat extracts created: {ncreated}')


def create_extract_cmems(filepath, options, sites, path_output):
    nc_sat = Dataset(filepath, 'r')
    # Retriving lat and long arrays
    if args.verbose:
        print('[INFO] Retrieving lat/long data...')
    var_lat, var_lon = get_lat_long_var_names(options)
    lat, lon = get_lat_long_arrays(nc_sat, var_lat, var_lon)

    # Retrieving global atribbutes
    if args.verbose:
        print('[INFO] Retrieving global attributes...')
    global_at = get_global_atrib(nc_sat, options)

    ncreated = 0
    # Working for each site, checking if there is in the image
    for site in sites:
        if args.verbose:
            print(f'[INFO]    Working for site: {site}')
        insitu_lat = sites[site]['latitude']
        insitu_lon = sites[site]['longitude']
        contain_flag = 0
        if cfs.contain_location(lat, lon, insitu_lat, insitu_lon) == 1:
            if lat.ndim == 1 and lon.ndim == 1:
                r = np.argmin(np.abs(lat - insitu_lat))
                c = np.argmin(np.abs(lon - insitu_lon))
            else:
                r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
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
            filename = filepath.split('/')[-1].replace('.', '_') + '_extract_' + site + '.nc'
            pdu = filepath.split('/')[-1]
            path_output_site = os.path.join(path_output, site)
            if not os.path.exists(path_output_site):
                os.mkdir(path_output_site)
            ofname = os.path.join(path_output_site, filename)
            global_at['station_name'] = site
            global_at['in_situ_lat'] = insitu_lat
            global_at['in_situ_lon'] = insitu_lon

            res = create_extract(ofname, pdu, options, nc_sat, global_at, lat, lon, r, c, None, None)
            if res:
                ncreated = ncreated + 1
                print(f'[INFO]    Extract file created: {ofname}')
            nc_sat.close()
        else:
            if args.verbose:
                print(f'[WARNING] Site {site} out of the image')

    return ncreated


def create_extract_cmems_multiple(ncpath, date, options, sites, path_output):
    if args.verbose:
        print(f'[INFO] NC Path: {ncpath}')
    strdate = date.strftime('%Y%j')
    band_list = ['400', '412_5', '442_5', '490', '510', '560', '620', '665', '673_75', '681_25', '708_75', '753_75',
                 '778_75']
    ncfiles = []
    for b in band_list:
        name = f'O{strdate}-rrs{b}-med-fr.nc'
        fname = os.path.join(ncpath, name)
        print(fname, os.path.exists(fname))
        ncfiles.append(fname)

    nc_sat = Dataset(ncfiles[0], 'r')

    # Retriving lat and long arrays
    if args.verbose:
        print('[INFO] Retrieving lat/long data...')
    var_lat, var_lon = get_lat_long_var_names(options)
    lat, lon = get_lat_long_arrays(nc_sat, var_lat, var_lon)

    # Retrieving global atribbutes
    if args.verbose:
        print('[INFO] Retrieving global attributes...')
    global_at = get_global_atrib(nc_sat, options)

    nc_sat.close()

    ncreated = 0
    # Working for each site, checking if there is in the image
    for site in sites:
        if args.verbose:
            print(f'[INFO]    Working for site: {site}')
        insitu_lat = sites[site]['latitude']
        insitu_lon = sites[site]['longitude']
        contain_flag = 0
        if cfs.contain_location(lat, lon, insitu_lat, insitu_lon) == 1:
            if lat.ndim == 1 and lon.ndim == 1:
                r = np.argmin(np.abs(lat - insitu_lat))
                c = np.argmin(np.abs(lon - insitu_lon))
            else:
                r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)
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
            # CMEMS2_O2021357 - rrs - med - fr_nc_extract_Venise.nc
            # filename = filepath.split('/')[-1].replace('.', '_') + '_extract_' + site + '.nc'
            # pdu = filepath.split('/')[-1]
            pdu = ncpath
            filename = f'CMEMS2_O{strdate}-rrs-med-fr_nc_extract_{site}.nc'

            path_output_site = os.path.join(path_output, site)
            if not os.path.exists(path_output_site):
                os.mkdir(path_output_site)
            ofname = os.path.join(path_output_site, filename)
            global_at['station_name'] = site
            global_at['in_situ_lat'] = insitu_lat
            global_at['in_situ_lon'] = insitu_lon
            res = create_extract_multiple(ofname, pdu, options, ncfiles, global_at, lat, lon, r, c)
            if res:
                ncreated = ncreated + 1
                print(f'[INFO]    Extract file created: {ofname}')

        else:
            if args.verbose:
                print(f'[WARNING] Site {site} out of the image')

    return ncreated


def get_cmems_product_day(path_source, org, datehere, dataset):
    path_day = path_source
    if org is not None:
        if org == 'YYYYjjj':
            yearstr = datehere.strftime('%Y')
            jjjstr = datehere.strftime('%j')
            path_day = os.path.join(path_source, yearstr, jjjstr)
    datefile = datehere.strftime('%Y%m%d')
    file = os.path.join(path_day, f'{datefile}_{dataset}.nc')
    if not os.path.exists(file):
        print(f'[WARNING] File: {file} does not exist. Skiping...')
        return None
    return file


def main():
    print('[INFO] Creating satellite extracts')

    if not args.config_file:
        return

    options = config_reader(args.config_file)

    if options.has_section('CSV_SELECTION') and options.has_option('CSV_SELECTION', 'path_csv') and options.has_option(
            'CSV_SELECTION', 'dataset'):
        path_csv = options['CSV_SELECTION']['path_csv']
        dataset = options['CSV_SELECTION']['dataset']
        if not os.path.exists(path_csv):
            print(f'[ERROR] Path csv {path_csv} was not found')
            return
        try:
            df = pd.read_csv(path_csv, ';')
        except:
            print(f'[ERROR] File {path_csv} is not a valid csv separated by ;')
            return
        col_date = 'date'
        col_lat = 'lat'
        col_lon = 'lon'
        format_date = '%Y-%m-%dT%H:%M'
        if options.has_option('CSV_SELECTION', 'col_date'):
            col_date = options['CSV_SELECTION']['col_date']
        if options.has_option('CSV_SELECTION', 'col_lat'):
            col_lat = options['CSV_SELECTION']['col_lat']
        if options.has_option('CSV_SELECTION', 'col_lon'):
            col_lon = options['CSV_SELECTION']['col_lon']
        if options.has_option('CSV_SELECTION', 'format_date'):
            format_date = options['CSV_SELECTION']['format_date']
        path_source, org, wce, time_start, time_stop = get_find_product_info(options)
        ncreated = 0
        for idx, row in df.iterrows():
            try:
                datestr = row[col_date].strip()
                datehere = dt.strptime(datestr, format_date)
                lathere = float(row[col_lat])
                lonhere = float(row[col_lon])
            except:
                print(f'[WARNING] Row {idx} is not valid. Date, latitute and/or longitude could not be parsed')
                continue
            fproduct = get_cmems_product_day(path_source, org, datehere, dataset)
            if fproduct is not None:
                created = launch_create_extract_station(fproduct, options, lathere, lonhere)
                if created:
                    ncreated = ncreated + 1
            # print(datehere, lathere, lonhere,fproduct)
        return

    # if options.has_option('file_path', 'path_skie') and options.has_option('file_path', 'path_skie_code'):
    #     path_skie = options['file_path']['path_skie']
    #     path_skie_code = options['file_path']['path_skie_code']
    #     if not os.path.exists(path_skie):
    #         print(f'[ERROR] Skie file: {path_skie} is not avialable')
    #         return
    #     if not os.path.exists(path_skie_code) or not os.path.isdir(path_skie_code):
    #         print(f'[ERROR] Skie code folder: {path_skie_code} is not avialable')
    #         return
    #     if args.verbose:
    #         print(f'[INFO] Path skie: {path_skie}')
    #         print(f'[INFO] Path skie: {path_skie_code}')
    #     sys.path.append(path_skie_code)
    #     try:
    #         from skie_csv import SKIE_CSV
    #     except ModuleNotFoundError:
    #         print(f'[ERROR] Skie module is not found')
    #         return
    #     skie_file = SKIE_CSV(path_skie)
    #     skie_file.start_list_dates()
    #     skie_file.extract_wl_colnames()
    #
    #     print('------------------------------')
    #     path_source, org, wce, time_start, time_stop = get_find_product_info(options)
    #     path_to_list = get_path_list_products(options)
    #     if path_to_list is not None:
    #         product_path_list = get_list_products(path_to_list, path_source, org, wce, time_start, time_stop)
    #         if len(product_path_list) == 0:
    #             if args.verbose:
    #                 print('-----------------')
    #             print(f'[WARNING] No valid datasets found on:  {path_source}')
    #             print('------------------------------')
    #             print('[INFO] COMPLETED. 0 sat extract files were created')
    #             return
    #         ncreated = 0
    #         for filepath in product_path_list:
    #             if args.verbose:
    #                 print('------------------------------')
    #                 pname = filepath.split('/')[-1]
    #                 print(f'[INFO]Extracting sat extract for product: {pname}')
    #             nhere = launch_create_extract_skie(filepath, skie_file, options)
    #             ncreated = ncreated + nhere
    #         print('------------------------------')
    #         print(f'COMPLETED. {ncreated} sat extract files were created')
    #     return

    if options.has_option('file_path', 'path_insitu') and options.has_option('file_path', 'path_insitu_code'):
        run_insitu_option(options)
        return

    if options.has_option('Time_and_sites_selection', 'time_list') and \
            options.has_option('file_path', 'cmems_product') and \
            options.has_option('file_path', 'cmems_dataset'):
        run_cmems_option(options)
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
