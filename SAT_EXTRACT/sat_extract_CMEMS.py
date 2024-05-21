import argparse
import configparser
import os
import sys

import pandas as pd
import pytz
from netCDF4 import Dataset
from datetime import datetime as dt
from datetime import timedelta
import numpy.ma as ma
import numpy as np
import subprocess

# import user defined functions from other .py
import sat_extract

code_home = os.path.dirname(os.path.dirname(sat_extract.__file__))
sys.path.append(code_home)
from SAT_EXTRACT.sat_extract import SatExtract
import COMMON.common_functions as cfs

parser = argparse.ArgumentParser(description="Create Match-up DataBase files (MDB) files.")
parser.add_argument("-d", "--download_sources", help="Download DOORS sources.", action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-p', "--product_file", help="Image file.")
parser.add_argument('-ac', "--atm_correction", help="Atmospheric correction algorithm (Default: C2RCC)",
                    choices=["C2RCC", "FUB"])

args = parser.parse_args()


def copy_extract(input_file, output_file):
    from netCDF4 import Dataset
    input_dataset = Dataset(input_file)
    ncout = Dataset(output_file, 'w', format='NETCDF4')

    # copy global attributes all at once via dictionary
    ncout.setncatts(input_dataset.__dict__)

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    for name, variable in input_dataset.variables.items():
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(input_dataset[name].__dict__)

        ncout[name][:] = input_dataset[name][:]

    # ncout.close()
    input_dataset.close()

    newEXTRACT = SatExtract(None)
    newEXTRACT.EXTRACT = ncout
    newEXTRACT.FILE_CREATED = True

    return newEXTRACT


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
                           extracts[extract]['irows'], None)
        if b:
            ncreated = ncreated + 1

    nc_sat.close()

    return ncreated


def launch_create_extract_station(filepath, options, insitu_lat, insitu_lon, insitu_time):
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
        # filename = filepath.split('/')[-1].replace('.', '_') + '_extract.nc'
        filename = filepath.split('/')[-1][:-3]
        filename = f'{filename}_{r}_{c}_extract.nc'
        pdu = filepath.split('/')[-1]

        ofname = os.path.join(path_output, filename)
        # global_at['station_name'] =
        global_at['in_situ_lat'] = insitu_lat
        global_at['in_situ_lon'] = insitu_lon
        insitu_time = insitu_time.timestamp()
        irows = [insitu_lat, insitu_lon, insitu_time]
        res = create_extract(ofname, pdu, options, nc_sat, global_at, lat, lon, r, c, 'STATION', irows)
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

    search_pattern = 'rrs_'
    wl_atrib = None
    if options.has_option('satellite_options', 'rrs_prefix'):
        search_pattern = options['satellite_options']['rrs_prefix']
        if not search_pattern.endswith('_'):
            search_pattern = f'{search_pattern}_'
    if options.has_option('satellite_options', 'wl_atrib'):
        wl_atrib = options['satellite_options']['wl_atrib']
    reflectance_bands, n_bands = get_reflectance_bands_info(nc_sat, search_pattern, wl_atrib)

    # print(reflectance_bands)
    # print(n_bands)

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

    if skie_file is not None and skie_file != 'STATION':
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
    if pdu is not None:
        newEXTRACT.create_pdu_variable(pdu, global_at['sensor'])

    # Rrs and wavelenghts
    satellite_Rrs = newEXTRACT.create_rrs_variable(global_at['sensor'])
    rbands = list(reflectance_bands.keys())
    wavelenghts = []
    for index in range(len(rbands)):
        rband = rbands[index]
        rbandvar = nc_sat.variables[rband]
        if rbandvar.ndim == 3:
            # bandarray = ma.array(rbandvar[:, :, :])
            # satellite_Rrs[0, index, :, :] = bandarray[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
            bandarray = ma.array(rbandvar[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])
            satellite_Rrs[0, index, :, :] = bandarray[:, :]
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
        if skie_file == 'STATION':
            insitu_lat_here = irows[0]
            insitu_lon_here = irows[1]
            insitu_time_here = irows[2]
            insitu_time = newEXTRACT.EXTRACT.createVariable('insitu_time', 'f8', ('satellite_id',), zlib=True,
                                                            complevel=6)
            insitu_time.units = "Seconds since 1970-1-1"
            insitu_time.description = 'In situ time in ISO 8601 format (UTC).'
            insitu_time[0] = insitu_time_here
            insitu_lat = newEXTRACT.EXTRACT.createVariable('insitu_latitude', 'f8', ('satellite_id',),
                                                           fill_value=-999,
                                                           zlib=True, complevel=6)
            insitu_lat.short_name = "latitude"
            insitu_lat.units = "degrees"
            insitu_lon = newEXTRACT.EXTRACT.createVariable('insitu_longitude', 'f8', ('satellite_id',),
                                                           fill_value=-999,
                                                           zlib=True, complevel=6)
            insitu_lon.short_name = "longitude"
            insitu_lon.units = "degrees"
            insitu_lat[0] = insitu_lat_here
            insitu_lon[0] = insitu_lon_here

        else:
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


def add_variable_single(newEXTRACT, extract, variable_list, variable_list_out, rrs_var_list):
    file = extract['file']
    limits = extract['limits']
    start_idx_y = limits[0]
    stop_idx_y = limits[1]
    start_idx_x = limits[2]
    stop_idx_x = limits[3]
    input_dataset = None
    if variable_list[0].strip() == '*':
        variable_list = []
        input_dataset = Dataset(file, 'r')
        for var in input_dataset.variables:
            if var in rrs_var_list:
                continue
            if input_dataset.variables[var].ndim == 3:
                variable_list.append(var)
        variable_list_out = variable_list

    if input_dataset is None:
        input_dataset = Dataset(file, 'r')

    for idx in range(len(variable_list)):
        variable_in = variable_list[idx]
        variable_out = f'satellite_{variable_list_out[idx]}'
        var_in = input_dataset.variables[variable_in]
        var_array = ma.array(var_in[:])
        var_array = np.array(var_array.filled(-999.0))

        if variable_out not in newEXTRACT.EXTRACT.variables:
            variable = newEXTRACT.create_2D_variable_general(variable_out, var_array, limits)
        else:
            variable = newEXTRACT.EXTRACT.variables[variable_out]
            variable[0, :, :] = var_array[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]

        for at in var_in.ncattrs():
            if at == '_FillValue':
                continue
            variable.setncattr(at, var_in.getncattr(at))

    input_dataset.close()

    return newEXTRACT


def add_variable_multiple(newEXTRACT, extract, variable_list, variable_list_out):
    list_files = extract['list_files']
    limits = extract['limits']
    start_idx_y = limits[0]
    stop_idx_y = limits[1]
    start_idx_x = limits[2]
    stop_idx_x = limits[3]
    nvar = len(list_files)
    for idx in range(nvar):
        file_in = list_files[idx]
        variable_in = variable_list[idx]
        variable_out = f'satellite_{variable_list_out[idx]}'
        nc_in = Dataset(file_in, 'r')
        var_in = nc_in.variables[variable_in]
        var_array = ma.array(var_in[:])
        var_array = np.array(var_array.filled(-999.0))

        if variable_out not in newEXTRACT.EXTRACT.variables:
            variable = newEXTRACT.create_2D_variable_general(variable_out, var_array, limits)
        else:
            variable = newEXTRACT.EXTRACT.variables[variable_out]
            variable[0, :, :] = var_array[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]

        for at in var_in.ncattrs():
            if at == '_FillValue':
                continue
            variable.setncattr(at, var_in.getncattr(at))
        nc_in.close()

    return newEXTRACT


def add_reflectance_single(newEXTRACT, extract, wl_list, var_list):
    if not 'satellite_bands' in newEXTRACT.EXTRACT.dimensions:
        print(f'[ERROR] Dimension satellite bands is not defined')
        return
    global_at = extract['1']['global_at']
    file = extract['file']
    limits = extract['limits']
    start_idx_y = limits[0]
    stop_idx_y = limits[1]
    start_idx_x = limits[2]
    stop_idx_x = limits[3]

    nwl = len(wl_list)

    if len(var_list) != nwl:
        return

    if 'satellite_Rrs' in newEXTRACT.EXTRACT.variables:
        satellite_Rrs = newEXTRACT.EXTRACT.variables['satellite_Rrs']
    else:
        satellite_Rrs = newEXTRACT.create_rrs_variable(global_at['sensor'])
    input_dataset = Dataset(file)
    for idx in range(nwl):
        var_name = var_list[idx]
        variable = input_dataset.variables[var_name]
        if variable.ndim == 3:
            bandarray = ma.array(variable[:, :, :])
            satellite_Rrs[0, idx, :, :] = bandarray[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
        elif variable.ndim == 2:
            bandarray = ma.array(variable[:, :])
            satellite_Rrs[0, idx, :, :] = bandarray[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]

    input_dataset.close()

    if 'satellite_bands' in newEXTRACT.EXTRACT.variables:
        satellite_bands = newEXTRACT.EXTRACT.variables['satellite_bands']
        satellite_bands[:] = wl_list
    else:
        newEXTRACT.create_satellite_bands_variable(wl_list)

    return newEXTRACT


def add_reflectance_multiple(newEXTRACT, extract, wl_list):
    if not 'satellite_bands' in newEXTRACT.EXTRACT.dimensions:
        print(f'[ERROR] Dimension satellite bands is not defined')
        return
    global_at = extract['1']['global_at']
    list_files = extract['list_files']
    limits = extract['limits']
    start_idx_y = limits[0]
    stop_idx_y = limits[1]
    start_idx_x = limits[2]
    stop_idx_x = limits[3]

    nwl = len(list_files)
    if 'satellite_Rrs' in newEXTRACT.EXTRACT.variables:
        satellite_Rrs = newEXTRACT.EXTRACT.variables['satellite_Rrs']
    else:
        satellite_Rrs = newEXTRACT.create_rrs_variable(global_at['sensor'])
    wavelengths = []
    for iwl in range(nwl):
        wl = wl_list[iwl]
        wavelengths.append(float(wl))
        input_dataset = Dataset(list_files[iwl])
        for name, variable in input_dataset.variables.items():
            wls = str(wl).replace('.', '_')
            ifind = name.find(wls)
            if ifind >= 0:
                if variable.ndim == 3:
                    bandarray = ma.array(variable[:, :, :])
                    satellite_Rrs[0, iwl, :, :] = bandarray[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
                elif variable.ndim == 2:
                    bandarray = ma.array(variable[:, :])
                    satellite_Rrs[0, iwl, :, :] = bandarray[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
        input_dataset.close()

    if 'satellite_bands' in newEXTRACT.EXTRACT.variables:
        satellite_bands = newEXTRACT.EXTRACT.variables['satellite_bands']
        satellite_bands[:] = wavelengths
    else:
        newEXTRACT.create_satellite_bands_variable(wavelengths)

    return newEXTRACT


def add_insitu_basic_info(newEXTRACT, extract, id, nid, csv_flag_menings):
    if 'insitu_id' not in newEXTRACT.EXTRACT.dimensions:
        newEXTRACT.create_dimension_insitu(nid)
    idx = id - 1
    ids = str(id)
    satellite_time = extract['satellite_time']
    insitu_time = extract[ids]['insitu_time']

    global_at = extract[ids]['global_at']
    insitu_lat = global_at['in_situ_lat']
    insitu_lon = global_at['in_situ_lon']

    if 'insitu_time' not in newEXTRACT.EXTRACT.variables:
        newEXTRACT.create_insitu_time_variable()
    if 'insitu_latitude' not in newEXTRACT.EXTRACT.variables and 'insitu_longitude' not in newEXTRACT.EXTRACT.variables:
        newEXTRACT.create_insitu_lat_long_variables()
    if 'time_difference' not in newEXTRACT.EXTRACT.variables:
        newEXTRACT.create_insitu_time_difference_variable()

    newEXTRACT.EXTRACT.variables['insitu_time'][0, idx] = float(insitu_time.timestamp())
    newEXTRACT.EXTRACT.variables['insitu_latitude'][0, idx] = insitu_lat
    newEXTRACT.EXTRACT.variables['insitu_longitude'][0, idx] = insitu_lon
    newEXTRACT.EXTRACT.variables['time_difference'][0, idx] = abs(insitu_time - satellite_time).total_seconds()

    if csv_flag_menings is None:
        return newEXTRACT

    for flag in csv_flag_menings:
        flag_meanings = csv_flag_menings[flag]
        if len(flag_meanings) > 32:
            continue
        flag_v = f'insitu_flag_{flag.lower()}'
        if flag_v not in newEXTRACT.EXTRACT.variables:
            newEXTRACT.create_insitu_flag_variable_version2(flag_v, flag_meanings)
        flag = global_at[flag]
        value_flag = 0
        if flag in flag_meanings:
            index_flag = flag_meanings.index(flag)
            value_flag = int(np.power(2, float(index_flag)))
        newEXTRACT.EXTRACT.variables[flag_v][0, idx] = value_flag

    return newEXTRACT


def start_extract(extract, ofname):
    global_at = extract['1']['global_at']
    satellite_time = extract['satellite_time']
    size_box = extract['size_box']
    n_bands = extract['n_bands']
    window = extract['limits']
    if 'list_files' in extract:
        nc_file = extract['list_files'][0]
    elif 'file' in extract:
        nc_file = extract['file']

    # print(global_at['in_situ_lat'], global_at['in_situ_lon'])

    newEXTRACT = SatExtract(ofname)
    if not newEXTRACT.FILE_CREATED:
        print(f'[ERROR] File {ofname} could not be created')
        return False

    if args.verbose:
        print(f'[INFO] Starting file: {ofname}')

    newEXTRACT.set_global_attributes(global_at)

    newEXTRACT.create_dimensions(size_box, n_bands)

    var_lat = extract['var_lat']
    var_lon = extract['var_lon']
    nc_sat = Dataset(nc_file, 'r')
    lat, long = get_lat_long_arrays(nc_sat, var_lat, var_lon)
    nc_sat.close()
    newEXTRACT.create_lat_long_variables(lat, long, window)
    newEXTRACT.create_satellite_time_variable(satellite_time)

    return newEXTRACT


def create_extract_multiple(ofname, pdu, options, nc_files, band_list, global_at, lat, long, r, c):
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
    newEXTRACT.EXTRACT.in_situ_lat = global_at['in_situ_lat']
    newEXTRACT.EXTRACT.in_situ_lon = global_at['in_situ_lon']

    newEXTRACT.create_dimensions(size_box, n_bands)

    newEXTRACT.create_lat_long_variables(lat, long, window)

    # Sat time start:  ,+9-2021-12-24T18:23:00.471Z
    nc_sat = Dataset(nc_files[0])
    if 'start_date' in nc_sat.ncattrs():
        sat_time = dt.strptime(nc_sat.start_date, '%Y-%m-%d').astimezone(pytz.utc)
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
    if pdu is not None:
        newEXTRACT.create_pdu_variable(pdu, global_at['sensor'])

    # Rrs and wavelenghts
    satellite_Rrs = newEXTRACT.create_rrs_variable(global_at['sensor'])
    wavelenghts = []
    for index in range(len(nc_files)):
        f = nc_files[index]
        nc_sat = Dataset(f)
        name_file = f.split('/')[-1]
        rband = name_file.split('-')[1].upper()
        wls = band_list[index].replace('_', '.')
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


def get_band_list(options):
    if options.has_option('satellite_options', 'band_list'):
        list_str = options['satellite_options']['band_list']
        band_list = [x.strip() for x in list_str.split(',')]
        return band_list
    else:
        band_list = ['400', '412_5', '442_5', '490', '510', '560', '620', '665', '673_75', '681_25', '708_75', '753_75',
                     '778_75']
        return band_list


def get_format_name(options):
    if options.has_option('satellite_options', 'format_name'):
        format_name = options['satellite_options']['format_name'].strip()
    else:
        format_name = 'O$DATE$-rrs$BAND$-med-fr.nc'
    if options.has_option('satellite_options', 'format_name_date'):
        format_name_date = options['satellite_options']['format_name_date'].strip()
    else:
        format_name_date = '%Y%j'

    return format_name, format_name_date


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


def add_insitu_global_atrib(at, site, latitude, longitude, other):
    at['site'] = site
    at['in_situ_lat'] = latitude
    at['in_situ_lon'] = longitude
    if other is not None:
        for key in other:
            at[key] = other[key]
    return at


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
                print(f'[WARNING] Attribute {atrib_name} is not available in dataset. Skiping...')

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
    time_start = None
    time_stop = None
    section = 'Time_and_sites_selection'
    if options.has_section(section):
        if options.has_option(section, 'time_start'):
            time_start = dt.strptime(options['Time_and_sites_selection']['time_start'], '%Y-%m-%d')
        if options.has_option(section, 'time_stop'):
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


def run_cmems_option_noreformat(options):
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

    if len(sites) > 1:
        print(f'[WARNING] Script implemented for only one site')
        return ncreated

    site = list(sites.keys())[0]
    print(site)
    product_name = options['file_path']['cmems_product']
    dataset_name = options['file_path']['cmems_dataset']
    path_source = options['file_path']['sat_source_dir']

    flist = options['Time_and_sites_selection']['time_list']
    ff = open(flist, 'r')
    for line in ff:
        strdate = line.strip()
        try:
            date = dt.strptime(strdate, '%Y-%m-%d')
            datestr = date.strftime('%Y%m%d')
            yyyy = date.strftime('%Y')
            jjj = date.strftime('%j')
            if args.verbose:
                print('----------------------------------')
                print(f'[INFO] Date: {strdate}')

            filename = f'{dataset_name}_extract_{datestr}_{site}.nc'
            ofname = os.path.join(path_output, filename)

            if os.path.exists(ofname):
                if args.verbose:
                    print(f'[INFO] Extract file for date: {strdate} already exist. Skipping...')
                continue

            path_nc = os.path.join(path_source, yyyy, jjj)
            if not os.path.exists(path_nc):
                print(f'[WARNING] Path nc files {path_nc} does not exist. Skipping...')
                continue
            nhere = create_extract_cmems_multiple(path_nc, date, options, sites, ofname)

            ncreated = ncreated + nhere
            # if args.verbose:
            #     print(f'[INFO] Removing file {filenc}')
            # os.remove(filenc)

        except:
            print('[ERROR] Error creating extract...')
            pass
    ff.close()
    if args.verbose:
        print('*********************************************************************************')
        print(f'[INFO] Extraction completed. Sat extracts created: {ncreated}')


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

    path_code_eistools = '/home/Luis.Gonzalezvilas/eistools'
    # path_code_eistools = '/store/COP2-OC-TAC/CODE/eistools'
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
                print(f'[WARNING] Path nc files {path_nc} does not exist. Skipping...')
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
            # pdu = filepath.split('/')[-1]
            pdu = None
            path_output_site = os.path.join(path_output, site)
            if not os.path.exists(path_output_site):
                os.mkdir(path_output_site)
            ofname = os.path.join(path_output_site, filename)
            global_at['site'] = site
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


def get_geo_info(options, file_nc, insitu_lat, insitu_lon):
    nc_sat = Dataset(file_nc, 'r')
    var_lat, var_lon = get_lat_long_var_names(options)
    lat, lon = get_lat_long_arrays(nc_sat, var_lat, var_lon)
    contain_flag = 0
    limits = None
    rc = None
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
            limits = [start_idx_y, stop_idx_y, start_idx_x, stop_idx_x]
            rc = [r, c]

    nc_sat.close()
    return limits, rc


def create_extract_cmems_multiple(ncpath, date, options, sites, ofname):
    if args.verbose:
        print(f'[INFO] NC Path: {ncpath}')

    band_list = get_band_list(options)
    format_name, format_name_date = get_format_name(options)
    strdate = date.strftime(format_name_date)
    if args.verbose:
        print(f'[INFO] Band list: {band_list}')
        print(f'[INFO] Format name: {format_name}')
        print(f'[INFO] Format name date: {format_name_date}')
    ncfiles = []
    for b in band_list:
        # name = f'O{strdate}-rrs{b}-med-fr.nc'
        bstr = f'{float(b):.2f}'
        if bstr.endswith('.00'):
            bstr = bstr[:-3]
        if bstr.find('.') > 0 and bstr.endswith('0'):
            bstr = bstr[:-1]
        bstr = bstr.replace('.', '_')
        name = format_name
        name = name.replace('$DATE$', strdate)
        name = name.replace('$BAND$', bstr)
        fname = os.path.join(ncpath, name)
        if args.verbose:
            print(f'[INFO] {fname} -> {os.path.exists(fname)}')
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
    # global_at = get_global_atrib(nc_sat, options)
    global_at = get_satellite_global_atrib_from_options(options)

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
            pdu = None  # not more implemented
            global_at['site'] = site
            global_at['in_situ_lat'] = insitu_lat
            global_at['in_situ_lon'] = insitu_lon
            res = create_extract_multiple(ofname, pdu, options, ncfiles, band_list, global_at, lat, lon, r, c)
            if res:
                ncreated = ncreated + 1
                print(f'[INFO]    Extract file created: {ofname}')

        else:
            if args.verbose:
                print(f'[WARNING] Site {site} out of the image')

    return ncreated


def get_cmems_product_day_dataset(path_source, org, datehere, dataset):
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


def get_cmems_product_day(path_source, org, datehere, dataset_name_file, dataset_name_format_date,cmems_download_options):
    path_day = path_source
    if org is not None:
        if org == 'YYYYjjj':
            yearstr = datehere.strftime('%Y')
            jjjstr = datehere.strftime('%j')
            path_year = os.path.join(path_source, yearstr)
            path_day = os.path.join(path_source, yearstr, jjjstr)
    datefile = datehere.strftime(dataset_name_format_date)
    namefile = dataset_name_file.replace('$DATE$', datefile)
    file = os.path.join(path_day, f'{namefile}')

    if not os.path.exists(file):
        ods = None
        if org == 'YYYYjjj':
            ods = '%Y/%j'
            path_year = os.path.join(path_source, yearstr)
            if not os.path.isdir(path_year):
                os.mkdir(path_year)
            if not os.path.isdir(path_day):
                os.mkdir(path_day)

        ##DONWLOAD CERTO SOURCES
        if namefile.find('CERTO') > 0:
            folder_cmd = 'OLCI'
            if namefile.find('OLCI') < 0:
                folder_cmd = 'MSI'
            cmd = f'wget --user=rsg_dump --password=yohlooHohw2Pa9ohv1Chi ftp://ftp.rsg.pml.ac.uk/DOORS_matchups/{folder_cmd}/{namefile} -O {file}'
            if args.verbose:
                print(f'[INFO] Trying download with cmd:')
                print(f'[INFO] {cmd}')

            proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            try:
                outs, errs = proc.communicate(timeout=1800)
            except subprocess.TimeoutExpired:
                proc.kill()
                outs, errs = proc.communicate()

        ##DOWNLOAD CMEMS SOURCES
        if cmems_download_options is not None:
            code_eistools = os.path.join(os.path.dirname(code_home),'eistools')
            if os.path.exists(code_eistools):
                sys.path.append(code_eistools)
                from cmems_lois import CMEMS_LOIS
                clois = CMEMS_LOIS(args.verbose)
                cmems_download_options['start_date'] = datehere
                cmems_download_options['end_date'] = datehere
                print(cmems_download_options)
                print(path_source)
                print(ods)
                clois.make_cmems_download(cmems_download_options, True, path_source, ods,True)

            else:
                print(f'[WARNING] Code {code_eistools} for downloading is not available')

    if os.path.exists(file) and os.stat(file).st_size == 0:
        os.remove(file)
        if org == 'YYYYjjj':
            if len(os.listdir(path_day)) == 0:
                os.rmdir(path_day)
            if len(os.listdir(path_year)) == 0:
                os.rmdir(path_year)

    if not os.path.exists(file):
        print(f'[WARNING] File: {file} does not exist for date: {datefile}. Skiping...')
        return None

    return file


def get_cmems_multiple_product_day(path_source, org, datehere, dataset_name_file, dataset_name_format_date,
                                   dataset_var_list):
    path_day = path_source
    if org is not None:
        if org == 'YYYYjjj':
            yearstr = datehere.strftime('%Y')
            jjjstr = datehere.strftime('%j')
            path_day = os.path.join(path_source, yearstr, jjjstr)
    strdate = datehere.strftime(dataset_name_format_date)
    ncfiles = []

    for var in dataset_var_list:
        name = dataset_name_file
        var = var.replace('.', '_')
        name = name.replace('$DATE$', strdate)
        name = name.replace('$BAND$', var)
        fname = os.path.join(path_day, name)
        if os.path.exists(fname):
            ncfiles.append(fname)
        else:
            print(f'[WARNING] File: {fname} was not found')

    if len(ncfiles) == len(dataset_var_list):
        return ncfiles
    else:
        return None


def get_satellite_time_from_global_attributes(fproduct):
    dataset = Dataset(fproduct)
    if 'time_coverage_start' in dataset.ncattrs() and 'time_coverage_end' in dataset.ncattrs():
        try:
            start_date = dt.strptime(dataset.time_coverage_start, '%d-%b-%Y %H:%M:%S.%f').replace(tzinfo=pytz.utc)
            end_date = dt.strptime(dataset.time_coverage_end, '%d-%b-%Y %H:%M:%S.%f').replace(tzinfo=pytz.utc)
            sat_time = (start_date.timestamp() + end_date.timestamp()) / 2
            sat_time = dt.utcfromtimestamp(sat_time).replace(tzinfo=pytz.utc)
            return sat_time

        except:
            return None
    dataset.close()
    return None


def download_doors_sources(options):
    section = 'doors_download'
    if not options.has_section(section):
        print(f'[ERROR] Section {section} is not included in  the config file')
        return
    compulsory_options = ['sat_source_dir', 'path_csv', 'dataset_name_file', 'dataset_name_format_date']
    options_dict = {}
    for coption in compulsory_options:
        if not options.has_option(section, coption):
            print(f'[ERROR] Option {coption} is compulsory in section {section}')
            return
        options_dict[coption] = options[section][coption]

    if not os.path.isdir(options_dict['sat_source_dir']):
        print(f'[ERROR] {options_dict["sat_source_dir"]} is not a valid directory')
        return
    path_csv = options_dict['path_csv']
    if not os.path.isfile(path_csv):
        print(f'[ERROR] {path_csv} is not a valid file')
        return
    file_out = os.path.join(os.path.dirname(path_csv), os.path.basename(path_csv)[:-4] + '_sources.csv')
    df = pd.read_csv(path_csv, sep=';')
    nrows = len(df.index)
    df = df.assign(source=[''] * nrows)
    df = df.assign(satellite_time=[''] * nrows)
    for index, row in df.iterrows():
        datehere = dt.strptime(row['Date'], '%Y-%m-%d')
        print(f'[INFO] Date: {datehere}')
        sfile = get_cmems_product_day(options_dict['sat_source_dir'], 'YYYYjjj', datehere,
                                      options_dict['dataset_name_file'], options_dict['dataset_name_format_date'],None)
        if sfile is not None:
            df.loc[index, 'source'] = sfile
            satellite_time = get_satellite_time_from_global_attributes(sfile)
            if satellite_time is not None:
                df.loc[index, 'satellite_time'] = satellite_time.strftime('%Y-%m-%d %H:%M:%S')
    df.to_csv(file_out, sep=';')
    print(f'[INFO] Completed')


def get_cmems_download_options(options):
    section = 'satellite_options'
    options_download = ['download_product', 'download_dataset', 'download_endpoint', 'download_bucket', 'download_tag']
    keys_download = ['product', 'dataset', 'endpoint', 'bucket', 'tag']
    keys_present = [False, False, True, False, False]  ##endpoint is not required, always True
    cmems_donwload_options = {'start_date': None, 'end_date': None}
    for idx, op in enumerate(options_download):
        if options.has_option(section, op):
            cmems_donwload_options[keys_download[idx]] = options[section][op]
            keys_present[idx] = True
        else:
            cmems_donwload_options[keys_download[idx]] = None

    if keys_present.count(True) == 5:
        return cmems_donwload_options
    else:
        return None

def test2():
    print('test2')
    import sat_extract
    import sys
    code_home = os.path.dirname(os.path.dirname(sat_extract.__file__))
    sys.path.append(code_home)
    code_eistools = os.path.join(os.path.dirname(code_home), 'eistools')
    if os.path.exists(code_eistools):
        sys.path.append(code_eistools)
        from cmems_lois import CMEMS_LOIS
        from datetime import datetime as dt
        clois = CMEMS_LOIS(args.verbose)
        cmems_download_options = {
            'start_date': dt(2024,1,1),
            'end_date': dt(2024,1,1),
            'product':'OCEANCOLOUR_GLO_BGC_L3_MY_009_107',
            'dataset':'c3s_obs-oc_glo_bgc-reflectance_my_l3-multi-4km_P1D',
            'bucket':'mdl-native-16',
            'endpoint': None,
            'tag':'202303'
        }
        ods = '%Y/%j'
        clois.make_cmems_download(cmems_download_options, True, '/store3/DOORS', ods, True)
    return True



def main():
    # if test2():
    #     return
    print('[INFO] Creating satellite extracts')
    if not args.config_file:
        return
    if not os.path.exists(args.config_file):
        print(f'[ERROR] File {args.config_file} does not exist')
        return

    options = config_reader(args.config_file)

    if args.download_sources:
        download_doors_sources(options)
        return

    path_output = get_output_path(options)
    if path_output is None:
        print(f'ERROR: {path_output} is not valid')
        return

    if options.has_section('CSV_SELECTION') and options.has_option('CSV_SELECTION', 'path_csv'):
        path_csv = options['CSV_SELECTION']['path_csv']
        if not os.path.exists(path_csv):
            print(f'[ERROR] Path csv: {path_csv} does not exist')
            return
        with open(path_csv) as f:
            first_line = f.readline().strip()
        path_csv_out = f'{path_csv[:-4]}_out.csv'
        fcsv_out = open(path_csv_out, 'w')
        fcsv_out.write(f'{first_line};Extract;Index')

        use_single_file = False
        n_bands = 0
        if options.has_option('CSV_SELECTION', 'use_single_file'):  ##SINGLE FILE SELECTION
            usf = options['CSV_SELECTION']['use_single_file']
            if usf.strip().lower() == 'true' or usf.strip() == '1':
                use_single_file = True

        dataset_name_file = options['CSV_SELECTION']['dataset_name_file']
        dataset_name_format_date = options['CSV_SELECTION']['dataset_name_format_date']
        s = options['CSV_SELECTION']['dataset_var_list']
        dataset_var_list = [x.strip() for x in s.split(',')]
        dataset_var_list_out = dataset_var_list
        if options.has_option('CSV_SELECTION', 'dataset_var_list_out'):
            s = options['CSV_SELECTION']['dataset_var_list_out']
            dataset_var_list_out = [x.strip() for x in s.split(',')]

        if use_single_file:
            rrs_list = []
            rrs_var_list = []
            is_reflectance = False
            if options.has_option('CSV_SELECTION', 'var_rrs_list') and options.has_option('CSV_SELECTION',
                                                                                          'var_rrs_format'):

                var_rrs_list = options['CSV_SELECTION']['var_rrs_list'].strip()
                var_rrs_format = options['CSV_SELECTION']['var_rrs_format'].strip()
                is_reflectance = True
                for r in var_rrs_list.split(','):
                    try:
                        rs = r.strip().replace('.', '_')
                        rrs_here = float(r.strip())
                        rrs_list.append(rrs_here)
                        var_here = var_rrs_format.replace('$BAND$', rs)
                        rrs_var_list.append(var_here)
                    except:
                        is_reflectance = False
                        break
                if is_reflectance:
                    n_bands = len(rrs_list)
        else:
            is_reflectance = True
            for var in dataset_var_list:
                try:
                    float(var)
                except:
                    is_reflectance = False
            if is_reflectance:
                n_bands = len(dataset_var_list)

        if not os.path.exists(path_csv):
            print(f'[ERROR] Path csv {path_csv} was not found')
            return
        try:
            df = pd.read_csv(path_csv, sep=';')
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
        col_time = None
        if options.has_option('CSV_SELECTION', 'col_time'):
            col_time = options['CSV_SELECTION']['col_time']
            format_time = '%H:%M:%S'
        if options.has_option('CSV_SELECTION', 'format_time'):
            format_time = options['CSV_SELECTION']['format_time']
        csv_flags = None
        csv_flags_meanings = None
        if options.has_option('CSV_SELECTION', 'csv_flags'):
            csv_flags_str = options['CSV_SELECTION']['csv_flags']
            csv_flags = [x.strip() for x in csv_flags_str.split(',')]

        path_source, org, wce, time_start, time_stop = get_find_product_info(options)
        size_box = get_box_size(options)
        var_lat, var_lon = get_lat_long_var_names(options)

        extract_list = {}
        if csv_flags is not None:
            csv_flags_meanings = {}

        for idx, row in df.iterrows():
            list = [str(x).strip() for x in row.to_list()]
            line_orig = ";".join(list)
            try:
                datehere = None
                if col_time is not None:
                    try:
                        datetimerow = f'{row[col_date].strip()}T{row[col_time].strip()}'
                        if np.isreal(datetimerow):
                            datetimerow = f'{datetimerow:.0f}'
                        format_datetime = f'{format_date}T{format_time}'
                        datehere = dt.strptime(datetimerow, format_datetime).replace(tzinfo=pytz.utc)
                    except:
                        pass
                if datehere is None:
                    datetimerow = row[col_date]
                    if np.isreal(datetimerow):
                        datetimerow = f'{datetimerow:.0f}'
                    format_datetime = format_date
                    datehere = dt.strptime(datetimerow, format_datetime).replace(tzinfo=pytz.utc)
                    datehere = datehere.replace(hour=12, minute=0, second=0).replace(tzinfo=pytz.utc)
            except:
                print(f'[WARNING] Row {idx} is not valid. Date/Time could not be parsed. Skipping...')
                fcsv_out.write('\n')
                fcsv_out.write(f'{line_orig};NaN;-1')
                continue
            lathere = float(row[col_lat])
            lonhere = float(row[col_lon])
            if np.isnan(lathere) or np.isnan(lonhere):
                print(f'[WARNING] Row {idx} is not valid. Latitude or longitude could not be parsed. Skipping...')
                fcsv_out.write('\n')
                fcsv_out.write(f'{line_orig};NaN;-1')
                continue
            if csv_flags is not None:
                for f in csv_flags:
                    val = row[f].strip()
                    val = val.replace(' ', '_')
                    if f not in csv_flags_meanings.keys():
                        csv_flags_meanings[f] = [val]
                    else:
                        lflags = csv_flags_meanings[f]
                        if val not in lflags:
                            lflags.append(val)
                            csv_flags_meanings[f] = lflags

            if use_single_file:
                cmems_download_options = get_cmems_download_options(options)
                fproduct = get_cmems_product_day(path_source, org, datehere, dataset_name_file,
                                                 dataset_name_format_date,cmems_download_options)
                if fproduct is not None:

                    limits, rc = get_geo_info(options, fproduct, lathere, lonhere)

                    if limits is not None:
                        global_at = get_satellite_global_atrib_from_options(options)
                        datehere_str = datehere.strftime('%Y%m%d')
                        site = f'{get_satellite_ref(global_at)}_{datehere_str}_{rc[0]}_{rc[1]}'
                        # print(f'[INFO] Site: {site}')
                        other = None
                        if csv_flags is not None:
                            other = {}
                            for f in csv_flags:
                                val = row[f].strip()
                                val = val.replace(' ', '_')
                                other[f] = val

                        global_at = add_insitu_global_atrib(global_at, site, lathere, lonhere, other)

                        ofname = os.path.join(path_output, f'extract_{site}.nc')
                        cmems_time = '11:00'
                        if options.has_option('satellite_options', 'satellite_time'):
                            cmems_time = options['satellite_options']['satellite_time'].strip()

                        satellite_time = get_satellite_time_from_global_attributes(fproduct)
                        if satellite_time is None:
                            try:
                                satellite_time = dt.strptime(f'{datehere_str}T{cmems_time}', '%Y%m%dT%H:%M').replace(
                                    tzinfo=pytz.utc)
                            except:
                                print(f'{cmems_time} is not a valid satellite time option. Skipping')
                                fcsv_out.write('\n')
                                fcsv_out.write(f'{line_orig};NaN;-1')
                                continue

                        if site not in extract_list.keys():
                            extract_list[site] = {
                                'ninsitu': 1,
                                'satellite_time': satellite_time,
                                'ofname': ofname,
                                'limits': limits,
                                'file': fproduct,
                                'size_box': size_box,
                                'n_bands': n_bands,
                                'var_lat': var_lat,
                                'var_lon': var_lon,
                                '1': {
                                    'insitu_time': datehere,
                                    'global_at': global_at,
                                }
                            }
                            fcsv_out.write('\n')
                            fcsv_out.write(f'{line_orig};{site}.nc;0')
                        else:
                            idx = extract_list[site]['ninsitu'] + 1
                            fcsv_out.write('\n')
                            fcsv_out.write(f'{line_orig};{site}.nc;{idx - 1}')
                            idxs = str(idx)
                            extract_list[site]['ninsitu'] = idx
                            extract_list[site][idxs] = {
                                'insitu_time': datehere,
                                'global_at': global_at,
                            }
                    else:
                        fcsv_out.write('\n')
                        fcsv_out.write(f'{line_orig};NaN;-1')

                else:
                    print(f'[WARNING] Files not found for date: {datehere.strftime("%Y-%m-%d")}. Skipping...')
                    fcsv_out.write('\n')
                    fcsv_out.write(f'{line_orig};NaN;-1')
            else:

                list_files = get_cmems_multiple_product_day(path_source, org, datehere, dataset_name_file,
                                                            dataset_name_format_date, dataset_var_list)

                if list_files is not None:

                    limits, rc = get_geo_info(options, list_files[0], lathere, lonhere)
                    if limits is not None:
                        global_at = get_satellite_global_atrib_from_options(options)
                        datehere_str = datehere.strftime('%Y%m%d')
                        site = f'{get_satellite_ref(global_at)}_{datehere_str}_{rc[0]}_{rc[1]}'
                        # print(f'[INFO] Site: {site}')
                        other = None
                        if csv_flags is not None:
                            other = {}
                            for f in csv_flags:
                                val = row[f].strip()
                                val = val.replace(' ', '_')
                                other[f] = val

                        global_at = add_insitu_global_atrib(global_at, site, lathere, lonhere, other)

                        ofname = os.path.join(path_output, f'extract_{site}.nc')
                        cmems_time = '11:00'
                        if options.has_option('satellite_options', 'satellite_time'):
                            cmems_time = options['satellite_options']['satellite_time'].strip()
                        try:
                            satellite_time = dt.strptime(f'{datehere_str}T{cmems_time}', '%Y%m%dT%H:%M').astimezone(
                                pytz.utc)
                            print(satellite_time)
                        except:
                            print(f'{cmems_time} is not a valid satellite time option. Skipping')
                            fcsv_out.write('\n')
                            fcsv_out.write(f'{line_orig};NaN;-1')
                            continue

                        if site not in extract_list.keys():
                            extract_list[site] = {
                                'ninsitu': 1,
                                'satellite_time': satellite_time,
                                'ofname': ofname,
                                'limits': limits,
                                'list_files': list_files,
                                'size_box': size_box,
                                'n_bands': n_bands,
                                'var_lat': var_lat,
                                'var_lon': var_lon,
                                '1': {
                                    'insitu_time': datehere,
                                    'global_at': global_at,
                                }
                            }
                            fcsv_out.write('\n')
                            fcsv_out.write(f'{line_orig};{site}.nc;0')
                        else:
                            idx = extract_list[site]['ninsitu'] + 1
                            fcsv_out.write('\n')
                            fcsv_out.write(f'{line_orig};{site}.nc;{idx - 1}')
                            idxs = str(idx)
                            extract_list[site]['ninsitu'] = idx
                            extract_list[site][idxs] = {
                                'insitu_time': datehere,
                                'global_at': global_at,
                            }
                    else:
                        fcsv_out.write('\n')
                        fcsv_out.write(f'{line_orig};NaN;-1')
                else:
                    print(f'[WARNING] Files not found for date: {datehere.strftime("%Y-%m-%d")}. Skipping...')
                    fcsv_out.write('\n')
                    fcsv_out.write(f'{line_orig};NaN;-1')

        fcsv_out.close()
        if use_single_file:
            for site in extract_list:
                extract = extract_list[site]
                ofname = extract['ofname']
                ofname_temp = f'{ofname[:-3]}_temp.nc'
                if os.path.exists(ofname):
                    newExtract = copy_extract(ofname, ofname_temp)
                else:
                    newExtract = start_extract(extract, ofname)

                if is_reflectance:
                    newExtract = add_reflectance_single(newExtract, extract, rrs_list, rrs_var_list)

                newExtract = add_variable_single(newExtract, extract, dataset_var_list, dataset_var_list_out,
                                                 rrs_var_list)

                nidx = 50

                newExtract = add_insitu_basic_info(newExtract, extract, 1, nidx, csv_flags_meanings)

                nhere = extract['ninsitu']
                if nhere > 1:
                    for idx in range(2, nhere + 1):
                        newExtract = add_insitu_basic_info(newExtract, extract, idx, nidx, csv_flags_meanings)
                newExtract.close_file()

                if os.path.exists(ofname_temp):
                    os.rename(ofname_temp, ofname)


        else:  ##MULTIPLE
            for site in extract_list:
                extract = extract_list[site]
                ofname = extract['ofname']
                ofname_temp = f'{ofname[:-3]}_temp.nc'
                if os.path.exists(ofname):
                    newExtract = copy_extract(ofname, ofname_temp)
                else:
                    newExtract = start_extract(extract, ofname)

                if is_reflectance:
                    newExtract = add_reflectance_multiple(newExtract, extract, dataset_var_list)
                else:
                    newExtract = add_variable_multiple(newExtract, extract, dataset_var_list, dataset_var_list_out)
                nidx = 50

                newExtract = add_insitu_basic_info(newExtract, extract, 1, nidx, csv_flags_meanings)

                nhere = extract['ninsitu']
                if nhere > 1:
                    for idx in range(2, nhere + 1):
                        newExtract = add_insitu_basic_info(newExtract, extract, idx, nidx, csv_flags_meanings)
                newExtract.close_file()

                if os.path.exists(ofname_temp):
                    os.rename(ofname_temp, ofname)

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
        run_cmems_option_noreformat(options)
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
