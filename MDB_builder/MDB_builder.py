#!/usr/bin/env python3
# coding: utf-8
"""
Created on Tue Jun 16 12:02:40 2020
Create Matchup Data Base (MDB) file (NetCDF4 file)
@author: javier.concha

Based on EUMETSAT MDB_Builder module (https://ocdb.readthedocs.io/en/latest/ocdb-MDB-user-manual.html)

Run as:

python MDB_builder.py -c path_to_config_file

"""
import shutil

"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
# %% imports
import os
import sys
import subprocess

from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from datetime import datetime
from datetime import timedelta
import pandas as pd
import configparser
from MDB_extra import MDBExtra

# import user defined functions from other .py
code_home = os.path.abspath('../')
sys.path.append(code_home)

import BRDF.brdf_olci as brdf
import COMMON.common_functions as cfs

os.environ['QT_QPA_PLATFORM'] = 'offscreen'  # to avoid error "QXcbConnection: Could not connect to display"
path2ncrcat = '/usr/bin/ncrcat'

import argparse

parser = argparse.ArgumentParser(
    description="Create Match-up DataBase files (MDB) files from satellite extracts and in situ L2 files.")
parser.add_argument("-d", "--debug", help="Debugging mode.", action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-sd', "--startdate", help="The Start Date - format YYYY-MM-DD ")
parser.add_argument('-ed', "--enddate", help="The End Date - format YYYY-MM-DD ")
parser.add_argument('-site', "--sitename", help="Site name.", choices=['VEIT', 'BEFR', 'BSBE'])
parser.add_argument('-ins', "--insitu", help="Satellite sensor name.", choices=['PANTHYR', 'HYPERNETS'])  # ,'HYPSTAR'])
parser.add_argument('-pi', "--path_to_ins", help="Path to in situ sources.")
parser.add_argument('-sat', "--satellite", help="Satellite sensor name.", choices=['OLCI', 'MSI'])
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-ps', "--path_to_sat", help="Path to satellite extracts.")
parser.add_argument('-o', "--output", help="Path to output")
parser.add_argument('-res', "--resolution", help="Resolution OL_2: WRR or WFR (for OLCI)")
parser.add_argument('-nl', "--nolist", help="Do not create satellite and in situ lists.", action="store_true")
parser.add_argument('-nd', "--nodelfiles", help="Do not delete temp files.", action="store_true")

args = parser.parse_args()


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


# user defined functions
def create_list_products(path_source, path_out, wce, res_str, type_product):
    path_to_list = f'{path_out}/file_{type_product}_{res_str}_list.txt'
    if not args.nolist:
        if wce == '*':
            cmd = f'find {path_source} |sort|uniq> {path_to_list}'
        else:
            cmd = f'find {path_source} -name {wce}|sort|uniq> {path_to_list}'
        if args.debug:
            print(f'CMD: {cmd}')
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(err)

    return path_to_list


def copy_nc(ifile, ofile):
    with Dataset(ifile) as src:
        dst = Dataset(ofile, 'w', format='NETCDF4')
        if args.debug:
            print(f'debug: copying {ifile} -----> {ofile}')

            # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)

        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))

        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            dst.createVariable(name, variable.datatype, variable.dimensions)

            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)

            dst[name][:] = src[name][:]
    return dst


def create_variable_from_other_variable(src, dst, name, variable):
    dst.createVariable(name, variable.datatype, variable.dimensions)

    # copy variable attributes all at once via dictionary
    dst[name].setncatts(src[name].__dict__)

    dst[name][:] = src[name][:]

    return dst


def clean_lists(path_out, res_str):
    list_path = f'{path_out}/OL_2_{res_str}_list.txt'
    sort_uniq(list_path, path_out)

    list_path = f'{path_out}/OL_2_{res_str}_list.txt'
    sort_uniq(list_path, path_out)

    if res_str == 'WFR':
        res_L1_str = 'EFR'
    elif res_str == 'WRR':
        res_L1_str = 'ERR'

    list_path = f'{path_out}/OL_1_{res_L1_str}_list.txt'
    sort_uniq(list_path, path_out)


def sort_uniq(list_path, path_out):
    if os.path.exists(list_path):
        cmd = f'sort {list_path}|uniq > {path_out}/temp_list.txt && {path_out}/temp_list.txt {list_path}'
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(err)
        os.remove(f'{path_out}/temp_list.txt')


def create_insitu_list_daily(path_to_insitu_list, date_str):
    # create list in situ per date YYYYMMDD
    YYYYMMDD_str = date_str[:8]
    path_to_list = f'{path_to_insitu_list[:-4]}_{YYYYMMDD_str}.txt'
    cmd = f'cat {path_to_insitu_list}|grep {YYYYMMDD_str}|sort|uniq> {path_to_list}'
    prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    out, err = prog.communicate()
    if err:
        print(err)

    return path_to_list


def add_insitu(extract_path, ofile, path_to_list_daily, datetime_str, time_window, ins_sensor):
    # print(f'Satellite time {datetime_str}')
    date_format = '%Y%m%dT%H%M%S'
    satellite_datetime = datetime.strptime(datetime_str, date_format)

    # to append to nc file
    if args.debug:
        print('debug MDB_builder Line 191: Creating copy of MDB file')

    new_MDB = copy_nc(extract_path, ofile)

    # add time window diff
    if args.debug:
        print('debug MDB_builder Line 197: Adding time window attribute')
    new_MDB.time_diff = f'{time_window * 60 * 60}'  # in seconds

    # create in situ dimensions
    if args.debug:
        print('debug MDB_builder Line 202: Creating in situ dimensions')

    new_MDB.createDimension('insitu_id', 30)
    new_MDB.createDimension('insitu_original_bands', 1603)

    # create variable
    if args.debug:
        print('debug MDB_builder Line 210: Creating in situ and time difference variables')

    insitu_time = new_MDB.createVariable('insitu_time', 'f4', ('satellite_id', 'insitu_id',), zlib=True, complevel=6)
    insitu_time.units = "Seconds since 1970-1-1"
    insitu_time.description = 'In situ time in ISO 8601 format (UTC).'

    insitu_filename = new_MDB.createVariable('insitu_filename', 'S2', ('satellite_id', 'insitu_id'), zlib=True,
                                             complevel=6)
    insitu_filename.description = 'In situ filename.'

    insitu_filepath = new_MDB.createVariable('insitu_filepath', 'S2', ('satellite_id', 'insitu_id'), zlib=True,
                                             complevel=6)
    insitu_filepath.description = 'In situ file path.'

    insitu_original_bands = new_MDB.createVariable('insitu_original_bands', 'f4', ('insitu_original_bands'),
                                                   fill_value=-999, zlib=True, complevel=6)
    insitu_original_bands.description = 'In situ bands in nm.'

    insitu_Rrs = new_MDB.createVariable('insitu_Rrs', 'f4', ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                        fill_value=-999, zlib=True, complevel=6)
    insitu_Rrs.description = 'In situ Rrs'

    time_difference = new_MDB.createVariable('time_difference', 'f4', ('satellite_id', 'insitu_id'), fill_value=-999,
                                             zlib=True, complevel=6)
    time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
    time_difference.units = "seconds"

    insitu_idx = 0
    # extract in situ data
    if args.debug:
        print('debug MDB_builder Line 240: Starting extraction...')
    with open(path_to_list_daily) as file:
        for idx, line in enumerate(file):
            # HYPERNETS ex: HYPERNETS_W_BEFR_L2A_REF_202103151201_202103231711_v1.1.nc
            # PANTHYR ex: AAOT_20190923_019046_20200923_104022_AZI_225_data.csv
            # time in ISO 8601 format (UTC). Ex: "2020-01-05T09:27:27.934965Z"
            if ins_sensor == 'PANTHYR':
                date_str = os.path.basename(line[:-1]).split('_')[3]
                time_str = os.path.basename(line[:-1]).split('_')[4]
            elif ins_sensor == 'HYPERNETS':
                date_str = os.path.basename(line[:-1]).split('_')[5][0:8]
                time_str = os.path.basename(line[:-1]).split('_')[5][-4:] + '00'

            YYYY_str = date_str[:4]
            MM_str = date_str[4:6]
            DD_str = date_str[6:8]
            HH_str = time_str[:2]
            mm_str = time_str[2:4]
            ss_str = time_str[4:6]
            insitu_datetime_str = f'{YYYY_str}-{MM_str}-{DD_str}T{HH_str}:{mm_str}:{ss_str}Z'
            insitu_datetime = datetime(int(YYYY_str), int(MM_str), int(DD_str), int(HH_str), int(mm_str), int(ss_str))

            time_diff = (insitu_datetime - satellite_datetime).total_seconds() / (60 * 60)

            if args.debug:
                print(f'debug MDB_builder Line 264. Time diff: {time_diff} Time Window {time_window}')

            if np.abs(time_diff) <= time_window:
                if args.debug:
                    print(f'debug: In situ found for: {line}')
                insitu_time[0, insitu_idx] = float(datetime.strptime(insitu_datetime_str,
                                                                     "%Y-%m-%dT%H:%M:%SZ").timestamp())  # Ex: 2021-02-24T11:31:00Z
                insitu_filename[0, insitu_idx] = os.path.basename(line[:-1])
                insitu_filepath[0, insitu_idx] = line[:-1]
                time_difference[0, insitu_idx] = float(time_diff) * 60 * 60  # in seconds

                if ins_sensor == 'PANTHYR':
                    # get data from csv using pandas
                    data = pd.read_csv(line[:-1], parse_dates=['timestamp'])

                    if insitu_idx == 0:
                        wl0 = data['wavelength'].tolist()
                        insitu_original_bands[:] = wl0
                    insitu_rhow_vec = [x for x, in data['rhow'][:]]
                    insitu_RrsArray = ma.array(insitu_rhow_vec).transpose() / np.pi
                    # insitu_Rrs[0, :, insitu_idx] = [ma.array(insitu_rhow_vec).transpose()]
                    insitu_Rrs[0, :, insitu_idx] = [insitu_RrsArray]
                    insitu_idx += 1
                    # print(rhow0)
                elif ins_sensor == 'HYPERNETS':
                    nc_ins = Dataset(line[:-1], 'r')
                    if insitu_idx == 0:
                        insitu_original_bands[:] = nc_ins.variables['wavelength'][:].tolist()
                        # ins_water_leaving_radiance = nc_ins.variables['water_leaving_radiance'][:]

                    insitu_rhow_vec = [x for x, in nc_ins.variables['reflectance'][:]]

                    insitu_RrsArray = ma.array(insitu_rhow_vec).transpose() / np.pi
                    # insitu_Rrs[0, :, insitu_idx] = [ma.array(insitu_rhow_vec).transpose()]
                    insitu_Rrs[0, :, insitu_idx] = [insitu_RrsArray]
                    insitu_idx += 1
                    nc_ins.close()
    new_MDB.close()
    if insitu_idx == 0:
        if os.path.exists(ofile):
            os.remove(ofile)
        if args.debug:
            print('Not in situ measurements within the time window. MDB file deleted!')
        return False
    else:
        return True


def add_insitu_aeronet(extract_path, ofile, areader, satellite_datetime, time_window, mdb_secondary):
    satellite_date = satellite_datetime.strftime('%Y-%m-%d')
    check_date = areader.prepare_data_fordate(satellite_date)
    if not check_date:
        print(f'[WARNING] No in situ spectra were found for date: {satellite_date}')
        return

    time_list = areader.extract_time_list()
    nspectra = len(time_list)
    nspectra_within_timewindow = 0
    spectra_within_timewindow = np.zeros(nspectra, dtype=bool)
    time_dif_seconds = np.zeros(nspectra, dtype=float)
    time_window_seconds = time_window * 3600
    for i in range(nspectra):
        time_dif_seconds[i] = float(abs((time_list[i] - satellite_datetime).total_seconds()))
        if time_dif_seconds[i] < time_window_seconds:
            spectra_within_timewindow[i] = True
            nspectra_within_timewindow = nspectra_within_timewindow + 1

    if args.verbose:
        print(f'[INFO] Number of in-situ spectra on date {satellite_date}:{nspectra}')
        print(f'[INFO] --> within the time window: {nspectra_within_timewindow}')

    if nspectra_within_timewindow == 0:
        print(f'[WARNING] No in situ spectra were found for date: {satellite_date} within the specified time window')
        return

    rrs = areader.extract_rrs(False)
    exactwl = areader.extract_spectral_data('Exact_Wavelengths', False)

    # check nc secondary
    ncsecondary = None
    if mdb_secondary is not None:
        dinput = Dataset(extract_path)
        sat_time = datetime.fromtimestamp(float(dinput.variables['satellite_time'][0]))
        # print(extract_path, sat_time, '=================================================================')
        lat_array = ma.array(dinput.variables['satellite_latitude'][:])
        lon_array = ma.array(dinput.variables['satellite_longitude'][:])

        ncsecondary = mdb_secondary.get_extradaset(dinput.satellite, dinput.platform, sat_time, lat_array, lon_array)

        dinput.close()
        if ncsecondary is None:
            print(f'[WARNING] No extra data found for date: {satellite_date}')
            return False

    # to append to nc file
    new_MDB = copy_nc(extract_path, ofile)

    # add varibles from mdb secondary
    if not ncsecondary is None:
        for varname in mdb_secondary.variables:
            var_secondary = ncsecondary.variables[varname]
            new_MDB = create_variable_from_other_variable(ncsecondary, new_MDB, varname, var_secondary)

    # add time window diff
    new_MDB.time_diff = f'{time_window_seconds}'  # in seconds

    n_insitu_bands = areader.nwl

    # create in situ dimensions
    new_MDB.createDimension('insitu_id', 50)
    new_MDB.createDimension('insitu_original_bands', n_insitu_bands)
    # new_MDB.createDimension('insitu_Rrs_bands', None)

    # create variable
    insitu_time = new_MDB.createVariable('insitu_time', 'f8', ('satellite_id', 'insitu_id',), zlib=True, complevel=6)
    insitu_time.units = "Seconds since 1970-1-1"
    insitu_time.description = 'In situ time in ISO 8601 format (UTC).'

    insitu_original_bands = new_MDB.createVariable('insitu_original_bands', 'f4', ('insitu_original_bands'),
                                                   fill_value=-999, zlib=True, complevel=6)
    insitu_original_bands.description = 'In situ bands in nm'

    insitu_Rrs = new_MDB.createVariable('insitu_Rrs', 'f4', ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                        fill_value=-999, zlib=True, complevel=6)
    insitu_Rrs.description = 'In situ Rrs'

    insitu_exact_wavelenghts = new_MDB.createVariable('insitu_exact_wavelenghts', 'f4',
                                                      ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                                      fill_value=-999, zlib=True, complevel=6)
    insitu_exact_wavelenghts.description = 'In situ bands in nm'

    time_difference = new_MDB.createVariable('time_difference', 'f4', ('satellite_id', 'insitu_id'), fill_value=-999,
                                             zlib=True, complevel=6)
    time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
    time_difference.units = "seconds"

    insitu_original_bands[:] = areader.dataset['Nominal_Wavelenghts'][:]

    insitu_idx = 0
    for i in range(nspectra):
        if not spectra_within_timewindow[i]:
            continue

        insitu_time[0, insitu_idx] = float(time_list[i].timestamp())
        time_difference[0, insitu_idx] = time_dif_seconds[i]

        rrs_here = np.array(rrs[i][:])
        rrs_here[rrs[i][:].mask] = -999
        insitu_Rrs[0, :, insitu_idx] = rrs_here

        wl_here = np.array(exactwl[i][:])
        wl_here[exactwl[i][:].mask] = -999
        insitu_exact_wavelenghts[0, :, insitu_idx] = wl_here

        insitu_idx = insitu_idx + 1

    new_MDB.close()
    if insitu_idx == 0:
        if os.path.exists(ofile):
            os.remove(ofile)
            if args.verbose:
                print('Not in situ measurements within the time window. MDB file deleted!')
        return False
    else:
        return True


def add_insitu_resto(extract_path, ofile, insitu_dataset, time_list, satellite_datetime, time_window):
    satellite_date = satellite_datetime.strftime('%Y-%m-%d')
    insitu_start_time = datetime.strptime(insitu_dataset.start_date, '%Y-%m-%d %H:%M')
    insitu_end_time = datetime.strptime(insitu_dataset.end_date, '%Y-%m-%d %H:%M')
    if satellite_datetime < insitu_start_time or satellite_datetime > insitu_end_time:
        print(f'[WARNING] No in situ spectra were found for date: {satellite_datetime}')
        return False
    # dtref = datetime(1970, 1, 1, 0, 0, 0).replace(microsecond=0)
    # start_search = satellite_datetime
    # sec = float((satellite_datetime. - self.dtref).total_seconds())

    nspectra = len(time_list)
    nspectra_within_timewindow = 0
    spectra_within_timewindow = np.zeros(nspectra, dtype=bool)
    time_dif_seconds = np.zeros(nspectra, dtype=float)
    time_window_seconds = time_window * 3600
    for i in range(nspectra):
        time_dif_seconds[i] = float(abs((time_list[i] - satellite_datetime).total_seconds()))
        if time_dif_seconds[i] < time_window_seconds:
            spectra_within_timewindow[i] = True
            nspectra_within_timewindow = nspectra_within_timewindow + 1

    if args.verbose:
        print(f'[INFO] Number of in-situ spectra on date {satellite_date}:{nspectra}')
        print(f'[INFO] --> within the time window: {nspectra_within_timewindow}')

    if nspectra_within_timewindow == 0:
        print(f'[WARNING] No in situ spectra were found for date: {satellite_date} within the specified time window')
        return False

    # sat extract nc file
    new_MDB = copy_nc(extract_path, ofile)

    # add time window diff
    new_MDB.time_diff = f'{time_window_seconds}'  # in seconds

    # n insitu bands
    nw = insitu_dataset.variables['Nominal_Wavelenghts'][:]
    n_insitu_bands = len(nw)

    # create in situ dimensions
    new_MDB.createDimension('insitu_id', 50)
    new_MDB.createDimension('insitu_original_bands', n_insitu_bands)

    # create variables
    insitu_time = new_MDB.createVariable('insitu_time', 'f8', ('satellite_id', 'insitu_id',), zlib=True, complevel=6)
    insitu_time.units = "Seconds since 1970-01-01"
    insitu_time.description = 'In situ time in ISO 8601 format (UTC).'

    insitu_original_bands = new_MDB.createVariable('insitu_original_bands', 'f4', ('insitu_original_bands'),
                                                   fill_value=-999, zlib=True, complevel=6)
    insitu_original_bands.description = 'In situ bands in nm'

    insitu_Rrs = new_MDB.createVariable('insitu_Rrs', 'f4', ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                        fill_value=-999, zlib=True, complevel=6)
    insitu_Rrs.description = 'In situ Rrs'

    # insitu_exact_wavelenghts = new_MDB.createVariable('insitu_exact_wavelenghts', 'f4',
    #                                                   ('satellite_id', 'insitu_original_bands', 'insitu_id'),
    #                                                   fill_value=-999, zlib=True, complevel=6)
    # insitu_exact_wavelenghts.description = 'In situ bands in nm'

    time_difference = new_MDB.createVariable('time_difference', 'f4', ('satellite_id', 'insitu_id'), fill_value=-999,
                                             zlib=True, complevel=6)
    time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
    time_difference.units = "seconds"

    insitu_original_bands[:] = nw

    insitu_idx = 0
    for i in range(nspectra):
        if not spectra_within_timewindow[i]:
            continue

        insitu_time[0, insitu_idx] = float(time_list[i].timestamp())
        time_difference[0, insitu_idx] = time_dif_seconds[i]

        rrs_here = np.ma.array(insitu_dataset.variables['RRS'][i][:])
        rrs_here[rrs_here.mask] = -999
        insitu_Rrs[0, :, insitu_idx] = rrs_here

        # wl_here = np.array(exactwl[i][:])
        # wl_here[exactwl[i][:].mask] = -999
        # insitu_exact_wavelenghts[0, :, insitu_idx] = wl_here

        insitu_idx = insitu_idx + 1

    new_MDB.close()

    return True


def add_insitu_meda(extract_path, ofile, path_to_list_daily, datetime_str, time_window, ins_sensor):
    # print(f'Satellite time {datetime_str}')
    date_format = '%Y%m%dT%H%M%S'
    satellite_datetime = datetime.strptime(datetime_str, date_format)

    # to append to nc file
    if args.debug:
        print('debug MDB_builder Line 191: Creating copy of MDB file')

    new_MDB = copy_nc(extract_path, ofile)

    # add time window diff
    if args.debug:
        print('debug MDB_builder Line 197: Adding time window attribute')
    new_MDB.time_diff = f'{time_window * 60 * 60}'  # in seconds

    # create in situ dimensions
    if args.debug:
        print('debug MDB_builder Line 202: Creating in situ dimensions')

    new_MDB.createDimension('insitu_id', 30)
    new_MDB.createDimension('insitu_original_bands', 7)

    # create variable
    if args.debug:
        print('debug MDB_builder Line 210: Creating in situ and time difference variables')

    insitu_time = new_MDB.createVariable('insitu_time', 'f4', ('satellite_id', 'insitu_id',), zlib=True, complevel=6)
    insitu_time.units = "Seconds since 1970-1-1"
    insitu_time.description = 'In situ time in ISO 8601 format (UTC).'

    insitu_filename = new_MDB.createVariable('insitu_filename', 'S2', ('satellite_id', 'insitu_id'), zlib=True,
                                             complevel=6)
    insitu_filename.description = 'In situ filename.'

    insitu_filepath = new_MDB.createVariable('insitu_filepath', 'S2', ('satellite_id', 'insitu_id'), zlib=True,
                                             complevel=6)
    insitu_filepath.description = 'In situ file path.'

    insitu_original_bands = new_MDB.createVariable('insitu_original_bands', 'f4', ('insitu_original_bands'),
                                                   fill_value=-999, zlib=True, complevel=6)
    insitu_original_bands.description = 'In situ bands in nm.'

    insitu_Rrs = new_MDB.createVariable('insitu_Rrs', 'f4', ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                        fill_value=-999, zlib=True, complevel=6)
    insitu_Rrs.description = 'In situ Rrs'

    time_difference = new_MDB.createVariable('time_difference', 'f4', ('satellite_id', 'insitu_id'), fill_value=-999,
                                             zlib=True, complevel=6)
    time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
    time_difference.units = "seconds"

    insitu_idx = 0
    # extract in situ data
    if args.debug:
        print('debug MDB_builder Line 240: Starting extraction...')
    with open(path_to_list_daily) as file:
        for idx, line in enumerate(file):

            ins_path = line[:-1]
            ins_filename = ins_path.split('/')[-1]
            nc_ins = Dataset(ins_path, 'r')
            ins_date = datetime.strptime(ins_filename.split('_')[3],'%y%m%d').replace(hour=0,minute=0,seconds=0,microsecond=0)

            ins_hours = nc_ins.variables['timetag'][:]
            for ihour in range(len(ins_hours)):
                ins_hour = ins_hours[idx]
                insitu_datetime = ins_date + timedelta(hours=ins_hour)
                time_diff = (insitu_datetime - satellite_datetime).total_seconds() / (60 * 60)
                if args.debug:
                    print(f'debug MDB_builder Line 592. Time diff: {time_diff} Time Window {time_window}')
                if np.abs(time_diff) <= time_window:
                    insitu_time[0, insitu_idx] = float(insitu_datetime.timestamp())  # Ex: 2021-02-24T11:31:00Z
                    insitu_filename[0, insitu_idx] = os.path.basename(line[:-1])
                    insitu_filepath[0, insitu_idx] = line[:-1]
                    time_difference[0, insitu_idx] = float(time_diff) * 60 * 60  # in seconds
                    insitu_RrsArray  = np.array(nc_ins.variables['rrs'][ihour,:])
                    insitu_Rrs[0, :, insitu_idx] = [insitu_RrsArray]
                    insitu_idx += 1
            nc_ins.close()
    new_MDB.close()
    if insitu_idx == 0:
        if os.path.exists(ofile):
            os.remove(ofile)
        if args.debug:
            print('Not in situ measurements within the time window. MDB file deleted!')
        return False
    else:
        return True


def get_simple_fileextracts_list(path_extract, wce, datetime_start, datetime_stop):
    list_files = []
    for f in os.listdir(path_extract):
        path_file = os.path.join(path_extract, f)
        include = False
        if wce is not None and f.find(wce) >= 0 and f.endswith('.nc'):
            include = True
        if wce is None and f.endswith('.nc'):
            include = True
        if datetime_start is not None and datetime_stop is not None and include:
            include = False
            nc_extract = Dataset(path_file)
            datehere = datetime.fromtimestamp(float(nc_extract.variables['satellite_time'][0]))
            if datetime_start <= datehere <= datetime_stop:
                include = True
            nc_extract.close()

        # if include:
        #     dataset_here = Dataset(path_file)
        #     if 'satellite_longitude' not in dataset_here.variables:
        #         print(f'[WARNING] file structure is not correct for: {f}')
        #         include = False
        #     dataset_here.close()
        if include:
            if args.verbose:
                print(f'[INFO] Appending file: {f}')
            list_files.append(path_file)
    return list_files


def get_start_end_date_from_options(options):
    if options['Time_and_sites_selection']['time_start']:
        datetime_start = datetime.strptime(options['Time_and_sites_selection']['time_start'], '%Y-%m-%d')
    else:
        datetime_start = datetime.strptime('2000-01-01', '%Y-%m-%d')
    if options['Time_and_sites_selection']['time_stop']:
        datetime_end = datetime.strptime(options['Time_and_sites_selection']['time_stop'], '%Y-%m-%d') + timedelta(
            seconds=59, minutes=59, hours=23)
    else:
        datetime_end = datetime.today()

    return datetime_start, datetime_end


def make_simple_builder(options, path_extract, path_out):
    wce = None
    satellite = 'SATELLITE'
    platform = 'PLATFORM'
    if options.has_option('satellite_options', 'satellite') and options.has_option('satellite_options', 'platform'):
        satellite = options['satellite_options']['satellite']
        platform = options['satellite_options']['platform']
        if platform.find(
                ',') > 0:  # if platform is a coma-separated list (e.g. 'A,B'), it's not included in the wildcard
            platform = platform.replace(',', '_')
            wce = f'{satellite}'
        else:
            wce = f'{satellite}{platform}'

    ac = 'AC'
    if options.has_option('satellite_options', 'ac'):
        ac = options['satellite_options']['ac']
    sensor = 'SENSOR'
    if options.has_option('satellite_options', 'sensor'):
        sensor = options['satellite_options']['sensor']

    site_type = options['Time_and_sites_selection']['insitu_type']

    datetime_start, datetime_stop = get_start_end_date_from_options(options)
    datetime_start_str = datetime_start.strftime('%Y%m%d')
    datetime_stop_str = datetime_stop.strftime('%Y%m%d')
    name_out = f'MDB_{satellite.upper()}{platform.upper()}_{sensor.upper()}_{ac.upper()}_{site_type.upper()}_{datetime_start_str}_{datetime_stop_str}.nc'
    ncout_file = os.path.join(path_out, name_out)
    if args.verbose:
        print(f'[INFO] Out file: {ncout_file}')

    list_files = get_simple_fileextracts_list(path_extract, wce, datetime_start, datetime_stop)

    if len(list_files) == 0:
        print(f'[WARNING] No sat extract files were found. Please review params in config file')
        return

    if len(list_files) > 100:
        if args.verbose:
            print(f'[INFO] Preparing contatenation of {len(list_files)} files...')
        list_files_tmp = []
        for icent in range(0, len(list_files), 100):
            if args.verbose:
                print(f'[INFO] Concatening: {icent} / {len(list_files)}')
            indextmp = int(icent / 100)
            list_files_here = list_files[icent:icent + 100]
            ncout_file_tmp = os.path.join(path_out, f'Temp_{indextmp}.nc')
            list_files_tmp.append(ncout_file_tmp)
            list_files_here.append(ncout_file_tmp)
            # concatenation
            cmd = [f"ncrcat -O -h"] + list_files_here
            cmd = " ".join(cmd)
            prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            out, err = prog.communicate()
            if err:
                print(f'[ERROR]{err}')
        list_files_tmp.append(ncout_file)
        cmd = [f"ncrcat -O -h"] + list_files_tmp
        cmd = " ".join(cmd)
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]{err}')

        [os.remove(f) for f in list_files_tmp[:-1]]


    else:
        list_files.append(ncout_file)
        # concatenation
        cmd = [f"ncrcat -O -h"] + list_files
        cmd = " ".join(cmd)
        if args.debug:
            print(f'CMD="{cmd}"')
        # os.system(cmd)
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]{err}')

    print(f'Concatenated file created: {ncout_file}')


def concatenate_nc_impl(list_files, path_out, ncout_file):
    if len(list_files) == 0:
        print(f'[WARNING] No sat extract files were found. Please review params in config file')
        return

    nfiles_ref = 100
    if len(list_files) > nfiles_ref:
        if args.verbose:
            print(f'[INFO] Preparing contatenation of {len(list_files)} files...')
        list_files_tmp = []
        for icent in range(0, len(list_files), nfiles_ref):
            if args.verbose:
                print(f'[INFO] Concatening: {icent} / {len(list_files)}')
            indextmp = int(icent / nfiles_ref)
            list_files_here = list_files[icent:icent + nfiles_ref]
            # if icent == 10:
            #     for tal in list_files_here:
            #         print(tal)
            #     print(
            #         '????????????????????????????????????????????????????????????????????????????????????????????????')
            ncout_file_tmp = os.path.join(path_out, f'Temp_{indextmp}.nc')
            list_files_tmp.append(ncout_file_tmp)
            list_files_here.append(ncout_file_tmp)
            # concatenation
            cmd = [f"ncrcat -O -h"] + list_files_here
            cmd = " ".join(cmd)
            prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            out, err = prog.communicate()
            if err:
                print(f'[ERROR]{err}')
        list_files_tmp.append(ncout_file)
        cmd = [f"ncrcat -O -h"] + list_files_tmp
        cmd = " ".join(cmd)
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]{err}')

        # [os.remove(f) for f in list_files_tmp[:-1]]
        if not args.nodelfiles:
            [os.remove(f) for f in list_files]

    else:
        list_files.append(ncout_file)
        # concatenation
        cmd = [f"ncrcat -O -h"] + list_files
        cmd = " ".join(cmd)
        if args.debug:
            print(f'CMD="{cmd}"')
        # os.system(cmd)
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]{err}')

        if not args.nodelfiles:
            [os.remove(f) for f in list_files[:-1]]
    print(f'Concatenated file created: {ncout_file}')


def check_single_mdbfile(mdbfile):
    valid = True
    # sizegood = 6 #multi
    sizegood = 11  # olci
    nc = Dataset(mdbfile)
    if nc.dimensions['satellite_bands'].size != sizegood:
        valid = False
    nc.close()
    return valid


def check_single_mdbfile_exist(prename, postname, list_mdbfiles_pathout):
    file_prev = None
    for name in list_mdbfiles_pathout:
        if name.startswith(prename) and name.endswith(postname):
            file_prev = name
            break
    return file_prev


def get_time_list_from_resto_dataset(insitu_dataset):
    time_list = []
    time_array = np.array(insitu_dataset.variables['Time'][:])
    for time in time_array:
        timehere = datetime.fromtimestamp(float(time))
        time_list.append(timehere)
    return time_list


def check():
    # base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/Irbe_Lighthouse/WFR/extracts'
    # for name in os.listdir(base):
    #     if not name.endswith('nc'):
    #         continue
    #     extract_path = os.path.join(base, name)
    #     dinput = Dataset(extract_path)
    #     sat_time = datetime.fromtimestamp(float(dinput.variables['satellite_time'][0]))
    #     print(name, sat_time, '=================================================================')
    # return True
    base = '/mnt/c/DATA_LUIS/OCTAC_WORK/MED_MATCHUPS/EXTRACTS/OLCI'
    pathextracts = os.path.join(base, 'Casablanca_Platform')
    for name in os.listdir(base):
        # if  name.endswith('extract_Casablanca_Platform.nc'):
        #     print(name)
        if not name.endswith('AERONET_Casablanca_Platform.nc'):
            continue
        ltal = name.split('_')
        datestr = ltal[3]
        datehere = datetime.strptime(datestr, '%Y%m%dT%H%M%S')
        daten = datehere.strftime('%Y%j')
        namen = f'CMEMS2_O{daten}-rrs-med-fr_nc_extract_Casablanca_Platform.nc'
        fextract = os.path.join(pathextracts, namen)
        fextractout = os.path.join(base, namen)
        # print(fextract,fextractout)
        os.replace(fextract, fextractout)
        print(namen)
    return True


# #############################
# %%
def main():
    print('[INFO] Creating MDB files!')
    ##CALLING ONLY FOR PRE-TESTING###
    # b = check()
    # if b:
    #     return
    #################################
    if args.debug:
        print('[DEBUG] Entering Debugging Mode:')

    # load config file (if exists)
    if args.config_file:
        if os.path.isfile(args.config_file):
            options = config_reader(args.config_file)

    # path to output
    if args.output:
        path_out = args.output
    elif args.config_file:
        if options['file_path']['output_dir']:
            path_out = options['file_path']['output_dir']
    if not os.path.isdir(path_out):
        os.mkdir(path_out)
    if args.verbose:
        print(f'[INFO] Path to output: {path_out}')

    list_mdbfiles_pathout = []
    for name in os.listdir(path_out):
        list_mdbfiles_pathout.append(name)

    # path to satellite extract
    if args.path_to_sat:
        satellite_path_source = args.path_to_sat
    elif args.config_file:
        if options['file_path']['sat_extract_dir']:
            satellite_path_source = options['file_path']['sat_extract_dir']
    if args.verbose:
        print(f'[INFO] Path to satellite sources: {satellite_path_source}')

    ##SYKE and INSITU MODE, SIMPLE BUILDER WITHOUT ADDING IN SITU REFLECTANCE DATA
    if options.has_option('Time_and_sites_selection', 'insitu_type'):
        if options['Time_and_sites_selection']['insitu_type'] == 'SYKE':
            print(f'[INFO] Entering simple concatenation mode...')
            make_simple_builder(options, satellite_path_source, path_out)
            return
        if options['Time_and_sites_selection']['insitu_type'] == 'INSITU':
            print(f'[INFO] Entering simple concatenation mode...')
            make_simple_builder(options, satellite_path_source, path_out)
            return

    # satellite, platform, sensor
    if args.satellite:
        sat_sensor = args.satellite
    elif args.config_file:
        sat_sensor = options['satellite_options']['sensor']
        sat_satellite = options['satellite_options']['satellite']
        sat_platform = options['satellite_options']['platform']
    if args.verbose:
        print(f'[INFO] Satellite: {sat_satellite.upper()}')
        print(f'[INFO] Satellite sensor: {sat_sensor.upper()}')
        print(f'[INFO] Satellite platform: {sat_platform.upper()}')

    # resolution
    if args.resolution:
        res = args.resolution
    elif args.config_file:
        res = options['satellite_options']['resolution']

    # atmospheric correction
    atm_corr = 'STANDARD'
    if args.config_file:
        if options.has_option('satellite_options', 'ac'):
            atm_corr = options['satellite_options']['ac']

    # in situ site
    if args.sitename:
        station_name = args.sitename
    elif args.config_file:
        station_name = options['Time_and_sites_selection']['sites']
    in_situ_lat, in_situ_lon = cfs.get_lat_lon_ins(station_name)  # in situ location based on the station name
    if args.verbose:
        print(f'[INFO] Station name: {station_name} with lat: {in_situ_lat}, long: {in_situ_lon}')

    # wild card expression for searching extracts
    if sat_sensor.upper() == 'OLCI' and atm_corr == 'STANDARD':
        wce = f'"{sat_satellite}{sat_platform}*OL_2_{res}*{station_name}*"'  # wild card expression
    else:
        wce = f'"{sat_satellite}{sat_platform}*nc"'

    if args.verbose:
        print(f'[INFO] Satellite extract Wild Card Expression: {wce}')

    # searching for file extracts
    path_to_satellite_list = create_list_products(satellite_path_source, path_out, wce, res, atm_corr)
    if args.verbose:
        print(f'[INFO] Path to satellite extract list: {path_to_satellite_list}')

    # in situ sensor: PANTHYR, HYPERNETS, AERONET, RESTO
    if args.insitu:
        ins_sensor = args.insitu
    elif args.config_file:
        if options['Time_and_sites_selection']['insitu_type']:
            ins_sensor = options['Time_and_sites_selection']['insitu_type']

    # wild card expression for in situ data
    if ins_sensor == 'PANTHYR':
        wce = f'"*AZI_270_data.csv"'  # wild card expression
    elif ins_sensor == 'HYPERNETS':  # 'HYPSTAR':
        wce = f'"HYPERNETS_W_*{station_name}*L2A_REF*_v1.2.nc"'
    elif ins_sensor == 'AERONET':
        wce = f'*{station_name}*'
    elif ins_sensor == 'RESTO':
        wce = f'RESTO_{station_name}'
    elif ins_sensor == 'MEDA':
        wce = f'meda_lam_opt_*_L1v1.nc'
    if args.debug:
        print(f'[DEBUG] In Situ Wild Card Expression: {wce}')
    # in situ path source
    if args.path_to_ins:
        insitu_path_source = args.path_to_ins
    elif args.config_file:
        if options['file_path']['ins_source_dir']:
            insitu_path_source = options['file_path']['ins_source_dir']
    if args.verbose:
        print(f'[INFO] Path to in situ data: {insitu_path_source}')
    # sarching for in situ data
    if ins_sensor != 'RESTO':
        path_to_insitu_list = create_list_products(insitu_path_source, path_out, wce, res, 'insitu')
    if args.debug:
        print(f'[DEBUG] Path to in situ list: {path_to_insitu_list}')
    areader = None
    if ins_sensor == 'AERONET':
        f = open(path_to_insitu_list)
        file_aeronet = f.readline()[:-1]
        f.close()
        aeronet_path = '/home/lois/PycharmProjects/aeronet'
        sys.path.append(aeronet_path)
        from base.anet_nc_reader import AERONETReader
        areader = AERONETReader(file_aeronet)
        if args.debug:
            print(f'[DEBUG] Path to AERONET NC file: {file_aeronet}')

    path_to_insitu = None
    if ins_sensor == 'RESTO':
        if args.config_file:
            if options.has_option('Time_and_sites_selection', 'name_file'):
                path_to_insitu = os.path.join(insitu_path_source, options['Time_and_sites_selection']['name_file'])
        if path_to_insitu is None:
            print(f'[ERROR] RESTO in situ path {path_to_insitu} was not defined')
            return
        if not os.path.exists(path_to_insitu):
            print(f'[ERROR] RESTO in situ path {path_to_insitu} does not exist')
            return
        if args.verbose:
            print(f'[INFO] Path to RESTO file: {path_to_insitu}')

    # time dif between in situ and sat data
    time_window = 3  # in hours (+- hours)
    if args.config_file:
        if options.has_option('Time_and_sites_selection', 'time_window'):
            time_window = int(options['Time_and_sites_selection']['time_window'])
    if args.verbose:
        print(f'[INFO] Time window: {time_window} hours')

    file_list = []  # for the concatenation later

    # defining start and end date
    if args.startdate:
        datetime_start = datetime.strptime(args.startdate, '%Y-%m-%d')
    elif args.config_file:
        if options['Time_and_sites_selection']['time_start']:
            datetime_start = datetime.strptime(options['Time_and_sites_selection']['time_start'], '%Y-%m-%d')
        else:
            datetime_start = datetime.strptime('2000-01-01', '%Y-%m-%d')
    else:
        datetime_start = datetime.strptime('2000-01-01', '%Y-%m-%d')
    if args.enddate:
        datetime_end = datetime.strptime(args.enddate, '%Y-%m-%d') + timedelta(seconds=59, minutes=59, hours=23)
    elif args.config_file:
        if options['Time_and_sites_selection']['time_stop']:
            datetime_end = datetime.strptime(options['Time_and_sites_selection']['time_stop'], '%Y-%m-%d') + timedelta(
                seconds=59, minutes=59, hours=23)
        else:
            datetime_end = datetime.today()
    else:
        datetime_end = datetime.today()

    if args.verbose:
        print(f'[INFO] Start date: {datetime_start}')
        print(f'[INFO] End date: {datetime_end}')

    # defining secondary extracts (if available)
    mdb_secondary = None
    if args.config_file:
        if options.has_option('file_path', 'sat_extract_secondary') and options.has_option('file_path',
                                                                                           'sat_extract_secondary_variables'):
            path_secondary = options['file_path']['sat_extract_secondary'].strip()
            if os.path.exists(path_secondary):
                if args.verbose:
                    print(f'[INFO] Started MDB secondary from path: {path_secondary}')
                mdb_secondary = MDBExtra(path_secondary, False)
                var_names = options['file_path']['sat_extract_secondary_variables'].split(',')
                for idx in range(len(var_names)):
                    var_names[idx] = var_names[idx].strip()
                mdb_secondary.variables = var_names

    with open(path_to_satellite_list, 'r') as file:
        for cnt, line in enumerate(file):
            extract_path = line[:-1]
            # extract date time info
            sensor_str = extract_path.split('/')[-1].split('_')[0]
            valid_extract = True
            if atm_corr == 'ACOLITE':
                res_str = 'EFR'
                lpath = extract_path.split('/')[-1].split('_')
                datetime_here = datetime(int(lpath[2]), int(lpath[3]), int(lpath[4]), int(lpath[5]), int(lpath[6]),
                                         int(lpath[7]))
                datetime_str = datetime_here.strftime('%Y%m%dT%H%M%S')
            elif atm_corr == 'CCI' or atm_corr == 'MULTI' or atm_corr == 'OLCI-L3':
                res_str = res
                nc_sat = Dataset(extract_path)
                if not 'satellite_time' in nc_sat.variables:
                    valid_extract = False
                    nc_sat.close()
                else:
                    datetime_here = datetime.fromtimestamp(float(nc_sat.variables['satellite_time'][0]))
                    if options.has_option('Time_and_sites_selection', 'time_sat_default'):
                        hm = options['Time_and_sites_selection']['time_sat_default']
                        dhm = datetime.strptime(hm, '%H:%M')
                        datetime_here = datetime_here.replace(hour=dhm.hour, minute=dhm.minute)
                    else:
                        datetime_here = datetime_here.replace(hour=11)
                    datetime_str = datetime_here.strftime('%Y%m%dT%H%M%S')
                    nc_sat.close()
                    sensor_str = atm_corr
            else:
                res_str = extract_path.split('/')[-1].split('_')[3]
                datetime_str = extract_path.split('/')[-1].split('_')[7]

            if not valid_extract:
                if args.verbose:
                    print('-----------------')
                    print(f'[WARNING] Extract {extract_path}is not valid. Skipping...')
                continue
            if args.verbose:
                print('-----------------')
                print(f'[INFO] Date: {datetime_str} Satellite/Platform: {sensor_str} Resolution: {res_str}')
                print(f'[INFO] Extract_path: {extract_path}')
            date_format = '%Y%m%dT%H%M%S'
            satellite_datetime = datetime.strptime(datetime_str, date_format)
            datetime_creation = datetime.today().strftime(date_format)
            if datetime_start <= satellite_datetime <= datetime_end:
                try:
                    if ins_sensor == 'AERONET':
                        path_to_list_daily = None
                        prefilename = f'MDB_{sensor_str}_{res_str}_{datetime_str}'
                        postfilename = f'{ins_sensor}_{station_name}.nc'
                        # print(prefilename, postfilename)
                        filename_prev = check_single_mdbfile_exist(prefilename, postfilename, list_mdbfiles_pathout)
                        if not filename_prev is None:
                            ofile = os.path.join(path_out, filename_prev)
                            if args.verbose:
                                print(f'[INFO] File already created: {ofile}')
                            if check_single_mdbfile(ofile):
                                file_list.append(ofile)  # for ncrcat later
                            else:
                                print(f'[WARNING] File {ofile} is not valid')
                                filecopy = os.path.join('/mnt/c/DATA_LUIS/OCTAC_WORK/MED_MATCHUPS/EXTRACTS/OLCI',
                                                        filename_prev)
                                shutil.copy(ofile, filecopy)
                        else:
                            filename = f'{prefilename}_{datetime_creation}_{postfilename}'
                            ofile = os.path.join(path_out, filename)
                            print(f'[INFO] Creating file from extract: {extract_path}')
                            if add_insitu_aeronet(extract_path, ofile, areader, satellite_datetime, time_window,
                                                  mdb_secondary):
                                print(f'[INFO] File created: {ofile}')
                                file_list.append(ofile)  # for ncrcat later

                    elif ins_sensor == 'RESTO':
                        path_to_list_daily = None
                        insitu_dataset = Dataset(path_to_insitu)
                        resto_time_list = get_time_list_from_resto_dataset(insitu_dataset)
                        prefilename = f'MDB_{sensor_str}_{res_str}_{datetime_str}'
                        postfilename = f'{ins_sensor}_{station_name}.nc'
                        # print(prefilename, postfilename)
                        filename_prev = check_single_mdbfile_exist(prefilename, postfilename, list_mdbfiles_pathout)
                        if not filename_prev is None:
                            ofile = os.path.join(path_out, filename_prev)
                            if args.verbose:
                                print(f'[INFO] File already created: {ofile}')
                            if check_single_mdbfile(ofile):
                                file_list.append(ofile)  # for ncrcat later
                        else:
                            filename = f'{prefilename}_{datetime_creation}_{postfilename}'
                            ofile = os.path.join(path_out, filename)
                            print(f'[INFO] Creating file from extract: {extract_path}')
                            if add_insitu_resto(extract_path, ofile, insitu_dataset, resto_time_list,
                                                satellite_datetime, time_window):
                                print(f'[INFO] File created: {ofile}')
                                file_list.append(ofile)  # for ncrcat later
                    elif ins_sensor == 'MEDA':
                        prefilename = f'MDB_{sensor_str}_{res_str}_{datetime_str}'
                        postfilename = f'{ins_sensor}_{station_name}.nc'
                        # print(prefilename, postfilename)
                        filename_prev = check_single_mdbfile_exist(prefilename, postfilename, list_mdbfiles_pathout)
                        if not filename_prev is None:
                            ofile = os.path.join(path_out, filename_prev)
                            if args.verbose:
                                print(f'[INFO] File already created: {ofile}')
                            if check_single_mdbfile(ofile):
                                file_list.append(ofile)  # for ncrcat later
                        else:
                            filename = f'{prefilename}_{datetime_creation}_{postfilename}'
                            ofile = os.path.join(path_out, filename)
                            print(f'[INFO] Creating file from extract (MEDA): {extract_path}')
                            path_to_list_daily = create_insitu_list_daily(path_to_insitu_list, datetime_str)
                            if add_insitu_meda(extract_path, ofile, path_to_list_daily, datetime_str, time_window,
                                          ins_sensor):
                                print(f'[INFO] file created: {ofile}')
                                file_list.append(ofile)  # for ncrcat later


                    else:
                        path_to_list_daily = create_insitu_list_daily(path_to_insitu_list, datetime_str)
                        if not os.stat(path_to_list_daily).st_size == 0:  # no PANTHYR data or not for that angle
                            filename = f'MDB_{sensor_str}_{res_str}_{datetime_str}_{datetime_creation}_{ins_sensor}_{station_name}.nc'
                            ofile = os.path.join(path_out, filename)
                            if os.path.exists(ofile):
                                if args.verbose:
                                    print(f'[INFO] File already created: {ofile}')
                                file_list.append(ofile)  # for ncrcat later
                            else:
                                if add_insitu(extract_path, ofile, path_to_list_daily, datetime_str, time_window,
                                              ins_sensor):
                                    print(f'[INFO] file created: {ofile}')
                                    file_list.append(ofile)  # for ncrcat later

                        else:
                            if args.verbose:
                                print('[WARNING] No in situ measurements found!')

                # except:
                except Exception as e:
                    if args.verbose:
                        print(f'[ERROR] Exception: {e}')
                    pass

                if path_to_list_daily is not None and os.path.exists(path_to_list_daily):
                    os.remove(path_to_list_daily)
            else:
                if args.verbose:
                    print('[WARNING] Out of time frame.')

    level_prod = 'L2'

    # calling subprocess for concatanating ncdf files # # ncrcat -h MDB_S3*.nc outcat2.nN

    ncout_file = os.path.join(path_out,
                              f'MDB_{sat_satellite}{sat_platform}_{sat_sensor.upper()}_{res_str}_{atm_corr}_{level_prod}_{ins_sensor}_{station_name}.nc')

    concatenate_nc_impl(file_list, path_out, ncout_file)
    # file_list.append(ncout_file)
    # # concatenation
    # cmd = [f"ncrcat -O -h"] + file_list
    # cmd = " ".join(cmd)
    # if args.verbose:
    #     print(f'CMD="{cmd}"')
    # os.system(cmd)
    # if not args.nodelfiles:
    #     [os.remove(f) for f in file_list[:-1]]

    print(f'[INFO]Concatenated file created: {ncout_file}')


# %%
if __name__ == '__main__':
    main()
