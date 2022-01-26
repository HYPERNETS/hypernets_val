#!/usr/bin/env python3
# coding: utf-8
"""
Created on Tue Jul 8 12:02:40 2021
Create extract from OLCI data as a NetCDF4 file
@author: javier.concha

Based on EUMETSAT MDB_Builder module (https://ocdb.readthedocs.io/en/latest/ocdb-MDB-user-manual.html)

Run as:

python MDB_builder.py -c path_to_config_file

"""
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
import math
from datetime import datetime
from datetime import timedelta
import configparser

# import user defined functions from other .py
code_home = os.path.abspath('../')
sys.path.append(code_home)

import BRDF.brdf_olci as brdf
import COMMON.common_functions as cfs

os.environ['QT_QPA_PLATFORM'] = 'offscreen'  # to avoid error "QXcbConnection: Could not connect to display"
path2ncrcat = '/opt/local/bin/ncrcat'

import argparse

parser = argparse.ArgumentParser(description="Create Match-up DataBase files (MDB) files.")
parser.add_argument("-d", "--debug", help="Debugging mode.", action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-sd', "--startdate", help="The Start Date - format YYYY-MM-DD ")
parser.add_argument('-ed', "--enddate", help="The End Date - format YYYY-MM-DD ")
parser.add_argument('-site', "--sitename", help="Site name.", choices=['VEIT', 'BEFR', 'BSBE'])
parser.add_argument('-sat', "--satellite", help="Satellite sensor name.", choices=['OLCI', 'MSI'])
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-ps', "--path_to_sat", help="Path to satellite sources.")
parser.add_argument('-o', "--output", help="Path to output")
parser.add_argument('-res', "--resolution", help="Resolution OL_2: WRR or WFR (for OLCI)")
parser.add_argument('-nl', "--nolist", help="Do not create satellite and in situ lists.", action="store_true")

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
def create_list_products(path_source, path_out, wce, res_str, type_product, dt_start, dt_end, org):
    if args.nolist:
        return None

    path_to_list = f'{path_out}/file_{type_product}_{res_str}_list.txt'

    if os.path.exists(path_to_list):
        print('Deleting previous file list...')
        cmd = f'rm {path_to_list}'
        print(cmd)
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(err)

    if org is None:
        cmd = f'find {path_source} -name {wce}|sort|uniq> {path_to_list}'
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(err)

    if org == 'yyyy_jjj':
        dt = dt_start
        while dt <= dt_end:
            print('Date: ', dt)
            year = dt.strftime('%Y')
            jday = dt.strftime('%j')
            path_day = os.path.join(path_source, year, jday)
            cmd = f'find {path_day} -name {wce}|sort|uniq>> {path_to_list}'
            prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            out, err = prog.communicate()
            if err:
                print(err)
            dt = dt + timedelta(hours=24)

    return path_to_list


def get_val_from_tie_point_grid(yPoint, xPoint, ySubsampling, xSubsampling, dataset):
    grid_height = dataset.shape[0]
    grid_width = dataset.shape[1]
    fi = (xPoint + 0.5) / xSubsampling
    fj = (yPoint + 0.5) / ySubsampling
    i0 = floor_and_crop(fi, 0, grid_width - 2)
    j0 = floor_and_crop(fj, 0, grid_height - 2)
    i1 = i0 + 1
    j1 = j0 + 1
    wi = fi - i0
    wj = fj - j0
    x00 = dataset[j0, i0]
    x10 = dataset[j0, i1]
    x01 = dataset[j1, i0]
    x11 = dataset[j1, i1]
    val = x00 + (wi * (x10 - x00)) + (wj * (x01 - x00)) + (wi * wj * (x11 + x00 - x01 - x10))
    return val


def floor_and_crop(v, minV, maxV):
    rv = math.floor(v)
    if rv < minV:
        return minV
    if rv > maxV:
        return maxV
    return rv


def extract_wind_and_angles(path_source, in_situ_lat, in_situ_lon):  # for OLCI
    # from Tie-Points grid (a coarser grid)
    filepah = os.path.join(path_source, 'tie_geo_coordinates.nc')
    nc_sat = Dataset(filepah, 'r')
    tie_lon = nc_sat.variables['longitude'][:]
    tie_lat = nc_sat.variables['latitude'][:]
    nc_sat.close()

    filepah = os.path.join(path_source, 'tie_meteo.nc')
    nc_sat = Dataset(filepah, 'r')
    horizontal_wind = nc_sat.variables['horizontal_wind'][:]
    nc_sat.close()

    filepah = os.path.join(path_source, 'tie_geometries.nc')
    nc_sat = Dataset(filepah, 'r')
    SZA = nc_sat.variables['SZA'][:]
    SAA = nc_sat.variables['SAA'][:]
    OZA = nc_sat.variables['OZA'][:]
    OAA = nc_sat.variables['OAA'][:]
    nc_sat.close()

    r, c = cfs.find_row_column_from_lat_lon(tie_lat, tie_lon, in_situ_lat, in_situ_lon)

    ws0 = horizontal_wind[r, c, 0]
    ws1 = horizontal_wind[r, c, 1]
    sza = SZA[r, c]
    saa = SAA[r, c]
    vza = OZA[r, c]
    vaa = OAA[r, c]

    return ws0, ws1, sza, saa, vza, vaa


def launch_create_extract(in_situ_sites, size_box, path_source, res_str, make_brdf):
    for site in in_situ_sites:
        try:
            in_situ_lat = float(in_situ_sites[site]['latitude'])
            in_situ_lon = float(in_situ_sites[site]['longitude'])
            path_output = in_situ_sites[site]['path_out']
            if not os.path.exists(path_output):
                os.mkdir(path_output)
            print(f'Creating extract for site: {site}')
            extract_path = create_extract(size_box, site, path_source, path_output, in_situ_lat, in_situ_lon, res_str,
                                          make_brdf)
            if not extract_path is None:
                print(f'file created: {extract_path}')
        except Exception as e:
            if args.debug:
                print(f'Exception: {e}')
                pass


def create_extract(size_box, station_name, path_source, path_output, in_situ_lat, in_situ_lon, res_str, make_brdf):
    if args.verbose:
        print(f'Creating extract for {station_name} from {path_source}')
    ofname = None
    # extract IFP-OL-2 version
    with open(os.path.join(path_source, 'xfdumanifest.xml'), 'r', encoding="utf-8") as read_obj:
        check_version = False
        for line in read_obj:
            if 'IPF-OL-2' in line and check_version == False:
                IPF_OL_2_version = line.split('"')[3]
                proc_version_str = f'IPF-OL-2 version {IPF_OL_2_version}'
                if args.verbose:
                    print(proc_version_str)
                check_version = True
                pass

    # open nc file
    coordinates_filename = 'geo_coordinates.nc'

    rhow_0400p00_filename = 'Oa01_reflectance.nc'
    rhow_0412p50_filename = 'Oa02_reflectance.nc'
    rhow_0442p50_filename = 'Oa03_reflectance.nc'
    rhow_0490p00_filename = 'Oa04_reflectance.nc'
    rhow_0510p00_filename = 'Oa05_reflectance.nc'
    rhow_0560p00_filename = 'Oa06_reflectance.nc'
    rhow_0620p00_filename = 'Oa07_reflectance.nc'
    rhow_0665p00_filename = 'Oa08_reflectance.nc'
    rhow_0673p75_filename = 'Oa09_reflectance.nc'
    rhow_0681p25_filename = 'Oa10_reflectance.nc'
    rhow_0708p75_filename = 'Oa11_reflectance.nc'
    rhow_0753p75_filename = 'Oa12_reflectance.nc'
    rhow_0778p75_filename = 'Oa16_reflectance.nc'
    rhow_0865p00_filename = 'Oa17_reflectance.nc'
    rhow_0885p00_filename = 'Oa18_reflectance.nc'
    rhow_1020p50_filename = 'Oa21_reflectance.nc'

    AOT_0865p50_filename = 'w_aer.nc'
    WQSF_filename = 'wqsf.nc'

    # reading latitude and longitude
    filepah = os.path.join(path_source, coordinates_filename)
    nc_sat = Dataset(filepah, 'r')
    lat = nc_sat.variables['latitude'][:, :]
    lon = nc_sat.variables['longitude'][:, :]
    contain_flag = cfs.contain_location(lat, lon, in_situ_lat, in_situ_lon)

    if contain_flag == 1:
        if not args.verbose:
            print('-----------------')
        r, c = cfs.find_row_column_from_lat_lon(lat, lon, in_situ_lat, in_situ_lon)
        # coordenadas son y, x

        start_idx_x = (c - int(size_box / 2))
        stop_idx_x = (c + int(size_box / 2) + 1)
        start_idx_y = (r - int(size_box / 2))
        stop_idx_y = (r + int(size_box / 2) + 1)

        if start_idx_y >= 0 and (stop_idx_y + 1) < lat.shape[0] and start_idx_x >= 0 and (stop_idx_x + 1) < lat.shape[
            1]:

            # print(r, c, in_situ_lat, in_situ_lon)
            # read variables
            if args.verbose:
                print('Reading variables...')

            filepah = os.path.join(path_source, rhow_0400p00_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0400p00 = nc_sat.variables['Oa01_reflectance'][:]
            satellite_start_time = nc_sat.start_time  # reading here start time
            # satellite_stop_time = nc_sat.stop_time

            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0412p50_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0412p50 = nc_sat.variables['Oa02_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0442p50_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0442p50 = nc_sat.variables['Oa03_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0490p00_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0490p00 = nc_sat.variables['Oa04_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0510p00_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0510p00 = nc_sat.variables['Oa05_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0560p00_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0560p00 = nc_sat.variables['Oa06_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0620p00_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0620p00 = nc_sat.variables['Oa07_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0665p00_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0665p00 = nc_sat.variables['Oa08_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0673p75_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0673p75 = nc_sat.variables['Oa09_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0681p25_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0681p25 = nc_sat.variables['Oa10_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0708p75_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0708p75 = nc_sat.variables['Oa11_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0753p75_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0753p75 = nc_sat.variables['Oa12_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0778p75_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0778p75 = nc_sat.variables['Oa16_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0865p00_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0865p00 = nc_sat.variables['Oa17_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_0885p00_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_0885p00 = nc_sat.variables['Oa18_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, rhow_1020p50_filename)
            nc_sat = Dataset(filepah, 'r')
            rhow_1020p50 = nc_sat.variables['Oa21_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, AOT_0865p50_filename)
            nc_sat = Dataset(filepah, 'r')
            AOT_0865p50 = nc_sat.variables['T865'][:]
            nc_sat.close()

            filepah = os.path.join(path_source, WQSF_filename)
            nc_sat = Dataset(filepah, 'r')
            WQSF = nc_sat.variables['WQSF'][:]
            WQSF_flag_masks = nc_sat.variables['WQSF'].flag_masks
            WQSF_flag_meanings = nc_sat.variables['WQSF'].flag_meanings
            nc_sat.close()

            # Calculate BRDF (it uses CHL_OC4ME_extract)
            if make_brdf:
                if args.verbose:
                    print('Making BRDF...')
                filepah = os.path.join(path_source, 'chl_oc4me.nc')
                nc_sat1 = Dataset(filepah, 'r')
                CHL_OC4ME = nc_sat1.variables['CHL_OC4ME'][:]
                nc_sat1.close()
                CHL_OC4ME_extract = ma.array(CHL_OC4ME[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])
                CHL_OC4ME_extract[~CHL_OC4ME_extract.mask] = ma.power(10, CHL_OC4ME_extract[~CHL_OC4ME_extract.mask])
                mask_chl = np.copy(CHL_OC4ME_extract.mask)
                CHL_OC4ME_extract[mask_chl] = 1  ##temporal value
                ws0, ws1, sza, saa, vza, vaa = extract_wind_and_angles(path_source, in_situ_lat, in_situ_lon)
                BRDF0 = np.full(CHL_OC4ME_extract.shape, np.nan)
                BRDF1 = np.full(CHL_OC4ME_extract.shape, np.nan)
                BRDF2 = np.full(CHL_OC4ME_extract.shape, np.nan)
                BRDF3 = np.full(CHL_OC4ME_extract.shape, np.nan)
                BRDF4 = np.full(CHL_OC4ME_extract.shape, np.nan)
                BRDF5 = np.full(CHL_OC4ME_extract.shape, np.nan)
                BRDF6 = np.full(CHL_OC4ME_extract.shape, np.nan)
                for ind0 in range(CHL_OC4ME_extract.shape[0]):
                    for ind1 in range(CHL_OC4ME_extract.shape[1]):
                        chl = CHL_OC4ME_extract[ind0, ind1]
                        # 412.5, 442.5, 490, 510, 560, 620, 660 bands
                        # 0      1      2    3    4    5    6   brdf index
                        # 412.5  442.5  490  510  560  620  665 OLCI bands
                        # 02     03     04   05   06   07   08  OLCI band names in L2
                        brdf_coeffs = brdf.brdf(ws0, ws1, chl, sza, saa, vza, vaa)
                        BRDF0[ind0, ind1] = brdf_coeffs[0, 0]
                        BRDF1[ind0, ind1] = brdf_coeffs[0, 1]
                        BRDF2[ind0, ind1] = brdf_coeffs[0, 2]
                        BRDF3[ind0, ind1] = brdf_coeffs[0, 3]
                        BRDF4[ind0, ind1] = brdf_coeffs[0, 4]
                        BRDF5[ind0, ind1] = brdf_coeffs[0, 5]
                        BRDF6[ind0, ind1] = brdf_coeffs[0, 6]
                BRDF0[mask_chl] = -999
                BRDF1[mask_chl] = -999
                BRDF2[mask_chl] = -999
                BRDF3[mask_chl] = -999
                BRDF4[mask_chl] = -999
                BRDF5[mask_chl] = -999
                BRDF6[mask_chl] = -999
                CHL_OC4ME_extract[mask_chl] = -999  # fill value

            # %% Save extract as netCDF4 file
            filename = path_source.split('/')[-1].replace('.', '_') + '_extract_' + station_name + '.nc'
            ofname = os.path.join(path_output, filename)

            satellite = filename[0:2]
            platform = filename[2]
            sensor = 'olci'

            if os.path.exists(ofname):
                os.remove(ofname)

            if args.verbose:
                print('Starting new extract...')

            new_EXTRACT = Dataset(ofname, 'w', format='NETCDF4')

            # Atributes
            new_EXTRACT.creation_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
            new_EXTRACT.satellite = satellite
            new_EXTRACT.platform = platform
            new_EXTRACT.sensor = sensor
            new_EXTRACT.description = f'{satellite}{platform} {sensor.upper()} {res_str} L2 extract'
            # new_EXTRACT.satellite_start_time = nc_sat.start_time
            # new_EXTRACT.satellite_stop_time = nc_sat.stop_time    
            # new_EXTRACT.satellite_PDU = path_source.split('/')[-1]
            # new_EXTRACT.satellite_path_source = path_source
            new_EXTRACT.satellite_aco_processor = 'Atmospheric Correction processor: xxx'
            new_EXTRACT.satellite_proc_version = proc_version_str

            # new_EXTRACT.datapolicy = 'Notice to users: Add data policy'
            # new_EXTRACT.insitu_sensor_processor_version = '0.0'
            new_EXTRACT.insitu_site_name = station_name

            new_EXTRACT.insitu_lat = in_situ_lat
            new_EXTRACT.insitu_lon = in_situ_lon

            if make_brdf:
                new_EXTRACT.satellite_ws0 = ws0
                new_EXTRACT.satellite_ws1 = ws1
                new_EXTRACT.satellite_SZA_center_pixel = sza
                new_EXTRACT.satellite_SAA_center_pixel = saa
                new_EXTRACT.satellite_VZA_center_pixel = vza
                new_EXTRACT.satellite_VAA_center_pixel = vaa

            # dimensions
            new_EXTRACT.createDimension('satellite_id', None)
            new_EXTRACT.createDimension('rows', size_box)
            new_EXTRACT.createDimension('columns', size_box)
            new_EXTRACT.createDimension('satellite_bands', 16)
            if make_brdf:
                new_EXTRACT.createDimension('satellite_BRDF_bands', 7)

            # geometry variables
            if args.verbose:
                print('Creating geometry variables...')
            filepah = os.path.join(path_source, 'tie_geometries.nc')
            nc_sat = Dataset(filepah, 'r')
            xsubsampling = nc_sat.getncattr('ac_subsampling_factor')
            ysubsampling = nc_sat.getncattr('al_subsampling_factor')
            SZA = nc_sat.variables['SZA'][:]
            SAA = nc_sat.variables['SAA'][:]
            OZA = nc_sat.variables['OZA'][:]
            OAA = nc_sat.variables['OAA'][:]
            nc_sat.close()
            SZAO = new_EXTRACT.createVariable('SZA', 'f4', ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                              zlib=True, complevel=6)
            SZAO.units = 'degress'
            SZAO.long_name = 'Sun Zenith Angle'

            SAAO = new_EXTRACT.createVariable('SAA', 'f4', ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                              zlib=True, complevel=6)
            SAAO.units = 'degress'
            SAAO.long_name = 'Sun Azimuth Angle'

            OZAO = new_EXTRACT.createVariable('OZA', 'f4', ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                              zlib=True, complevel=6)
            OZAO.units = 'degress'
            OZAO.long_name = 'Observation Zenith Angle'

            OAAO = new_EXTRACT.createVariable('OAA', 'f4', ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                              zlib=True, complevel=6)
            OAAO.units = 'degress'
            OAAO.long_name = 'Observation Azimuth Angle'

            for yy in range(size_box):
                for xx in range(size_box):
                    yPos = start_idx_y + yy
                    xPos = start_idx_x + xx
                    SZAO[0, yy, xx] = get_val_from_tie_point_grid(yPos, xPos, ysubsampling, xsubsampling, SZA)
                    SAAO[0, yy, xx] = get_val_from_tie_point_grid(yPos, xPos, ysubsampling, xsubsampling, SAA)
                    OZAO[0, yy, xx] = get_val_from_tie_point_grid(yPos, xPos, ysubsampling, xsubsampling, OZA)
                    OAAO[0, yy, xx] = get_val_from_tie_point_grid(yPos, xPos, ysubsampling, xsubsampling, OAA)

            # time
            if args.verbose:
                print('Adding time and PDU...')
            satellite_time = new_EXTRACT.createVariable('satellite_time', 'f4', ('satellite_id'), fill_value=-999,
                                                        zlib=True, complevel=6)
            satellite_time[0] = float(datetime.strptime(satellite_start_time, "%Y-%m-%dT%H:%M:%S.%fZ").timestamp())
            satellite_time.units = "Seconds since 1970-1-1"

            # PDU
            satellite_PDU = new_EXTRACT.createVariable('satellite_PDU', 'S2', ('satellite_id'), zlib=True,
                                                       complevel=6)  # string
            satellite_PDU[0] = path_source.split('/')[-1]
            satellite_PDU.long_name = "OLCI source PDU name"

            # latitude
            if args.verbose:
                print('Adding latitude and longitude...')
            satellite_latitude = new_EXTRACT.createVariable('satellite_latitude', 'f8',
                                                            ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                                            zlib=True, complevel=6)
            satellite_latitude[0, :, :] = [lat[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]]
            satellite_latitude.short_name = 'latitude'

            # longitude
            satellite_longitude = new_EXTRACT.createVariable('satellite_longitude', 'f8',
                                                             ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                                             zlib=True, complevel=6)
            satellite_longitude[0, :, :] = [lon[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]]
            satellite_longitude.short_name = 'longitude'

            # Variable satellite_bands (wavelenghts for Rrs)
            if args.verbose:
                print('Adding satellite bands and reflectances...')
            satellite_bands = new_EXTRACT.createVariable('satellite_bands', 'f4', ('satellite_bands'), fill_value=-999,
                                                         zlib=True, complevel=6)
            satellite_bands[:] = [0400.00, 0412.50, 0442.50, 0490.00, 0510.00, 0560.00, 0620.00, 0665.00, 0673.75,
                                  0681.25, 0708.75, 0753.75, 0778.75, 0865.00, 0885.00, 1020.50]
            satellite_bands.units = 'nm'

            # Variable satellite_Rrs (NOT BRDF-corrected remote sensing reflectance)
            satellite_Rrs = new_EXTRACT.createVariable('satellite_Rrs', 'f4',
                                                       ('satellite_id', 'satellite_bands', 'rows', 'columns'),
                                                       fill_value=-999, zlib=True, complevel=6)
            satellite_Rrs[0, 0, :, :] = ma.array(rhow_0400p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 1, :, :] = ma.array(rhow_0412p50[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 2, :, :] = ma.array(rhow_0442p50[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 3, :, :] = ma.array(rhow_0490p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 4, :, :] = ma.array(rhow_0510p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 5, :, :] = ma.array(rhow_0560p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 6, :, :] = ma.array(rhow_0620p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 7, :, :] = ma.array(rhow_0665p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 8, :, :] = ma.array(rhow_0673p75[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 9, :, :] = ma.array(rhow_0681p25[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 10, :, :] = ma.array(rhow_0708p75[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 11, :, :] = ma.array(rhow_0753p75[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 12, :, :] = ma.array(rhow_0778p75[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 13, :, :] = ma.array(rhow_0865p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 14, :, :] = ma.array(rhow_0885p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs[0, 15, :, :] = ma.array(rhow_1020p50[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) / np.pi
            satellite_Rrs.short_name = 'Satellite Rrs'
            satellite_Rrs.long_name = "Above water Remote Sensing Reflectance for OLCI acquisition"
            satellite_Rrs.units = "sr-1"

            # BRDF-corrected varibles: satellite_BRDF_bancs, satellite_BRDF_Rrs, satellite_BRDF_fQ, chl_oc4me
            if make_brdf:
                # double satellite_BRDF_bands     (satellite_BRDF_bands) ;
                satellite_BRDF_bands = new_EXTRACT.createVariable('satellite_BRDF_bands', 'f4',
                                                                  ('satellite_BRDF_bands'),
                                                                  fill_value=-999, zlib=True, complevel=6)
                satellite_BRDF_bands[:] = [412.50, 442.50, 490.00, 510.00, 560.00, 620.00, 665.00]
                satellite_BRDF_bands.units = 'nm'

                satellite_BRDF_Rrs = new_EXTRACT.createVariable('satellite_BRDF_Rrs', 'f4',
                                                                ('satellite_id', 'satellite_BRDF_bands', 'rows',
                                                                 'columns'),
                                                                fill_value=-999, zlib=True, complevel=6)
                BRDF_Rrs = ma.array(rhow_0412p50[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x] * BRDF0) / np.pi
                BRDF_Rrs[mask_chl] = -999
                satellite_BRDF_Rrs[0, 0, :, :] = BRDF_Rrs
                BRDF_Rrs = ma.array(rhow_0442p50[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x] * BRDF1) / np.pi
                BRDF_Rrs[mask_chl] = -999
                satellite_BRDF_Rrs[0, 1, :, :] = BRDF_Rrs
                BRDF_Rrs = ma.array(rhow_0490p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x] * BRDF2) / np.pi
                BRDF_Rrs[mask_chl] = -999
                satellite_BRDF_Rrs[0, 2, :, :] = BRDF_Rrs
                BRDF_Rrs = ma.array(rhow_0510p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x] * BRDF3) / np.pi
                BRDF_Rrs[mask_chl] = -999
                satellite_BRDF_Rrs[0, 3, :, :] = BRDF_Rrs
                BRDF_Rrs = ma.array(rhow_0560p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x] * BRDF4) / np.pi
                BRDF_Rrs[mask_chl] = -999
                satellite_BRDF_Rrs[0, 4, :, :] = BRDF_Rrs
                BRDF_Rrs = ma.array(rhow_0620p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x] * BRDF5) / np.pi
                BRDF_Rrs[mask_chl] = -999
                satellite_BRDF_Rrs[0, 5, :, :] = BRDF_Rrs
                BRDF_Rrs = ma.array(rhow_0665p00[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x] * BRDF6) / np.pi
                BRDF_Rrs[mask_chl] = -999
                satellite_BRDF_Rrs[0, 6, :, :] = BRDF_Rrs
                satellite_BRDF_Rrs.description = 'Satellite Rrs BRDF-corrected'
                satellite_BRDF_Rrs.short_name = 'Satellite Rrs.'
                satellite_BRDF_Rrs.long_name = "Above water Remote Sensing Reflectance for OLCI acquisition with BRDF correction applied";
                satellite_BRDF_Rrs.units = "sr-1"
                satellite_BRDF_fQ = new_EXTRACT.createVariable('satellite_BRDF_fQ', 'f4',
                                                               ('satellite_id', 'satellite_BRDF_bands', 'rows',
                                                                'columns'),
                                                               fill_value=-999, zlib=True, complevel=6)
                satellite_BRDF_fQ[0, 0, :, :] = [ma.array(BRDF0)]
                satellite_BRDF_fQ[0, 1, :, :] = [ma.array(BRDF1)]
                satellite_BRDF_fQ[0, 2, :, :] = [ma.array(BRDF2)]
                satellite_BRDF_fQ[0, 3, :, :] = [ma.array(BRDF3)]
                satellite_BRDF_fQ[0, 4, :, :] = [ma.array(BRDF4)]
                satellite_BRDF_fQ[0, 5, :, :] = [ma.array(BRDF5)]
                satellite_BRDF_fQ[0, 6, :, :] = [ma.array(BRDF6)]
                satellite_BRDF_fQ.description = 'Satellite BRDF fQ coefficients'
                satellite_chl_oc4me = new_EXTRACT.createVariable('chl_oc4me', 'f4', ('satellite_id', 'rows', 'columns'),
                                                                 fill_value=-999, zlib=True, complevel=6)
                satellite_chl_oc4me[0, :, :] = ma.array(CHL_OC4ME_extract)
                satellite_chl_oc4me.description = 'Satellite Chlorophyll-a concentration from OC4ME.'

            # AOT
            satellite_AOT_0865p50_box = new_EXTRACT.createVariable('satellite_AOT_0865p50', 'f4',
                                                                   ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                                                   zlib=True, complevel=6)
            satellite_AOT_0865p50_box[0, :, :] = ma.array(AOT_0865p50[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])
            satellite_AOT_0865p50_box.description = 'Satellite Aerosol optical thickness'

            # WQSF: Quality Flags
            satellite_WQSF = new_EXTRACT.createVariable('satellite_WQSF', 'f4', ('satellite_id', 'rows', 'columns'),
                                                        fill_value=-999, zlib=True, complevel=6)
            satellite_WQSF[0, :, :] = [ma.array(WQSF[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])]
            satellite_WQSF.description = 'Satellite Level 2 WATER Product, Classification, Quality and Science Flags Data Set'
            satellite_WQSF.flag_masks = WQSF_flag_masks
            satellite_WQSF.flag_meanings = WQSF_flag_meanings

            new_EXTRACT.close()
        else:
            if args.verbose:
                print('Warning: Index out of bound!')
    else:
        if args.verbose:
            print('Warning: File does NOT contains the in situ location!')

    return ofname


#############################
# %%
def main():
    print('Creating satellite extracts.')

    if args.debug:
        print('Entering Debugging Mode:')
    # load config file
    if args.config_file:
        if os.path.isfile(args.config_file) == True:
            options = config_reader(args.config_file)
    # else:
    #     print(args.config_file + ' does not exist. Please provide a valid config file path')
    #     sys.exit()    

    # path to satellite source
    if args.path_to_sat:
        satellite_path_source = args.path_to_sat
    elif args.config_file:
        if not options.has_option('file_path', 'sat_source_dir'):
            print(
                f'Error: Section file_path, Option sat_source_dir is not available in the config file: {args.config_file}')
            return
        if options['file_path']['sat_source_dir']:
            satellite_path_source = options['file_path']['sat_source_dir']

    if not os.path.exists(satellite_path_source):
        print(f'ERROR path: {satellite_path_source} does not exit')
        return

    if args.verbose:
        print(f'Path to satellite sources: {satellite_path_source}')

    # temporay path
    tmp_path = None
    if args.config_file:
        if options.has_option('file_path', 'tmp_dir') and options['file_path']['tmp_dir']:
            tmp_path = options['file_path']['tmp_dir']
            tmp_path_del = os.path.join(tmp_path, '*')
            cmd = f'rm -r {tmp_path_del}'
            print(cmd)
            prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            out, err = prog.communicate()
            if err:
                print(err)

    # create list of sat granules
    if not args.config_file:
        if args.resolution == 'WRR':
            res = 'WRR'
        else:
            res = 'WFR'
    else:
        res = 'WFR'
        if options.has_option('satellite_options', 'resolution'):
            res = options['satellite_options']['resolution']

    # path to output
    if args.output:
        path_out = args.output
    elif args.config_file:
        if options['file_path']['output_dir']:
            path_out = options['file_path']['output_dir']

    if args.verbose:
        print(f'Path to output: {path_out}')
    if not os.path.isdir(path_out):
        os.mkdir(path_out)

    # in situ sites
    in_situ_sites = {}
    if args.sitename:
        station_name = args.sitename
        in_situ_lat, in_situ_lon = cfs.get_lat_lon_ins(station_name)  # in situ location based on the station name
        in_situ_sites[station_name] = {
            'latitude': in_situ_lat,
            'longitude': in_situ_lon,
            'path_out': os.path.join(path_out, station_name)
        }
    elif args.config_file:
        if not options.has_option('Time_and_sites_selection', 'sites'):
            print(
                f'Error: Section Time_and_sites_selection, Option sites is not available in the config file: {args.config_file}')
            return

        if not options.has_option('Time_and_sites_selection', 'sites_file') or \
                (options.has_option('Time_and_sites_selection', 'sites_file') and not
                options['Time_and_sites_selection']['sites_file']):
            site_list = options['Time_and_sites_selection']['sites'].split(',')
            site_list = [s.strip() for s in site_list]
            in_situ_sites = cfs.get_sites_from_list(site_list, path_out)
        else:
            site_list = None
            region_list = None
            if options['Time_and_sites_selection']['sites']:
                site_list = options['Time_and_sites_selection']['sites'].split(',')
                site_list = [s.strip() for s in site_list]
            if not options['Time_and_sites_selection']['sites'] and \
                    (options.has_option('Time_and_sites_selection', 'sites_region') and
                     options['Time_and_sites_selection']['sites_region']):
                region_list = options['Time_and_sites_selection']['sites_region'].split(',')
                region_list = [r.strip() for r in region_list]
            in_situ_sites = cfs.get_sites_from_file(options['Time_and_sites_selection']['sites_file'], site_list,
                                                    region_list, path_out)

    if args.verbose:
        for site in in_situ_sites:
            lath = in_situ_sites[site]['latitude']
            lonh = in_situ_sites[site]['longitude']
            print(f'station_name: {site} with lat: {lath}, lon: {lonh}')

    if args.config_file:
        datetime_start = datetime.strptime('2000-01-01', '%Y-%m-%d')
        datetime_end = datetime.today()
        if options.has_option('Time_and_sites_selection', 'time_start'):
            try:
                datetime_start = datetime.strptime(options['Time_and_sites_selection']['time_start'], '%Y-%m-%d')
            except:
                print(f'WARNING: time_start format is not valid. Usind dafult value: 2000-01-01')

        if options.has_option('Time_and_sites_selection', 'time_end'):
            try:
                datetime_end = datetime.strptime(options['Time_and_sites_selection']['time_stop'],
                                                 '%Y-%m-%d') + timedelta(seconds=59, minutes=59, hours=23)
            except:
                print(f'WARNING: time_end format is not valid. Usind dafult value: today')
    else:
        if args.startdate:
            datetime_start = datetime.strptime(args.startdate, '%Y-%m-%d')
        else:
            datetime_start = datetime.strptime('2000-01-01', '%Y-%m-%d')
        if args.enddate:
            datetime_end = datetime.strptime(args.enddate, '%Y-%m-%d') + timedelta(seconds=59, minutes=59, hours=23)
        else:
            datetime_end = datetime.today()

    wce = f'"*OL_2_{res}*SEN3*"'  # wild card expression

    org = None
    if options.has_option('file_path', 'sat_source_dir_organization') and \
            options['file_path']['sat_source_dir_organization']:
        org = options['file_path']['sat_source_dir_organization']

    path_to_satellite_list = create_list_products(satellite_path_source, path_out, wce, res, 'satellite',
                                                  datetime_start, datetime_end, org)

    if args.verbose:
        print(f'Satellite List: {path_to_satellite_list}')

    # if os.path.exists(f'{path_out}/OL_2_{res}_list.txt'):
    #     os.remove(f'{path_out}/OL_2_{res}_list.txt')

    # create extract and save it in internal folder
    size_box = 25
    if args.config_file:
        if options.has_option('satellite_options', 'extract_size'):
            size_box = int(options['satellite_options']['extract_size'])

    if args.verbose:
        print(f'Start date: {datetime_start}')
        print(f'End date: {datetime_end}')

    make_brdf = False
    if args.config_file:
        if options.has_option('satellite_options', 'brdf'):
            make_brdf = True if options['satellite_options']['brdf'] == 'T' else False

    day_ref = -1
    with open(path_to_satellite_list, 'r') as file:
        for cnt, line in enumerate(file):
            if args.verbose:
                print('-----------------')
            path_to_sat_source = line[:-1]
            if path_to_sat_source.endswith('.zip'):
                if not tmp_path is None and os.path.exists(tmp_path):
                    cmd = f'unzip -o {path_to_sat_source} -d {tmp_path}'
                    prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
                    out, err = prog.communicate()
                    if err:
                        print(err)
                    namesat = path_to_sat_source.split('/')[-1]
                    path_to_sat_source = os.path.join(tmp_path, namesat[0:namesat.find('.zip')])
                    if not os.path.exists(path_to_sat_source):
                        continue
                else:
                    continue
            # extract date time info
            sensor_str = path_to_sat_source.split('/')[-1].split('_')[0]
            res_str = path_to_sat_source.split('/')[-1].split('_')[3]
            datetime_str = path_to_sat_source.split('/')[-1].split('_')[7]
            if args.verbose:
                print(f'{datetime_str} {sensor_str} {res_str}')
            date_format = '%Y%m%dT%H%M%S'
            satellite_datetime = datetime.strptime(datetime_str, date_format)
            day_of_year = int(satellite_datetime.strftime('%j'))

            if datetime_start <= satellite_datetime <= datetime_end:
                launch_create_extract(in_situ_sites, size_box, path_to_sat_source, res_str, make_brdf)
                if day_of_year != day_ref and not tmp_path is None:
                    if day_ref != -1:
                        print('Deleting temporary files...')
                        tmp_path_del = os.path.join(tmp_path, '*')
                        cmd = f'rm -r {tmp_path_del}'
                        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
                        out, err = prog.communicate()
                        if err:
                            print(err)
                    day_ref = day_of_year


            else:
                if args.verbose:
                    print('Warning: Out of time frame.')


# %%
if __name__ == '__main__':
    main()
