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
import zipfile

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
import pandas as pd
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
from COMMON.check_geo import CHECK_GEO
from SAT_EXTRACT.sat_extract import SatExtract

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
parser.add_argument('-nl', "--nolist",
                    help="Do not create initial satellite lists, checking day by day allowing download",
                    action="store_true")
parser.add_argument('-adownload', "--allow_download", help="Allow download", action="store_true")

args = parser.parse_args()


class SatExtractOLCI(SatExtract):

    def __init__(self, ofname, variable_list):
        if ofname is not None:
            SatExtract.__init__(self, ofname)
        if variable_list is not None:
            self.variable_list = variable_list
        else:
            self.variable_list = {}
        self.sensor = 'olci'
        self.n_bands = 16
        self.wl_list = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.75, 681.25, 708.75, 753.75, 778.75, 865, 885,
                        1020.5]

    def set_variable_list(self, path_source, extra_bands):
        self.variable_list = {}
        reflectance_bands = {
            'Oa01': 0,
            'Oa02': 1,
            'Oa03': 2,
            'Oa04': 3,
            'Oa05': 4,
            'Oa06': 5,
            'Oa07': 6,
            'Oa08': 7,
            'Oa09': 8,
            'Oa10': 9,
            'Oa11': 10,
            'Oa12': 11,
            'Oa16': 12,
            'Oa17': 13,
            'Oa18': 14,
            'Oa21': 15,
        }
        other_bands = {
            'T865': {
                'band': 'satellite_AOT_0865P50',
                'description': 'Satellite Aerosol optical thickness'
            },
            'WQSF': {
                'band': 'satellite_WQSF',
                'description': 'Satellite Level 2 WATER Product, Classification, Quality and Science Flags Data Set'
            }

        }
        if extra_bands is not None:
            for key in extra_bands:
                other_bands[key] = {
                    'band': f'satellite_{key}',
                    'description': ''
                }
        extra_bands = other_bands


        for name in os.listdir(path_source):
            if not name.endswith('nc'):
                continue
            if name.startswith('tie'):
                continue
            fnc = os.path.join(path_source, name)
            dataset = Dataset(fnc)
            for name_var, variable in dataset.variables.items():
                dims_names = variable.get_dims()
                if len(dims_names) != 2:
                    continue
                if dims_names[0].name == 'rows' and dims_names[1].name == 'columns':
                    apply = 0
                    index_rrs = -1
                    extract_band = ''
                    desc = ''
                    if name_var.endswith('_reflectance'):
                        apply = 2
                        band = name_var.split('_')[0]
                        index_rrs = reflectance_bands[band]
                        extract_band = 'satellite_Rrs'  # also defined using apply==2
                    if name_var in extra_bands.keys():
                        apply = 1
                        extract_band = extra_bands[name_var]['band']
                        desc = extra_bands[name_var]['description']
                    self.variable_list[name_var] = {
                        'path': name,
                        'apply': apply,
                        'index_rrs': index_rrs,
                        'extract_band': extract_band,
                        'description': desc
                    }
            dataset.close()

    def set_variable_data(self, path_source, window):
        for name_var in self.variable_list:

            if self.variable_list[name_var]['apply'] == 0:
                continue
            file_path = os.path.join(path_source, self.variable_list[name_var]['path'])
            if not os.path.exists(file_path):
                print(f'[WARNING] Path to variable: {file_path} not found. Data will not be available')
                continue

            if self.variable_list[name_var]['apply'] == 2:  ##satellite rrs band, it was already created
                index_rrs = self.variable_list[name_var]['index_rrs']
                if args.verbose:
                    print(f'[INFO] Creating rrs band: {index_rrs}')
                if 'satellite_Rrs' not in self.EXTRACT.variables:
                    return False
                extract_variable = self.EXTRACT.variables['satellite_Rrs']
                dataset = Dataset(file_path, 'r')
                rrs_array = dataset.variables[name_var][:]
                dataset.close()
                extract_variable[0, index_rrs, :, :] = ma.array(
                    rrs_array[window[0]:window[1], window[2]:window[3]]) / np.pi

            if self.variable_list[name_var]['apply'] == 1:  ##other bands, not previously created
                extract_band = self.variable_list[name_var]['extract_band']
                description = self.variable_list[name_var]['description']
                dataset = Dataset(file_path, 'r')
                var_atts = dataset.variables[name_var].ncattrs()
                if len(description)==0 and 'long_name' in var_atts:
                    description = dataset.variables[name_var].long_name
                if 'flag_masks' in var_atts and 'flag_meanings' in var_atts:
                    if args.verbose:
                        print(f'[INFO] Creating flag band: {extract_band}')
                    flag_band = dataset.variables[name_var]
                    self.create_flag_variable(extract_band, flag_band, description, flag_band.flag_masks,
                                              flag_band.flag_meanings, window)
                else:
                    if args.verbose:
                        print(f'[INFO] Creating band: {extract_band}')
                    var_array = dataset.variables[name_var]
                    extract_var = self.create_2D_variable_general(extract_band, var_array, window)
                    extract_var.description = description




                dataset.close()

        return True

    def set_geometry_data(self, path_source, size_box, window):
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



        start_idx_y = window[0]
        start_idx_x = window[1]
        for yy in range(size_box):
            for xx in range(size_box):
                yPos = start_idx_y + yy
                xPos = start_idx_x + xx
                self.EXTRACT.variables['satellite_SZA'][0, yy, xx] = get_val_from_tie_point_grid(yPos, xPos, ysubsampling, xsubsampling, SZA)
                self.EXTRACT.variables['satellite_SAA'][0, yy, xx] = get_val_from_tie_point_grid(yPos, xPos, ysubsampling, xsubsampling, SAA)
                self.EXTRACT.variables['satellite_OZA'][0, yy, xx] = get_val_from_tie_point_grid(yPos, xPos, ysubsampling, xsubsampling, OZA)
                self.EXTRACT.variables['satellite_OAA'][0, yy, xx] = get_val_from_tie_point_grid(yPos, xPos, ysubsampling, xsubsampling, OAA)

    def get_version(self, path_source):
        # extract IFP-OL-2 version
        proc_version_str = 'IPF-OL-2'
        with open(os.path.join(path_source, 'xfdumanifest.xml'), 'r', encoding="utf-8") as read_obj:
            check_version = False
            for line in read_obj:
                if 'IPF-OL-2' in line and check_version == False:
                    IPF_OL_2_version = line.split('"')[3]
                    proc_version_str = f'IPF-OL-2 version {IPF_OL_2_version}'
                    if args.verbose:
                        print(f'[INFO] Version: {proc_version_str}')
                    check_version = True
                    pass
        return proc_version_str

    def get_global_atrib(self, path_source, options):
        proc_version_str = self.get_version(path_source)
        filename = path_source.split('/')[-1].replace('.', '_')
        satellite = filename[0:2]
        platform = filename[2]

        at = {'sensor': self.sensor, 'satellite': satellite, 'platform': platform, 'res': '',
              'aco_processor': 'STANDARD', 'proc_version': proc_version_str, 'station_name': '', 'in_situ_lat': -999,
              'in_situ_lon': -999}

        if options is not None:
            at['station_name'] = options['station_name']
            at['in_situ_lat'] = options['in_situ_lat']
            at['in_situ_lon'] = options['in_situ_lon']
            at['res'] = options['resolution']

        return at

    def get_sat_time(self, path_source):
        filepah = os.path.join(path_source, 'Oa01_reflectance.nc')
        nc_sat = Dataset(filepah, 'r')
        satellite_start_time = datetime.strptime(nc_sat.start_time, "%Y-%m-%dT%H:%M:%S.%fZ")
        nc_sat.close()
        return satellite_start_time


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


def check_wce(name, wce):
    check_w = True
    if wce is None:
        return check_w
    wce = wce.replace("\"", "")
    wces = wce.strip().split('*')
    for s in wces:
        if not s:
            continue
        if name.find(s) < 0:
            check_w = False
    return check_w


def get_dates_and_platform_from_file_name(name):
    from datetime import datetime as dt
    platform = 'S3'
    start_date = dt.now()
    end_date = dt.now()
    try:
        platform = name.split('_')[0]
        start_date = dt.strptime(name.split('_')[7], '%Y%m%dT%H%M%S')
        end_date = dt.strptime(name.split('_')[8], '%Y%m%dT%H%M%S')
    except:
        pass
    return platform, start_date, end_date


def check_products_to_download(info_edac, products, info_path):
    products_to_download = []
    for p in products:
        pname = str(p)
        check_info_path = False

        for name in info_path:
            if info_path[name]['platform'] != info_edac[pname]['platform']:
                continue
            if info_path[name]['start_date'] >= info_edac[pname]['start_date'] and info_path[name]['end_date'] <= \
                    info_edac[pname]['end_date']:
                check_info_path = True
        if not check_info_path:
            products_to_download.append(p)

    return products_to_download


def create_extracts_day_by_day(date_list, path_source, site_name, insitu_lat, insitu_lon, wce, unzip_path, path_out,
                               size_box, make_brdf):
    ncreated = 0

    for idx in range(len(date_list)):
        fproducts, iszipped = get_olci_products_day_download(date_list[idx], path_source, insitu_lat, insitu_lon, wce,
                                                             unzip_path)
        nproducts = len(fproducts)
        if nproducts == 0:
            if args.verbose:
                print(f'[WARNING] No products found for {date_list[idx]}. Skipping day...')
            continue

        for id in range(len(fproducts)):
            path_product = fproducts[id]
            if args.verbose:
                print('---------------------------------------------------')
                print(f'DATE: {date_list[idx]}')
                print(f'GRANULE: {path_product}')
            res_str = path_product.split('/')[-1].split('_')[3]

            ofname = create_extract(size_box, site_name, path_product, path_out, insitu_lat, insitu_lon, res_str,
                                    make_brdf, None)
            if ofname is not None:
                if args.verbose:
                    print(f'Sat extract {ofname} was created')
                ncreated = ncreated + 1

    print('------------------------------')
    print(f'COMPLETED. {ncreated} sat extract files were created')


def get_olci_products_day_download(date, path_source, insitu_lat, insitu_lon, wce, unzip_path):
    year = date.strftime('%Y')
    jday = date.strftime('%j')
    path_year = os.path.join(path_source, year)
    path_day = os.path.join(path_year, jday)
    info_path = {}
    ##final variables
    fproducts = []
    iszipped = []
    if os.path.exists(path_day):
        for name in os.listdir(path_day):
            if not check_wce(name, wce):
                continue
            prod_path = os.path.join(path_day, name)
            if name.endswith('.SEN3') or zipfile.is_zipfile(prod_path):
                if args.verbose:
                    print(f'[INFO] Adding already available granule {name} to the path list')
                platform, start_date, end_date = get_dates_and_platform_from_file_name(name)
                info_path[prod_path] = {
                    'platform': platform,
                    'start_date': start_date,
                    'end_date': end_date,
                }

    edac = check_donwload()
    if edac is not None:
        info_edac, products = check_list_products_eumetsat(edac, date, insitu_lat, insitu_lon)
        if len(info_edac) == 0:
            products_to_download = []  ##no products to be downloaded
        else:
            if len(info_path) == 0:  ##all the products in info_edac must me donwloaded
                products_to_download = products
            else:
                products_to_download = check_products_to_download(info_edac, products, info_path)

        if len(products_to_download) > 0:
            if args.verbose:
                print(f'[INFO] {len(products_to_download)} product(s) not found. Starting download...')
            if not os.path.exists(path_day):
                if not os.path.exists(path_year):
                    os.mkdir(path_year)
                os.mkdir(path_day)
            path_products, path_unavailable = launch_download(edac, date, path_day, insitu_lat, insitu_lon,
                                                              products_to_download)

            if path_products is None and info_path is None:
                print(f'[WARNING] No granules were found or available for downloading for {date}.')
                return fproducts, iszipped
            if len(path_products) > 0:
                for path_product in path_products:
                    name = path_product.split('/')[-1]
                    if args.verbose:
                        print(f'[INFO] Adding downloaded granule {name} to the path list')
                    platform, start_date, end_date = get_dates_and_platform_from_file_name(name)
                    info_path[path_product] = {
                        'platform': platform,
                        'start_date': start_date,
                        'end_date': end_date,
                    }

    if len(info_path) == 0:
        if edac is None:
            postmessage = ' Download is not enabled.'
        else:
            postmessage = ' Granules could not be downloaded.'
        print(f'[WARNING] No granules were found for {date}.{postmessage}')
        return fproducts, iszipped

    if args.verbose:
        print('[INFO] Checking granules...')

    cgeo = CHECK_GEO()
    for prod_path in info_path:
        name = prod_path.split('/')[-1]
        do_zip_here = False
        contain_flag = 0
        if name.endswith('.SEN3'):
            path_prod_u = prod_path
            cgeo.start_polygon_from_prod_manifest_file(path_prod_u)
            contain_flag = cgeo.check_point_lat_lon(insitu_lat, insitu_lon)
        elif zipfile.is_zipfile(prod_path):
            cgeo.start_polygon_image_from_zip_manifest_file(prod_path)
            contain_flag = cgeo.check_point_lat_lon(insitu_lat, insitu_lon)
            if contain_flag == 1:
                do_zip_here = True
                path_prod_t = prod_path.split('/')[-1][0:-4]
                path_prod_u = os.path.join(unzip_path, path_prod_t)
                if os.path.exists(path_prod_u):
                    do_zip_here = False
        if args.verbose:
            print(f'[INFO] Checking product {name} contain flag: {contain_flag}')
        if do_zip_here:
            with zipfile.ZipFile(prod_path, 'r') as zprod:
                if args.verbose:
                    print(f'[INFO] Unziping {name} to {unzip_path}')
                zprod.extractall(path=unzip_path)
        if contain_flag == 1:
            fproducts.append(path_prod_u)
            iszipped.append(do_zip_here)

    return fproducts, iszipped


def check_list_products_eumetsat(edac, date, insitu_lat, insitu_lon):
    date_str = date.strftime('%Y-%m-%d')
    products, product_names, collection_id = edac.search_olci_by_point(date_str, 'FR', 'L2', insitu_lat, insitu_lon, -1,
                                                                       -1, 'NT')

    info = {}
    if len(product_names) == 0:
        return info, None

    for name in product_names:
        platform, start_date, end_date = get_dates_and_platform_from_file_name(name)
        info[name] = {
            'platform': platform,
            'start_date': start_date,
            'end_date': end_date,
            'available': False
        }
    return info, products


def launch_download(edac, date, path_day, insitu_lat, insitu_lon, products):
    if args.verbose:
        print(f'[INFO] Launching download for day: {date}')

    product_names = None
    if products is None:
        date_str = date.strftime('%Y-%m-%d')
        products, product_names, collection_id = edac.search_olci_by_point(date_str, 'FR', 'L2', insitu_lat, insitu_lon,
                                                                           -1, -1)

    if products is None:
        return None

    if products is not None and product_names is None:
        product_names = [str(p) for p in products]
    edac.download_product_from_product_list(products, path_day, False)

    path_products = []
    path_unavailable = []
    for name in product_names:
        paths_here = [os.path.join(path_day, name), os.path.join(path_day, f'{name}.zip'),
                      os.path.join(path_day, f'{name}.tar')]
        available = False
        for path_here in paths_here:
            if os.path.exists(path_here):
                path_products.append(path_here)
                available = True
        if not available:
            path_unavailable.append(path_here)

    return path_products, path_unavailable


def check_donwload():
    import sat_extract
    code_home = os.path.dirname(os.path.dirname(os.path.dirname(sat_extract.__file__)))
    code_download = os.path.join(code_home, 'cnrdownload')
    if os.path.exists(code_download):
        sys.path.append(code_download)
        try:
            from eumdac_lois import EUMDAC_LOIS
            edac = EUMDAC_LOIS(True, None)
        except:
            print(f'[WARNING] Error loading package eumdac_lois. Download is not enabled')
            return None

        return edac
    else:
        print(f'[WARNING] Package {code_download} is not available. Download is not enabled')
        return None


# user defined functions
def create_list_products(path_source, path_out, wce, res_str, date_list, org):
    if args.nolist:
        return None

    path_to_list = f'{path_out}/file_satellite_{res_str}_list.txt'

    if os.path.exists(path_to_list):
        print('Deleting previous file list...')
        cmd = f'rm {path_to_list}'
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

    if org == 'yyyy/jjj':
        for dt in date_list:
            print('Date: ', dt)
            year = dt.strftime('%Y')
            jday = dt.strftime('%j')
            path_day = os.path.join(path_source, year, jday)
            cmd = f'find {path_day} -name {wce}|sort|uniq>> {path_to_list}'
            prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            out, err = prog.communicate()
            if err:
                print(err)

    if args.verbose:
        print(f'Satellite List: {path_to_list}')

    return path_to_list


def get_date_list_from_start_end_date(dt_start, dt_end):
    date_list = []
    dt = dt_start
    while dt <= dt_end:
        date_list.append(dt)
        dt = dt + timedelta(hours=24)
    return date_list


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
            datehere = datetime.strptime(dateherestr, '%Y-%m-%d')
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


def get_params_time(options):
    date_list = None
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

        if options.has_option('Time_and_sites_selection', 'time_list_file'):
            time_list_file = options['Time_and_sites_selection']['time_list_file']
            if time_list_file:
                date_list, datetime_start, datetime_end = get_date_list_from_file(time_list_file, datetime_start,
                                                                                  datetime_end)

    else:
        if args.startdate:
            datetime_start = datetime.strptime(args.startdate, '%Y-%m-%d')
        else:
            datetime_start = datetime.strptime('2000-01-01', '%Y-%m-%d')
        if args.enddate:
            datetime_end = datetime.strptime(args.enddate, '%Y-%m-%d') + timedelta(seconds=59, minutes=59, hours=23)
        else:
            datetime_end = datetime.today()

    if date_list is None:
        date_list = get_date_list_from_start_end_date(datetime_start, datetime_end)
    ndates = len(date_list)
    if args.verbose:
        print(f'[INFO] Start date: {datetime_start}')
        print(f'[INFO] End date: {datetime_end}')
        print(f'[INFO] # of dates: {ndates}')

    return datetime_start, datetime_end, date_list


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


def launch_create_extract_syke(filepath, skie_file, options):
    ncreated = 0

    path_output = get_output_path(options)
    if path_output is None:
        print(f'ERROR: {path_output} is not valid')
        return ncreated

    # Start dataset
    # if args.verbose:
    #     if args.verbose:
    #         print('[INFO] Starting dataset...')
    # nc_sat = Dataset(filepath, 'r')

    # Retriving lat and long arrays
    if args.verbose:
        print('[INFO] Retrieving lat/long data...')
    lat, lon = get_lat_long_arrays(filepath)

    sat_time = get_sat_time(filepath)
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

        b = create_extract_syke(ofname, pdu, options, filepath, global_at, lat, lon, r, c, skie_file,
                                extracts[extract]['irows'])

        if b:
            ncreated = ncreated + 1

    return ncreated


def create_extract_syke(ofname, pdu, options, filepath, global_at, lat, long, r, c, skie_file, irows):
    size_box = get_box_size(options)
    start_idx_x = (c - int(size_box / 2))
    stop_idx_x = (c + int(size_box / 2) + 1)
    start_idx_y = (r - int(size_box / 2))
    stop_idx_y = (r + int(size_box / 2) + 1)
    window = [start_idx_y, stop_idx_y, start_idx_x, stop_idx_x]

    reflectance_bands, n_bands = get_reflectance_bands_info(filepath)
    if n_bands == 0:
        print('[ERROR] reflectance bands are not defined')
        return False
    flag_band_name = 'wqsf.nc'

    newEXTRACT = SatExtract(ofname)
    if not newEXTRACT.FILE_CREATED:
        print(f'[ERROR] File {ofname} could not be created')
        return False

    if args.verbose:
        print(f'[INFO]  Creating file: {ofname}')
    newEXTRACT.set_global_attributes(global_at)
    if skie_file is not None:
        newEXTRACT.create_dimensions_incluidinginsitu(size_box, n_bands, skie_file.get_n_bands(), 30)
    else:
        newEXTRACT.create_dimensions(size_box, n_bands)

    newEXTRACT.create_lat_long_variables(lat, long, window)

    # Sat time start:  ,+9-2021-12-24T18:23:00.471Z
    sat_start_time = get_sat_time(filepath)
    newEXTRACT.create_satellite_time_variable(sat_start_time)

    # pdu variable
    newEXTRACT.create_pdu_variable(pdu, global_at['sensor'])

    # Rrs and wavelenghts
    satellite_Rrs = newEXTRACT.create_rrs_variable(global_at['sensor'])
    wavelenghts = []
    index = 0
    for rband in reflectance_bands:
        path_band = reflectance_bands[rband]['file_path']
        nc_sat = Dataset(path_band)
        # nc_sat.variables['Oa01_reflectance'][:]
        bandarray = ma.array(nc_sat.variables[rband][:, :])
        satellite_Rrs[0, index, :, :] = bandarray[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x] / np.pi
        nc_sat.close()
        wl = reflectance_bands[rband]['wavelenght']
        wavelenghts.append(wl)
        index = index + 1
    newEXTRACT.create_satellite_bands_variable(wavelenghts)

    # flags
    flag_file = os.path.join(filepath, flag_band_name)
    nc_sat = Dataset(flag_file)
    flag_band = nc_sat.variables['WQSF']
    newEXTRACT.create_flag_variable('satellite_WQSF', flag_band, flag_band.long_name, flag_band.flag_masks,
                                    flag_band.flag_meanings, window)
    nc_sat.close()

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
            sat_time_here = sat_start_time
            time_diff = float(abs((insitu_time_here - sat_time_here).total_seconds()))
            insitu_timedif_var[0, index_var] = time_diff
            latp, lonp = skie_file.get_lat_lon_at_subdb(irow)
            insitu_lat[0, index_var] = latp
            insitu_lon[0, index_var] = lonp
            index_var = index_var + 1

    newEXTRACT.close_file()

    return True


def get_reflectance_bands_info(path_source):
    reflectance_bands = {}
    bands_ints = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 17, 18, 21]
    wl_list = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.75, 681.25, 708.75, 753.75, 778.75, 865, 885, 1020.5]
    nbands = len(bands_ints)
    for index in range(nbands):
        v = bands_ints[index]
        wl = wl_list[index]
        band_name = f'Oa{v:02d}_reflectance'
        band_path = f'{band_name}.nc'
        file_path = os.path.join(path_source, band_path)
        reflectance_bands[band_name] = {
            'wavelenght': wl,
            'file_path': file_path
        }

    return reflectance_bands, nbands


def get_lat_long_arrays(path_source):
    coordinates_filename = 'geo_coordinates.nc'
    filepah = os.path.join(path_source, coordinates_filename)
    nc_sat = Dataset(filepah, 'r')
    lat = nc_sat.variables['latitude'][:, :]
    lon = nc_sat.variables['longitude'][:, :]

    return lat, lon


def check_location_source(path_source, in_situ_lat, in_situ_lon, size_box):
    lat, lon = get_lat_long_arrays(path_source)

    contain_flag = cfs.contain_location(lat, lon, in_situ_lat, in_situ_lon)
    r = -1
    c = -1
    if contain_flag == 1:
        r, c = cfs.find_row_column_from_lat_lon(lat, lon, in_situ_lat, in_situ_lon)

        start_idx_y, stop_idx_y, start_idx_x, stop_idx_x = get_window_limits(r, c, size_box)

        if start_idx_y >= 0 and (stop_idx_y + 1) < lat.shape[0] and start_idx_x >= 0 and (stop_idx_x + 1) < lat.shape[
            1]:
            pass
        else:
            r = -1
            c = -1

    return r, c


def get_window_limits(r, c, size_box):
    start_idx_y = (r - int(size_box / 2))
    stop_idx_y = (r + int(size_box / 2) + 1)
    start_idx_x = (c - int(size_box / 2))
    stop_idx_x = (c + int(size_box / 2) + 1)

    return start_idx_y, stop_idx_y, start_idx_x, stop_idx_x


def get_sat_time(path_source):
    filepah = os.path.join(path_source, 'Oa01_reflectance.nc')
    nc_sat = Dataset(filepah, 'r')
    satellite_start_time = datetime.strptime(nc_sat.start_time, "%Y-%m-%dT%H:%M:%S.%fZ")
    nc_sat.close()
    return satellite_start_time


def get_global_atrib(path_source):
    filename = path_source.split('/')[-1].replace('.', '_')
    satellite = filename[0:2]
    platform = filename[2]
    sensor = 'olci'
    at = {'sensor': sensor, 'satellite': satellite, 'platform': platform, 'res': '',
          'aco_processor': 'STANDARD', 'proc_version': '', 'station_name': '', 'in_situ_lat': -999, 'in_situ_lon': -999}

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


def get_olci_products_day(path_source, unzip_path, org, wce, lathere, lonhere, datehere):
    path_search = path_source
    if wce is not None:
        wce = wce.replace('*', '')
    if org is not None:
        if org == 'YYYYmmdd':
            path_search = os.path.join(path_source, datehere.strftime('%Y%m%d'))
        if org == 'YYYYjjj':
            path_search = os.path.join(path_source, datehere.strftime('%Y'), datehere.strftime('%j'))

    fproducts = []
    iszipped = []
    cgeo = CHECK_GEO()

    for name in os.listdir(path_search):
        if wce is not None:
            if name.find(wce) < 0:
                continue
        prod_path = os.path.join(path_search, name)
        do_zip_here = False
        contain_flag = 0
        if name.endswith('.SEN3'):
            path_prod_u = prod_path
            cgeo.start_polygon_from_prod_manifest_file(path_prod_u)
            contain_flag = cgeo.check_point_lat_lon(lathere, lonhere)
        elif zipfile.is_zipfile(prod_path):
            cgeo.start_polygon_image_from_zip_manifest_file(prod_path)
            contain_flag = cgeo.check_point_lat_lon(lathere, lonhere)
            if contain_flag == 1:
                do_zip_here = True
                path_prod_u = prod_path.split('/')[-1][0:-4]
                path_prod_u = os.path.join(unzip_path, path_prod_u)
        if do_zip_here:
            with zipfile.ZipFile(prod_path, 'r') as zprod:
                if args.verbose:
                    print(f'[INFO] Unziping {name} to {unzip_path}')
                zprod.extractall(path=unzip_path)
        if contain_flag == 1:
            fproducts.append(path_prod_u)
            iszipped.append(do_zip_here)

    return fproducts, iszipped


def launch_create_extract(in_situ_sites, size_box, path_source, res_str, make_brdf):
    for site in in_situ_sites:
        try:
            in_situ_lat = float(in_situ_sites[site]['latitude'])
            in_situ_lon = float(in_situ_sites[site]['longitude'])
            path_output = in_situ_sites[site]['path_out']
            if not os.path.exists(path_output):
                os.mkdir(path_output)
            if args.verbose:
                print(f'Creating extract for site: {site}')
            extract_path = create_extract(size_box, site, path_source, path_output, in_situ_lat, in_situ_lon, res_str,
                                          make_brdf, None)
            if not extract_path is None:
                print(f'file created: {extract_path}')
        except Exception as e:
            if args.verbose:
                print(f'Exception: {e}')
                pass


def create_extractv2(path_product, path_output, options):
    size_box = options['size_box']
    station_name = options['station_name']
    in_situ_lat = options['in_situ_lat']
    in_situ_lon = options['in_situ_lon']
    res_str = options['resolution']
    make_brdf = options['make_brdf']
    insitu_info = options['insitu_info']
    variable_list = options['variable_list']
    if args.verbose:
        print(f'[INFO] Creating extract v.2 for {station_name} from {path_product}')

    if args.verbose:
        print(f'[INFO] Checking in situ  location -> lat: {in_situ_lat}; lon: {in_situ_lon}')
    r, c = check_location_source(path_product, in_situ_lat, in_situ_lon, size_box)
    if r == -1 and c == -1:
        if args.verbose:
            print('f[WARNING] File does NOT contains the in situ location! Skipping...')
        return

    start_idx_y, stop_idx_y, start_idx_x, stop_idx_x = get_window_limits(r, c, size_box)
    window = [start_idx_y, stop_idx_y, start_idx_x, stop_idx_x]
    if args.verbose:
        print(f'[INFO] Getting extraction window: {stop_idx_y}-{stop_idx_y}:{start_idx_x}-{stop_idx_x}')

    filename = path_product.split('/')[-1].replace('.', '_') + '_extract_' + station_name + '.nc'
    ofname = os.path.join(path_output, filename)
    if args.verbose:
        print(f'[INFO] Starting extract file: {ofname}')
    newExtract = SatExtractOLCI(ofname, variable_list)
    newExtract.create_dimensions(size_box, newExtract.n_bands)
    newExtract.set_global_attributes(newExtract.get_global_atrib(path_product, options))
    lat, lon = get_lat_long_arrays(path_product)
    newExtract.create_lat_long_variables(lat, lon, window)
    newExtract.create_satellite_time_variable(newExtract.get_sat_time(path_product))
    newExtract.create_satellite_bands_variable(newExtract.wl_list)
    newExtract.create_rrs_variable(newExtract.sensor)
    newExtract.create_geometry_variables()
    b = newExtract.set_variable_data(path_product, window)
    newExtract.set_geometry_data(path_product,size_box,window)
    newExtract.close_file()

    if not b:
        os.remove(ofname)
        return None

    return ofname

    # ofname = None
    # solci = SatExtractOLCI(ofname)
    # solci.set_version(path_source)

    # ofname = None
    # proc_version_str = check_version(path_source)
    # variable_list = get_variable_list(path_source,extra_bands)
    # for var in variable_list:
    #     print(var,variable_list[var]['apply'],variable_list[var]['path'])


def get_variable_list(path_source, extra_bands):
    variable_list = {}
    reflectance_bands = {
        'Oa01': 0,
        'Oa02': 1,
        'Oa03': 2,
        'Oa04': 3,
        'Oa05': 4,
        'Oa06': 5,
        'Oa07': 6,
        'Oa08': 7,
        'Oa09': 8,
        'Oa10': 9,
        'Oa11': 10,
        'Oa12': 11,
        'Oa16': 12,
        'Oa17': 13,
        'Oa18': 14,
        'Oa21': 15,
    }
    other_bands = {
        'T865': {
            'band': 'satellite_AOT_0865P50',
            'description': 'Satellite Aerosol optical thickness'
        },
        'WQSF': {
            'band': 'satellite_WQSF',
            'description': 'Satellite Level 2 WATER Product, Classification, Quality and Science Flags Data Set'
        }

    }
    if extra_bands is None:
        extra_bands = other_bands

    for name in os.listdir(path_source):
        if not name.endswith('nc'):
            continue
        if name.startswith('tie'):
            continue
        fnc = os.path.join(path_source, name)
        dataset = Dataset(fnc)
        for name_var, variable in dataset.variables.items():
            dims_names = variable.get_dims()
            if len(dims_names) != 2:
                continue
            if dims_names[0].name == 'rows' and dims_names[1].name == 'columns':
                apply = 0
                index_rrs = -1
                if name_var.endswith('_reflectance'):
                    apply = 2
                    band = name_var.split('_')[0]
                    index_rrs = reflectance_bands[band]
                variable_list[name_var] = {
                    'path': fnc,
                    'apply': apply,
                    'index_rrs': index_rrs
                }
        dataset.close()
    return variable_list


def check_version(path_source):
    # extract IFP-OL-2 version
    proc_version_str = 'IPF-OL-2'
    with open(os.path.join(path_source, 'xfdumanifest.xml'), 'r', encoding="utf-8") as read_obj:
        check_version = False
        for line in read_obj:
            if 'IPF-OL-2' in line and check_version == False:
                IPF_OL_2_version = line.split('"')[3]
                proc_version_str = f'IPF-OL-2 version {IPF_OL_2_version}'
                if args.verbose:
                    print(f'[INFO] {proc_version_str}')
                check_version = True
                pass
    return proc_version_str


def create_extract(size_box, station_name, path_source, path_output, in_situ_lat, in_situ_lon, res_str, make_brdf,
                   insitu_info):
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
                    print(f'[INFO] {proc_version_str}')
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
            new_EXTRACT.satellite_aco_processor = 'STANDARD'
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
            SZAO = new_EXTRACT.createVariable('satellite_SZA', 'f4', ('satellite_id', 'rows', 'columns'),
                                              fill_value=-999,
                                              zlib=True, complevel=6)
            SZAO.units = 'degress'
            SZAO.long_name = 'Sun Zenith Angle'

            SAAO = new_EXTRACT.createVariable('satellite_SAA', 'f4', ('satellite_id', 'rows', 'columns'),
                                              fill_value=-999,
                                              zlib=True, complevel=6)
            SAAO.units = 'degress'
            SAAO.long_name = 'Sun Azimuth Angle'

            OZAO = new_EXTRACT.createVariable('satellite_OZA', 'f4', ('satellite_id', 'rows', 'columns'),
                                              fill_value=-999,
                                              zlib=True, complevel=6)
            OZAO.units = 'degress'
            OZAO.long_name = 'Observation Zenith Angle'

            OAAO = new_EXTRACT.createVariable('satellite_OAA', 'f4', ('satellite_id', 'rows', 'columns'),
                                              fill_value=-999,
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
            satellite_time = new_EXTRACT.createVariable('satellite_time', 'f8', ('satellite_id'), fill_value=-999,
                                                        zlib=True, complevel=6)
            satellite_time[0] = float(datetime.strptime(satellite_start_time, "%Y-%m-%dT%H:%M:%S.%fZ").timestamp())
            satellite_time.units = "Seconds since 1970-1-1"

            # PDU- DEPRECATED NOT INCLUDED IN THE NEW VERSION
            # satellite_PDU = new_EXTRACT.createVariable('satellite_PDU', 'S2', ('satellite_id'), zlib=True,
            #                                            complevel=6)  # string
            # satellite_PDU[0] = path_source.split('/')[-1]
            # satellite_PDU.long_name = "OLCI source PDU name"

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
            if args.verbose:
                print('Adding AOT...')
            satellite_AOT_0865p50_box = new_EXTRACT.createVariable('satellite_AOT_0865p50', 'f4',
                                                                   ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                                                   zlib=True, complevel=6)
            satellite_AOT_0865p50_box[0, :, :] = ma.array(AOT_0865p50[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])
            satellite_AOT_0865p50_box.description = 'Satellite Aerosol optical thickness'

            # WQSF: Quality Flags
            if args.verbose:
                print('Adding WQSF flags...')
            satellite_WQSF = new_EXTRACT.createVariable('satellite_WQSF', 'f4', ('satellite_id', 'rows', 'columns'),
                                                        fill_value=-999, zlib=True, complevel=6)
            satellite_WQSF[0, :, :] = [ma.array(WQSF[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])]
            satellite_WQSF.description = 'Satellite Level 2 WATER Product, Classification, Quality and Science Flags Data Set'
            satellite_WQSF.flag_masks = WQSF_flag_masks
            satellite_WQSF.flag_meanings = WQSF_flag_meanings

            # other bands

            # in situ info
            if insitu_info is not None:
                if args.verbose:
                    print('Adding in situ info...')
                insitu_lat_here = insitu_info[0]
                insitu_lon_here = insitu_info[1]
                insitu_time_here = insitu_info[2]
                insitu_time = new_EXTRACT.createVariable('insitu_time', 'f8', ('satellite_id',), zlib=True,
                                                         complevel=6)
                insitu_time.units = "Seconds since 1970-1-1"
                insitu_time.description = 'In situ time in ISO 8601 format (UTC).'
                insitu_time[0] = insitu_time_here
                insitu_lat = new_EXTRACT.createVariable('insitu_latitude', 'f8', ('satellite_id',),
                                                        fill_value=-999,
                                                        zlib=True, complevel=6)
                insitu_lat.short_name = "latitude"
                insitu_lat.units = "degrees"
                insitu_lon = new_EXTRACT.createVariable('insitu_longitude', 'f8', ('satellite_id',),
                                                        fill_value=-999,
                                                        zlib=True, complevel=6)
                insitu_lon.short_name = "longitude"
                insitu_lon.units = "degrees"
                insitu_lat[0] = insitu_lat_here
                insitu_lon[0] = insitu_lon_here

            new_EXTRACT.close()
        else:
            if args.verbose:
                print('Warning: Index out of bound!')
    else:
        if args.verbose:
            print('Warning: File does NOT contains the in situ location!')

    return ofname


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
    if options.has_option('Time_and_sites_selection', 'time_start'):
        time_start = datetime.strptime(options['Time_and_sites_selection']['time_start'], '%Y-%m-%d')
    if options.has_option('Time_and_sites_selection', 'time_stop'):
        time_stop = datetime.strptime(options['Time_and_sites_selection']['time_stop'], '%Y-%m-%d') + timedelta(
            hours=24)

    return path_source, org, wce, time_start, time_stop


def get_path_list_products(options):
    path_out = get_output_path(options)
    if path_out is None:
        return None
    # path_to_list = f'{path_out}/file_{platform}_list.txt'
    path_to_list = f'{path_out}/file_list.txt'
    return path_to_list


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

            STATUS_STR = 'NO PASSED'
            if check_product(path_to_sat_source, time_start, time_stop):
                product_path_list.append(path_to_sat_source)
                STATUS_STR = 'OK'

            if args.verbose:
                print('-----------------')
                pname = path_to_sat_source.split('/')[-1]
                print(f'[INFO ]Check: {pname}. STATUS: {STATUS_STR}')

    return product_path_list


def check_product(filepath, time_start, time_stop):
    if not os.path.isdir(filepath):
        print(f'[WARNING] {filepath} is not a valid dataset. Skipping...')
        return False
    filename = filepath.split('/')[-1]
    if not filename.startswith('S3') or not filename.endswith('SEN3'):
        print(f'[WARNING] {filepath} is not a valid dataset. Skipping...')
        return False
    # try:
    #     nc_sat = Dataset(filepath, 'r')
    # except OSError:
    #     if args.verbose:
    #         print(f'[WARNING] {filepath} is not a valid dataset. Skipping...')
    #     return False

    # checkTime = False
    # time_sat = get_sat_time(filepath)

    # if 'start_time' in nc_sat.ncattrs():
    #     time_sat = dt.strptime(nc_sat.start_time, '%Y-%m-%d %H:%M:%S')
    #     checkTime = True
    # checkTime = True
    check_res = True
    # if not checkTime and check_res:
    #     if args.verbose:
    #         print(f'[WARNING] Attribute time_coverage_start is not available in the dataset. Skipping...')
    #     check_res = False
    # if check_res:
    #     if time_sat < time_start or time_sat > time_stop:
    #         if args.verbose:
    #             print('[WARNING] Product is out of the temporal coverage. Skipping...')
    #         check_res = False

    return check_res


def get_basic_options_from_file_config(options):
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

    # make_download
    make_download = False
    if options.has_option('satellite_options', 'download'):
        if options['satellite_options']['download'].upper() == 'T' or options['satellite_options'][
            'download'].upper() == 'TRUE':
            make_download = True

    #Ohter bands
    extra_bands = None
    if options.has_option('satellite_options','extra_bands'):
        str_val = options['satellite_options']['extra_bands']
        extra_bands = [x.strip() for x in str_val.split(',')]

    basic_options = {
        'satellite_path_source': satellite_path_source,
        'path_out': path_out,
        'tmp_path': tmp_path,
        'size_box': size_box,
        'make_brdf': make_brdf,
        'org': org,
        'resolution': res,
        'make_download': make_download,
        'extra_bands': extra_bands
    }
    return basic_options


def get_basic_options_from_arguments():
    ##parameters with default values
    size_box = 25
    make_brdf = False
    make_download = False
    org = None

    res = 'WFR'
    if args.resolution == 'WRR':
        res = 'WRR'

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

    basic_options = {
        'satellite_path_source': satellite_path_source,
        'path_out': path_out,
        'tmp_path': tmp_path,
        'size_box': size_box,
        'make_brdf': make_brdf,
        'org': org,
        'resolution': res,
        'make_download': make_download
    }
    return basic_options


def create_dir(path):
    try:
        os.mkdir(path)
        return path
    except:
        return None


def get_insitu_sites(options, path_out):
    in_situ_sites = {}
    if args.sitename:
        station_name = args.sitename
        in_situ_lat, in_situ_lon = cfs.get_lat_lon_ins(station_name)  # in situ location based on the station name
        path_out_site = path_out
        if path_out.find(station_name) == -1:
            path_out_site = os.path.join(path_out, station_name)
        in_situ_sites[station_name] = {
            'latitude': in_situ_lat,
            'longitude': in_situ_lon,
            'path_out': path_out_site
        }
    elif args.config_file:
        if not options.has_option('Time_and_sites_selection', 'sites'):
            print(
                f'Error: Section Time_and_sites_selection, Option sites is not available in the config file: {args.config_file}')
            return None
        site_list_option = options['Time_and_sites_selection']['sites'].strip()
        if not site_list_option:
            print(
                f'Error: Section Time_and_sites_selection, Option sites is not available in the config file: {args.config_file}')
            return None
        site_list = site_list_option.split(',')
        site_list = [s.strip() for s in site_list]
        if not options.has_option('Time_and_sites_selection', 'sites_file'):
            in_situ_sites = cfs.get_sites_from_list(site_list, path_out)
        else:
            site_file_option = options['Time_and_sites_selection']['sites_file'].strip()
            if not site_list_option:
                in_situ_sites = cfs.get_sites_from_list(site_list, path_out)
            else:
                region_list = None
                if options.has_option('Time_and_sites_selection', 'sites_region'):
                    region_list_option = options['Time_and_sites_selection']['sites_region'].strip()
                    if region_list_option:
                        region_list = region_list_option.split(',')
                        region_list = [r.strip() for r in region_list]
                in_situ_sites = cfs.get_sites_from_file(site_file_option, site_list, region_list, path_out)

    if args.verbose:
        for site in in_situ_sites:
            lath = in_situ_sites[site]['latitude']
            lonh = in_situ_sites[site]['longitude']
            print(f'[INFO] Station: station_name: {site}  Lat: {lath} Lon: {lonh}')
    return in_situ_sites


#############################
# %%
def main():
    print('Creating satellite extracts.')

    if args.debug:
        print('Entering Debugging Mode:')

    # Getting basic_options from config_file or arguments
    basic_options = None
    if args.config_file:
        if os.path.isfile(args.config_file):
            options = config_reader(args.config_file)
            basic_options = get_basic_options_from_file_config(options)
    if basic_options is None:
        basic_options = get_basic_options_from_arguments()
    if basic_options is None:
        return
    if args.verbose:
        for option in basic_options:
            print(f'[INFO] {option}:{basic_options[option]}')

    size_box = basic_options['size_box']
    res = basic_options['resolution']
    path_out = basic_options['path_out']
    make_brdf = basic_options['make_brdf']
    satellite_path_source = basic_options['satellite_path_source']
    tmp_path = basic_options['tmp_path']
    org = basic_options['org']
    wce = f'"*OL_2_{res}*SEN3*"'  # wild card expression

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


    ##MOVING TARGET CSV SELECTION
    if options.has_section('CSV_SELECTION') and options.has_option('CSV_SELECTION', 'path_csv'):
        path_csv = options['CSV_SELECTION']['path_csv']
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
        tmp_path = path_source
        if args.config_file:
            if options.has_option('file_path', 'tmp_dir') and options['file_path']['tmp_dir']:
                tmp_path = options['file_path']['tmp_dir']
        ncreated = 0
        for idx, row in df.iterrows():
            try:
                datestr = row[col_date].strip()
                datehere = datetime.strptime(datestr, format_date)
                lathere = float(row[col_lat])
                lonhere = float(row[col_lon])
            except:
                print(f'[WARNING] Row {idx} is not valid. Date, latitute and/or longitude could not be parsed')
                continue
            fproducts, iszipped = get_olci_products_day(path_source, tmp_path, org, wce, lathere, lonhere, datehere)
            nproducts = len(fproducts)
            if nproducts == 0:
                if args.verbose:
                    print(f'[WARNING] No products found for {datehere}')
                if args.allow_download:  ##make download in needed
                    date_list = [datehere]
                    site = f'site_{idx}'
                    create_extracts_day_by_day(date_list, satellite_path_source, site, lathere, lonhere, wce,
                                               tmp_path,
                                               path_out, size_box, make_brdf)
                    continue
                continue
            insitu_time = datehere.timestamp()
            insitu_info = [lathere, lonhere, insitu_time]
            variable_list = None

            for id in range(len(fproducts)):
                path_product = fproducts[id]
                if args.verbose:
                    print('---------------------------------------------------')
                    print(f'[INFO] DATE: {datehere}')
                    print(f'[INFO] GRANULE: {path_product}')
                res_str = path_product.split('/')[-1].split('_')[3]
                basic_options['resolution'] = res_str
                ids = f'{idx}_{id}'
                basic_options['station_name'] = ids
                basic_options['in_situ_lat'] = lathere
                basic_options['in_situ_lon'] = lonhere
                basic_options['insitu_info'] = insitu_info
                extra_bands = basic_options['extra_bands']
                if variable_list is None:
                    solci = SatExtractOLCI(None, None)
                    solci.set_variable_list(path_product, extra_bands)
                    variable_list = solci.variable_list
                basic_options['variable_list'] = variable_list
                # variable_list = options['variable_list']

                ofname = create_extractv2(path_product, path_out, basic_options)
                if ofname is not None:
                    if args.verbose:
                        print(f'Sat extract {ofname} was created')
                    ncreated = ncreated + 1
        print('------------------------------')
        print(f'COMPLETED. {ncreated} sat extract files were created')
        return

    # in situ sites
    in_situ_sites = get_insitu_sites(options, path_out)

    # time options
    datetime_start, datetime_end, date_list = get_params_time(options)

    # make_download if needed
    if args.nolist:
        for site in in_situ_sites:
            insitu_lat = in_situ_sites[site]['latitude']
            insitu_lon = in_situ_sites[site]['longitude']
            create_extracts_day_by_day(date_list, satellite_path_source, site, insitu_lat, insitu_lon, wce, tmp_path,
                                       path_out, size_box, make_brdf)
        return

    # satellite list
    path_to_satellite_list = create_list_products(satellite_path_source, path_out, wce, res, date_list, org)

    if args.verbose:
        print(f'Start date: {datetime_start}')
        print(f'End date: {datetime_end}')

    day_ref = -1
    with open(path_to_satellite_list, 'r') as file:
        for cnt, line in enumerate(file):
            if args.verbose:
                print('-----------------')
            path_to_sat_source = line[:-1]
            iszipped = False
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
                    iszipped = True
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
                if day_of_year != day_ref and not tmp_path is None and iszipped:
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
