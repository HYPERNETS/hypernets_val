#!/usr/bin/env python
# coding: utf-8

"""
Create extract from MSI data as a NetCDF4 file
@author: clemence.goyens

Based on EUMETSAT MDB_Builder module (https://ocdb.readthedocs.io/en/latest/ocdb-MDB-user-manual.html)

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


import os
import sys
import argparse
import configparser
from netCDF4 import Dataset
import h5py
import numpy as np
import numpy.ma as ma
from datetime import datetime
from os import listdir
from os.path import isfile, join

code_home = os.path.abspath('../')
print(code_home)
sys.path.append(code_home)
sys.path.append("/home/cgoyens/BackUpThinkpadClem/Projects/HYPERNETS/SatValidation/hypernets_val/")

import COMMON.common_functions as cfs

parser = argparse.ArgumentParser(description="Create Match-up DataBase files (MDB) files.")
parser.add_argument("-d", "--debug", help="Debugging mode.", action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-p', "--product_file", help="Image file.")

args = parser.parse_args()
#

# Names,LAT,LON, type,coordinator
# Names,LAT,LON, type,coordinator
def get_lat_lon_ins(station_name, S2_timestart=None):
    if (station_name == 'ZEEBRUGGE' or station_name == 'M1BE' or station_name == 'MOW1'): # Water RBINS
        Latitude=51.362
        Longitude=3.12
    elif station_name == 'Thornton': # Water RBINS
        Latitude=51.5325
        Longitude=2.95528
    elif (station_name == 'Blankaart' or station_name == 'BSBE' or station_name == 'BNBE'): # Water RBINS
        Latitude=50.988277
        Longitude=2.830315
    elif station_name == 'MESHURO':   # Water LOV
        Latitude=43.32
        Longitude=4.86666666666667
    elif station_name == 'MAFR':   # Water LOV
        if S2_timestart<datetime.strptime("2022/09/27", '%Y/%m/%d'):
            Latitude=45.543777
            Longitude=-1.042159
        else:
            Latitude = 45.54907109791979
            Longitude = -1.0418157401551864
    elif (station_name == 'BERRE' or station_name == 'BEFR'):    # Water LOV
        Latitude=43.4423106
        Longitude=5.0971775
    elif station_name == 'LPAR': # Water IAFE
        Latitude=-34.818
        Longitude=-57.8959
    elif station_name == 'CHASCOMUS':  # Water IAFE
        Latitude=-35.582747
        Longitude=-58.018314
    elif station_name == 'AAOT':   # Water CNR
        Latitude=45.314242
        Longitude=12.508312
    elif (station_name == 'GAIT' or station_name == 'LAKE_GARDA'): # Water CNR
        Latitude=45.57694
        Longitude=10.57944
    elif (station_name == 'Venise' or station_name == 'Venise_PANTHYR' or station_name == 'VEIT'): # Adriatic Sea
        Latitude=45.313900
        Longitude=12.508300
    elif station_name == 'LAMPEDUSA':  # Water CNR
        Latitude=35.49344
        Longitude=-12.46773
    elif station_name == 'TARTU_temp_sites': # Water tartu
        Latitude=58.26708
        Longitude=26.47054
    #elif (station_name == 'Berre' or station_name == 'BEFR'): # Adriatic Sea
    #    Latitude=43.4484
    #   Longitude=5.1012
    else:
        print('ERROR: station not found: '+station_name)
    return Latitude, Longitude


def qcgen(pixel_classif_flags):
    rows = np.shape(pixel_classif_flags)[0]
    columns = np.shape(pixel_classif_flags)[1]

    # Flags settings
    idepix_flag = np.zeros([rows, columns])

    # IDEPIX Flags
    idepix_flag_masks = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072,
                         262144, 524288, 1048576]
    idepix_flag_meanings = ["IDEPIX_INVALID", "IDEPIX_CLOUD", "IDEPIX_CLOUD_AMBIGUOUS",
                            "IDEPIX_CLOUD_SURE", "IDEPIX_CLOUD_BUFFER", "IDEPIX_CLOUD_SHADOW",
                            "IDEPIX_SNOW_ICE", "IDEPIX_BRIGHT", "IDEPIX_WHITE",
                            "IDEPIX_COASTLINE", "IDEPIX_LAND", "IDEPIX_CIRRUS_SURE",
                            #                           "IDEPIX_COASTLINE", "IDEPIX_CIRRUS_SURE",
                            "IDEPIX_CIRRUS_AMBIGUOUS", "IDEPIX_CLEAR_LAND", "IDEPIX_CLEAR_WATER",
                            "IDEPIX_WATER",
                            "IDEPIX_BRIGHTWHITE", "IDEPIX_VEG_RISK", "IDEPIX_MOUNTAIN_SHADOW",
                            "IDEPIX_POTENTIAL_SHADOW", "IDEPIX_CLUSTERED_CLOUD_SHADOW"]
    idepix_flag_list = ["IDEPIX_INVALID", "IDEPIX_CLOUD", "IDEPIX_CLOUD_AMBIGUOUS", "IDEPIX_CLOUD_SURE",
                        "IDEPIX_CLOUD_BUFFER", "IDEPIX_CLOUD_SHADOW", "IDEPIX_SNOW_ICE",
                        "IDEPIX_COASTLINE", "IDEPIX_LAND", "IDEPIX_CIRRUS_SURE"]  # flags to be applied
    #     idepix_flag_list = ["IDEPIX_VEG_RISK"]  # flags to be applied
    TMP = [idepix_flag_meanings.index(x) if x in idepix_flag_list else 9999999 for x in idepix_flag_meanings]
    INDEX = np.where(np.array(TMP) < 9999999)[0]
    for d in INDEX:
        TEST = idepix_flag_masks[d] & pixel_classif_flags
        idepix_flag[TEST == idepix_flag_masks[d]] = 2
    return idepix_flag, idepix_flag_masks, idepix_flag_meanings


def extract_wind_and_angles(path_source, in_situ_lat, in_situ_lon, file):
    # from Tie-Points grid (a coarser grid)
    filepah = os.path.join(path_source, file)
    nc_sat = Dataset(filepah, 'r')
    tie_lon = nc_sat.variables['lon'][:]
    tie_lat = nc_sat.variables['lat'][:]
    SZA = nc_sat.variables['SZA'][:]
    nc_sat.close()

    r, c = cfs.find_row_column_from_lat_lon(tie_lat, tie_lon, in_situ_lat, in_situ_lon)

    ws0 = 0  # horizontal_wind[r,c,0]
    ws1 = 0  # horizontal_wind[r,c,1]
    sza = SZA[r, c]
    saa = 0#SAA[r, c]
    vza = 0#OZA[r, c]
    vaa = 0#OAA[r, c]

    return ws0, ws1, sza, saa, vza, vaa


def create_extract_MSI(file, options):


    if options.has_option('satellite_options', 'extract_size'):
        try:
            size_box = int(options['satellite_options']['extract_size'])
        except:
            size_box = 25
    else:
        size_box = 25


    #size_box=options["satellite_options"]['extract_size']



    station_name = options["Time_and_sites_selection"]['sites']
    path_output= options["file_path"]['output_dir']
    print(path_output)
    res_str= options["satellite_options"]['resolution']
    proc_version_str= options["satellite_options"]['proc_version']
    processor = options["satellite_options"]['processor']
    path_source=options["file_path"]['sat_source_dir']


    #  if args.verbose:
    print(f'Creating extract for {station_name} from {path_source}')
    filepath = os.path.join(path_source, file)

    try:
        nc_sat = Dataset(filepath)
    except OSError :
        print(f'[ERROR] Unknown file format')
        return None

    lat = nc_sat.variables['lat'][:]
    lon = nc_sat.variables['lon'][:]

    if len(lat.shape)==1:
        lat=np.array(lat)
    if len(lon.shape)==1:
        lon=np.array(lon)

    s2time = [datetime.strptime(nc_sat.__dict__["time_coverage_start"],
                                         '%Y-%m-%dT%H:%M:%S') if "time_coverage_start" in nc_sat.__dict__
              else datetime.strptime(nc_sat.__dict__["start_date"], '%d-%b-%Y %H:%M:%S.%f')]
    print(s2time)

    in_situ_lat, in_situ_lon = get_lat_lon_ins(station_name, s2time[0])

    val_list = file.split('_')
    sat_time = None
    for v in val_list:
        try:
            sat_time = datetime.strptime(v, '%Y%m%dT%H%M%S')
            break
        except ValueError:
            continue

    start_date = sat_time
    ## to be done later
    contain_flag = cfs.contain_location(lat,lon,in_situ_lat,in_situ_lon)
    print(contain_flag)
    if contain_flag==1:
        #         if not args.verbose:
        #             print('-----------------')
        r, c = cfs.find_row_column_from_lat_lon(lat, lon, in_situ_lat, in_situ_lon)

        start_idx_x = (c - int(size_box / 2))
        stop_idx_x = (c + int(size_box / 2) + 1)
        start_idx_y = (r - int(size_box / 2))
        stop_idx_y = (r + int(size_box / 2) + 1)
        print(r)
        print(c)

        if r >= 0 and r + 1 < lat.shape[0] and c >= 0 and c + 1 < lon.shape[0]:


            if processor.lower() in ['acolite', 'ac']:
                try:
                    print('acolite')
                    rhow_443 = nc_sat.variables['Rrs443_a'][:]
                    rhow_492 = nc_sat.variables['Rrs492_a'][:]
                    rhow_560 = nc_sat.variables['Rrs560_a'][:]
                    rhow_665 = nc_sat.variables['Rrs665_a'][:]
                    rhow_704 = nc_sat.variables['Rrs704_a'][:]
                    rhow_740 = nc_sat.variables['Rrs740_a'][:]
                    rhow_783 = nc_sat.variables['Rrs783_a'][:]
                    rhow_865 = nc_sat.variables['Rrs865_a'][:]
                    nbands = 8
                except KeyError:
                    return
            if processor.lower() in ['c2rcc']:
                try:
                    print('c2rcc')
                    rhow_443 = nc_sat.variables['Rrs443_c'][:, :]
                    rhow_492 = nc_sat.variables['Rrs492_c'][:, :]
                    rhow_560 = nc_sat.variables['Rrs560_c'][:, :]
                    rhow_665 = nc_sat.variables['Rrs665_c'][:, :]
                    rhow_704 = nc_sat.variables['Rrs704_c'][:, :]
                    rhow_740 = nc_sat.variables['Rrs740_c'][:, :]
                    rhow_783 = nc_sat.variables['Rrs783_c'][:, :]
                    rhow_865 = nc_sat.variables['Rrs865_c'][:, :]
                    nbands = 8
                except KeyError:
                    return

            # add idepix flags
            pixel_classif_flags = nc_sat.variables['pixel_classif_flags'][:, :]
            idepix_flag, idepix_flag_mask, idepix_flag_meanings = qcgen(pixel_classif_flags)
            #             loc = np.where(idepix_flag>0)

            WQSF = nc_sat.variables['pixel_classif_flags'][:]
            WQSF_flag_masks = idepix_flag_mask
            WQSF_flag_meanings = idepix_flag_meanings
            nc_sat.close()

            # %% Calculate BRDF
            ws0, ws1, sza, saa, vza, vaa = extract_wind_and_angles(path_source, in_situ_lat, in_situ_lon, file)

            # %% Save extract as netCDF4 file
            filebase="{}/{}".format(path_output, file).split('__ACI')
            ofname = filebase[0] + '_extract_' + station_name + '_' + processor + filebase[1]
            #ofname = os.path.join(path_output, filename)
            print(ofname)

            if os.path.exists(ofname):
                os.remove(ofname)

            new_EXTRACT = Dataset(ofname, 'w', format='NETCDF4')
            new_EXTRACT.creation_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")

            # MAKE A DEF WITH READ METADATA AND EXTRACT THE REQUIRED METADATA OR FROM FILENAME

            satellite = file[0:2]
            platform = file[2]
            sensor = file[4:7]

            # read the metadata and add the SPACECRAFT_NAME as satellite
            # '    Level-1C_DataStrip_ID:General_Info:Datatake_Info:SPACECRAFT_NAME: Sentinel-2A',

            new_EXTRACT.satellite = satellite
            new_EXTRACT.platform = platform

            new_EXTRACT.sensor = sensor
            new_EXTRACT.description = f'{satellite}{platform} {sensor.upper()} {res_str} L-1C extract'
            # new_EXTRACT.satellite_PDU = path_source.split('/')[-1]
            # new_EXTRACT.satellite_path_source = path_source
            new_EXTRACT.satellite_aco_processor = 'Atmospheric Correction processor: ' + processor

            # MAKE A DEF WITH READ METADATA AND EXTRACT THE REQUIRED METADATA
            # read the metadata and add the SPACECRAFT_NAME as satellite
            new_EXTRACT.satellite_proc_version = proc_version_str

            # new_EXTRACT.datapolicy = 'Notice to users: Add data policy'
            # new_EXTRACT.insitu_sensor_processor_version = '0.0'
            new_EXTRACT.insitu_site_name = station_name

            new_EXTRACT.insitu_lat = in_situ_lat
            new_EXTRACT.insitu_lon = in_situ_lon

            #new_EXTRACT.satellite_ws0 = ws0
            #new_EXTRACT.satellite_ws1 = ws1
            new_EXTRACT.satellite_SZA_center_pixel = sza
            #new_EXTRACT.satellite_SAA_center_pixel = saa
            #new_EXTRACT.satellite_VZA_center_pixel = vza
            #new_EXTRACT.satellite_VAA_center_pixel = vaa

            # dimensions
            new_EXTRACT.createDimension('satellite_id', None)
            new_EXTRACT.createDimension('rows', size_box)
            new_EXTRACT.createDimension('columns', size_box)
            new_EXTRACT.createDimension('satellite_bands', nbands)

            filepah = os.path.join(path_source, file)
            nc_sat = Dataset(filepah, 'r')
            tie_lon = nc_sat.variables['lon'][:]
            tie_lat = nc_sat.variables['lat'][:]
            SZA = nc_sat.variables['SZA'][:]
            nc_sat.close()

            # variables
            satellite_SZA = new_EXTRACT.createVariable('SZA', 'f4', ('satellite_id', 'rows', 'columns'),
                                                       fill_value=-999, zlib=True, complevel=6)
            satellite_SZA[0, :, :] = [SZA[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]]
            satellite_SZA.long_name = 'Sun Zenith Angle'
            satellite_SZA.units = 'degrees'

            # satellite_SAA = new_EXTRACT.createVariable('SAA', 'f4', ('satellite_id', 'rows', 'columns'),
            #                                            fill_value=-999, zlib=True, complevel=6)
            # satellite_SAA[0, :, :] = [SAA[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]]
            # satellite_SAA.long_name = 'Sun Azimuth Angle'
            # satellite_SAA.units = 'degrees'
            #
            # satellite_OZA = new_EXTRACT.createVariable('OZA', 'f4', ('satellite_id', 'rows', 'columns'),
            #                                            fill_value=-999, zlib=True, complevel=6)
            # satellite_OZA[0, :, :] = [OZA[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]]
            # satellite_OZA.long_name = 'Observation Zenith Angle'
            # satellite_OZA.units = 'degrees'
            #
            # satellite_OAA = new_EXTRACT.createVariable('OAA', 'f4', ('satellite_id', 'rows', 'columns'),
            #                                            fill_value=-999, zlib=True, complevel=6)
            # satellite_OAA[0, :, :] = [OAA[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]]
            # satellite_OAA.long_name = 'Observation Azimuth Angle'
            # satellite_OAA.units = 'degrees'

            satellite_time = new_EXTRACT.createVariable('satellite_time', 'f4', ('satellite_id'), fill_value=-999,
                                                        zlib=True, complevel=6)
            # ValueError: time data '2021-03-03T10:48:50.106616+00:00' does not match format '%d-%b-%Y %H:%M:%S.%f'
            satellite_time[0] = float(start_date.timestamp())
            #float(datetime.strptime(start_date, "%Y-%m-%dT%H:%M:%S.%f+00:00").timestamp())
            satellite_time.units = "Seconds since 1970-1-1"

#            satellite_PDU = new_EXTRACT.createVariable('satellite_PDU',str, ('satellite_id'),fill_value=-999, zlib=True,complevel=6)
#            satellite_PDU[0] = "{}/{}".format(path_source, file).split('__ACI.nc')[0]
#            satellite_PDU.long_name = "S2A source PDU name"

            satellite_latitude = new_EXTRACT.createVariable('satellite_latitude', 'f4',
                                                            ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                                            zlib=True, complevel=6)

            if len(lat.shape)==2:
                satellite_latitude[0, :, :] = [lat[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]]
                satellite_latitude.short_name = 'latitude'
            else:
                satellite_latitude[0, :, 0] = [lat[start_idx_y:stop_idx_y]]
                satellite_latitude.short_name = 'latitude'

            satellite_longitude = new_EXTRACT.createVariable('satellite_longitude', 'f4',
                                                             ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                                             zlib=True, complevel=6)
            if len(lon.shape)==2:
                satellite_longitude[0, :, :] = [lon[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]]
                satellite_longitude.short_name = 'longitude'
            else:
                satellite_longitude[0, 0, :] = [lon[start_idx_x:stop_idx_x]]
                satellite_longitude.short_name = 'longitude'

            S2Awl = [442.7, 492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 945.1, 1373.5, 1613.7,
                     2202.4]  # sentinel2 wavelengths original

            S2Awl_c2rcc = [442.7, 492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 864.7]  # no 832.8 wavelength

            S2Awl_ac = [442.7, 492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 1613.7,
                        2202.4]  # no 945.1 wavelength

            # double satellite_bands          (satellite_bands) ;
            satellite_bands = new_EXTRACT.createVariable('satellite_bands', 'f4', ('satellite_bands'), fill_value=-999,
                                                         zlib=True, complevel=6)
            satellite_bands.units = 'nm'
            for i in range(nbands):
                if nbands == 8:
                    satellite_bands[i] = S2Awl_c2rcc[i]
                else:
                    if nbands == 11:
                        satellite_bands[i] = S2Awl_ac[i]  # must check the rhow wavelength values from ac and c2rcc

            #             S2Awl=[442.7, 492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 945.1, 1373.5, 1613.7, 2202.4]

            #             S2Bwl=[442.3, 492.1, 559, 665, 703.8, 739.1, 779.7, 833, 864, 943.2, 1376.9, 1610.4, 2185.7]

            # NOT BRDF-corrected
            satellite_Rrs = new_EXTRACT.createVariable('satellite_Rrs', 'f4',
                                                       ('satellite_id', 'satellite_bands', 'rows', 'columns'),
                                                       fill_value=-999, zlib=True, complevel=6)
            if nbands == 8:
                satellite_Rrs[0, 0, :, :] = [ma.array(rhow_443[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) ]
                satellite_Rrs[0, 1, :, :] = [ma.array(rhow_492[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) ]
                satellite_Rrs[0, 2, :, :] = [ma.array(rhow_560[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) ]
                satellite_Rrs[0, 3, :, :] = [ma.array(rhow_665[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) ]
                satellite_Rrs[0, 4, :, :] = [ma.array(rhow_704[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) ]
                satellite_Rrs[0, 5, :, :] = [ma.array(rhow_740[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) ]
                satellite_Rrs[0, 6, :, :] = [ma.array(rhow_783[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) ]
                satellite_Rrs[0, 7, :, :] = [ma.array(rhow_865[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]) ]

            satellite_Rrs.short_name = 'Satellite Rrs.'
            satellite_Rrs.long_name = "Above water Remote Sensing Reflectance for S2-MSI acquisition without BRDF correction applied";
            satellite_Rrs.units = "sr-1";

            satellite_WQSF = new_EXTRACT.createVariable('satellite_WQSF', 'f4', ('satellite_id', 'rows', 'columns'),
                                                        fill_value=-999, zlib=True, complevel=6)
            satellite_WQSF[0, :, :] = [ma.array(WQSF[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])]
            satellite_WQSF.description = 'Satellite Level 1C WATER Product, Classification, Quality and Science Flags Data Set'
            satellite_WQSF.flag_masks = WQSF_flag_masks
            satellite_WQSF.flag_meanings = WQSF_flag_meanings

            new_EXTRACT.close()
            print('Extract created!')
            return ofname

        else:
            print('Warning: Index out of bound!')

    else:
        if args.verbose:
            print('Warning: File does NOT contains the in situ location!')


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



def get_find_product_info(options):
    path_source = options['file_path']['sat_source_dir']
    org = None
    if options.has_option('file_path', 'sat_source_dir_organization') and options['file_path'][
        'sat_source_dir_organization']:
        org = options['file_path']['sat_source_dir_organization']
    wce = '*'
    if options.has_option('file_path', 'wce') and options['file_path']['wce']:
        wce = options['file_path']['wce']

    return path_source, org, wce

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

def main():
    print('[INFO]Creating satellite extracts')

    if not args.config_file:
        return
    options = config_reader(args.config_file)
    print('------------------------------')
    path_source, org, wce = get_find_product_info(options)

    product_path_list = [f for f in listdir(path_source) if isfile(join(path_source, f))]

    if len(product_path_list) == 0:
        if args.verbose:
            print('-----------------')
        print(f'[WARNING] No valid datasets found on:  {path_source}')
        print('------------------------------')
        print('[INFO] COMPLETED. 0 sat extract files were created')
        return
    ncreated = []
    for filepath in product_path_list:
        print(f'[INFO]Extracting sat extract for product: {filepath}')
        nhere = create_extract_MSI(filepath, options)
        ncreated.append(nhere)
    print('------------------------------')
    print(f'COMPLETED. {len(ncreated)} sat extract files were created')


if __name__ == '__main__':
    main()