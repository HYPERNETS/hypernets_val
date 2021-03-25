#!/usr/bin/env python3
# coding: utf-8
"""
Created on Tue Jun 16 12:02:40 2020
Create Matchup Data Base (MDB) file (NetCDF4 file)
@author: javier.concha

Based on EUMETSAT MDB_Builder module (https://ocdb.readthedocs.io/en/latest/ocdb-MDB-user-manual.html)

Run as:

python MDB_builder.py -p PATH_TO_SAT_SOURCES

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
"""
Sat:
- filename (path to zip or unzip?)
- datetime (1)
- xsize
- ysize
- bands
dimensions: (xsize,ysize,bands)

In situ:
- station name OR lat, lon in situ
- lat in situ
- lon in situ
- time window
- datetime (N)
- bands
dimensions: (bands,datetime)

from $ ncdump -h MDB_S3A_OLCI_L2_AERONET_Venise.nc:
dimensions:
    satellite_id = UNLIMITED ; // (456 currently)
    rows = 25 ;
    columns = 25 ;
    wind_vectors = 2 ;
    satellite_detectors = 3700 ;
    satellite_BRDF_Bands = 11 ;
    satellite_bands = 16 ;
    insitu_id = 24 ;
    insitu_original_bands = 29 ;
    insitu_Rrs_bands = 16 ;

(base) 203$ ncdump -h extract_Venise.nc 
netcdf extract_Venise {
dimensions:
    satellite_id = UNLIMITED ; // (1 currently)
    size_box_x = 25 ;
    size_box_y = 25 ;   

"""
#%% imports
import os
import sys
import subprocess

from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from datetime import datetime
from datetime import timedelta
import pandas as pd

# import user defined functions from other .py
code_home = os.path.abspath('../')
sys.path.append(code_home)

import BRDF.brdf_olci as brdf
import COMMON.common_functions as cfs

os.environ['QT_QPA_PLATFORM']='offscreen' # to avoid error "QXcbConnection: Could not connect to display"

import argparse
parser = argparse.ArgumentParser(description="Create Match-up DataBase files (MDB) files.")
parser.add_argument("-d", "--debug", help="Debugging mode.",action="store_true")
parser.add_argument("-t", "--test", help="Test mode.",action="store_true")
parser.add_argument('-sd', "--startdate", help="The Start Date - format YYYY-MM-DD ")
parser.add_argument('-ed', "--enddate", help="The End Date - format YYYY-MM-DD ")
parser.add_argument('-site', "--sitename", help="Site name.",required=True,choices=['VEIT','BEFR','BSBE'])
parser.add_argument('-ins', "--insitu", help="Satellite sensor name.",required=True,choices=['PANTHYR', 'HYPERNETS']) # ,'HYPSTAR'])
parser.add_argument('-pi', "--path_to_ins", help="Path to in situ sources.",required=True)
parser.add_argument('-sat', "--satellite", help="Satellite sensor name.",choices=['OLCI', 'MSI'])
parser.add_argument('-ps', "--path_to_sat", help="Path to satellite sources.")
parser.add_argument('-o', "--output", help="Path to output")
parser.add_argument('-res', "--resolution", help="Resolution OL_2: WRR or WFR (for OLCI)")
parser.add_argument('-nl', "--nolist", help="Do not create satellite and in situ lists.",action="store_true")

args = parser.parse_args()

# user defined functions
def create_list_products(path_source,path_out,wce,res_str,type_product):
    path_to_list = f'{path_out}/file_{type_product}_{res_str}_list.txt'
    if not args.nolist:
        cmd = f'find {path_source} -name {wce}|sort|uniq> {path_to_list}'
        prog = subprocess.Popen(cmd, shell=True,stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(err)  
    
    return path_to_list

def extract_wind_and_angles(path_source,in_situ_lat,in_situ_lon):
    # from Tie-Points grid (a coarser grid)
    filepah = os.path.join(path_source,'tie_geo_coordinates.nc')
    nc_f0 = Dataset(filepah,'r')
    tie_lon = nc_f0.variables['longitude'][:]
    tie_lat = nc_f0.variables['latitude'][:]
    
    filepah = os.path.join(path_source,'tie_meteo.nc')
    nc_f0 = Dataset(filepah,'r')
    horizontal_wind = nc_f0.variables['horizontal_wind'][:]
    nc_f0.close()
    
    filepah = os.path.join(path_source,'tie_geometries.nc')
    nc_f1 = Dataset(filepah,'r')
    SZA = nc_f1.variables['SZA'][:]    
    SAA = nc_f1.variables['SAA'][:]  
    OZA = nc_f1.variables['OZA'][:]  
    OAA = nc_f1.variables['OAA'][:] 
    nc_f1.close()        
      
    r, c = cfs.find_row_column_from_lat_lon(tie_lat,tie_lon,in_situ_lat,in_situ_lon)
    
    ws0 = horizontal_wind[r,c,0]
    ws1 = horizontal_wind[r,c,1]  
    sza = SZA[r,c]
    saa = SAA[r,c]
    vza = OZA[r,c]
    vaa = OAA[r,c]

    return ws0, ws1, sza, saa, vza, vaa

def create_extract(size_box,station_name,path_source,path_output,in_situ_lat,in_situ_lon,res_str,insitu_sensor):
    # extract IFP-OL-2 version
    with open(os.path.join(path_source,'xfdumanifest.xml'),'r', encoding="utf-8") as read_obj:
        check_version = False
        for line in read_obj:
            if 'IPF-OL-2' in line and check_version == False:
                IPF_OL_2_version = line.split('"')[3]
                proc_version_str = f'IPF-OL-2 version {IPF_OL_2_version}'
                print(proc_version_str)
                check_version = True
                pass

    
    #% open nc file
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
    
    filepah = os.path.join(path_source,coordinates_filename)
    nc_f0 = Dataset(filepah,'r')
    
    lat = nc_f0.variables['latitude'][:,:]
    lon = nc_f0.variables['longitude'][:,:]
    
    
    
    contain_flag = cfs.contain_location(lat,lon,in_situ_lat,in_situ_lon)
    
    if contain_flag:
        if not args.debug:
            print('-----------------')
        r, c = cfs.find_row_column_from_lat_lon(lat,lon,in_situ_lat,in_situ_lon)
        
        start_idx_x = (r-int(size_box/2))
        stop_idx_x = (r+int(size_box/2)+1)
        start_idx_y = (c-int(size_box/2))
        stop_idx_y = (c+int(size_box/2)+1)
    
        
        if r>=0 and r+1<lat.shape[0] and c>=0 and c+1<lat.shape[1]:
            # read nc file
            filepah = os.path.join(path_source,rhow_0400p00_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0400p00 = nc_f1.variables['Oa01_reflectance'][:]
            nc_f1.close()

            filepah = os.path.join(path_source,rhow_0412p50_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0412p50 = nc_f1.variables['Oa02_reflectance'][:]
            nc_f1.close()
            
            filepah = os.path.join(path_source,rhow_0442p50_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0442p50 = nc_f1.variables['Oa03_reflectance'][:]
            nc_f1.close()
            
            filepah = os.path.join(path_source,rhow_0490p00_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0490p00 = nc_f1.variables['Oa04_reflectance'][:]
            nc_f1.close()
            
            filepah = os.path.join(path_source,rhow_0510p00_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0510p00 = nc_f1.variables['Oa05_reflectance'][:]
            nc_f1.close()
            
            filepah = os.path.join(path_source,rhow_0560p00_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0560p00 = nc_f1.variables['Oa06_reflectance'][:]
            nc_f1.close()
    
            filepah = os.path.join(path_source,rhow_0620p00_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0620p00 = nc_f1.variables['Oa07_reflectance'][:]
            nc_f1.close()
            
            filepah = os.path.join(path_source,rhow_0665p00_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0665p00 = nc_f1.variables['Oa08_reflectance'][:]
            nc_f1.close()
            
            filepah = os.path.join(path_source,rhow_0673p75_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0673p75 = nc_f1.variables['Oa09_reflectance'][:]
            nc_f1.close()

            filepah = os.path.join(path_source,rhow_0681p25_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0681p25 = nc_f1.variables['Oa10_reflectance'][:]
            nc_f1.close()

            filepah = os.path.join(path_source,rhow_0708p75_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0708p75 = nc_f1.variables['Oa11_reflectance'][:]
            nc_f1.close() 

            filepah = os.path.join(path_source,rhow_0753p75_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0753p75 = nc_f1.variables['Oa12_reflectance'][:]
            nc_f1.close()

            filepah = os.path.join(path_source,rhow_0778p75_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0778p75 = nc_f1.variables['Oa16_reflectance'][:]
            nc_f1.close()          
            
            filepah = os.path.join(path_source,rhow_0865p00_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0865p00 = nc_f1.variables['Oa17_reflectance'][:]
            nc_f1.close()

            filepah = os.path.join(path_source,rhow_0885p00_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_0885p00 = nc_f1.variables['Oa18_reflectance'][:]
            nc_f1.close()
            
            filepah = os.path.join(path_source,rhow_1020p50_filename)
            nc_f1 = Dataset(filepah,'r')
            rhow_1020p50 = nc_f1.variables['Oa21_reflectance'][:]
            nc_f1.close()
            
            filepah = os.path.join(path_source,AOT_0865p50_filename)
            nc_f1 = Dataset(filepah,'r')
            AOT_0865p50 = nc_f1.variables['T865'][:]
            nc_f1.close()
    
            filepah = os.path.join(path_source,WQSF_filename)
            nc_f1 = Dataset(filepah,'r')
            WQSF = nc_f1.variables['WQSF'][:]
            WQSF_flag_masks = nc_f1.variables['WQSF'].flag_masks
            WQSF_flag_meanings = nc_f1.variables['WQSF'].flag_meanings
            nc_f1.close()

            #%% Calculate BRDF
            ws0, ws1, sza, saa, vza, vaa = extract_wind_and_angles(path_source,in_situ_lat,in_situ_lon)
            
            filepah = os.path.join(path_source,'chl_oc4me.nc')
            nc_f11 = Dataset(filepah,'r')
            CHL_OC4ME = nc_f11.variables['CHL_OC4ME'][:]
            nc_f11.close()
            CHL_OC4ME_extract = ma.array(CHL_OC4ME[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])
            BRDF0 = np.full(CHL_OC4ME_extract.shape,np.nan)
            BRDF1 = np.full(CHL_OC4ME_extract.shape,np.nan)
            BRDF2 = np.full(CHL_OC4ME_extract.shape,np.nan)
            BRDF3 = np.full(CHL_OC4ME_extract.shape,np.nan)
            BRDF4 = np.full(CHL_OC4ME_extract.shape,np.nan)
            BRDF5 = np.full(CHL_OC4ME_extract.shape,np.nan)
            BRDF6 = np.full(CHL_OC4ME_extract.shape,np.nan)
            for ind0 in range(CHL_OC4ME_extract.shape[0]):
                for ind1 in range(CHL_OC4ME_extract.shape[1]):
                    chl = CHL_OC4ME_extract[ind0,ind1]
                    # 412.5, 442.5, 490, 510, 560, 620, 660 bands
                    # 0      1      2    3    4    5    6   brdf index
                    # 412.5  442.5  490  510  560  620  665 OLCI bands
                    # 02     03     04   05   06   07   08  OLCI band names in L2
                    brdf_coeffs = brdf.brdf(ws0, ws1, chl, sza, saa, vza, vaa)
                    BRDF0[ind0,ind1] = brdf_coeffs[0,0]
                    BRDF1[ind0,ind1] = brdf_coeffs[0,1]
                    BRDF2[ind0,ind1] = brdf_coeffs[0,2]
                    BRDF3[ind0,ind1] = brdf_coeffs[0,3]
                    BRDF4[ind0,ind1] = brdf_coeffs[0,4]
                    BRDF5[ind0,ind1] = brdf_coeffs[0,5]
                    BRDF6[ind0,ind1] = brdf_coeffs[0,6]

            #%% Save extract as netCDF4 file
            path_out = os.path.join(path_output,'EXTRACTS')
            filename = path_source.split('/')[-1].replace('.','_')+'_extract_'+station_name+'.nc'
            ofname = os.path.join(path_out,filename)
            
            print(filename)

            satellite = filename[0:2]
            platform = filename[2]
            sensor = 'olci'
            
            if os.path.exists(ofname):
              os.remove(ofname)
            
            fmb = Dataset(ofname, 'w', format='NETCDF4')
            fmb.MDB_software_version = '0.0'
            fmb.creation_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")            
            fmb.satellite = satellite
            fmb.platform = platform
            fmb.sensor = sensor
            fmb.description = f'{satellite}{platform} {sensor.upper()} {res_str} L2 - {insitu_sensor} Matchup Data Base'
            fmb.satellite_start_time = nc_f0.start_time
            fmb.satellite_stop_time = nc_f0.stop_time    
            fmb.satellite_PDU = path_source.split('/')[-1]
            fmb.satellite_path_source = path_source
            fmb.satellite_aco_processor = 'Atmospheric Correction processor: xxx'
            fmb.satellite_proc_version = proc_version_str
            fmb.satellite_WQSF_flag_masks = WQSF_flag_masks
            fmb.satellite_WQSF_flag_meanings = WQSF_flag_meanings

            fmb.datapolicy = 'Notice to users: Add data policy'
            fmb.insitu_sensor_processor_version = '0.0'
            fmb.insitu_site_name = station_name

            fmb.insitu_lat = in_situ_lat
            fmb.insitu_lon = in_situ_lon

            fmb.satellite_ws0 = ws0
            fmb.satellite_ws1 = ws1
            fmb.satellite_SZA_center_pixel = sza
            fmb.satellite_SAA_center_pixel = saa
            fmb.satellite_VZA_center_pixel = vza
            fmb.satellite_VAA_center_pixel = vaa
            
            # dimensions
            fmb.createDimension('satellite_size_box_x', size_box)
            fmb.createDimension('satellite_size_box_y', size_box)
            fmb.createDimension('satellite_bands', 16)
            fmb.createDimension('satellite_BRDF_bands', 7)
            
            # variables  
            # satellite_SZA = fmb.createVariable('satellite_SZA', 'f4', ('satellite_size_box_x','satellite_size_box_y'), fill_value=-999, zlib=True, complevel=6)
            # satellite_SZA[:] = [SZA[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]]
            # satellite_SZA.long_name = 'Sun Zenith Angle'
            # satellite_SZA.long_name = 'Sun Zenith Angle'    
            # satellite_SZA.units = 'degrees'
# SZA
# SAA
# OZA
# OAA

            satellite_latitude = fmb.createVariable('satellite_latitude',  'f4', ('satellite_size_box_x','satellite_size_box_y'), fill_value=-999, zlib=True, complevel=6) 
            satellite_latitude[:,:] = [lat[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]]
            
            satellite_longitude = fmb.createVariable('satellite_longitude',  'f4', ('satellite_size_box_x','satellite_size_box_y'), fill_value=-999, zlib=True, complevel=6)
            satellite_longitude[:,:] = [lon[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]]

            # double satellite_bands          (satellite_bands) ;
            satellite_bands = fmb.createVariable('satellite_bands',  'f4', ('satellite_bands'), fill_value=-999, zlib=True, complevel=6) 
            satellite_bands[:] = [0400.00,0412.50,0442.50,0490.00,0510.00,0560.00,0620.00,0665.00,0673.75,0681.25,0708.75,0753.75,0778.75,0865.00,0885.00,1020.50]

            # double satellite_BRDF_Bands     (satellite_BRDF_Bands) ;
            satellite_BRDF_bands = fmb.createVariable('satellite_BRDF_bands',  'f4', ('satellite_BRDF_bands'), fill_value=-999, zlib=True, complevel=6) 
            satellite_BRDF_bands[:] = [412.50,442.50,490.00,510.00,560.00,620.00,665.00]
    
            # NOT BRDF-corrected
            satellite_rhow = fmb.createVariable('satellite_rhow', 'f4', ('satellite_bands','satellite_size_box_x','satellite_size_box_y'), fill_value=-999, zlib=True, complevel=6)
            satellite_rhow[0,:,:] = [ma.array(rhow_0400p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[1,:,:] = [ma.array(rhow_0412p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[2,:,:] = [ma.array(rhow_0442p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[3,:,:] = [ma.array(rhow_0490p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[4,:,:] = [ma.array(rhow_0510p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[5,:,:] = [ma.array(rhow_0560p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[6,:,:] = [ma.array(rhow_0620p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[7,:,:] = [ma.array(rhow_0665p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[8,:,:] = [ma.array(rhow_0673p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[9,:,:] = [ma.array(rhow_0681p25[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[10,:,:] = [ma.array(rhow_0708p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[11,:,:] = [ma.array(rhow_0753p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[12,:,:] = [ma.array(rhow_0778p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[13,:,:] = [ma.array(rhow_0865p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[14,:,:] = [ma.array(rhow_0885p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[15,:,:] = [ma.array(rhow_1020p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow.description = 'Satellite rhow.'
            
            # BRDF-corrected
            satellite_BRDF_rhow = fmb.createVariable('satellite_BRDF_rhow', 'f4', ('satellite_BRDF_bands','satellite_size_box_x','satellite_size_box_y'), fill_value=-999, zlib=True, complevel=6)
            satellite_BRDF_rhow[0,:,:] = [ma.array(rhow_0412p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF0)]
            satellite_BRDF_rhow[1,:,:] = [ma.array(rhow_0442p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF1)]
            satellite_BRDF_rhow[2,:,:] = [ma.array(rhow_0490p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF2)]
            satellite_BRDF_rhow[3,:,:] = [ma.array(rhow_0510p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF3)]
            satellite_BRDF_rhow[4,:,:] = [ma.array(rhow_0560p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF4)]
            satellite_BRDF_rhow[5,:,:] = [ma.array(rhow_0620p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF5)]
            satellite_BRDF_rhow[6,:,:] = [ma.array(rhow_0665p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF6)]
            satellite_BRDF_rhow.description = 'Satellite rhow BRDF-corrected'
            
            satellite_AOT_0865p50_box = fmb.createVariable('satellite_AOT_0865p50', 'f4', ('satellite_size_box_x','satellite_size_box_y'), fill_value=-999, zlib=True, complevel=6)
            satellite_AOT_0865p50_box[:,:] = [ma.array(AOT_0865p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_AOT_0865p50_box.description = 'Satellite Aerosol optical thickness'
    
            satellite_WQSF = fmb.createVariable('satellite_WQSF', 'f4', ('satellite_size_box_x','satellite_size_box_y'), fill_value=-999, zlib=True, complevel=6)
            satellite_WQSF[:,:] = [ma.array(WQSF[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_WQSF.description = 'Satellite Level 2 WATER Product, Classification, Quality and Science Flags Data Set'
            
            satellite_BRDF_fQ = fmb.createVariable('satellite_BRDF_fQ', 'f4', ('satellite_BRDF_bands','satellite_size_box_x','satellite_size_box_y'), fill_value=-999, zlib=True, complevel=6)
            satellite_BRDF_fQ[0,:,:] = [ma.array(BRDF0)]
            satellite_BRDF_fQ[1,:,:] = [ma.array(BRDF1)]
            satellite_BRDF_fQ[2,:,:] = [ma.array(BRDF2)]
            satellite_BRDF_fQ[3,:,:] = [ma.array(BRDF3)]
            satellite_BRDF_fQ[4,:,:] = [ma.array(BRDF4)]
            satellite_BRDF_fQ[5,:,:] = [ma.array(BRDF5)]
            satellite_BRDF_fQ[6,:,:] = [ma.array(BRDF6)]
            satellite_BRDF_fQ.description = 'Satellite BRDF fQ coefficients'

            satellite_chl_oc4me = fmb.createVariable('chl_oc4me', 'f4', ('satellite_size_box_x','satellite_size_box_y'), fill_value=-999, zlib=True, complevel=6)
            satellite_chl_oc4me[:,:] = [ma.array(CHL_OC4ME_extract)]
            satellite_chl_oc4me.description = 'Satellite Chlorophyll-a concentration from OC4ME.'

            fmb.close()
            # print('Extract created!')
                
        else:
            print('Index out of bound!')
    else:
        if args.debug:
            print('File does NOT contains the in situ location!')
    
    return ofname

def copy_nc(ifile,ofile):
    with Dataset(ifile) as src:
        dst = Dataset(ofile, "w")
        
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

def add_OL_12_to_list(path_source,path_output,res_str):
    # create OL_1 and OL_2 and OL_2 trimmed lists
    if res_str == 'WFR':
        res_L1_str = 'EFR' 
    elif res_str == 'WRR':
        res_L1_str = 'ERR'
        
    cmd = f'cat {path_source}/xfdumanifest.xml | grep OL_2_{res_str}|grep -v product|cut -d '+"'"+'"'+"'"+f' -f2>> {path_output}/OL_2_{res_str}_list.txt'  
    prog = subprocess.Popen(cmd, shell=True,stderr=subprocess.PIPE)
    out, err = prog.communicate()
    if err:
        print(err)  
    # elif args.debug:
    #     print('Run:')
    #     print(cmd)
    
    cmd = f'echo {path_source.split("/")[-1]}>> {path_output}/OL_2_{res_str}_list.txt'
    prog = subprocess.Popen(cmd, shell=True,stderr=subprocess.PIPE)
    out, err = prog.communicate()
    if err:
        print(err) 
    # elif args.debug:
    #     print('Run:')
    #     print(cmd)
        
    cmd = f'cat {path_source}/xfdumanifest.xml | grep OL_1_{res_L1_str}|cut -d '+"'"+'"'+"'"+f' -f2>> {path_output}/OL_1_{res_L1_str}_list.txt'  
    prog = subprocess.Popen(cmd, shell=True,stderr=subprocess.PIPE)
    out, err = prog.communicate()
    if err:
        print(err) 
    # elif args.debug:
    #     print('Run:')
    #     print(cmd)
        
def clean_lists(path_out,res_str):
    list_path = f'{path_out}/OL_2_{res_str}_list.txt'
    sort_uniq(list_path,path_out)
    
    list_path = f'{path_out}/OL_2_{res_str}_list.txt'
    sort_uniq(list_path,path_out)
    
    if res_str == 'WFR':
        res_L1_str = 'EFR'
    elif res_str == 'WRR':
        res_L1_str = 'ERR'

    list_path = f'{path_out}/OL_1_{res_L1_str}_list.txt'
    sort_uniq(list_path,path_out)
    
def sort_uniq(list_path,path_out): 
    if os.path.exists(list_path):
        cmd = f'sort {list_path}|uniq > {path_out}/temp_list.txt && {path_out}/temp_list.txt {list_path}'  
        prog = subprocess.Popen(cmd, shell=True,stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(err) 
        os.remove(f'{path_out}/temp_list.txt')   
   
def create_insitu_list_daily(path_to_insitu_list,date_str):
    # create list in situ per date YYYYMMDD
    YYYYMMDD_str = date_str[:8]
    path_to_list = f'{path_to_insitu_list[:-4]}_{YYYYMMDD_str}.txt'
    cmd = f'cat {path_to_insitu_list}|grep {YYYYMMDD_str}|sort|uniq> {path_to_list}'
    prog = subprocess.Popen(cmd, shell=True,stderr=subprocess.PIPE)
    out, err = prog.communicate()
    if err:
        print(err)
        
    return path_to_list
 
    
def add_insitu(extract_path,ofile,path_to_list_daily,datetime_str,time_window):
    # print(f'Satellite time {datetime_str}')
    date_format = '%Y%m%dT%H%M%S'
    satellite_datetime = datetime.strptime(datetime_str, date_format)
      
    # to append to nc file
    nc_f0 = copy_nc(extract_path,ofile)

    # add time window diff
    nc_f0.time_diff = f'{time_window*60*60}' # in seconds
    
    # create in situ dimensions
    nc_f0.createDimension('insitu_id', None)
    nc_f0.createDimension('insitu_bands', None)
    
    # create variable 
    insitu_time=nc_f0.createVariable('insitu_time', 'S2', ('insitu_id'), zlib=True, complevel=6)
    insitu_time.description  = 'In situ time in ISO 8601 format (UTC).'
    
    insitu_filename=nc_f0.createVariable('insitu_filename', 'S2', ('insitu_id'), zlib=True, complevel=6)
    insitu_filename.description  = 'In situ filename.'
    
    insitu_filepath=nc_f0.createVariable('insitu_filepath', 'S2', ('insitu_id'), zlib=True, complevel=6)
    insitu_filepath.description  = 'In situ file path.'
    
    insitu_bands=nc_f0.createVariable('insitu_bands', 'f4', ('insitu_bands'), fill_value=-999, zlib=True, complevel=6)
    insitu_bands.description  = 'In situ bands in nm.'
    
    insitu_rhow=nc_f0.createVariable('insitu_rhow', 'f4', ('insitu_bands','insitu_id'), fill_value=-999, zlib=True, complevel=6)
    insitu_rhow.description  = 'In situ rhow.'

    time_difference=nc_f0.createVariable('time_difference', 'f4', ('insitu_bands'), fill_value=-999, zlib=True, complevel=6)
    time_difference.long_name = "Absolute time difference between satellite acquisition and AERONET acquisition"
    time_difference.units = "seconds"
     
    insitu_idx = 0
    # extract in situ data
    with open(path_to_list_daily) as file:
        for idx, line in enumerate(file):
            # HYPERNETS ex: HYPERNETS_W_BEFR_L2A_REF_202103151201_202103231711_v1.1.nc
            # PANTHYR ex: AAOT_20190923_019046_20200923_104022_AZI_225_data.csv
            # time in ISO 8601 format (UTC). Ex: "2020-01-05T09:27:27.934965Z"
            if args.insitu == 'PANTHYR':             
                date_str = os.path.basename(line[:-1]).split('_')[3]
                time_str = os.path.basename(line[:-1]).split('_')[4]
            elif args.insitu == 'HYPERNETS':
                date_str = os.path.basename(line[:-1]).split('_')[5][0:8]
                time_str = os.path.basename(line[:-1]).split('_')[5][-4:]+'00'            
                     
            YYYY_str = date_str[:4]
            MM_str = date_str[4:6]
            DD_str = date_str[6:8]
            HH_str = time_str[:2]
            mm_str = time_str[2:4]
            ss_str = time_str[4:6]
            insitu_datetime_str = f'{YYYY_str}-{MM_str}-{DD_str}T{HH_str}:{mm_str}:{ss_str}Z'
            insitu_datetime = datetime(int(YYYY_str),int(MM_str),int(DD_str),int(HH_str),int(mm_str),int(ss_str))
            
            time_diff = (insitu_datetime-satellite_datetime).total_seconds()/(60*60)
            if np.abs(time_diff) <= time_window:
                insitu_time[insitu_idx] = insitu_datetime_str
                insitu_filename[insitu_idx] = os.path.basename(line[:-1])
                insitu_filepath[insitu_idx] = line[:-1]
                time_difference[insitu_idx] = float(time_diff)*60*60 # in seconds

                if args.insitu == 'PANTHYR':            
                    # get data from csv using pandas
                    data = pd.read_csv(line[:-1],parse_dates=['timestamp'])   
                    
                    if insitu_idx == 0:
                        wl0 = data['wavelength'].tolist()
                        insitu_bands[:] = wl0
                    insitu_rhow[:,insitu_idx] =  data['rhow'].tolist()
                    insitu_idx += 1
                        # print(rhow0)
                elif args.insitu == 'HYPERNETS':
                    nc_fi = Dataset(line[:-1],'r')
                    if insitu_idx == 0:
                        insitu_bands[:] = nc_fi.variables['wavelength'][:].tolist()
                        ins_water_leaving_radiance = nc_fi.variables['water_leaving_radiance'][:]
                    
                    insitu_rhow[:,insitu_idx] =  nc_fi.variables['reflectance'][:].tolist()
                    insitu_idx += 1
                    nc_fi.close()
    if insitu_idx == 0:
        if os.path.exists(ofile):
            os.remove(ofile)
        if args.debug:
            print('Not in situ measurements within the time window. MDB file deleted!')
        return False
    else:
        return True
    
    nc_f0.close()
    
# #############################
#%%
def main():
    if sys.platform == 'linux':
        satellite_path_source1 = '/dst04-data1/OC/OLCI'
        satellite_path_source2 = 'trimmed_sources/'
        satellite_path_source = os.path.join(satellite_path_source1,satellite_path_source2)    
        path_out = '/home/Javier.Concha/MDB_py/ODATA'
    elif sys.platform == 'darwin':
        # satellite_path_source = os.path.join(code_home,'IDATA','SAT','S3A')
        satellite_path_source1 = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/Images/OLCI/'
        satellite_path_source2 = 'trimmed_sources/'
        satellite_path_source = os.path.join(satellite_path_source1,satellite_path_source2)  
        path_out = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/HYPERNETS/HYPERNETS_D7p2/MDB_py/ODATA'
    else:
        print('Error: host flag is not either mac or vm')
    print('Main Code!')

    if args.path_to_sat:
       satellite_path_source = args.path_to_sat 
    print(f'Path to satellite sources: {satellite_path_source}')
    
    insitu_path_source = args.path_to_ins
    print(f'Path to in situ {args.insitu} sources: {insitu_path_source}')
    
    if args.debug:
        print('Entering Debugging Mode:')

    if args.test:
        path_out = path_out + '/test/'
        
    if args.output:
        path_out = os.path.join(path_out,args.output)
        if not os.path.isdir(path_out):
            os.mkdir(path_out)
        if not os.path.isdir(os.path.join(path_out,'EXTRACTS')):
             os.mkdir(os.path.join(path_out,'EXTRACTS'))
    
    # look for in situ data within t hours
    # save nc file
    
    station_name = args.sitename
    
    
    # in situ location based on the station name
    in_situ_lat, in_situ_lon = cfs.get_lat_lon_ins(station_name)
    if args.debug:
        print(f'station_name: {station_name} with lat: {in_situ_lat}, lon: {in_situ_lon}')
        
    
    # create list of sat granules
    if args.resolution == 'WRR':
        res = 'WRR'
    else:
        res = 'WFR'
        
    wce = f'"*OL_2_{res}*SEN3"' # wild card expression
    path_to_satellite_list = create_list_products(satellite_path_source,path_out,wce,res,'satellite')
    
    if os.path.exists(f'{path_out}/OL_2_{res}_list.txt'):
        os.remove(f'{path_out}/OL_2_{res}_list.txt')

    if os.path.exists(f'{path_out}/OL_2_{res}_list.txt'):
        os.remove(f'{path_out}/OL_2_{res}_list.txt')

    if res == 'WFR':
        res_L1_str = 'EFR'
    elif res == 'WRR':
        res_L1_str = 'ERR'
        
    if os.path.exists(f'{path_out}/OL_1_{res_L1_str}_list.txt'):
        os.remove(f'{path_out}/OL_1_{res_L1_str}_list.txt')   
    
    # create list of in situ files
    if args.insitu == 'PANTHYR':
        wce = f'"*AZI_270_data.csv"' # wild card expression
    elif args.insitu == 'HYPERNETS': # 'HYPSTAR':
        wce =  f'"HYPERNETS_W_*{args.sitename}*L2A_REF*_v1.1.nc"'
    
    path_to_insitu_list = create_list_products(insitu_path_source,path_out,wce,res,'insitu')
    
    # create extract and save it in internal folder
    size_box = 25
    insitu_sensor = args.insitu
    # in situ 
    
    time_window = 3 # in hours (+- hours)
    
    if args.startdate:
        datetime_start = datetime.strptime(args.startdate, '%Y-%m-%d')
    else:
        datetime_start = datetime.strptime('2000-01-01', '%Y-%m-%d')
    if args.enddate:
        datetime_end = datetime.strptime(args.enddate, '%Y-%m-%d') + timedelta(seconds=59,minutes=59,hours=23)
    else:
        datetime_end = datetime.today()      
         
    if args.debug:
        print(f'Start date: {datetime_start}')  
        print(f'End date: {datetime_end}')
    
    with open(path_to_satellite_list,'r') as file:
            for cnt, line in enumerate(file):
                path_to_sat_source = line[:-1]
                # extract date time info
                sensor_str = path_to_sat_source.split('/')[-1].split('_')[0]
                res_str = path_to_sat_source.split('/')[-1].split('_')[3]
                datetime_str = path_to_sat_source.split('/')[-1].split('_')[7]
                if args.debug:
                    print('-----------------')
                    print(f'{datetime_str} {sensor_str} {res_str}')                
                date_format = '%Y%m%dT%H%M%S'
                satellite_datetime = datetime.strptime(datetime_str, date_format)
                datetime_creation = datetime.today().strftime(date_format)
                if satellite_datetime >= datetime_start and satellite_datetime <= datetime_end:
                    try:
                        path_to_list_daily = create_insitu_list_daily(path_to_insitu_list,datetime_str)
                        if not os.stat(path_to_list_daily).st_size == 0: # no PANTHYR data or not for that angle
                            extract_path = \
                                create_extract(size_box,station_name,path_to_sat_source,path_out,in_situ_lat,in_situ_lon,res_str,insitu_sensor)
                            
                            if args.output:
                                ofile = os.path.join(path_out,f'MDB_{sensor_str}_{res_str}_{datetime_str}_{datetime_creation}_{insitu_sensor}_{station_name}.nc')
                            else:
                                ofile = os.path.join(path_out,f'MDBs',f'MDB_{sensor_str}_{res_str}_{datetime_str}_{datetime_creation}_{insitu_sensor}_{station_name}.nc')
                
                            if add_insitu(extract_path,ofile,path_to_list_daily,datetime_str,time_window):
                                add_OL_12_to_list(path_to_sat_source,path_out,res_str)
                                print(f'file created: {ofile}')
                        else:
                            if args.debug:
                                print('No in situ measurements found!')
                    
                    # except:
                    except Exception as e:
                        if args.debug:
                            print(f'Exception: {e}')
                        pass
                    
                    if os.path.exists(path_to_list_daily):
                            os.remove(path_to_list_daily)
                else:
                    if args.debug:
                        print('Out of time frame.')
# %%
if __name__ == '__main__':
    main()        
