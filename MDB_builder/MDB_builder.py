#!/usr/bin/env python3
# coding: utf-8
"""
Created on Tue Jun 16 12:02:40 2020
Create Matchup Data Set MDB file (NetCDF4 file)
@author: javier.concha

Run as:

python MDB_builder.py -sat S3A -ins AERONET -idir IDATA/AERONET -odir ODATA/ 

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


import shutil

import zipfile
import subprocess

from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import argparse
from datetime import datetime

# import user defined functions from other .py
code_home = os.path.abspath('../')
sys.path.append(code_home)

import BRDF.brdf_olci as brdf
import COMMON.common_functions as cfs

os.environ['QT_QPA_PLATFORM']='offscreen' # to avoid error "QXcbConnection: Could not connect to display"

# user defined functions
def create_list_products(path_source,path_out,wce):
    cmd = f'find {path_source} -name {wce}> {path_out}/file_list.txt'
    prog = subprocess.Popen(cmd, shell=True,stderr=subprocess.PIPE)
    out, err = prog.communicate()
    if err:
        print(err)  
    path_to_list = f'{path_out}/file_list.txt'
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

def extract_box(size_box,station_name,path_source,path_output,in_situ_lat,in_situ_lon):
    #%
  
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
            
            if os.path.exists(ofname):
              os.remove(ofname)
            
            fmb = Dataset(ofname, 'w', format='NETCDF4')
            fmb.description = 'OLCI NxN extract'
            fmb.start_time = nc_f0.start_time
            fmb.stop_time = nc_f0.stop_time    
            fmb.input_file = path_source
            
            fmb.createDimension('satellite_id', None)
            fmb.createDimension('size_box_x', size_box)
            fmb.createDimension('size_box_y', size_box)
            
            row_center = fmb.createVariable('row_center',  'f4', ('satellite_id'), fill_value=-999, zlib=True, complevel=6) 
            row_center[:] = r
            row_center.description = 'row index to the original L2 file'
            col_center = fmb.createVariable('col_center',  'f4', ('satellite_id'), fill_value=-999, zlib=True, complevel=6) 
            col_center[:] = c
            col_center.description = 'column index to the original L2 file'
            
            lat_insitu = fmb.createVariable('lat_insitu',  'f4', ('satellite_id'), fill_value=-999, zlib=True, complevel=6) 
            lat_insitu[:] = in_situ_lat
            lat_insitu.description = 'latitude of in situ measurement'
            
            lon_insitu = fmb.createVariable('lon_insitu',  'f4', ('satellite_id'), fill_value=-999, zlib=True, complevel=6) 
            lon_insitu[:] = in_situ_lon
            lon_insitu.description = 'longitude of in situ measurement'
            
            ws0_value = fmb.createVariable('ws0_value',  'f4', ('satellite_id'), fill_value=-999, zlib=True, complevel=6) 
            ws0_value[:] = ws0
            ws1_value = fmb.createVariable('ws1_value',  'f4', ('satellite_id'), fill_value=-999, zlib=True, complevel=6) 
            ws1_value[:] = ws1
            sza_value = fmb.createVariable('sza_value',  'f4', ('satellite_id'), fill_value=-999, zlib=True, complevel=6) 
            sza_value[:] = sza
            saa_value = fmb.createVariable('saa_value',  'f4', ('satellite_id'), fill_value=-999, zlib=True, complevel=6) 
            saa_value[:] = saa
            vza_value = fmb.createVariable('vza_value',  'f4', ('satellite_id'), fill_value=-999, zlib=True, complevel=6) 
            vza_value[:] = vza
            vaa_value = fmb.createVariable('vaa_value',  'f4', ('satellite_id'), fill_value=-999, zlib=True, complevel=6) 
            vaa_value[:] = vaa
            
            latitude = fmb.createVariable('latitude',  'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6) 
            latitude[:,:] = [lat[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]]
            
            longitude = fmb.createVariable('longitude',  'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            longitude[:,:] = [lon[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]]
    
            # NOT BRDF-corrected
            rhow_0400p00_box=fmb.createVariable('rhow_0400p00', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0400p00_box[:,:] = [ma.array(rhow_0400p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0400p00_box.description = 'rhow(0400.00) NOT brdf-corrected'

            rhow_0412p50_box=fmb.createVariable('rhow_0412p50', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0412p50_box[:,:] = [ma.array(rhow_0412p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0412p50_box.description = 'rhow(0412.50) NOT brdf-corrected'
            
            rhow_0442p50_box=fmb.createVariable('rhow_0442p50', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0442p50_box[:,:] = [ma.array(rhow_0442p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0442p50_box.description = 'rhow(0442.50) NOT brdf-corrected'
            
            rhow_0490p00_box=fmb.createVariable('rhow_0490p00', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0490p00_box[:,:] = [ma.array(rhow_0490p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0490p00_box.description = 'rhow(0490.00) NOT brdf-corrected'
            
            rhow_0510p00_box=fmb.createVariable('rhow_0510p00', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0510p00_box[:,:] = [ma.array(rhow_0510p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0510p00_box.description = 'rhow(0510.00) NOT brdf-corrected'
            
            rhow_0560p00_box=fmb.createVariable('rhow_0560p00', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0560p00_box[:,:] = [ma.array(rhow_0560p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0560p00_box.description = 'rhow(0560.00) NOT brdf-corrected'
    
            rhow_0620p00_box=fmb.createVariable('rhow_0620p00', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0620p00_box[:,:] = [ma.array(rhow_0620p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0620p00_box.description = 'rhow(0620.00) NOT brdf-corrected'
            
            rhow_0665p00_box=fmb.createVariable('rhow_0665p00', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0665p00_box[:,:] = [ma.array(rhow_0665p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0665p00_box.description = 'rhow(0665.00) NOT brdf-corrected'
    
            rhow_0673p75_box=fmb.createVariable('rhow_0673p75', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0673p75_box[:,:] = [ma.array(rhow_0673p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0673p75_box.description = 'rhow(0673.75) NOT brdf-corrected'

            rhow_0681p25_box=fmb.createVariable('rhow_0681p25', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0681p25_box[:,:] = [ma.array(rhow_0681p25[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0681p25_box.description = 'rhow(0681.25) NOT brdf-corrected'

            rhow_0708p75_box=fmb.createVariable('rhow_0708p75', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0708p75_box[:,:] = [ma.array(rhow_0708p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0708p75_box.description = 'rhow(0708.75) NOT brdf-corrected'

            rhow_0753p75_box=fmb.createVariable('rhow_0753p75', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0753p75_box[:,:] = [ma.array(rhow_0753p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0753p75_box.description = 'rhow(0753.75) NOT brdf-corrected'

            rhow_0778p75_box=fmb.createVariable('rhow_0778p75', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0778p75_box[:,:] = [ma.array(rhow_0778p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0778p75_box.description = 'rhow(0778.75) NOT brdf-corrected'
            
            rhow_0865p00_box=fmb.createVariable('rhow_0865p00', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0865p00_box[:,:] = [ma.array(rhow_0865p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0865p00_box.description = 'rhow(0865.00) NOT brdf-corrected'

            rhow_0885p00_box=fmb.createVariable('rhow_0885p00', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0885p00_box[:,:] = [ma.array(rhow_0885p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_0885p00_box.description = 'rhow(0885.00) NOT brdf-corrected'
            
            rhow_1020p50_box=fmb.createVariable('rhow_1020p50', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_1020p50_box[:,:] = [ma.array(rhow_1020p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            rhow_1020p50_box.description = 'rhow(1020.50) NOT brdf-corrected'
            
            # BRDF-corrected
            rhow_0412p50_fq=fmb.createVariable('rhow_0412p50_fq', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0412p50_fq[:,:] = [ma.array(rhow_0412p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF0)]
            rhow_0412p50_fq.description = 'rhow(0412.50) brdf-corrected'
            
            rhow_0442p50_fq=fmb.createVariable('rhow_0442p50_fq', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0442p50_fq[:,:] = [ma.array(rhow_0442p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF1)]
            rhow_0442p50_fq.description = 'rhow(0442.50) brdf-corrected'
            
            rhow_0490p00_fq=fmb.createVariable('rhow_0490p00_fq', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0490p00_fq[:,:] = [ma.array(rhow_0490p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF2)]
            rhow_0490p00_fq.description = 'rhow(0490.00) brdf-corrected'
            
            rhow_0510p00_fq=fmb.createVariable('rhow_0510p00_fq', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0510p00_fq[:,:] = [ma.array(rhow_0510p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF3)]
            rhow_0510p00_fq.description = 'rhow(0510.00) brdf-corrected'
            
            rhow_0560p00_fq=fmb.createVariable('rhow_0560p00_fq', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0560p00_fq[:,:] = [ma.array(rhow_0560p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF4)]
            rhow_0560p00_fq.description = 'rhow(0560.00) brdf-corrected'
    
            rhow_0620p00_fq=fmb.createVariable('rhow_0620p00_fq', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0620p00_fq[:,:] = [ma.array(rhow_0620p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF5)]
            rhow_0620p00_fq.description = 'rhow(0620.00) brdf-corrected'
            
            rhow_0665p00_fq=fmb.createVariable('rhow_0665p00_fq', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            rhow_0665p00_fq[:,:] = [ma.array(rhow_0665p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF6)]
            rhow_0665p00_fq.description = 'rhow(0665.00) brdf-corrected'
            
            AOT_0865p50_box=fmb.createVariable('AOT_0865p50', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            AOT_0865p50_box[:,:] = [ma.array(AOT_0865p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            AOT_0865p50_box.description = 'Aerosol optical thickness'
    
            WQSF_box=fmb.createVariable('WQSF', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            WQSF_box[:,:] = [ma.array(WQSF[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            WQSF_box.description = 'OLCI Level 2 WATER Product, Classification, Quality and Science Flags Data Set'
            
            fq_0 = fmb.createVariable('fq_0', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            fq_0[:,:] = [ma.array(BRDF0)]
    
            fq_1 = fmb.createVariable('fq_1', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            fq_1[:,:] = [ma.array(BRDF1)]
    
            fq_2 = fmb.createVariable('fq_2', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            fq_2[:,:] = [ma.array(BRDF2)]
    
            fq_3 = fmb.createVariable('fq_3', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            fq_3[:,:] = [ma.array(BRDF3)]
    
            fq_4 = fmb.createVariable('fq_4', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            fq_4[:,:] = [ma.array(BRDF4)]
    
            fq_5 = fmb.createVariable('fq_5', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            fq_5[:,:] = [ma.array(BRDF5)]
    
            fq_6 = fmb.createVariable('fq_6', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            fq_6[:,:] = [ma.array(BRDF6)]
    
            chl_oc4me = fmb.createVariable('chl_oc4me', 'f4', ('satellite_id','size_box_x','size_box_y'), fill_value=-999, zlib=True, complevel=6)
            chl_oc4me[:,:] = [ma.array(CHL_OC4ME_extract)]
            
            fmb.close()
            print('Extract created!')
        else:
            print('Index out of bound!')
    else:
        print('File does NOT contains the in situ location!')
    return ofname

# #############################
# def add_insitu(ofname,satellite_id):
#     # in situ
#     path = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/HYPERNETS_D7p2/MDB_py/IDATA/INSITU/AERONET'   
# #    filename = station_name+'_20V3_20190927_20200110.nc'
#     filename = 'Venise_20_201601001_201612031.nc'
#     # filename = station_name+'_20V3_20180622_20180822.nc'
#     filename_insitu = os.path.join(path,filename)
#     if not os.path.exists(filename_insitu):
#         print('File does not exist')
        
#     nc_f0 = Dataset(filename_insitu,'r')

#     Time = nc_f0.variables['Time'][:]
#     Level = nc_f0.variables['Level'][:]
#     Julian_day = nc_f0.variables['Julian_day'][:]
#     Exact_wavelengths = nc_f0.variables['Exact_wavelengths'][:]
#     Lwn_fonQ = nc_f0.variables['Lwn_fonQ'][:]

#     nc_f0.close()

#     day_vec    = np.array([float(Time[i].replace(' ',':').split(':')[0]) for i in range(0,len(Time))])
#     month_vec  = np.array([float(Time[i].replace(' ',':').split(':')[1]) for i in range(0,len(Time))])
#     year_vec   = np.array([float(Time[i].replace(' ',':').split(':')[2]) for i in range(0,len(Time))])
#     hour_vec   = np.array([float(Time[i].replace(' ',':').split(':')[3]) for i in range(0,len(Time))])
#     minute_vec = np.array([float(Time[i].replace(' ',':').split(':')[4]) for i in range(0,len(Time))])
#     second_vec = np.array([float(Time[i].replace(' ',':').split(':')[5]) for i in range(0,len(Time))])

#     Julian_day_vec =np.array([float(Julian_day[i]) for i in range(0,len(Time))])
#     date_format = "%d:%m:%Y %H:%M:%S"
#     ins_time = np.array([datetime.strptime(Time[i], date_format) for i in range(0,len(Time))])

#     doy_vec = np.array([int(float(Julian_day[i])) for i in range(0,len(Time))])


#     # openin MDB created
#     year_str = ofname.split('/')[-3]
#     doy_str = ofname.split('/')[-2]       
#     nc_f1 = Dataset(ofname,'r')
    
#     date_format = "%Y-%m-%dT%H:%M:%S.%fZ" 
#     sat_start_time = datetime.strptime(nc_f1.start_time, date_format)
#     sat_stop_time = datetime.strptime(nc_f1.stop_time, date_format)
    
#     delta_time = 3# float in hours       
#     time_diff = ins_time - sat_stop_time
#     dt_hour = [i.total_seconds()/(60*60) for i in time_diff] # time diffence between in situ measurements and sat in hours
#     idx_min = np.argmin(np.abs(dt_hour))
#     matchup_idx_vec = np.abs(dt_hour) <= delta_time 

#     nday = sum(matchup_idx_vec)
#     print(str(nday)+' matchups per '+year_str+' '+doy_str)
#     if nday >=1:
#         print('----------------------------')
#         print('line '+str(cnt))
#         print('--Zibordi et al. 2009')
#         print(str(nday)+' matchups per '+year_str+' '+doy_str)

#         print(mu_Lwn_0412p50_fq_ins_zi.append(Lwn_fonQ[idx_min,3])) # 412,
#         print(mu_Lwn_0442p50_fq_ins_zi.append(Lwn_fonQ[idx_min,5])) # 441.8
#         print(mu_Lwn_0490p00_fq_ins_zi.append(Lwn_fonQ[idx_min,6])) # 488.5
#         if Exact_wavelengths[idx_min,13] != -999:
#             idx_560 = 13
#         elif Exact_wavelengths[idx_min,12] != -999:
#             idx_560 = 12
#         else: 
#             idx_560 = 11
#         print(mu_Lwn_0560p00_fq_ins_zi.append(Lwn_fonQ[idx_min,idx_560])) # 551,
#         print(mu_Lwn_0665p00_fq_ins_zi.append(Lwn_fonQ[idx_min,15])) # 667.9  

#%%
# def main():
print('Main Code!')


# look for in situ data within t hours
# save nc file

path_source = os.path.join(code_home,'IDATA','SAT','S3A')
path_source1 = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/Images/OLCI/'
path_source2 = 'trimmed_sources_NEW_Venice/'
path_source = os.path.join(path_source1,path_source2)
path_out = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/HYPERNETS_D7p2/MDB_py/ODATA'

station_name = 'Venise'

# in situ location based on the station name
lat_in_situ, lon_in_situ = cfs.get_lat_lon_ins(station_name)

# create list of sat granules
res = 'WRR'
wce = f'"*OL_2_{res}*trim*"' # wild card expression
path_to_list = create_list_products(path_source,path_out,wce)

# create extract and save it in internal folder
size_box = 25

with open(path_to_list,'r') as file:
        for cnt, line in enumerate(file):
            print('------------------')
            path_to_sat_source = line[:-1]
            ofname = extract_box(size_box,station_name,path_to_sat_source,path_out,lat_in_situ,lon_in_situ)


        
#%%
# if __name__ == '__main__':
#     main()        
