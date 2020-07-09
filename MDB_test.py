#!/usr/bin/env python3
# coding: utf-8
"""
Created on Thu May  7 21:28:30 2020

@author: javier.concha
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
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

# User defined functions
sys.path.insert(0,'/Users/javier.concha/Desktop/Javier/2019_ROMA/CNR_Research/HYPERNETS_Validation_Protocols/python_scripts')
import common_functions

#%%
path_main = '/Users/javier.concha/Desktop/Javier/2019_ROMA/CNR_Research/HYPERNETS_D7p2/data/AERONET-OC'
filename = 'MDB_S3A_OLCI_L2_AERONET_Venise.nc'
path_file = os.path.join(path_main,filename)
nc_f0 = Dataset(path_file,'r')

# satellite_id:long_name = "Index with acquisition time (central pixel)" ;
# satellite_id:units = "Seconds since 1970-1-1" ;
satellite_id = nc_f0.dimensions['satellite_id']

# double satellite_OAA(satellite_id, rows, columns) ;
# double satellite_OZA(satellite_id, rows, columns) ;
# double satellite_SAA(satellite_id, rows, columns) ;
# double satellite_SZA(satellite_id, rows, columns) ;
# double central_latitude(satellite_id) ;
# double central_longitude(satellite_id) ;
# double satellite_latitude(satellite_id, rows, columns) ;
# double satellite_longitude(satellite_id, rows, columns) ;
# double satellite_altitude(satellite_id, rows, columns) ;
# double satellite_solar_flux(satellite_id, satellite_bands, columns) ;
# double satellite_FWHM(satellite_id, satellite_bands, columns) ;
# double satellite_lambda0(satellite_id, satellite_bands, columns) ;
# string insitu_site_name(satellite_id) ;
# double insitu_time(satellite_id, insitu_id) ;
# int64  insitu_available_measurements(satellite_id) ;
# double time_difference(satellite_id, insitu_id) ;
# double insitu_Day_of_Year_fraction(satellite_id, insitu_id) ;
# double insitu_Oa01_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa02_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa03_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa04_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa05_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa06_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa07_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa08_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa09_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa10_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa11_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa12_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa16_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa17_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa18_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa21_Rrs(satellite_id, insitu_id) ;
# double insitu_Oa01_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa02_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa03_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa04_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa05_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa06_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa07_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa08_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa09_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa10_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa11_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa12_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa16_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa17_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa18_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Oa21_Rrs_applied_shift(satellite_id, insitu_id) ;
# double insitu_Chlorophyll(satellite_id, insitu_id) ;
# double insitu_SAM(satellite_id, insitu_id) ;
# double insitu_Lw(satellite_id, insitu_original_bands, insitu_id) ;
# double insitu_Lw_Q(satellite_id, insitu_original_bands, insitu_id) ;
# double insitu_Lwn(satellite_id, insitu_original_bands, insitu_id) ;
# double insitu_Lwn_fQ(satellite_id, insitu_original_bands, insitu_id) ;
# insitu_original_bands = 29 ;
# double insitu_original_bands(insitu_original_bands) ;
#         insitu_original_bands:units = "nm" ;
#         insitu_original_bands:long_name = "Nominal AERONET bands central wavelengths" ;
# double insitu_Exact_Wavelengths(satellite_id, insitu_original_bands, insitu_id) ;
# double insitu_F0(satellite_id, insitu_original_bands, insitu_id) ;
# double insitu_Solar_Zenith_Angle(satellite_id, insitu_original_bands, insitu_id) ;
# double insitu_Solar_Azimuth_Angle(satellite_id, insitu_original_bands, insitu_id) ;

# satellites
# 0400p00,0412p50,0442p50,0490p00,0510p00,0560p00,0620p00,0665p00,0673p75,0681p25,0708p75,0753p75,0778p75,0865p00,0885p00,1020p00

satellite_PDU = nc_f0.variables['satellite_PDU'][:]
satellite_bands  = nc_f0.variables['satellite_bands'][:]
satellite_Oa01_Rrs  = nc_f0.variables['satellite_Oa01_Rrs'][:] # 0
satellite_Oa02_Rrs  = nc_f0.variables['satellite_Oa02_Rrs'][:] # 1
satellite_Oa03_Rrs  = nc_f0.variables['satellite_Oa03_Rrs'][:] # 2
satellite_Oa04_Rrs  = nc_f0.variables['satellite_Oa04_Rrs'][:] # 3
satellite_Oa05_Rrs  = nc_f0.variables['satellite_Oa05_Rrs'][:] # 4
satellite_Oa06_Rrs  = nc_f0.variables['satellite_Oa06_Rrs'][:] # 5
satellite_Oa07_Rrs  = nc_f0.variables['satellite_Oa07_Rrs'][:] # 6
satellite_Oa08_Rrs  = nc_f0.variables['satellite_Oa08_Rrs'][:] # 7
satellite_Oa09_Rrs  = nc_f0.variables['satellite_Oa09_Rrs'][:] # 8
satellite_Oa10_Rrs  = nc_f0.variables['satellite_Oa10_Rrs'][:] # 9
satellite_Oa11_Rrs  = nc_f0.variables['satellite_Oa11_Rrs'][:] # 10
satellite_Oa12_Rrs  = nc_f0.variables['satellite_Oa12_Rrs'][:] # 11
satellite_Oa16_Rrs  = nc_f0.variables['satellite_Oa16_Rrs'][:] # 12
satellite_Oa17_Rrs  = nc_f0.variables['satellite_Oa17_Rrs'][:] # 13
satellite_Oa18_Rrs  = nc_f0.variables['satellite_Oa18_Rrs'][:] # 14
satellite_Oa21_Rrs  = nc_f0.variables['satellite_Oa21_Rrs'][:] # 15

satellite_OAA = nc_f0.variables['satellite_OAA'][:]
satellite_OZA = nc_f0.variables['satellite_OZA'][:]
satellite_SAA = nc_f0.variables['satellite_SAA'][:]
satellite_SZA = nc_f0.variables['satellite_SZA'][:]

central_latitude = nc_f0.variables['central_latitude'][:]
central_longitude = nc_f0.variables['central_longitude'][:]

satellite_latitude = nc_f0.variables['satellite_latitude'][:]
satellite_longitude = nc_f0.variables['satellite_longitude'][:]
satellite_altitude = nc_f0.variables['satellite_altitude'][:]
satellite_solar_flux = nc_f0.variables['satellite_solar_flux'][:]
satellite_FWHM = nc_f0.variables['satellite_FWHM'][:]
satellite_lambda0 = nc_f0.variables['satellite_lambda0'][:]
insitu_site_name = nc_f0.variables['insitu_site_name'][:]

# in situ
insitu_time = nc_f0.variables['insitu_time'][:] # insitu_time(satellite_id, insitu_id)
insitu_available_measurements = nc_f0.variables['insitu_available_measurements'][:]
time_difference = nc_f0.variables['time_difference'][:]
insitu_Day_of_Year_fraction = nc_f0.variables['insitu_Day_of_Year_fraction'][:]
insitu_Oa01_Rrs = nc_f0.variables['insitu_Oa01_Rrs'][:] # double insitu_Oa01_Rrs(satellite_id, insitu_id) ;
insitu_Oa02_Rrs = nc_f0.variables['insitu_Oa02_Rrs'][:]
insitu_Oa03_Rrs = nc_f0.variables['insitu_Oa03_Rrs'][:]
insitu_Oa04_Rrs = nc_f0.variables['insitu_Oa04_Rrs'][:]
insitu_Oa05_Rrs = nc_f0.variables['insitu_Oa05_Rrs'][:]
insitu_Oa06_Rrs = nc_f0.variables['insitu_Oa06_Rrs'][:]
insitu_Oa07_Rrs = nc_f0.variables['insitu_Oa07_Rrs'][:]
insitu_Oa08_Rrs = nc_f0.variables['insitu_Oa08_Rrs'][:]
insitu_Oa09_Rrs = nc_f0.variables['insitu_Oa09_Rrs'][:]
insitu_Oa10_Rrs = nc_f0.variables['insitu_Oa10_Rrs'][:]
insitu_Oa11_Rrs = nc_f0.variables['insitu_Oa11_Rrs'][:]
insitu_Oa12_Rrs = nc_f0.variables['insitu_Oa12_Rrs'][:]
insitu_Oa16_Rrs = nc_f0.variables['insitu_Oa16_Rrs'][:]
insitu_Oa17_Rrs = nc_f0.variables['insitu_Oa17_Rrs'][:]
insitu_Oa18_Rrs = nc_f0.variables['insitu_Oa18_Rrs'][:]
insitu_Oa21_Rrs = nc_f0.variables['insitu_Oa21_Rrs'][:]
insitu_Oa01_Rrs_applied_shift = nc_f0.variables['insitu_Oa01_Rrs_applied_shift'][:] # double insitu_Oa01_Rrs_applied_shift(satellite_id, insitu_id) ;
insitu_Oa02_Rrs_applied_shift = nc_f0.variables['insitu_Oa02_Rrs_applied_shift'][:]
insitu_Oa03_Rrs_applied_shift = nc_f0.variables['insitu_Oa03_Rrs_applied_shift'][:]
insitu_Oa04_Rrs_applied_shift = nc_f0.variables['insitu_Oa04_Rrs_applied_shift'][:]
insitu_Oa05_Rrs_applied_shift = nc_f0.variables['insitu_Oa05_Rrs_applied_shift'][:]
insitu_Oa06_Rrs_applied_shift = nc_f0.variables['insitu_Oa06_Rrs_applied_shift'][:]
insitu_Oa07_Rrs_applied_shift = nc_f0.variables['insitu_Oa07_Rrs_applied_shift'][:]
insitu_Oa08_Rrs_applied_shift = nc_f0.variables['insitu_Oa08_Rrs_applied_shift'][:]
insitu_Oa09_Rrs_applied_shift = nc_f0.variables['insitu_Oa09_Rrs_applied_shift'][:]
insitu_Oa10_Rrs_applied_shift = nc_f0.variables['insitu_Oa10_Rrs_applied_shift'][:]
insitu_Oa11_Rrs_applied_shift = nc_f0.variables['insitu_Oa11_Rrs_applied_shift'][:]
insitu_Oa12_Rrs_applied_shift = nc_f0.variables['insitu_Oa12_Rrs_applied_shift'][:]
insitu_Oa16_Rrs_applied_shift = nc_f0.variables['insitu_Oa16_Rrs_applied_shift'][:]
insitu_Oa17_Rrs_applied_shift = nc_f0.variables['insitu_Oa17_Rrs_applied_shift'][:]
insitu_Oa18_Rrs_applied_shift = nc_f0.variables['insitu_Oa18_Rrs_applied_shift'][:]
insitu_Oa21_Rrs_applied_shift = nc_f0.variables['insitu_Oa21_Rrs_applied_shift'][:]
insitu_Chlorophyll = nc_f0.variables['insitu_Chlorophyll'][:]
insitu_SAM = nc_f0.variables['insitu_SAM'][:]
insitu_Lw = nc_f0.variables['insitu_Lw'][:]
insitu_Lw_Q = nc_f0.variables['insitu_Lw_Q'][:]
insitu_Lwn = nc_f0.variables['insitu_Lwn'][:]
insitu_Lwn_fQ = nc_f0.variables['insitu_Lwn_fQ'][:]
insitu_original_bands = nc_f0.variables['insitu_original_bands'][:]
insitu_Exact_Wavelengths = nc_f0.variables['insitu_Exact_Wavelengths'][:]
insitu_F0 = nc_f0.variables['insitu_F0'][:]
insitu_Solar_Zenith_Angle = nc_f0.variables['insitu_Solar_Zenith_Angle'][:]
insitu_Solar_Azimuth_Angle = nc_f0.variables['insitu_Solar_Azimuth_Angle'][:]

#%%


for sat_idx in range(100):
    lat = satellite_latitude[sat_idx,:,:]
    lon = satellite_longitude[sat_idx,:,:]
    lat0 = central_latitude[sat_idx]
    lon0 = central_longitude[sat_idx]
    r,c = common_functions.find_row_column_from_lat_lon(lat,lon,lat0,lon0)
    
    wl_idx = 0
    wl = satellite_bands[wl_idx]
    
    
    
    size_box = 3
    extract = common_functions.extract_box(satellite_Oa01_Rrs[sat_idx,:,:],r,c,size_box)
    
    Rrs_0400p00_extract = common_functions.extract_box(satellite_Oa01_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0412p50_extract = common_functions.extract_box(satellite_Oa02_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0442p50_extract = common_functions.extract_box(satellite_Oa03_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0490p00_extract = common_functions.extract_box(satellite_Oa04_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0510p00_extract = common_functions.extract_box(satellite_Oa05_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0560p00_extract = common_functions.extract_box(satellite_Oa06_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0620p00_extract = common_functions.extract_box(satellite_Oa07_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0665p00_extract = common_functions.extract_box(satellite_Oa08_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0673p75_extract = common_functions.extract_box(satellite_Oa09_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0681p25_extract = common_functions.extract_box(satellite_Oa10_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0708p75_extract = common_functions.extract_box(satellite_Oa11_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0753p75_extract = common_functions.extract_box(satellite_Oa12_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0778p75_extract = common_functions.extract_box(satellite_Oa16_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0865p00_extract = common_functions.extract_box(satellite_Oa17_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_0885p00_extract = common_functions.extract_box(satellite_Oa18_Rrs[sat_idx,:,:],r,c,size_box)
    Rrs_1020p00_extract = common_functions.extract_box(satellite_Oa21_Rrs[sat_idx,:,:],r,c,size_box)

    sat_spectrum = np.array([Rrs_0400p00_extract.mean()
        ,Rrs_0412p50_extract.mean()
        ,Rrs_0442p50_extract.mean()
        ,Rrs_0490p00_extract.mean()
        ,Rrs_0510p00_extract.mean()
        ,Rrs_0560p00_extract.mean()
        ,Rrs_0620p00_extract.mean()
        ,Rrs_0665p00_extract.mean()
        ,Rrs_0673p75_extract.mean()
        ,Rrs_0681p25_extract.mean()
        ,Rrs_0708p75_extract.mean()
        ,Rrs_0753p75_extract.mean()
        ,Rrs_0778p75_extract.mean()
        ,Rrs_0865p00_extract.mean()
        ,Rrs_0885p00_extract.mean()
        ,Rrs_1020p00_extract.mean()])
    
    plt.figure()
    plt.title(satellite_PDU[sat_idx][:31])
    plt.plot(satellite_bands,sat_spectrum,'-*')

    for ins_idx, line in enumerate(insitu_time[sat_idx,:]):
        ins_spectrum = np.array([insitu_Oa01_Rrs[sat_idx,ins_idx]
            ,insitu_Oa02_Rrs[sat_idx,ins_idx]
            ,insitu_Oa03_Rrs[sat_idx,ins_idx]
            ,insitu_Oa04_Rrs[sat_idx,ins_idx]
            ,insitu_Oa05_Rrs[sat_idx,ins_idx]
            ,insitu_Oa06_Rrs[sat_idx,ins_idx]
            ,insitu_Oa07_Rrs[sat_idx,ins_idx]
            ,insitu_Oa08_Rrs[sat_idx,ins_idx]
            ,insitu_Oa09_Rrs[sat_idx,ins_idx]
            ,insitu_Oa10_Rrs[sat_idx,ins_idx]
            ,insitu_Oa11_Rrs[sat_idx,ins_idx]
            ,insitu_Oa12_Rrs[sat_idx,ins_idx]
            ,insitu_Oa16_Rrs[sat_idx,ins_idx]
            ,insitu_Oa17_Rrs[sat_idx,ins_idx]
            ,insitu_Oa18_Rrs[sat_idx,ins_idx]
            ,insitu_Oa21_Rrs[sat_idx,ins_idx]])
        if not np.all(np.isnan(ins_spectrum)):
            print('---------------')
            print(f'sat_idx: {sat_idx};ins_idx: {ins_idx}')
            print(insitu_Oa01_Rrs[sat_idx,ins_idx])
            print(insitu_Oa02_Rrs[sat_idx,ins_idx])
            print(insitu_Oa03_Rrs[sat_idx,ins_idx])
            print(insitu_Oa04_Rrs[sat_idx,ins_idx])
            print(insitu_Oa05_Rrs[sat_idx,ins_idx])
            print(insitu_Oa06_Rrs[sat_idx,ins_idx])
            print(insitu_Oa07_Rrs[sat_idx,ins_idx])
            print(insitu_Oa08_Rrs[sat_idx,ins_idx])
            print(insitu_Oa09_Rrs[sat_idx,ins_idx])
            print(insitu_Oa10_Rrs[sat_idx,ins_idx])
            print(insitu_Oa11_Rrs[sat_idx,ins_idx])
            print(insitu_Oa12_Rrs[sat_idx,ins_idx])
            print(insitu_Oa16_Rrs[sat_idx,ins_idx])
            print(insitu_Oa17_Rrs[sat_idx,ins_idx])
            print(insitu_Oa18_Rrs[sat_idx,ins_idx])
            print(insitu_Oa21_Rrs[sat_idx,ins_idx])
            plt.plot(satellite_bands,ins_spectrum)
