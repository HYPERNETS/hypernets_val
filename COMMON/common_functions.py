#!/usr/bin/env python3
# coding: utf-8
"""
Created on Thu Sep 19 12:00:01 2019
common functions
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
import  numpy as np
import os
from netCDF4 import Dataset
import datetime
#%%
def get_lat_lon_ins(station_name):
    if station_name == 'Galata_Platform': # Black Sea
        Latitude=43.044624
        Longitude=28.193190
    elif station_name == 'Gustav_Dalen_Tower': # Baltic Sea
        Latitude=58.594170
        Longitude=17.466830
    elif station_name == 'Helsinki_Lighthouse': # Baltic Sea
        Latitude=59.948970
        Longitude=24.926360
    elif station_name == 'Lake_Erie':
        Latitude=41.825600
        Longitude=-83.193600
    elif station_name == 'LISCO':
        Latitude=40.954517
        Longitude=-73.341767
    elif station_name == 'Palgrunden':
        Latitude=58.755333
        Longitude=13.151500
    elif station_name == 'Thornton_C-power':
        Latitude=51.532500
        Longitude=2.955278
    elif station_name == 'USC_SEAPRISM':
        Latitude=33.563710
        Longitude=-118.117820
    elif station_name == 'USC_SEAPRISM_2':
        Latitude=33.563710
        Longitude=-118.117820
    elif (station_name == 'Venise' or station_name == 'Venise_PANTHYR'): # Adriatic Sea
        Latitude=45.313900
        Longitude=12.508300
    elif station_name == 'WaveCIS_Site_CSI_6':
        Latitude=28.866667
        Longitude=-90.483333
    elif station_name == 'Gloria': # Black Sea
        Latitude=44.599970
        Longitude=29.359670
    else:
        print('ERROR: station not found: '+station_name)
    return Latitude, Longitude

#%    root mean squared error
def rmse(predictions, targets):
    return np.sqrt(((np.asarray(predictions) - np.asarray(targets)) ** 2).mean())    

#%% Get FO from Thuiller
def get_F0(wl,path_main):
    path_to_file = os.path.join(path_main,'Thuillier_F0.nc')
    nc_f0 = Dataset(path_to_file,'r')
    Wavelength = nc_f0.variables['Wavelength'][:]
    F0 = nc_f0.variables['F0'][:]
    F0_val = np.interp(wl,Wavelength,F0)
    return F0_val

def convert_DMS_to_decimal(DD,MM,SS,cardinal):
    D = DD + (MM/60) + (SS/3600)
    if (cardinal == 'S') or (cardinal == 'W'):
        D = -D
    return D

def contain_location(lat,lon,in_situ_lat,in_situ_lon):    
    if in_situ_lat >= lat.min()  and in_situ_lat <= lat.max() and in_situ_lon >= lon.min() and in_situ_lon <= lon.max():
        contain_flag = 1
    else:
        contain_flag = 0
        
    return contain_flag

def find_row_column_from_lat_lon(lat,lon,lat0,lon0):
    #% closest squared distance
    #% lat and lon are arrays of MxN
    #% lat0 and lon0 is the coordinates of one point
    if contain_location(lat,lon,lat0,lon0):
        dist_squared = (lat-lat0)**2 + (lon-lon0)**2
        r, c = np.unravel_index(np.argmin(dist_squared),lon.shape) # index to the closest in the latitude and longitude arrays
    else:
        print('Warning: Location not contained in the file!!!')
        r = np.nan
        c = np.nan
    return r, c

def extract_box(granule,r,c,size_box):
    offset = int(size_box/2)
    extract = granule[r-offset:r+offset+1,c-offset:c+offset+1]
    return extract
    
def doy_from_YYYYMMDD(date):
    adate = datetime.datetime.strptime(date,"%Y%m%d")
    day_of_year = adate.timetuple().tm_yday
    return day_of_year
#%%
def main():
    print('common_functions.py loaded!')
#%%
if __name__ == '__main__':
    main() 