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
import configparser

# import user defined functions from other .py
code_home = os.path.abspath('../')
sys.path.append(code_home)

import BRDF.brdf_olci as brdf
import COMMON.common_functions as cfs

os.environ['QT_QPA_PLATFORM']='offscreen' # to avoid error "QXcbConnection: Could not connect to display"
path2ncrcat = '/opt/local/bin/ncrcat'

import argparse
parser = argparse.ArgumentParser(description="Create Match-up DataBase files (MDB) files from satellite extracts and in situ L2 files.")
parser.add_argument("-d", "--debug", help="Debugging mode.",action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.",action="store_true")
parser.add_argument('-sd', "--startdate", help="The Start Date - format YYYY-MM-DD ")
parser.add_argument('-ed', "--enddate", help="The End Date - format YYYY-MM-DD ")
parser.add_argument('-site', "--sitename", help="Site name.",choices=['VEIT','BEFR','BSBE'])
parser.add_argument('-ins', "--insitu", help="Satellite sensor name.",choices=['PANTHYR', 'HYPERNETS']) # ,'HYPSTAR'])
parser.add_argument('-pi', "--path_to_ins", help="Path to in situ sources.")
parser.add_argument('-sat', "--satellite", help="Satellite sensor name.",choices=['OLCI', 'MSI'])
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-ps', "--path_to_sat", help="Path to satellite extracts.")
parser.add_argument('-o', "--output", help="Path to output")
parser.add_argument('-res', "--resolution", help="Resolution OL_2: WRR or WFR (for OLCI)")
parser.add_argument('-nl', "--nolist", help="Do not create satellite and in situ lists.",action="store_true")
parser.add_argument('-nd', "--nodelfiles", help="Do not delete temp files.",action="store_true")


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
def create_list_products(path_source,path_out,wce,res_str,type_product):
    path_to_list = f'{path_out}/file_{type_product}_{res_str}_list.txt'
    if not args.nolist:
        cmd = f'find {path_source} -name {wce}|sort|uniq> {path_to_list}'
        if args.debug:
            print(f'CMD: {cmd}')
        prog = subprocess.Popen(cmd, shell=True,stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(err)  
    
    return path_to_list

def copy_nc(ifile,ofile):
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
 
    
def add_insitu(extract_path,ofile,path_to_list_daily,datetime_str,time_window,ins_sensor):
    # print(f'Satellite time {datetime_str}')
    date_format = '%Y%m%dT%H%M%S'
    satellite_datetime = datetime.strptime(datetime_str, date_format)
      
    # to append to nc file
    new_MDB = copy_nc(extract_path,ofile)

    # add time window diff
    new_MDB.time_diff = f'{time_window*60*60}' # in seconds
    
    # create in situ dimensions
    new_MDB.createDimension('insitu_id', 30)
    new_MDB.createDimension('insitu_original_bands', 1605)
    # new_MDB.createDimension('insitu_Rrs_bands', None)
    
    # create variable 
    insitu_time=new_MDB.createVariable('insitu_time', 'f4', ('satellite_id','insitu_id',), zlib=True, complevel=6)
    insitu_time.units = "Seconds since 1970-1-1"
    insitu_time.description  = 'In situ time in ISO 8601 format (UTC).'
    
    insitu_filename=new_MDB.createVariable('insitu_filename', 'S2', ('satellite_id','insitu_id'), zlib=True, complevel=6)
    insitu_filename.description  = 'In situ filename.'
    
    insitu_filepath=new_MDB.createVariable('insitu_filepath', 'S2', ('satellite_id','insitu_id'), zlib=True, complevel=6)
    insitu_filepath.description  = 'In situ file path.'
    
    insitu_original_bands=new_MDB.createVariable('insitu_original_bands', 'f4', ('insitu_original_bands'), fill_value=-999, zlib=True, complevel=6)
    insitu_original_bands.description  = 'In situ bands in nm.'
    
    insitu_rhow=new_MDB.createVariable('insitu_rhow', 'f4', ('satellite_id','insitu_original_bands','insitu_id'), fill_value=-999, zlib=True, complevel=6)
    insitu_rhow.description  = 'In situ rhow.'

    time_difference=new_MDB.createVariable('time_difference', 'f4', ('satellite_id','insitu_id'), fill_value=-999, zlib=True, complevel=6)
    time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
    time_difference.units = "seconds"
     
    insitu_idx = 0
    # extract in situ data
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
                if args.debug:
                    print(f'debug: In situ found for: {line}')
                insitu_time[0,insitu_idx] = float(datetime.strptime(insitu_datetime_str,"%Y-%m-%dT%H:%M:%SZ").timestamp()) # Ex: 2021-02-24T11:31:00Z
                insitu_filename[0,insitu_idx] = os.path.basename(line[:-1])
                insitu_filepath[0,insitu_idx] = line[:-1]
                time_difference[0,insitu_idx] = float(time_diff)*60*60 # in seconds



                if ins_sensor == 'PANTHYR':            
                    # get data from csv using pandas
                    data = pd.read_csv(line[:-1],parse_dates=['timestamp'])   
                    
                    if insitu_idx == 0:
                        wl0 = data['wavelength'].tolist()
                        insitu_original_bands[:] = wl0
                    insitu_rhow_vec = [x for x, in data['rhow'][:]] 
                    insitu_rhow[0,:,insitu_idx] =  [ma.array(insitu_rhow_vec).transpose()]

                    insitu_idx += 1
                        # print(rhow0)
                elif ins_sensor == 'HYPERNETS':
                    nc_ins = Dataset(line[:-1],'r')
                    if insitu_idx == 0:
                        insitu_original_bands[:] = nc_ins.variables['wavelength'][:].tolist()
                        ins_water_leaving_radiance = nc_ins.variables['water_leaving_radiance'][:]

                    insitu_rhow_vec = [x for x, in nc_ins.variables['reflectance'][:]]
                    insitu_rhow[0,:,insitu_idx] =  [ma.array(insitu_rhow_vec).transpose()]
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
    
# #############################
#%%
def main():
    print('Creating MDB files!')
    if args.debug:
        print('Entering Debugging Mode:')
    # load config file

    if args.config_file:
        if os.path.isfile(args.config_file):
            options = config_reader(args.config_file)
    # else:
    #     print(args.config_file + ' does not exist. Please provide a valid config file path')
    #     sys.exit()    

    # path to satellite extract
    if args.path_to_sat:
        satellite_path_source = args.path_to_sat
    elif args.config_file:
        if options['file_path']['sat_extract_dir']:
            satellite_path_source = options['file_path']['sat_extract_dir']
    if args.verbose:
        print(f'Path to satellite sources: {satellite_path_source}')
        
    if args.satellite:
        sat_sensor = args.satellite
    elif args.config_file:
        sat_sensor = options['satellite_options']['sensor']
        sat_satellite = options['satellite_options']['satellite']
        sat_platform = options['satellite_options']['platform']
    if args.verbose:
        print(f'Satellite: {sat_satellite.upper()}')  
        print(f'Satellite sensor: {sat_sensor.upper()}')   
        print(f'Satellite platform: {sat_platform.upper()}') 
        
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
        
    # create list of sat extracts
    if args.resolution:
        res = args.resolution
    elif args.config_file:
        res = options['satellite_options']['resolution']

    
    if os.path.exists(f'{path_out}/OL_2_{res}_list.txt'):
        os.remove(f'{path_out}/OL_2_{res}_list.txt')

    if os.path.exists(f'{path_out}/OL_2_{res}_list.txt'):
        os.remove(f'{path_out}/OL_2_{res}_list.txt')

    # in situ site
    if args.sitename:
        station_name = args.sitename
    elif args.config_file:
        station_name = options['Time_and_sites_selection']['sites']
    in_situ_lat, in_situ_lon = cfs.get_lat_lon_ins(station_name) # in situ location based on the station name
    if args.verbose:
        print(f'station_name: {station_name} with lat: {in_situ_lat}, lon: {in_situ_lon}')
        
    if sat_sensor.upper() == 'OLCI':
        wce = f'"{sat_satellite}{sat_platform}*OL_2_{res}*{station_name}*"' # wild card expression
    if args.debug:
        print(f'Satellite extract Wild Card Expression: {wce}')


    path_to_satellite_list = create_list_products(satellite_path_source,path_out,wce,res,'satellite')
    if args.debug:
        print(f'Path to satellite extract list: {path_to_satellite_list}')        
    
    # create list of in situ files
    if args.insitu:
        ins_sensor = args.insitu
    elif args.config_file:
        if options['Time_and_sites_selection']['insitu_type']:
            ins_sensor = options['Time_and_sites_selection']['insitu_type']
    
    if ins_sensor == 'PANTHYR':
        wce = f'"*AZI_270_data.csv"' # wild card expression
    elif ins_sensor == 'HYPERNETS': # 'HYPSTAR':
        wce =  f'"HYPERNETS_W_*{station_name}*L2A_REF*_v1.2.nc"'
    if args.debug:
        print(f'In Situ Wild Card Expression: {wce}')
        
    if args.path_to_ins:
        insitu_path_source = args.path_to_ins
    elif args.config_file:
        if options['file_path']['ins_source_dir']:
            insitu_path_source = options['file_path']['ins_source_dir']
    if args.verbose:
        print(f'Path to in situ data: {insitu_path_source}')
    
    path_to_insitu_list = create_list_products(insitu_path_source,path_out,wce,res,'insitu')
    if args.debug:
        print(f'Path to in situ list: {path_to_insitu_list}')
    
    # in situ
    time_window = 3 # in hours (+- hours)

    file_list = [] # for the concatenation later
    
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
        datetime_end = datetime.strptime(args.enddate, '%Y-%m-%d') + timedelta(seconds=59,minutes=59,hours=23)
    elif args.config_file:
        if options['Time_and_sites_selection']['time_stop']:
            datetime_end = datetime.strptime(options['Time_and_sites_selection']['time_stop'], '%Y-%m-%d') + timedelta(seconds=59,minutes=59,hours=23)
        else:
            datetime_end = datetime.today()
    else:
        datetime_end = datetime.today()      
         
    if args.verbose:
        print(f'Start date: {datetime_start}')  
        print(f'End date: {datetime_end}')
    
    with open(path_to_satellite_list,'r') as file:
            for cnt, line in enumerate(file):
                extract_path = line[:-1]
                # extract date time info
                sensor_str = extract_path.split('/')[-1].split('_')[0]
                res_str = extract_path.split('/')[-1].split('_')[3]
                datetime_str = extract_path.split('/')[-1].split('_')[7]
                if args.verbose:
                    print('-----------------')
                    print(f'{datetime_str} {sensor_str} {res_str}')  
                if args.debug:
                    print(f'extract_path: {extract_path}')
                date_format = '%Y%m%dT%H%M%S'
                satellite_datetime = datetime.strptime(datetime_str, date_format)
                datetime_creation = datetime.today().strftime(date_format)
                if datetime_start <= satellite_datetime <= datetime_end:
                    try:
                        path_to_list_daily = create_insitu_list_daily(path_to_insitu_list,datetime_str)
                        if not os.stat(path_to_list_daily).st_size == 0: # no PANTHYR data or not for that angle
                            filename = f'MDB_{sensor_str}_{res_str}_{datetime_str}_{datetime_creation}_{ins_sensor}_{station_name}.nc'
                            ofile = os.path.join(path_out,filename)
                
                            if add_insitu(extract_path,ofile,path_to_list_daily,datetime_str,time_window,ins_sensor):
                                print(f'file created: {ofile}')
                                file_list.append(ofile) # for ncrcat later

                        else:
                            if args.verbose:
                                print('No in situ measurements found!')
                    
                    # except:
                    except Exception as e:
                        if args.debug:
                            print(f'Exception: {e}')
                        pass
                    
                    if os.path.exists(path_to_list_daily):
                            os.remove(path_to_list_daily)
                else:
                    if args.verbose:
                        print('Out of time frame.')
   
    level_prod = 'L2'

    #calling subprocess for concatanating ncdf files # # ncrcat -h MDB_S3*.nc outcat2.nc
    ncout_file = os.path.join(path_out,f'MDB_{sat_satellite}{sat_platform}_{sat_sensor.upper()}_{res_str}_{level_prod}_{ins_sensor}_{station_name}.nc')
    file_list.append(ncout_file)

    # concatenation
    cmd = [f"ncrcat -O -h"] + file_list
    cmd  = " ".join(cmd)
    if args.debug:
        print(f'CMD="{cmd}"')
    os.system(cmd)

    if not args.nodelfiles:
        [os.remove(f) for f in file_list[:-1]]
    
    print(f'Concatenated file created: {ncout_file}')
    # llll = subprocess.Popen(cmd, shell=True)
    # out, err = llll.communicate()
    # if err:
    #     print ('ncrcat process failed!')
    #     print(err)  
    # else:
    #     print (f'ncrcat process successful!\nfile created: {ncout_file}')
# %%
if __name__ == '__main__':
    main()        
