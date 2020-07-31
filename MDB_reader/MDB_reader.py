#!/usr/bin/env python3
# coding: utf-8
"""
Created on Mon Jul 20 22:49:47 2020

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
import subprocess
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from datetime import datetime
import pandas
from matplotlib import pyplot as plt
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
from scipy import stats

import configparser

# User Defined Functions
code_home = os.path.abspath('../')
sys.path.append(code_home)
import COMMON.common_functions as cfs
import COMMON.apply_OLCI_flags as apply_OLCI_flags

#%%
# create list of sat granules
def create_list_MDBs(path_to_source,path_out,wce,type_product):
    path_to_list = f'{path_out}/file_{type_product}_list.txt'
    cmd = f'find {path_to_source} -name {wce}|sort|uniq> {path_to_list}'
    prog = subprocess.Popen(cmd, shell=True,stderr=subprocess.PIPE)
    out, err = prog.communicate()
    if err:
        print(err)  
    
    return path_to_list

def plot_scatter(x,y,str1,path_out,prot_name,sensor_name,station_vec,min_val,max_val): 

    # replace nan in y (sat data)
    x = np.array(x)
    y = np.array(y)
    station_vec = np.array(station_vec)

    x = x[~np.isnan(y)] # it is assumed that only sat data could be nan
    station_vec = station_vec[~np.isnan(y)]
    y = y[~np.isnan(y)]


    rmse_val = np.nan
    mean_abs_rel_diff = np.nan
    mean_rel_diff = np.nan
    r_value = np.nan
    rmse_val_Venise = np.nan
    mean_abs_rel_diff_Venise = np.nan
    mean_rel_diff_Venise = np.nan
    r_value_Venise = np.nan
    rmse_val_Gloria = np.nan
    mean_abs_rel_diff_Gloria = np.nan
    mean_rel_diff_Gloria = np.nan
    r_value_Gloria = np.nan
    rmse_val_Galata_Platform = np.nan
    mean_abs_rel_diff_Galata_Platform = np.nan
    mean_rel_diff_Galata_Platform = np.nan
    r_value_Galata_Platform = np.nan
    rmse_val_Helsinki_Lighthouse = np.nan
    mean_abs_rel_diff_Helsinki_Lighthouse = np.nan
    mean_rel_diff_Helsinki_Lighthouse = np.nan
    r_value_Helsinki_Lighthouse = np.nan
    rmse_val_Gustav_Dalen_Tower = np.nan
    mean_abs_rel_diff_Gustav_Dalen_Tower = np.nan
    mean_rel_diff_Gustav_Dalen_Tower = np.nan
    r_value_Gustav_Dalen_Tower  = np.nan

    count_Venise = 0
    count_Gloria = 0
    count_Galata_Platform = 0
    count_Helsinki_Lighthouse = 0
    count_Gustav_Dalen_Tower = 0

    plt.figure()
    #plt.errorbar(x, y, xerr=e_x, yerr=e_y, fmt='or')
    for cnt, line in enumerate(y):
        if station_vec[cnt] == 'Venise_PAN':
            mrk_color = 'r'
            count_Venise = count_Venise+1
        elif station_vec[cnt] == 'Gloria':
            mrk_color = 'g'
            count_Gloria = count_Gloria+1
        elif station_vec[cnt] == 'Galata_Platform':
            mrk_color = 'b'
            count_Galata_Platform = count_Galata_Platform+1
        elif station_vec[cnt] == 'Helsinki_Lighthouse':
            mrk_color = 'm'
            count_Helsinki_Lighthouse = count_Helsinki_Lighthouse+1
        elif station_vec[cnt] == 'Gustav_Dalen_Tower':
            mrk_color = 'c'
            count_Gustav_Dalen_Tower = count_Gustav_Dalen_Tower+1

        if prot_name == 'ba':
            mrk_style = 'x'
        elif  prot_name == 'zi':
            mrk_style = '+' 

        plt.plot(x[cnt], y[cnt],color=mrk_color,marker=mrk_style)
    plt.axis([min_val, max_val, min_val, max_val])
    plt.gca().set_aspect('equal', adjustable='box')
    # plot 1:1 line
    xmin, xmax = plt.gca().get_xlim()
    ymin, ymax = plt.gca().get_ylim()
    plt.plot([xmin,xmax],[ymin, ymax],'--k')
    
    # Generated linear fit
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    line = slope*np.array([xmin,xmax],dtype=np.float32)+intercept
    plt.plot([xmin,xmax], line)
    # plt.legend(['1:1','Regression Line'])
    plt.xlabel(r'$\rho^{PANTHYR}_{W}$',fontsize=12)
    plt.ylabel(r'$\rho^{'+sensor_name+'}_{W}$',fontsize=12)
    if (xmin<0 or ymin<0):
        plt.plot([xmin,xmax],[0, 0],'--k',linewidth = 0.7)  
    
    # stats
    N = len(x)
    
    ref_obs = np.asarray(x)
    sat_obs = np.asarray(y)
    rmse_val = cfs.rmse(sat_obs,ref_obs)

            # the mean of relative (signed) percent differences
    rel_diff = 100*(ref_obs-sat_obs)/ref_obs
    mean_rel_diff = np.mean(rel_diff)
        
        #  the mean of absolute (unsigned) percent differences
    mean_abs_rel_diff = np.mean(np.abs(rel_diff))

    cond_station = np.asarray(station_vec)=='Venise_PAN'
    if sum(cond_station):
        ref_obs_Venise = ref_obs[cond_station]
        sat_obs_Venise = sat_obs[cond_station]
        slope_Venise, intercept_Venise, r_value_Venise, p_value_Venise, std_err_Venise = stats.linregress(ref_obs_Venise,sat_obs_Venise)
        rmse_val_Venise = cfs.rmse(sat_obs_Venise,ref_obs_Venise)
        rel_diff_Venise = 100*(ref_obs_Venise-sat_obs_Venise)/ref_obs_Venise
        mean_rel_diff_Venise = np.mean(rel_diff_Venise)
        mean_abs_rel_diff_Venise = np.mean(np.abs(rel_diff_Venise))
    
        cond_station = np.asarray(station_vec)=='Gloria'
    if sum(cond_station):    
        ref_obs_Gloria = ref_obs[cond_station]
        sat_obs_Gloria = sat_obs[cond_station]
        slope_Gloria, intercept_Gloria, r_value_Gloria, p_value_Gloria, std_err_Gloria = stats.linregress(ref_obs_Gloria,sat_obs_Gloria)
        rmse_val_Gloria = cfs.rmse(sat_obs_Gloria,ref_obs_Gloria)
        rel_diff_Gloria = 100*(ref_obs_Gloria-sat_obs_Gloria)/ref_obs_Gloria
        mean_rel_diff_Gloria = np.mean(rel_diff_Gloria)
        mean_abs_rel_diff_Gloria = np.mean(np.abs(rel_diff_Gloria))
        
        cond_station = np.asarray(station_vec)=='Galata_Platform'
    if sum(cond_station):    
        ref_obs_Galata_Platform = ref_obs[cond_station]
        sat_obs_Galata_Platform = sat_obs[cond_station]
        slope_Galata_Platform, intercept_Galata_Platform, r_value_Galata_Platform, p_value_Galata_Platform, std_err_Galata_Platform = stats.linregress(ref_obs_Galata_Platform,sat_obs_Galata_Platform)
        rmse_val_Galata_Platform = cfs.rmse(sat_obs_Galata_Platform,ref_obs_Galata_Platform)
        rel_diff_Galata_Platform = 100*(ref_obs_Galata_Platform-sat_obs_Galata_Platform)/ref_obs_Galata_Platform
        mean_rel_diff_Galata_Platform = np.mean(rel_diff_Galata_Platform)
        mean_abs_rel_diff_Galata_Platform = np.mean(np.abs(rel_diff_Galata_Platform))
        
        cond_station = np.asarray(station_vec)=='Helsinki_Lighthouse'
    if sum(cond_station):    
        ref_obs_Helsinki_Lighthouse = ref_obs[cond_station]
        sat_obs_Helsinki_Lighthouse = sat_obs[cond_station]
        slope_Helsinki_Lighthouse, intercept_Helsinki_Lighthouse, r_value_Helsinki_Lighthouse, p_value_Helsinki_Lighthouse, std_err_Helsinki_Lighthouse = stats.linregress(ref_obs_Helsinki_Lighthouse,sat_obs_Helsinki_Lighthouse)
        rmse_val_Helsinki_Lighthouse = cfs.rmse(sat_obs_Helsinki_Lighthouse,ref_obs_Helsinki_Lighthouse)
        rel_diff_Helsinki_Lighthouse = 100*(ref_obs_Helsinki_Lighthouse-sat_obs_Helsinki_Lighthouse)/ref_obs_Helsinki_Lighthouse
        mean_rel_diff_Helsinki_Lighthouse = np.mean(rel_diff_Helsinki_Lighthouse)
        mean_abs_rel_diff_Helsinki_Lighthouse = np.mean(np.abs(rel_diff_Helsinki_Lighthouse))
        
        cond_station = np.asarray(station_vec)=='Gustav_Dalen_Tower'
    if sum(cond_station):    
        ref_obs_Gustav_Dalen_Tower = ref_obs[cond_station]
        sat_obs_Gustav_Dalen_Tower = sat_obs[cond_station]
        slope_Gustav_Dalen_Tower, intercept_Gustav_Dalen_Tower, r_value_Gustav_Dalen_Tower, p_value_Gustav_Dalen_Tower, std_err_Gustav_Dalen_Tower = stats.linregress(ref_obs_Gustav_Dalen_Tower,sat_obs_Gustav_Dalen_Tower)
        rmse_val_Gustav_Dalen_Tower = cfs.rmse(sat_obs_Gustav_Dalen_Tower,ref_obs_Gustav_Dalen_Tower)
        rel_diff_Gustav_Dalen_Tower = 100*(ref_obs_Gustav_Dalen_Tower-sat_obs_Gustav_Dalen_Tower)/ref_obs_Gustav_Dalen_Tower
        mean_rel_diff_Gustav_Dalen_Tower = np.mean(rel_diff_Gustav_Dalen_Tower)
        mean_abs_rel_diff_Gustav_Dalen_Tower = np.mean(np.abs(rel_diff_Gustav_Dalen_Tower))
    

    str2 = str1
    # to print without .0
    if str1[-2:]=='.0':
        str2 = str2[:-2]
        
    
    str0 = '{}nm\nN={:d}\nrmse={:,.4f}\nMAPD={:,.0f}%\nMPD={:,.0f}%\n$r^2$={:,.2f}'\
    .format(str2,\
            N,\
            rmse_val,\
            mean_abs_rel_diff,\
            mean_rel_diff,\
            r_value**2)
        
    plt.text(0.05, 0.65, str0,horizontalalignment='left', fontsize=12,transform=plt.gca().transAxes)

    if prot_name == 'ba':
        prot_name_str = 'BW06'
    elif  prot_name == 'zi':
        prot_name_str = 'ZMB18' 

    plt.title(prot_name_str)    
    
    ofname = sensor_name+'_scatter_matchups_'+str1.replace(".","p")+'_'+prot_name+'.pdf'
    ofname = os.path.join(path_out,'source',ofname)
    
    plt.savefig(ofname, dpi=300)

    # latex table
    if str1 == '412.5':
        print('proto & nm & N & rmse & MAPD & MPD & $r^2$\n')
    str_table = '{} & {} & {:d} & {:,.4f} & {:,.1f} & {:,.1f} & {:,.2f}\\\\'\
    .format(prot_name_str,\
            str2,\
            N,\
            rmse_val,\
            mean_abs_rel_diff,\
            mean_rel_diff,\
            r_value**2)

    print(str_table)
 
    print('count_Venise: '+str(count_Venise))
    # print('count_Gloria: '+str(count_Gloria))
    # print('count_Galata_Platform: '+str(count_Galata_Platform))
    # print('count_Helsinki_Lighthouse: '+str(count_Helsinki_Lighthouse))
    # print('count_Gustav_Dalen_Tower: '+str(count_Gustav_Dalen_Tower))

    # plt.show()   
    return rmse_val, mean_abs_rel_diff, mean_rel_diff, r_value**2,\
        rmse_val_Venise, mean_abs_rel_diff_Venise, mean_rel_diff_Venise, r_value_Venise**2,\
        rmse_val_Gloria, mean_abs_rel_diff_Gloria, mean_rel_diff_Gloria, r_value_Gloria**2,\
        rmse_val_Galata_Platform, mean_abs_rel_diff_Galata_Platform, mean_rel_diff_Galata_Platform, r_value_Galata_Platform**2,\
        rmse_val_Helsinki_Lighthouse, mean_abs_rel_diff_Helsinki_Lighthouse, mean_rel_diff_Helsinki_Lighthouse, r_value_Helsinki_Lighthouse**2,\
        rmse_val_Gustav_Dalen_Tower, mean_abs_rel_diff_Gustav_Dalen_Tower, mean_rel_diff_Gustav_Dalen_Tower, r_value_Gustav_Dalen_Tower**2
#%%

# class PANTHYR_class(object):
#     def __init__(self,options):

def config_reader(FILEconfig):
    """
    Reads and checks configuration file
    Args:
        FILEconfig(str): configuration file path
    
    Return:
        Config object
    """
    Config = configparser.ConfigParser()
    Config.read(FILEconfig)
    input_path = str(Config['file_path']['input_directory'])
    output_path = str(Config['file_path']['output_directory'])
    if input_path == '' or input_path == ' ' or output_path == '' or output_path == ' ':
        print ('please provide input and output directories path')
        sys.exit()
    if Config['Time_and_sites_selection']['insitu_type'] not in (['PANTHYR']):
        print ("please select a valid in situ type: 'PANTHYR'")
        sys.exit()
    if Config['satellite_options']['platform'] not in (['A','B']):
        print ("please select a valid name for platform: A or B")
        sys.exit()
    if Config['Time_and_sites_selection']['insitu_type'] == 'AERONET' and Config['Time_and_sites_selection']['sites'] in (['',' ']):
        Config['Time_and_sites_selection']['sites'] = 'ALL'
    if Config['Filtering_options']['flags'] in [' ','']:
        Config['Filtering_options']['flags'] = 'None'
    if int(Config['satellite_options']['window_size']) > int(Config['satellite_options']['miniprods_size']) or int(Config['satellite_options']['window_size']) % 2 == 0:
        print ('windows_size must be an odd number between 1 and %d' %int(Config['satellite_options']['miniprods_size']))
        sys.exit()
    return Config
# def main():
    """business logic for when running this module as the primary one!"""
print('Main Code!')
#%%
path_main = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/HYPERNETS_D7p2/MDB_py/'

# config_file = file_config_parse.config_file
config_file = os.path.join(path_main,'MDB_reader','config_file_OLCI_PANTHYR.ini')

if os.path.isfile (config_file) == True:
    Config = config_reader(config_file)
else:
    print (config_file + ' does not exist. Please provide a valid config file path')
    sys.exit()

input_directory = Config['file_path']['input_directory']
path_main2 = 'ODATA/MDBs'
path_to_source = os.path.join(path_main,path_main2)
output_directory = Config['file_path']['output_directory']

'''
date_list.txt created as:
% cat file_list_local.txt|cut -d _ -f4|sort|uniq>date_list.txt
''' 
# PANTHYR Data
insitu_sensor = Config['Time_and_sites_selection']['insitu_type']
satellite_sensor = Config['satellite_options']['satellite'] # A and B

# create list of MDBs
type_product = 'MDB'
res = Config['satellite_options']['resolution']
wce = f'"{type_product}*{satellite_sensor}*{res}*{insitu_sensor}*.nc"' # wild card expression
path_to_list = create_list_MDBs(path_to_source,output_directory,wce,type_product)

with open(path_to_list) as file:
    for idx, line in enumerate(file):
        print(idx)
        MDBfile = line[:-1]
        print(MDBfile)


   
# #%%
# if __name__ == '__main__':
#     main()