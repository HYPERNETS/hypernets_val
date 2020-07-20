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
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from datetime import datetime
import pandas
from matplotlib import pyplot as plt
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
from scipy import stats
# to import apply_flags_OLCI.py
import sys
import subprocess
path_main = '/Users/javier.concha/Desktop/Javier/2019_ROMA/CNR_Research/HYPERNETS_Validation_Protocols/python_scripts/'
sys.path.insert(0,path_main)
import apply_flags_OLCI as OLCI_flags

# User Defined Functions
import common_functions

#%%
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
    rmse_val = common_functions.rmse(sat_obs,ref_obs)

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
        rmse_val_Venise = common_functions.rmse(sat_obs_Venise,ref_obs_Venise)
        rel_diff_Venise = 100*(ref_obs_Venise-sat_obs_Venise)/ref_obs_Venise
        mean_rel_diff_Venise = np.mean(rel_diff_Venise)
        mean_abs_rel_diff_Venise = np.mean(np.abs(rel_diff_Venise))
    
        cond_station = np.asarray(station_vec)=='Gloria'
    if sum(cond_station):    
        ref_obs_Gloria = ref_obs[cond_station]
        sat_obs_Gloria = sat_obs[cond_station]
        slope_Gloria, intercept_Gloria, r_value_Gloria, p_value_Gloria, std_err_Gloria = stats.linregress(ref_obs_Gloria,sat_obs_Gloria)
        rmse_val_Gloria = common_functions.rmse(sat_obs_Gloria,ref_obs_Gloria)
        rel_diff_Gloria = 100*(ref_obs_Gloria-sat_obs_Gloria)/ref_obs_Gloria
        mean_rel_diff_Gloria = np.mean(rel_diff_Gloria)
        mean_abs_rel_diff_Gloria = np.mean(np.abs(rel_diff_Gloria))
        
        cond_station = np.asarray(station_vec)=='Galata_Platform'
    if sum(cond_station):    
        ref_obs_Galata_Platform = ref_obs[cond_station]
        sat_obs_Galata_Platform = sat_obs[cond_station]
        slope_Galata_Platform, intercept_Galata_Platform, r_value_Galata_Platform, p_value_Galata_Platform, std_err_Galata_Platform = stats.linregress(ref_obs_Galata_Platform,sat_obs_Galata_Platform)
        rmse_val_Galata_Platform = common_functions.rmse(sat_obs_Galata_Platform,ref_obs_Galata_Platform)
        rel_diff_Galata_Platform = 100*(ref_obs_Galata_Platform-sat_obs_Galata_Platform)/ref_obs_Galata_Platform
        mean_rel_diff_Galata_Platform = np.mean(rel_diff_Galata_Platform)
        mean_abs_rel_diff_Galata_Platform = np.mean(np.abs(rel_diff_Galata_Platform))
        
        cond_station = np.asarray(station_vec)=='Helsinki_Lighthouse'
    if sum(cond_station):    
        ref_obs_Helsinki_Lighthouse = ref_obs[cond_station]
        sat_obs_Helsinki_Lighthouse = sat_obs[cond_station]
        slope_Helsinki_Lighthouse, intercept_Helsinki_Lighthouse, r_value_Helsinki_Lighthouse, p_value_Helsinki_Lighthouse, std_err_Helsinki_Lighthouse = stats.linregress(ref_obs_Helsinki_Lighthouse,sat_obs_Helsinki_Lighthouse)
        rmse_val_Helsinki_Lighthouse = common_functions.rmse(sat_obs_Helsinki_Lighthouse,ref_obs_Helsinki_Lighthouse)
        rel_diff_Helsinki_Lighthouse = 100*(ref_obs_Helsinki_Lighthouse-sat_obs_Helsinki_Lighthouse)/ref_obs_Helsinki_Lighthouse
        mean_rel_diff_Helsinki_Lighthouse = np.mean(rel_diff_Helsinki_Lighthouse)
        mean_abs_rel_diff_Helsinki_Lighthouse = np.mean(np.abs(rel_diff_Helsinki_Lighthouse))
        
        cond_station = np.asarray(station_vec)=='Gustav_Dalen_Tower'
    if sum(cond_station):    
        ref_obs_Gustav_Dalen_Tower = ref_obs[cond_station]
        sat_obs_Gustav_Dalen_Tower = sat_obs[cond_station]
        slope_Gustav_Dalen_Tower, intercept_Gustav_Dalen_Tower, r_value_Gustav_Dalen_Tower, p_value_Gustav_Dalen_Tower, std_err_Gustav_Dalen_Tower = stats.linregress(ref_obs_Gustav_Dalen_Tower,sat_obs_Gustav_Dalen_Tower)
        rmse_val_Gustav_Dalen_Tower = common_functions.rmse(sat_obs_Gustav_Dalen_Tower,ref_obs_Gustav_Dalen_Tower)
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
# def main():
    """business logic for when running this module as the primary one!"""
print('Main Code!')
#%%

path_main = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/PANTHYR/AAOT/'
path_out = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/PANTHYR/Figures/'

'''
date_list.txt created as:
% cat file_list_local.txt|cut -d _ -f4|sort|uniq>date_list.txt
''' 
# PANTHYR Data
list_name = 'file_list_PANTHYR.txt'
path_data = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/PANTHYR/AAOT/data'

# create list in txt file with file starting with "S3A_OL_2_WFR____"
# cmd = 'find '+path_data+'/20* -name "*data.csv"|sort|uniq>'+path_data+'/'+list_name
cmd = 'find '+path_data+'/20* -name "*AZI_270_data.csv"|sort|uniq>'+path_data+'/'+list_name # only with AZI 270
#        print(cmd)
# New process, connected to the Python interpreter through pipes:
prog = subprocess.Popen(cmd, shell=True,stderr=subprocess.PIPE)
out, err = prog.communicate()
if not err:
    path_to_list = os.path.join(path_data,list_name)
    year_vec_PAN   = []
    month_vec_PAN  = []
    day_vec_PAN    = []
    hour_vec_PAN   = []
    minute_vec_PAN = []
    second_vec_PAN = []
    ins_time_PAN   = []
    path_vec_PAN   = []
    azi_vec_PAN    = []
    
    with open(path_to_list,'r') as file0:
        for cnt0, line0 in enumerate(file0):
            date_str = line0.split('_')[5]
            time_str = line0.split('_')[6]
            year_vec_PAN.append(float(date_str[0:4]))
            month_vec_PAN.append(float(date_str[4:6])) 
            day_vec_PAN.append(float(date_str[6:8]))   
            hour_vec_PAN.append(float(time_str[0:2]))  
            minute_vec_PAN.append(float(time_str[2:4]))
            second_vec_PAN.append(float(time_str[4:6]))
            ins_time_PAN.append(datetime(int(date_str[0:4]),int(date_str[4:6]),\
                int(date_str[6:8]),int(time_str[0:2]),int(time_str[2:4]),int(time_str[4:6])))
            path_vec_PAN.append(line0)
            azi_vec_PAN.append(float(line0.split('_')[8]))
    
    year_vec_PAN   = np.array(year_vec_PAN)
    month_vec_PAN  = np.array(month_vec_PAN)
    day_vec_PAN    = np.array(day_vec_PAN)
    hour_vec_PAN   = np.array(hour_vec_PAN)
    minute_vec_PAN = np.array(minute_vec_PAN)
    second_vec_PAN = np.array(second_vec_PAN)
    ins_time_PAN = np.array(ins_time_PAN)
    azi_vec_PAN = np.array(azi_vec_PAN)

# AERONET-OC Data
year_vec_AOC       = []   
month_vec_AOC      = []  
day_vec_AOC        = []    
hour_vec_AOC       = []   
minute_vec_AOC     = [] 
second_vec_AOC     = [] 
Julian_day_vec_AOC = [] 
ins_time_AOC       = [] 
doy_vec_AOC        = [] 

filename = 'Venise_20V3_20190927_20200206.nc'
# filename = station_name+'_20V3_20180622_20180822.nc'
path = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/HYPERNETS_Validation_Protocols/python_scripts/netcdf_file'
filename_insitu = os.path.join(path,filename)
if not os.path.exists(filename_insitu):
    print('File does not exist')
    
nc_f0 = Dataset(filename_insitu,'r')

Time = nc_f0.variables['Time'][:]
Level = nc_f0.variables['Level'][:]
Julian_day = nc_f0.variables['Julian_day'][:]
Exact_wavelengths = nc_f0.variables['Exact_wavelengths'][:]
Lwn_fonQ = nc_f0.variables['Lwn_fonQ'][:]

nc_f0.close()

year_vec_AOC   = np.array([float(Time[i].replace(' ',':').split(':')[2]) for i in range(0,len(Time))])
month_vec_AOC  = np.array([float(Time[i].replace(' ',':').split(':')[1]) for i in range(0,len(Time))])
day_vec_AOC    = np.array([float(Time[i].replace(' ',':').split(':')[0]) for i in range(0,len(Time))])
hour_vec_AOC   = np.array([float(Time[i].replace(' ',':').split(':')[3]) for i in range(0,len(Time))])
minute_vec_AOC = np.array([float(Time[i].replace(' ',':').split(':')[4]) for i in range(0,len(Time))])
second_vec_AOC = np.array([float(Time[i].replace(' ',':').split(':')[5]) for i in range(0,len(Time))])

Julian_day_vec_AOC =np.array([float(Julian_day[i]) for i in range(0,len(Time))])
date_format = "%d:%m:%Y %H:%M:%S"
ins_time_AOC = np.array([datetime.strptime(Time[i], date_format) for i in range(0,len(Time))])

doy_vec_AOC = np.array([int(float(Julian_day[i])) for i in range(0,len(Time))])

#%% matchups

# extract.nc is created by create_extract.py
# Solar spectral irradiance F0 in uW/cm^2/nm
path_to_Thuillier_FO = '/Users/javier.concha/Desktop/Javier/2019_ROMA/CNR_Research/HYPERNETS_Validation_Protocols/python_scripts/'
F0_0412p50 = common_functions.get_F0(412.5,path_to_Thuillier_FO)  
F0_0442p50 = common_functions.get_F0(442.5,path_to_Thuillier_FO)
F0_0490p00 = common_functions.get_F0(490.0,path_to_Thuillier_FO)
F0_0560p00 = common_functions.get_F0(560.0,path_to_Thuillier_FO)
F0_0665p00 = common_functions.get_F0(665.0,path_to_Thuillier_FO)  

# from AERONET-OC: Lwn in [mW/(cm^2 sr um)]

# Zibordi: initialization
matchups_PAN_rhow_0400p00_ins_zi = []
matchups_PAN_rhow_0412p50_fq_ins_zi = []
matchups_PAN_rhow_0442p50_fq_ins_zi = []
matchups_PAN_rhow_0490p00_fq_ins_zi = []
matchups_PAN_rhow_0510p00_fq_ins_zi = []
matchups_PAN_rhow_0560p00_fq_ins_zi = []
matchups_PAN_rhow_0620p00_fq_ins_zi = []
matchups_PAN_rhow_0665p00_fq_ins_zi = []
matchups_PAN_rhow_0673p75_ins_zi = []
matchups_PAN_rhow_0681p25_ins_zi = []
matchups_PAN_rhow_0708p75_ins_zi = []
matchups_PAN_rhow_0753p75_ins_zi = []
matchups_PAN_rhow_0778p75_ins_zi = []
matchups_PAN_rhow_0865p00_ins_zi = []
matchups_PAN_rhow_0885p00_ins_zi = []
matchups_PAN_rhow_1020p50_ins_zi = []


matchups_PAN_rhow_0400p00_sat_zi = []
matchups_PAN_rhow_0412p50_fq_sat_zi = []
matchups_PAN_rhow_0442p50_fq_sat_zi = []
matchups_PAN_rhow_0490p00_fq_sat_zi = []
matchups_PAN_rhow_0510p00_fq_sat_zi = []
matchups_PAN_rhow_0560p00_fq_sat_zi = []
matchups_PAN_rhow_0620p00_fq_sat_zi = []
matchups_PAN_rhow_0665p00_fq_sat_zi = []
matchups_PAN_rhow_0673p75_sat_zi = []
matchups_PAN_rhow_0681p25_sat_zi = []
matchups_PAN_rhow_0708p75_sat_zi = []
matchups_PAN_rhow_0753p75_sat_zi = []
matchups_PAN_rhow_0778p75_sat_zi = []
matchups_PAN_rhow_0865p00_sat_zi = []
matchups_PAN_rhow_0885p00_sat_zi = []
matchups_PAN_rhow_1020p50_sat_zi = []

matchups_PAN_rhow_0400p00_ins_zi_station = []
matchups_PAN_rhow_0412p50_fq_ins_zi_station = []
matchups_PAN_rhow_0442p50_fq_ins_zi_station = []
matchups_PAN_rhow_0490p00_fq_ins_zi_station = []
matchups_PAN_rhow_0510p00_fq_ins_zi_station = []
matchups_PAN_rhow_0560p00_fq_ins_zi_station = []
matchups_PAN_rhow_0620p00_fq_ins_zi_station = []
matchups_PAN_rhow_0665p00_fq_ins_zi_station = []
matchups_PAN_rhow_0673p75_ins_zi_station = []
matchups_PAN_rhow_0681p25_ins_zi_station = []
matchups_PAN_rhow_0708p75_ins_zi_station = []
matchups_PAN_rhow_0753p75_ins_zi_station = []
matchups_PAN_rhow_0778p75_ins_zi_station = []
matchups_PAN_rhow_0865p00_ins_zi_station = []
matchups_PAN_rhow_0885p00_ins_zi_station = []
matchups_PAN_rhow_1020p50_ins_zi_station = []

matchups_PAN_rhow_0400p00_sat_zi_stop_time = []
matchups_PAN_rhow_0412p50_fq_sat_zi_stop_time = []
matchups_PAN_rhow_0442p50_fq_sat_zi_stop_time = []
matchups_PAN_rhow_0490p00_fq_sat_zi_stop_time = []
matchups_PAN_rhow_0510p00_fq_sat_zi_stop_time = []
matchups_PAN_rhow_0560p00_fq_sat_zi_stop_time = []
matchups_PAN_rhow_0620p00_fq_sat_zi_stop_time = []
matchups_PAN_rhow_0665p00_fq_sat_zi_stop_time = []
matchups_PAN_rhow_0673p75_sat_zi_stop_time = []
matchups_PAN_rhow_0681p25_sat_zi_stop_time = []  
matchups_PAN_rhow_0708p75_sat_zi_stop_time = []
matchups_PAN_rhow_0753p75_sat_zi_stop_time = []
matchups_PAN_rhow_0778p75_sat_zi_stop_time = []
matchups_PAN_rhow_0865p00_sat_zi_stop_time = []
matchups_PAN_rhow_0885p00_sat_zi_stop_time = []  
matchups_PAN_rhow_1020p50_sat_zi_stop_time = []   

matchups_PAN_rhow_0400p00_ins_zi_time = []
matchups_PAN_rhow_0412p50_fq_ins_zi_time = []
matchups_PAN_rhow_0442p50_fq_ins_zi_time = []
matchups_PAN_rhow_0490p00_fq_ins_zi_time = []
matchups_PAN_rhow_0510p00_fq_ins_zi_time = []
matchups_PAN_rhow_0560p00_fq_ins_zi_time = []
matchups_PAN_rhow_0620p00_fq_ins_zi_time = []
matchups_PAN_rhow_0665p00_fq_ins_zi_time = []
matchups_PAN_rhow_0673p75_ins_zi_time = []
matchups_PAN_rhow_0681p25_ins_zi_time = [] 
matchups_PAN_rhow_0708p75_ins_zi_time = []
matchups_PAN_rhow_0753p75_ins_zi_time = []
matchups_PAN_rhow_0778p75_ins_zi_time = []
matchups_PAN_rhow_0865p00_ins_zi_time = []
matchups_PAN_rhow_0885p00_ins_zi_time = [] 
matchups_PAN_rhow_1020p50_ins_zi_time = []

# Bailey and Werdell: initialization
matchups_PAN_rhow_0400p00_ins_ba = []
matchups_PAN_rhow_0412p50_fq_ins_ba = []
matchups_PAN_rhow_0442p50_fq_ins_ba = []
matchups_PAN_rhow_0490p00_fq_ins_ba = []
matchups_PAN_rhow_0510p00_fq_ins_ba = []
matchups_PAN_rhow_0560p00_fq_ins_ba = []
matchups_PAN_rhow_0620p00_fq_ins_ba = []
matchups_PAN_rhow_0665p00_fq_ins_ba = []
matchups_PAN_rhow_0673p75_ins_ba = []
matchups_PAN_rhow_0681p25_ins_ba = []
matchups_PAN_rhow_0708p75_ins_ba = []
matchups_PAN_rhow_0753p75_ins_ba = []
matchups_PAN_rhow_0778p75_ins_ba = []
matchups_PAN_rhow_0865p00_ins_ba = []
matchups_PAN_rhow_0885p00_ins_ba = []
matchups_PAN_rhow_1020p50_ins_ba = []

matchups_PAN_rhow_0400p00_sat_ba = []
matchups_PAN_rhow_0412p50_fq_sat_ba = []
matchups_PAN_rhow_0442p50_fq_sat_ba = []
matchups_PAN_rhow_0490p00_fq_sat_ba = []
matchups_PAN_rhow_0510p00_fq_sat_ba = []
matchups_PAN_rhow_0560p00_fq_sat_ba = []
matchups_PAN_rhow_0620p00_fq_sat_ba = []
matchups_PAN_rhow_0665p00_fq_sat_ba = []
matchups_PAN_rhow_0673p75_sat_ba = []
matchups_PAN_rhow_0681p25_sat_ba = []
matchups_PAN_rhow_0708p75_sat_ba = []
matchups_PAN_rhow_0753p75_sat_ba = []
matchups_PAN_rhow_0778p75_sat_ba = []
matchups_PAN_rhow_0865p00_sat_ba = []
matchups_PAN_rhow_0885p00_sat_ba = []
matchups_PAN_rhow_1020p50_sat_ba = []

matchups_PAN_rhow_0400p00_ins_ba_station = []
matchups_PAN_rhow_0412p50_fq_ins_ba_station = []
matchups_PAN_rhow_0442p50_fq_ins_ba_station = []
matchups_PAN_rhow_0490p00_fq_ins_ba_station = []
matchups_PAN_rhow_0510p00_fq_ins_ba_station = [] 
matchups_PAN_rhow_0560p00_fq_ins_ba_station = []
matchups_PAN_rhow_0620p00_fq_ins_ba_station = []
matchups_PAN_rhow_0665p00_fq_ins_ba_station = []
matchups_PAN_rhow_0673p75_ins_ba_station = []
matchups_PAN_rhow_0681p25_ins_ba_station = [] 
matchups_PAN_rhow_0708p75_ins_ba_station = []
matchups_PAN_rhow_0753p75_ins_ba_station = []
matchups_PAN_rhow_0778p75_ins_ba_station = []
matchups_PAN_rhow_0865p00_ins_ba_station = []
matchups_PAN_rhow_0885p00_ins_ba_station = [] 
matchups_PAN_rhow_1020p50_ins_ba_station = [] 

matchups_PAN_rhow_0400p00_sat_ba_stop_time = []
matchups_PAN_rhow_0412p50_fq_sat_ba_stop_time = []
matchups_PAN_rhow_0442p50_fq_sat_ba_stop_time = []
matchups_PAN_rhow_0490p00_fq_sat_ba_stop_time = []
matchups_PAN_rhow_0510p00_fq_sat_ba_stop_time = []
matchups_PAN_rhow_0560p00_fq_sat_ba_stop_time = []
matchups_PAN_rhow_0620p00_fq_sat_ba_stop_time = []
matchups_PAN_rhow_0665p00_fq_sat_ba_stop_time = []
matchups_PAN_rhow_0673p75_sat_ba_stop_time = []
matchups_PAN_rhow_0681p25_sat_ba_stop_time = []
matchups_PAN_rhow_0708p75_sat_ba_stop_time = []
matchups_PAN_rhow_0753p75_sat_ba_stop_time = []
matchups_PAN_rhow_0778p75_sat_ba_stop_time = []
matchups_PAN_rhow_0865p00_sat_ba_stop_time = []
matchups_PAN_rhow_0885p00_sat_ba_stop_time = []
matchups_PAN_rhow_1020p50_sat_ba_stop_time = []

matchups_PAN_rhow_0400p00_ins_ba_time = []
matchups_PAN_rhow_0412p50_fq_ins_ba_time = []
matchups_PAN_rhow_0442p50_fq_ins_ba_time = []
matchups_PAN_rhow_0490p00_fq_ins_ba_time = []
matchups_PAN_rhow_0510p00_fq_ins_ba_time = []   
matchups_PAN_rhow_0560p00_fq_ins_ba_time = []
matchups_PAN_rhow_0620p00_fq_ins_ba_time = []
matchups_PAN_rhow_0665p00_fq_ins_ba_time = []
matchups_PAN_rhow_0673p75_ins_ba_time = []
matchups_PAN_rhow_0681p25_ins_ba_time = []
matchups_PAN_rhow_0708p75_ins_ba_time = []
matchups_PAN_rhow_0753p75_ins_ba_time = []
matchups_PAN_rhow_0778p75_ins_ba_time = []
matchups_PAN_rhow_0865p00_ins_ba_time = []
matchups_PAN_rhow_0885p00_ins_ba_time = []
matchups_PAN_rhow_1020p50_ins_ba_time = []   

station_name = 'Venise_PAN'
# path_to_OLCI = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/PANTHYR/OLCI'
path_to_OLCI ='/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/HYPERNETS_Validation_Protocols/python_scripts/data/output'

# create list in txt file with file starting with "S3A_OL_2_WFR____"
list_name_OLCI_PANTHYR = 'file_list_Venise_PANTHYR.txt'
cmd = 'find '+path_to_OLCI+'/20* -name "*Venise_PANTHYR.nc"|sort|uniq>'+path_to_OLCI+'/'+list_name_OLCI_PANTHYR
#        print(cmd)
# New process, connected to the Python interpreter through pipes:
prog = subprocess.Popen(cmd, shell=True,stderr=subprocess.PIPE)
out, err = prog.communicate()
if err:
    print('ERROR: in '+cmd)

path_to_list = os.path.join(path_to_OLCI,list_name_OLCI_PANTHYR)
with open(path_to_list,'r') as file:
    for cnt, line in enumerate(file): 

        # print('----------------------------')
        # print('line '+str(cnt))
        year_str = line.split('/')[-3]
        doy_str = line.split('/')[-2]  
        path_to_extract = os.path.join(path_to_OLCI,line[:-1])     
        nc_f1 = Dataset(path_to_extract,'r')
        
        date_format = "%Y-%m-%dT%H:%M:%S.%fZ" 
        sat_start_time = datetime.strptime(nc_f1.start_time, date_format)
        sat_stop_time = datetime.strptime(nc_f1.stop_time, date_format)
        
        # from sat
        rhow_0400p00 = nc_f1.variables['rhow_0400p00'][:]
        rhow_0412p50_fq = nc_f1.variables['rhow_0412p50_fq'][:]
        rhow_0442p50_fq = nc_f1.variables['rhow_0442p50_fq'][:]
        rhow_0490p00_fq = nc_f1.variables['rhow_0490p00_fq'][:]
        rhow_0510p00_fq = nc_f1.variables['rhow_0510p00_fq'][:]
        rhow_0560p00_fq = nc_f1.variables['rhow_0560p00_fq'][:]
        rhow_0620p00_fq = nc_f1.variables['rhow_0620p00_fq'][:]
        rhow_0665p00_fq = nc_f1.variables['rhow_0665p00_fq'][:]
        rhow_0673p75 = nc_f1.variables['rhow_0673p75'][:]
        rhow_0681p25 = nc_f1.variables['rhow_0681p25'][:]
        rhow_0708p75 = nc_f1.variables['rhow_0708p75'][:]
        rhow_0753p75 = nc_f1.variables['rhow_0753p75'][:]
        rhow_0778p75 = nc_f1.variables['rhow_0778p75'][:]
        rhow_0865p00 = nc_f1.variables['rhow_0865p00'][:]
        rhow_0885p00 = nc_f1.variables['rhow_0885p00'][:]
        rhow_1020p50 = nc_f1.variables['rhow_1020p50'][:]
                
        WQSF = nc_f1.variables['WQSF'][:]
        AOT_0865p50 = nc_f1.variables['AOT_0865p50'][:]
        sza = nc_f1.variables['sza_value'][:]
        vza = nc_f1.variables['vza_value'][:]
        
#############################################
        # Zibordi et al. 2018
        delta_time = 2# float in hours       
        time_diff = ins_time_PAN - sat_stop_time
        dt_hour = [i.total_seconds()/(60*60) for i in time_diff] # time diffence between in situ measurements and sat in hours
        idx_min_PAN = np.argmin(np.abs(dt_hour))
        idxs_min_PAN = np.where(np.abs(dt_hour) == np.abs(dt_hour).min())
        idxs_allday_PAN = np.where((year_vec_PAN == sat_stop_time.year) & (month_vec_PAN == sat_stop_time.month) & (day_vec_PAN == sat_stop_time.day))
        matchup_idx_vec_PAN = np.abs(dt_hour) <= delta_time
        nday_PAN_zi = sum(matchup_idx_vec_PAN) 
        
        #########################################
        # plot in situ data
        # plot PANTHYR data
        plt.figure()
        scatter_legend = []
        # plot all PANTHYR data for the day
        for idx in range(len(idxs_allday_PAN[0])):
            # extract in situ data
            name_base = path_vec_PAN[idxs_allday_PAN[0][idx]][:-9]
            data = pandas.read_csv(os.path.join(path_data,name_base+'data.csv'),parse_dates=['timestamp']) 
            meta = pandas.read_csv(os.path.join(path_data,name_base+'meta.csv'))      
            wl0 = data['wavelength']
            rhow0 = data['rhow']

            if azi_vec_PAN[idxs_allday_PAN[0][idx]] == 135:
                linestyle = ':'
            elif azi_vec_PAN[idxs_allday_PAN[0][idx]] == 225:
                linestyle = '--' 
            elif azi_vec_PAN[idxs_allday_PAN[0][idx]] == 270:
                linestyle = '-' 

            plt.plot(wl0,rhow0,color='lightgray',linestyle = linestyle,label='_nolegend_') 

        # plot closest PANTHYR data
        for idx in range(len(idxs_min_PAN[0])):
            # extract in situ data
            name_base = path_vec_PAN[idxs_min_PAN[0][idx]][:-9]
            data = pandas.read_csv(os.path.join(path_data,name_base+'data.csv'),parse_dates=['timestamp']) 
            meta = pandas.read_csv(os.path.join(path_data,name_base+'meta.csv'))      
            wl = data['wavelength']
            rhow = data['rhow']

            if azi_vec_PAN[idxs_min_PAN[0][idx]] == 135:
                linestyle = ':'
            elif azi_vec_PAN[idxs_min_PAN[0][idx]] == 225:
                linestyle = '--' 
            elif azi_vec_PAN[idxs_min_PAN[0][idx]] == 270:
                linestyle = '-' 

            plt.plot(wl,rhow,'b',linestyle = linestyle)
            scatter_legend.append('PANTHYR '+str(ins_time_PAN[idxs_min_PAN[0][idx]])[11:-3]+' ('+str(int(azi_vec_PAN[idxs_min_PAN[0][idx]]))+')')
            # plt.plot([wl_0412p50_PAN,wl_0442p50_PAN,wl_0490p00_PAN,\
            #     wl_0560p00_PAN,wl_0665p00_PAN],\
            #     [rhow_0412p50_PAN,rhow_0442p50_PAN,rhow_0490p00_PAN,rhow_0560p00_PAN,rhow_0665p00_PAN],'ob')

        
        
        # Matchup with AERONET-OC. Warning: Only when there is matchup with PANTHYR
        time_diff = ins_time_AOC - sat_stop_time
        delta_time = 3
        dt_hour = [i.total_seconds()/(60*60) for i in time_diff] # time diffence between in situ measurements and sat in hours
        idx_min_AOC = np.argmin(np.abs(dt_hour))
        matchup_idx_vec_AOC = np.abs(dt_hour) <= delta_time 
        nday_AOC = sum(matchup_idx_vec_AOC) 
        if nday_AOC >=1:

            if Lwn_fonQ[idx_min_AOC,0] != -999:
                rhow_AOC_0340p00 = Lwn_fonQ[idx_min_AOC,0]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,0],path_to_Thuillier_FO)
            else:
                rhow_AOC_0340p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,1] != -999:
                rhow_AOC_0380p00 = Lwn_fonQ[idx_min_AOC,1]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,1],path_to_Thuillier_FO)
            else:
                rhow_AOC_0380p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,2] != -999:
                rhow_AOC_0400p00 = Lwn_fonQ[idx_min_AOC,2]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,2],path_to_Thuillier_FO)
            else:
                rhow_AOC_0400p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,3] != -999:
                rhow_AOC_0412p00 = Lwn_fonQ[idx_min_AOC,3]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,3],path_to_Thuillier_FO)
            else:
                rhow_AOC_0412p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,4] != -999:
                rhow_AOC_0440p00 = Lwn_fonQ[idx_min_AOC,4]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,4],path_to_Thuillier_FO)
            else:
                rhow_AOC_0440p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,5] != -999:
                rhow_AOC_0443p00 = Lwn_fonQ[idx_min_AOC,5]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,5],path_to_Thuillier_FO)
            else:
                rhow_AOC_0443p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,6] != -999:
                rhow_AOC_0490p00 = Lwn_fonQ[idx_min_AOC,6]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,6],path_to_Thuillier_FO)
            else:
                rhow_AOC_0490p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,7] != -999:
                rhow_AOC_0500p00 = Lwn_fonQ[idx_min_AOC,7]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,7],path_to_Thuillier_FO)
            else:
                rhow_AOC_0500p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,8] != -999:
                rhow_AOC_0510p00 = Lwn_fonQ[idx_min_AOC,8]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,8],path_to_Thuillier_FO)
            else:
                rhow_AOC_0510p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,9] != -999:
                rhow_AOC_0531p00 = Lwn_fonQ[idx_min_AOC,9]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,9],path_to_Thuillier_FO)
            else:
                rhow_AOC_0531p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,10] != -999:
                rhow_AOC_0532p00 = Lwn_fonQ[idx_min_AOC,10]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,10],path_to_Thuillier_FO)
            else:
                rhow_AOC_0532p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,11] != -999:
                rhow_AOC_0551p00 = Lwn_fonQ[idx_min_AOC,11]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,11],path_to_Thuillier_FO)
            else:
                rhow_AOC_0551p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,12] != -999:
                rhow_AOC_0555p00 = Lwn_fonQ[idx_min_AOC,12]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,12],path_to_Thuillier_FO)
            else:
                rhow_AOC_0555p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,13] != -999:
                rhow_AOC_0560p00 = Lwn_fonQ[idx_min_AOC,13]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,13],path_to_Thuillier_FO)
            else:
                rhow_AOC_0560p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,14] != -999:
                rhow_AOC_0620p00 = Lwn_fonQ[idx_min_AOC,14]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,14],path_to_Thuillier_FO)
            else:
                rhow_AOC_0620p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,15] != -999:
                rhow_AOC_0667p00 = Lwn_fonQ[idx_min_AOC,15]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,15],path_to_Thuillier_FO)
            else:
                rhow_AOC_0667p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,16] != -999:
                rhow_AOC_0675p00 = Lwn_fonQ[idx_min_AOC,16]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,16],path_to_Thuillier_FO)
            else:
                rhow_AOC_0675p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,17] != -999:
                rhow_AOC_0681p00 = Lwn_fonQ[idx_min_AOC,17]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,17],path_to_Thuillier_FO)
            else:
                rhow_AOC_0681p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,18] != -999:
                rhow_AOC_0709p00 = Lwn_fonQ[idx_min_AOC,18]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,18],path_to_Thuillier_FO)
            else:
                rhow_AOC_0709p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,19] != -999:
                rhow_AOC_0779p00 = Lwn_fonQ[idx_min_AOC,19]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,19],path_to_Thuillier_FO)
            else:
                rhow_AOC_0779p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,20] != -999:
                rhow_AOC_0865p00 = Lwn_fonQ[idx_min_AOC,20]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,20],path_to_Thuillier_FO)
            else:
                rhow_AOC_0865p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,21] != -999:
                rhow_AOC_0870p00 = Lwn_fonQ[idx_min_AOC,21]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,21],path_to_Thuillier_FO)
            else:
                rhow_AOC_0870p00 = np.nan
            if Lwn_fonQ[idx_min_AOC,22] != -999:
                rhow_AOC_1020p00 = Lwn_fonQ[idx_min_AOC,22]*np.pi/common_functions.get_F0(Exact_wavelengths[idx_min_AOC,22],path_to_Thuillier_FO)
            else:
                rhow_AOC_1020p00 = np.nan
            
            plt.plot([\
                Exact_wavelengths[idx_min_AOC,0],\
                Exact_wavelengths[idx_min_AOC,1],\
                Exact_wavelengths[idx_min_AOC,2],\
                Exact_wavelengths[idx_min_AOC,3],\
                Exact_wavelengths[idx_min_AOC,4],\
                Exact_wavelengths[idx_min_AOC,5],\
                Exact_wavelengths[idx_min_AOC,6],\
                Exact_wavelengths[idx_min_AOC,7],\
                Exact_wavelengths[idx_min_AOC,8],\
                Exact_wavelengths[idx_min_AOC,9],\
                Exact_wavelengths[idx_min_AOC,10],\
                Exact_wavelengths[idx_min_AOC,11],\
                Exact_wavelengths[idx_min_AOC,12],\
                Exact_wavelengths[idx_min_AOC,13],\
                Exact_wavelengths[idx_min_AOC,14],\
                Exact_wavelengths[idx_min_AOC,15],\
                Exact_wavelengths[idx_min_AOC,16],\
                Exact_wavelengths[idx_min_AOC,17],\
                Exact_wavelengths[idx_min_AOC,18],\
                Exact_wavelengths[idx_min_AOC,19],\
                Exact_wavelengths[idx_min_AOC,20],\
                Exact_wavelengths[idx_min_AOC,21],\
                Exact_wavelengths[idx_min_AOC,22]],[\
                rhow_AOC_0340p00,\
                rhow_AOC_0380p00,\
                rhow_AOC_0400p00,\
                rhow_AOC_0412p00,\
                rhow_AOC_0440p00,\
                rhow_AOC_0443p00,\
                rhow_AOC_0490p00,\
                rhow_AOC_0500p00,\
                rhow_AOC_0510p00,\
                rhow_AOC_0531p00,\
                rhow_AOC_0532p00,\
                rhow_AOC_0551p00,\
                rhow_AOC_0555p00,\
                rhow_AOC_0560p00,\
                rhow_AOC_0620p00,\
                rhow_AOC_0667p00,\
                rhow_AOC_0675p00,\
                rhow_AOC_0681p00,\
                rhow_AOC_0709p00,\
                rhow_AOC_0779p00,\
                rhow_AOC_0865p00,\
                rhow_AOC_0870p00,\
                rhow_AOC_1020p00],'ob',mfc='none')
            scatter_legend.append('AERONET-OC '+str(ins_time_AOC[idx_min_AOC])[11:-3])  
        #########################################
        # END plot in situ data
#############################################
        # CONTINUE Zibordi et al. 2018
        rhow_0400p00_PAN = rhow[np.argmin(np.abs(wl-400.00))] 
        rhow_0412p50_PAN = rhow[np.argmin(np.abs(wl-412.50))] 
        rhow_0442p50_PAN = rhow[np.argmin(np.abs(wl-442.50))] 
        rhow_0490p00_PAN = rhow[np.argmin(np.abs(wl-490.00))] 
        rhow_0510p00_PAN = rhow[np.argmin(np.abs(wl-510.00))] 
        rhow_0560p00_PAN = rhow[np.argmin(np.abs(wl-560.00))] 
        rhow_0620p00_PAN = rhow[np.argmin(np.abs(wl-620.00))] 
        rhow_0665p00_PAN = rhow[np.argmin(np.abs(wl-665.00))] 
        rhow_0673p75_PAN = rhow[np.argmin(np.abs(wl-673.75))] 
        rhow_0681p25_PAN = rhow[np.argmin(np.abs(wl-681.25))]
        rhow_0708p75_PAN = rhow[np.argmin(np.abs(wl-708.75))] 
        rhow_0753p75_PAN = rhow[np.argmin(np.abs(wl-753.75))] 
        rhow_0778p75_PAN = rhow[np.argmin(np.abs(wl-778.75))] 
        rhow_0865p00_PAN = rhow[np.argmin(np.abs(wl-865.00))] 
        rhow_0885p00_PAN = rhow[np.argmin(np.abs(wl-885.00))]

        wl_0400p00_PAN = wl[np.argmin(np.abs(wl-400.00))] 
        wl_0412p50_PAN = wl[np.argmin(np.abs(wl-412.50))] 
        wl_0442p50_PAN = wl[np.argmin(np.abs(wl-442.50))] 
        wl_0490p00_PAN = wl[np.argmin(np.abs(wl-490.00))] 
        wl_0510p00_PAN = wl[np.argmin(np.abs(wl-510.00))]
        wl_0560p00_PAN = wl[np.argmin(np.abs(wl-560.00))] 
        wl_0620p00_PAN = wl[np.argmin(np.abs(wl-620.00))] 
        wl_0665p00_PAN = wl[np.argmin(np.abs(wl-665.00))] 
        wl_0673p75_PAN = wl[np.argmin(np.abs(wl-673.75))] 
        wl_0681p25_PAN = wl[np.argmin(np.abs(wl-681.25))]  
        wl_0708p75_PAN = wl[np.argmin(np.abs(wl-708.75))] 
        wl_0753p75_PAN = wl[np.argmin(np.abs(wl-753.75))] 
        wl_0778p75_PAN = wl[np.argmin(np.abs(wl-778.75))] 
        wl_0865p00_PAN = wl[np.argmin(np.abs(wl-865.00))] 
        wl_0885p00_PAN = wl[np.argmin(np.abs(wl-885.00))] 

        if nday_PAN_zi >=1:
            print('----------------------------')
            print('line '+str(cnt))
            print('--Zibordi et al. 2018')
            print(str(nday_PAN_zi)+' matchups per '+year_str+' '+doy_str)
            print(sat_stop_time)
            print(data['timestamp'][0])
#            print(rhow[idx_min_PAN,:])
#            print(Exact_wavelengths[idx_min_PAN,:])
            
            center_px = int(len(rhow_0412p50_fq)/2 + 0.5)
            size_box = 3
            start_idx_x = int(center_px-int(size_box/2))
            stop_idx_x = int(center_px+int(size_box/2)+1)
            start_idx_y = int(center_px-int(size_box/2))
            stop_idx_y = int(center_px+int(size_box/2)+1)
            rhow_0400p00_box = rhow_0400p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0412p50_fq_box = rhow_0412p50_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0442p50_fq_box = rhow_0442p50_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0490p00_fq_box = rhow_0490p00_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0510p00_fq_box = rhow_0510p00_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0560p00_fq_box = rhow_0560p00_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0620p00_fq_box = rhow_0620p00_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0665p00_fq_box = rhow_0665p00_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0673p75_box = rhow_0673p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0681p25_box = rhow_0681p25[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0708p75_box = rhow_0708p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0753p75_box = rhow_0753p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0778p75_box = rhow_0778p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0865p00_box = rhow_0865p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0885p00_box = rhow_0885p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_1020p50_box = rhow_1020p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]

            AOT_0865p50_box = AOT_0865p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]

            # plot sat data
                # numpy.all(a, axis=None, out=None, keepdims=<no value>)[source]
                # Test whether all array elements along a given axis evaluate to True.
                # When all .mask are True, i.e. all invalid, .all() is True. Or if at least one is valid -> .all() is False
                # Therefore, if at least one band extract has a valid number, then plot it. 
                # This is plot the extract before applying filtering criteria.
            if not rhow_0400p00_box.mask.all() or not rhow_0412p50_fq_box.mask.all()\
                or not rhow_0442p50_fq_box.mask.all() or not rhow_0490p00_fq_box.mask.all()\
                or not rhow_0510p00_fq_box.mask.all() or not rhow_0560p00_fq_box.mask.all()\
                or not rhow_0620p00_fq_box.mask.all() or not rhow_0665p00_fq_box.mask.all()\
                or not rhow_0673p75_box.mask.all() or not rhow_0681p25_box.mask.all()\
                or not rhow_0708p75_box.mask.all() or not rhow_0753p75_box.mask.all()\
                or not rhow_0778p75_box.mask.all() or not rhow_0865p00_box.mask.all()\
                or not rhow_0885p00_box.mask.all() or not rhow_1020p50_box.mask.all():
                # plot Zibordi
                plt.plot([400.00,412.50,442.50,490.00,510.00,560.00,620.00,665.00,673.75,681.25,708.75,753.75,778.75,865.00,885.00,1020.5],\
                    [rhow_0400p00_box.mean(),rhow_0412p50_fq_box.mean(),rhow_0442p50_fq_box.mean(),rhow_0490p00_fq_box.mean(),\
                    rhow_0510p00_fq_box.mean(),rhow_0560p00_fq_box.mean(),rhow_0620p00_fq_box.mean(),rhow_0665p00_fq_box.mean(),\
                    rhow_0673p75_box.mean(),rhow_0681p25_box.mean(),rhow_0708p75_box.mean(),rhow_0753p75_box.mean(),\
                    rhow_0778p75_box.mean(),rhow_0865p00_box.mean(),rhow_0885p00_box.mean(),rhow_1020p50_box.mean()],'+r',markersize=8,linewidth=3)
                scatter_legend.append('OLCI: ZMB18') 
            
            flags_mask = OLCI_flags.create_mask(WQSF[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])
                # Ex. flags_mask when all valid:
                # [[0 0 0 0 0]
                #  [0 0 0 0 0]
                #  [0 0 0 0 0]
                #  [0 0 0 0 0]
                #  [0 0 0 0 0]]
            print('flags_mask:')
            print(flags_mask)

            # from AERONET-OC V3 file
            # 0         1         2         3         4         5         6         7         8         9         10        11        12        13        14        15        16        17        18        19        20        21        22 
            # Lw[340nm],Lw[380nm],Lw[400nm],Lw[412nm],Lw[440nm],Lw[443nm],Lw[490nm],Lw[500nm],Lw[510nm],Lw[531nm],Lw[532nm],Lw[551nm],Lw[555nm],Lw[560nm],Lw[620nm],Lw[667nm],Lw[675nm],Lw[681nm],Lw[709nm],Lw[779nm],Lw[865nm],Lw[870nm],Lw[1020nm]    
            # -999,     -999,     -999,     412,      -999,     441.8,    488.5,    -999,     -999,     -999,     530.3,    551,      -999,     -999,     -999,     667.9,    -999,     -999,     -999,     -999,     -999,     870.8,    1020.5,
                
            if sza<=70 and vza<=56 and not flags_mask.any(): # if any of the pixels if flagged, Fails validation criteria because all have to be valid in Zibordi 2018
                Lwn_560 = rhow_0560p00_fq_box
                Lwn_560_CV = np.abs(Lwn_560.std()/Lwn_560.mean())    
                
                AOT_0865p50_CV = np.abs(AOT_0865p50_box.std()/AOT_0865p50_box.mean())

                rhow_0400p00_sat_passed_plot = np.nan
                rhow_0412p50_sat_passed_plot = np.nan
                rhow_0442p50_sat_passed_plot = np.nan
                rhow_0490p00_sat_passed_plot = np.nan
                rhow_0510p00_sat_passed_plot = np.nan
                rhow_0560p00_sat_passed_plot = np.nan
                rhow_0620p00_sat_passed_plot = np.nan
                rhow_0665p00_sat_passed_plot = np.nan
                rhow_0673p75_sat_passed_plot = np.nan
                rhow_0681p25_sat_passed_plot = np.nan
                rhow_0708p75_sat_passed_plot = np.nan
                rhow_0753p75_sat_passed_plot = np.nan
                rhow_0778p75_sat_passed_plot = np.nan
                rhow_0865p00_sat_passed_plot = np.nan
                rhow_0885p00_sat_passed_plot = np.nan                
                
                if Lwn_560_CV <= 0.2 and AOT_0865p50_CV <= 0.2:
                    # if any is invalid, do not calculated matchup
                    if not ((rhow_0400p00_box.mask.any() or np.isnan(rhow_0400p00_box).any())\
                        or (rhow_0412p50_fq_box.mask.any() or np.isnan(rhow_0412p50_fq_box).any())\
                        or (rhow_0442p50_fq_box.mask.any() or np.isnan(rhow_0442p50_fq_box).any())\
                        or (rhow_0490p00_fq_box.mask.any() or np.isnan(rhow_0490p00_fq_box).any())\
                        or (rhow_0510p00_fq_box.mask.any() or np.isnan(rhow_0510p00_fq_box).any())\
                        or (rhow_0560p00_fq_box.mask.any() or np.isnan(rhow_0560p00_fq_box).any())\
                        or (rhow_0620p00_fq_box.mask.any() or np.isnan(rhow_0620p00_fq_box).any())\
                        or (rhow_0665p00_fq_box.mask.any() or np.isnan(rhow_0665p00_fq_box).any())\
                        or (rhow_0673p75_box.mask.any() or np.isnan(rhow_0673p75_box).any())\
                        or (rhow_0681p25_box.mask.any() or np.isnan(rhow_0681p25_box).any())\
                        or (rhow_0708p75_box.mask.any() or np.isnan(rhow_0708p75_box).any())\
                        or (rhow_0753p75_box.mask.any() or np.isnan(rhow_0753p75_box).any())\
                        or (rhow_0778p75_box.mask.any() or np.isnan(rhow_0778p75_box).any())\
                        or (rhow_0865p00_box.mask.any() or np.isnan(rhow_0865p00_box).any())\
                        or (rhow_0885p00_box.mask.any() or np.isnan(rhow_0885p00_box).any())):
  
                    # Rrs 0400p00
                    # print('400.0')
                    # if not (rhow_0400p00_box.mask.any() == True or np.isnan(rhow_0400p00_box).any() == True):
                        rhow_0400p00_sat_passed_plot = rhow_0400p00_box.mean()
                    #     print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0400p00_sat_zi.append(rhow_0400p00_box.mean())
                        matchups_PAN_rhow_0400p00_ins_zi.append(rhow_0400p00_PAN) # 412,
                        matchups_PAN_rhow_0400p00_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0400p00_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0400p00_ins_zi_time.append(ins_time_PAN[idx_min_PAN])

                    # Rrs 0412p50
                    # print('412.5')
                    # if not (rhow_0412p50_fq_box.mask.any() == True or np.isnan(rhow_0412p50_fq_box).any() == True):
                        rhow_0412p50_sat_passed_plot = rhow_0412p50_fq_box.mean()
                    #     print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0412p50_fq_sat_zi.append(rhow_0412p50_fq_box.mean())
                        matchups_PAN_rhow_0412p50_fq_ins_zi.append(rhow_0412p50_PAN) # 412,
                        matchups_PAN_rhow_0412p50_fq_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0412p50_fq_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0412p50_fq_ins_zi_time.append(ins_time_PAN[idx_min_PAN])
                        
                    # Rrs 0442p50
                    # print('442.5')
                    # if not (rhow_0442p50_fq_box.mask.any() == True or np.isnan(rhow_0442p50_fq_box).any() == True):
                        rhow_0442p50_sat_passed_plot = rhow_0442p50_fq_box.mean()
                        # print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0442p50_fq_sat_zi.append(rhow_0442p50_fq_box.mean())
                        matchups_PAN_rhow_0442p50_fq_ins_zi.append(rhow_0442p50_PAN) # 441.8
                        matchups_PAN_rhow_0442p50_fq_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0442p50_fq_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0442p50_fq_ins_zi_time.append(ins_time_PAN[idx_min_PAN])
                        
                    # Rrs 0490p00
                    # print('490.0')
                    # if not (rhow_0490p00_fq_box.mask.any() == True or np.isnan(rhow_0490p00_fq_box).any() == True):
                        rhow_0490p00_sat_passed_plot = rhow_0490p00_fq_box.mean()
                        # print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0490p00_fq_sat_zi.append(rhow_0490p00_fq_box.mean())
                        matchups_PAN_rhow_0490p00_fq_ins_zi.append(rhow_0490p00_PAN) # 488.5
                        matchups_PAN_rhow_0490p00_fq_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0490p00_fq_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0490p00_fq_ins_zi_time.append(ins_time_PAN[idx_min_PAN])

                    # Rrs 0510p00
                    # print('510.0')
                    # if not (rhow_0510p00_fq_box.mask.any() == True or np.isnan(rhow_0510p00_fq_box).any() == True):
                        rhow_0510p00_sat_passed_plot = rhow_0510p00_fq_box.mean()
                        # print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0510p00_fq_sat_zi.append(rhow_0510p00_fq_box.mean())
                        matchups_PAN_rhow_0510p00_fq_ins_zi.append(rhow_0510p00_PAN) # 488.5
                        matchups_PAN_rhow_0510p00_fq_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0510p00_fq_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0510p00_fq_ins_zi_time.append(ins_time_PAN[idx_min_PAN])    
                        
                    # Rrs 0560p00
                    # print('560.0')
                    # if not (rhow_0560p00_fq_box.mask.any() == True or np.isnan(rhow_0560p00_fq_box).any() == True):
                        rhow_0560p00_sat_passed_plot = rhow_0560p00_fq_box.mean()
                        # print('At least one element in sat product is invalid!')
                        matchups_PAN_rhow_0560p00_fq_sat_zi.append(rhow_0560p00_fq_box.mean())
                        matchups_PAN_rhow_0560p00_fq_ins_zi.append(rhow_0560p00_PAN) # 551,
                        matchups_PAN_rhow_0560p00_fq_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0560p00_fq_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0560p00_fq_ins_zi_time.append(ins_time_PAN[idx_min_PAN])

                    # Rrs 0620p00
                    # print('620.0')
                    # if not (rhow_0620p00_fq_box.mask.any() == True or np.isnan(rhow_0620p00_fq_box).any() == True):
                        rhow_0620p00_sat_passed_plot = rhow_0620p00_fq_box.mean()
                        # print('At least one element in sat product is invalid!')
                        matchups_PAN_rhow_0620p00_fq_sat_zi.append(rhow_0620p00_fq_box.mean())
                        matchups_PAN_rhow_0620p00_fq_ins_zi.append(rhow_0620p00_PAN) # 551,
                        matchups_PAN_rhow_0620p00_fq_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0620p00_fq_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0620p00_fq_ins_zi_time.append(ins_time_PAN[idx_min_PAN])    
                        
                    # Rrs 0665p00
                    # print('665.0')
                    # if not (rhow_0665p00_fq_box.mask.any() == True or np.isnan(rhow_0665p00_fq_box).any() == True):
                        rhow_0665p00_sat_passed_plot = rhow_0665p00_fq_box.mean()
                        # print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0665p00_fq_sat_zi.append(rhow_0665p00_fq_box.mean())
                        matchups_PAN_rhow_0665p00_fq_ins_zi.append(rhow_0665p00_PAN) # 667.9    
                        matchups_PAN_rhow_0665p00_fq_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0665p00_fq_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0665p00_fq_ins_zi_time.append(ins_time_PAN[idx_min_PAN])

                    # Rrs 0673p75
                    # print('673.75')
                    # if not (rhow_0673p75_box.mask.any() == True or np.isnan(rhow_0673p75_box).any() == True):
                        rhow_0673p75_sat_passed_plot = rhow_0673p75_box.mean()
                        # print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0673p75_sat_zi.append(rhow_0673p75_box.mean())
                        matchups_PAN_rhow_0673p75_ins_zi.append(rhow_0673p75_PAN) # 667.9    
                        matchups_PAN_rhow_0673p75_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0673p75_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0673p75_ins_zi_time.append(ins_time_PAN[idx_min_PAN])  

                    # Rrs 0681p25
                    # print('681.25')
                    # if not (rhow_0681p25_box.mask.any() == True or np.isnan(rhow_0681p25_box).any() == True):
                        rhow_0681p25_sat_passed_plot = rhow_0681p25_box.mean()
                        # print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0681p25_sat_zi.append(rhow_0681p25_box.mean())
                        matchups_PAN_rhow_0681p25_ins_zi.append(rhow_0681p25_PAN) # 667.9    
                        matchups_PAN_rhow_0681p25_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0681p25_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0681p25_ins_zi_time.append(ins_time_PAN[idx_min_PAN]) 

                    # Rrs 0708p75
                    # print('708.75')
                    # if not (rhow_0708p75_box.mask.any() == True or np.isnan(rhow_0708p75_box).any() == True):
                        rhow_0708p75_sat_passed_plot = rhow_0708p75_box.mean()
                        # print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0708p75_sat_zi.append(rhow_0708p75_box.mean())
                        matchups_PAN_rhow_0708p75_ins_zi.append(rhow_0708p75_PAN) # 667.9    
                        matchups_PAN_rhow_0708p75_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0708p75_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0708p75_ins_zi_time.append(ins_time_PAN[idx_min_PAN])     

                    # Rrs 0753p75
                    # print('753.75')
                    # if not (rhow_0753p75_box.mask.any() == True or np.isnan(rhow_0753p75_box).any() == True):
                        rhow_0753p75_sat_passed_plot = rhow_0753p75_box.mean()
                        # print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0753p75_sat_zi.append(rhow_0753p75_box.mean())
                        matchups_PAN_rhow_0753p75_ins_zi.append(rhow_0753p75_PAN) # 667.9    
                        matchups_PAN_rhow_0753p75_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0753p75_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0753p75_ins_zi_time.append(ins_time_PAN[idx_min_PAN])  

                    # Rrs 0778p75
                    # print('778.75')
                    # if not (rhow_0778p75_box.mask.any() == True or np.isnan(rhow_0778p75_box).any() == True):
                        rhow_0778p75_sat_passed_plot = rhow_0778p75_box.mean()
                        # print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0778p75_sat_zi.append(rhow_0778p75_box.mean())
                        matchups_PAN_rhow_0778p75_ins_zi.append(rhow_0778p75_PAN) # 667.9    
                        matchups_PAN_rhow_0778p75_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0778p75_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0778p75_ins_zi_time.append(ins_time_PAN[idx_min_PAN])

                    # Rrs 0865p00
                    # print('865.0')
                    # if not (rhow_0865p00_box.mask.any() == True or np.isnan(rhow_0865p00_box).any() == True):
                        rhow_0865p00_sat_passed_plot = rhow_0865p00_box.mean()
                        # print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0865p00_sat_zi.append(rhow_0865p00_box.mean())
                        matchups_PAN_rhow_0865p00_ins_zi.append(rhow_0865p00_PAN) # 667.9    
                        matchups_PAN_rhow_0865p00_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0865p00_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0865p00_ins_zi_time.append(ins_time_PAN[idx_min_PAN]) 

                    # Rrs 0885p00
                    # print('885.0')
                    # if not (rhow_0885p00_box.mask.any() == True or np.isnan(rhow_0885p00_box).any() == True):
                        rhow_0885p00_sat_passed_plot = rhow_0885p00_box.mean()
                        # print('At least one element in sat product is invalid!')
                    # else:
                        matchups_PAN_rhow_0885p00_sat_zi.append(rhow_0885p00_box.mean())
                        matchups_PAN_rhow_0885p00_ins_zi.append(rhow_0885p00_PAN) # 667.9    
                        matchups_PAN_rhow_0885p00_ins_zi_station.append(station_name)
                        matchups_PAN_rhow_0885p00_sat_zi_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0885p00_ins_zi_time.append(ins_time_PAN[idx_min_PAN])   

                        plt.plot([400.00,412.50,442.50,490.00,510.00,560.00,620.00,665.00,673.75,681.25,708.75,753.75,778.75,865.00,885.00],\
                            [rhow_0400p00_sat_passed_plot,rhow_0412p50_sat_passed_plot,rhow_0442p50_sat_passed_plot,\
                            rhow_0490p00_sat_passed_plot,rhow_0510p00_sat_passed_plot,rhow_0560p00_sat_passed_plot,\
                            rhow_0620p00_sat_passed_plot,rhow_0665p00_sat_passed_plot,rhow_0673p75_sat_passed_plot,\
                            rhow_0681p25_sat_passed_plot,rhow_0708p75_sat_passed_plot,rhow_0753p75_sat_passed_plot,\
                            rhow_0778p75_sat_passed_plot,rhow_0865p00_sat_passed_plot,rhow_0885p00_sat_passed_plot],'+k')
                        scatter_legend.append('OLCI: ZMB18 (passed)')       


        #         else:
        #             print('CV exceeds criteria: CV[Lwn(560)]='+str(Lwn_560_CV)+'; CV[AOT(865.5)]='+str(AOT_0865p50_CV))
        #     else:
        #         print('Angles exceeds criteria: sza='+str(sza)+'; vza='+str(vza)+'; OR some pixels are flagged!')
        # else:
        #     print('Not matchups per '+year_str+' '+doy_str)

        plt.xlabel('Wavelength (nm)',fontsize=12)
        plt.ylabel(r'$\rho_{W}$',fontsize=12)
        plt.title('OLCI Time: '+str(sat_stop_time)[:-10])

#############################################
        # Bailey and Werdell 2006 
        delta_time = 3# float in hours       
        time_diff = ins_time_PAN - sat_stop_time
        dt_hour = [i.total_seconds()/(60*60) for i in time_diff] # time diffence between in situ measurements and sat in hours
        idx_min_PAN = np.argmin(np.abs(dt_hour))
        matchup_idx_vec_PAN = np.abs(dt_hour) <= delta_time 

        nday_PAN_ba = sum(matchup_idx_vec_PAN)
        if nday_PAN_ba >=1:
            print('--Bailey and Werdell 2006')
            print(str(nday_PAN_ba)+' matchups per '+year_str+' '+doy_str)
#           print(Lwn_fonQ[idx_min_PAN,:])
#           print(Exact_wavelengths[idx_min_PAN,:])
            
            
            center_px = int(len(rhow_0412p50_fq)/2 + 0.5)
            size_box = 5
            NTP = size_box*size_box # Number Total Pixels, excluding land pixels, Bailey and Werdell 2006
            start_idx_x = int(center_px-int(size_box/2))
            stop_idx_x = int(center_px+int(size_box/2)+1)
            start_idx_y = int(center_px-int(size_box/2))
            stop_idx_y = int(center_px+int(size_box/2)+1)
            rhow_0400p00_box = rhow_0400p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0412p50_fq_box = rhow_0412p50_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0442p50_fq_box = rhow_0442p50_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0490p00_fq_box = rhow_0490p00_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0510p00_fq_box = rhow_0510p00_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0560p00_fq_box = rhow_0560p00_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0620p00_fq_box = rhow_0620p00_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0665p00_fq_box = rhow_0665p00_fq[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0673p75_box = rhow_0673p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0681p25_box = rhow_0681p25[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0708p75_box = rhow_0708p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0753p75_box = rhow_0753p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0778p75_box = rhow_0778p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0865p00_box = rhow_0865p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_0885p00_box = rhow_0885p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            rhow_1020p50_box = rhow_1020p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]

            AOT_0865p50_box = AOT_0865p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]
            
            print('rhow_0412p50_fq_box:')
            print(rhow_0412p50_fq_box)
            print('rhow_0412p50_fq_box.mask:')
            print(rhow_0412p50_fq_box.mask)
            
            flags_mask = OLCI_flags.create_mask(WQSF[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])
            print('flags_mask:')
            print(flags_mask)
            # plot sat data
                # numpy.all(a, axis=None, out=None, keepdims=<no value>)[source]
                # Test whether all array elements along a given axis evaluate to True.
                # When all .mask are True, i.e. all invalid, .all() is True. Or if at least one is valid -> .all() is False
                # Therefore, if at least one band extract has a valid number, then plot it. 
                # This is plot the extract before applying filtering criteria.
            if not rhow_0400p00_box.mask.all() or not rhow_0412p50_fq_box.mask.all()\
                or not rhow_0442p50_fq_box.mask.all() or not rhow_0490p00_fq_box.mask.all()\
                or not rhow_0510p00_fq_box.mask.all() or not rhow_0560p00_fq_box.mask.all()\
                or not rhow_0620p00_fq_box.mask.all() or not rhow_0665p00_fq_box.mask.all()\
                or not rhow_0673p75_box.mask.all() or not rhow_0681p25_box.mask.all()\
                or not rhow_0708p75_box.mask.all() or not rhow_0753p75_box.mask.all()\
                or not rhow_0778p75_box.mask.all() or not rhow_0865p00_box.mask.all()\
                or not rhow_0885p00_box.mask.all() or not rhow_1020p50_box.mask.all():
                # plot BW06
                plt.plot([400.00,412.50,442.50,490.00,510.00,560.00,620.00,665.00,673.75,681.25,708.75,753.75,778.75,865.00,885.00,1020.5],\
                    [rhow_0400p00_box.mean(),rhow_0412p50_fq_box.mean(),rhow_0442p50_fq_box.mean(),rhow_0490p00_fq_box.mean(),\
                    rhow_0510p00_fq_box.mean(),rhow_0560p00_fq_box.mean(),rhow_0620p00_fq_box.mean(),rhow_0665p00_fq_box.mean(),\
                    rhow_0673p75_box.mean(),rhow_0681p25_box.mean(),rhow_0708p75_box.mean(),rhow_0753p75_box.mean(),\
                    rhow_0778p75_box.mean(),rhow_0865p00_box.mean(),rhow_0885p00_box.mean(),rhow_1020p50_box.mean()],'xr',markersize=8,linewidth=3)
                scatter_legend.append('OLCI: BW06')
            
            NGP = np.count_nonzero(flags_mask == 0) # Number Good Pixels, Bailey and Werdell 2006
            
            if sza<=75 and vza<=60 and NGP>NTP/2+1:

                # if nan, change mask
                rhow_0400p00_box = ma.masked_invalid(rhow_0400p00_box)
                rhow_0412p50_fq_box = ma.masked_invalid(rhow_0412p50_fq_box)
                rhow_0442p50_fq_box = ma.masked_invalid(rhow_0442p50_fq_box)
                rhow_0490p00_fq_box = ma.masked_invalid(rhow_0490p00_fq_box)
                rhow_0510p00_fq_box = ma.masked_invalid(rhow_0510p00_fq_box)
                rhow_0560p00_fq_box = ma.masked_invalid(rhow_0560p00_fq_box)
                rhow_0620p00_fq_box = ma.masked_invalid(rhow_0620p00_fq_box)
                rhow_0665p00_fq_box = ma.masked_invalid(rhow_0665p00_fq_box)
                rhow_0673p75_box = ma.masked_invalid(rhow_0673p75_box)
                rhow_0681p25_box = ma.masked_invalid(rhow_0681p25_box)
                rhow_0708p75_box = ma.masked_invalid(rhow_0708p75_box)
                rhow_0753p75_box = ma.masked_invalid(rhow_0753p75_box)
                rhow_0778p75_box = ma.masked_invalid(rhow_0778p75_box)
                rhow_0865p00_box = ma.masked_invalid(rhow_0865p00_box)
                rhow_0885p00_box = ma.masked_invalid(rhow_0885p00_box)
                rhow_1020p50_box = ma.masked_invalid(rhow_1020p50_box)
            
                AOT_0865p50_box = ma.masked_invalid(AOT_0865p50_box)

                NGP_rhow_0400p00 = np.count_nonzero(rhow_0400p00_box.mask == 0)
                NGP_rhow_0412p50 = np.count_nonzero(rhow_0412p50_fq_box.mask == 0)
                NGP_rhow_0442p50 = np.count_nonzero(rhow_0442p50_fq_box.mask == 0)
                NGP_rhow_0490p00 = np.count_nonzero(rhow_0490p00_fq_box.mask == 0)
                NGP_rhow_0510p00 = np.count_nonzero(rhow_0510p00_fq_box.mask == 0)
                NGP_rhow_0560p00 = np.count_nonzero(rhow_0560p00_fq_box.mask == 0)
                NGP_rhow_0620p00 = np.count_nonzero(rhow_0620p00_fq_box.mask == 0)
                NGP_rhow_0665p00 = np.count_nonzero(rhow_0665p00_fq_box.mask == 0)
                NGP_rhow_0673p75 = np.count_nonzero(rhow_0673p75_box.mask == 0)
                NGP_rhow_0681p25 = np.count_nonzero(rhow_0681p25_box.mask == 0)
                NGP_rhow_0708p75 = np.count_nonzero(rhow_0708p75_box.mask == 0)
                NGP_rhow_0753p75 = np.count_nonzero(rhow_0753p75_box.mask == 0)
                NGP_rhow_0778p75 = np.count_nonzero(rhow_0778p75_box.mask == 0)
                NGP_rhow_0865p00 = np.count_nonzero(rhow_0865p00_box.mask == 0)
                NGP_rhow_0885p00 = np.count_nonzero(rhow_0885p00_box.mask == 0)
                NGP_rhow_1020p50 = np.count_nonzero(rhow_1020p50_box.mask == 0)

                NGP_AOT_0865p50 = np.count_nonzero(AOT_0865p50_box.mask == 0)

                mean_unfiltered_rhow_0400p00 = rhow_0400p00_box.mean()
                mean_unfiltered_rhow_0412p50 = rhow_0412p50_fq_box.mean()
                mean_unfiltered_rhow_0442p50 = rhow_0442p50_fq_box.mean()
                mean_unfiltered_rhow_0490p00 = rhow_0490p00_fq_box.mean()
                mean_unfiltered_rhow_0510p00 = rhow_0510p00_fq_box.mean()
                mean_unfiltered_rhow_0560p00 = rhow_0560p00_fq_box.mean()
                mean_unfiltered_rhow_0620p00 = rhow_0620p00_fq_box.mean()
                mean_unfiltered_rhow_0665p00 = rhow_0665p00_fq_box.mean()
                mean_unfiltered_rhow_0673p75 = rhow_0673p75_box.mean()
                mean_unfiltered_rhow_0681p25 = rhow_0681p25_box.mean()
                mean_unfiltered_rhow_0708p75 = rhow_0708p75_box.mean()
                mean_unfiltered_rhow_0753p75 = rhow_0753p75_box.mean()
                mean_unfiltered_rhow_0778p75 = rhow_0778p75_box.mean()
                mean_unfiltered_rhow_0865p00 = rhow_0865p00_box.mean()
                mean_unfiltered_rhow_0885p00 = rhow_0885p00_box.mean()
                mean_unfiltered_rhow_1020p50 = rhow_1020p50_box.mean()
                
                mean_unfiltered_AOT_0865p50 = AOT_0865p50_box.mean()

                std_unfiltered_rhow_0400p00 = rhow_0400p00_box.std()
                std_unfiltered_rhow_0412p50 = rhow_0412p50_fq_box.std()
                std_unfiltered_rhow_0442p50 = rhow_0442p50_fq_box.std()
                std_unfiltered_rhow_0490p00 = rhow_0490p00_fq_box.std()
                std_unfiltered_rhow_0510p00 = rhow_0510p00_fq_box.std()
                std_unfiltered_rhow_0560p00 = rhow_0560p00_fq_box.std()
                std_unfiltered_rhow_0620p00 = rhow_0620p00_fq_box.std()
                std_unfiltered_rhow_0665p00 = rhow_0665p00_fq_box.std()
                std_unfiltered_rhow_0673p75 = rhow_0673p75_box.std()
                std_unfiltered_rhow_0681p25 = rhow_0681p25_box.std()
                std_unfiltered_rhow_0708p75 = rhow_0708p75_box.std()
                std_unfiltered_rhow_0753p75 = rhow_0753p75_box.std()
                std_unfiltered_rhow_0778p75 = rhow_0778p75_box.std()
                std_unfiltered_rhow_0865p00 = rhow_0865p00_box.std()
                std_unfiltered_rhow_0885p00 = rhow_0885p00_box.std()
                std_unfiltered_rhow_1020p50 = rhow_1020p50_box.std()
                
                std_unfiltered_AOT_0865p50 = AOT_0865p50_box.std()

                # mask values that are not within +/- 1.5*std of mean\ 
                rhow_0400p00_box = ma.masked_outside(rhow_0400p00_box,mean_unfiltered_rhow_0400p00\
                    -1.5*std_unfiltered_rhow_0400p00\
                    , mean_unfiltered_rhow_0400p00\
                    +1.5*std_unfiltered_rhow_0412p50)
                rhow_0412p50_fq_box = ma.masked_outside(rhow_0412p50_fq_box,mean_unfiltered_rhow_0412p50\
                    -1.5*std_unfiltered_rhow_0412p50\
                    , mean_unfiltered_rhow_0412p50\
                    +1.5*std_unfiltered_rhow_0412p50)
                rhow_0442p50_fq_box = ma.masked_outside(rhow_0442p50_fq_box,mean_unfiltered_rhow_0442p50\
                    -1.5*std_unfiltered_rhow_0442p50\
                    , mean_unfiltered_rhow_0442p50\
                    +1.5*std_unfiltered_rhow_0442p50)
                rhow_0490p00_fq_box = ma.masked_outside(rhow_0490p00_fq_box,mean_unfiltered_rhow_0490p00\
                    -1.5*std_unfiltered_rhow_0490p00\
                    , mean_unfiltered_rhow_0490p00\
                    +1.5*std_unfiltered_rhow_0490p00)
                rhow_0510p00_fq_box = ma.masked_outside(rhow_0510p00_fq_box,mean_unfiltered_rhow_0510p00\
                    -1.5*std_unfiltered_rhow_0510p00\
                    , mean_unfiltered_rhow_0510p00\
                    +1.5*std_unfiltered_rhow_0510p00)
                rhow_0560p00_fq_box = ma.masked_outside(rhow_0560p00_fq_box,mean_unfiltered_rhow_0560p00\
                    -1.5*std_unfiltered_rhow_0560p00\
                    , mean_unfiltered_rhow_0560p00\
                    +1.5*std_unfiltered_rhow_0560p00)
                rhow_0620p00_fq_box = ma.masked_outside(rhow_0620p00_fq_box,mean_unfiltered_rhow_0620p00\
                    -1.5*std_unfiltered_rhow_0620p00\
                    , mean_unfiltered_rhow_0620p00\
                    +1.5*std_unfiltered_rhow_0620p00)
                rhow_0665p00_fq_box = ma.masked_outside(rhow_0665p00_fq_box,mean_unfiltered_rhow_0665p00\
                    -1.5*std_unfiltered_rhow_0665p00\
                    , mean_unfiltered_rhow_0665p00\
                    +1.5*std_unfiltered_rhow_0665p00)
                rhow_0673p75_box = ma.masked_outside(rhow_0673p75_box,mean_unfiltered_rhow_0673p75\
                    -1.5*std_unfiltered_rhow_0673p75\
                    , mean_unfiltered_rhow_0673p75\
                    +1.5*std_unfiltered_rhow_0673p75)
                rhow_0708p75_box = ma.masked_outside(rhow_0708p75_box,mean_unfiltered_rhow_0708p75\
                    -1.5*std_unfiltered_rhow_0708p75\
                    , mean_unfiltered_rhow_0708p75\
                    +1.5*std_unfiltered_rhow_0708p75)
                rhow_0753p75_box = ma.masked_outside(rhow_0753p75_box,mean_unfiltered_rhow_0753p75\
                    -1.5*std_unfiltered_rhow_0753p75\
                    , mean_unfiltered_rhow_0753p75\
                    +1.5*std_unfiltered_rhow_0753p75)
                rhow_0778p75_box = ma.masked_outside(rhow_0778p75_box,mean_unfiltered_rhow_0778p75\
                    -1.5*std_unfiltered_rhow_0778p75\
                    , mean_unfiltered_rhow_0778p75\
                    +1.5*std_unfiltered_rhow_0778p75)
                rhow_0865p00_box = ma.masked_outside(rhow_0865p00_box,mean_unfiltered_rhow_0865p00\
                    -1.5*std_unfiltered_rhow_0865p00\
                    , mean_unfiltered_rhow_0865p00\
                    +1.5*std_unfiltered_rhow_0865p00)
                rhow_0885p00_box = ma.masked_outside(rhow_0885p00_box,mean_unfiltered_rhow_0885p00\
                    -1.5*std_unfiltered_rhow_0885p00\
                    , mean_unfiltered_rhow_0885p00\
                    +1.5*std_unfiltered_rhow_0885p00)
                rhow_1020p50_box = ma.masked_outside(rhow_1020p50_box,mean_unfiltered_rhow_1020p50\
                    -1.5*std_unfiltered_rhow_1020p50\
                    , mean_unfiltered_rhow_1020p50\
                    +1.5*std_unfiltered_rhow_1020p50)

                AOT_0865p50_box = ma.masked_outside(AOT_0865p50_box,mean_unfiltered_AOT_0865p50\
                    -1.5*std_unfiltered_AOT_0865p50\
                    , mean_unfiltered_AOT_0865p50\
                    +1.5*std_unfiltered_AOT_0865p50)

                mean_filtered_rhow_0400p00 = rhow_0400p00_box.mean()
                mean_filtered_rhow_0412p50 = rhow_0412p50_fq_box.mean()
                mean_filtered_rhow_0442p50 = rhow_0442p50_fq_box.mean()
                mean_filtered_rhow_0490p00 = rhow_0490p00_fq_box.mean()
                mean_filtered_rhow_0510p00 = rhow_0510p00_fq_box.mean()
                mean_filtered_rhow_0560p00 = rhow_0560p00_fq_box.mean()
                mean_filtered_rhow_0620p00 = rhow_0620p00_fq_box.mean()
                mean_filtered_rhow_0665p00 = rhow_0665p00_fq_box.mean()
                mean_filtered_rhow_0673p75 = rhow_0673p75_box.mean()
                mean_filtered_rhow_0681p25 = rhow_0681p25_box.mean()
                mean_filtered_rhow_0708p75 = rhow_0708p75_box.mean()
                mean_filtered_rhow_0753p75 = rhow_0753p75_box.mean()
                mean_filtered_rhow_0778p75 = rhow_0778p75_box.mean()
                mean_filtered_rhow_0865p00 = rhow_0865p00_box.mean()
                mean_filtered_rhow_0885p00 = rhow_0885p00_box.mean()
                mean_filtered_rhow_1020p50 = rhow_1020p50_box.mean()

                mean_filtered_AOT_0865p50  = AOT_0865p50_box.mean()

                std_filtered_rhow_0400p00 = rhow_0400p00_box.std()
                std_filtered_rhow_0412p50 = rhow_0412p50_fq_box.std()
                std_filtered_rhow_0442p50 = rhow_0442p50_fq_box.std()
                std_filtered_rhow_0490p00 = rhow_0490p00_fq_box.std()
                std_filtered_rhow_0510p00 = rhow_0510p00_fq_box.std()
                std_filtered_rhow_0560p00 = rhow_0560p00_fq_box.std()
                std_filtered_rhow_0620p00 = rhow_0620p00_fq_box.std()
                std_filtered_rhow_0665p00 = rhow_0665p00_fq_box.std()
                std_filtered_rhow_0673p75 = rhow_0673p75_box.std()
                std_filtered_rhow_0681p25 = rhow_0681p25_box.std()
                std_filtered_rhow_0708p75 = rhow_0708p75_box.std()
                std_filtered_rhow_0753p75 = rhow_0753p75_box.std()
                std_filtered_rhow_0778p75 = rhow_0778p75_box.std()
                std_filtered_rhow_0865p00 = rhow_0865p00_box.std()
                std_filtered_rhow_0885p00 = rhow_0885p00_box.std()
                std_filtered_rhow_1020p50 = rhow_1020p50_box.std()
                
                std_filtered_AOT_0865p50  = AOT_0865p50_box.std()

                CV_filtered_rhow_0400p00 = std_filtered_rhow_0400p00/mean_filtered_rhow_0400p00
                CV_filtered_rhow_0412p50 = std_filtered_rhow_0412p50/mean_filtered_rhow_0412p50
                CV_filtered_rhow_0442p50 = std_filtered_rhow_0442p50/mean_filtered_rhow_0442p50
                CV_filtered_rhow_0490p00 = std_filtered_rhow_0490p00/mean_filtered_rhow_0490p00
                CV_filtered_rhow_0510p00 = std_filtered_rhow_0510p00/mean_filtered_rhow_0510p00
                CV_filtered_rhow_0560p00 = std_filtered_rhow_0560p00/mean_filtered_rhow_0560p00
                CV_filtered_rhow_0620p00 = std_filtered_rhow_0620p00/mean_filtered_rhow_0620p00
                CV_filtered_rhow_0665p00 = std_filtered_rhow_0665p00/mean_filtered_rhow_0665p00
                CV_filtered_rhow_0673p75 = std_filtered_rhow_0673p75/mean_filtered_rhow_0673p75
                CV_filtered_rhow_0681p25 = std_filtered_rhow_0681p25/mean_filtered_rhow_0681p25
                CV_filtered_rhow_0708p75 = std_filtered_rhow_0708p75/mean_filtered_rhow_0708p75
                CV_filtered_rhow_0753p75 = std_filtered_rhow_0753p75/mean_filtered_rhow_0753p75
                CV_filtered_rhow_0778p75 = std_filtered_rhow_0778p75/mean_filtered_rhow_0778p75
                CV_filtered_rhow_0865p00 = std_filtered_rhow_0865p00/mean_filtered_rhow_0865p00
                CV_filtered_rhow_0885p00 = std_filtered_rhow_0885p00/mean_filtered_rhow_0885p00
                CV_filtered_rhow_1020p50 = std_filtered_rhow_1020p50/mean_filtered_rhow_1020p50
                
                CV_filtered_AOT_0865p50  = std_filtered_AOT_0865p50/mean_filtered_AOT_0865p50  
                
                CVs = [CV_filtered_rhow_0412p50, CV_filtered_rhow_0442p50,\
                                     CV_filtered_rhow_0490p00, CV_filtered_rhow_0560p00,\
                                     CV_filtered_AOT_0865p50]
                print(CVs)
                MedianCV = np.nanmedian(np.abs(CVs))

                print('Median CV='+str(MedianCV))

                rhow_0400p00_sat_passed_plot = np.nan
                rhow_0412p50_sat_passed_plot = np.nan
                rhow_0442p50_sat_passed_plot = np.nan
                rhow_0490p00_sat_passed_plot = np.nan
                rhow_0510p00_sat_passed_plot = np.nan
                rhow_0560p00_sat_passed_plot = np.nan
                rhow_0620p00_sat_passed_plot = np.nan
                rhow_0665p00_sat_passed_plot = np.nan
                rhow_0673p75_sat_passed_plot = np.nan
                rhow_0681p25_sat_passed_plot = np.nan
                rhow_0708p75_sat_passed_plot = np.nan
                rhow_0753p75_sat_passed_plot = np.nan
                rhow_0778p75_sat_passed_plot = np.nan
                rhow_0865p00_sat_passed_plot = np.nan
                rhow_0885p00_sat_passed_plot = np.nan

                if MedianCV <= 0.15:
                    # Rrs 0400p00
                    # print('400.0')
                    if not NGP_rhow_0400p00<NTP/2+1:
                        rhow_0400p00_sat_passed_plot = mean_filtered_rhow_0400p00
                        # print('Exceeded: NGP_rhow_0400p00='+str(NGP_rhow_0400p00))
                    # else:
                        matchups_PAN_rhow_0400p00_sat_ba.append(mean_filtered_rhow_0400p00)
                        matchups_PAN_rhow_0400p00_ins_ba.append(rhow_0400p00_PAN) # 412,
                        matchups_PAN_rhow_0400p00_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0400p00_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0400p00_ins_ba_time.append(ins_time_PAN[idx_min_PAN])

                    # Rrs 0412p50
                    # print('412.5')
                    if not NGP_rhow_0412p50<NTP/2+1:
                        rhow_0412p50_sat_passed_plot = mean_filtered_rhow_0412p50
                        # print('Exceeded: NGP_rhow_0412p50='+str(NGP_rhow_0412p50))
                    # else:
                        matchups_PAN_rhow_0412p50_fq_sat_ba.append(mean_filtered_rhow_0412p50)
                        matchups_PAN_rhow_0412p50_fq_ins_ba.append(rhow_0412p50_PAN) # 412,
                        matchups_PAN_rhow_0412p50_fq_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0412p50_fq_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0412p50_fq_ins_ba_time.append(ins_time_PAN[idx_min_PAN])
                        
                    # Rrs 0442p50
                    # print('442.5')
                    if not NGP_rhow_0442p50<NTP/2+1:
                        rhow_0442p50_sat_passed_plot = mean_filtered_rhow_0442p50
                        # print('Exceeded: NGP_rhow_0442p50='+str(NGP_rhow_0442p50))
                    # else:
                        matchups_PAN_rhow_0442p50_fq_sat_ba.append(mean_filtered_rhow_0442p50)
                        matchups_PAN_rhow_0442p50_fq_ins_ba.append(rhow_0442p50_PAN) # 441.8
                        matchups_PAN_rhow_0442p50_fq_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0442p50_fq_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0442p50_fq_ins_ba_time.append(ins_time_PAN[idx_min_PAN])
                        
                    # Rrs 0490p00
                    # print('490.0')
                    if not NGP_rhow_0490p00<NTP/2+1:
                        rhow_0490p00_sat_passed_plot = mean_filtered_rhow_0490p00
                        # print('Exceeded: NGP_rhow_0490p00='+str(NGP_rhow_0490p00))
                    # else:
                        matchups_PAN_rhow_0490p00_fq_sat_ba.append(mean_filtered_rhow_0490p00)
                        matchups_PAN_rhow_0490p00_fq_ins_ba.append(rhow_0490p00_PAN) # 488.5
                        matchups_PAN_rhow_0490p00_fq_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0490p00_fq_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0490p00_fq_ins_ba_time.append(ins_time_PAN[idx_min_PAN])

                    # Rrs 0510p00
                    # print('510.0')
                    if not NGP_rhow_0510p00<NTP/2+1:
                        rhow_0510p00_sat_passed_plot = mean_filtered_rhow_0510p00
                        # print('Exceeded: NGP_rhow_0510p00='+str(NGP_rhow_0510p00))
                    # else:
                        matchups_PAN_rhow_0510p00_fq_sat_ba.append(mean_filtered_rhow_0510p00)
                        matchups_PAN_rhow_0510p00_fq_ins_ba.append(rhow_0510p00_PAN) # 488.5
                        matchups_PAN_rhow_0510p00_fq_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0510p00_fq_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0510p00_fq_ins_ba_time.append(ins_time_PAN[idx_min_PAN])    
                        
                    # Rrs 0560p00
                    # print('560.0')
                    if not NGP_rhow_0560p00<NTP/2+1:
                        rhow_0560p00_sat_passed_plot = mean_filtered_rhow_0560p00
                        # print('Exceeded: NGP_rhow_0560p00='+str(NGP_rhow_0560p00))
                    # else:
                        matchups_PAN_rhow_0560p00_fq_sat_ba.append(mean_filtered_rhow_0560p00)
                        matchups_PAN_rhow_0560p00_fq_ins_ba.append(rhow_0560p00_PAN) # 551,
                        matchups_PAN_rhow_0560p00_fq_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0560p00_fq_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0560p00_fq_ins_ba_time.append(ins_time_PAN[idx_min_PAN])

                    # Rrs 0620p00
                    # print('620.0')
                    if not NGP_rhow_0620p00<NTP/2+1:
                        rhow_0620p00_sat_passed_plot = mean_filtered_rhow_0620p00
                        # print('Exceeded: NGP_rhow_0620p00='+str(NGP_rhow_0620p00))
                    # else:
                        matchups_PAN_rhow_0620p00_fq_sat_ba.append(mean_filtered_rhow_0620p00)
                        matchups_PAN_rhow_0620p00_fq_ins_ba.append(rhow_0620p00_PAN) # 551,
                        matchups_PAN_rhow_0620p00_fq_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0620p00_fq_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0620p00_fq_ins_ba_time.append(ins_time_PAN[idx_min_PAN])
                        
                    # Rrs 0665p00
                    # print('665.0')
                    if not NGP_rhow_0665p00<NTP/2+1:
                        rhow_0665p00_sat_passed_plot = mean_filtered_rhow_0665p00
                        # print('Exceeded: NGP_rhow_0665p00='+str(NGP_rhow_0665p00))
                    # else:
                        matchups_PAN_rhow_0665p00_fq_sat_ba.append(mean_filtered_rhow_0665p00)
                        matchups_PAN_rhow_0665p00_fq_ins_ba.append(rhow_0665p00_PAN) # 667.9    
                        matchups_PAN_rhow_0665p00_fq_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0665p00_fq_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0665p00_fq_ins_ba_time.append(ins_time_PAN[idx_min_PAN])

                    # Rrs 0673p75
                    # print('673.75')
                    if not NGP_rhow_0673p75<NTP/2+1:
                        rhow_0673p75_sat_passed_plot = mean_filtered_rhow_0673p75
                        # print('Exceeded: NGP_rhow_0673p75='+str(NGP_rhow_0673p75))
                    # else:
                        matchups_PAN_rhow_0673p75_sat_ba.append(mean_filtered_rhow_0673p75)
                        matchups_PAN_rhow_0673p75_ins_ba.append(rhow_0673p75_PAN) # 412,
                        matchups_PAN_rhow_0673p75_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0673p75_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0673p75_ins_ba_time.append(ins_time_PAN[idx_min_PAN])
                    
                    # Rrs 0681p25
                    # print('681.25')
                    if not NGP_rhow_0681p25<NTP/2+1:
                        rhow_0681p25_sat_passed_plot = mean_filtered_rhow_0681p25
                        # print('Exceeded: NGP_rhow_0681p25='+str(NGP_rhow_0681p25))
                    # else:
                        matchups_PAN_rhow_0681p25_sat_ba.append(mean_filtered_rhow_0681p25)
                        matchups_PAN_rhow_0681p25_ins_ba.append(rhow_0681p25_PAN) # 412,
                        matchups_PAN_rhow_0681p25_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0681p25_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0681p25_ins_ba_time.append(ins_time_PAN[idx_min_PAN])

                    # Rrs 0708p75
                    # print('708.75')
                    if not NGP_rhow_0708p75<NTP/2+1:
                        rhow_0708p75_sat_passed_plot = mean_filtered_rhow_0708p75
                        # print('Exceeded: NGP_rhow_0708p75='+str(NGP_rhow_0708p75))
                    # else:
                        matchups_PAN_rhow_0708p75_sat_ba.append(mean_filtered_rhow_0708p75)
                        matchups_PAN_rhow_0708p75_ins_ba.append(rhow_0708p75_PAN) # 412,
                        matchups_PAN_rhow_0708p75_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0708p75_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0708p75_ins_ba_time.append(ins_time_PAN[idx_min_PAN])
                        
                    # Rrs 0753p75
                    # print('753.75')
                    if not NGP_rhow_0753p75<NTP/2+1:
                        rhow_0753p75_sat_passed_plot = mean_filtered_rhow_0753p75
                        # print('Exceeded: NGP_rhow_0753p75='+str(NGP_rhow_0753p75))
                    # else:
                        matchups_PAN_rhow_0753p75_sat_ba.append(mean_filtered_rhow_0753p75)
                        matchups_PAN_rhow_0753p75_ins_ba.append(rhow_0753p75_PAN) # 412,
                        matchups_PAN_rhow_0753p75_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0753p75_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0753p75_ins_ba_time.append(ins_time_PAN[idx_min_PAN])
                        
                    # Rrs 0778p75
                    # print('778.75')
                    if not NGP_rhow_0778p75<NTP/2+1:
                        rhow_0778p75_sat_passed_plot = mean_filtered_rhow_0778p75
                        # print('Exceeded: NGP_rhow_0778p75='+str(NGP_rhow_0778p75))
                    # else:
                        matchups_PAN_rhow_0778p75_sat_ba.append(mean_filtered_rhow_0778p75)
                        matchups_PAN_rhow_0778p75_ins_ba.append(rhow_0778p75_PAN) # 412,
                        matchups_PAN_rhow_0778p75_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0778p75_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0778p75_ins_ba_time.append(ins_time_PAN[idx_min_PAN])
                        
                    # Rrs 0865p00
                    # print('865.0')
                    if not NGP_rhow_0865p00<NTP/2+1:
                        rhow_0865p00_sat_passed_plot = mean_filtered_rhow_0865p00
                        # print('Exceeded: NGP_rhow_0865p00='+str(NGP_rhow_0865p00))
                    # else:
                        matchups_PAN_rhow_0865p00_sat_ba.append(mean_filtered_rhow_0865p00)
                        matchups_PAN_rhow_0865p00_ins_ba.append(rhow_0865p00_PAN) # 412,
                        matchups_PAN_rhow_0865p00_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0865p00_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0865p00_ins_ba_time.append(ins_time_PAN[idx_min_PAN])
                        
                    # Rrs 0885p00
                    # print('885.0')
                    if not NGP_rhow_0885p00<NTP/2+1:
                        rhow_0885p00_sat_passed_plot = mean_filtered_rhow_0885p00
                        # print('Exceeded: NGP_rhow_0885p00='+str(NGP_rhow_0885p00))
                    # else:
                        matchups_PAN_rhow_0885p00_sat_ba.append(mean_filtered_rhow_0885p00)
                        matchups_PAN_rhow_0885p00_ins_ba.append(rhow_0885p00_PAN) # 412,
                        matchups_PAN_rhow_0885p00_ins_ba_station.append(station_name)
                        matchups_PAN_rhow_0885p00_sat_ba_stop_time.append(sat_stop_time)
                        matchups_PAN_rhow_0885p00_ins_ba_time.append(ins_time_PAN[idx_min_PAN])
                    
                    # if any of the bands passes the criteria, plot it
                    if (not NGP_rhow_0400p00<NTP/2+1 or not NGP_rhow_0412p50<NTP/2+1 or not NGP_rhow_0442p50<NTP/2+1\
                        or not NGP_rhow_0490p00<NTP/2+1 or not NGP_rhow_0510p00<NTP/2+1 or not NGP_rhow_0560p00<NTP/2+1\
                        or not NGP_rhow_0620p00<NTP/2+1 or not NGP_rhow_0665p00<NTP/2+1 or not NGP_rhow_0673p75<NTP/2+1\
                        or not NGP_rhow_0681p25<NTP/2+1 or not NGP_rhow_0708p75<NTP/2+1 or not NGP_rhow_0753p75<NTP/2+1\
                        or not NGP_rhow_0778p75<NTP/2+1 or not NGP_rhow_0865p00<NTP/2+1 or not NGP_rhow_0885p00<NTP/2+1):
                        plt.plot([400.00,412.50,442.50,490.00,510.00,560.00,620.00,665.00,673.75,681.25,708.75,753.75,778.75,865.00,885.00],\
                            [rhow_0400p00_sat_passed_plot,rhow_0412p50_sat_passed_plot,rhow_0442p50_sat_passed_plot,\
                            rhow_0490p00_sat_passed_plot,rhow_0510p00_sat_passed_plot,rhow_0560p00_sat_passed_plot,\
                            rhow_0620p00_sat_passed_plot,rhow_0665p00_sat_passed_plot,rhow_0673p75_sat_passed_plot,\
                            rhow_0681p25_sat_passed_plot,rhow_0708p75_sat_passed_plot,rhow_0753p75_sat_passed_plot,\
                            rhow_0778p75_sat_passed_plot,rhow_0865p00_sat_passed_plot,rhow_0885p00_sat_passed_plot],'xk')
                        scatter_legend.append('OLCI: BW06 (passed)')        
        #         else:
        #             print('Median CV exceeds criteria: Median[CV]='+str(MedianCV))
        #     else:
        #         print('Angles exceeds criteria: sza='+str(sza)+'; vza='+str(vza)+'; OR NGP='+str(NGP)+'< NTP/2+1='+str(NTP/2+1)+'!')
        # else:
        #     print('Not matchups per '+year_str+' '+doy_str)   


        plt.legend(scatter_legend,fontsize=12)
        ofname = year_str+'_'+doy_str+'_'+'_spectrum.pdf'
        ofname = os.path.join(path_out,'source',ofname)
        plt.savefig(ofname, dpi=300)


#%% plots   

prot_name = 'zi'
sensor_name = 'OLCI'
rmse_val_0400p00_zi, mean_abs_rel_diff_0400p00_zi, mean_rel_diff_0400p00_zi, r_sqr_0400p00_zi,\
rmse_val_0400p00_zi_Venise,mean_abs_rel_diff_0400p00_zi_Venise, mean_rel_diff_0400p00_zi_Venise, r_sqr_0400p00_zi_Venise,\
rmse_val_0400p00_zi_Gloria,mean_abs_rel_diff_0400p00_zi_Gloria, mean_rel_diff_0400p00_zi_Gloria, r_sqr_0400p00_zi_Gloria,\
rmse_val_0400p00_zi_Galata_Platform,mean_abs_rel_diff_0400p00_zi_Galata_Platform, mean_rel_diff_0400p00_zi_Galata_Platform, r_sqr_0400p00_zi_Galata_Platform,\
rmse_val_0400p00_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0400p00_zi_Helsinki_Lighthouse, mean_rel_diff_0400p00_zi_Helsinki_Lighthouse, r_sqr_0400p00_zi_Helsinki_Lighthouse,\
rmse_val_0400p00_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0400p00_zi_Gustav_Dalen_Tower, mean_rel_diff_0400p00_zi_Gustav_Dalen_Tower, r_sqr_0400p00_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0400p00_ins_zi,matchups_PAN_rhow_0400p00_sat_zi,'400.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0400p00_ins_zi_station,min_val=-0.01,max_val=0.1)

rmse_val_0412p50_zi, mean_abs_rel_diff_0412p50_zi, mean_rel_diff_0412p50_zi, r_sqr_0412p50_zi,\
rmse_val_0412p50_zi_Venise,mean_abs_rel_diff_0412p50_zi_Venise, mean_rel_diff_0412p50_zi_Venise, r_sqr_0412p50_zi_Venise,\
rmse_val_0412p50_zi_Gloria,mean_abs_rel_diff_0412p50_zi_Gloria, mean_rel_diff_0412p50_zi_Gloria, r_sqr_0412p50_zi_Gloria,\
rmse_val_0412p50_zi_Galata_Platform,mean_abs_rel_diff_0412p50_zi_Galata_Platform, mean_rel_diff_0412p50_zi_Galata_Platform, r_sqr_0412p50_zi_Galata_Platform,\
rmse_val_0412p50_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0412p50_zi_Helsinki_Lighthouse, mean_rel_diff_0412p50_zi_Helsinki_Lighthouse, r_sqr_0412p50_zi_Helsinki_Lighthouse,\
rmse_val_0412p50_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0412p50_zi_Gustav_Dalen_Tower, mean_rel_diff_0412p50_zi_Gustav_Dalen_Tower, r_sqr_0412p50_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0412p50_fq_ins_zi,matchups_PAN_rhow_0412p50_fq_sat_zi,'412.5',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0412p50_fq_ins_zi_station,min_val=-0.01,max_val=0.1)

rmse_val_0442p50_zi, mean_abs_rel_diff_0442p50_zi, mean_rel_diff_0442p50_zi, r_sqr_0442p50_zi,\
rmse_val_0442p50_zi_Venise,mean_abs_rel_diff_0442p50_zi_Venise, mean_rel_diff_0442p50_zi_Venise, r_sqr_0442p50_zi_Venise,\
rmse_val_0442p50_zi_Gloria,mean_abs_rel_diff_0442p50_zi_Gloria, mean_rel_diff_0442p50_zi_Gloria, r_sqr_0442p50_zi_Gloria,\
rmse_val_0442p50_zi_Galata_Platform,mean_abs_rel_diff_0442p50_zi_Galata_Platform, mean_rel_diff_0442p50_zi_Galata_Platform, r_sqr_0442p50_zi_Galata_Platform,\
rmse_val_0442p50_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0442p50_zi_Helsinki_Lighthouse, mean_rel_diff_0442p50_zi_Helsinki_Lighthouse, r_sqr_0442p50_zi_Helsinki_Lighthouse,\
rmse_val_0442p50_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0442p50_zi_Gustav_Dalen_Tower, mean_rel_diff_0442p50_zi_Gustav_Dalen_Tower, r_sqr_0442p50_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0442p50_fq_ins_zi,matchups_PAN_rhow_0442p50_fq_sat_zi,'442.5',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0442p50_fq_ins_zi_station,min_val=-0.01,max_val=0.1)

rmse_val_0490p00_zi, mean_abs_rel_diff_0490p00_zi, mean_rel_diff_0490p00_zi, r_sqr_0490p00_zi,\
rmse_val_0490p00_zi_Venise,mean_abs_rel_diff_0490p00_zi_Venise, mean_rel_diff_0490p00_zi_Venise, r_sqr_0490p00_zi_Venise,\
rmse_val_0490p00_zi_Gloria,mean_abs_rel_diff_0490p00_zi_Gloria, mean_rel_diff_0490p00_zi_Gloria, r_sqr_0490p00_zi_Gloria,\
rmse_val_0490p00_zi_Galata_Platform,mean_abs_rel_diff_0490p00_zi_Galata_Platform, mean_rel_diff_0490p00_zi_Galata_Platform, r_sqr_0490p00_zi_Galata_Platform,\
rmse_val_0490p00_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0490p00_zi_Helsinki_Lighthouse, mean_rel_diff_0490p00_zi_Helsinki_Lighthouse, r_sqr_0490p00_zi_Helsinki_Lighthouse,\
rmse_val_0490p00_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0490p00_zi_Gustav_Dalen_Tower, mean_rel_diff_0490p00_zi_Gustav_Dalen_Tower, r_sqr_0490p00_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0490p00_fq_ins_zi,matchups_PAN_rhow_0490p00_fq_sat_zi,'490.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0490p00_fq_ins_zi_station,min_val=-0.01,max_val=0.12)

rmse_val_0510p00_zi, mean_abs_rel_diff_0510p00_zi, mean_rel_diff_0510p00_zi, r_sqr_0510p00_zi,\
rmse_val_0510p00_zi_Venise,mean_abs_rel_diff_0510p00_zi_Venise, mean_rel_diff_0510p00_zi_Venise, r_sqr_0510p00_zi_Venise,\
rmse_val_0510p00_zi_Gloria,mean_abs_rel_diff_0510p00_zi_Gloria, mean_rel_diff_0510p00_zi_Gloria, r_sqr_0510p00_zi_Gloria,\
rmse_val_0510p00_zi_Galata_Platform,mean_abs_rel_diff_0510p00_zi_Galata_Platform, mean_rel_diff_0510p00_zi_Galata_Platform, r_sqr_0510p00_zi_Galata_Platform,\
rmse_val_0510p00_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0510p00_zi_Helsinki_Lighthouse, mean_rel_diff_0510p00_zi_Helsinki_Lighthouse, r_sqr_0510p00_zi_Helsinki_Lighthouse,\
rmse_val_0510p00_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0510p00_zi_Gustav_Dalen_Tower, mean_rel_diff_0510p00_zi_Gustav_Dalen_Tower, r_sqr_0510p00_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0510p00_fq_ins_zi,matchups_PAN_rhow_0510p00_fq_sat_zi,'510.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0510p00_fq_ins_zi_station,min_val=-0.01,max_val=0.12)

rmse_val_0560p00_zi, mean_abs_rel_diff_0560p00_zi, mean_rel_diff_0560p00_zi, r_sqr_0560p00_zi,\
rmse_val_0560p00_zi_Venise,mean_abs_rel_diff_0560p00_zi_Venise, mean_rel_diff_0560p00_zi_Venise, r_sqr_0560p00_zi_Venise,\
rmse_val_0560p00_zi_Gloria,mean_abs_rel_diff_0560p00_zi_Gloria, mean_rel_diff_0560p00_zi_Gloria, r_sqr_0560p00_zi_Gloria,\
rmse_val_0560p00_zi_Galata_Platform,mean_abs_rel_diff_0560p00_zi_Galata_Platform, mean_rel_diff_0560p00_zi_Galata_Platform, r_sqr_0560p00_zi_Galata_Platform,\
rmse_val_0560p00_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0560p00_zi_Helsinki_Lighthouse, mean_rel_diff_0560p00_zi_Helsinki_Lighthouse, r_sqr_0560p00_zi_Helsinki_Lighthouse,\
rmse_val_0560p00_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0560p00_zi_Gustav_Dalen_Tower, mean_rel_diff_0560p00_zi_Gustav_Dalen_Tower, r_sqr_0560p00_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0560p00_fq_ins_zi,matchups_PAN_rhow_0560p00_fq_sat_zi,'560.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0560p00_fq_ins_zi_station,min_val=-0.01,max_val=0.12)

rmse_val_0620p00_zi, mean_abs_rel_diff_0620p00_zi, mean_rel_diff_0620p00_zi, r_sqr_0620p00_zi,\
rmse_val_0620p00_zi_Venise,mean_abs_rel_diff_0620p00_zi_Venise, mean_rel_diff_0620p00_zi_Venise, r_sqr_0620p00_zi_Venise,\
rmse_val_0620p00_zi_Gloria,mean_abs_rel_diff_0620p00_zi_Gloria, mean_rel_diff_0620p00_zi_Gloria, r_sqr_0620p00_zi_Gloria,\
rmse_val_0620p00_zi_Galata_Platform,mean_abs_rel_diff_0620p00_zi_Galata_Platform, mean_rel_diff_0620p00_zi_Galata_Platform, r_sqr_0620p00_zi_Galata_Platform,\
rmse_val_0620p00_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0620p00_zi_Helsinki_Lighthouse, mean_rel_diff_0620p00_zi_Helsinki_Lighthouse, r_sqr_0620p00_zi_Helsinki_Lighthouse,\
rmse_val_0620p00_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0620p00_zi_Gustav_Dalen_Tower, mean_rel_diff_0620p00_zi_Gustav_Dalen_Tower, r_sqr_0620p00_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0620p00_fq_ins_zi,matchups_PAN_rhow_0620p00_fq_sat_zi,'620.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0620p00_fq_ins_zi_station,min_val=-0.01,max_val=0.08)

rmse_val_0665p00_zi, mean_abs_rel_diff_0665p00_zi, mean_rel_diff_0665p00_zi, r_sqr_0665p00_zi,\
rmse_val_0665p00_zi_Venise,mean_abs_rel_diff_0665p00_zi_Venise, mean_rel_diff_0665p00_zi_Venise, r_sqr_0665p00_zi_Venise,\
rmse_val_0665p00_zi_Gloria,mean_abs_rel_diff_0665p00_zi_Gloria, mean_rel_diff_0665p00_zi_Gloria, r_sqr_0665p00_zi_Gloria,\
rmse_val_0665p00_zi_Galata_Platform,mean_abs_rel_diff_0665p00_zi_Galata_Platform, mean_rel_diff_0665p00_zi_Galata_Platform, r_sqr_0665p00_zi_Galata_Platform,\
rmse_val_0665p00_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0665p00_zi_Helsinki_Lighthouse, mean_rel_diff_0665p00_zi_Helsinki_Lighthouse, r_sqr_0665p00_zi_Helsinki_Lighthouse,\
rmse_val_0665p00_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0665p00_zi_Gustav_Dalen_Tower, mean_rel_diff_0665p00_zi_Gustav_Dalen_Tower, r_sqr_0665p00_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0665p00_fq_ins_zi,matchups_PAN_rhow_0665p00_fq_sat_zi,'665.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0665p00_fq_ins_zi_station,min_val=-0.01,max_val=0.06)

rmse_val_0673p75_zi, mean_abs_rel_diff_0673p75_zi, mean_rel_diff_0673p75_zi, r_sqr_0673p75_zi,\
rmse_val_0673p75_zi_Venise,mean_abs_rel_diff_0673p75_zi_Venise, mean_rel_diff_0673p75_zi_Venise, r_sqr_0673p75_zi_Venise,\
rmse_val_0673p75_zi_Gloria,mean_abs_rel_diff_0673p75_zi_Gloria, mean_rel_diff_0673p75_zi_Gloria, r_sqr_0673p75_zi_Gloria,\
rmse_val_0673p75_zi_Galata_Platform,mean_abs_rel_diff_0673p75_zi_Galata_Platform, mean_rel_diff_0673p75_zi_Galata_Platform, r_sqr_0673p75_zi_Galata_Platform,\
rmse_val_0673p75_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0673p75_zi_Helsinki_Lighthouse, mean_rel_diff_0673p75_zi_Helsinki_Lighthouse, r_sqr_0673p75_zi_Helsinki_Lighthouse,\
rmse_val_0673p75_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0673p75_zi_Gustav_Dalen_Tower, mean_rel_diff_0673p75_zi_Gustav_Dalen_Tower, r_sqr_0673p75_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0673p75_ins_zi,matchups_PAN_rhow_0673p75_sat_zi,'673.75',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0673p75_ins_zi_station,min_val=-0.01,max_val=0.06)

rmse_val_0681p25_zi, mean_abs_rel_diff_0681p25_zi, mean_rel_diff_0681p25_zi, r_sqr_0681p25_zi,\
rmse_val_0681p25_zi_Venise,mean_abs_rel_diff_0681p25_zi_Venise, mean_rel_diff_0681p25_zi_Venise, r_sqr_0681p25_zi_Venise,\
rmse_val_0681p25_zi_Gloria,mean_abs_rel_diff_0681p25_zi_Gloria, mean_rel_diff_0681p25_zi_Gloria, r_sqr_0681p25_zi_Gloria,\
rmse_val_0681p25_zi_Galata_Platform,mean_abs_rel_diff_0681p25_zi_Galata_Platform, mean_rel_diff_0681p25_zi_Galata_Platform, r_sqr_0681p25_zi_Galata_Platform,\
rmse_val_0681p25_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0681p25_zi_Helsinki_Lighthouse, mean_rel_diff_0681p25_zi_Helsinki_Lighthouse, r_sqr_0681p25_zi_Helsinki_Lighthouse,\
rmse_val_0681p25_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0681p25_zi_Gustav_Dalen_Tower, mean_rel_diff_0681p25_zi_Gustav_Dalen_Tower, r_sqr_0681p25_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0681p25_ins_zi,matchups_PAN_rhow_0681p25_sat_zi,'681.25',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0681p25_ins_zi_station,min_val=-0.01,max_val=0.06)

rmse_val_0708p75_zi, mean_abs_rel_diff_0708p75_zi, mean_rel_diff_0708p75_zi, r_sqr_0708p75_zi,\
rmse_val_0708p75_zi_Venise,mean_abs_rel_diff_0708p75_zi_Venise, mean_rel_diff_0708p75_zi_Venise, r_sqr_0708p75_zi_Venise,\
rmse_val_0708p75_zi_Gloria,mean_abs_rel_diff_0708p75_zi_Gloria, mean_rel_diff_0708p75_zi_Gloria, r_sqr_0708p75_zi_Gloria,\
rmse_val_0708p75_zi_Galata_Platform,mean_abs_rel_diff_0708p75_zi_Galata_Platform, mean_rel_diff_0708p75_zi_Galata_Platform, r_sqr_0708p75_zi_Galata_Platform,\
rmse_val_0708p75_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0708p75_zi_Helsinki_Lighthouse, mean_rel_diff_0708p75_zi_Helsinki_Lighthouse, r_sqr_0708p75_zi_Helsinki_Lighthouse,\
rmse_val_0708p75_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0708p75_zi_Gustav_Dalen_Tower, mean_rel_diff_0708p75_zi_Gustav_Dalen_Tower, r_sqr_0708p75_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0708p75_ins_zi,matchups_PAN_rhow_0708p75_sat_zi,'708.75',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0708p75_ins_zi_station,min_val=-0.01,max_val=0.04)

rmse_val_0753p75_zi, mean_abs_rel_diff_0753p75_zi, mean_rel_diff_0753p75_zi, r_sqr_0753p75_zi,\
rmse_val_0753p75_zi_Venise,mean_abs_rel_diff_0753p75_zi_Venise, mean_rel_diff_0753p75_zi_Venise, r_sqr_0753p75_zi_Venise,\
rmse_val_0753p75_zi_Gloria,mean_abs_rel_diff_0753p75_zi_Gloria, mean_rel_diff_0753p75_zi_Gloria, r_sqr_0753p75_zi_Gloria,\
rmse_val_0753p75_zi_Galata_Platform,mean_abs_rel_diff_0753p75_zi_Galata_Platform, mean_rel_diff_0753p75_zi_Galata_Platform, r_sqr_0753p75_zi_Galata_Platform,\
rmse_val_0753p75_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0753p75_zi_Helsinki_Lighthouse, mean_rel_diff_0753p75_zi_Helsinki_Lighthouse, r_sqr_0753p75_zi_Helsinki_Lighthouse,\
rmse_val_0753p75_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0753p75_zi_Gustav_Dalen_Tower, mean_rel_diff_0753p75_zi_Gustav_Dalen_Tower, r_sqr_0753p75_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0753p75_ins_zi,matchups_PAN_rhow_0753p75_sat_zi,'753.75',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0753p75_ins_zi_station,min_val=-0.001,max_val=0.009)

rmse_val_0778p75_zi, mean_abs_rel_diff_0778p75_zi, mean_rel_diff_0778p75_zi, r_sqr_0778p75_zi,\
rmse_val_0778p75_zi_Venise,mean_abs_rel_diff_0778p75_zi_Venise, mean_rel_diff_0778p75_zi_Venise, r_sqr_0778p75_zi_Venise,\
rmse_val_0778p75_zi_Gloria,mean_abs_rel_diff_0778p75_zi_Gloria, mean_rel_diff_0778p75_zi_Gloria, r_sqr_0778p75_zi_Gloria,\
rmse_val_0778p75_zi_Galata_Platform,mean_abs_rel_diff_0778p75_zi_Galata_Platform, mean_rel_diff_0778p75_zi_Galata_Platform, r_sqr_0778p75_zi_Galata_Platform,\
rmse_val_0778p75_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0778p75_zi_Helsinki_Lighthouse, mean_rel_diff_0778p75_zi_Helsinki_Lighthouse, r_sqr_0778p75_zi_Helsinki_Lighthouse,\
rmse_val_0778p75_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0778p75_zi_Gustav_Dalen_Tower, mean_rel_diff_0778p75_zi_Gustav_Dalen_Tower, r_sqr_0778p75_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0778p75_ins_zi,matchups_PAN_rhow_0778p75_sat_zi,'778.75',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0778p75_ins_zi_station,min_val=-0.0002,max_val=0.009)

rmse_val_0865p00_zi, mean_abs_rel_diff_0865p00_zi, mean_rel_diff_0865p00_zi, r_sqr_0865p00_zi,\
rmse_val_0865p00_zi_Venise,mean_abs_rel_diff_0865p00_zi_Venise, mean_rel_diff_0865p00_zi_Venise, r_sqr_0865p00_zi_Venise,\
rmse_val_0865p00_zi_Gloria,mean_abs_rel_diff_0865p00_zi_Gloria, mean_rel_diff_0865p00_zi_Gloria, r_sqr_0865p00_zi_Gloria,\
rmse_val_0865p00_zi_Galata_Platform,mean_abs_rel_diff_0865p00_zi_Galata_Platform, mean_rel_diff_0865p00_zi_Galata_Platform, r_sqr_0865p00_zi_Galata_Platform,\
rmse_val_0865p00_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0865p00_zi_Helsinki_Lighthouse, mean_rel_diff_0865p00_zi_Helsinki_Lighthouse, r_sqr_0865p00_zi_Helsinki_Lighthouse,\
rmse_val_0865p00_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0865p00_zi_Gustav_Dalen_Tower, mean_rel_diff_0865p00_zi_Gustav_Dalen_Tower, r_sqr_0865p00_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0865p00_ins_zi,matchups_PAN_rhow_0865p00_sat_zi,'865.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0865p00_ins_zi_station,min_val=-0.001,max_val=0.006)

rmse_val_0885p00_zi, mean_abs_rel_diff_0885p00_zi, mean_rel_diff_0885p00_zi, r_sqr_0885p00_zi,\
rmse_val_0885p00_zi_Venise,mean_abs_rel_diff_0885p00_zi_Venise, mean_rel_diff_0885p00_zi_Venise, r_sqr_0885p00_zi_Venise,\
rmse_val_0885p00_zi_Gloria,mean_abs_rel_diff_0885p00_zi_Gloria, mean_rel_diff_0885p00_zi_Gloria, r_sqr_0885p00_zi_Gloria,\
rmse_val_0885p00_zi_Galata_Platform,mean_abs_rel_diff_0885p00_zi_Galata_Platform, mean_rel_diff_0885p00_zi_Galata_Platform, r_sqr_0885p00_zi_Galata_Platform,\
rmse_val_0885p00_zi_Helsinki_Lighthouse,mean_abs_rel_diff_0885p00_zi_Helsinki_Lighthouse, mean_rel_diff_0885p00_zi_Helsinki_Lighthouse, r_sqr_0885p00_zi_Helsinki_Lighthouse,\
rmse_val_0885p00_zi_Gustav_Dalen_Tower,mean_abs_rel_diff_0885p00_zi_Gustav_Dalen_Tower, mean_rel_diff_0885p00_zi_Gustav_Dalen_Tower, r_sqr_0885p00_zi_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0885p00_ins_zi,matchups_PAN_rhow_0885p00_sat_zi,'885.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0885p00_ins_zi_station,min_val=-0.001,max_val=0.005)

#% plots  
prot_name = 'ba' 
sensor_name = 'OLCI'
rmse_val_0400p00_ba, mean_abs_rel_diff_0400p00_ba, mean_rel_diff_0400p00_ba, r_sqr_0400p00_ba,\
rmse_val_0400p00_ba_Venise,mean_abs_rel_diff_0400p00_ba_Venise, mean_rel_diff_0400p00_ba_Venise, r_sqr_0400p00_ba_Venise,\
rmse_val_0400p00_ba_Gloria,mean_abs_rel_diff_0400p00_ba_Gloria, mean_rel_diff_0400p00_ba_Gloria, r_sqr_0400p00_ba_Gloria,\
rmse_val_0400p00_ba_Galata_Platform,mean_abs_rel_diff_0400p00_ba_Galata_Platform, mean_rel_diff_0400p00_ba_Galata_Platform, r_sqr_0400p00_ba_Galata_Platform,\
rmse_val_0400p00_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0400p00_ba_Helsinki_Lighthouse, mean_rel_diff_0400p00_ba_Helsinki_Lighthouse, r_sqr_0400p00_ba_Helsinki_Lighthouse,\
rmse_val_0400p00_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0400p00_ba_Gustav_Dalen_Tower, mean_rel_diff_0400p00_ba_Gustav_Dalen_Tower, r_sqr_0400p00_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0400p00_ins_ba,matchups_PAN_rhow_0400p00_sat_ba,'400.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0400p00_ins_ba_station,min_val=-0.01,max_val=0.1)

rmse_val_0412p50_ba, mean_abs_rel_diff_0412p50_ba, mean_rel_diff_0412p50_ba, r_sqr_0412p50_ba,\
rmse_val_0412p50_ba_Venise,mean_abs_rel_diff_0412p50_ba_Venise, mean_rel_diff_0412p50_ba_Venise, r_sqr_0412p50_ba_Venise,\
rmse_val_0412p50_ba_Gloria,mean_abs_rel_diff_0412p50_ba_Gloria, mean_rel_diff_0412p50_ba_Gloria, r_sqr_0412p50_ba_Gloria,\
rmse_val_0412p50_ba_Galata_Platform,mean_abs_rel_diff_0412p50_ba_Galata_Platform, mean_rel_diff_0412p50_ba_Galata_Platform, r_sqr_0412p50_ba_Galata_Platform,\
rmse_val_0412p50_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0412p50_ba_Helsinki_Lighthouse, mean_rel_diff_0412p50_ba_Helsinki_Lighthouse, r_sqr_0412p50_ba_Helsinki_Lighthouse,\
rmse_val_0412p50_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0412p50_ba_Gustav_Dalen_Tower, mean_rel_diff_0412p50_ba_Gustav_Dalen_Tower, r_sqr_0412p50_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0412p50_fq_ins_ba,matchups_PAN_rhow_0412p50_fq_sat_ba,'412.5',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0412p50_fq_ins_ba_station,min_val=-0.01,max_val=0.1)

rmse_val_0442p50_ba, mean_abs_rel_diff_0442p50_ba, mean_rel_diff_0442p50_ba, r_sqr_0442p50_ba,\
rmse_val_0442p50_ba_Venise,mean_abs_rel_diff_0442p50_ba_Venise, mean_rel_diff_0442p50_ba_Venise, r_sqr_0442p50_ba_Venise,\
rmse_val_0442p50_ba_Gloria,mean_abs_rel_diff_0442p50_ba_Gloria, mean_rel_diff_0442p50_ba_Gloria, r_sqr_0442p50_ba_Gloria,\
rmse_val_0442p50_ba_Galata_Platform,mean_abs_rel_diff_0442p50_ba_Galata_Platform, mean_rel_diff_0442p50_ba_Galata_Platform, r_sqr_0442p50_ba_Galata_Platform,\
rmse_val_0442p50_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0442p50_ba_Helsinki_Lighthouse, mean_rel_diff_0442p50_ba_Helsinki_Lighthouse, r_sqr_0442p50_ba_Helsinki_Lighthouse,\
rmse_val_0442p50_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0442p50_ba_Gustav_Dalen_Tower, mean_rel_diff_0442p50_ba_Gustav_Dalen_Tower, r_sqr_0442p50_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0442p50_fq_ins_ba,matchups_PAN_rhow_0442p50_fq_sat_ba,'442.5',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0442p50_fq_ins_ba_station,min_val=-0.01,max_val=0.1)

rmse_val_0490p00_ba, mean_abs_rel_diff_0490p00_ba, mean_rel_diff_0490p00_ba, r_sqr_0490p00_ba,\
rmse_val_0490p00_ba_Venise,mean_abs_rel_diff_0490p00_ba_Venise, mean_rel_diff_0490p00_ba_Venise, r_sqr_0490p00_ba_Venise,\
rmse_val_0490p00_ba_Gloria,mean_abs_rel_diff_0490p00_ba_Gloria, mean_rel_diff_0490p00_ba_Gloria, r_sqr_0490p00_ba_Gloria,\
rmse_val_0490p00_ba_Galata_Platform,mean_abs_rel_diff_0490p00_ba_Galata_Platform, mean_rel_diff_0490p00_ba_Galata_Platform, r_sqr_0490p00_ba_Galata_Platform,\
rmse_val_0490p00_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0490p00_ba_Helsinki_Lighthouse, mean_rel_diff_0490p00_ba_Helsinki_Lighthouse, r_sqr_0490p00_ba_Helsinki_Lighthouse,\
rmse_val_0490p00_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0490p00_ba_Gustav_Dalen_Tower, mean_rel_diff_0490p00_ba_Gustav_Dalen_Tower, r_sqr_0490p00_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0490p00_fq_ins_ba,matchups_PAN_rhow_0490p00_fq_sat_ba,'490.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0490p00_fq_ins_ba_station,min_val=-0.01,max_val=0.12)

rmse_val_0510p00_ba, mean_abs_rel_diff_0510p00_ba, mean_rel_diff_0510p00_ba, r_sqr_0510p00_ba,\
rmse_val_0510p00_ba_Venise,mean_abs_rel_diff_0510p00_ba_Venise, mean_rel_diff_0510p00_ba_Venise, r_sqr_0510p00_ba_Venise,\
rmse_val_0510p00_ba_Gloria,mean_abs_rel_diff_0510p00_ba_Gloria, mean_rel_diff_0510p00_ba_Gloria, r_sqr_0510p00_ba_Gloria,\
rmse_val_0510p00_ba_Galata_Platform,mean_abs_rel_diff_0510p00_ba_Galata_Platform, mean_rel_diff_0510p00_ba_Galata_Platform, r_sqr_0510p00_ba_Galata_Platform,\
rmse_val_0510p00_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0510p00_ba_Helsinki_Lighthouse, mean_rel_diff_0510p00_ba_Helsinki_Lighthouse, r_sqr_0510p00_ba_Helsinki_Lighthouse,\
rmse_val_0510p00_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0510p00_ba_Gustav_Dalen_Tower, mean_rel_diff_0510p00_ba_Gustav_Dalen_Tower, r_sqr_0510p00_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0510p00_fq_ins_ba,matchups_PAN_rhow_0510p00_fq_sat_ba,'510.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0510p00_fq_ins_ba_station,min_val=-0.01,max_val=0.12)

rmse_val_0560p00_ba, mean_abs_rel_diff_0560p00_ba, mean_rel_diff_0560p00_ba, r_sqr_0560p00_ba,\
rmse_val_0560p00_ba_Venise,mean_abs_rel_diff_0560p00_ba_Venise, mean_rel_diff_0560p00_ba_Venise, r_sqr_0560p00_ba_Venise,\
rmse_val_0560p00_ba_Gloria,mean_abs_rel_diff_0560p00_ba_Gloria, mean_rel_diff_0560p00_ba_Gloria, r_sqr_0560p00_ba_Gloria,\
rmse_val_0560p00_ba_Galata_Platform,mean_abs_rel_diff_0560p00_ba_Galata_Platform, mean_rel_diff_0560p00_ba_Galata_Platform, r_sqr_0560p00_ba_Galata_Platform,\
rmse_val_0560p00_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0560p00_ba_Helsinki_Lighthouse, mean_rel_diff_0560p00_ba_Helsinki_Lighthouse, r_sqr_0560p00_ba_Helsinki_Lighthouse,\
rmse_val_0560p00_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0560p00_ba_Gustav_Dalen_Tower, mean_rel_diff_0560p00_ba_Gustav_Dalen_Tower, r_sqr_0560p00_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0560p00_fq_ins_ba,matchups_PAN_rhow_0560p00_fq_sat_ba,'560.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0560p00_fq_ins_ba_station,min_val=-0.01,max_val=0.12)

rmse_val_0620p00_ba, mean_abs_rel_diff_0620p00_ba, mean_rel_diff_0620p00_ba, r_sqr_0620p00_ba,\
rmse_val_0620p00_ba_Venise,mean_abs_rel_diff_0620p00_ba_Venise, mean_rel_diff_0620p00_ba_Venise, r_sqr_0620p00_ba_Venise,\
rmse_val_0620p00_ba_Gloria,mean_abs_rel_diff_0620p00_ba_Gloria, mean_rel_diff_0620p00_ba_Gloria, r_sqr_0620p00_ba_Gloria,\
rmse_val_0620p00_ba_Galata_Platform,mean_abs_rel_diff_0620p00_ba_Galata_Platform, mean_rel_diff_0620p00_ba_Galata_Platform, r_sqr_0620p00_ba_Galata_Platform,\
rmse_val_0620p00_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0620p00_ba_Helsinki_Lighthouse, mean_rel_diff_0620p00_ba_Helsinki_Lighthouse, r_sqr_0620p00_ba_Helsinki_Lighthouse,\
rmse_val_0620p00_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0620p00_ba_Gustav_Dalen_Tower, mean_rel_diff_0620p00_ba_Gustav_Dalen_Tower, r_sqr_0620p00_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0620p00_fq_ins_ba,matchups_PAN_rhow_0620p00_fq_sat_ba,'620.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0620p00_fq_ins_ba_station,min_val=-0.01,max_val=0.08)

rmse_val_0665p00_ba, mean_abs_rel_diff_0665p00_ba, mean_rel_diff_0665p00_ba, r_sqr_0665p00_ba,\
rmse_val_0665p00_ba_Venise,mean_abs_rel_diff_0665p00_ba_Venise, mean_rel_diff_0665p00_ba_Venise, r_sqr_0665p00_ba_Venise,\
rmse_val_0665p00_ba_Gloria,mean_abs_rel_diff_0665p00_ba_Gloria, mean_rel_diff_0665p00_ba_Gloria, r_sqr_0665p00_ba_Gloria,\
rmse_val_0665p00_ba_Galata_Platform,mean_abs_rel_diff_0665p00_ba_Galata_Platform, mean_rel_diff_0665p00_ba_Galata_Platform, r_sqr_0665p00_ba_Galata_Platform,\
rmse_val_0665p00_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0665p00_ba_Helsinki_Lighthouse, mean_rel_diff_0665p00_ba_Helsinki_Lighthouse, r_sqr_0665p00_ba_Helsinki_Lighthouse,\
rmse_val_0665p00_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0665p00_ba_Gustav_Dalen_Tower, mean_rel_diff_0665p00_ba_Gustav_Dalen_Tower, r_sqr_0665p00_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0665p00_fq_ins_ba,matchups_PAN_rhow_0665p00_fq_sat_ba,'665.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0665p00_fq_ins_ba_station,min_val=-0.01,max_val=0.06)

rmse_val_0673p75_ba, mean_abs_rel_diff_0673p75_ba, mean_rel_diff_0673p75_ba, r_sqr_0673p75_ba,\
rmse_val_0673p75_ba_Venise,mean_abs_rel_diff_0673p75_ba_Venise, mean_rel_diff_0673p75_ba_Venise, r_sqr_0673p75_ba_Venise,\
rmse_val_0673p75_ba_Gloria,mean_abs_rel_diff_0673p75_ba_Gloria, mean_rel_diff_0673p75_ba_Gloria, r_sqr_0673p75_ba_Gloria,\
rmse_val_0673p75_ba_Galata_Platform,mean_abs_rel_diff_0673p75_ba_Galata_Platform, mean_rel_diff_0673p75_ba_Galata_Platform, r_sqr_0673p75_ba_Galata_Platform,\
rmse_val_0673p75_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0673p75_ba_Helsinki_Lighthouse, mean_rel_diff_0673p75_ba_Helsinki_Lighthouse, r_sqr_0673p75_ba_Helsinki_Lighthouse,\
rmse_val_0673p75_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0673p75_ba_Gustav_Dalen_Tower, mean_rel_diff_0673p75_ba_Gustav_Dalen_Tower, r_sqr_0673p75_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0673p75_ins_ba,matchups_PAN_rhow_0673p75_sat_ba,'673.75',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0673p75_ins_ba_station,min_val=-0.01,max_val=0.06)

rmse_val_0681p25_ba, mean_abs_rel_diff_0681p25_ba, mean_rel_diff_0681p25_ba, r_sqr_0681p25_ba,\
rmse_val_0681p25_ba_Venise,mean_abs_rel_diff_0681p25_ba_Venise, mean_rel_diff_0681p25_ba_Venise, r_sqr_0681p25_ba_Venise,\
rmse_val_0681p25_ba_Gloria,mean_abs_rel_diff_0681p25_ba_Gloria, mean_rel_diff_0681p25_ba_Gloria, r_sqr_0681p25_ba_Gloria,\
rmse_val_0681p25_ba_Galata_Platform,mean_abs_rel_diff_0681p25_ba_Galata_Platform, mean_rel_diff_0681p25_ba_Galata_Platform, r_sqr_0681p25_ba_Galata_Platform,\
rmse_val_0681p25_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0681p25_ba_Helsinki_Lighthouse, mean_rel_diff_0681p25_ba_Helsinki_Lighthouse, r_sqr_0681p25_ba_Helsinki_Lighthouse,\
rmse_val_0681p25_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0681p25_ba_Gustav_Dalen_Tower, mean_rel_diff_0681p25_ba_Gustav_Dalen_Tower, r_sqr_0681p25_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0681p25_ins_ba,matchups_PAN_rhow_0681p25_sat_ba,'681.25',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0681p25_ins_ba_station,min_val=-0.01,max_val=0.06)

rmse_val_0708p75_ba, mean_abs_rel_diff_0708p75_ba, mean_rel_diff_0708p75_ba, r_sqr_0708p75_ba,\
rmse_val_0708p75_ba_Venise,mean_abs_rel_diff_0708p75_ba_Venise, mean_rel_diff_0708p75_ba_Venise, r_sqr_0708p75_ba_Venise,\
rmse_val_0708p75_ba_Gloria,mean_abs_rel_diff_0708p75_ba_Gloria, mean_rel_diff_0708p75_ba_Gloria, r_sqr_0708p75_ba_Gloria,\
rmse_val_0708p75_ba_Galata_Platform,mean_abs_rel_diff_0708p75_ba_Galata_Platform, mean_rel_diff_0708p75_ba_Galata_Platform, r_sqr_0708p75_ba_Galata_Platform,\
rmse_val_0708p75_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0708p75_ba_Helsinki_Lighthouse, mean_rel_diff_0708p75_ba_Helsinki_Lighthouse, r_sqr_0708p75_ba_Helsinki_Lighthouse,\
rmse_val_0708p75_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0708p75_ba_Gustav_Dalen_Tower, mean_rel_diff_0708p75_ba_Gustav_Dalen_Tower, r_sqr_0708p75_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0708p75_ins_ba,matchups_PAN_rhow_0708p75_sat_ba,'708.75',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0708p75_ins_ba_station,min_val=-0.01,max_val=0.04)

rmse_val_0753p75_ba, mean_abs_rel_diff_0753p75_ba, mean_rel_diff_0753p75_ba, r_sqr_0753p75_ba,\
rmse_val_0753p75_ba_Venise,mean_abs_rel_diff_0753p75_ba_Venise, mean_rel_diff_0753p75_ba_Venise, r_sqr_0753p75_ba_Venise,\
rmse_val_0753p75_ba_Gloria,mean_abs_rel_diff_0753p75_ba_Gloria, mean_rel_diff_0753p75_ba_Gloria, r_sqr_0753p75_ba_Gloria,\
rmse_val_0753p75_ba_Galata_Platform,mean_abs_rel_diff_0753p75_ba_Galata_Platform, mean_rel_diff_0753p75_ba_Galata_Platform, r_sqr_0753p75_ba_Galata_Platform,\
rmse_val_0753p75_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0753p75_ba_Helsinki_Lighthouse, mean_rel_diff_0753p75_ba_Helsinki_Lighthouse, r_sqr_0753p75_ba_Helsinki_Lighthouse,\
rmse_val_0753p75_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0753p75_ba_Gustav_Dalen_Tower, mean_rel_diff_0753p75_ba_Gustav_Dalen_Tower, r_sqr_0753p75_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0753p75_ins_ba,matchups_PAN_rhow_0753p75_sat_ba,'753.75',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0753p75_ins_ba_station,min_val=-0.001,max_val=0.009)

rmse_val_0778p75_ba, mean_abs_rel_diff_0778p75_ba, mean_rel_diff_0778p75_ba, r_sqr_0778p75_ba,\
rmse_val_0778p75_ba_Venise,mean_abs_rel_diff_0778p75_ba_Venise, mean_rel_diff_0778p75_ba_Venise, r_sqr_0778p75_ba_Venise,\
rmse_val_0778p75_ba_Gloria,mean_abs_rel_diff_0778p75_ba_Gloria, mean_rel_diff_0778p75_ba_Gloria, r_sqr_0778p75_ba_Gloria,\
rmse_val_0778p75_ba_Galata_Platform,mean_abs_rel_diff_0778p75_ba_Galata_Platform, mean_rel_diff_0778p75_ba_Galata_Platform, r_sqr_0778p75_ba_Galata_Platform,\
rmse_val_0778p75_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0778p75_ba_Helsinki_Lighthouse, mean_rel_diff_0778p75_ba_Helsinki_Lighthouse, r_sqr_0778p75_ba_Helsinki_Lighthouse,\
rmse_val_0778p75_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0778p75_ba_Gustav_Dalen_Tower, mean_rel_diff_0778p75_ba_Gustav_Dalen_Tower, r_sqr_0778p75_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0778p75_ins_ba,matchups_PAN_rhow_0778p75_sat_ba,'778.75',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0778p75_ins_ba_station,min_val=-0.001,max_val=0.009)

rmse_val_0865p00_ba, mean_abs_rel_diff_0865p00_ba, mean_rel_diff_0865p00_ba, r_sqr_0865p00_ba,\
rmse_val_0865p00_ba_Venise,mean_abs_rel_diff_0865p00_ba_Venise, mean_rel_diff_0865p00_ba_Venise, r_sqr_0865p00_ba_Venise,\
rmse_val_0865p00_ba_Gloria,mean_abs_rel_diff_0865p00_ba_Gloria, mean_rel_diff_0865p00_ba_Gloria, r_sqr_0865p00_ba_Gloria,\
rmse_val_0865p00_ba_Galata_Platform,mean_abs_rel_diff_0865p00_ba_Galata_Platform, mean_rel_diff_0865p00_ba_Galata_Platform, r_sqr_0865p00_ba_Galata_Platform,\
rmse_val_0865p00_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0865p00_ba_Helsinki_Lighthouse, mean_rel_diff_0865p00_ba_Helsinki_Lighthouse, r_sqr_0865p00_ba_Helsinki_Lighthouse,\
rmse_val_0865p00_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0865p00_ba_Gustav_Dalen_Tower, mean_rel_diff_0865p00_ba_Gustav_Dalen_Tower, r_sqr_0865p00_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0865p00_ins_ba,matchups_PAN_rhow_0865p00_sat_ba,'865.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0865p00_ins_ba_station,min_val=-0.001,max_val=0.006)

rmse_val_0885p00_ba, mean_abs_rel_diff_0885p00_ba, mean_rel_diff_0885p00_ba, r_sqr_0885p00_ba,\
rmse_val_0885p00_ba_Venise,mean_abs_rel_diff_0885p00_ba_Venise, mean_rel_diff_0885p00_ba_Venise, r_sqr_0885p00_ba_Venise,\
rmse_val_0885p00_ba_Gloria,mean_abs_rel_diff_0885p00_ba_Gloria, mean_rel_diff_0885p00_ba_Gloria, r_sqr_0885p00_ba_Gloria,\
rmse_val_0885p00_ba_Galata_Platform,mean_abs_rel_diff_0885p00_ba_Galata_Platform, mean_rel_diff_0885p00_ba_Galata_Platform, r_sqr_0885p00_ba_Galata_Platform,\
rmse_val_0885p00_ba_Helsinki_Lighthouse,mean_abs_rel_diff_0885p00_ba_Helsinki_Lighthouse, mean_rel_diff_0885p00_ba_Helsinki_Lighthouse, r_sqr_0885p00_ba_Helsinki_Lighthouse,\
rmse_val_0885p00_ba_Gustav_Dalen_Tower,mean_abs_rel_diff_0885p00_ba_Gustav_Dalen_Tower, mean_rel_diff_0885p00_ba_Gustav_Dalen_Tower, r_sqr_0885p00_ba_Gustav_Dalen_Tower\
= plot_scatter(\
    matchups_PAN_rhow_0885p00_ins_ba,matchups_PAN_rhow_0885p00_sat_ba,'885.0',path_out,prot_name,sensor_name,\
    matchups_PAN_rhow_0885p00_ins_ba_station,min_val=-0.001,max_val=0.005)

# #%%
# if __name__ == '__main__':
#     main()