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
import datetime
import time
import calendar
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
# import COMMON.apply_OLCI_flags as apply_OLCI_flags
import COMMON.Class_Flags_OLCI as flag

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
    Reads and checks configuration file for the validation
    Args:
        FILEconfig(str): configuration file path
    
    Return:
        options object
    """
    options = configparser.ConfigParser()
    options.read(FILEconfig)
    input_path = str(options['file_path']['input_directory'])
    output_path = str(options['file_path']['output_directory'])
    if input_path == '' or input_path == ' ' or output_path == '' or output_path == ' ':
        print ('please provide input and output directories path')
        sys.exit()
    if options['insitu_options']['sensor'] not in (['PANTHYR']):
        print ("please select a valid in situ type: 'PANTHYR'")
        sys.exit()
    if options['satellite_options']['platform'] not in (['A','B']):
        print ("please select a valid name for platform: A or B")
        sys.exit()
    if options['insitu_options']['sensor'] == 'AERONET' and options['Time_and_sites_selection']['sites'] in (['',' ']):
        options['Time_and_sites_selection']['sites'] = 'ALL'
    if options['Filtering_options']['flags'] in [' ','']:
        options['Filtering_options']['flags'] = 'None'
    if int(options['satellite_options']['window_size']) > int(options['satellite_options']['miniprods_size']) or int(options['satellite_options']['window_size']) % 2 == 0:
        print ('windows_size must be an odd number between 1 and %d' %int(options['satellite_options']['miniprods_size']))
        sys.exit()
    return options


class PANTHYR_class(object):
    def __init__(self,options):
        from MDB_builder.MDBs_config_class import MDBs_config
        config_MDB = MDBs_config()
        config_MDB.defaults_olci()
        config_MDB.defaults_PANTHYR()


        self.satellite_Nbands = config_MDB.n_sensor_bands #No. OLCI bands
        self.insitu_Nbands = config_MDB.Rrs_Nbands_PANTHYR
        self.RRS = np.ma.empty((0,self.satellite_Nbands))
        self.RRS_std = np.ma.empty((0,self.satellite_Nbands))
        self.RRS_panthyr = np.ma.empty((0,self.insitu_Nbands))
        self.RRS_panthyr_mask_shift = np.ma.empty((0,self.insitu_Nbands))
        self.panthyr_min = np.ma.empty((0,self.insitu_Nbands))
        self.panthyr_num = np.ma.empty((0,self.insitu_Nbands))
        self.panthyr_max = np.ma.empty((0,self.insitu_Nbands))
        self.date = np.ma.empty((0,1))
        self.site = np.ma.empty((0,1))
        self.pdu = np.ma.empty((0,1))
        self.Panthyrdate = np.ma.empty((0,1))
        # self.Aerolevel = np.ma.empty((0,1))

        self.MDB_PANTHYR_reader(options)


    def MDB_PANTHYR_reader(self,options):
        """

        Finds and reads AERONET MDBs and filters data as defined in the configuration file.
        Args:
            options (dict): options from config. file
                                                     
        Used for internal purpose and should not be called directly.

        """
        startDate = datetime.datetime.strptime(options['Time_and_sites_selection']['time_start'],'%Y-%m-%d')
        endDate = datetime.datetime.strptime(options['Time_and_sites_selection']['time_stop'],'%Y-%m-%d')
        startDate = calendar.timegm(startDate.timetuple())
        endDate = calendar.timegm(endDate.timetuple()) + 86400.0
        dim_window = int(options['satellite_options']['window_size'])
        dim_total = int(options['satellite_options']['miniprods_size'])
        central_pixel = int(np.floor(dim_total / 2))
        n_bands = self.satellite_Nbands
        delta_t = int(options['Filtering_options']['time_difference'])
        #walk trhough input directory looking for MDBs
        input_list = []
        site_list = options['Time_and_sites_selection']['sites']
        site_list = site_list.replace(" ", "")
        site_list = str.split(site_list,',')
        # files_shortlisted = [('MDB_%s%s_%s_L2_AERONET_%s.nc'%(str(options['satellite']).upper().replace(" ", ""),str(options['platform']).replace(" ", ""),str(options['sensor']).upper().replace(" ", ""),s)) for s in site_list]
        # for filename in os.listdir(options['input_directory']):
        #     if filename in files_shortlisted:            
        #         input_list.append(filename)

        path_to_source = options['file_path']['input_directory']
        output_directory = options['file_path']['output_directory']
        '''
        date_list.txt created as:
        % cat file_list_local.txt|cut -d _ -f4|sort|uniq>date_list.txt
        ''' 
        # PANTHYR Data
        insitu_sensor = options['insitu_options']['sensor']
        satellite_sensor = options['satellite_options']['satellite'] # A and B

        # create list of MDBs
        type_product = 'MDB'
        res = options['satellite_options']['resolution']
        wce = f'"{type_product}*{satellite_sensor}*{res}*{insitu_sensor}*.nc"' # wild card expression
        path_to_list = create_list_MDBs(path_to_source,output_directory,wce,type_product)

        with open(path_to_list) as file:
            for idx, line in enumerate(file):
                # print(idx)
                MDBfile_path = line[:-1]
                input_list.append(MDBfile_path)        
                
        #definition of flags
        flag_list = str(options['Filtering_options']['flags'])
        flag_list = flag_list.replace(" ", "")
        flag_list = str.split(flag_list,',')
        flagging = flag.Class_Flags_OLCI()
        #extracting data from MDB
        MDB_index = 0
        check = 0

        for MDBpath in input_list:
            # open each MDB
            print('------------------')
            print(MDBpath.split('/')[-1])
            nc = Dataset(os.path.join(MDBpath))
            
            #check if it is an AERONET MDB
            if nc.satellite+nc.platform == str(options['satellite_options']['satellite']).replace(" ", "").upper()+str(options['satellite_options']['platform']).replace(" ", ""):
                #import current MDB variables
                insitu_bands = nc.variables['insitu_bands'][:] # insitu_bands(insitu_bands)
                satellite_bands = nc.variables['satellite_bands'][:]

                self.insitu_bands = insitu_bands
                self.satellite_bands = satellite_bands
                nc.satellite_stop_time
                curr_satellite_date = nc.satellite_stop_time
#                 curr_level = nc.variables[options['insitu_prefix']+'level'][:]
                curr_insitu_date = nc.variables['insitu_time'][:] # insitu_time(insitu_id)
                curr_site = nc.insitu_site_name

                time_difference = nc.variables['time_difference'][:]
                index_time = np.argmin(np.abs(time_difference))

                OZA = nc.satellite_VZA_center_pixel
                SZA = nc.satellite_SZA_center_pixel
                # print(f'OZA: {OZA:.1f}; SZA: {SZA:.1f}')

                #flags mask
                
                satellite_WQSF = nc.variables['satellite_WQSF'][central_pixel-int(np.floor(dim_window/2)):central_pixel + int(np.floor(dim_window / 2)) + 1,central_pixel
                                                                        - int(np.floor(dim_window / 2)):central_pixel + int(np.floor(dim_window / 2)) + 1]
                if (str(flag_list[0])) != 'None':
                    mask = flagging.Mask(satellite_WQSF,flag_list)
                    mask [np.where(mask != 0)] = 1
                else: 
                    mask = np.full(satellite_WQSF.shape,0,dtype=int)
                #NTP = Nuber of Total non-LAND Pixels
                land = flagging.Mask(satellite_WQSF,(['LAND']))
                inland_w = flagging.Mask(satellite_WQSF,(['INLAND_WATER']))
                land[np.where(inland_w > 0)] = 0
                NTP = np.power(dim_window,2) - np.sum(land,axis=(0,1))  

                print(f'NTP: {NTP:.0f}')              

                #extracting date and level (masked through selected OLCI flags)

                if time_difference[index_time] < delta_t*60*60\
                    and OZA <= float(options['Filtering_options']['sensor_zenith_max'])\
                    and SZA <= float(options['Filtering_options']['sun_zenith_max']):
                    print(f'time difference: {time_difference[index_time]}, within delta_t: {delta_t}')

                    satellite_rhow = nc.variables['satellite_rhow'][:]
                    for sat_band_index in range(0,satellite_rhow.shape[0]):
                        wl = satellite_bands[sat_band_index]
                        ins_band_index = np.argmin(np.abs(wl-insitu_bands))
                        print(f'Closest in situ band to sat band {wl}: {insitu_bands[ins_band_index]} nm')
    
    
                        curr_rrs = satellite_rhow[sat_band_index,central_pixel - int(np.floor(dim_window / 2)):central_pixel + int(np.floor(dim_window / 2) + 1),central_pixel 
                                                                            - int(np.floor(dim_window / 2)):central_pixel + int(np.floor(dim_window / 2)) + 1]
                        curr_rrs_mean = curr_rrs.mean()
                        curr_rrs_stdv = curr_rrs.std()
    
                        
    
                        if float(satellite_bands[sat_band_index]) == 560:
                            curr_rrs_CV = curr_rrs_stdv/curr_rrs_mean
                            print(f'Extract mean: {curr_rrs_mean:.4f}; stdev: {curr_rrs_stdv:.4f}, CV: {100*curr_rrs_CV:.4f} %')
                        else: 
                            print(f'Extract mean: {curr_rrs_mean:.4f}; stdev: {curr_rrs_stdv:.4f}')

                # for var in nc.variables:
                #     if 'satellite_rhow' in var: # float satellite_rhow(satellite_bands, satellite_size_box_x, satellite_size_box_y)
                        
                #         var[band_index,:,:]
                #         print(var)

                #         band_index += 1
                
                #mask by selected sites and period
                # mask_site = np.ones(nc.variables[options['sat_prefix'] +'time'][:].shape)
                # for selected_site in site_list:
                #     mask_site [np.where(curr_site == selected_site)] = 0
                # if options['sites'] == 'ALL':
                #     mask_site = mask_site * 0
#                 mask_site [np.where(curr_satellite_date > endDate)] = 1
#                 mask_site [np.where(curr_satellite_date < startDate)] = 1
#                 #mask by quality level for AERONET data
#                 mask_level = np.zeros(curr_level.shape)
#                 mask_level[np.ma.where(curr_level < float(options['aeronet_level_min']))] = 1
#                 #QA_insitu = nc.variables[options['insitu_prefix']+'QA'][:]
#                 #if int(options['qa']) == 0:
#                 #    pass
#                 #else:
#                 #    QA_mask = QA_insitu >= float(options['qa'])
#                 #    mask_level[~QA_mask] = 1
#                 curr_site = np.ma.masked_array(curr_site,mask_site)
#                 curr_pdu = np.ma.masked_array(nc.variables[options['sat_prefix'] +'PDU'][:],mask_site)
#                 curr_satellite_date = np.ma.masked_array(curr_satellite_date,mask_site)
#                 curr_mb = np.ma.empty((nc.variables[options['sat_prefix'] +'time'].size,0))  
#                 curr_mb_std = np.ma.empty((nc.variables[options['sat_prefix'] +'time'].size,0))
#                 curr_mbAER = np.ma.empty((nc.variables[options['sat_prefix'] +'time'].size,0))
#                 curr_mbAERmin = np.ma.empty((nc.variables[options['sat_prefix'] +'time'].size,0))
#                 curr_mbAERmax = np.ma.empty((nc.variables[options['sat_prefix'] +'time'].size,0))
#                 curr_mbAERnum = np.ma.empty((nc.variables[options['sat_prefix'] +'time'].size,0))
#                 #time difference mask for AERONET data
#                 time_diff = nc.variables['time_difference'][:]
#                 mask_time = np.zeros(time_diff.shape)
#                 mask_time [np.ma.where(time_diff > delta_t)] = 1
#                 time_diff = np.ma.masked_array(time_diff,mask_time)
#                 #angles mask
#                 OZA = nc.variables[options['sat_prefix'] +'OZA'][:,central_pixel - int(np.floor(dim_window/2)) : central_pixel+ int(np.floor(dim_window / 2)) + 1,central_pixel
#                                                                  - int(np.floor(dim_window / 2)) : central_pixel + int(np.floor(dim_window / 2)) + 1]
#                 OZA = np.ma.masked_array(OZA,OZA > float(options['sensor_zenith_max']))
#                 SZA = nc.variables[options['sat_prefix'] +'SZA'][:,central_pixel-int(np.floor(dim_window/2)):central_pixel + int(np.floor(dim_window / 2)) + 1,central_pixel
#                                                                  - int(np.floor(dim_window / 2)) : central_pixel + int(np.floor(dim_window / 2)) + 1]
#                 SZA = np.ma.masked_array(SZA,SZA > float(options['sun_zenith_max']))
#                 #flags mask
#                 if (str(flag_list[0])) != 'None':
#                     mask = flagging.Mask(nc.variables['satellite_WQSF'][central_pixel-int(np.floor(dim_window/2)):central_pixel + int(np.floor(dim_window / 2)) + 1,central_pixel
#                                                                                      - int(np.floor(dim_window / 2)):central_pixel + int(np.floor(dim_window / 2)) + 1],(flag_list))
#                     mask [np.where(mask != 0)] = 1
#                 else: 
#                     mask = np.full(nc.variables['satellite_WQSF'][central_pixel-int(np.floor(dim_window/2)):central_pixel + int(np.floor(dim_window / 2)) + 1,central_pixel
#                                                                                - int(np.floor(dim_window / 2)):central_pixel + int(np.floor(dim_window / 2)) + 1].shape,0,dtype=int)
#                 #NTP = Nuber of Total non-LAND Pixels
#                 land = flagging.Mask(nc.variables['satellite_WQSF'][central_pixel-int(np.floor(dim_window/2)):central_pixel + int(np.floor(dim_window / 2)) + 1,central_pixel
#                                                                                  - int(np.floor(dim_window / 2)):central_pixel + int(np.floor(dim_window / 2)) + 1],(['LAND']))
#                 inland_w = flagging.Mask(nc.variables['satellite_WQSF'][central_pixel-int(np.floor(dim_window/2)):central_pixel + int(np.floor(dim_window / 2)) + 1,central_pixel
#                                                                                      - int(np.floor(dim_window / 2)):central_pixel + int(np.floor(dim_window / 2)) + 1],(['INLAND_WATER']))
#                 land[np.where(inland_w > 0)] = 0
#                 NTP = np.power(dim_window,2) - np.sum(land,axis=(1,2))
#                 #extracting date and level (masked through selected OLCI flags)
#                 band_index = 0
#                 cv_412 = np.ma.empty((curr_satellite_date.size,0))
#                 curr_insitu_date = np.ma.masked_array(curr_insitu_date,mask_time)
#                 curr_insitu_date = np.ma.masked_array(curr_insitu_date,mask_level)
#                 curr_level = np.ma.masked_array(curr_level,mask_level)
#                 curr_level = np.ma.masked_array(curr_level,mask_time)
#                 #take the first not filtered (by quality level) spectrum 
#                 validAERindex = np.ma.zeros(nc.dimensions[options['sat_prefix'] +'id'].size,dtype = int)
#                 validAeroDate = np.ma.empty(nc.dimensions[options['sat_prefix'] +'id'].size)
#                 validAeroLevel = np.ma.empty(nc.dimensions[options['sat_prefix'] +'id'].size)
#                 for pos in range (0, nc.dimensions[options['sat_prefix'] +'id'].size):
#                         try:
#                             validAERindex [pos] = int(np.min(np.where(curr_insitu_date[pos,:].mask == 0)))
#                         except ValueError:
#                             validAERindex [pos] = np.int(0)
#                         validAeroDate [pos] = curr_insitu_date[pos,validAERindex[pos]]
#                         validAeroLevel [pos] = curr_level[pos,validAERindex[pos]]

#                 #extracting Rrs from dim_window X dim_window area
#                 for var in nc.variables:
#                     if options['sat_prefix'] in var and ('Rrs' in var or options['sat_prefix'] +'T865' in var):
#                         curr_rrs = np.ma.masked_array(nc.variables[var][:,central_pixel - int(np.floor(dim_window / 2)):central_pixel + int(np.floor(dim_window / 2) + 1),central_pixel 
#                                                                         - int(np.floor(dim_window / 2)):central_pixel + int(np.floor(dim_window / 2)) + 1],np.ndarray.astype((mask),int))
#                         #retrieve Rrs values without BRDF applied
#                         if (options['brdf'] == 'False' or options['brdf'] == 'F' 
#                             or options['brdf'] == 'FALSE') and band_index < options['brdf_n_bands'] and var != options['sat_prefix'] +'T865':
#                             curr_brdf = nc.variables[options['sat_prefix'] +'BRDF'][:,central_pixel  - int(np.floor(dim_window/2)):central_pixel 
#                                                                  + int(np.floor(dim_window/2)+1),central_pixel  - int(np.floor(dim_window/2)):central_pixel 
#                                                                  + int(np.floor(dim_window/2))+1,band_index]
#                             curr_brdf[curr_brdf.mask] = 1
#                             curr_rrs = np.ma.divide(curr_rrs,curr_brdf)
#                         numberValid = np.ma.count(curr_rrs,axis = (1,2))
                        
#                         curr_rrs = np.ma.masked_array(curr_rrs,OZA.mask)
#                         curr_rrs = np.ma.masked_array(curr_rrs,SZA.mask)                        
#                         curr_mean = np.ma.mean(curr_rrs, axis = (1,2))
#                         curr_std = np.ma.std(curr_rrs, axis = (1,2))
#                         #calculate and filter by Mean+1.5*std and recalculate mean 
#                         unfilMean = np.ma.repeat(curr_mean[:,np.newaxis,np.newaxis],dim_window,axis=1)
#                         unfilMean = np.ma.repeat(unfilMean,dim_window,axis=2)
#                         unfilstd = np.ma.repeat(curr_std[:,np.newaxis,np.newaxis],dim_window,axis=1)
#                         unfilstd = np.ma.repeat(unfilstd,dim_window,axis=2)
#                         tresh_up = unfilMean + 1.5 * unfilstd
#                         tresh_bottom = unfilMean - 1.5 * unfilstd
#                         mask_outlier = np.ma.zeros(curr_rrs.shape)
#                         mask_outlier[np.where(curr_rrs - tresh_up > 0)] = 1 
#                         mask_outlier[np.where(tresh_bottom - curr_rrs > 0)] = 1
#                         if options['outliers'] == 'False' or options['outliers'] == 'F' or options['outliers'] == 'FALSE':
#                             mask_outlier = mask_outlier * 0
#                         curr_rrs_filtered = np.ma.masked_array(curr_rrs,mask_outlier)
#                         curr_mean = np.ma.median(curr_rrs_filtered, axis = (1,2))
#                         curr_mean = np.ma.masked_array(curr_mean,mask_site)
#                         curr_std = np.ma.std(curr_rrs_filtered, axis = (1,2))
#                         curr_std = np.ma.masked_array(curr_std,mask_site)
#                         #mask out whole matchup if #of valid pixel < defined trheshold * NTP
#                         curr_mean = np.ma.masked_array(curr_mean,numberValid < float(options['valid_min_pixel']) * NTP)
#                         curr_std = np.ma.masked_array(curr_std,numberValid < float(options['valid_min_pixel']) * NTP) 
#                         curr_cv = np.ma.divide(curr_std,curr_mean)
#                         if (band_index < n_bands and satellite_bands[band_index] == 560):
#                             #name is 412 but is done on 560
#                             cv_412 = np.ma.append(cv_412,curr_cv[:,np.newaxis],axis = 1)
#                         if var != options['sat_prefix'] +'T865':
#                             curr_mean = curr_mean.reshape(curr_mean.size,1)
#                             curr_std = curr_std.reshape(curr_std.size,1)
#                             curr_mb = np.ma.append(curr_mb,curr_mean,axis =1)
#                             curr_mb_std = np.ma.append(curr_mb_std,curr_std,axis =1)
#                             band_index += 1
#                     if options['insitu_prefix'] in var and 'Rrs' in var:
#                         if 'bands' in var or 'applied_shift' in var:
#                             pass
#                         else:
#                             #read aeronet values and filter for time_difference and quality level. 
#                             #Take the closest valid data using validAERindex 
#                             currAER = np.ma.masked_array(nc.variables[var][:],mask_time)
#                             curr_shift = np.ma.masked_array(nc.variables[var + '_applied_shift'][:]).astype(int)
#                             currAER = np.ma.masked_array(currAER,mask_level)
#                             if options['bandshifting'] in (['FALSE','False','F']) :
#                                 currAER = np.ma.masked_array(currAER,curr_shift)
#                             currAERmin = np.ma.min(currAER,axis = 1)
#                             currAERmax = np.ma.max(currAER,axis = 1)
#                             currAERnum = np.ma.count(currAER,axis = 1)
#                             validAER = np.ma.empty(currAER[:,0].shape)
#                             for pos in range (0, nc.dimensions[options['sat_prefix'] +'id'].size):
#                                 validAER [pos] = currAER[pos,validAERindex[pos]]
#                             curr_mbAER = np.ma.append(curr_mbAER,np.ma.reshape(validAER,(currAER[:,0].size,1)),axis = 1)
#                             curr_mbAERmin = np.ma.append(curr_mbAERmin,np.ma.reshape(currAERmin,(currAERmin.size,1)),axis = 1)
#                             curr_mbAERmax = np.ma.append(curr_mbAERmax,np.ma.reshape(currAERmax,(currAERmax.size,1)),axis = 1)
#                             curr_mbAERnum = np.ma.append(curr_mbAERnum,np.ma.reshape(currAERnum,(currAERnum.size,1)),axis = 1)
#                 #filter out whole matchups if coefficient of variation at 560 > threshold
#                 cv_412 = np.median(cv_412,axis = 1)
#                 cv_412 = np.abs(np.repeat(np.reshape(cv_412,(cv_412.size,1)),self.n_bands,axis = 1))
#                 cv_412 = np.ma.masked_greater(cv_412,float(options ['cv_max']))
#                 if options['cv_'].lower() == 'false':
#                     cv_412.mask = cv_412.mask * 0
#                 curr_mb = np.ma.masked_array(curr_mb,cv_412.mask)
#                 curr_mb = curr_mb.reshape(curr_mb.shape[0],n_bands)
#                 curr_mb_std = np.ma.masked_array(curr_mb_std,cv_412.mask)
#                 curr_mb_std = curr_mb_std.reshape(curr_mb_std.shape[0],n_bands)
#                 curr_mbAER = curr_mbAER.reshape(curr_mbAER.shape[0],insitu_bands.size)
#                 curr_mbAERmin = curr_mbAERmin.reshape(curr_mbAERmin.shape[0],insitu_bands.size)
#                 curr_mbAERmax = curr_mbAERmax.reshape(curr_mbAERmax.shape[0],insitu_bands.size)
#                 curr_mbAERnum = curr_mbAERnum.reshape(curr_mbAERnum.shape[0],insitu_bands.size)
#                 self.RRS = np.ma.append(self.RRS,curr_mb,axis = 0)
#                 self.RRS_std = np.ma.append(self.RRS_std,curr_mb_std,axis = 0) 
#                 self.RRS_aeronet = np.ma.append(self.RRS_aeronet,curr_mbAER,axis = 0)
#                 self.aeronet_max = np.ma.append(self.aeronet_max,curr_mbAERmax,axis = 0)
#                 self.aeronet_min = np.ma.append(self.aeronet_min,curr_mbAERmin,axis = 0)
#                 self.aeronet_num = np.ma.append(self.aeronet_num,curr_mbAERnum,axis = 0)
#                 curr_site = np.ma.masked_array(curr_site,curr_mb[:,0].mask)
#                 curr_satellite_date = np.ma.masked_array(curr_satellite_date,curr_mb[:,0].mask)
#                 validAeroDate = np.ma.masked_array(validAeroDate,curr_mb[:,0].mask)
#                 validAeroLevel = np.ma.masked_array(validAeroLevel,curr_mb[:,0].mask)
#                 self.date = np.ma.append(self.date,curr_satellite_date)
#                 self.site = np.ma.append(self.site,curr_site)
#                 self.pdu = np.ma.append(self.pdu,curr_pdu)
#                 self.Aerodate = np.ma.append(self.Aerodate,validAeroDate)
#                 self.Aerolevel = np.ma.append(self.Aerolevel,validAeroLevel)
#                 check = 1
#             else:
#                 pass
#         if check == 1:
#                 nc.close()
#                 return 1
#         else:
#              return 0         
# #%%                
# def main():
"""business logic for when running this module as the primary one!"""
print('Main Code!')

path_main = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/HYPERNETS_D7p2/MDB_py/'

# Read config. file
# config_file = file_config_parse.config_file
config_file = os.path.join(path_main,'MDB_reader','config_file_OLCI_PANTHYR.ini')

if os.path.isfile(config_file) == True:
    options = config_reader(config_file)
else:
    print(config_file + ' does not exist. Please provide a valid config file path')
    sys.exit()

if options['insitu_options']['sensor'] == 'PANTHYR':
    #read AERONET MDB
    dataTOplot = PANTHYR_class(options)
    # if dataTOplot.RRS.size == 0:
    #     print ('no PANTHYR MDBs found!')
    # elif np.sum(~dataTOplot.date.mask * ~dataTOplot.Aerodate.mask) == 0:
    #     print ('no valid matchups found')
    # else:
    #     #plot data and save stats
    #     [df_data,df_overall,header] = plot_matchups(dataTOplot,options)
    #     write_csv_stat(df_data,df_overall,header,options)   


#%%
# if __name__ == '__main__':
#     main()