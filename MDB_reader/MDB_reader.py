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

# for plotting
color_dict = dict({\
 '400.00':'LightBlue',\
 '412.50':'DeepSkyBlue',\
 '442.50':'DodgerBlue',\
 '490.00':'Blue',\
 '510.00':'ForestGreen',\
 '560.00':'Green',\
 '620.00':'LightCoral',\
 '665.00':'Red',\
 '673.75':'Crimson',\
 '681.25':'FireBrick',\
 '708.75':'Silver',\
 '753.75':'Gray',\
 '778.75':'DimGray',\
 '865.00':'SlateGray',\
 '885.00':'DarkSlateGray',\
'1020.50':'Black'})

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
        
        N_MUs = 0 # number of MUs

        r_s = central_pixel - int(np.floor(dim_window/2))       # starting row
        r_e = central_pixel + int(np.floor(dim_window/2)) + 1   # ending row
        c_s = central_pixel - int(np.floor(dim_window/2))       # starting col
        c_e = central_pixel + int(np.floor(dim_window/2)) + 1   # ending col

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
        satellite_sensor = options['satellite_options']['satellite'] # S3
        platform = options['satellite_options']['platform']

        # create list of MDBs
        type_product = 'MDB'
        res = options['satellite_options']['resolution']
        wce = f'"{type_product}*{satellite_sensor}{platform}*{res}*{insitu_sensor}*.nc"' # wild card expression
        path_to_list = create_list_MDBs(path_to_source,output_directory,wce,type_product)
        
        plt.figure()

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
        # MDB_index = 0
        # check = 0

        for MDBpath in input_list:
            # open each MDB
            print('------------------')
            print(MDBpath.split('/')[-1])
            check = 0
            try:
                nc = Dataset(os.path.join(MDBpath))
                
            # except:
            except Exception as e:
                print(f'Exception: {e}')
                pass    
            
            #check if it is an AERONET MDB
            if nc.satellite+nc.platform == str(options['satellite_options']['satellite']).replace(" ", "").upper()+str(options['satellite_options']['platform']).replace(" ", "")\
                and (options['satellite_options']['proc_version'] in nc.satellite_proc_version):
                #import current MDB variables
                insitu_bands = nc.variables['insitu_bands'][:] # insitu_bands(insitu_bands)
                satellite_bands = nc.variables['satellite_bands'][:]

                satellite_BRDF_bands = nc.variables['satellite_BRDF_bands'][:]
                satellite_BRDF_bands_list = list(satellite_BRDF_bands)

                self.insitu_bands = insitu_bands
                self.satellite_bands = satellite_bands
                nc.satellite_stop_time
                # curr_satellite_date = nc.satellite_stop_time
#                 curr_level = nc.variables[options['insitu_prefix']+'level'][:]
                # curr_insitu_date = nc.variables['insitu_time'][:] # insitu_time(insitu_id)
                curr_site = nc.insitu_site_name

                time_difference = nc.variables['time_difference'][:]
                ins_time_index = np.argmin(np.abs(time_difference))

                OZA = nc.satellite_VZA_center_pixel
                SZA = nc.satellite_SZA_center_pixel
                print(f'OZA: {OZA:.1f}; SZA: {SZA:.1f}')
                
                sat_proc_version_str = nc.satellite_proc_version

                #flags mask
                satellite_WQSF = nc.variables['satellite_WQSF'][r_s:r_e,c_s:c_e]
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

                if time_difference[ins_time_index] < delta_t*60*60\
                    and OZA <= float(options['Filtering_options']['sensor_zenith_max'])\
                    and SZA <= float(options['Filtering_options']['sun_zenith_max']):
                    print(f'time difference: {time_difference[ins_time_index]}, within delta_t: {delta_t}')

                    satellite_rhow = nc.variables['satellite_rhow'][:]
                    satellite_BRDF_rhow = nc.variables['satellite_BRDF_rhow'][:]
                    insitu_rrs = nc.variables['insitu_rhow'][:]
                    insitu_rrs = insitu_rrs/np.pi # transform from rhow to Rrs

                    curr_ins_rrs = []
                    curr_sat_rrs_mean = []
                    curr_bands = []
                    
                    print(f'size insitu_rrs: {insitu_rrs.shape}')

                    for sat_band_index in range(0,satellite_rhow.shape[0]):
                        wl = satellite_bands[sat_band_index]
                        if options['satellite_options']['BRDF'] == 'T' and (wl in satellite_BRDF_bands_list):
                            print(f'Band {wl} in BRDF bands.')
                            sat_BRDF_band_index = np.argmin(np.abs(wl-satellite_BRDF_bands))
                            # satellite extract
                            curr_sat_box = satellite_BRDF_rhow[sat_BRDF_band_index,r_s:r_e,c_s:c_e]
                        else:
                            # satellite extract
                            curr_sat_box = satellite_rhow[sat_band_index,r_s:r_e,c_s:c_e]

                        ins_band_index = np.argmin(np.abs(wl-insitu_bands))
                        print(f'Closest in situ band to sat band {wl}: {insitu_bands[ins_band_index]} nm')

                        # print(curr_sat_box)
                        curr_sat_box = np.ma.masked_where(curr_sat_box == 65535,curr_sat_box)
                        # print(curr_sat_box)
                        numberValid = np.ma.count(curr_sat_box)
                         

                        #calculate and filter by Mean+1.5*std and recalculate mean 
                        unfilMean = curr_sat_box.mean()
                        unfilstd = curr_sat_box.std()
                        tresh_up = unfilMean + 1.5 * unfilstd
                        tresh_bottom = unfilMean - 1.5 * unfilstd
                        mask_outlier = np.ma.zeros(curr_sat_box.shape)
                        mask_outlier[np.where(curr_sat_box - tresh_up > 0)] = 1 
                        mask_outlier[np.where(tresh_bottom - curr_sat_box > 0)] = 1
                        if options['Filtering_options']['outliers'] == 'False' or options['Filtering_options']['outliers'] == 'F' or options['Filtering_options']['outliers'] == 'FALSE':
                            mask_outlier = mask_outlier * 0
                        # print(mask_outlier)
                        curr_sat_box_filtered = np.ma.masked_array(curr_sat_box,mask_outlier)
                        curr_sat_box_mean = np.ma.mean(curr_sat_box_filtered)
                        curr_sat_box_std = np.ma.std(curr_sat_box_filtered)

                        #mask out whole matchup if #of valid pixel < defined trheshold * NTP
                        if numberValid >= (float(options['Filtering_options']['valid_min_pixel']) * NTP):
                            # in situ
                            curr_ins_rrs.append(insitu_rrs[ins_band_index,ins_time_index])
                            curr_sat_rrs_mean.append(curr_sat_box_mean/np.pi) # transform rhow to Rrs
                            curr_bands.append(wl)
                            print(f'in situ: {insitu_rrs[ins_band_index,ins_time_index]}; sat: {curr_sat_box_mean}')
                            
                            if wl == 560:
                                #name is 412 but is done on 560
                                curr_sat_box_cv_412 = curr_sat_box_std/curr_sat_box_mean
                                print(f'cv_412: {curr_sat_box_cv_412}')
                                check = 1
                                
                        else:
                            print(f'Not included: NTP= {NTP:.0f}; numberValid= {numberValid}')

                    if check and curr_sat_box_cv_412 <= float(options['Filtering_options']['cv_max']):
                        N_MUs += 1
                        for sat_band_index in range(len(curr_bands)):
                            if curr_bands[sat_band_index] in satellite_BRDF_bands_list and options['satellite_options']['BRDF'] == 'T':
                                mfc = 'Gray'
                                lw = 1.5
                            else:
                                mfc = None
                                lw = None
                            if options['plot_options']['to_plot'] == 'rhow':
                                plt.scatter(curr_ins_rrs[sat_band_index]*np.pi,curr_sat_rrs_mean[sat_band_index]*np.pi,\
                                    c=color_dict[f'{curr_bands[sat_band_index]:.2f}'],edgecolors=mfc,linewidths=lw)
                            
                            elif options['plot_options']['to_plot'] == 'Rrs':
                                plt.scatter(curr_ins_rrs[sat_band_index],curr_sat_rrs_mean[sat_band_index],\
                                    c=color_dict[f'{curr_bands[sat_band_index]:.2f}'],edgecolors=mfc,linewidths=lw)
                            


        plt.gca().set_aspect('equal', adjustable='box')

        if options['plot_options']['to_plot'] == 'rhow':
            plt.xlabel(r'PANTHYR $\rho_{W}$',fontsize=12)
            plt.ylabel(r'OLCI $\rho_{W}$',fontsize=12)    
        elif options['plot_options']['to_plot'] == 'Rrs':
            plt.xlabel(r'PANTHYR $R_{rs}$',fontsize=12)
            plt.ylabel(r'OLCI $R_{rs}$',fontsize=12)               


        plt.legend(['400.00', '412.50', '442.50', '490.00', '510.00', '560.00', '620.00', '665.00',\
                    '673.75', '681.25', '708.75', '753.75', '778.75', '865.00', '885.00','1020.50'],
                   loc='upper left',\
                   bbox_to_anchor=(1.001, 1))
        xmin, xmax = plt.gca().get_xlim()
        ymin, ymax = plt.gca().get_ylim()
        xmin = np.min([xmin,ymin])
        xmax = np.max([xmax,ymax])
        print(f'xmin={xmin}; xmax={xmax}')
        plt.plot([xmin,xmax],[xmin,xmax],'--k')

        ofname = f'{satellite_sensor}{platform}_{res}_{insitu_sensor}'
        plt.title(ofname+ f'; N = {N_MUs}; {sat_proc_version_str}')
        if options['satellite_options']['BRDF'] == 'T':
            brdf_str = '_BRDF'
        else:
            brdf_str = ''
        ofname = ofname +'_'+ sat_proc_version_str.replace(' ','_').replace('.','p')+brdf_str+ '.pdf'
        ofname = os.path.join(output_directory,ofname)
        if 'T' in options['plot_options']['save_plot']:
            plt.savefig(ofname,dpi=300)
        print(ofname)
        print(f'N match-ups: {N_MUs}')
      
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