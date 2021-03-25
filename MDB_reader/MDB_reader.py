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
import pandas as pd
from matplotlib import pyplot as plt
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
from scipy import stats

import configparser

# User Defined Functions
code_home = os.path.abspath('../')
sys.path.append(code_home)
import COMMON.common_functions as cfs
# import COMMON.apply_OLCI_flags as apply_OLCI_flags
import COMMON.Class_Flags_OLCI as flag

debug = True

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

plot_lims_Rrs = dict({\
 '400.00':[-0.002,0.015],\
 '412.50':[-0.002,0.018],\
 '442.50':[-0.002,0.030],\
 '490.00':[-0.002,0.030],\
 '510.00':[-0.002,0.030],\
 '560.00':[-0.002,0.030],\
 '620.00':[-0.002,0.015],\
 '665.00':[-0.002,0.01],\
 '673.75':[-0.002,0.01],\
 '681.25':[-0.002,0.01],\
 '708.75':[-0.002,0.006],\
 '753.75':[-0.0002,0.0016],\
 '778.75':[-0.0002,0.0016],\
 '865.00':[-0.0001,0.0010],\
 '885.00':[-0.0002,0.0008],\
'1020.50':[-0.0003,0.0007]})
    
plot_lims_LWN = dict({\
 '400.00':[-0.4,4.0],\
 '412.50':[-0.4,4.0],\
 '442.50':[-0.2,6.0],\
 '490.00':[-0.2,8.0],\
 '510.00':[-0.1,8.0],\
 '560.00':[-0.1,8.0],\
 '620.00':[-0.1,2.5],\
 '665.00':[-0.1,3.0],\
 '673.75':[-0.1,1.6],\
 '681.25':[-0.1,1.6],\
 '708.75':[-0.1,0.9],\
 '753.75':[-0.02,0.36],\
 '778.75':[-0.02,0.4],\
 '865.00':[-0.01,0.09],\
 '885.00':[-0.02,0.08],\
'1020.50':[-0.04,0.1]})    

olci_band_list  = ['400.00', '412.50', '442.50', '490.00', '510.00', '560.00', '620.00', '665.00',\
                    '673.75', '681.25', '708.75', '753.75', '778.75', '865.00', '885.00','1020.50']
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

def plot_scatter(x,y,sat_band,ins_band,path_out,prot_name,sensor_name,\
                 station,sat_proc_version_str,satellite_sensor,\
                 platform,res,insitu_sensor,brdf_str,options,df): 

    # replace nan in y (sat data)
    x = np.array(x)
    y = np.array(y)
    # station_vec = np.array(station_vec)

    x = x[~np.isnan(y)] # it is assumed that only sat data could be nan
    # station_vec = station_vec[~np.isnan(y)]
    y = y[~np.isnan(y)]

    if options['plot_options']['to_plot'] == 'LWN':
        x = x*cfs.get_F0(ins_band)
        y = y*cfs.get_F0(sat_band)   


    count_Venise = 0
    
    # if curr_bands[sat_band_index] in satellite_BRDF_bands_list and options['satellite_options']['BRDF'] == 'T':
    #     mfc = 'Gray'
    #     lw = 1.5
    # else:
    #     mfc = None
    #     lw = None
    
    plt.figure()
    #plt.errorbar(x, y, xerr=e_x, yerr=e_y, fmt='or')
    for cnt, line in enumerate(y):
        if station == 'Venise_PAN':
            mrk_color = 'r'
            count_Venise = count_Venise+1


        if prot_name == 'ba':
            mrk_style = 'o'
        elif  prot_name == 'zi':
            mrk_style = '+' 

        plt.plot(x[cnt], y[cnt],color=mrk_color,marker=mrk_style)
        
    # plot_lims_Rrs     
    if options['plot_options']['to_plot'] == 'Rrs':    
        plt.xlim(plot_lims_Rrs[f'{sat_band:0.2f}'])
        plt.ylim(plot_lims_Rrs[f'{sat_band:0.2f}'])
    elif options['plot_options']['to_plot'] == 'LWN':
        plt.xlim(plot_lims_LWN[f'{sat_band:0.2f}'])
        plt.ylim(plot_lims_LWN[f'{sat_band:0.2f}'])
        
    plt.gca().set_aspect('equal', adjustable='box')
    # plot 1:1 line
    xmin, xmax = plt.gca().get_xlim()
    ymin, ymax = plt.gca().get_ylim()
    plt.plot([xmin,xmax],[ymin, ymax],'--k')
    
    # Generated linear fit
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    # line = slope*np.array([xmin,xmax],dtype=np.float32)+intercept
    # plt.plot([xmin,xmax], line)
    # plt.legend(['1:1','Regression Line'])
    
    if options['plot_options']['to_plot'] == 'rhow':
        plt.xlabel(f'{insitu_sensor} '+r'$\rho_{W}$',fontsize=12)
        plt.ylabel(r'OLCI $\rho_{W}$',fontsize=12)    
    elif options['plot_options']['to_plot'] == 'Rrs':
        plt.xlabel(f'{insitu_sensor} '+r'$R_{rs}$',fontsize=12)
        plt.ylabel(r'OLCI $R_{rs}$',fontsize=12) 
    elif options['plot_options']['to_plot'] == 'LWN':
        plt.xlabel(f'{insitu_sensor} '+r'$L_{WN}$'+f'({str(ins_band)})',fontsize=12)
        plt.ylabel(r'OLCI $L_{WN}$'+f'({str(sat_band)})',fontsize=12)     
    
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

    sat_band_str = str(sat_band)
    str2 = sat_band_str
    # # to print without .0
    # if str1[-2:]=='.0':
    #     str2 = str2[:-2]
        
    
    str0 = '{}\nN={:d}\nRMSD={:,.4f}\nMAPD={:,.0f}%\nMPD={:,.0f}%\n$r^2$={:,.2f}'\
    .format(sat_band_str,\
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

    ofname = f'{satellite_sensor}{platform}_{res}_{insitu_sensor}'
    plt.title(ofname+ f'; {sat_proc_version_str}')  
    
    ofname = f'{satellite_sensor}{platform}_{res}_{insitu_sensor}_{sat_proc_version_str[-5:].replace(".","p")}_{sat_band_str.replace(".","p")}.pdf'
    ofname = os.path.join(path_out,ofname)
    
    plt.savefig(ofname, dpi=300)

    # latex table
    if sat_band_str == '412.5':
        print('proto & nm & N & RMSD & MAPD & MPD & $r^2$\n')
    str_table = '{} & {} & {:d} & {:,.4f} & {:,.1f} & {:,.1f} & {:,.2f}\\\\'\
    .format(prot_name_str,\
            str2,\
            N,\
            rmse_val,\
            mean_abs_rel_diff,\
            mean_rel_diff,\
            r_value**2)

    print(str_table)
    df = df.append({'proc':sat_proc_version_str,'sat':satellite_sensor,\
                 'platform':platform,'wl':sat_band_str,'N':N,\
                    'rmse':rmse_val,'MAPD':mean_abs_rel_diff,'MPD':mean_rel_diff,\
                    'r^2':r_value**2},ignore_index=True) 
    
 
    print('count_Venise: '+str(count_Venise))
    
    return df
    # print('count_Gloria: '+str(count_Gloria))
    # print('count_Galata_Platform: '+str(count_Galata_Platform))
    # print('count_Helsinki_Lighthouse: '+str(count_Helsinki_Lighthouse))
    # print('count_Gustav_Dalen_Tower: '+str(count_Gustav_Dalen_Tower))

    # plt.show()   
    # return rmse_val, mean_abs_rel_diff, mean_rel_diff, r_value**2,\
    #     rmse_val_Venise, mean_abs_rel_diff_Venise, mean_rel_diff_Venise, r_value_Venise**2,\
    #     rmse_val_Gloria, mean_abs_rel_diff_Gloria, mean_rel_diff_Gloria, r_value_Gloria**2,\
    #     rmse_val_Galata_Platform, mean_abs_rel_diff_Galata_Platform, mean_rel_diff_Galata_Platform, r_value_Galata_Platform**2,\
    #     rmse_val_Helsinki_Lighthouse, mean_abs_rel_diff_Helsinki_Lighthouse, mean_rel_diff_Helsinki_Lighthouse, r_value_Helsinki_Lighthouse**2,\
        # rmse_val_Gustav_Dalen_Tower, mean_abs_rel_diff_Gustav_Dalen_Tower, mean_rel_diff_Gustav_Dalen_Tower, r_value_Gustav_Dalen_Tower**2
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
    if options['insitu_options']['sensor'] not in (['PANTHYR','HYPERNETS']):
        print ("please select a valid in situ type: PANTHYR or HYPERNETS ")
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
        
        columns = ['datetime','PDU','insitu_filename','OZA','SZA','bands','MU','outlier','version','ws',\
            'LWN_ins_400.00','LWN_ins_412.50','LWN_ins_442.50','LWN_ins_490.00','LWN_ins_510.00','LWN_ins_560.00','LWN_ins_620.00',\
            'LWN_ins_665.00','LWN_ins_673.75','LWN_ins_681.25','LWN_ins_708.75','LWN_ins_753.75','LWN_ins_778.75','LWN_ins_865.00',\
            'LWN_ins_885.00','LWN_sat_400.00','LWN_sat_412.50','LWN_sat_442.50','LWN_sat_490.00','LWN_sat_510.00',\
            'LWN_sat_560.00','LWN_sat_620.00','LWN_sat_665.00','LWN_sat_673.75','LWN_sat_681.25','LWN_sat_708.75','LWN_sat_753.75',\
            'LWN_sat_778.75','LWN_sat_865.00','LWN_sat_885.00','LWN_sat_1020.50',
            'Rrs_ins_400.00','Rrs_ins_412.50','Rrs_ins_442.50','Rrs_ins_490.00','Rrs_ins_510.00','Rrs_ins_560.00','Rrs_ins_620.00',\
            'Rrs_ins_665.00','Rrs_ins_673.75','Rrs_ins_681.25','Rrs_ins_708.75','Rrs_ins_753.75','Rrs_ins_778.75','Rrs_ins_865.00',\
            'Rrs_ins_885.00','Rrs_ins_1020.50','Rrs_sat_400.00','Rrs_sat_412.50','Rrs_sat_442.50','Rrs_sat_490.00','Rrs_sat_510.00',\
            'Rrs_sat_560.00','Rrs_sat_620.00','Rrs_sat_665.00','Rrs_sat_673.75','Rrs_sat_681.25','Rrs_sat_708.75','Rrs_sat_753.75',\
            'Rrs_sat_778.75','Rrs_sat_865.00','Rrs_sat_885.00','Rrs_sat_1020.50']
        df_matchups = pd.DataFrame(columns=columns)
        
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
        # In situ data
        insitu_sensor = options['insitu_options']['sensor']
        satellite_sensor = options['satellite_options']['satellite'] # S3
        platform = options['satellite_options']['platform']

        # create list of MDBs
        type_product = 'MDB'
        res = options['satellite_options']['resolution']
        wce = f'"{type_product}*{satellite_sensor}{platform}*{res}*{insitu_sensor}*.nc"' # wild card expression
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
        #extracting data from MDB
        # MDB_index = 0
        # check = 0
        
        matchups_sat_array = np.empty([0,len(olci_band_list)]) # dim: MU number, band number
        matchups_ins_array = np.empty([0,len(olci_band_list)])
        
        columns = ['version','PDU','ins: LWN(753.75)','sat: LWN(753.75)']
        df_outliers = pd.DataFrame(columns=columns)
        
        date_list = []
        outliers_ins_753p75 = []
        outliers_ins_753p75_values = []
        outliers_sat_753p75_values = []
        
        
        for MDBpath in input_list:
            # open each MDB
            print('------------------')
            print(MDBpath.split('/')[-1])
            check_cv560 = False
            try:
                nc = Dataset(os.path.join(MDBpath))
                
            # except:
            except Exception as e:
                print(f'Exception: {e}')
                pass    
            
            
            #check if it is an AERONET MDB
            if nc.satellite+nc.platform == str(options['satellite_options']['satellite']).replace(" ", "").upper()+str(options['satellite_options']['platform']).replace(" ", "")\
                and (options['satellite_options']['proc_version'] in nc.satellite_proc_version)\
                    and MDBpath.split('_')[-4][:8] not in date_list: # check repeated
                #import current MDB variables
                insitu_bands = nc.variables['insitu_bands'][:] # insitu_bands(insitu_bands)
                satellite_bands = nc.variables['satellite_bands'][:]

                satellite_BRDF_bands = nc.variables['satellite_BRDF_bands'][:]
                satellite_BRDF_bands_list = list(satellite_BRDF_bands)

                self.insitu_bands = insitu_bands
                self.satellite_bands = satellite_bands
                curr_satellite_date = nc.satellite_stop_time
#                 curr_level = nc.variables[options['insitu_prefix']+'level'][:]
                # curr_insitu_date = nc.variables['insitu_time'][:] # insitu_time(insitu_id)
                curr_site = nc.insitu_site_name

                time_difference = nc.variables['time_difference'][:]
                ins_time_index = np.argmin(np.abs(time_difference))

                OZA = nc.satellite_VZA_center_pixel
                SZA = nc.satellite_SZA_center_pixel
                ws0 = nc.satellite_ws0
                ws1 = nc.satellite_ws1
                ws = np.sqrt((float(ws0)**2)+(float(ws1)**2))

                
                sat_proc_version_str = nc.satellite_proc_version

                #flags mask
                flagging = flag.Class_Flags_OLCI(nc.satellite_WQSF_flag_masks,nc.satellite_WQSF_flag_meanings)
                satellite_WQSF = nc.variables['satellite_WQSF'][r_s:r_e,c_s:c_e]
                if (str(flag_list[0])) != 'None':
                    mask = flagging.Mask(satellite_WQSF,flag_list)
                    mask [np.where(mask != 0)] = 1
                else: 
                    mask = np.full(satellite_WQSF.shape,0,dtype=int)
                print(mask)
                NGP = np.power(dim_window,2) - np.sum(mask) # Number Good Pixels excluding flagged pixels
                #NTP = Number Total non-LAND Pixels
                land = flagging.Mask(satellite_WQSF,(['LAND']))
                inland_w = flagging.Mask(satellite_WQSF,(['INLAND_WATER']))
                land[np.where(inland_w > 0)] = 0
                NTP = np.power(dim_window,2) - np.sum(land,axis=(0,1))

                # conditions minimum valid pixels
                if float(options['Filtering_options']['valid_min_pixel']) == 1: # 100% like Zibordi
                    cond_min_valid_pxs = NGP == np.power(dim_window,2)
                else:
                    cond_min_valid_pxs = NGP > (float(options['Filtering_options']['valid_min_pixel']) * NTP + 1)

                if debug:
                    print(f'OZA: {OZA:.1f}; SZA: {SZA:.1f}; ws: {ws:.1f}; NGP: {NGP}')
                if time_difference[ins_time_index] < delta_t*60*60\
                    and OZA <= float(options['Filtering_options']['sensor_zenith_max'])\
                    and SZA <= float(options['Filtering_options']['sun_zenith_max'])\
                    and cond_min_valid_pxs:
                    # print(f'time difference: {time_difference[ins_time_index]}, within delta_t: {delta_t}')
                    
                    date_list.append(MDBpath.split('_')[-4][:8]) # to avoid multiple MUs per day
                    
                    satellite_rhow = nc.variables['satellite_rhow'][:]
                    satellite_BRDF_rhow = nc.variables['satellite_BRDF_rhow'][:]
                    insitu_rrs = nc.variables['insitu_rhow'][:]
                    insitu_rrs = insitu_rrs/np.pi # transform from rhow to Rrs
                    insitu_finame_list = nc.variables['insitu_filename'][:]

                    curr_ins_rrs = []
                    curr_sat_rrs_mean = []
                    curr_bands = []
                    
                    # print(f'size insitu_rrs: {insitu_rrs.shape}')
                    outlier_flag = False
                    for sat_band_index in range(0,satellite_rhow.shape[0]):
                        wl = satellite_bands[sat_band_index]
                        if options['satellite_options']['BRDF'] == 'T' and (wl in satellite_BRDF_bands_list):
                            if debug:
                                print(f'Band {wl} in BRDF bands.')
                            sat_BRDF_band_index = np.argmin(np.abs(wl-satellite_BRDF_bands))
                            # satellite extract
                            curr_sat_box = satellite_BRDF_rhow[sat_BRDF_band_index,r_s:r_e,c_s:c_e]
                        else:
                            # satellite extract
                            curr_sat_box = satellite_rhow[sat_band_index,r_s:r_e,c_s:c_e]

                        ins_band_index = np.argmin(np.abs(wl-insitu_bands))
                        if debug:
                            print(f'Closest in situ band to sat band {wl}: {insitu_bands[ins_band_index]} nm')

                        curr_sat_box = np.ma.masked_where(curr_sat_box == 65535,curr_sat_box)
                        curr_sat_box = np.ma.masked_invalid(curr_sat_box)
                        curr_sat_box = np.ma.masked_array(curr_sat_box,mask)

                        numberValid = np.ma.count(curr_sat_box)

                        # conditions minimum valid pixels. Note: condition repeated because some masking (e.g. invalid pixels) was made
                        # after the first time this condition was asked
                        if float(options['Filtering_options']['valid_min_pixel']) == 1: # 100% like Zibordi
                            cond_min_valid_pxs = numberValid == np.power(dim_window,2)
                        else:
                            cond_min_valid_pxs = numberValid > (float(options['Filtering_options']['valid_min_pixel']) * NTP + 1)

                        #mask out whole matchup if #of valid pixel < defined trheshold * NTP
                        if cond_min_valid_pxs:

                            #calculate and filter by Mean+1.5*std and recalculate mean 
                            unfilMean = curr_sat_box.mean()
                            unfilstd = curr_sat_box.std()
                            tresh_up = unfilMean + 1.5 * unfilstd
                            tresh_bottom = unfilMean - 1.5 * unfilstd
                            mask_outlier = np.ma.zeros(curr_sat_box.shape)
                            mask_outlier[np.where(curr_sat_box > tresh_up)] = 1 
                            mask_outlier[np.where(curr_sat_box < tresh_bottom)] = 1
                            if options['Filtering_options']['outliers'] in ['False','F','FALSE']:
                                mask_outlier = mask_outlier * 0
                            curr_sat_box_filtered = np.ma.masked_array(curr_sat_box,mask_outlier)
    
                            curr_sat_box_mean = np.ma.mean(curr_sat_box_filtered)
                            curr_sat_box_std = np.ma.std(curr_sat_box_filtered)


                            # in situ. curr_ins_rrs is a list with N=number of bands elements
                            ins_value = insitu_rrs[ins_band_index,ins_time_index]
                            insitu_filename = insitu_finame_list[ins_time_index]
                            curr_ins_rrs.append(ins_value)
                            curr_sat_rrs_mean.append(curr_sat_box_mean/np.pi) # transform rhow to Rrs
                            curr_bands.append(wl)
                            if debug:
                                print(f'in situ: {insitu_rrs[ins_band_index,ins_time_index]}; sat: {curr_sat_box_mean}')
                            
                            if wl == 560:
                                #name is 412 but is done on 560
                                curr_sat_box_cv_560 = curr_sat_box_std/curr_sat_box_mean
                                print(f'cv_560: {curr_sat_box_cv_560}')
                                check_cv560 = True
                            
                            # to detect outliers from the in situ data
                            # if wl == 753.75 \
                            #     and cfs.get_F0(753.75)*ins_value>0.05 \
                            #         and curr_sat_box_cv_560 <= float(options['Filtering_options']['cv_max']): # a MU
                            #     print('in situ 753.75 band > 0.05')
                            #     print(nc.satellite_PDU)
                            #     outlier_flag = True
                            #     print(outlier_flag)
                            #     outliers_ins_753p75.append(nc.satellite_PDU)
                            #     outliers_ins_753p75_values.append(cfs.get_F0(753.75)*ins_value)
                            #     outliers_sat_753p75_values.append(cfs.get_F0(753.75)*curr_sat_box_mean/np.pi)
                            #     df_outliers =  df_outliers.append({'version':sat_proc_version_str,\
                            #                         'PDU':nc.satellite_PDU,\
                            #                         'ins: LWN(753.75)':cfs.get_F0(753.75)*ins_value,\
                            #                         'sat: LWN(753.75)':cfs.get_F0(753.75)*curr_sat_box_mean/np.pi},ignore_index=True)
                                
                        else:
                            print(f'Not included: NTP= {NTP:.0f}; NGP={NGP}; numberValid= {numberValid}')
                    
            # scatter plots  
            MU_flag = False
            if check_cv560 and curr_sat_box_cv_560 <= float(options['Filtering_options']['cv_max']):
                N_MUs += 1
                MU_flag = True
                matchups_ins_array = np.append(matchups_ins_array,[curr_ins_rrs],axis=0)
                matchups_sat_array = np.append(matchups_sat_array,[curr_sat_rrs_mean],axis=0)
                # plotting all bands
                for sat_band_index in range(len(curr_bands)-1):
                    if curr_bands[sat_band_index] in satellite_BRDF_bands_list and options['satellite_options']['BRDF'] == 'T':
                        mfc = 'Gray'
                        lw = 1.5
                    else:
                        mfc = None
                        lw = None
                        
                    ins_band_index = np.argmin(np.abs(curr_bands[sat_band_index]-insitu_bands))
                    ins_band = insitu_bands[ins_band_index]
                    sat_band = curr_bands[sat_band_index]
                    if options['plot_options']['to_plot'] == 'rhow':
                        plt.scatter(curr_ins_rrs[sat_band_index]*np.pi,curr_sat_rrs_mean[sat_band_index]*np.pi,\
                            c=color_dict[f'{curr_bands[sat_band_index]:.2f}'],edgecolors=mfc,linewidths=lw)
                    elif options['plot_options']['to_plot'] == 'Rrs':
                        plt.scatter(curr_ins_rrs[sat_band_index],curr_sat_rrs_mean[sat_band_index],\
                            c=color_dict[f'{curr_bands[sat_band_index]:.2f}'],edgecolors=mfc,linewidths=lw)
                    elif options['plot_options']['to_plot'] == 'LWN':
                        plt.scatter(curr_ins_rrs[sat_band_index]*cfs.get_F0(ins_band),\
                            curr_sat_rrs_mean[sat_band_index]*cfs.get_F0(sat_band),\
                            c=color_dict[f'{curr_bands[sat_band_index]:.2f}'],edgecolors=mfc,linewidths=lw)

                df_matchups = df_matchups.append(\
                            {'datetime':curr_satellite_date,\
                            'PDU':nc.satellite_PDU,\
                            'insitu_filename': insitu_filename,\
                            'OZA':OZA,'SZA':SZA,\
                            'LWN_ins_400.00':cfs.get_F0(400.00)*curr_ins_rrs[0],\
                            'LWN_ins_412.50':cfs.get_F0(412.50)*curr_ins_rrs[1],\
                            'LWN_ins_442.50':cfs.get_F0(442.50)*curr_ins_rrs[2],\
                            'LWN_ins_490.00':cfs.get_F0(490.00)*curr_ins_rrs[3],\
                            'LWN_ins_510.00':cfs.get_F0(510.00)*curr_ins_rrs[4],\
                            'LWN_ins_560.00':cfs.get_F0(560.00)*curr_ins_rrs[5],\
                            'LWN_ins_620.00':cfs.get_F0(620.00)*curr_ins_rrs[6],\
                            'LWN_ins_665.00':cfs.get_F0(665.00)*curr_ins_rrs[7],\
                            'LWN_ins_673.75':cfs.get_F0(673.75)*curr_ins_rrs[8],\
                            'LWN_ins_681.25':cfs.get_F0(681.25)*curr_ins_rrs[9],\
                            'LWN_ins_708.75':cfs.get_F0(708.75)*curr_ins_rrs[10],\
                            'LWN_ins_753.75':cfs.get_F0(753.75)*curr_ins_rrs[11],\
                            'LWN_ins_778.75':cfs.get_F0(778.75)*curr_ins_rrs[12],\
                            'LWN_ins_865.00':cfs.get_F0(865.00)*curr_ins_rrs[13],\
                            'LWN_ins_885.00':cfs.get_F0(885.00)*curr_ins_rrs[14],\
                            'LWN_ins_1020.50':cfs.get_F0(1020.50)*curr_ins_rrs[15],\
                            'LWN_sat_400.00':cfs.get_F0(400.00)*curr_sat_rrs_mean[0],\
                            'LWN_sat_412.50':cfs.get_F0(412.50)*curr_sat_rrs_mean[1],\
                            'LWN_sat_442.50':cfs.get_F0(442.50)*curr_sat_rrs_mean[2],\
                            'LWN_sat_490.00':cfs.get_F0(490.00)*curr_sat_rrs_mean[3],\
                            'LWN_sat_510.00':cfs.get_F0(510.00)*curr_sat_rrs_mean[4],\
                            'LWN_sat_560.00':cfs.get_F0(560.00)*curr_sat_rrs_mean[5],\
                            'LWN_sat_620.00':cfs.get_F0(620.00)*curr_sat_rrs_mean[6],\
                            'LWN_sat_665.00':cfs.get_F0(665.00)*curr_sat_rrs_mean[7],\
                            'LWN_sat_673.75':cfs.get_F0(673.75)*curr_sat_rrs_mean[8],\
                            'LWN_sat_681.25':cfs.get_F0(681.25)*curr_sat_rrs_mean[9],\
                            'LWN_sat_708.75':cfs.get_F0(708.75)*curr_sat_rrs_mean[10],\
                            'LWN_sat_753.75':cfs.get_F0(753.75)*curr_sat_rrs_mean[11],\
                            'LWN_sat_778.75':cfs.get_F0(778.75)*curr_sat_rrs_mean[12],\
                            'LWN_sat_865.00':cfs.get_F0(865.00)*curr_sat_rrs_mean[13],\
                            'LWN_sat_885.00':cfs.get_F0(885.00)*curr_sat_rrs_mean[14],\
                            'LWN_sat_1020.50':cfs.get_F0(1020.50)*curr_sat_rrs_mean[15],\
                            'Rrs_ins_400.00':curr_ins_rrs[0],\
                            'Rrs_ins_412.50':curr_ins_rrs[1],\
                            'Rrs_ins_442.50':curr_ins_rrs[2],\
                            'Rrs_ins_490.00':curr_ins_rrs[3],\
                            'Rrs_ins_510.00':curr_ins_rrs[4],\
                            'Rrs_ins_560.00':curr_ins_rrs[5],\
                            'Rrs_ins_620.00':curr_ins_rrs[6],\
                            'Rrs_ins_665.00':curr_ins_rrs[7],\
                            'Rrs_ins_673.75':curr_ins_rrs[8],\
                            'Rrs_ins_681.25':curr_ins_rrs[9],\
                            'Rrs_ins_708.75':curr_ins_rrs[10],\
                            'Rrs_ins_753.75':curr_ins_rrs[11],\
                            'Rrs_ins_778.75':curr_ins_rrs[12],\
                            'Rrs_ins_865.00':curr_ins_rrs[13],\
                            'Rrs_ins_885.00':curr_ins_rrs[14],\
                            'Rrs_ins_1020.50':curr_ins_rrs[15],\
                            'Rrs_sat_400.00':curr_sat_rrs_mean[0],\
                            'Rrs_sat_412.50':curr_sat_rrs_mean[1],\
                            'Rrs_sat_442.50':curr_sat_rrs_mean[2],\
                            'Rrs_sat_490.00':curr_sat_rrs_mean[3],\
                            'Rrs_sat_510.00':curr_sat_rrs_mean[4],\
                            'Rrs_sat_560.00':curr_sat_rrs_mean[5],\
                            'Rrs_sat_620.00':curr_sat_rrs_mean[6],\
                            'Rrs_sat_665.00':curr_sat_rrs_mean[7],\
                            'Rrs_sat_673.75':curr_sat_rrs_mean[8],\
                            'Rrs_sat_681.25':curr_sat_rrs_mean[9],\
                            'Rrs_sat_708.75':curr_sat_rrs_mean[10],\
                            'Rrs_sat_753.75':curr_sat_rrs_mean[11],\
                            'Rrs_sat_778.75':curr_sat_rrs_mean[12],\
                            'Rrs_sat_865.00':curr_sat_rrs_mean[13],\
                            'Rrs_sat_885.00':curr_sat_rrs_mean[14],\
                            'Rrs_sat_1020.50':curr_sat_rrs_mean[15],\
                            'bands':curr_bands,\
                            'MU':MU_flag,\
                            'outlier':outlier_flag,\
                            'version':sat_proc_version_str,\
                            'ws':ws},ignore_index=True)
            
        # scatter plot all bands
        plt.gca().set_aspect('equal', adjustable='box')

        if options['plot_options']['to_plot'] == 'rhow':
            plt.xlabel(f'{insitu_sensor} '+r'$\rho_{W}$',fontsize=12)
            plt.ylabel(r'OLCI $\rho_{W}$',fontsize=12)    
        elif options['plot_options']['to_plot'] == 'Rrs':
            plt.xlabel(f'{insitu_sensor} '+r'$R_{rs}$',fontsize=12)
            plt.ylabel(r'OLCI $R_{rs}$',fontsize=12) 
        elif options['plot_options']['to_plot'] == 'LWN':
            plt.xlabel(f'{insitu_sensor} '+r'$L_{WN}$',fontsize=12)
            plt.ylabel(r'OLCI $L_{WN}$',fontsize=12)               


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
            brdf_str = 'BRDF'
        else:
            brdf_str = ''
        ofname = ofname +'_'+ sat_proc_version_str.replace(' ','_').replace('.','p')+f'_{brdf_str}'+ '.pdf'
        ofname = os.path.join(output_directory,ofname)
        if 'T' in options['plot_options']['save_plot']:
            plt.savefig(ofname,dpi=300)
        print(ofname)
        print(f'N match-ups: {N_MUs}')
        
        columns = ['proc','sat',\
                 'platform','wl','N','rmse','MAPD','MPD','r^2']
        df = pd.DataFrame(columns=columns)
        
        # plot by band
        for sat_band_index in range(len(curr_bands)):
            x = matchups_ins_array[:,sat_band_index]
            y = matchups_sat_array[:,sat_band_index]
            ins_band_index = np.argmin(np.abs(curr_bands[sat_band_index]-insitu_bands))
            ins_band = insitu_bands[ins_band_index]
            sat_band = curr_bands[sat_band_index]
            df = plot_scatter(x,y,sat_band,ins_band,output_directory,'ba','OLCI','Venise_PAN',sat_proc_version_str,\
                         satellite_sensor,platform,res,insitu_sensor,brdf_str,options,df)

        print(df.to_latex(index=False,float_format="{:0.2f}".format))
        ofname_csv = os.path.join(output_directory,f'metrics_{satellite_sensor}{platform}_{res}_{sat_proc_version_str[-5:].replace(".","p")}.csv')
        print(ofname_csv)
        df.to_csv(ofname_csv,index=False)
        
        # outliers 
        plt.figure()
        plt.scatter(outliers_ins_753p75_values,outliers_sat_753p75_values)
        for i, txt in enumerate(outliers_ins_753p75):
            plt.annotate(txt[:31], (outliers_ins_753p75_values[i],outliers_sat_753p75_values[i]),rotation=40)
        plt.xlim(plot_lims_LWN['753.75'])
        plt.ylim(plot_lims_LWN['753.75'])
        plt.gca().set_aspect('equal', adjustable='box')
        
        print(df_outliers.to_csv(index=False))
        
        self.df_matchups = df_matchups
# # #%%                
# def main():
#     """business logic for when running this module as the primary one!"""
#     print('Main Code!')

path_main = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/HYPERNETS/HYPERNETS_D7p2/MDB_py/'

# Read config. file
# config_file = file_config_parse.config_file
sat_proc_version = '7.00' # '6.13' 0r '7.00'
if sat_proc_version == '6.13':
    config_file = os.path.join(path_main,'MDB_reader','config_file_OLCI_PANTHYR.ini')
elif sat_proc_version == '7.00':
    config_file = os.path.join(path_main,'MDB_reader','config_file_OLCI_HYPERNETS.ini')
elif sat_proc_version == '7.01':
    config_file = os.path.join(path_main,'MDB_reader','config_file_OLCI_PANTHYR_07.01.ini')
    
if os.path.isfile(config_file) == True:
    options = config_reader(config_file)
else:
    print(config_file + ' does not exist. Please provide a valid config file path')
    sys.exit()

if options['insitu_options']['sensor'] == 'HYPERNETS':
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




#%% Plot scatter
bands = np.array(olci_band_list[:-1],dtype=float)
df_MU = dataTOplot.df_matchups

# ins_data = df_MU['LWN_ins_753.75']
# sat_data = df_MU['LWN_sat_753.75']
# plt.figure()
# plt.scatter(ins_data,sat_data)
# plt.gca().set_aspect('equal', adjustable='box')

#% Plot spectra
outlier_flag = False

df_outliers = df_MU.loc[(df_MU['outlier'] == outlier_flag)]
LWN_outliers_ins = df_outliers.loc[:,'LWN_ins_400.00':'LWN_ins_885.00'].to_numpy()
LWN_outliers_sat = df_outliers.loc[:,'LWN_sat_400.00':'LWN_sat_885.00'].to_numpy()
if outlier_flag:
    path_out = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/OCTAC/EUMETSAT/Figures/OUTLIERS'
else:
    path_out = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/OCTAC/EUMETSAT/Figures/'
    
# sat_bands = ['400.00', '412.50', '442.50', '490.00', '510.00', '560.00', '620.00', '665.00',\
#                     '673.75', '681.25', '708.75', '753.75', '778.75', '865.00', '885.00']
sat_bands = ['442.50','490.00','560.00','665.00', '753.75', '778.75']
    
for idx in range(len(LWN_outliers_ins)):
    plt.figure(figsize = (20,10))
    plt.subplot(2,2,4)
    plt.plot(bands,LWN_outliers_sat[idx],'k.-')
    plt.plot(bands,LWN_outliers_ins[idx],'b.--')
    plt.ylim([-0.5,7])
    plt.xlim([400,885])
    PDU = df_outliers.iloc[idx]['PDU']
    ws = df_outliers.iloc[idx]['ws']
    plt.suptitle(PDU+'\n'+f'IPF-OL-2 version: {sat_proc_version}; wind speed = {ws:.1f}')
    plt.legend([PDU[:3],'PANTHYR'])
    plt.ylabel(r'$L_{WN}$',fontsize=14)
    plt.xlabel('Wavelength (nm)',fontsize=14)
    
    for idx2 in range(len(sat_bands)):
        ins_data = df_MU[f'LWN_ins_{sat_bands[idx2]}']
        sat_data = df_MU[f'LWN_sat_{sat_bands[idx2]}']
        plt.subplot(2,4,idx2+1)
        plt.scatter(ins_data,sat_data)
        plt.scatter(df_outliers.iloc[idx][f'LWN_ins_{sat_bands[idx2]}'],df_outliers.iloc[idx][f'LWN_sat_{sat_bands[idx2]}'])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlim(plot_lims_LWN[f'{sat_bands[idx2]}'])
        plt.ylim(plot_lims_LWN[f'{sat_bands[idx2]}']) 
        plt.title(f'{sat_bands[idx2]} nm', x=0.5, y=0.9)
        plt.xlabel(r'PANTHYR $L_{WN}$',fontsize=10)
        plt.ylabel(r'OLCI $L_{WN}$',fontsize=10)  

        # plot 1:1 line
        xmin, xmax = plt.gca().get_xlim()
        ymin, ymax = plt.gca().get_ylim()
        plt.plot([xmin,xmax],[ymin, ymax],'--k')
        
        ofname = os.path.join(path_out,f'{PDU[:31]}_{sat_proc_version}_spectra.pdf')
        plt.savefig(ofname)
#%% Plot spectra
path_insitu = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/HYPERNETS/HYPERNETS_D7p2/MDB_py/IDATA/INSITU/HYPERNETS/BEFR/'
bands = np.array(olci_band_list,dtype=float)
df_MU = dataTOplot.df_matchups

# ins_data = df_MU['LWN_ins_753.75']
# sat_data = df_MU['LWN_sat_753.75']
# plt.figure()
# plt.scatter(ins_data,sat_data)
# plt.gca().set_aspect('equal', adjustable='box')

#% Plot spectra
Rrs_ins = df_MU.loc[:,'Rrs_ins_400.00':'Rrs_ins_1020.50'].to_numpy()
Rrs_sat = df_MU.loc[:,'Rrs_sat_400.00':'Rrs_sat_1020.50'].to_numpy()
path_out = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/HYPERNETS/HYPERNETS_D7p2/MDB_py/ODATA/Figures/HYPERNETS'
    
# sat_bands = ['400.00', '412.50', '442.50', '490.00', '510.00', '560.00', '620.00', '665.00',\
#                     '673.75', '681.25', '708.75', '753.75', '778.75', '865.00', '885.00']
# sat_bands = ['442.50','490.00','560.00','665.00', '753.75', '778.75']
    
for idx in range(len(Rrs_ins)):
    plt.figure(figsize = (20,10))
    plt.plot(bands,Rrs_sat[idx],'k.-')
    plt.plot(bands,Rrs_ins[idx],'b.--')
    insitu_filename = df_MU['insitu_filename'][idx]
    nc_fi = Dataset(os.path.join(path_insitu,insitu_filename),'r')
    ins_wavelength = nc_fi.variables['wavelength'][:]
    ins_reflectance = nc_fi.variables['reflectance'][:]
    plt.plot(ins_wavelength,ins_reflectance/np.pi)
    # plt.ylim([-0.5,7])
    # plt.xlim([400,885])
    PDU = df_outliers.iloc[idx]['PDU']
    ws = df_outliers.iloc[idx]['ws']
    plt.suptitle(PDU+'\n'+f'IPF-OL-2 version: {sat_proc_version}; wind speed = {ws:.1f}'+f'\n{insitu_filename}')
    plt.legend([PDU[:3],'HYPERNETS spectrally sampled','HYPERNETS'],fontsize=16)
    plt.ylabel(r'$R_{rs}$',fontsize=18)
    plt.xlabel('Wavelength (nm)',fontsize=18)
         
    ofname = os.path.join(path_out,f'{PDU[:31]}_{sat_proc_version}_spectra.pdf')
    plt.savefig(ofname)
        
#%% all
bands = np.array(olci_band_list,dtype=float)
plt.figure(figsize = (20,10))
# sat
plt.subplot(2,1,1)
df_outliers = df_MU.loc[(df_MU['outlier'] == False)]
Rrs_outliers_sat = df_outliers.loc[:,'LWN_sat_400.00':'LWN_sat_885.00'].to_numpy()
for idx in range(len(LWN_outliers_sat)):
    plt.plot(bands,LWN_outliers_sat[idx],'k.-')
df_outliers = df_MU.loc[(df_MU['outlier'] == True)]
LWN_outliers_sat = df_outliers.loc[:,'LWN_sat_400.00':'LWN_sat_885.00'].to_numpy()
for idx in range(len(LWN_outliers_sat)):
    plt.plot(bands,LWN_outliers_sat[idx],'r.-')
plt.ylim([-0.5,7])
plt.xlim([400,885])
plt.title(PDU[:3]+f' IPF-OL-2 version: {sat_proc_version}', x=0.5, y=0.9)
plt.ylabel(r'$L_{WN}$',fontsize=14)
# plt.xlabel('Wavelength (nm)',fontsize=14)

# ins
plt.subplot(2,1,2)
df_outliers = df_MU.loc[(df_MU['outlier'] == False)]
LWN_outliers_ins = df_outliers.loc[:,'LWN_ins_400.00':'LWN_ins_885.00'].to_numpy()
for idx in range(len(LWN_outliers_ins)):
    plt.plot(bands,LWN_outliers_ins[idx],'k.-')
df_outliers = df_MU.loc[(df_MU['outlier'] == True)]
LWN_outliers_ins = df_outliers.loc[:,'LWN_ins_400.00':'LWN_ins_885.00'].to_numpy()
for idx in range(len(LWN_outliers_ins)):
    plt.plot(bands,LWN_outliers_ins[idx],'r.-')    
plt.ylim([-0.5,7])
plt.xlim([400,885])
plt.title('PANTHYR', x=0.5, y=0.9)
plt.ylabel(r'$L_{WN}$',fontsize=14)
plt.xlabel('Wavelength (nm)',fontsize=14)

plt.suptitle('(*) in red: Outliers identified from the 778.5 nm band for IPF-OL-2 version 7.00.')

path_out = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/OCTAC/EUMETSAT/Figures/OUTLIERS'
ofname = os.path.join(path_out,f'{PDU[:3]}_{sat_proc_version}_spectra.pdf')
plt.savefig(ofname)

#%% plot 865/665
plt.figure(figsize = (20,10))
plt.subplot(1,2,1)
rhow_778_75 = df_MU['LWN_sat_778.75']*np.pi/cfs.get_F0(778.75)
rhow_665 = df_MU['LWN_sat_665.00']*np.pi/cfs.get_F0(665.00)
plt.scatter(rhow_778_75,rhow_665, facecolors='none', edgecolors='b')
# plt.gca().set_aspect('equal', adjustable='box')
plt.xlim([0,0.01])
plt.ylim([0,0.06])
plt.xlabel(f'{PDU[:3]}'+r' $\rho_{W}(778.75)$',fontsize=14)
plt.ylabel(f'{PDU[:3]}'+r' $\rho_{W}(665)$',fontsize=14)
x = np.array([0,0.01])
y = 4.228*x # from Kevin
plt.plot(x,y,'k--')


plt.subplot(1,2,2)
rhow_778_75 = df_MU['LWN_ins_778.75']*np.pi/cfs.get_F0(778.75)
rhow_665 = df_MU['LWN_ins_665.00']*np.pi/cfs.get_F0(665.00)
plt.scatter(rhow_778_75,rhow_665, facecolors='none', edgecolors='b')
# plt.gca().set_aspect('equal', adjustable='box')
plt.xlim([0,0.01])
plt.ylim([0,0.06])
plt.xlabel(r'PANTHYR $\rho_{W}(778.75)$',fontsize=14)
plt.ylabel(r'PANTHYR $\rho_{W}(665)$',fontsize=14)
plt.plot(x,y,'k--')
plt.suptitle(f'{PDU[:3]}; IPF-OL-2 version: {sat_proc_version}',fontsize=14)
ofname = os.path.join(path_out,f'{PDU[:3]}_{sat_proc_version}_779-665.pdf')
plt.savefig(ofname)
#%%
# if __name__ == '__main__':
#     main()