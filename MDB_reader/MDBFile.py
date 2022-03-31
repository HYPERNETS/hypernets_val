import os
import sys
from datetime import datetime
from datetime import timedelta

from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pandas as pd

code_home = os.path.abspath('../')
sys.path.append(code_home)
import COMMON.Class_Flags_OLCI as flag

DIMENSION_NAMES = ['satellite_id', 'rows', 'columns', 'satellite_bands', 'insitu_id', 'insitu_original_bands']
VAR_NAMES = ['satellite_time', 'satellite_latitude', 'satellite_longitude', 'satellite_bands',
             'satellite_Rrs', 'insitu_time', 'insitu_original_bands', 'insitu_Rrs', 'time_difference']

ATRIB_NAMES = ['creation_time', 'satellite', 'platform', 'sensor', 'description', 'satellite_aco_processor',
               'satellite_proc_version', 'insitu_site_name', 'insitu_lat', 'insitu_lon', 'time_diff']


class MDBFile:

    def __init__(self, file_path):
        global DIMENSION_NAMES
        global VAR_NAMES
        global ATRIB_NAMES
        self.file_path = file_path
        self.VALID = True
        file_name = file_path.split('/')[-1]
        print(f'[INFO]Starting MDBFile: {file_name}')
        try:
            self.nc = Dataset(file_path)
            self.variables = self.nc.variables
            self.dimensions = self.nc.dimensions
            self.flag_band_name = 'satellite_WQSF'
            self.VALID = self.check_structure()
        # except:
        except Exception as e:
            self.VALID = False
            print(f'Exception starting MDB File: {e}')

        if not self.VALID:
            print(f'MDB File: {file_path} is not valid')

        if self.VALID:
            self.n_mu_total = len(self.dimensions['satellite_id'])
            print('Total mu: ', self.n_mu_total)
            self.sat_times = []
            for st in self.variables['satellite_time']:
                self.sat_times.append(datetime(1970, 1, 1) + timedelta(seconds=int(st)))
            self.start_date = self.sat_times[0]
            self.end_date = self.sat_times[-1]
            self.insitu_bands = self.nc.variables['insitu_original_bands'][:]  # insitu_bands(insitu_bands)
            self.satellite_bands = self.nc.variables['satellite_bands'][:]
            self.n_insitu_bands = len(self.insitu_bands)
            self.n_satellite_bands = len(self.satellite_bands)
            self.info = {}
            atribs = self.nc.ncattrs()
            for at in atribs:
                self.info[at] = self.nc.getncattr(at)
            self.wlref = self.satellite_bands
            self.wlref_sat_indices = list(range(len(self.satellite_bands)))
            # print(self.info)
            # print(self.flag_band_name)

        # Sat filtering options
        self.flag_list = ''
        self.window_size = 3
        self.valid_min_pixels = 1
        self.delta_t = 7200
        self.set_default_sat_filtering_options()

        # Variables defining a specific MU. To load a MU, uses load_mu_data
        self.index_mu = 0
        self.mu_sat_time = []
        self.mu_insitu_time = []
        self.ins_time_index = -1
        self.satellite_rrs = []
        self.insitu_rrs = []
        self.mu_curr_sat_rrs_mean = []
        self.mu_curr_ins_rrs = []
        self.mu_valid_bands = []

        # Variables to make validation...
        col_names = ['Index', 'Index_MU', 'Index_Band', 'Sat_Time', 'Ins_Time', 'Time_Diff', 'Wavelenght', 'Ins_Rrs',
                     'Sat_Rrs', 'Valid']
        self.df_validation = pd.DataFrame(columns=col_names)
        self.df_validation_valid = pd.DataFrame(columns=col_names)
        self.mu_dates = {}

    # Checking atrib, var and dimensions names
    def check_structure(self):
        check_var = True
        check_dim = True
        check_atrib = True
        for var in VAR_NAMES:
            if var not in self.variables:
                print(var)
                check_var = False
        for dim in DIMENSION_NAMES:
            if dim not in self.dimensions:
                check_dim = False
        for atrib in ATRIB_NAMES:
            if atrib not in self.nc.ncattrs():
                check_atrib = False

        if check_var == False or check_dim == False or check_atrib == False:
            return False

        if not self.flag_band_name in self.variables:
            if self.nc.satellite_aco_processor == 'Polymer' and 'satellite_bitmask' in self.variables:
                self.flag_band_name = 'satellite_bitmask'

            if self.nc.satellite_aco_processor == 'C2RCC' and 'satellite_c2rcc_flags' in self.variables:
                self.flag_band_name = 'satellite_c2rcc_flags'

        return True

    # Set default sat filtering options
    def set_default_sat_filtering_options(self):
        flag_list = 'CLOUD,CLOUD_AMBIGUOUS,CLOUD_MARGIN,INVALID,COSMETIC,SATURATED,SUSPECT,HISOLZEN,HIGHGLINT,SNOW_ICE,AC_FAIL,WHITECAPS,RWNEG_O2,RWNEG_O3,RWNEG_O4,RWNEG_O5,RWNEG_O6,RWNEG_O7,RWNEG_O8'
        flag_list = flag_list.replace(" ", "")
        self.flag_list = str.split(flag_list, ',')
        self.window_size = 3
        self.valid_min_pixels = 1
        self.delta_t = 7200

    def get_dimensions(self):
        # Dimensions
        central_r = int(np.floor(len(self.dimensions['rows']) / 2))
        central_c = int(np.floor(len(self.dimensions['columns']) / 2))
        r_s = central_r - int(np.floor(self.window_size / 2))  # starting row
        r_e = central_r + int(np.floor(self.window_size / 2)) + 1  # ending row
        c_s = central_c - int(np.floor(self.window_size / 2))  # starting col
        c_e = central_c + int(np.floor(self.window_size / 2)) + 1  # ending col
        return central_r, central_c, r_s, r_e, c_s, c_e

    def retrieve_ins_info_mu(self, index_mu):
        time_difference = self.variables['time_difference'][index_mu]
        ins_time_index = np.argmin(np.abs(time_difference))
        ins_time = self.variables['insitu_time'][index_mu][ins_time_index]
        mu_insitu_time = datetime.fromtimestamp(
            int(ins_time))  # datetime(1970, 1, 1) + timedelta(seconds=int(ins_time))
        time_condition = False
        if np.abs(time_difference[ins_time_index]) < self.delta_t:
            time_condition = True
        return ins_time_index, mu_insitu_time, time_condition

    def get_flag_mask(self, index_mu):
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()

        if self.flag_band_name == 'satellite_bitmask':  # POLYMER BIT MASK
            flagging = flag.Class_Flags_Polymer(self.variables[self.flag_band_name].flag_masks,
                                                self.variables[self.flag_band_name].flag_meanings)
            satellite_flag_band = self.variables[self.flag_band_name][index_mu, r_s:r_e, c_s:c_e]
            mask = flagging.MaskGeneral(satellite_flag_band)
            mask[np.where(mask != 0)] = 1
            NGP = np.power(self.window_size, 2) - np.sum(mask)  # Number Good Pixels excluding flagged pixels
            land = flagging.Mask(satellite_flag_band, (['LAND']))
            inland_w = flagging.Mask(satellite_flag_band, (['CASE2']))
            land[np.where(inland_w > 0)] = 0
            NTP = np.power(self.window_size, 2) - np.sum(land, axis=(0, 1))
        elif self.flag_band_name == 'satellite_c2rcc_flags':  # C2RCC FLAGS
            # flagging = flag.Class_Flags_OLCI(self.variables[self.flag_band_name].flag_masks,
            #                                  self.variables[self.flag_band_name].flag_meanings)
            satellite_flag_band = self.variables[self.flag_band_name][index_mu, r_s:r_e, c_s:c_e]
            valuePE = np.uint64(2147483648)
            mask = np.ones(satellite_flag_band.shape, dtype=np.uint64)
            mask[satellite_flag_band == valuePE] = 0
            NGP = np.power(self.window_size, 2) - np.sum(mask)  # Number Good Pixels excluding flagged pixels
            NTP = np.power(self.window_size, 2)
        else:
            flagging = flag.Class_Flags_OLCI(self.variables[self.flag_band_name].flag_masks,
                                             self.variables[self.flag_band_name].flag_meanings)
            satellite_flag_band = self.variables[self.flag_band_name][index_mu, r_s:r_e, c_s:c_e]
            if (str(self.flag_list[0])) != 'None':
                mask = flagging.Mask(satellite_flag_band, self.flag_list)
                mask[np.where(mask != 0)] = 1
            else:
                mask = np.full(satellite_flag_band.shape, 0, dtype=int)
            NGP = np.power(self.window_size, 2) - np.sum(mask)  # Number Good Pixels excluding flagged pixels
            land = flagging.Mask(satellite_flag_band, (['LAND']))
            inland_w = flagging.Mask(satellite_flag_band, (['INLAND_WATER']))
            land[np.where(inland_w > 0)] = 0
            NTP = np.power(self.window_size, 2) - np.sum(land, axis=(0, 1))

        # conditions minimum valid pixels
        if float(self.valid_min_pixels) == 1:  # 100% like Zibordi
            cond_min_valid_pxs = NGP == np.power(self.window_size, 2)
        else:
            cond_min_valid_pxs = NGP > (float(self.valid_min_pixels) * NTP + 1)

        return mask, cond_min_valid_pxs, NGP, NTP

    # Funcion to load data from a specific MU
    def load_mu_data(self, index_mu):
        if not self.VALID:
            return False

        if index_mu < 0 or index_mu >= self.n_mu_total:
            print('Not valid index_mu')
            return False

        # Index match-up
        self.index_mu = index_mu

        # Sat and instrument rrs
        self.insitu_rrs = self.variables['insitu_Rrs'][index_mu]
        # self.insitu_rrs = self.insitu_rrs / np.pi  # transform from rhow to Rr
        self.satellite_rrs = self.variables['satellite_Rrs'][index_mu]

        # Sat and instrument time
        self.mu_sat_time = self.sat_times[index_mu]
        self.ins_time_index, self.mu_insitu_time, time_condition = self.retrieve_ins_info_mu(index_mu)

        # Dimensions
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()

        # Flag mask
        mask, cond_min_valid_pxs, NGP, NTP = self.get_flag_mask(index_mu)

        # TEMPORAL
        filtro_temporal = True
        if time_condition and cond_min_valid_pxs:
            for iref in range(len(self.wlref_sat_indices)):  # range(0, self.n_satellite_bands):
                sat_band_index = self.wlref_sat_indices[iref]

                wl = self.satellite_bands[sat_band_index]
                # print('------------------------------', sat_band_index, wl)
                # print(sat_band_index, wl)
                ins_band_index = np.argmin(np.abs(wl - self.insitu_bands))
                if self.insitu_rrs.mask[ins_band_index, self.ins_time_index]:
                    filtro_temporal = False
                if np.isnan(self.insitu_rrs[ins_band_index, self.ins_time_index]):
                    filtro_temporal = False
                    continue
                val = self.insitu_rrs[ins_band_index, self.ins_time_index]
                if val > 0.008 and sat_band_index == 6:
                    filtro_temporal = False
                if val > 0.004 and sat_band_index == 7:
                    filtro_temporal = False
                if val < 0:
                    filtro_temporal = False

        # Getting spectra for comparison
        self.mu_valid_bands = [False] * len(self.wlref_sat_indices)

        if time_condition and cond_min_valid_pxs and filtro_temporal:
            self.mu_curr_ins_rrs = []
            self.mu_curr_sat_rrs_mean = []
            for iref in range(len(self.wlref_sat_indices)):
                sat_band_index = self.wlref_sat_indices[iref]
                wl = self.satellite_bands[sat_band_index]
                curr_sat_box = self.satellite_rrs[sat_band_index, r_s:r_e, c_s:c_e]
                ins_band_index = np.argmin(np.abs(wl - self.insitu_bands))
                difwl = abs(wl - self.insitu_bands[ins_band_index])
                check_dif_wl = True
                check_ins_value = True
                if difwl > 10:
                    check_dif_wl = False

                if self.insitu_rrs.mask[ins_band_index, self.ins_time_index]:
                    check_ins_value = False
                if np.isnan(self.insitu_rrs[ins_band_index, self.ins_time_index]):
                    check_ins_value = False
                if not np.isnan(self.insitu_rrs[ins_band_index, self.ins_time_index]) and self.insitu_rrs[
                    ins_band_index, self.ins_time_index] < 0:
                    check_ins_value = False

                curr_sat_box = np.ma.masked_where(curr_sat_box == -999, curr_sat_box)
                curr_sat_box = np.ma.masked_invalid(curr_sat_box)
                curr_sat_box = np.ma.masked_array(curr_sat_box, mask)
                numberValid = np.ma.count(curr_sat_box)
                if float(self.valid_min_pixels) == 1:  # 100% like Zibordi
                    cond_min_valid_pxs = numberValid == np.power(self.window_size, 2)
                else:
                    cond_min_valid_pxs = numberValid > (float(self.valid_min_pixels) * NTP + 1)

                if cond_min_valid_pxs and check_dif_wl and check_ins_value:
                    curr_sat_box_mean = curr_sat_box.mean()
                    self.mu_curr_sat_rrs_mean.append(curr_sat_box_mean)
                    ins_value = self.insitu_rrs[ins_band_index, self.ins_time_index]
                    self.mu_curr_ins_rrs.append(ins_value)
                    self.mu_valid_bands[iref] = True

            # print(self.mu_curr_sat_rrs_mean)
            # print(self.mu_curr_ins_rrs)
        if sum(self.mu_valid_bands) == len(self.wlref_sat_indices):
            return True
        else:
            return False

    def prepare_df_validation(self):
        # ['Index','Index_MU','Index_Band','Sat_Time','Ins_Time','Time_Diff','Wavelenght','Ins_Rrs','Sat_Rrs','Valid']
        index_tot = 0
        index_valid_tot = 0
        nmu_valid = 0
        for index_mu in range(self.n_mu_total):
            mu_valid = self.load_mu_data(index_mu)
            if mu_valid:
                nmu_valid = nmu_valid + 1
            sdate = self.mu_sat_time.strftime('%Y-%m-%d')
            if not sdate in self.mu_dates.keys():
                self.mu_dates[sdate] = {
                    'satellite': self.info['satellite'].upper(),
                    'platform': self.info['platform'].upper(),
                    'sensor': self.info['sensor'].upper(),
                    'site': self.info['insitu_site_name'].upper(),
                    'ac': self.info['satellite_aco_processor'].upper(),
                    'mu_valid': mu_valid,
                }
                if self.mu_dates[sdate]['ac'] == 'ATMOSPHERIC CORRECTION PROCESSOR: XXX':
                    self.mu_dates[sdate]['ac'] = 'STANDARD'
            else:
                print('[WARNING] A single MDB file should not contain more than one match-up in a specific date')

            # print(f'Match-up # {index_mu} Valid: {mu_valid}')
            index_valid = 0
            # for sat_band_index in range(0, self.n_satellite_bands):
            for iref in range(len(self.wlref_sat_indices)):  # range(0, self.n_satellite_bands):
                sat_band_index = self.wlref_sat_indices[iref]

                time_diff = round(abs((self.mu_sat_time - self.mu_insitu_time).total_seconds() / 3600), 2)
                row = {
                    'Index': [index_tot],
                    'Index_MU': [index_mu],
                    'Index_Band': [sat_band_index],
                    'Sat_Time': [self.mu_sat_time.strftime('%Y-%m-%d %H:%M')],
                    'Ins_Time': [self.mu_insitu_time.strftime('%Y-%m-%d %H:%M')],
                    'Time_Diff': [time_diff],
                    'Wavelenght': [self.wlref[iref]],
                    'Ins_Rrs': [-999],
                    'Sat_Rrs': [-999],
                    'Valid': [mu_valid]  # [self.mu_valid_bands[sat_band_index]]
                }
                if mu_valid:  # self.mu_valid_bands[sat_band_index]:
                    row['Ins_Rrs'] = [self.mu_curr_ins_rrs[index_valid]]
                    row['Sat_Rrs'] = [self.mu_curr_sat_rrs_mean[index_valid]]
                    index_valid = index_valid + 1

                self.df_validation = pd.concat([self.df_validation, pd.DataFrame.from_dict(row)], ignore_index=True)
                index_tot = index_tot + 1

                if mu_valid:
                    row['Index'] = index_valid_tot
                    self.df_validation_valid = pd.concat([self.df_validation_valid, pd.DataFrame.from_dict(row)],
                                                         ignore_index=True)
                    index_valid_tot = index_valid_tot + 1

        print(f'# total match-ups: {self.n_mu_total} Valid: {nmu_valid}')

    ##PLOT FUNCTIONS
    def plot_spectra(self, hfig):
        print('Plot spectra...')
        legend = [f"S3{self.nc.platform}-{self.mu_sat_time.strftime('%d-%m-%Y %H:%M')}",
                  f"Hypstar-{self.mu_insitu_time.strftime('%d-%m-%Y %H:%M')}"]
        row_names = []
        for w in self.satellite_bands:
            wn = f'{w}'
            row_names.append(wn)
        ydata = np.array([self.mu_curr_sat_rrs_mean, self.mu_curr_ins_rrs])
        df = pd.DataFrame(np.transpose(ydata), columns=legend, index=row_names)
        title = f"Match-up #{self.index_mu} {self.mu_sat_time.strftime('%d-%m-%Y')}"

        if hfig is None:
            hfig = plt.figure()
        df.plot(lw=2, marker='.', markersize=10)
        rindex = list(range(16))
        plt.xticks(rindex, row_names, rotation=45, fontsize=12)
        plt.xlabel('Wavelength(nm)', fontsize=14)
        plt.ylabel('R$_{rs}$[1/sr]', fontsize=14)
        plt.title(title, fontsize=16)
        plt.gcf().subplots_adjust(bottom=0.18)
        plt.gcf().subplots_adjust(left=0.20)
        plt.gcf().canvas.draw()

        return hfig

        #
        # ofname = f"Spectra_S3{platform}_{sat_time.strftime('%Y%m%d')}.jpg"
        # ofname = os.path.join(path_out, ofname)
        # print('Name out:', ofname)
        #
        # plt.savefig(ofname, dpi=300)
        # plt.close()
        #
        # return ofname

    def get_file_name_base(self):
        sat = self.info['satellite']
        platform = self.info['platform']
        res = 'WFR'
        insitu_sensor = 'AERONET'

        insitu_site = self.info['insitu_site_name']
        ofname = f'{sat}{platform}_{res}_{insitu_sensor}_{insitu_site}'
        return ofname

    def get_file_name(self, wl):
        dates = self.start_date.strftime('%Y%m%d') + '_' + self.end_date.strftime('%Y%m%d')
        if wl is None:
            file_name = self.get_file_name_base() + f'_{dates}'
        else:
            file_name = self.get_file_name_base() + f'_{wl}_{dates}'
        return file_name

    def get_file_name_param(self, param):
        dates = self.start_date.strftime('%Y%m%d') + '_' + self.end_date.strftime('%Y%m%d')
        file_name = self.get_file_name_base() + f'_{dates}_{param}'
        return file_name

    def get_title(self):
        dates = self.start_date.strftime('%Y-%m-%d') + ' to ' + self.end_date.strftime('%Y-%m-%d')
        title = self.get_file_name_base() + f' {dates}'
        return title

    def set_wl_ref(self, wllist):
        self.wlref = wllist
        self.wlref_sat_indices = []
        for wl in wllist:
            index = np.argmin(np.abs(wl - self.satellite_bands))
            self.wlref_sat_indices.append(index)
