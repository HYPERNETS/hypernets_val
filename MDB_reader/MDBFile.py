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

from QC_INSITU import QC_INSITU
from QC_SAT import QC_SAT

DIMENSION_NAMES = ['satellite_id', 'rows', 'columns', 'satellite_bands', 'insitu_id', 'insitu_original_bands']
VAR_NAMES = ['satellite_time', 'satellite_latitude', 'satellite_longitude', 'satellite_bands',
             'satellite_Rrs', 'insitu_time', 'insitu_original_bands', 'insitu_Rrs', 'time_difference']

ATRIB_NAMES = ['creation_time', 'satellite', 'platform', 'sensor', 'description', 'satellite_aco_processor',
               'satellite_proc_version', 'insitu_site_name', 'insitu_lat', 'insitu_lon']


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
            print('[INFO ]Total mu: ', self.n_mu_total)
            self.sat_times = []

            for st in self.variables['satellite_time']:
                self.sat_times.append(datetime.fromtimestamp(float(st)))
            # for idx in range(len(self.nc.variables['satellite_PDU'])):
            #     pdu = self.nc.variables['satellite_PDU'][idx]
            #     time_pdu = self.get_sat_time_from_fname(pdu)
            #     # self.sat_times.append(self.get_sat_time_from_fname(pdu))
            #     st = self.variables['satellite_time'][idx]
            #     # dt_1 = datetime(1970, 1, 1) + timedelta(seconds=int(st))
            #     # dt_2 = datetime(1970, 1, 1) + timedelta(seconds=float(st))
            #     dt_3 = datetime.fromtimestamp(float(st))
            #     # print(time_pdu,dt_1,dt_2,dt_3)
            #     self.sat_times.append(dt_3)

            # for st in self.variables['satellite_time']:
            #     #self.sat_times.append(datetime(1970, 1, 1) + timedelta(seconds=int(st)))
            #     self.sat_times.append(datetime.fromtimestamp(float(st)))
            self.start_date = self.sat_times[0]
            self.end_date = self.sat_times[-1]
            self.insitu_bands = self.nc.variables['insitu_original_bands'][:]  # insitu_bands(insitu_bands)
            self.satellite_bands = self.nc.variables['satellite_bands'][:]
            self.n_insitu_bands = len(self.insitu_bands)
            self.n_satellite_bands = len(self.satellite_bands)
            self.info = {}
            atribs = self.nc.ncattrs()
            for at in atribs:
                if at == 'satellite_aco_processor':
                    self.info[at] = self.nc.getncattr(at).upper()
                else:
                    self.info[at] = self.nc.getncattr(at)

            if self.info['satellite_aco_processor'] == 'ATMOSPHERIC CORRECTION PROCESSOR: XXX':
                self.info['satellite_aco_processor'] = 'STANDARD'

            if self.info['satellite_aco_processor'] == 'CLIMATE CHANGE INITIATIVE - EUROPEAN SPACE AGENCY':
                self.info['satellite_aco_processor'] = 'CCI'

            self.wlref = self.satellite_bands
            self.wlref_sat_indices = list(range(len(self.satellite_bands)))
            # self.set_wl_ref_insitu_indices()

        self.delta_t = 7200

        # Variables defining a specific MU. To load a MU, uses load_mu_data
        self.index_mu = 0
        self.mu_sat_time = []
        self.mu_insitu_time = []
        self.ins_time_index = -1
        self.satellite_rrs = []
        self.insitu_rrs = []
        self.mu_curr_sat_rrs_mean = []
        self.mu_curr_ins_rrs = []

        # QUALITY CONTROL

        self.only_complete_spectra_valid = False
        self.qc_insitu = QC_INSITU(self.variables['insitu_Rrs'], self.variables['insitu_original_bands'])
        self.qc_insitu.time_max = self.delta_t
        self.qc_insitu.set_wllist_using_wlref(self.wlref)

        if self.nc.satellite_aco_processor == 'ACOLITE' or self.nc.satellite_aco_processor == 'Climate Change Initiative - European Space Agency':
            self.qc_sat = QC_SAT(self.variables['satellite_Rrs'], self.satellite_bands, None,
                                 self.info['satellite_aco_processor'])
        elif len(self.nc.satellite_aco_processor) == 0:
            self.qc_sat = QC_SAT(self.variables['satellite_Rrs'], self.satellite_bands, None, 'Climate Change Initiative - European Space Agency')
        else:
            self.qc_sat = QC_SAT(self.variables['satellite_Rrs'], self.satellite_bands,
                                 self.variables[self.flag_band_name], self.info['satellite_aco_processor'])

        # Variables to make validation...
        self.col_names = ['Index', 'Index_MU', 'Index_Band', 'Sat_Time', 'Ins_Time', 'Time_Diff', 'Wavelenght',
                          'Ins_Rrs',
                          'Sat_Rrs', 'Valid']
        self.df_validation = None
        self.df_validation_valid = None
        self.mu_dates = {}

    # Checking atrib, var and dimensions names
    def check_structure(self):
        check_var = True
        check_dim = True
        check_atrib = True
        for var in VAR_NAMES:
            if var not in self.variables:
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
            if self.nc.satellite_aco_processor.upper() == 'POLYMER' and 'satellite_bitmask' in self.variables:
                self.flag_band_name = 'satellite_bitmask'

            if self.nc.satellite_aco_processor.upper() == 'C2RCC' and 'satellite_c2rcc_flags' in self.variables:
                self.flag_band_name = 'satellite_c2rcc_flags'

            if self.nc.satellite_aco_processor.upper() == 'FUB' and 'satellite_quality_flags':
                self.flag_band_name = 'satellite_quality_flags'

            if self.nc.satellite_aco_processor.upper() == 'ACOLITE':
                self.flag_band_name = 'NONE'

        return True

    # Set qc sat filtering options
    def set_default_filtering_options(self):
        print(self.qc_sat.info_flag)
        # flag_list = 'CLOUD,CLOUD_AMBIGUOUS,CLOUD_MARGIN,INVALID,COSMETIC,SATURATED,SUSPECT,HISOLZEN,HIGHGLINT,SNOW_ICE,AC_FAIL,WHITECAPS,RWNEG_O2,RWNEG_O3,RWNEG_O4,RWNEG_O5,RWNEG_O6,RWNEG_O7,RWNEG_O8'
        # if self.nc.satellite_aco_processor == 'FUB':
        #     flag_list = 'land'
        #     # flag_list = 'land,coastline,fresh_inland_water,bright,straylight_risk,invalid,cosmetic,duplicated,sun_glint_risk,dubious,saturated_Oa01,saturated_Oa02,saturated_Oa03,saturated_Oa04,saturated_Oa05,saturated_Oa06,saturated_Oa07,saturated_Oa08,saturated_Oa09,saturated_Oa10,saturated_Oa11,saturated_Oa12,saturated_Oa13,saturated_Oa14,saturated_Oa15,saturated_Oa16,saturated_Oa17,saturated_Oa18,saturated_Oa19,saturated_Oa20,saturated_Oa21'
        #
        # flag_list = flag_list.replace(" ", "")
        # self.flag_list = str.split(flag_list, ',')
        # self.window_size = 3
        # self.valid_min_pixels = 1
        # self.delta_t = 7200

    def get_dimensions(self):
        # Dimensions
        central_r = int(np.floor(len(self.dimensions['rows']) / 2))
        central_c = int(np.floor(len(self.dimensions['columns']) / 2))
        r_s = central_r - int(np.floor(self.window_size / 2))  # starting row
        r_e = central_r + int(np.floor(self.window_size / 2)) + 1  # ending row
        c_s = central_c - int(np.floor(self.window_size / 2))  # starting col
        c_e = central_c + int(np.floor(self.window_size / 2)) + 1  # ending col
        return central_r, central_c, r_s, r_e, c_s, c_e

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

        elif self.flag_band_name == 'NONE':  # ACOLITE
            index_1020 = np.argmin(np.abs(1020 - self.satellite_bands))
            bandref = self.variables['satellite_Rrs'][index_mu, index_1020, r_s:r_e, c_s:c_e]
            mask = np.zeros(bandref.shape, dtype=np.uint64)
            mask[bandref > 0.05] = 1
            NTP = np.power(self.window_size, 2) - np.sum(mask)
            NGP = np.power(self.window_size, 2)
        else:  # STANDARD/FUB
            flagging = flag.Class_Flags_OLCI(self.variables[self.flag_band_name].flag_masks,
                                             self.variables[self.flag_band_name].flag_meanings)
            satellite_flag_band = self.variables[self.flag_band_name][index_mu, r_s:r_e, c_s:c_e]
            if (str(self.flag_list[0])) != 'None':
                mask = flagging.Mask(satellite_flag_band, self.flag_list)
                mask[np.where(mask != 0)] = 1
            else:
                mask = np.full(satellite_flag_band.shape, 0, dtype=int)

            if self.nc.satellite_aco_processor == 'ACOLITE':
                index_1020 = np.argmin(np.abs(412 - self.satellite_bands))
                bandref = self.variables['satellite_Rrs'][index_mu, index_1020, r_s:r_e, c_s:c_e]
                mask[bandref > 0.015] = 1

            NGP = np.power(self.window_size, 2) - np.sum(mask)  # Number Good Pixels excluding flagged pixels

            if self.nc.satellite_aco_processor == 'FUB':
                land = flagging.Mask(satellite_flag_band, (['land']))
                inland_w = flagging.Mask(satellite_flag_band, (['fresh_inland_water']))
            else:
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

    def retrieve_ins_info_mu_spectra(self, index_mu):
        time_difference_prev = self.variables['time_difference'][index_mu]
        time_difference = np.ma.copy(time_difference_prev)
        times_here = self.variables['insitu_time'][index_mu]
        sat_time_here = self.sat_times[index_mu]
        for idx in range(len(times_here)):
            itime = times_here[idx]
            if not np.ma.is_masked(itime):
                # insitu_time_here = datetime(1970, 1, 1) + timedelta(seconds=int(itime))
                insitu_time_here = datetime.fromtimestamp(float(itime))
                time_diff_here = abs((sat_time_here - insitu_time_here).total_seconds())
                time_difference[idx] = time_diff_here

        # time_difference = time_difference_prev

        if 'insitu_exact_wavelenghts' in self.variables:
            exact_wl = self.variables['insitu_exact_wavelenghts'][index_mu]
        else:
            exact_wl = self.variables['insitu_original_bands']

        ins_time_index, time_condition, valid_insitu, spectrum_complete, rrs_values = \
            self.qc_insitu.get_finalspectrum_mu(index_mu, time_difference, exact_wl, self.wlref)

        if time_condition and valid_insitu:
            ins_time = self.variables['insitu_time'][index_mu][ins_time_index]
            mu_insitu_time = datetime.fromtimestamp(int(ins_time))
        else:  ##aunque los datos sean invalidos (time dif>max time dif), obtenemos el mu_insitu_time como referencia
            ins_time_index = np.argmin(np.abs(time_difference))
            ins_time = self.variables['insitu_time'][index_mu][ins_time_index]
            mu_insitu_time = datetime.fromtimestamp(int(ins_time))

        return ins_time_index, mu_insitu_time, time_condition, valid_insitu, spectrum_complete, rrs_values

    # def retrieve_ins_info_mu(self, index_mu):
    #     time_difference = self.variables['time_difference'][index_mu]
    #     ins_time_index = np.argmin(np.abs(time_difference))
    #     ins_time = self.variables['insitu_time'][index_mu][ins_time_index]
    #     mu_insitu_time = datetime.fromtimestamp(
    #         int(ins_time))  # datetime(1970, 1, 1) + timedelta(seconds=int(ins_time))
    #     time_condition = False
    #     if np.abs(time_difference[ins_time_index]) < self.delta_t:
    #         time_condition = True
    #     return ins_time_index, mu_insitu_time, time_condition

    def load_mu_datav2(self, index_mu):

        is_mu_valid = False
        load_info = {
            'status': '',
            'spectrum_complete': False,
            'valid_bands': []
        }

        if not self.VALID:
            load_info['status'] = -1  # 'NO VALID MDB FILE'
            return is_mu_valid, load_info

        if index_mu < 0 or index_mu >= self.n_mu_total:
            load_info['status'] = -2  # f'NO VALID MATCH-UP INDEX:{index_mu}'
            return is_mu_valid, load_info

        # Index match-up
        self.index_mu = index_mu

        # Sat and instrument rrs
        self.insitu_rrs = self.variables['insitu_Rrs'][index_mu]
        self.satellite_rrs = self.variables['satellite_Rrs'][index_mu]

        # Sat and instrument time
        self.mu_sat_time = self.sat_times[index_mu]
        # THIS STEP IS NOW DONE BEFORE PREPARING DF FOR VALIDATION
        # if self.info['satellite_aco_processor'] == 'CCI':
        #     self.mu_sat_time = self.mu_sat_time.replace(hour=11)

        # self.ins_time_index, self.mu_insitu_time, time_condition = self.retrieve_ins_info_mu(index_mu)
        self.ins_time_index, self.mu_insitu_time, time_condition, valid_insitu, spectrum_complete, rrs_ins_values = \
            self.retrieve_ins_info_mu_spectra(index_mu)

        load_info['spectrum_complete'] = spectrum_complete

        # if not time_condition:
        #     load_info['status'] = -3  # f'IN SITU DATA OUT OF TIME WINDOW'
        #     return is_mu_valid, load_info

        if not valid_insitu:
            load_info['status'] = -4  # f'INVALID INSITU DATA'
            return is_mu_valid, load_info

        if not spectrum_complete and self.qc_insitu.only_complete_spectra:
            load_info['status'] = -5  # f'INCOMPLETE IN SITU SPECTRUM'
            return is_mu_valid, load_info

        cond_min_pixels, cond_stats, valid_mu, sat_values = self.qc_sat.get_match_up_values(index_mu)
        if not valid_mu:
            load_info['status'] = -6  # f'NO VALID SAT DATA'
            return is_mu_valid, load_info

        # Getting spectra for comparison
        mu_valid_bands = [False] * len(self.wlref_sat_indices)
        self.mu_curr_ins_rrs = []
        self.mu_curr_sat_rrs_mean = []
        for iref in range(len(self.wlref_sat_indices)):
            # sat_band_index = self.wlref_sat_indices[iref]
            check_ins_value = True
            if rrs_ins_values.mask[iref]:
                check_ins_value = False
            if check_ins_value:
                self.mu_curr_sat_rrs_mean.append(sat_values[iref])
                self.mu_curr_ins_rrs.append(rrs_ins_values[iref])
                mu_valid_bands[iref] = True

        if sum(mu_valid_bands) == 0:
            load_info['status'] = -6  # f'NO VALID IN SITU DATA'
            return is_mu_valid, load_info

        load_info['status'] = 1  # f'OK'
        is_mu_valid = True
        if sum(mu_valid_bands) == len(self.wlref_sat_indices):
            load_info['spectrum_complete'] = True
            return is_mu_valid, load_info
        else:
            load_info['spectrum_complete'] = False
            load_info['valid_bands'] = mu_valid_bands
            return is_mu_valid, load_info

    # Funcion to load data from a specific MU
    # def load_mu_data(self, index_mu):
    #     if not self.VALID:
    #         return False
    #
    #     if index_mu < 0 or index_mu >= self.n_mu_total:
    #         print('Not valid index_mu')
    #         return False
    #
    #     # Index match-up
    #     self.index_mu = index_mu
    #
    #     # Sat and instrument rrs
    #     self.insitu_rrs = self.variables['insitu_Rrs'][index_mu]
    #     # self.insitu_rrs = self.insitu_rrs * np.pi  # transform from rhow to Rr
    #
    #     self.satellite_rrs = self.variables['satellite_Rrs'][index_mu]
    #
    #     # Sat and instrument time
    #     self.mu_sat_time = self.sat_times[index_mu]
    #     self.ins_time_index, self.mu_insitu_time, time_condition = self.retrieve_ins_info_mu(index_mu)
    #
    #     # Dimensions
    #     central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
    #
    #     # Flag mask
    #     mask, cond_min_valid_pxs, NGP, NTP = self.get_flag_mask(index_mu)

    # TEMPORAL
    # filtro_temporal = True
    # if time_condition and cond_min_valid_pxs:
    #     for iref in range(len(self.wlref_sat_indices)):  # range(0, self.n_satellite_bands):
    #         sat_band_index = self.wlref_sat_indices[iref]
    #         wl = self.satellite_bands[sat_band_index]
    #         # print(sat_band_index, wl)
    #         ins_band_index = np.argmin(np.abs(wl - self.insitu_bands))
    #         # print('------------------------------', sat_band_index, wl, ins_band_index,
    #         #       self.insitu_rrs[ins_band_index, self.ins_time_index])
    #
    #         if self.insitu_rrs.mask[ins_band_index, self.ins_time_index]:
    #             #print('me lo evalua aqui...')
    #             filtro_temporal = False
    #         if np.isnan(self.insitu_rrs[ins_band_index, self.ins_time_index]):
    #             #print('o mas bien es aqui...')
    #             filtro_temporal = False
    #             continue
    #         val = self.insitu_rrs[ins_band_index, self.ins_time_index]
    #         if val > 0.004 and sat_band_index == 2:
    #             filtro_temporal = False
    #         if val > 0.01 and sat_band_index == 4:
    #             filtro_temporal = False
    #         # if val > 0.008 and sat_band_index == 6:
    #         #     filtro_temporal = False
    #         # if val > 0.004 and sat_band_index == 7:
    #         #     filtro_temporal = False
    #         if val < 0:
    #             filtro_temporal = False

    # Getting spectra for comparison
    # self.mu_valid_bands = [False] * len(self.wlref_sat_indices)
    #
    # if time_condition and cond_min_valid_pxs:  # and filtro_temporal:
    #     self.mu_curr_ins_rrs = []
    #     self.mu_curr_sat_rrs_mean = []
    #     for iref in range(len(self.wlref_sat_indices)):
    #         sat_band_index = self.wlref_sat_indices[iref]
    #         wl = self.satellite_bands[sat_band_index]
    #         curr_sat_box = self.satellite_rrs[sat_band_index, r_s:r_e, c_s:c_e]
    #         ins_band_index = np.argmin(np.abs(wl - self.insitu_bands))
    #         difwl = abs(wl - self.insitu_bands[ins_band_index])
    #         check_dif_wl = True
    #         check_ins_value = True
    #         if difwl > 10:
    #             check_dif_wl = False
    #
    #         if self.insitu_rrs.mask[ins_band_index, self.ins_time_index]:
    #             check_ins_value = False
    #         if np.isnan(self.insitu_rrs[ins_band_index, self.ins_time_index]):
    #             check_ins_value = False
    #         if not np.isnan(self.insitu_rrs[ins_band_index, self.ins_time_index]) and self.insitu_rrs[
    #             ins_band_index, self.ins_time_index] < 0:
    #             check_ins_value = False
    #
    #         curr_sat_box = np.ma.masked_where(curr_sat_box == -999, curr_sat_box)
    #         curr_sat_box = np.ma.masked_invalid(curr_sat_box)
    #         curr_sat_box = np.ma.masked_array(curr_sat_box, mask)
    #         numberValid = np.ma.count(curr_sat_box)
    #         if float(self.valid_min_pixels) == 1:  # 100% like Zibordi
    #             cond_min_valid_pxs = numberValid == np.power(self.window_size, 2)
    #         else:
    #             cond_min_valid_pxs = numberValid > (float(self.valid_min_pixels) * NTP + 1)

    # if not check_ins_value:
    #     print(f'[WARNING] Wavelenght:  {wl} is not valid. Max diff: {difwl}')

    #         if cond_min_valid_pxs and check_dif_wl and check_ins_value:
    #             curr_sat_box_mean = curr_sat_box.mean()
    #             self.mu_curr_sat_rrs_mean.append(curr_sat_box_mean)
    #             ins_value = self.insitu_rrs[ins_band_index, self.ins_time_index]
    #             self.mu_curr_ins_rrs.append(ins_value)
    #             self.mu_valid_bands[iref] = True
    #             if iref == 1 and curr_sat_box_mean >= 0.0075:
    #                 self.mu_valid_bands[iref] = False
    #
    #     # print(self.mu_curr_sat_rrs_mean)
    #     # print(self.mu_curr_ins_rrs)
    # if sum(self.mu_valid_bands) == len(self.wlref_sat_indices):
    #     return True
    # else:
    #     return False

    def get_mu_key(self):
        sdate = self.mu_sat_time.strftime('%Y-%m-%d')
        idate = self.mu_insitu_time.strftime('%H:%M:%S')
        mukey = f'{sdate}_{idate}'
        return mukey

    def set_hour_sat_time(self, hour, minute):
        for index_mu in range(self.n_mu_total):
            sat_time_prev = self.sat_times[index_mu]
            sat_time_new = sat_time_prev.replace(hour=hour, minute=minute)
            self.sat_times[index_mu] = sat_time_new

    def prepare_df_validation(self):
        print('[INFO] Preparing DF for validation...')
        nbands = len(self.wlref)
        ntot = nbands * self.n_mu_total
        self.df_validation = pd.DataFrame(columns=self.col_names, index=list(range(ntot)))
        # ['Index','Index_MU','Index_Band','Sat_Time','Ins_Time','Time_Diff','Wavelenght','Ins_Rrs','Sat_Rrs','Valid']
        index_tot = 0
        # index_valid_tot = 0
        nmu_valid = 0
        nmu_valid_complete = 0
        for index_mu in range(self.n_mu_total):
            if index_mu % 100 == 0:
                print(f'[INFO] MU: {index_mu} of {self.n_mu_total}')
            # print(f'[INFO] MU: {index_mu} of {self.n_mu_total}')
            mu_valid, info_mu = self.load_mu_datav2(index_mu)

            # if not mu_valid:
            #     status = info_mu['status']
            #     print(f'[WARNING] MU: {index_mu} no valid: {status}')

            # invalid = [2,94,136,140,159]
            # if index_mu in invalid:
            #     mu_valid = False
            # if index_mu==94:
            #     mu_valid = False

            if mu_valid:
                nmu_valid = nmu_valid + 1

            spectrum_complete = info_mu['spectrum_complete']
            mu_valid_bands = info_mu['valid_bands']
            # n_good_bands = len(self.wlref_sat_indices)
            # if spectrum_complete:
            #     nmu_valid_complete = nmu_valid_complete + 1
            # else:
            #     n_good_bands = sum(mu_valid_bands)

            mukey = self.get_mu_key()
            time_diff = round(abs((self.mu_sat_time - self.mu_insitu_time).total_seconds() / 3600), 2)
            if not mukey in self.mu_dates.keys():
                self.mu_dates[mukey] = {
                    'Sat_Time': self.mu_sat_time.strftime('%Y-%m-%d %H:%M'),
                    'Ins_Time': self.mu_insitu_time.strftime('%Y-%m-%d %H:%M'),
                    'Time_Diff': time_diff,
                    'satellite': self.info['satellite'].upper(),
                    'platform': self.info['platform'].upper(),
                    'sensor': self.info['sensor'].upper(),
                    'site': self.info['insitu_site_name'].upper(),
                    'ac': self.info['satellite_aco_processor'].upper(),
                    'mu_valid': mu_valid,
                    'spectrum_complete': spectrum_complete,
                    'n_good_bands': 0
                }
                if self.mu_dates[mukey]['ac'] == 'ATMOSPHERIC CORRECTION PROCESSOR: XXX':
                    self.mu_dates[mukey]['ac'] = 'STANDARD'
            else:
                print('[WARNING] A single MDB file should not contain more than one match-up in a specific time/date')

            index_valid = 0
            n_good_bands = 0

            for iref in range(len(self.wlref_sat_indices)):

                sat_band_index = self.wlref_sat_indices[iref]
                valid_here = mu_valid
                if not spectrum_complete and mu_valid:
                    valid_here = mu_valid_bands[iref]
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
                    'Valid': [valid_here]  # [self.mu_valid_bands[sat_band_index]]
                }
                if valid_here:  # self.mu_valid_bands[sat_band_index]:
                    n_good_bands = n_good_bands + 1
                    row['Ins_Rrs'] = [self.mu_curr_ins_rrs[index_valid]]
                    row['Sat_Rrs'] = [self.mu_curr_sat_rrs_mean[index_valid]]
                    index_valid = index_valid + 1

                self.df_validation.iloc[index_tot] = pd.DataFrame.from_dict(row)

                index_tot = index_tot + 1

            self.mu_dates[mukey]['n_good_bands'] = n_good_bands
            spectrum_complete = n_good_bands == len(self.wlref_sat_indices)
            self.mu_dates[mukey]['spectrum_complete'] = spectrum_complete
            if spectrum_complete:
                nmu_valid_complete = nmu_valid_complete + 1

        self.df_validation_valid = self.df_validation[self.df_validation['Valid']][:]

        print(
            f'[INFO]# total match-ups: {self.n_mu_total} Valid: {nmu_valid}  With complete spectrum: {nmu_valid_complete}')
        return nmu_valid

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
        insitu_sensor = 'HYPSTAR'

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
        # self.set_wl_ref_insitu_indices()

    def set_wlsatrange_aswlref(self, wlmin, wlmax):
        self.wlref = []
        self.wlref_sat_indices = []
        for index in range(len(self.satellite_bands)):
            wl = self.satellite_bands[index]
            if wlmin <= wl <= wlmax:
                self.wlref.append(wl)
                self.wlref_sat_indices.append(index)

    def set_wlsatlist_aswlref(self, wlsatlist):
        wllist = self.qc_sat.get_wl_sat_list_from_wlreflist(wlsatlist)
        if len(wllist) == len(wlsatlist):
            self.set_wl_ref(wllist)

        # self.set_wl_ref_insitu_indices()

    # def set_wl_ref_insitu_indices(self):
    #     self.wlref_insitu_indices = []
    #     for wl in self.wlref:
    #         ins_index = np.argmin(np.abs(wl - self.insitu_bands))
    #         self.wlref_insitu_indices.append(ins_index)

    def get_nearestinsituwl_atsatwl(self, maxdiff, wlmin, wlmax):
        if wlmin is None:
            wlmin = min(self.satellite_bands)
        if wlmax is None:
            wlmax = max(self.satellite_bands)
        wllist = []
        for wl in self.satellite_bands:
            if wlmin <= wl <= wlmax:
                ins_index = np.argmin(np.abs(wl - self.insitu_bands))
                ins_wl = self.insitu_bands[ins_index]
                difwl = abs(wl - ins_wl)
                print(wl, ins_wl)
                if difwl <= maxdiff:
                    wllist.append(ins_wl)
        return wllist

    def get_sat_time_from_fname(self, fname):
        val_list = fname.split('_')
        sat_time = None
        for v in val_list:
            try:
                sat_time = datetime.strptime(v, '%Y%m%dT%H%M%S')
                break
            except ValueError:
                continue
        return sat_time

    def close(self):
        if self.VALID:
            self.nc.close()
