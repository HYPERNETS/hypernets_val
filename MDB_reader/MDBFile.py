import math
import os
import sys
from datetime import datetime
from datetime import timedelta

from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import pytz

code_home = os.path.dirname(os.path.dirname(__file__))
sys.path.append(code_home)
import COMMON.Class_Flags_OLCI as flag
from COMMON import common_functions as cfs

from QC_INSITU import QC_INSITU
from QC_SAT import QC_SAT

# from skimage import exposure
from matplotlib.colors import ListedColormap

DIMENSION_NAMES = ['satellite_id', 'rows', 'columns', 'satellite_bands', 'insitu_id', 'insitu_original_bands']
VAR_NAMES = ['satellite_time', 'satellite_latitude', 'satellite_longitude', 'satellite_bands',
             'satellite_Rrs', 'insitu_time', 'insitu_original_bands', 'insitu_Rrs', 'time_difference']

ATRIB_NAMES = ['creation_time', 'satellite', 'platform', 'sensor', 'description', 'satellite_aco_processor',
               'satellite_proc_version', 'site', 'in_situ_lat', 'in_situ_lon']


class MDBFile:

    def __init__(self, file_path):
        global DIMENSION_NAMES
        global VAR_NAMES
        global ATRIB_NAMES
        self.file_path = file_path
        self.VALID = True
        file_name = file_path.split('/')[-1]
        print(f'[INFO] Starting MDBFile: {file_name}')
        try:
            self.nc = Dataset(file_path)
            self.variables = self.nc.variables
            self.dimensions = self.nc.dimensions
            self.flag_band_name = 'satellite_WQSF'
            if self.check_structure() == 0:
                self.VALID = False

        except Exception as e:
            self.VALID = False
            print(f'[ERROR] Exception starting NetCDF file: {e}')
            return


        if not self.VALID:
            print(f'[ERROR] MDB File: {file_path} is not a valid MDB file with the correct structure.')
            return

        if self.VALID:
            self.n_mu_total = len(self.dimensions['satellite_id'])
            self.n_insitu_day = len(self.dimensions['insitu_id'])
            print('[INFO] Total mu: ', self.n_mu_total)
            self.sat_times = []
            for st in self.variables['satellite_time']:
                sat_time_here = datetime.utcfromtimestamp(float(st))
                if sat_time_here.hour == 0 and sat_time_here.minute == 0 and sat_time_here.second == 0:
                    sat_time_here = sat_time_here.replace(hour=11)
                self.sat_times.append(sat_time_here)
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
                if file_name.find('CCIALL') > 0:
                    self.info['satellite_aco_processor'] = 'CCIALL'
                else:
                    self.info['satellite_aco_processor'] = 'CCI'

            self.info['res'] = 'WFR'
            self.info['insitu_sensor'] = 'HYPSTAR'

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
        elif len(self.nc.satellite_aco_processor) == 0 or self.nc.satellite_aco_processor == 'CMEMS':
            self.qc_sat = QC_SAT(self.variables['satellite_Rrs'], self.satellite_bands, None,
                                 'Climate Change Initiative - European Space Agency')
        else:
            self.qc_sat = QC_SAT(self.variables['satellite_Rrs'], self.satellite_bands,
                                 self.variables[self.flag_band_name], self.info['satellite_aco_processor'])

        self.window_size = self.qc_sat.window_size

        # Variables to make validation...
        self.col_names = ['Index', 'Index_MU', 'Index_Band', 'Sat_Time', 'Ins_Time', 'Time_Diff', 'Wavelenght',
                          'Ins_Rrs',
                          'Sat_Rrs', 'Valid', 'status']
        self.df_validation = None
        self.df_validation_valid = None
        self.mu_dates = {}

        self.df_mu = None
        self.col_names_mu = ['Index_MU', 'Sat_Time', 'Ins_Time', 'Time_Diff', 'satellite', 'platform', 'sensor', 'site',
                             'ac', 'mu_valid', 'mu_insitu_id', 'spectrum_complete', 'n_good_bands', 'Ins_Lat',
                             'Ins_Lon', 'status']

        # variables controlling display images
        self.rgb_bands = [665, 560, 490]

        self.PI_DIVIDED = False
        self.var_mu_valid = 'mu_valid'
        self.allow_repeated = False

    def check_repeated(self):
        sat_times = np.array(self.sat_times)
        for idx in range(1, self.n_mu_total):
            sat_time_idx = sat_times[idx].replace(hour=0, minute=0, second=0, microsecond=0)
            for idu in range(idx):
                if sat_time_idx == sat_times[idu].replace(hour=0, minute=0, second=0, microsecond=0):
                    print(f'[WARNING] There are repeated satellite dates.')
                    return False
        return True

    def check_wl(self):
        if 'mu_wavelength' not in self.variables:
            return True

        checkw = True
        mu_wl = np.unique(np.array(self.variables['mu_wavelength'][:]))
        wl_list = list(self.satellite_bands)
        for mwl in mu_wl:
            if mwl not in wl_list:
                print(f'[WARNING] {mwl} is not available in satellite bands variable')
                checkw = False

        return checkw

    def correct_wl(self):
        print('[INFO] Correcting wavelengths... ')
        mu_wl = np.array(self.variables['mu_wavelength'][:])
        mu_wl_end = mu_wl.copy()
        wl_list = list(self.satellite_bands)
        nreassigned = 0
        for idx in range(len(mu_wl)):
            mwl = mu_wl[idx]
            if mwl in wl_list:
                continue
            diffwl = np.abs(mwl - self.satellite_bands)
            imin = np.argmin(diffwl)
            if diffwl[imin] <= 1:
                mu_wl_end[idx] = self.satellite_bands[imin]
                nreassigned = nreassigned + 1
                # print(f'[INFO] Wavelength: {mwl} set to {mu_wl_end[idx]}')
            else:
                print(f'[ERROR] Wavelength could not be reassigned')
                return

        print(f'[INFO] Number of reassigned data points: {nreassigned}')

        print(f'[INFO] Creating temporary file with corrected wavelengths...')
        file_temp = os.path.join(os.path.dirname(self.file_path), 'Temp.nc')
        if os.path.exists(file_temp):
            os.remove(file_temp)

        ncout = Dataset(file_temp, 'w', format='NETCDF4')

        # copy global attributes all at once via dictionary
        ncout.setncatts(self.nc.__dict__)

        # copy dimensions (satellite_id is defined as unlimited, so we do not need to change it)
        for name, dimension in self.nc.dimensions.items():
            ncout.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))

        for name, variable in self.nc.variables.items():
            fill_value = None
            if '_FillValue' in variable.ncattrs():
                fill_value = variable._FillValue
            ncout.createVariable(name, variable.datatype, variable.dimensions, zlib=True, fill_value=fill_value,
                                 shuffle=True, complevel=6)
            ncout[name].setncatts(self.nc[name].__dict__)

            if name == 'mu_wavelength':
                ncout[name][:] = mu_wl_end[:]
            else:
                ncout[name][:] = self.nc[name][:]
        ncout.close()

        self.nc.close()
        print(f'[INFO] Final file: {self.file_path}')
        os.rename(file_temp, self.file_path)
        print(f'[INFO] Completed')

    def remove_repeated(self):

        sat_times = np.array(self.sat_times)

        idx_included = [True] * self.n_mu_total
        n_excluded = 0
        print(f'[INFO] Searching repeated ids...')

        for idx in range(1, self.n_mu_total):
            sat_time_idx = sat_times[idx].replace(hour=0, minute=0, second=0, microsecond=0)
            for idu in range(idx):
                if sat_time_idx == sat_times[idu].replace(hour=0, minute=0, second=0, microsecond=0):
                    idx_included[idx] = False
                    n_excluded = n_excluded + 1
                    print(
                        f'[INFO] Repeated satellite id: {idx} with time: {sat_times[idx]} is equal to {idu} -> {sat_times[idu]}')

        print(f'[INFO] Number of repeated ids: {n_excluded}')

        indices = np.array(idx_included, dtype=np.bool)

        print(f'[INFO] Creating temporary file without repeated ids...')
        file_temp = os.path.join(os.path.dirname(self.file_path), 'Temp.nc')
        if os.path.exists(file_temp):
            os.remove(file_temp)

        ncout = Dataset(file_temp, 'w', format='NETCDF4')

        # copy global attributes all at once via dictionary
        ncout.setncatts(self.nc.__dict__)

        # copy dimensions (satellite_id is defined as unlimited, so we do not need to change it)
        for name, dimension in self.nc.dimensions.items():
            ncout.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))

        for name, variable in self.nc.variables.items():
            fill_value = None
            if '_FillValue' in variable.ncattrs():
                fill_value = variable._FillValue
            ncout.createVariable(name, variable.datatype, variable.dimensions, zlib=True, fill_value=fill_value,
                                 shuffle=True, complevel=6)
            ncout[name].setncatts(self.nc[name].__dict__)

            if 'satellite_id' in variable.dimensions:
                var_new_array = np.array(variable[indices == True])
                nnew = var_new_array.shape[0]
                # print('Longitud de nueva variable: ',name,len(indices[0]))
                for idx in range(nnew):
                    if name == 'satellite_PDU':
                        ncout.variables[name][idx] = str(var_new_array[idx])
                    else:
                        ncout.variables[name][idx] = var_new_array[idx]
            else:
                ncout[name][:] = self.nc[name][:]
        ncout.close()

        self.nc.close()
        print(f'[INFO] Final file: {self.file_path}')
        os.rename(file_temp, self.file_path)
        print(f'[INFO] Completed')

    def create_file_with_flag_bands(self, file_out):

        ncout = Dataset(file_out, 'w', format='NETCDF4')

        # copy global attributes all at once via dictionary
        ncout.setncatts(self.nc.__dict__)

        # copy dimensions (
        for name, dimension in self.nc.dimensions.items():
            ncout.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        # copy variables
        for name, variable in self.nc.variables.items():
            fill_value = None
            if '_FillValue' in variable.ncattrs():
                fill_value = variable._FillValue
            ncout.createVariable(name, variable.datatype, variable.dimensions, zlib=True, fill_value=fill_value,
                                 shuffle=True, complevel=6)
            ncout[name].setncatts(self.nc[name].__dict__)
            ncout[name][:] = self.nc[name][:]
        ncout.close()

        # adding flag variables
        name_flag_names = ['flag_platform', 'flag_sensor', 'flag_ac', 'flag_insitu']

        self.nc.close()

    def check_structure(self):
        
        for var in VAR_NAMES:
            if var not in self.variables:
                print(f'[ERROR] Variable: {var} is not available in MDB file')
                return 0
        for dim in DIMENSION_NAMES:
            if dim not in self.dimensions:
                print(f'[ERROR] Dimension: {dim} is not available in MDB file')
                return 0

        check_atrib = True
        for atrib in ATRIB_NAMES:
            # if atrib == 'site' or atrib == 'in_situ_lat' or atrib == 'in_situ_lon':
            #     continue
            if atrib not in self.nc.ncattrs():
                print(f'[WARNING] Attribute: {atrib} is not available in MDB file')
                check_atrib = False


        if not self.flag_band_name in self.variables:
            if self.nc.satellite_aco_processor.upper() == 'POLYMER' and 'satellite_bitmask' in self.variables:
                self.flag_band_name = 'satellite_bitmask'

            if self.nc.satellite_aco_processor.upper() == 'C2RCC' and 'satellite_c2rcc_flags' in self.variables:
                self.flag_band_name = 'satellite_c2rcc_flags'

            if self.nc.satellite_aco_processor.upper() == 'FUB' and 'satellite_quality_flags':
                self.flag_band_name = 'satellite_quality_flags'

            if self.nc.satellite_aco_processor.upper() == 'ACOLITE':
                self.flag_band_name = 'NONE'

        if check_atrib:
            return 1
        else:
            return 2

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

        # CHECK IF EXIST SPATIAL INDEX
        if 'insitu_spatial_index' in self.variables:
            self.qc_insitu.spatial_index_array = self.variables['insitu_spatial_index'][index_mu]
        else:
            self.qc_insitu.spatial_index_array = None

        for idx in range(len(times_here)):
            itime = times_here[idx]

            if not np.ma.is_masked(itime) and not np.isnan(itime):
                insitu_time_here = datetime.utcfromtimestamp(float(itime))
                time_diff_here = abs((sat_time_here - insitu_time_here).total_seconds())
                time_difference[idx] = time_diff_here
                diff_time = abs(time_difference_prev[idx] - time_diff_here)
                if diff_time >= 150:
                    print('[WARNING] Unexplained time difference. Please review sat and insitu times are in utc: ')
                    ssource = self.variables['satellite_PDU'][index_mu]
                    isource = self.variables['insitu_filename'][index_mu][idx]
                    print(f'--> Satellite time in file: {sat_time_here} Name file: {ssource}')
                    print(f'--> In situ time in file: {insitu_time_here} Name file: {isource}')
                    print(f'--> Time difference in file: {time_difference_prev[idx]} Recomputed: {time_diff_here}')

            else:
                time_difference[idx] = np.ma.masked

        if 'insitu_exact_wavelenghts' in self.variables:
            exact_wl = self.variables['insitu_exact_wavelenghts'][index_mu]
        else:
            exact_wl = self.variables['insitu_original_bands']

        ins_time_index, time_condition, valid_insitu, spectrum_complete, rrs_values = \
            self.qc_insitu.get_finalspectrum_mu(index_mu, time_difference, exact_wl, self.wlref)

        # print(ins_time_index, time_condition,valid_insitu)

        if time_condition and valid_insitu:
            ins_time = self.variables['insitu_time'][index_mu][ins_time_index]
            mu_insitu_time = datetime.utcfromtimestamp(float(ins_time))
        else:  ##aunque los datos sean invalidos (time dif>max time dif), obtenemos el mu_insitu_time como referencia

            ins_time_index = np.argmin(np.abs(time_difference))
            ins_time = self.variables['insitu_time'][index_mu][ins_time_index]
            # s_index = spatial_index[ins_time_index]
            # print('Me llega aqui con index_mu igual a: ', index_mu, 'ins time: ',ins_time, 's_index must be zero: ',s_index)
            if np.ma.is_masked(ins_time):  ##all the ins situ time is masked,we can do mu_insitu_time==mu_sat_time
                mu_insitu_time = sat_time_here
            else:
                mu_insitu_time = datetime.utcfromtimestamp(float(ins_time))

        # if index_mu==0:
        #     print('----------->', mu_insitu_time)

        return ins_time_index, mu_insitu_time, time_condition, valid_insitu, spectrum_complete, rrs_values

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
        name_variable_insitu = 'insitu_Rrs'
        if not self.qc_insitu.apply_nir_correction:
            name_variable_insitu = 'insitu_Rrs_nosc'

        self.insitu_rrs = self.variables[name_variable_insitu][index_mu]

        self.satellite_rrs = self.variables['satellite_Rrs'][index_mu]

        # Sat and instrument time
        self.mu_sat_time = self.sat_times[index_mu]

        # THIS STEP IS NOW DONE BEFORE PREPARING DF FOR VALIDATION
        # if self.info['satellite_aco_processor'] == 'CCI':
        #     self.mu_sat_time = self.mu_sat_time.replace(hour=11)

        # print(index_mu)
        # if index_mu==95 or index_mu==181:
        #     self.mu_curr_ins_rrs = []
        #     self.mu_curr_sat_rrs_mean = []
        #     load_info['status'] = -3  # f'IN SITU DATA OUT OF TIME WINDOW'
        #     return is_mu_valid, load_info

        self.ins_time_index, self.mu_insitu_time, time_condition, valid_insitu, spectrum_complete, rrs_ins_values = \
            self.retrieve_ins_info_mu_spectra(index_mu)

        # print(self.valid_insitu)

        if rrs_ins_values is not None and self.PI_DIVIDED:
            rrs_ins_values = rrs_ins_values / np.pi

        load_info['spectrum_complete'] = spectrum_complete

        if not time_condition:
            self.mu_curr_ins_rrs = []
            self.mu_curr_sat_rrs_mean = []
            load_info['status'] = -3  # f'IN SITU DATA OUT OF TIME WINDOW'

            # print('---out of time error------------')
            # print(self.ins_time_index)
            # print(self.mu_insitu_time, ' versus ', self.mu_sat_time)
            # print(self.qc_insitu.time_max)

            # print('???????????????????????')

            return is_mu_valid, load_info

        if not valid_insitu:
            self.mu_curr_ins_rrs = []
            self.mu_curr_sat_rrs_mean = []
            load_info['status'] = -4  # f'INVALID INSITU DATA'
            return is_mu_valid, load_info

        if not spectrum_complete and self.qc_insitu.only_complete_spectra:
            self.mu_curr_ins_rrs = []
            self.mu_curr_sat_rrs_mean = []
            load_info['status'] = -5  # f'INCOMPLETE IN SITU SPECTRUM'
            return is_mu_valid, load_info

        cond_min_pixels, cond_stats, valid_mu, sat_values = self.qc_sat.get_match_up_values(index_mu)

        if not valid_mu:
            self.mu_curr_ins_rrs = []
            self.mu_curr_sat_rrs_mean = []
            load_info['status'] = -6  # f'NO VALID SAT DATA'
            return is_mu_valid, load_info

        # Getting spectra for comparison
        mu_valid_bands = [False] * len(self.wlref_sat_indices)
        self.mu_curr_ins_rrs = []
        self.mu_curr_sat_rrs_mean = []
        for iref in range(len(self.wlref_sat_indices)):
            if not np.ma.is_masked(rrs_ins_values[iref]):
                iband = self.wlref_sat_indices[iref]
                sat_value_here = sat_values[iref]
                if self.qc_sat.pi_multiplied[iband]:
                    sat_value_here = sat_value_here * np.pi
                if self.qc_sat.pi_divided[iband]:
                    sat_value_here = sat_value_here / np.pi

                self.mu_curr_sat_rrs_mean.append(sat_value_here)
                self.mu_curr_ins_rrs.append(rrs_ins_values[iref])
                mu_valid_bands[iref] = True

        if sum(mu_valid_bands) == 0:
            load_info['status'] = -7  # f'NO VALID IN SITU DATA'
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

    def save_flag_images(self, path_img):
        print('[INFO] Saving FLAG images...')
        for index_mu in range(self.n_mu_total):
            flagging = flag.Class_Flags_OLCI(self.variables[self.flag_band_name].flag_masks,
                                             self.variables[self.flag_band_name].flag_meanings)
            satellite_flag_band = np.array(self.variables[self.flag_band_name][index_mu])

            mask = np.zeros((25, 25))

            # clouds
            flag_clouds = ['CLOUD', 'CLOUD_AMBIGUOUS', 'CLOUD_MARGIN']
            mask_clouds = flagging.Mask(satellite_flag_band, flag_clouds)
            mask[np.where(mask_clouds != 0)] = 1

            # other
            flag_others = ['INVALID', 'COSMETIC', 'SATURATED', 'SUSPECT', 'HISOLZEN', 'HIGHGLINT', 'SNOW_ICE',
                           'AC_FAIL', 'WHITECAPS']
            mask_others = flagging.Mask(satellite_flag_band, flag_others)
            mask[np.where(mask_others != 0)] = 2

            # negative reflectance
            flag_rwneg = ['RWNEG_O2', 'RWNEG_O3', 'RWNEG_O4', 'RWNEG_O5', 'RWNEG_O6', 'RWNEG_O7', 'RWNEG_O8']
            mask_rwneg = flagging.Mask(satellite_flag_band, flag_rwneg)
            mask[np.where(mask_rwneg != 0)] = 3

            # land
            flag_land = ['LAND', 'COASTLINE']
            mask_land = flagging.Mask(satellite_flag_band, flag_land)
            mask[np.where(mask_land != 0)] = 4

            mu_sat_time = self.sat_times[index_mu]
            mu_sat_time_str = mu_sat_time.strftime('%Y%m%d_%H%M')
            file_name = f'SATFLAG_{mu_sat_time_str}.png'
            file_img = os.path.join(path_img, file_name)
            print(file_img)
            extent = [0, 25, 0, 25]
            cmap_here = ListedColormap(['b', 'w', 'red', 'orange', 'maroon'])

            plt.imshow(mask, interpolation=None, extent=extent, cmap=cmap_here, vmin=0, vmax=5)

            plt.hlines(11, 11, 14, colors=['r'])
            plt.hlines(14, 11, 14, colors=['r'])
            plt.vlines(11, 11, 14, colors=['r'])
            plt.vlines(14, 11, 14, colors=['r'])
            plt.colorbar()
            plt.savefig(file_img)
            plt.close()

    def save_rgb_images(self, path_img):
        print('[INFO] Saving RGB images...')
        for index_mu in range(self.n_mu_total):
            satellite_rrs = self.variables['satellite_Rrs'][index_mu]

            indexred = np.argmin(np.abs(self.rgb_bands[0] - self.satellite_bands))
            indexgreen = np.argmin(np.abs(self.rgb_bands[1] - self.satellite_bands))
            indexblue = np.argmin(np.abs(self.rgb_bands[2] - self.satellite_bands))

            print(indexred, indexgreen, indexblue)

            band_red = np.ma.array(satellite_rrs[indexred][:][:])
            band_green = np.ma.array(satellite_rrs[indexgreen][:][:])
            band_blue = np.ma.array(satellite_rrs[indexblue][:][:])

            rgb = np.ma.zeros((25, 25, 3))

            rgb[:, :, 0] = self.scale_array(band_red, None, None)
            rgb[:, :, 1] = self.scale_array(band_green, None, None)
            rgb[:, :, 2] = self.scale_array(band_blue, None, None)
            # rgbi = rgb * 255
            # rgbi = rgb.astype(int)

            # rgb[rgb.mask] = 1
            # rgb = exposure.equalize_adapthist(rgb)

            mu_sat_time = self.sat_times[index_mu]
            mu_sat_time_str = mu_sat_time.strftime('%Y%m%d_%H%M')
            file_name = f'SATRGB_{mu_sat_time_str}.png'
            file_img = os.path.join(path_img, file_name)
            print(file_img)
            extent = [0, 25, 0, 25]
            plt.imshow(rgb, interpolation=None, extent=extent)
            # plt.xticks(np.arange(2, 25, 3))
            # plt.yticks(np.arange(2, 25, 3))
            # plt.grid()
            # plt.grid(b=True, which='major', color='#666666', linestyle='-')
            # plt.minorticks_on()
            # plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
            plt.hlines(11, 11, 14, colors=['r'])
            plt.hlines(14, 11, 14, colors=['r'])
            plt.vlines(11, 11, 14, colors=['r'])
            plt.vlines(14, 11, 14, colors=['r'])
            # plt.colorbar()
            plt.savefig(file_img)
            plt.close()

    def save_wl_images(self, path_img):
        print('[INFO] Saving RGB images...')
        for index_mu in range(self.n_mu_total):
            if index_mu != 15:
                continue
            satellite_rrs = self.variables['satellite_Rrs'][index_mu]
            for idx in range(14):
                wl = self.satellite_bands[idx]
                wls = str(wl).replace('.', '_')
                band = np.ma.array(satellite_rrs[idx][:][:])
                mu_sat_time = self.sat_times[index_mu]
                mu_sat_time_str = mu_sat_time.strftime('%Y%m%d_%H%M')
                file_name = f'SAT_{index_mu}_{wls}_{mu_sat_time_str}.png'
                file_img = os.path.join(path_img, file_name)
                print(file_img)
                extent = [0, 25, 0, 25]
                cmap = plt.colormaps['gray']
                plt.imshow(band, cmap=cmap, interpolation=None, extent=extent)

                # plt.hlines(11, 11, 14, colors=['r'])
                # plt.hlines(14, 11, 14, colors=['r'])
                # plt.vlines(11, 11, 14, colors=['r'])
                # plt.vlines(14, 11, 14, colors=['r'])
                # plt.colorbar()
                plt.savefig(file_img)
                plt.close()
        for index_mu in range(self.n_mu_total):
            name_bands = ['satellite_WQSF', 'satellite_AOT_0865p50', 'OAA', 'OZA']
            for name_band in name_bands:
                satellite_band = self.variables['satellite_WQSF'][index_mu]
                if index_mu != 15:
                    continue
                band = satellite_band[:][:]
                file_name = f'SAT_{index_mu}_{name_band}.png'
                file_img = os.path.join(path_img, file_name)
                print(file_img)
                extent = [0, 25, 0, 25]
                if name_band == 'satellite_WQSF':
                    cmap = plt.colormaps['Set1']
                else:
                    cmap = plt.colormaps['gray']
                plt.imshow(band, cmap=cmap, interpolation=None, extent=extent)
                plt.savefig(file_img)
                plt.close()

    def scale_array(self, array, min_value, max_value):
        # print('=====================================')
        # print(array)
        if min_value is None:
            min_value = np.min(array)
        if max_value is None:
            max_value = np.max(array)
        array_out = (array - min_value) / (max_value - min_value)

        return array_out

    def prepare_df_validation(self):
        print('[INFO] Preparing DF for validation...')
        nbands = len(self.wlref)
        ntot = nbands * self.n_mu_total
        self.df_validation = pd.DataFrame(columns=self.col_names, index=list(range(ntot)))
        # ['Index','Index_MU','Index_Band','Sat_Time','Ins_Time','Time_Diff','Wavelenght','Ins_Rrs','Sat_Rrs','Valid']
        index_tot = 0
        nmu_valid = 0
        nmu_valid_complete = 0
        status_error = [0] * 8
        for index_mu in range(self.n_mu_total):
            if index_mu % 100 == 0:
                print(f'[INFO] MU: {index_mu} of {self.n_mu_total}')

            mu_valid, info_mu = self.load_mu_datav2(index_mu)

            if info_mu['status'] < 0:
                spos = info_mu['status'] * (-1)
                status_error[spos] = status_error[spos] + 1

            # print(index_mu)
            # print(info_mu)
            # print('------------')

            if mu_valid:
                nmu_valid = nmu_valid + 1

            spectrum_complete = info_mu['spectrum_complete']
            mu_valid_bands = info_mu['valid_bands']

            mukey = self.get_mu_key()
            time_diff = round(abs((self.mu_sat_time - self.mu_insitu_time).total_seconds() / 3600), 2)
            ins_lat = -999.0
            ins_lon = -999.0
            if 'insitu_latitude' in self.variables and 'insitu_longitude' in self.variables:
                ins_lat = float(self.variables['insitu_latitude'][index_mu, self.ins_time_index])
                ins_lon = float(self.variables['insitu_longitude'][index_mu, self.ins_time_index])

            # if index_mu==0:
            #     print(self.mu_insitu_time,self.mu_sat_time,time_diff)
            # compability with old mdb
            key_site = 'site'
            if key_site not in self.info.keys() and 'insitu_site_name' in self.info.keys():
                key_site = 'insitu_site_name'

            if mukey in self.mu_dates.keys() and self.allow_repeated:
                index = 1
                while index >= 1:
                    mukey = f'{mukey}_{index}'
                    if mukey in self.mu_dates.keys():
                        index = index + 1
                    else:
                        break

            if not mukey in self.mu_dates.keys():
                self.mu_dates[mukey] = {
                    'Index_MU': index_mu,
                    'Sat_Time': self.mu_sat_time.strftime('%Y-%m-%d %H:%M'),
                    'Ins_Time': self.mu_insitu_time.strftime('%Y-%m-%d %H:%M'),
                    'Ins_Lat': ins_lat,
                    'Ins_Lon': ins_lon,
                    'Time_Diff': time_diff,
                    'satellite': self.info['satellite'].upper(),
                    'platform': self.info['platform'].upper(),
                    'sensor': self.info['sensor'].upper(),
                    'site': self.info[key_site].upper(),
                    'ac': self.info['satellite_aco_processor'].upper(),
                    'mu_valid': mu_valid,
                    'mu_insitu_id': self.ins_time_index,
                    'spectrum_complete': spectrum_complete,
                    'n_good_bands': 0,
                }
                if self.mu_dates[mukey]['ac'] == 'ATMOSPHERIC CORRECTION PROCESSOR: XXX':
                    self.mu_dates[mukey]['ac'] = 'STANDARD'
            else:
                print('[WARNING] A single MDB file should not contain more than one match-up in a specific time/date')
                print(f'[WARNING] REPEATED MU KEY: {mukey}')

            index_valid = 0
            n_good_bands = 0

            for iref in range(len(self.wlref_sat_indices)):

                sat_band_index = self.wlref_sat_indices[iref]
                valid_here = mu_valid

                if not spectrum_complete and mu_valid:
                    valid_here = mu_valid_bands[iref]

                # valid_here_int = 1
                # if not valid_here:
                #     valid_here_int = 0
                #     spos = info_mu['status'] * (-1)
                #     if spos<=3:
                #         valid_here_int = -1

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
                    'Valid': [valid_here],
                    'status': [info_mu['status']]
                }
                if valid_here:  # self.mu_valid_bands[sat_band_index]:
                    n_good_bands = n_good_bands + 1
                    row['Ins_Rrs'] = [self.mu_curr_ins_rrs[index_valid]]
                    row['Sat_Rrs'] = [self.mu_curr_sat_rrs_mean[index_valid]]
                    index_valid = index_valid + 1

                # print('Llega aqui-> ',index_tot,pd.__version__)
                # print(pd.DataFrame.from_dict(row))
                # print(self.df_validation.columns.values)
                row_here = pd.DataFrame.from_dict(row)

                # print(type(row_here))
                self.df_validation.iloc[index_tot] = row_here.iloc[0]
                # print(self.df_validation.index[index_tot])

                index_tot = index_tot + 1

            self.mu_dates[mukey]['n_good_bands'] = n_good_bands
            spectrum_complete = n_good_bands == len(self.wlref_sat_indices)
            self.mu_dates[mukey]['spectrum_complete'] = spectrum_complete
            if spectrum_complete:
                nmu_valid_complete = nmu_valid_complete + 1

        self.df_validation_valid = self.df_validation[self.df_validation['Valid']][:]

        self.prepare_df_mu()
        print(
            f'[INFO] #total match-ups: {self.n_mu_total} Valid: {nmu_valid}  With complete spectrum: {nmu_valid_complete}')

        valid_options = ['Satellite', 'Platform', 'Sensor', 'Site', 'Total', 'Valid', 'NO VALID MDB FILE',
                         'NO VALID MATCH-UP INDEX', 'OUT OF TIME WINDOW', 'INVALID IN SITU DATA',
                         'INCOMPLETE IN SITU SPECTRUM', 'INVALID SAT DATA', 'INVALID IN SITU VALUES']

        df_valid = pd.DataFrame(index=list(range(1)), columns=valid_options)
        df_valid.iloc[0].at['Satellite'] = self.nc.satellite
        df_valid.iloc[0].at['Platform'] = self.nc.platform
        df_valid.iloc[0].at['Sensor'] = self.nc.sensor
        key_site = 'site'
        if key_site not in self.nc.ncattrs() and 'insitu_site_name' in self.nc.ncattrs():
            key_site = 'insitu_site_name'
        print(key_site)
        df_valid.iloc[0].at['Site'] = self.nc.getncattr(key_site)
        df_valid.iloc[0].at['Total'] = self.n_mu_total
        df_valid.iloc[0].at['Valid'] = nmu_valid

        for idx in range(1, len(status_error)):
            idxl = idx + 5
            col_name = valid_options[idxl]
            sbase = f'# Invalid match-ups [{col_name}]: '
            df_valid.iloc[0].at[col_name] = status_error[idx]
            # if idx == 1:
            #     sbase = '# Invalid match-ups [NO VALID MDB FILE]: '
            # if idx == 2:
            #     sbase = '# Invalid match-ups [NO VALID MATCH-UP INDEX]: '
            # if idx == 3:
            #     sbase = '# Invalid match-ups [OUT OF TIME WINDOW]: '
            # if idx == 4:
            #         sbase = '# Invalid match-ups [INVALID IN SITU DATA]: '
            # if idx == 5:
            #         sbase = '# Invalid match-ups [INCOMPLETE IN SITU SPECTRUM]: '
            # if idx == 6:
            #         sbase = '# Invalid match-ups [INVALID SAT DATA]: '
            # if idx == 7:
            #         sbase = '# Invalid match-ups [INVALID IN SITU VALUES]: '
            if status_error[idx] > 0:
                print(f'[INFO] -> {sbase}{status_error[idx]}')

        return nmu_valid, df_valid

    def prepare_df_mu(self):
        self.df_mu = pd.DataFrame(columns=self.col_names_mu, index=list(range(len(self.mu_dates))))

        index_tot = 0
        for mukey in self.mu_dates:
            row = self.mu_dates[mukey]
            for c in row:
                self.df_mu.iloc[index_tot].at[c] = row[c]  # pd.DataFrame.from_dict(row)
            index_tot = index_tot + 1

    ##PLOT FUNCTIONS
    def plot_spectra(self, hfig):
        if len(self.mu_curr_sat_rrs_mean) == 0 or len(self.mu_curr_ins_rrs) == 0:
            return

        print('Plot spectra...')
        legend = [f"S3{self.nc.platform}-{self.mu_sat_time.strftime('%d-%m-%Y %H:%M')}",
                  f"Hypstar-{self.mu_insitu_time.strftime('%d-%m-%Y %H:%M')}"]
        row_names = []
        for w in self.wlref:
            wn = f'{w}'
            row_names.append(wn)
        ydata = np.array([self.mu_curr_sat_rrs_mean, self.mu_curr_ins_rrs])
        df = pd.DataFrame(np.transpose(ydata), columns=legend, index=row_names)
        title = f"Match-up #{self.index_mu} {self.mu_sat_time.strftime('%d-%m-%Y')}"

        if hfig is None:
            hfig = plt.figure()
        df.plot(lw=2, marker='.', markersize=10)
        rindex = list(range(len(self.wlref)))
        plt.xticks(rindex, row_names, rotation=45, fontsize=12)
        plt.xlabel('Wavelength(nm)', fontsize=14)
        plt.ylabel('R$_{rs}$[1/sr]', fontsize=14)
        plt.title(title, fontsize=16)
        plt.gcf().subplots_adjust(bottom=0.18)
        plt.gcf().subplots_adjust(left=0.20)
        plt.gcf().canvas.draw()

        ofname = f"Spectra_{self.index_mu}.jpg"
        path_base = os.path.dirname(self.file_path)
        name_mdb = os.path.basename(self.file_path)
        path_out = os.path.join(path_base, f'{name_mdb[:-3]}')
        if not os.path.exists(path_out):
            os.mkdir(path_out)
        path_out_mu = os.path.join(path_out, 'MUPlots')
        if not os.path.exists(path_out_mu):
            os.mkdir(path_out_mu)
        ofname = os.path.join(path_out_mu, ofname)
        print('Name out:', ofname)
        plt.savefig(ofname, dpi=300)
        plt.close()

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
        # res = 'WFR'
        res = self.info['res']
        # insitu_sensor = 'HYPSTAR'
        insitu_sensor = self.info['insitu_sensor']

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

    def check_bands(self, wllist, maxdiff):

        for wl in wllist:
            index = np.argmin(np.abs(wl - self.satellite_bands))
            diff = np.abs(wl - self.satellite_bands[index])
            if diff > maxdiff:
                print(f'[ERROR] Band {wl} is not avaiable in satellite band list:')
                print(f'[ERROR] -> {self.satellite_bands}')
                return -1
        if len(wllist) != len(self.satellite_bands):
            return 0

        return 1

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

    def get_insitu_wl(self):
        wl = np.array(self.nc.variables['insitu_original_bands'][:])
        return wl

    def get_sat_wl_as_strlist(self, wllist):
        wlstrlist = []
        wlsat_orig = np.array(self.nc.variables['satellite_bands'])
        for wl in wllist:
            isat = np.argmin(np.abs(wlsat_orig - wl))
            wlstr = f'{wlsat_orig[isat]:.2f}'
            wlstr = wlstr.replace('.00', '')
            wlstrlist.append(wlstr)
        return wlstrlist

    def get_mu_insitu_spectra(self, index_mu, scale_factor):
        import numpy.ma as ma
        noriginal_bands = len(self.dimensions['insitu_original_bands'])
        spectra_all = ma.zeros((self.n_insitu_day, noriginal_bands))
        var_insitu = self.nc.variables['insitu_Rrs']
        insitu_valid = np.array(self.nc.variables['insitu_valid'][index_mu])
        insitu_id = self.nc.variables['mu_insitu_id'][index_mu]
        insitu_valid[insitu_id] = 2
        for idx in range(self.n_insitu_day):
            spectra_here = ma.array(var_insitu[index_mu, :, idx]).transpose()
            if np.sum(np.isnan(spectra_here)) > 0:
                insitu_valid[idx] = -1
            else:
                if scale_factor is not None:
                    spectra_here = spectra_here * scale_factor
                spectra_all[idx, :] = spectra_here[:]

        spectra_selected = spectra_all[insitu_valid == 2]
        spectra_valid = spectra_all[insitu_valid == 1]
        spectra_invalid = spectra_all[insitu_valid == 0]

        return spectra_selected, spectra_valid, spectra_invalid

    def get_all_insitu_spectra(self, scale_factor, use_rhow, compute_stats):
        import numpy.ma as ma
        nspectra = self.n_mu_total * self.n_insitu_day
        noriginal_bands = len(self.dimensions['insitu_original_bands'])
        all_spectra = ma.zeros((nspectra, noriginal_bands))
        all_spectra_validity = np.zeros((nspectra,))
        var_insitu = self.nc.variables['insitu_Rrs']
        var_insitu_valid = self.nc.variables['insitu_valid']
        index_row = 0
        for index_mu in range(self.n_mu_total):
            if (index_mu % 500) == 0:
                print(f'[INFO] Getting spectra {index_mu} / {self.n_mu_total}')
            idx_selected = -1
            if self.nc.variables['mu_valid'][index_mu] == 1:
                idx_selected = self.nc.variables['mu_insitu_id'][index_mu]
            for idx in range(self.n_insitu_day):
                if np.ma.is_masked(self.nc.variables['insitu_time'][index_mu, idx]):
                    continue
                all_spectra_validity[index_row] = var_insitu_valid[index_mu, idx]
                if idx == idx_selected:
                    all_spectra_validity[index_row] = 2
                spectra_here = ma.array(var_insitu[index_mu, :, idx]).transpose()
                if use_rhow:
                    spectra_here = spectra_here * np.pi
                if scale_factor is not None:
                    spectra_here = spectra_here * scale_factor
                all_spectra[index_row, :] = spectra_here[:]
                index_row = index_row + 1

        all_spectra = all_spectra[0:index_row, :]
        all_spectra_validity = all_spectra_validity[0:index_row]

        spectra_stats = None
        if compute_stats:
            spectra_good = all_spectra[all_spectra_validity >= 1]
            import statistics as st
            spectra_avg = ma.mean(spectra_good, axis=0)
            spectra_std = ma.std(spectra_good, axis=0)
            indices_max = ma.argmax(spectra_good, axis=0)
            imax = st.mode(indices_max)
            spectra_max_real = spectra_good[imax, :]
            spectra_max = ma.max(spectra_good, axis=0)

            indices_min = ma.argmin(spectra_good, axis=0)
            imin = st.mode(indices_min)
            spectra_mim_real = spectra_good[imin, :]
            spectra_min = ma.min(spectra_good, axis=0)

            spectra_stats = {
                'avg': spectra_avg,
                'std': spectra_std,
                'spectra_min_real': spectra_mim_real,
                'spectra_max_real': spectra_max_real,
                'spectra_min': spectra_min,
                'spectra_max': spectra_max,
            }

        return all_spectra, all_spectra_validity, spectra_stats

    def get_all_insitu_valid_spectra(self, scale_factor):
        import numpy.ma as ma
        nspectra = self.n_mu_total * self.n_insitu_day
        noriginal_bands = len(self.dimensions['insitu_original_bands'])
        spectra_good = ma.zeros((nspectra, noriginal_bands))
        var_insitu = self.nc.variables['insitu_Rrs']
        var_insitu_valid = self.nc.variables['insitu_valid']
        index_row = 0
        for index_mu in range(self.n_mu_total):
            for idx in range(self.n_insitu_day):
                if var_insitu_valid[index_mu, idx] == 1:
                    spectra_here = ma.array(var_insitu[index_mu, :, idx]).transpose()
                    if scale_factor is not None:
                        spectra_here = spectra_here * scale_factor
                    spectra_good[index_row, :] = spectra_here[:]
                    index_row = index_row + 1

        spectra_good = spectra_good[0:index_row, :]

        import statistics as st
        spectra_avg = ma.mean(spectra_good, axis=0)
        spectra_std = ma.std(spectra_good, axis=0)
        indices_max = ma.argmax(spectra_good, axis=0)
        imax = st.mode(indices_max)
        spectra_max_real = spectra_good[imax, :]
        spectra_max = ma.max(spectra_good, axis=0)

        indices_min = ma.argmin(spectra_good, axis=0)
        imin = st.mode(indices_min)
        spectra_mim_real = spectra_good[imin, :]
        spectra_min = ma.min(spectra_good, axis=0)

        spectra_stats = {
            'avg': spectra_avg,
            'std': spectra_std,
            'spectra_min_real': spectra_mim_real,
            'spectra_max_real': spectra_max_real,
            'spectra_min': spectra_min,
            'spectra_max': spectra_max,
        }

        return spectra_good, spectra_stats

    def get_flag_insitu_valid_spectra(self, scale_factor, flag_name, flag_value):
        import numpy.ma as ma
        nspectra = self.n_mu_total * self.n_insitu_day
        noriginal_bands = len(self.dimensions['insitu_original_bands'])
        spectra_good = ma.zeros((nspectra, noriginal_bands))
        var_insitu = self.nc.variables['insitu_Rrs']
        var_insitu_valid = self.nc.variables['insitu_valid']
        var_flag = self.nc.variables[flag_name]
        index_row = 0
        for index_mu in range(self.n_mu_total):
            for idx in range(self.n_insitu_day):
                if var_insitu_valid[index_mu, idx] == 1 and var_flag[index_mu] == flag_value:
                    spectra_here = ma.array(var_insitu[index_mu, :, idx]).transpose()
                    if scale_factor is not None:
                        spectra_here = spectra_here * scale_factor
                    spectra_good[index_row, :] = spectra_here[:]
                    index_row = index_row + 1

        spectra_good = spectra_good[0:index_row, :]

        import statistics as st
        spectra_avg = ma.mean(spectra_good, axis=0)
        spectra_std = ma.std(spectra_good, axis=0)
        indices_max = ma.argmax(spectra_good, axis=0)
        imax = st.mode(indices_max)
        spectra_max_real = spectra_good[imax, :]
        spectra_max = ma.max(spectra_good, axis=0)

        indices_min = ma.argmin(spectra_good, axis=0)
        imin = st.mode(indices_min)
        spectra_mim_real = spectra_good[imin, :]
        spectra_min = ma.min(spectra_good, axis=0)

        spectra_stats = {
            'avg': spectra_avg,
            'std': spectra_std,
            'spectra_min_real': spectra_mim_real,
            'spectra_max_real': spectra_max_real,
            'spectra_min': spectra_min,
            'spectra_max': spectra_max,
        }

        return spectra_good, spectra_stats

    def get_spectra_stats(self, spectra_good):
        import statistics as st
        import numpy.ma as ma
        # for s in spectra_good:
        #     print(s)
        spectra_avg = ma.mean(spectra_good, axis=0)
        # print(spectra_avg)
        spectra_std = ma.std(spectra_good, axis=0)
        indices_max = ma.argmax(spectra_good, axis=0)
        imax = st.mode(indices_max)
        spectra_max_real = spectra_good[imax, :]
        spectra_max = ma.max(spectra_good, axis=0)

        indices_min = ma.argmin(spectra_good, axis=0)
        imin = st.mode(indices_min)
        spectra_mim_real = spectra_good[imin, :]
        spectra_min = ma.min(spectra_good, axis=0)

        spectra_median = ma.median(spectra_good, axis=0)

        num_nonmasked = spectra_median.count()
        nspectra = spectra_median.shape[0]
        spectra_p25 = np.percentile(spectra_good, 25, axis=0)
        spectra_p75 = np.percentile(spectra_good, 75, axis=0)
        if spectra_p25.count() < num_nonmasked:
            spectra_p25 = np.ma.zeros(spectra_median.shape)
            spectra_p75 = np.ma.zeros(spectra_median.shape)
            for idx in range(nspectra):
                array_here = spectra_good[:, idx]
                array_here = array_here[array_here.mask == False]
                if len(array_here) > 0:
                    spectra_p25[idx] = np.percentile(array_here, 25)
                    spectra_p75[idx] = np.percentile(array_here, 75)
                else:
                    spectra_p25[idx] = ma.masked
                    spectra_p75[idx] = ma.masked

        spectra_stats = {
            'avg': spectra_avg,
            'std': spectra_std,
            'spectra_min_real': spectra_mim_real,
            'spectra_max_real': spectra_max_real,
            'spectra_min': spectra_min,
            'spectra_max': spectra_max,
            'median': spectra_median,
            'p25': spectra_p25,
            'p75': spectra_p75
        }
        return spectra_stats

    def get_mu_spectra_insitu_and_sat_withwlvalues(self, index_mu, wlvalues, scale_factor, flag_name, flag_value,
                                                   flag_array):
        if self.nc.variables[self.var_mu_valid][index_mu] == 1:

            if flag_name is not None:
                if flag_array is not None:
                    if flag_array[index_mu] != flag_value:
                        return None, None
                else:
                    if self.nc.variables[flag_name][index_mu] != flag_value:
                        return None, None

            import netCDF4
            dfv = netCDF4.default_fillvals['f4']
            var_insitu = np.array(self.nc.variables['mu_ins_rrs'])
            var_satrrs = np.array(self.nc.variables['mu_sat_rrs'])
            var_wl = np.array(self.nc.variables['mu_wavelength'])
            var_satid = np.array(self.nc.variables['mu_satellite_id'])
            insitu_spectra = var_insitu[var_satid == index_mu]
            sat_spectra = var_satrrs[var_satid == index_mu]
            wl = var_wl[var_satid == index_mu]

            n_bands = len(wlvalues)
            import numpy.ma as ma
            insitu_spectra_end = ma.masked_equal(np.zeros(n_bands), 0)
            # print(insitu_spectra_end)
            sat_spectra_end = ma.masked_equal(np.zeros(n_bands), 0)
            # DOTAL = False
            for idx in range(n_bands):
                wlhere = wlvalues[idx]
                if wlhere in wl:
                    insitu_spectra_end[idx] = insitu_spectra[wl == wlhere]
                    sat_spectra_end[idx] = sat_spectra[wl == wlhere]
                    if insitu_spectra_end[idx] == dfv:
                        insitu_spectra_end[idx] = np.ma.masked
                    if sat_spectra_end[idx] == dfv:
                        sat_spectra_end[idx] = np.ma.masked
            if scale_factor is not None:
                insitu_spectra_end = insitu_spectra_end * scale_factor
                sat_spectra_end = sat_spectra_end * scale_factor

            return insitu_spectra_end, sat_spectra_end

        return None, None

    def get_mu_spectra_insitu_and_sat(self, index_mu, scale_factor):

        if self.nc.variables[self.var_mu_valid][index_mu] == 1:
            var_insitu = np.array(self.nc.variables['mu_ins_rrs'])
            var_satrrs = np.array(self.nc.variables['mu_sat_rrs'])
            var_wl = np.array(self.nc.variables['mu_wavelength'])
            var_satid = np.array(self.nc.variables['mu_satellite_id'])

            insitu_spectra = var_insitu[var_satid == index_mu]
            sat_spectra = var_satrrs[var_satid == index_mu]
            if scale_factor is not None:
                insitu_spectra = insitu_spectra * scale_factor
                sat_spectra = sat_spectra * scale_factor
            wl = var_wl[var_satid == index_mu]
            valid = np.ones(wl.shape, dtype=bool)
            wl = wl[valid]
            insitu_spectra = insitu_spectra[valid]
            sat_spectra = sat_spectra[valid]
            return wl, insitu_spectra, sat_spectra
        else:
            return None, None, None

    def get_all_spectra_insitu_sat_with_wlvalues_group(self, scale_factor, wlvalues, flag_name, flag_value, flag_array,
                                                       group_value, group_array):
        import numpy.ma as ma
        nspectra = self.n_mu_total
        n_bands = len(wlvalues)
        insitu_spectra_good = ma.zeros((nspectra, n_bands))
        sat_spectra_good = ma.zeros((nspectra, n_bands))
        index_here = 0
        for index_mu in range(self.n_mu_total):
            if group_array[index_mu] != group_value:
                continue
            # wl, insitu_spectrum, sat_spectrum = self.get_mu_spectra_insitu_and_sat(index_mu, scale_factor)
            insitu_spectrum, sat_spectrum = self.get_mu_spectra_insitu_and_sat_withwlvalues(index_mu, wlvalues,
                                                                                            scale_factor, flag_name,
                                                                                            flag_value, flag_array)
            if insitu_spectrum is not None and sat_spectrum is not None:
                insitu_spectra_good[index_here, :] = insitu_spectrum[:]
                sat_spectra_good[index_here, :] = sat_spectrum[:]
                index_here = index_here + 1
        sat_spectra_good = sat_spectra_good[0:index_here, :]
        insitu_spectra_good = insitu_spectra_good[0:index_here, :]
        sat_stats = self.get_spectra_stats(sat_spectra_good)
        insitu_stats = self.get_spectra_stats(insitu_spectra_good)

        return sat_stats, insitu_stats

    def get_all_spectra_insitu_sat_with_wlvalues(self, scale_factor, wlvalues, flag_name, flag_value, flag_array):
        import numpy.ma as ma
        nspectra = self.n_mu_total
        n_bands = len(wlvalues)
        insitu_spectra_good = ma.zeros((nspectra, n_bands))
        sat_spectra_good = ma.zeros((nspectra, n_bands))
        index_here = 0
        for index_mu in range(self.n_mu_total):
            # wl, insitu_spectrum, sat_spectrum = self.get_mu_spectra_insitu_and_sat(index_mu, scale_factor)
            insitu_spectrum, sat_spectrum = self.get_mu_spectra_insitu_and_sat_withwlvalues(index_mu, wlvalues,
                                                                                            scale_factor, flag_name,
                                                                                            flag_value, flag_array)
            if insitu_spectrum is not None and sat_spectrum is not None:
                insitu_spectra_good[index_here, :] = insitu_spectrum[:]
                sat_spectra_good[index_here, :] = sat_spectrum[:]
                index_here = index_here + 1
        sat_spectra_good = sat_spectra_good[0:index_here, :]
        insitu_spectra_good = insitu_spectra_good[0:index_here, :]
        sat_stats = self.get_spectra_stats(sat_spectra_good)
        insitu_stats = self.get_spectra_stats(insitu_spectra_good)

        return sat_stats, insitu_stats

    def get_all_spectra_insitu_sat(self, scale_factor):
        import numpy.ma as ma
        nspectra = self.n_mu_total
        # n_bands = len(np.unique(np.array(self.nc.variables['mu_wavelength'])))
        n_bands = len(np.array(self.nc.variables['satellite_bands']))
        # n_bands = 8

        insitu_spectra_good = ma.zeros((nspectra, n_bands))
        sat_spectra_good = ma.zeros((nspectra, n_bands))
        index_here = 0
        wavelength = None
        for index_mu in range(self.n_mu_total):
            wl, insitu_spectrum, sat_spectrum = self.get_mu_spectra_insitu_and_sat(index_mu, scale_factor)
            if insitu_spectrum is not None and sat_spectrum is not None:
                wavelength = wl
                insitu_spectra_good[index_here, :] = insitu_spectrum[:]
                sat_spectra_good[index_here, :] = sat_spectrum[:]
                index_here = index_here + 1
        sat_spectra_good = sat_spectra_good[0:index_here, :]
        insitu_spectra_good = insitu_spectra_good[0:index_here, :]
        sat_stats = self.get_spectra_stats(sat_spectra_good)
        insitu_stats = self.get_spectra_stats(insitu_spectra_good)

        return wavelength, sat_stats, insitu_stats

    def get_flag_value(self, flag_var, flag):
        flag_variable = self.variables[flag_var]
        if isinstance(flag_variable.flag_meanings, list):
            flag_meanings = ' '.join(flag_variable.flag_meanings)
        else:
            flag_meanings = flag_variable.flag_meanings
        flag_list = flag_meanings.upper().split(' ')
        flag_values = flag_variable.flag_values
        flag_value = None
        try:
            index = flag_list.index(flag.upper())
            flag_value = flag_values[index]
        except:
            pass
        return flag_value

    def get_flag_list_names(self, flag_var):
        flag_names = []
        if flag_var not in self.variables:
            return flag_names
        flag_values = np.unique(self.nc.variables[flag_var][:])
        all_flag_values = self.nc.variables[flag_var].flag_values.tolist()
        all_flag_names = self.nc.variables[flag_var].flag_meanings.split(" ")
        if isinstance(all_flag_values, int):
            all_flag_values = [all_flag_values]
        for flag_value in flag_values:
            idx = all_flag_values.index(flag_value)
            if idx >= 0:
                flag_names.append(all_flag_names[idx])
        return flag_names

    def get_flag_name(self, flag_var, index_mu):
        flag_value = self.nc.variables[flag_var][index_mu]
        all_flag_values = self.nc.variables[flag_var].flag_values.tolist()
        if isinstance(all_flag_values, int):
            all_flag_values = [all_flag_values]
        all_flag_names = self.nc.variables[flag_var].flag_meanings.split(" ")
        # print('==================')
        # print(flag_var)
        # print(type(all_flag_values),all_flag_values)
        # print(all_flag_names)
        # print(flag_value)
        # print('=====================')
        idx = all_flag_values.index(flag_value)
        if idx >= 0:
            return all_flag_names[idx]
        else:
            return None

    def get_flag_names_mu(self, index_mu):
        sat = self.get_flag_name('flag_satellite', index_mu)
        sensor = self.get_flag_name('flag_sensor', index_mu)
        satsensor = f'{sat.upper()}{sensor.upper()}'
        site = self.get_flag_name('flag_site', index_mu).upper()
        ac = self.get_flag_name('flag_ac', index_mu).upper()
        # from datetime import datetime as dt
        date_here = self.sat_times[index_mu].strftime('%Y%m%d')
        return date_here, satsensor, site, ac

    def get_times_mu(self, index_mu):
        ac = self.get_flag_name('flag_ac', index_mu).upper()
        sat_time = self.sat_times[index_mu].strftime('%Y%m%d')
        ins_time_array = np.array(self.variables['insitu_time'])
        inst_time_here = datetime.utcfromtimestamp(float(ins_time_array[index_mu][0]))
        ins_time = inst_time_here.strftime('%Y%m%d%H%M%S')
        return sat_time, ins_time, ac

    def obtain_mu_valid_with_ac(self, ac_ref):
        # step 1: obtain valid keys with ac
        key_valids = []
        mu_valid = np.array(self.nc.variables['mu_valid'][:])
        for index_mu in range(self.n_mu_total):
            datehere, satsensor, site, ac = self.get_flag_names_mu(index_mu)
            if mu_valid[index_mu] and ac == ac_ref:
                key_here = f'{datehere}{satsensor}{site}'
                key_valids.append(key_here)
        # step2: obtain mu with valid vales for specific ac
        mu_valid_ac_ref = np.zeros(self.n_mu_total, dtype=np.int)
        for index_mu in range(self.n_mu_total):
            datehere, satsensor, site, ac = self.get_flag_names_mu(index_mu)
            if ac == ac_ref:
                if mu_valid[index_mu]:
                    mu_valid_ac_ref[index_mu] = 1
            else:
                key_here = f'{datehere}{satsensor}{site}'
                if key_here in key_valids:
                    mu_valid_ac_ref[index_mu] = 1

        return mu_valid_ac_ref

    def obtain_info_common_mu(self):

        acs = self.get_flag_list_names('flag_ac')
        if len(acs) <= 1:
            print(f'[INFO] Common match-ups are set only if more than one AC processor is available')
            return None
        print(f'[INFO] Potential atmospheric correction processors: {acs}')
        info_common_mu = {}
        mu_valid = np.array(self.nc.variables['mu_valid'][:])
        print(f'[INFO] Checking MU info...')
        for index_mu in range(self.n_mu_total):
            datehere, satsensor, site, ac = self.get_flag_names_mu(index_mu)
            fill_acs = False
            if satsensor not in info_common_mu.keys():
                info_common_mu[satsensor] = {site: {datehere: {}}}
                fill_acs = True
            if satsensor in info_common_mu.keys() and site not in info_common_mu[satsensor]:
                info_common_mu[satsensor][site] = {datehere: {}}
                fill_acs = True
            if satsensor in info_common_mu.keys() and site in info_common_mu[satsensor] and datehere not in \
                    info_common_mu[satsensor][site]:
                info_common_mu[satsensor][site][datehere] = {}
                fill_acs = True
            if fill_acs:
                for ac_here in acs:
                    info_common_mu[satsensor][site][datehere][ac_here] = 0

            if mu_valid[index_mu]:
                info_common_mu[satsensor][site][datehere][ac] = 1

        print(f'[INFO] Obtaining valid common match-ups...')
        mu_valid_common = np.zeros(self.n_mu_total, dtype=np.int)
        n_mu_valid_common = 0
        for index_mu in range(self.n_mu_total):
            if mu_valid[index_mu] == 1:
                datehere, satsensor, site, ac = self.get_flag_names_mu(index_mu)
                sum = 0
                for ac_here in acs:
                    sum = sum + info_common_mu[satsensor][site][datehere][ac_here]

                if sum == len(acs):
                    mu_valid_common[index_mu] = 1
                    n_mu_valid_common = n_mu_valid_common + 1

        # print(f'[INFO]n_mu_valid_common')
        # print(len(acs))
        n_mu_valid_common = n_mu_valid_common / len(acs)
        print(f'[INFO] # valid commons match-ups: {n_mu_valid_common}')

        if n_mu_valid_common > 0:
            return mu_valid_common
        else:
            return None

    def obtain_info_common_mu_nosat(self):

        acs = self.get_flag_list_names('flag_ac')
        # acs = ['CCI','POLYMER']
        if len(acs) <= 1:
            print(f'[INFO] Common match-ups are set only if more than one AC processor is available')
            return None
        print(f'[INFO] Potential atmospheric correction processors: {acs}')
        info_common_mu = {}
        mu_valid = np.array(self.nc.variables['mu_valid'][:])
        print(f'[INFO] Checking MU info...')
        for index_mu in range(self.n_mu_total):
            datehere, satsensor, site, ac = self.get_flag_names_mu(index_mu)
            fill_acs = False
            if site not in info_common_mu:
                info_common_mu[site] = {datehere: {}}
                fill_acs = True
            if datehere not in info_common_mu[site]:
                info_common_mu[site][datehere] = {}
                fill_acs = True
            if fill_acs:
                for ac_here in acs:
                    info_common_mu[site][datehere][ac_here] = 0

            if mu_valid[index_mu]:
                info_common_mu[site][datehere][ac] = 1

        print(f'[INFO] Obtaining valid common match-ups...')
        mu_valid_common = np.zeros(self.n_mu_total, dtype=np.int)
        n_mu_valid_common = 0
        for index_mu in range(self.n_mu_total):
            if mu_valid[index_mu] == 1:
                datehere, satsensor, site, ac = self.get_flag_names_mu(index_mu)

                if ac not in acs:
                    mu_valid_common[index_mu] = 1
                    continue

                sum = 0
                for ac_here in acs:
                    sum = sum + info_common_mu[site][datehere][ac_here]

                if sum == len(acs):
                    mu_valid_common[index_mu] = 1
                    n_mu_valid_common = n_mu_valid_common + 1

        # print(f'[INFO]n_mu_valid_common')
        # print(len(acs))
        n_mu_valid_common = n_mu_valid_common / len(acs)
        print(f'[INFO] # valid commons match-ups: {n_mu_valid_common}')

        if n_mu_valid_common > 0:
            return mu_valid_common
        else:
            return None

    def obtain_info_common_mu_ins(self):
        acs = self.get_flag_list_names('flag_ac')
        if len(acs) <= 1:
            print(f'[INFO] Common match-ups are set only if more than one AC processor is available')
            return None
        print(f'[INFO] Potential atmospheric correction processors: {acs}')
        info_common_mu = {}
        mu_valid = np.array(self.nc.variables['mu_valid'][:])
        print(f'[INFO] Checking MU info...')
        for index_mu in range(self.n_mu_total):
            sat_time, ins_time, ac = self.get_times_mu(index_mu)
            key_time = f'{sat_time}_{ins_time}'
            if key_time not in info_common_mu.keys():
                info_common_mu[key_time] = {}
                for ac_here in acs:
                    info_common_mu[key_time][ac_here] = 0
            if mu_valid[index_mu]:
                info_common_mu[key_time][ac] = 1

        print(f'[INFO] Obtaining valid common match-ups...')
        mu_valid_common = np.zeros(self.n_mu_total, dtype=np.int)
        n_mu_valid_common = 0
        for index_mu in range(self.n_mu_total):
            if mu_valid[index_mu] == 1:
                sat_time, ins_time, ac = self.get_times_mu(index_mu)
                key_time = f'{sat_time}_{ins_time}'
                sum = 0
                for ac_here in acs:
                    sum = sum + info_common_mu[key_time][ac_here]
                if sum == len(acs):
                    mu_valid_common[index_mu] = 1
                    n_mu_valid_common = n_mu_valid_common + 1

        # print(f'[INFO]n_mu_valid_common')
        # print(len(acs))
        n_mu_valid_common = n_mu_valid_common / len(acs)
        print(f'[INFO] # valid commons match-ups: {n_mu_valid_common}')

        if n_mu_valid_common > 0:
            return mu_valid_common
        else:
            return None

    def analyze_sat_flags_advanced(self, flag_options, valid_array):
        flag_info = {}
        for flag_output in flag_options:
            flag_var_name = flag_options[flag_output]['flag_var_name']
            flag_ref = flag_options[flag_output]['flag_ref']
            flag_method = flag_options[flag_output]['flag_method']
            var_group_name = flag_options[flag_output]['var_group_name']
            var_group_val = flag_options[flag_output]['var_group_value']
            seriesid = flag_options[flag_output]['seriesid']
            seriesoutput = flag_options[flag_output]['seriesoutput']

            array, nmu, nmacro, nmacrow = self.analyze_sat_flags_advanced_impl(flag_var_name, flag_ref, flag_method,
                                                                               var_group_name, var_group_val,
                                                                               valid_array)

            parray = (array / self.n_mu_total) * 100
            central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()

            flag_info[flag_output] = {
                'array': array,
                'parray': parray,
                'ntotal': np.sum(array),
                'ntotalw': np.sum(array[r_s:r_e, c_s:c_e]),
                'nmu': nmu,
                'nmacro': nmacro,
                'nmacrow': nmacrow,
                'pmacro': float((nmacro / nmu) * 100),
                'pmacrow': float((nmacrow / nmu) * 100),
                'seriesid': seriesid,
                'seriesoutput': seriesoutput
            }

        return flag_info

    def analyze_sat_flags_advanced_impl(self, flag_var_name, flag_ref, flag_method, var_group_name, var_group_val,
                                        valid_array):
        flag_variable = self.nc.variables[flag_var_name]
        if isinstance(flag_variable.flag_meanings, list):
            flag_meanings = ' '.join(flag_variable.flag_meanings)
        else:
            flag_meanings = flag_variable.flag_meanings
            flag_meanings = flag_meanings.replace('  ', ' ')

        flag_masks = flag_variable.flag_masks
        # if np.min(flag_masks)>0 and np.max(flag_masks)<=np.iinfo(np.int64)
        # print(np.min(flag_masks),np.max(flag_masks))

        flag_masks = flag_masks.astype(np.uint64)

        # if flag_var_name=='satellite_bitmask':
        #     print('-------------------------->',flag_meanings)
        #     print('-------------------------->',type(flag_meanings))
        #     print('-------------------------->',flag_masks)

        flagging = flag.Class_Flags_OLCI(flag_masks, flag_meanings)
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        nrows = len(self.dimensions['rows'])
        ncols = len(self.dimensions['columns'])
        array = np.zeros((nrows, ncols))
        nmu = 0
        nmacro = 0
        nmacrow = 0

        if flag_method == 'eq':
            flag_ref_list = [flag_ref]

        if flag_method == 'list':
            flag_ref_list = []
            flag_ref_l = flag_ref.split(',')
            flag_all_list = [x.strip() for x in flag_meanings.split(' ')]
            for name in flag_ref_l:
                if name.strip() in flag_all_list:
                    flag_ref_list.append(name.strip())

        group_array = None
        if var_group_name is not None:
            group_array = np.array(self.variables[var_group_name][:])

        for index_mu in range(self.n_mu_total):
            if group_array is not None:
                if group_array[index_mu] != var_group_val:
                    continue
            if valid_array is not None:
                if valid_array[index_mu] == 0:
                    continue
            nmu = nmu + 1
            flag_array_mu = np.array(self.variables[flag_var_name][index_mu]).astype(np.uint64)

            mask = flagging.Mask(flag_array_mu, flag_ref_list)
            # if flag_var_name == 'satellite_c2rcc_flags' and flag_ref_list[0]=='Cloud_risk':
            #     print(index_mu,'--->',flag_array_mu[11,12],flag_array_mu[12,12],flag_array_mu[13,12])
            mask[np.where(mask != 0)] = 1
            if np.sum(mask) > 0:
                nmacro = nmacro + 1
            mask_window = mask[r_s:r_e, c_s:c_e]

            if np.sum(mask_window) > 0:
                nmacrow = nmacrow + 1
            # if flag_ref == 'sun_glint_risk' and np.sum(mask_window)==9:
            #     print(index_mu)

            array = array + mask

        return array, nmu, nmacro, nmacrow

    def analyze_sat_flags(self, flag_var_name, flag_list):
        flag_variable = self.nc.variables[flag_var_name]
        if isinstance(flag_variable.flag_meanings, list):
            flag_meanings = ' '.join(flag_variable.flag_meanings)
        else:
            flag_meanings = flag_variable.flag_meanings
        if flag_list is None:
            flag_list = flag_variable.flag_meanings.split(' ')
        flag_info = {}

        flag_masks = flag_variable.flag_masks
        if flag_masks.dtype == 'int64':
            flag_masks = flag_masks.astype(np.uint64)

        flagging = flag.Class_Flags_OLCI(flag_masks, flag_meanings)
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        for flag_name in flag_list:
            # flag_value = flag_values[idx]
            # flag_name = flag_list[idx]
            array = np.zeros((25, 25))
            nmacro = 0
            nmacrow = 0
            for index_mu in range(self.n_mu_total):
                flag_array_mu = np.array(self.variables[self.flag_band_name][index_mu])
                mask = flagging.Mask(flag_array_mu, [flag_name])
                mask[np.where(mask != 0)] = 1
                if np.sum(mask) > 0:
                    nmacro = nmacro + 1
                mask_window = mask[r_s:r_e, c_s:c_e]
                if np.sum(mask_window) > 0:
                    nmacrow = nmacrow + 1
                array = array + mask
            parray = (array / self.n_mu_total) * 100
            flag_info[flag_name] = {
                'array': array,
                'parray': parray,
                'ntotal': np.sum(array),
                'nmacro': nmacro,
                'nmacrow': nmacrow,
                'pmacro': float((nmacro / self.n_mu_total) * 100),
                'pmacrow': float((nmacrow / self.n_mu_total) * 100)
            }
        return flag_info

    def analyse_mu_temporal(self, onlyvalid, file_out):
        year_min = self.start_date.year
        year_max = self.end_date.year + 1
        year = list(range(year_min, year_max))
        year.reverse()
        month = list(range(1, 13))
        dfall_month = pd.DataFrame(index=year, columns=month, dtype=np.float)
        dfall_month[:] = 0
        mu_valid = np.array(self.nc.variables['mu_valid'])
        for idx in range(self.n_mu_total):
            date_here = self.sat_times[idx]
            year_here = date_here.year
            month_here = date_here.month
            if onlyvalid:
                if mu_valid[idx] == 1:
                    dfall_month.at[year_here, month_here] = dfall_month.at[year_here, month_here] + 1
            else:
                dfall_month.at[year_here, month_here] = dfall_month.at[year_here, month_here] + 1
        from matplotlib import cm
        import seaborn as sns
        h = plt.Figure()
        cmap = cm.get_cmap('RdYlBu_r')
        dfall_month_withnan = dfall_month
        dfall_month_withnan[dfall_month == 0] = np.nan
        sns.heatmap(dfall_month_withnan, annot=False, cmap=cmap, linewidths=1, linecolor='black')
        plt.xlabel('Month')
        plt.ylabel('Year')
        # plt.title(title)
        if file_out is not None:
            plt.savefig(file_out, dpi=300)
        plt.close(h)

    def analyse_mu_temporal_flag(self, onlyvalid, varvalid, name_flag_var, flag_list):
        year_min = self.start_date.year
        year_max = self.end_date.year + 1

        correct = 1
        if varvalid == 'mu_valid_common':  # number of common mu is divided by the number of ac
            if 'flag_ac' in self.nc.variables:
                flag_var_ac = self.nc.variables['flag_ac']
                correct = len(flag_var_ac.flag_meanings.split(' '))

        monthl = []
        for year in range(year_min, year_max):
            for month in range(1, 13):
                # if year == 2023 and month >= 4:
                #     continue
                date_here = datetime(year, month, 1)
                # print(date_here)
                monthl.append(date_here.strftime('%Y-%m'))

        print(monthl)

        flag_var = self.nc.variables[name_flag_var]
        flag_var_array = np.array(flag_var)
        if flag_list is None:
            flag_list = flag_var.flag_meanings.split(' ')
            flag_values = flag_var.flag_values
        else:
            all_flag_list = flag_var.flag_meanings.split(' ')
            all_flag_values = flag_var.flag_values
            flag_values = np.zeros(all_flag_values.shape)
            for idx in range(len(flag_list)):
                flag_here = flag_list[idx]
                index = all_flag_list.index(flag_here)
                flag_values[idx] = all_flag_values[index]

        dfall_month = pd.DataFrame(index=flag_list, columns=monthl, dtype=np.float)
        dfall_month[:] = 0
        total_mu = False
        if varvalid is not None and not onlyvalid:
            total_mu = True
        if varvalid is None:
            varvalid = 'mu_valid'
        mu_valid = np.array(self.nc.variables[varvalid])
        for idx in range(self.n_mu_total):
            date_here = self.sat_times[idx]
            date_here_str = date_here.strftime('%Y-%m')
            flag_value = flag_var_array[idx]
            nw = np.where(flag_values == flag_value)
            if len(nw[0]) == 0:
                continue
            flag_index = nw[0][0]
            flag_here = flag_list[flag_index]
            if onlyvalid:
                if mu_valid[idx] == 1:
                    dfall_month.at[flag_here, date_here_str] = dfall_month.at[flag_here, date_here_str] + (1 / correct)
            else:
                if total_mu:
                    if mu_valid[idx] >= 0:
                        dfall_month.at[flag_here, date_here_str] = dfall_month.at[flag_here, date_here_str] + (
                                1 / correct)
                else:
                    dfall_month.at[flag_here, date_here_str] = dfall_month.at[flag_here, date_here_str] + (1 / correct)
        return dfall_month
        # from matplotlib import cm
        # import seaborn as sns
        # h = plt.Figure()
        # cmap = cm.get_cmap('RdYlBu_r')
        # dfall_month_withnan = dfall_month
        # dfall_month_withnan[dfall_month == 0] = np.nan
        # sns.heatmap(dfall_month_withnan, annot=False, cmap=cmap, linewidths=1, linecolor='black')
        # plt.xlabel('Year-Month')
        # plt.ylabel('Site')
        # plt.gcf().tight_layout()
        # if file_out is not None:
        #     plt.savefig(file_out, dpi=300)
        #     file_out_csv = file_out.replace('.jpg','.csv')
        #     dfall_month_withnan.to_csv(file_out_csv,sep=';')
        # plt.close(h)

    # varvalid: mu_valid or mu_valid_common
    def analyse_mu_flag(self, name_flag_var, flag_list, varvalid):
        flag_var = self.nc.variables[name_flag_var]
        flag_var_array = np.array(flag_var)
        correct = 1
        ntot = 1
        if varvalid == 'mu_valid_common':  # number of common mu is divided by the number of ac
            if 'flag_ac' in self.nc.variables:
                flag_var_ac = self.nc.variables['flag_ac']
                flag_var_ac_array = np.array(flag_var_ac)
                flag_ac_values = flag_var_ac.flag_values
                correct = len(flag_var_ac.flag_meanings.split(' '))
                ntot = ntot + correct

        print(f'CORRECT VALUE {correct}')

        if flag_list is None:
            # flag_list = flag_var.flag_meanings.split(' ')
            flag_values = flag_var.flag_values
        else:
            all_flag_list = flag_var.flag_meanings.split(' ')
            all_flag_values = flag_var.flag_values
            flag_values = np.zeros(all_flag_values.shape)
            for idx in range(len(flag_list)):
                flag_here = flag_list[idx]
                index = all_flag_list.index(flag_here)
                flag_values[idx] = all_flag_values[index]

        mu_valid = np.array(self.nc.variables[varvalid])

        nflag = len(flag_values)
        n_values = np.zeros((nflag, 1))
        n_valid = np.zeros((nflag, ntot))

        for idx in range(self.n_mu_total):
            flag_value = flag_var_array[idx]

            nw = np.where(flag_values == flag_value)
            if len(nw[0]) == 0:
                continue
            flag_index = nw[0][0]
            # flag_index = np.where(flag_values == flag_value)[0][0]
            n_values[flag_index, 0] = n_values[flag_index, 0] + (1 / correct)
            if mu_valid[idx] == 1:
                n_valid[flag_index, 0] = n_valid[flag_index, 0] + (1 / correct)

        if ntot > 1:  ##data for each ac,
            mu_valid = np.array(self.nc.variables['mu_valid'])
            for idx in range(self.n_mu_total):
                flag_value = flag_var_array[idx]
                flag_index = np.where(flag_values == flag_value)[0][0]
                flag_ac_here = flag_var_ac_array[idx]
                flag_ac_index = np.where(flag_ac_values == flag_ac_here)[0][0] + 1
                if mu_valid[idx] == 1:
                    n_valid[flag_index, flag_ac_index] = n_valid[flag_index, flag_ac_index] + 1

        # print(n_values, n_valid)
        porc_values = (n_valid / n_values) * 100

        n_values = np.floor(n_values)

        return n_values, n_valid, porc_values

        # print(porc_values)
        # plt.figure()
        # plt.barh(flag_list, porc_values)
        # plt.grid(b=True, which='major', color='gray', linestyle='--')
        # plt.xlabel('% Valid Match-ups', fontsize=12)
        # plt.gcf().tight_layout()
        # if file_out is not None:
        #     plt.savefig(file_out, dpi=300)
        # plt.close()

    def get_full_array(self, namevariable):
        if self.VALID and namevariable in self.variables:
            return np.array(self.variables[namevariable][:])
        else:
            return None

    def get_dims_variable(self, namevariable):
        if self.VALID and namevariable in self.variables:
            return self.variables[namevariable].dimensions
        else:
            return None

    def get_fill_value(self, namevariable):
        fillValue = None
        if self.VALID and namevariable in self.variables:
            try:
                fillValue = self.variables[namevariable]._FillValue
            except:
                pass
        return fillValue

    def close(self):
        if self.VALID:
            self.nc.close()
