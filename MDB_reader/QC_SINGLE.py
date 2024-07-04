import os

import pytz
from netCDF4 import Dataset
import numpy as np
import COMMON.Class_Flags_OLCI as flag
import math


class QC_SINGLE:

    def __init__(self, path_mdb, satellite_variable, insitu_variable):
        self.path_mdb = path_mdb
        self.dataset = None
        if os.path.exists(path_mdb):
            try:
                self.dataset = Dataset(path_mdb)
            except:
                self.dataset = None
        self.satellite_variable = satellite_variable
        self.insitu_variable = insitu_variable

        self.window_size = 3
        self.min_valid_pixels = 9
        self.use_Bailey_Werdell = False
        self.apply_outliers = True
        self.outliers_info = {
            'central_stat': 'avg',
            'dispersion_stat': 'std',
            'factor': 1.5
        }
        self.time_diff_max = 120.0 * 60.0
        self.fix_time_sat = None

        self.nmu = -1
        self.NTP = self.window_size * self.window_size  # total number of pixels
        self.NTPW = self.NTP  # total number of water pixels (excluding land/inland waters), could vary with MU
        self.NVP = 0  # number of valid pixels (excluding flag pixels), varies with MU
        self.flag_mask = None  # mask
        self.ac_processor = None
        self.n_masked_invalid = 0

        self.info_flag = {}  ##masks based o flags
        self.th_masks = []  ##masks based on therholds
        self.check_statistics = []  ##macropixel filters based on statistics
        self.insitu_flags = {}  ##insitu flags
        self.insitu_th = {}  ##insitu threhlods

        stat_list = ['n_values', 'avg', 'std', 'median', 'min', 'max', 'CV']
        self.statistics = {
            'without_outliers': {x: 0 for x in stat_list},
            'with_outliers': {x: 0 for x in stat_list}
        }

        ##temporal sampling method: closest, only_single
        self.temporal_sampling_method = 'closest'

    def check_sat_insitu_variables(self):
        if self.dataset is None:
            print(
                f'[ERROR] Error setting Quality Control protocol: dataset is not defined, {self.path_mdb} is not a valid NetCDF file')
            return False
        if not self.satellite_variable in self.dataset.variables:
            print(
                f'[ERROR] Error setting Quality Control protocol: {self.satellite_variable} is not available in path: {self.path_mdb}')
            return False
        if not self.insitu_variable in self.dataset.variables:
            print(
                f'[ERROR] Error setting Quality Control protocol: {self.insitu_variable} is not available in path: {self.path_mdb}')
            return False
        var_sat = self.dataset.variables[self.satellite_variable]
        if len(var_sat.shape) == 3:
            self.nmu = var_sat.shape[0]
        else:
            print(f'[ERROR] Satellite variable {self.satellite_variable} should have 3 dimensions')
            return False
        var_ins = self.dataset.variables[self.insitu_variable]
        if len(var_ins.shape) == 2:
            if var_ins.shape[0] == self.nmu:
                print(f'[INFO] Number of potential match-ups set to: {self.nmu}')
            else:
                print(
                    f'[ERROR] Disagreement between the number of match-ups between the satellite and insitu variable ({self.nmu} vs. {var_ins.shape[0]})')
                return False
        else:
            print(f'[ERROR] In situ variable {self.insitu_variable} should have 2 dimensions')
            return False
        return True

    def set_window_size(self, wsize):
        self.window_size = wsize
        self.NTP = self.window_size * self.window_size
        self.NTPW = self.NTP
        if self.min_valid_pixels > self.NTP:
            self.min_valid_pixels = self.NTP

    def add_threshold_mask(self, band_name, value_th, type_th):
        th_mask = {
            'index_sat': -1,
            'band_name': band_name,
            'value_th': value_th,
            'type_th': type_th,
            'n_masked': 0
        }
        self.th_masks.append(th_mask)

    def add_bands_statistics(self, name_band, type_stat, value_th, type_th):
        if self.ncdataset is None:
            print('[WARNING] Statistics for bands could not be added as dataset was not defined')
            return
        if name_band not in self.ncdataset.variables:
            return
        var_band = self.ncdataset.variables[name_band]
        check_val = {
            'variable': var_band,
            'name_band': name_band,
            'type_stat': type_stat,
            'value_th': value_th,
            'type_th': type_th
        }
        self.check_statistics.append(check_val)

    def add_insitu_flag_expression(self, flag_band, flag_list, remove_spectra):
        if self.ncdataset is None:
            return

        if flag_band not in self.ncdataset.variables:
            return

        flag_variable = self.ncdataset.variables[flag_band]

        if flag_list == 'ALL':  # APPLY ALL THE FLAGS
            flag_list = flag_variable.flag_meanings.split()
        else:
            flag_list = [x.strip() for x in flag_list.split(',')]

        oflags = None
        try:
            from COMMON.Class_Flags_OLCI import Class_Flags_OLCI
            flag_values = [np.uint32(x.strip()) for x in flag_variable.flag_mask.split(',')]
            oflags = Class_Flags_OLCI(flag_values, flag_variable.flag_meanings)
        except:
            print(f'[WARNING] Flag class could not be defined for variable: {flag_band}')
            pass

        self.insitu_flags[flag_band] = {
            'variable': flag_variable,
            'flag_list': flag_list,
            'remove_spectra': remove_spectra,
            'oflags': oflags
        }

    def add_insitu_band_thersholds(self, band_name, type_th, value_min, value_max, isangle):
        if self.ncdataset is None:
            return

        if band_name not in self.ncdataset.variables:
            return

        band_variable = self.ncdataset.variables[band_name]

        self.insitu_th[band_variable] = {
            'variable': band_variable,
            'th_type': type_th,
            'value_min': value_min,
            'value_max': value_max,
            'isangle': isangle
        }

    def prepare_new_match_up(self):
        self.NTPW = self.NTP  # total number of water pixels (excluding land/inland waters), could vary with MU
        self.NVP = 0  # number of valid pixels (excluding flag pixels), varies with MU
        self.flag_mask = None  # mask based on flagging
        stat_list = ['n_values', 'avg', 'std', 'median', 'min', 'max', 'CV']
        self.statistics = {
            'without_outliers': {x: 0 for x in stat_list},
            'with_outliers': {x: 0 for x in stat_list}
        }

    ##main method for masking
    def compute_masks_and_check_roi(self, index_mu):
        land = self.compute_flag_masks(index_mu)
        self.compute_invalid_mask(index_mu)
        self.compute_th_masks(index_mu)
        self.NVP = self.NTP - np.sum(self.flag_mask)
        self.NTPW = self.NTP - np.sum(land, axis=(0, 1))
        min_valid_pixels = self.min_valid_pixels
        if self.use_Bailey_Werdell:
            min_valid_pixels = math.floor(0.50 * self.NTPW) + 1
        cond_min_pixels = False
        if self.NVP >= min_valid_pixels:
            cond_min_pixels = True

        return cond_min_pixels

    ##masks based on flags
    def compute_flag_masks(self, index_mu):
        flag_mask = np.zeros((self.window_size, self.window_size), dtype=np.uint64)
        land = np.zeros((self.window_size, self.window_size), dtype=np.uint64)
        for flag_band in self.info_flag.keys():
            flag_mask_here, land_here = self.compute_flag_mask_impl(index_mu, flag_band)
            if flag_mask_here is not None:
                flag_mask = np.add(flag_mask, flag_mask_here)
            if land_here is not None:
                land = np.add(land, land_here)
            self.info_flag[flag_band]['nflagged'] = np.sum(flag_mask)

        if self.flag_mask is None:
            self.flag_mask = flag_mask

        self.flag_mask[flag_mask > 0] = 1

        land[land > 0] = 1

        return land

    def compute_flag_mask_impl(self, index_mu, flag_band):
        land = None
        flag_mask = None
        if index_mu < 0 or index_mu >= self.nmu:
            return land, flag_mask
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        if flag_band not in self.info_flag.keys():
            return land, flag_mask
        satellite_flag = self.info_flag[flag_band]['variable']
        if satellite_flag is None:
            return land, flag_mask
        if not 'flag_meanings' in satellite_flag:
            return land, flag_mask
        # flag_meanings_ string separated by spaces or list, it's compulsolry
        flag_meanings = satellite_flag.flag_meanings
        if isinstance(flag_meanings, list):
            flag_meanings = ' '.join(flag_meanings)
        satellite_flag_band = satellite_flag[index_mu, r_s:r_e, c_s:c_e]
        # float32 is not allowed
        if str(satellite_flag.dtype) == 'float32':
            satellite_flag_band = satellite_flag_band.astype('uint64')

        ##flag_list_tobe_applied, comma-separared str or list
        flag_list_tobe_applied = self.info_flag[flag_band]['flag_list']
        if isinstance(flag_list_tobe_applied, str):
            flag_list_tobe_applied = [x.strip() for x in flag_list_tobe_applied.split(',')]

        if self.info_flag[flag_band]['ac_processor'] == 'POLYMER':
            flagging = flag.Class_Flags_Polymer(satellite_flag.flag_masks, flag_meanings)
            flag_mask = flagging.MaskGeneral(satellite_flag_band)
            flag_mask[np.where(flag_mask != 0)] = 1
        elif self.info_flag[flag_band]['ac_processor'] == 'IDEPIX':
            flagging = flag.Class_Flags_Idepix(satellite_flag.flag_masks, flag_meanings)
            flag_mask = flagging.Mask(satellite_flag_band, flag_list_tobe_applied)
            flag_mask[np.where(flag_mask != 0)] = 1
        else:
            ##we must be sure that flag_mask must be uint64
            flag_masks = satellite_flag.flag_masks.astype('uint64')
            flagging = flag.Class_Flags_OLCI(flag_masks, flag_meanings)
            flag_mask = flagging.Mask(satellite_flag_band, flag_list_tobe_applied)
            flag_mask[np.where(flag_mask != 0)] = 1

        flag_land = self.info_flag[flag_band]['flag_land']
        if flag_land is not None and flag_land.strip().lower() == 'none':
            flag_land = None
        flag_inlandwater = self.info_flag[flag_band]['flag_inlandwater']
        if flag_inlandwater is not None and flag_inlandwater.strip().lower() == 'none':
            flag_inlandwater = None

        if flag_land is not None:
            land = flagging.Mask(satellite_flag_band, ([flag_land]))
            land[np.where(land != 0)] = 1
            if flag_inlandwater is not None:
                inland_w = flagging.Mask(satellite_flag_band, ([flag_inlandwater]))
                land[np.where(inland_w != 0)] = 0

        return flag_mask, land

    def compute_invalid_mask(self, index_mu):
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        datahere = np.ma.squeeze(self.dataset[self.satellite_variable][index_mu, r_s:r_e, c_s:c_e])
        mask_invalid = np.zeros(datahere.shape, dtype=np.uint64)
        mask_invalid[datahere.mask] = 1
        self.n_masked_invalid = np.sum(mask_invalid)
        if self.flag_mask is None:
            self.flag_mask = mask_invalid
        self.flag_mask[mask_invalid > 0] = 1

    def compute_th_masks(self, index_mu):
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        mask_thershold = np.zeros((self.window_size, self.window_size), dtype=np.uint64)

        for idx in range(len(self.th_masks)):
            th_mask = self.th_masks[idx]
            if th_mask['index_sat'] >= 0:
                band_here = self.satellite_rrs[index_mu, th_mask['index_sat'], r_s:r_e, c_s:c_e]
            else:
                var_here = self.ncdataset.variables[th_mask['band_name']]
                band_here = var_here[index_mu, r_s:r_e, c_s:c_e]

            mask_thershold_here = np.zeros(band_here.shape, dtype=np.uint64)
            if th_mask['type_th'] == 'greater':
                mask_thershold_here[band_here > th_mask['value_th']] = 1
            elif th_mask['type_th'] == 'lower':
                mask_thershold_here[band_here < th_mask['value_th']] = 1
            n_masked = np.sum(mask_thershold_here)
            th_mask['n_masked'] = n_masked
            self.th_masks[idx] = th_mask
            mask_thershold = mask_thershold + mask_thershold_here

        if self.flag_mask is None:
            self.flag_mask = mask_thershold

        self.flag_mask[mask_thershold > 0] = 1

    def compute_statistics(self, index_mu):
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        data_here = self.dataset.variables[self.satellite_variable][index_mu, r_s:r_e, c_s:c_e]
        data_valid = data_here[self.flag_mask == 0]
        stats = self.compute_statistics_impl(self.statistics['without_outliers'], data_valid)
        self.statistics['without_outliers'] = stats
        if self.apply_outliers:
            cvalue = self.statistics['without_outliers'][self.outliers_info['central_stat']]
            dvalue = self.statistics['without_outliers'][self.outliers_info['dispersion_stat']]
            min_th = cvalue - (dvalue * self.outliers_info['factor'])
            max_th = cvalue + (dvalue * self.outliers_info['factor'])
            mask_outliers = np.zeros(data_valid.shape)
            mask_outliers[data_valid > max_th] = 1
            mask_outliers[data_valid < min_th] = 1
            n_outliers = np.sum(mask_outliers)
            if n_outliers > 0:
                data_valid = data_valid[mask_outliers == 0]
            stats = self.compute_statistics_impl(self.statistics['with_outliers'], data_valid)
            self.statistics['with_outliers'] = stats

        for check_stat_here in self.check_statistics:
            name_band = check_stat_here['name_band']
            var_here = check_stat_here['variable']
            var_here_array = var_here[index_mu, r_s:r_e, c_s:c_e]
            var_here_valid = var_here_array[~var_here_array.mask]
            if not name_band in self.statistics:
                self.statistics[name_band] = {
                    'n_values': 0,
                    'avg': 0,
                    'std': 0,
                    'median': 0,
                    'min': 0,
                    'max': 0,
                    'CV': 0
                }
            stats = self.compute_statistics_impl(self.statistics[name_band], var_here_valid)
            self.statistics[name_band] = stats

        return True

    def compute_statistics_impl(self, stats, array):
        if not np.any(array):
            stats['n_values'] = 0
            stats['avg'] = 0
            stats['std'] = 0
            stats['median'] = 0
            stats['min'] = 0
            stats['max'] = 0
            stats['CV'] = 0
        else:
            stats['n_values'] = len(array)
            stats['avg'] = np.mean(array)
            stats['std'] = np.std(array)
            stats['median'] = np.median(array)
            stats['min'] = np.min(array)
            stats['max'] = np.max(array)
            CV = (stats['std'] / abs(stats['avg'])) * 100
            stats['CV'] = CV
        return stats

    def do_check_statistics(self):
        CHECK = True
        for check_stat in self.check_statistics:
            name_band = check_stat['name_band']
            if name_band in self.statistics:
                val_here = self.statistics[name_band][check_stat['type_stat']]
                if check_stat['type_th'] == 'greater' and val_here > check_stat['value_th']:
                    CHECK = False
                if check_stat['type_th'] == 'lower' and val_here < check_stat['value_th']:
                    CHECK = False

        return CHECK

    def get_match_up_value(self, index_mu):
        self.prepare_new_match_up()
        cond_min_pixels = self.compute_masks_and_check_roi(index_mu)
        cond_stats = False
        valid_mu = False
        value = -999.0
        if cond_min_pixels:
            outliers_str = 'without_outliers'
            if self.apply_outliers:
                outliers_str = 'with_outliers'
            self.compute_statistics(index_mu)
            cond_stats = self.do_check_statistics()
            valid_mu = True if cond_stats else False
            if valid_mu:
                value = self.statistics[outliers_str][self.stat_value]

        return cond_min_pixels, cond_stats, valid_mu, value

    def get_ins_value(self,index_mu):
        data_ins_mu = self.dataset.variables[self.insitu_variable][index_mu]

        if 'satellite_time' in self.dataset.variables and 'insitu_time' in self.dataset.variables:
            from datetime import datetime as dt
            satellite_time = self.dataset.variables['satellite_time'][index_mu]
            insitu_time = self.dataset.variables['insitu_time'][index_mu]
            if self.fix_time_sat:
                satellite_time_day = dt.utcfromtimestamp(float(satellite_time)).strftime('%Y-%m-%d')
                satellite_time = dt.strptime(f'{satellite_time_day}T{self.fix_time_sat}','%Y-%m-%dT%H:%M').replace(tzinfo=pytz.utc).timestamp()
                #print('-->',dt.utcfromtimestamp(satellite_time))
            # ##temporal, change satellite time
            # satellite_time_r = dt.utcfromtimestamp(float(satellite_time))
            # satellite_time_r = satellite_time_r.replace(hour=13,tzinfo = pytz.UTC)
            # #print('estamos aqui...',satellite_time_r)
            # satellite_time = satellite_time_r.timestamp()
            # ##fin temporal

            nid = insitu_time.shape[0]
            satellite_time_n = np.ma.repeat(satellite_time,nid)
            time_diff_mu = np.ma.abs(insitu_time-satellite_time_n)
        else:
            time_diff_mu = np.ma.zeros(data_ins_mu.shape)

        time_diff_mu_copy = np.ma.masked_where(data_ins_mu.mask,time_diff_mu,True)
        n_non_masked = np.ma.count(time_diff_mu_copy)
        if n_non_masked>=1:
            if self.temporal_sampling_method=='closest':
                insitu_id = np.ma.argmin(time_diff_mu_copy)
                value = data_ins_mu[insitu_id]
                time_diff_value = time_diff_mu_copy[insitu_id]
            if self.temporal_sampling_method=='only_single':
                if n_non_masked==1:
                    insitu_id = np.ma.argmin(time_diff_mu_copy)
                    value = data_ins_mu[insitu_id]
                    time_diff_value = time_diff_mu_copy[insitu_id]
                else:
                    time_diff_value = -999.0
                    insitu_id = -999
                    value = -999.0
        else:
            time_diff_value = -999.0
            insitu_id = -999
            value = -999.0

        return time_diff_mu,insitu_id,time_diff_value,value

    def get_dimensions(self):
        # Dimensions
        nrows = self.dataset[self.satellite_variable].shape[1]
        ncols = self.dataset[self.satellite_variable].shape[2]
        central_r = int(np.floor(nrows / 2))
        central_c = int(np.floor(ncols / 2))
        r_s = central_r - int(np.floor(self.window_size / 2))  # starting row
        r_e = central_r + int(np.floor(self.window_size / 2)) + 1  # ending row
        c_s = central_c - int(np.floor(self.window_size / 2))  # starting col
        c_e = central_c + int(np.floor(self.window_size / 2)) + 1  # ending col
        return central_r, central_c, r_s, r_e, c_s, c_e

    def close_dataset(self):
        if self.dataset is not None: self.dataset.close()
