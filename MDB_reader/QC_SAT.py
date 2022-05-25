import numpy as np
import COMMON.Class_Flags_OLCI as flag
import math
import BSC_QAA.bsc_qaa_EUMETSAT as bsc_qaa
from QCBase import QCBase


class QC_SAT:

    def __init__(self, satellite_rrs, sat_bands, satellite_flag, ac_processor):
        self.name = ''
        self.satellite_rrs = satellite_rrs
        self.sat_bands = sat_bands
        self.nmu = self.satellite_rrs.shape[0]
        self.nbands = self.satellite_rrs.shape[1]

        self.stat_value = 'avg'

        self.window_size = 3
        self.min_valid_pixels = 9
        self.use_Bailey_Werdell = False  # if true, minimun number of valid pixels is established
        # as min_porc_valid_pixels (default 50%) +1 of NTPW
        # if false, minimun number of valid pixels==min_valid_pixels

        self.apply_outliers = True
        self.outliers_info = {
            'central_stat': 'avg',
            'dispersion_stat': 'std',
            'factor': 1.5
        }

        self.NTP = self.window_size * self.window_size  # total number of pixels
        self.NTPW = self.NTP  # total number of water pixels (excluding land/inland waters), could vary with MU
        self.NVP = 0  # number of valid pixels (excluding flag pixels), varies with MU
        self.flag_mask = None  # mask based on flagging

        self.info_flag = {}
        if satellite_flag is not None and ac_processor is not None:
            flag_list, flag_land, flag_inlandwater = self.get_flag_defaults(ac_processor)
            self.info_flag[satellite_flag.name] = {
                'variable': satellite_flag,
                'flag_list': flag_list,
                'flag_land': flag_land,
                'flag_inlandwater': flag_inlandwater,
                'ac_processor': ac_processor,
                'nflagged': 0,
                'flag_stats': None
            }

        self.th_masks = []
        self.check_statistics = []

        self.invalid_mask = {}
        for sat_index in range(self.nbands):
            sat_index_str = str(sat_index)
            self.invalid_mask[sat_index_str] = {
                'wavelength': self.sat_bands[sat_index],
                'apply_mask': True,
                'n_masked': 0,
                'ref': f'rrs_{self.sat_bands[sat_index]:.0f}_invalid'
            }

        self.statistics = {}
        for sat_index in range(self.nbands):
            sat_index_str = str(sat_index)
            stat_list = {
                'n_values': 0,
                'avg': 0,
                'std': 0,
                'median': 0,
                'min': 0,
                'max': 0,
                'CV': 0
            }
            self.statistics[sat_index_str] = {
                'wavelength': self.sat_bands[sat_index],
                'without_outliers': stat_list,
                'with_outliers': stat_list
            }
        # print(self.statistics[sat_index_str]['without_outliers'])

        self.max_diff_wl = 5

        self.apply_band_shifting = False
        self.wl_ref = None

    def set_window_size(self, wsize):
        self.window_size = wsize
        self.NTP = self.window_size * self.window_size
        self.NTPW = self.NTP
        if self.min_valid_pixels > self.NTP:
            self.min_valid_pixels = self.NTP

    def get_wl_sat_list_from_wlreflist(self, wlref):
        wllist = []
        for wl in wlref:
            index = self.get_index_sat_from_wlvalue(wl)
            if index >= 0:
                wllist.append(self.sat_bands[index])
        return wllist

    def prepare_new_match_up(self):
        self.NTPW = self.NTP  # total number of water pixels (excluding land/inland waters), could vary with MU
        self.NVP = 0  # number of valid pixels (excluding flag pixels), varies with MU
        self.flag_mask = None  # mask based on flagging
        self.statistics = {}
        for sat_index in range(self.nbands):
            sat_index_str = str(sat_index)
            stat_list = {
                'n_values': 0,
                'avg': 0,
                'std': 0,
                'median': 0,
                'min': 0,
                'max': 0,
                'CV': 0
            }
            self.statistics[sat_index_str] = {
                'wavelength': self.sat_bands[sat_index],
                'without_outliers': stat_list,
                'with_outliers': stat_list
            }

    def compute_statistics(self, index_mu):
        cond_min_pixels = self.compute_masks_and_check_roi(index_mu)
        if not cond_min_pixels:
            return False

        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        for sat_index in range(self.nbands):
            sat_index_str = str(sat_index)
            rrs_here = self.satellite_rrs[index_mu, sat_index, r_s:r_e, c_s:c_e]
            rrs_valid = rrs_here[self.flag_mask == 0]
            stats = self.compute_statistics_impl(self.statistics[sat_index_str]['without_outliers'], rrs_valid)
            self.statistics[sat_index_str]['without_outliers'] = stats
            # self.statistics[sat_index_str]['without_outliers']['avg'] = np.mean(rrs_valid)
            # self.statistics[sat_index_str]['without_outliers']['std'] = np.std(rrs_valid)

            if self.apply_outliers:
                cvalue = self.statistics[sat_index_str]['without_outliers'][self.outliers_info['central_stat']]
                dvalue = self.statistics[sat_index_str]['without_outliers'][self.outliers_info['dispersion_stat']]
                min_th = cvalue - (dvalue * self.outliers_info['factor'])
                max_th = cvalue + (dvalue * self.outliers_info['factor'])
                mask_outliers = np.zeros(rrs_valid.shape)
                mask_outliers[rrs_valid > max_th] = 1
                mask_outliers[rrs_valid < min_th] = 1
                n_outliers = np.sum(mask_outliers)
                # self.statistics[sat_index_str]['n_good'] = self.NVP - n_outliers
                if n_outliers > 0:
                    rrs_valid = rrs_valid[mask_outliers == 0]
                stats = self.compute_statistics_impl(self.statistics[sat_index_str]['with_outliers'], rrs_valid)
                self.statistics[sat_index_str]['with_outliers'] = stats

        return True

    def compute_statistics_impl(self, stats, array):
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
            index_sat = str(check_stat['index_sat'])
            # type_stat = check_stat['type_stat']

            if index_sat in self.statistics:
                outliers_str = 'with_outliers'
                if not check_stat['with_outliers']:
                    outliers_str = 'without_outliers'
                val_here = self.statistics[index_sat][outliers_str][check_stat['type_stat']]
                # print(val_here, check_stat['value_th'])
                if check_stat['type_th'] == 'greater' and val_here > check_stat['value_th']:
                    CHECK = False
                if check_stat['type_th'] == 'lower' and val_here < check_stat['value_th']:
                    CHECK = False

        return CHECK

    def get_match_up_values(self, index_mu):
        self.prepare_new_match_up()
        cond_min_pixels = self.compute_masks_and_check_roi(index_mu)

        cond_stats = False
        valid_mu = False

        wl_orig = []
        if self.wl_ref is None:
            indexes_bands = range(self.nbands)
            wl_orig = self.sat_bands
        else:
            indexes_bands = []
            for wl in self.wl_ref:
                index = self.get_index_sat_from_wlvalue(wl)
                wl_orig.append(self.sat_bands[index])
                indexes_bands.append(index)

        values = [0] * len(indexes_bands)

        if cond_min_pixels:
            outliers_str = 'without_outliers'
            if self.apply_outliers:
                outliers_str = 'with_outliers'
            self.compute_statistics(index_mu)

            cond_stats = self.do_check_statistics()
            if cond_stats:
                valid_mu = True
            for idx in range(len(indexes_bands)):
                sat_index = indexes_bands[idx]
                sat_index_str = str(sat_index)
                values[idx] = self.statistics[sat_index_str][outliers_str][self.stat_value]
            if self.apply_band_shifting:
                values = bsc_qaa.bsc_qaa(values, wl_orig, self.wl_ref)

        return cond_min_pixels, cond_stats, valid_mu, values

    def compute_masks_and_check_roi(self, index_mu):
        land = self.compute_flag_masks(index_mu)
        self.compute_invalid_masks(index_mu)
        self.compute_th_masks(index_mu)
        self.NVP = self.NTP - np.sum(self.flag_mask)
        self.NTPW = self.NTP - np.sum(land, axis=(0, 1))
        # print(f'[INFO] Index mu: {index_mu}')
        # print(f'[INFO] Number total of pixels: {self.NTP}')
        # print(f'[INFO] Water pixels: {self.NTPW}')
        # print(f'[INFO] Valid (no-flag) pixels: {self.NVP}')

        min_valid_pixels = self.min_valid_pixels
        if self.use_Bailey_Werdell:
            min_valid_pixels = math.floor(0.50 * self.NTPW) + 1

        cond_min_pixels = False
        if self.NVP >= min_valid_pixels:
            cond_min_pixels = True

        return cond_min_pixels

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

    def compute_flag_stats(self, index_mu):
        for flag_band in self.info_flag.keys():
            self.compute_flag_stats_impl(index_mu, flag_band)

    def compute_th_masks(self, index_mu):
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        mask_thershold = np.zeros((self.window_size, self.window_size), dtype=np.uint64)
        for idx in range(len(self.th_masks)):
            th_mask = self.th_masks[idx]
            rrs_here = self.satellite_rrs[index_mu, th_mask['index_sat'], r_s:r_e, c_s:c_e]
            mask_thershold_here = np.zeros(rrs_here.shape, dtype=np.uint64)
            if th_mask['type_th'] == 'greater':
                mask_thershold_here[rrs_here > th_mask['value_th']] = 1
            elif th_mask['type_th'] == 'lower':
                mask_thershold_here[rrs_here < th_mask['value_th']] = 1
            n_masked = np.sum(mask_thershold_here)
            th_mask['n_masked'] = n_masked
            self.th_masks[idx] = th_mask
            mask_thershold = mask_thershold + mask_thershold_here

        if self.flag_mask is None:
            self.flag_mask = mask_thershold

        self.flag_mask[mask_thershold > 0] = 1

    def compute_invalid_masks(self, index_mu):
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        mask_invalid = np.zeros((self.window_size, self.window_size), dtype=np.uint64)
        for sat_index in range(self.nbands):
            sat_index_str = str(sat_index)
            if self.invalid_mask[sat_index_str]['apply_mask']:
                rrshere = self.satellite_rrs[index_mu, sat_index, r_s:r_e, c_s:c_e]

                mask_invalid_here = np.zeros(rrshere.shape, dtype=np.uint64)
                mask_invalid_here[rrshere.mask] = 1
                n_masked = np.sum(mask_invalid_here)
                self.invalid_mask[sat_index_str]['n_masked'] = n_masked
                mask_invalid = mask_invalid + mask_invalid_here

        if self.flag_mask is None:
            self.flag_mask = mask_invalid

        self.flag_mask[mask_invalid > 0] = 1

        # for key in self.invalid_mask:
        #     print(self.invalid_mask[key]['ref'],self.invalid_mask[key]['n_masked'])

    ##ADDING QUALITY CONTROL PROTOCOLS------------------------------
    # Add a thershold mask.
    # index_sat: index band (if -1, index_sat is obtained from wl_sat)
    # wl_sat: wavelength (used for computing index_sat)
    # value_th: threshold
    # type_th:[greater, lower]
    def add_theshold_mask(self, index_sat, wl_sat, value_th, type_th):
        if index_sat == -1:
            index_sat = self.get_index_sat_from_wlvalue(wl_sat)
        if index_sat < 0:
            return
        if index_sat >= self.nbands:
            return

        th_mask = {
            'index_sat': index_sat,
            'value_th': value_th,
            'type_th': type_th,
            'n_masked': 0
        }
        self.th_masks.append(th_mask)

    def add_band_statistics(self, index_sat, wl_sat, type_stat, with_outliers, value_th, type_th):
        if index_sat == -1:
            index_sat = self.get_index_sat_from_wlvalue(wl_sat)
        if index_sat < 0:
            return
        if index_sat >= self.nbands:
            return

        check_val = {
            'index_sat': index_sat,
            'type_stat': type_stat,
            'with_outliers': with_outliers,
            'value_th': value_th,
            'type_th': type_th
        }
        self.check_statistics.append(check_val)

    ##IMPLEMENTATIONS-----------------------------------
    def compute_flag_stats_impl(self, index_mu, flag_band):
        if index_mu < 0 or index_mu >= self.nmu:
            return
        if flag_band not in self.info_flag.keys():
            return
        if self.info_flag[flag_band]['flag_stats'] is not None:
            return
        satellite_flag = self.info_flag[flag_band]['variable']
        if satellite_flag is None:
            return

        if self.info_flag[flag_band]['ac_processor'] == 'POLYMER':
            flagging = flag.Class_Flags_Polymer(satellite_flag.flag_masks, satellite_flag.flag_meanings)
        else:
            flagging = flag.Class_Flags_OLCI(satellite_flag.flag_masks, satellite_flag.flag_meanings)
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        satellite_flag_band = satellite_flag[index_mu, r_s:r_e, c_s:c_e]
        flag_list_here = str.split(satellite_flag.flag_meanings, ' ')
        flag_stats = {}
        for flag_here in flag_list_here:
            flag_ref = f'{flag_band}.{flag_here}'

            mask_here = flagging.Mask(satellite_flag_band, ([flag_here]))
            mask_here[np.where(mask_here != 0)] = 1
            nflagg_here = np.sum(mask_here)
            flag_stats[flag_ref] = nflagg_here

        self.info_flag[flag_band]['flag_stats'] = flag_stats

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
            flag_mask = np.zeros((self.window_size, self.window_size), dtype=np.uint64)
            return flag_mask, land

        satellite_flag_band = satellite_flag[index_mu, r_s:r_e, c_s:c_e]
        if self.info_flag[flag_band]['ac_processor'] == 'POLYMER':
            flagging = flag.Class_Flags_Polymer(satellite_flag.flag_masks, satellite_flag.flag_meanings)
            flag_mask = flagging.MaskGeneral(satellite_flag_band)
            flag_mask[np.where(flag_mask != 0)] = 1
        else:
            flagging = flag.Class_Flags_OLCI(satellite_flag.flag_masks, satellite_flag.flag_meanings)
            if self.info_flag[flag_band]['ac_processor'] == 'C2RCC':  # C2RCC FLAGS
                if index_mu == 0:
                    print(satellite_flag.flag_meanings)
                    print(satellite_flag.flag_masks)
                valuePE = np.uint64(2147483648)
                flag_mask = np.ones(satellite_flag_band.shape, dtype=np.uint64)
                flag_mask[satellite_flag_band == valuePE] = 0
            else:
                flag_mask = flagging.Mask(satellite_flag_band, self.info_flag[flag_band]['flag_list'])
                flag_mask[np.where(flag_mask != 0)] = 1

        flag_land = self.info_flag[flag_band]['flag_land']
        flag_inlandwater = self.info_flag[flag_band]['flag_inlandwater']

        if flag_land is not None:
            land = flagging.Mask(satellite_flag_band, ([flag_land]))
            land[np.where(land != 0)] = 1
            if flag_inlandwater is not None:
                inland_w = flagging.Mask(satellite_flag_band, ([flag_inlandwater]))
                land[np.where(inland_w != 0)] = 0

        return flag_mask, land

    # self.invalid_mask[sat_index_str] = {
    #     'wavelength': self.invalid_mask[sat_index],
    #     'apply_mask': True,
    #     'n_masked': 0
    # }

    ##UTILITIES/DEFAUTLS-----------------------------------
    def get_index_sat_from_wlvalue(self, wl_sat):
        index_sat = np.argmin(np.abs(wl_sat - self.sat_bands))
        if np.abs(wl_sat - self.sat_bands[index_sat]) > self.max_diff_wl:
            index_sat = -1
        return index_sat

    def get_dimensions(self):
        # Dimensions
        nrows = self.satellite_rrs.shape[2]
        ncols = self.satellite_rrs.shape[3]
        central_r = int(np.floor(nrows / 2))
        central_c = int(np.floor(ncols / 2))
        r_s = central_r - int(np.floor(self.window_size / 2))  # starting row
        r_e = central_r + int(np.floor(self.window_size / 2)) + 1  # ending row
        c_s = central_c - int(np.floor(self.window_size / 2))  # starting col
        c_e = central_c + int(np.floor(self.window_size / 2)) + 1  # ending col
        return central_r, central_c, r_s, r_e, c_s, c_e

    def get_flag_defaults(self, ac_processor):
        flag_list = None
        flag_land = None
        flag_inlandwaters = None
        if ac_processor == 'STANDARD':
            flag_list = 'CLOUD,CLOUD_AMBIGUOUS,CLOUD_MARGIN,INVALID,COSMETIC,SATURATED,SUSPECT,HISOLZEN,HIGHGLINT,SNOW_ICE,AC_FAIL,WHITECAPS,RWNEG_O2,RWNEG_O3,RWNEG_O4,RWNEG_O5,RWNEG_O6,RWNEG_O7,RWNEG_O8'
            flag_land = 'LAND'
            flag_inlandwaters = 'INLAND_WATER'
        if ac_processor == 'POLYMER':
            # flag_list = 'LAND,CLOUD_BASE,L1_INVALID,NEGATIVE_BB,OUT_OF_BOUNDS,EXCEPTION,THICK_AEROSOL,HIGH_AIR_MASS,EXTERNAL_MASK'
            flag_list = 'LAND,CLOUD_BASE'
            flag_land = 'LAND'

        if ac_processor == 'C2RCC':
            flag_list = 'Rtosa_OOS, Rtosa_OOR, Rhow_OOR, Cloud_risk, Iop_OOR, Apig_at_max, Adet_at_max, Agelb_at_max, Bpart_at_max, Bwit_at_max, Apig_at_min, Adet_at_min, Agelb_at_min, Bpart_at_min, Bwit_at_min, Rhow_OOS, Kd489_OOR,Kdmin_OOR, Kd489_at_max, Kdmin_at_max'
            # flag_list =  ''
        if ac_processor == 'FUB':
            # flag_list = 'land,coastline,fresh_inland_water,bright,straylight_risk,invalid,cosmetic,duplicated,sun_glint_risk,dubious,saturated_Oa01,saturated_Oa02,saturated_Oa03,saturated_Oa04,saturated_Oa05,saturated_Oa06,saturated_Oa07,saturated_Oa08,saturated_Oa09,saturated_Oa10,saturated_Oa11,saturated_Oa12,saturated_Oa13,saturated_Oa14,saturated_Oa15,saturated_Oa16'
            flag_list = 'land,coastline,fresh_inland_water,bright,straylight_risk,invalid,cosmetic,sun_glint_risk,dubious,saturated_Oa01,saturated_Oa02,saturated_Oa03,saturated_Oa04,saturated_Oa05,saturated_Oa06,saturated_Oa07,saturated_Oa08,saturated_Oa09,saturated_Oa10,saturated_Oa11,saturated_Oa12,saturated_Oa13,saturated_Oa14,saturated_Oa15,saturated_Oa16'
            flag_land = 'land'
            flag_inlandwaters = 'fresh_inland_water'
            # flag_list = 'land,coastline,fresh_inland_water,bright,straylight_risk,invalid,cosmetic,duplicated,sun_glint_risk,dubious,saturated_Oa01,saturated_Oa02,saturated_Oa03,saturated_Oa04,saturated_Oa05,saturated_Oa06,saturated_Oa07,saturated_Oa08,saturated_Oa09,saturated_Oa10,saturated_Oa11,saturated_Oa12,saturated_Oa13,saturated_Oa14,saturated_Oa15,saturated_Oa16,saturated_Oa17,saturated_Oa18,saturated_Oa19,saturated_Oa20,saturated_Oa21'

        if flag_list is not None:
            flag_list = flag_list.replace(" ", "")
            flag_list = str.split(flag_list, ',')

        return flag_list, flag_land, flag_inlandwaters

    # eumetsat_defults: windows_size should be 3 (min_valid_pixels==9) o 5 (use_Bailey_Werdell=True)
    def set_eumetsat_defaults(self, window_size):
        self.stat_value = 'avg'

        self.window_size = window_size
        if window_size == 3:
            self.min_valid_pixels = 9
            self.use_Bailey_Werdell = False
        if window_size == 9:
            self.use_Bailey_Werdell = True

        self.apply_outliers = True
        self.outliers_info = {
            'central_stat': 'avg',
            'dispersion_stat': 'std',
            'factor': 1.5
        }

        self.add_band_statistics(-1, 560, 'CV', True, 20, 'greater')

    def set_qc_from_qcbase(self, qcbase):

        self.name = qcbase.sat_name
        self.stat_value = qcbase.sat_stat_value
        self.window_size = qcbase.sat_window_size
        self.min_valid_pixels = qcbase.sat_min_valid_pixels
        self.use_Bailey_Werdell = qcbase.sat_use_Bailey_Werdell
        self.max_diff_wl = qcbase.sat_max_diff_wl
        self.apply_band_shifting = qcbase.sat_apply_band_shifting

        self.apply_outliers = qcbase.sat_apply_outliers
        self.outliers_info = qcbase.sat_outliers_info

        if len(qcbase.sat_th_mask) > 0:
            for thm in qcbase.sat_th_mask:
                self.add_theshold_mask(thm['index_sat'], thm['wl_sat'], thm['value_th'], thm['type_th'])

        if len(qcbase.sat_check_statistics) > 0:
            for cst in qcbase.sat_check_statistics:
                self.add_band_statistics(cst['index_sat'], cst['wl_sat'], cst['type_stat'], cst['with_outliers'],
                                         cst['value_th'], cst['type_th'])

        if len(qcbase.sat_info_flag) > 0:
            for info_f in qcbase.sat_info_flag:
                name_flag = info_f['sat_flag_name']
                if name_flag is None:
                    name_flag = self.info_flag.keys()[0]
                if name_flag is self.info_flag.keys():
                    if info_f['flag_list'] is not None:
                        self.info_flag[name_flag]['flag_list'] = info_f['flag_list']
                    if info_f['flag_land'] is not None:
                        self.info_flag[name_flag]['flag_land'] = info_f['flag_land']
                    if info_f['flag_inlandwater'] is not None:
                        self.info_flag[name_flag]['flag_inlandwater'] = info_f['flag_land']
