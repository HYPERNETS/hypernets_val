class QCBase:
    def __init__(self):
        # Satellite quality control
        self.sat_name = ''
        self.sat_stat_value = 'avg'
        self.sat_window_size = 3
        self.sat_min_valid_pixels = 9
        self.sat_use_Bailey_Werdell = False

        self.sat_max_diff_wl = 5
        self.sat_apply_band_shifting = False

        self.sat_apply_outliers = True
        self.sat_outliers_info = {
            'central_stat': 'avg',
            'dispersion_stat': 'std',
            'factor': 1.5
        }
        self.sat_info_flag = []
        self.sat_th_mask = []
        self.sat_check_statistics = []

    def set_sat_outliers_info(self, central_stat, dispersion_stat, factor):
        self.sat_apply_outliers = True
        self.sat_outliers_info['central_stat'] = central_stat
        self.sat_outliers_info['dispersion_stat'] = dispersion_stat
        self.sat_outliers_info['factor'] = factor

    # sat_flag_name could be None to use default
    def add_sat_info_flag(self, sat_flag_name, flag_list, flag_land, flag_inlandwater):
        sat_info_flag_here = {
            'sat_flag_name': sat_flag_name,
            'flag_list': flag_list,
            'flag_land': flag_land,
            'flag_inlandwater': flag_inlandwater
        }
        self.sat_info_flag.append(sat_info_flag_here)

    def add_sat_theshold_mask(self, index_sat, wl_sat, value_th, type_th):
        th_mask_here = {
            'index_sat': index_sat,
            'wl_sat': wl_sat,
            'value_th': value_th,
            'type_th': type_th
        }
        self.sat_th_mask.append(th_mask_here)

    def add_sat_band_statistics(self, index_sat, wl_sat, type_stat, with_outliers, value_th, type_th):
        check_val_here = {
            'index_sat': index_sat,
            'wl_sat': wl_sat,
            'type_stat': type_stat,
            'with_outliers': with_outliers,
            'value_th': value_th,
            'type_th': type_th
        }
        self.sat_check_statistics.append(check_val_here)

    def set_sat_eumetsat_defaults(self, window_size):
        self.sat_stat_value = 'avg'
        self.sat_window_size = window_size
        if window_size == 3:
            self.sat_min_valid_pixels = 9
            self.sat_use_Bailey_Werdell = False
        if window_size == 9:
            self.sat_use_Bailey_Werdell = True
        self.set_sat_outliers_info('avg', 'std', 1.5)
        self.add_sat_band_statistics(-1, 560, 'CV', True, 20, 'greater')
