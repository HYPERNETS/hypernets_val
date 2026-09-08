import math
import numpy as np
import COMMON.flag_functions as ffs
import COMMON.Class_Flags_OLCI as flag
import BSC_QAA.bsc_qaa_EUMETSAT as bsc_qaa





class QC_SAT:

    def __init__(self, dataset,versose=False):

        self.dataset = dataset
        self.wl_list = None
        self.indices_wl = None
        self.verbose = versose

        ##basic information
        self.bands_variable = 'satellite_bands'
        self.spectral_variable = 'satellite_Rrs'
        self.spectral_variable_unc = 'satellite_Rrs_unc'
        self.window_size = 3
        self.min_valid_pixels = 9
        self.min_valid_porc = 50.0

        self.use_min_valid_porc = False  # if true, minimum number of valid pixels is established as min_porc_valid_pixels (default 50%) +1 of NP
        self.min_valid_porc_with_valid = True #if valid, NP is the NTWP (number of water pixels), if false, NP, is the NTP (number of total pixels)

        self.stat_value = 'avg'

        self.outliers_info = {'apply': True, 'central_stat': 'avg','dispersion_stat': 'std','factor': 1.5}

        self.max_diff_wl = 5.0

        self.flag_land = [None]*2
        self.flag_inland_water = [None]*2

        self.filter_invalid = 'all'
        self.filter_flag = []
        self.filter_spectral_th = []
        self.filter_var_th = []
        self.filter_macropixel_spectral = []
        self.filter_macropixel_var = []

        self.potential_stat_values = ['avg', 'std', 'min', 'max', 'median']
        self.central_stat_values = ['avg', 'median']
        self.dispersion_stat_values = ['std', 'iqr']
        self.th_types = ['greater', 'gt', 'gte', 'lower', 'lt', 'lte']
        self.macropixel_spatial_stats = ['avg', 'median', 'std', 'iqr', 'min', 'max', 'CV']
        self.macropixel_spectral_stats = ['all', 'any', 'avg', 'median', 'std', 'iqr', 'min', 'max', 'CV']
        self.filter_invalid_options  = ['all','any']

        ##dimensions
        self.n_bands = None ##defined for wl_list
        self.n_rows = None
        self.n_cols = None
        self.n_mu = None

        self.is_valid = True

    def set_basic_info(self,options_config):
        self.bands_variable = options_config['bands_variable']
        self.spectral_variable = options_config['spectral_variable']
        self.spectral_variable_unc = options_config['spectral_variable_unc']
        self.window_size = options_config['window_size']
        self.min_valid_pixels = options_config['min_valid_pixels']
        self.min_valid_porc = options_config['min_valid_porc']
        self.use_min_valid_porc = options_config['use_min_valid_porc']
        self.min_valid_porc_with_valid = options_config['min_valid_porc_with_valid']

        self.stat_value = options_config['stat_value']
        self.outliers_info['apply'] = options_config['apply_outliers']
        self.max_diff_wl  = options_config['max_diff_wl']
        if self.outliers_info['apply']:
            o_info = options_config['outliers_info']
            self.outliers_info['central_stat'] = o_info[0]
            self.outliers_info['dispersion_stat'] = o_info[1]
            self.outliers_info['factor'] = o_info[2]

        self.flag_land = options_config['flag_land']
        self.flag_inland_water = options_config['flag_inland_water']

        self.filter_invalid = options_config['filter_invalid']

    def set_wl_list(self,options_config):
        if not self.bands_variable in self.dataset.variables:
            print(f'[ERROR][QC_SAT] {self.bands_variable} variable is not available in the NetCDF dataset')
            return

        original_bands = self.dataset.variables[self.bands_variable][:]
        wl_list = options_config['wl_list']
        wl_min = options_config['wl_min']
        wl_max = options_config['wl_max']
        if wl_list is not None and (wl_min is not None or wl_max is not None):
            print(f'[WARNING] As wl_list is given, wl_min and wl_max values are not used for setting the wavelength list')


        if wl_list is None:
            if wl_min is not None or wl_max is not None:
                if wl_min is None:
                    wl_min = np.min(original_bands)
                if wl_max  is None:
                    wl_max = np.max(original_bands)
                if wl_max<wl_min:
                    print(f'[ERROR][QC_SAT] wl_max ({wl_max}) must be greater or equal than wl_min ({wl_min})')
                    return
                wl_list = original_bands[(original_bands>=wl_min) & (original_bands<=wl_max)]
                print(f'[INFO][QC_SAT] wl_list set to {len(wl_list)} bands between {np.min(wl_list)} and {np.max(wl_list)}')
            elif wl_min is None and wl_max is None:
                wl_list = original_bands
                print(f'[INFO][QC_SAT] wl_list set to the original satellite bands with {len(wl_list)} bands between {np.min(wl_list)} and {np.max(wl_list)}')


        wl_list = np.array(wl_list)

        n_original_bands = len(original_bands)
        n_list = len(wl_list)

        wl_list_m = np.repeat(wl_list,n_original_bands).reshape((n_list,n_original_bands))
        wl_original_m = np.tile(original_bands,n_list).reshape((n_list,n_original_bands))

        diff_wl = np.abs(wl_list_m-wl_original_m)
        min_diff_wl = np.min(diff_wl,axis=1)
        if np.max(min_diff_wl)>=self.max_diff_wl:
            print(f'[ERROR][QC_SAT] Some bands given in the parameter wl_list are not available as wavelength difference with the satellite original bands is greater than the allowed maximum of {self.max_diff_wl} nm')
            no_valid_bands = wl_list[min_diff_wl>self.max_diff_wl]
            print(f'Please review the following bands: ')
            for no_valid_band in no_valid_bands:
                print(f'  {no_valid_band} nm')
            print(f'Or you can also modify the allowed  maximum difference using the max_diff_wl parameter in QC_SAT')
            return


        self.indices_wl = np.argmin(diff_wl,axis=1)
        self.wl_list = original_bands[self.indices_wl]




        # self.name = ''
        # self.satellite_rrs = satellite_rrs
        # self.satellite_rrs_unc  = None
        # self.sat_bands = sat_bands
        #
        # self.qc_sat_results = {}
        #
        # self.pi_multiplied = [False] * len(self.sat_bands)
        # self.pi_divided = [False] * len(self.sat_bands)
        # self.nmu = self.satellite_rrs.shape[0]
        # self.nbands = self.satellite_rrs.shape[1]
        #
        # self.stat_value = 'avg'
        #
        # self.window_size = 3
        # self.min_valid_pixels = 9

        #
        # self.apply_outliers = True
        # self.outliers_info = {
        #     'central_stat': 'avg',
        #     'dispersion_stat': 'std',
        #     'factor': 1.5
        # }
        #
        # self.NTP = self.window_size * self.window_size  # total number of pixels
        # self.NTPW = self.NTP  # total number of water pixels (excluding land/inland waters), could vary with MU
        # self.NVP = 0  # number of valid pixels (excluding flag pixels), varies with MU
        # self.flag_mask = None  # mask based on flagging
        # self.ac_processor = ac_processor
        #
        # self.info_flag = {}
        # if satellite_flag is not None and ac_processor is not None:
        #     flag_list, flag_land, flag_inlandwater = self.get_flag_defaults(ac_processor)
        #     self.info_flag[satellite_flag.name] = {
        #         'variable': satellite_flag,
        #         'flag_list': flag_list,
        #         'flag_land': flag_land,
        #         'flag_inlandwater': flag_inlandwater,
        #         'ac_processor': ac_processor,
        #         'nflagged': 0,
        #         'flag_stats': None
        #     }
        #
        # self.th_masks = []
        # self.check_statistics = []
        #
        # ##values for statistis no rrs
        # self.check_statistics_norrs = []
        # self.statistics_norrs = {}
        # self.ncdataset = None
        #
        # self.invalid_mask = {}
        # for sat_index in range(self.nbands):
        #     sat_index_str = str(sat_index)
        #     self.invalid_mask[sat_index_str] = {
        #         'wavelength': self.sat_bands[sat_index],
        #         'apply_mask': True,
        #         'n_masked': 0,
        #         'ref': f'rrs_{self.sat_bands[sat_index]:.0f}_invalid'
        #     }
        #
        # self.statistics = {}
        # self.statistics_unc = {}
        # for sat_index in range(self.nbands):
        #     sat_index_str = str(sat_index)
        #     stat_list = {
        #         'n_values': 0,
        #         'avg': 0,
        #         'std': 0,
        #         'median': 0,
        #         'min': 0,
        #         'max': 0,
        #         'CV': 0
        #     }
        #     self.statistics[sat_index_str] = {
        #         'wavelength': self.sat_bands[sat_index],
        #         'without_outliers': stat_list,
        #         'with_outliers': stat_list
        #     }
        #     self.statistics_unc[sat_index_str] = {
        #         'wavelength': self.sat_bands[sat_index],
        #         'without_outliers': stat_list,
        #         'with_outliers': stat_list
        #     }
        # # print(self.statistics[sat_index_str]['without_outliers'])
        #
        # self.max_diff_wl = 10
        #
        # self.apply_band_shifting = False
        # self.wl_ref = None
        # self.mu_invalid_list = []
        #
        # self.indices_valid_bands = None
        #
        # self.apply_olci_gains = False
        # self.olci_gains_s3a = {
        #     '400': 0.97546,
        #     '412.5': 0.97406,
        #     '442.5': 0.97492,
        #     '490': 0.9689,
        #     '510': 0.97184,
        #     '560': 0.97571,
        #     '620': 0.98001,
        #     '665': 0.97834,
        #     '673.75': 0.9786,
        #     '681.25': 0.97908,
        #     '708.75': 0.98013,
        #     '753.75': 0.98552,
        #     '778.75': 0.98772,
        #     '865': 0.986,
        #     '885': 0.98657,
        #     '1020': 0.91316
        # }
        # self.olci_gains_s3b = {
        #     '400': 0.99458,
        #     '412.5': 0.9901,
        #     '442.5': 0.99221,
        #     '490': 0.9862,
        #     '510': 0.98898,
        #     '560': 0.99114,
        #     '620': 0.99769,
        #     '665': 0.99684,
        #     '673.75': 0.99716,
        #     '681.25': 0.99802,
        #     '708.75': 0.99782,
        #     '753.75': 1.00163,
        #     '778.75': 1.00259,
        #     '865': 1,
        #     '885': 1.00089,
        #     '1020': 0.94064
        # }

    def set_filter_flag(self,options_config, key_values = None):
        self.filter_flag = get_filter_list(options_config,'filter_flag_',key_values=key_values)

    def set_filter_spectral_th(self,options_config,key_values = None):
        self.filter_spectral_th = get_filter_list(options_config,'filter_spectral_th_',key_values=key_values)

    def set_filter_var_th(self,options_config,key_values = None):
        self.filter_var_th = get_filter_list(options_config,'filter_var_th_',key_values=key_values)

    def set_filter_macropixel_spectral(self,options_config,key_values = None):
        self.filter_macropixel_spectral = get_filter_list(options_config,'filter_macropixel_spectral_',key_values=key_values)

    def set_filter_macropixel_var(self,options_config,key_values = None):
        self.filter_macropixel_var = get_filter_list(options_config,'filter_macropixel_var_',key_values=key_values)

    def check_parameters(self,potential_stat_values=None):
        check_qc = self.check_spectral_variable(self.spectral_variable,self.bands_variable)
        n_rows, n_cols = -1, -1
        if check_qc:
            n_rows = self.dataset.variables[self.spectral_variable].shape[2]
            n_cols = self.dataset.variables[self.spectral_variable].shape[3]

        if self.wl_list is None or self.indices_wl is None:
            check_qc = False

        check_w = self.check_window_size(self.window_size,n_rows,n_cols)
        if check_w:
            n_all = self.window_size * self.window_size
            if not self.use_min_valid_porc and self.min_valid_pixels > n_all:
                print(f'[ERROR][QC_SAT] min_valid_pixels {min_valid_pixels} should be lower or equal to {n_all} for a window size of {self.window_size} * {self.window_size} pixels')
                check_qc = False
        else:
            check_qc = False

        if self.use_min_valid_porc and (self.min_valid_porc<0 or self.min_valid_porc>100):
            print(f'[ERROR][QC_SAT] min_valid_porc {self.min_valid_porc} should be between 0 and 100')
            check_qc = False

        if self.filter_invalid not in self.filter_invalid_options:
            print(f'[ERROR][QC_SAT] filter_invalid {self.filter_invalid} should be one of {self.filter_invalid_options}')
            check_qc = False

        if potential_stat_values is None:
            potential_stat_values = self.potential_stat_values

        if self.stat_value not in potential_stat_values:
            print(f'[ERROR][QC_SAT] Stat value ({self.stat_value}) is not valid, should be one of {potential_stat_values}')
            check_qc = False

        if self.outliers_info['apply']:
            if self.outliers_info['central_stat'] not in self.central_stat_values:
                print(f'[ERROR][QC_SAT] central_stat in outliers_info should be one of {self.central_stat_values}')
                check_qc = False
            if self.outliers_info['dispersion_stat'] not in self.dispersion_stat_values:
                print(f'[ERROR][QC_SAT] dispersion_stat in outliers_info should be one of {self.dispersion_stat_values}')
                check_qc = False
            factor_str = self.outliers_info['factor']
            try:
                self.outliers_info['factor'] = float(factor_str)
            except Exception as ex:
                print(f'[ERROR][QC_SAT] Factor in outliers_info {factor_str} must be a float. Exception: {ex}')
                self.outliers_info['factor'] = None
                check_qc = False

        if self.flag_land is not None:
            if not self.check_flag_list(self.flag_land[0],[self.flag_land[1]]):
                check_qc = False

        if self.flag_inland_water is not None:
            if not self.check_flag_list(self.flag_inland_water[0],[self.flag_inland_water[1]]):
                check_qc = False

        if len(self.filter_flag)>0:
            for idx in range(len(self.filter_flag)):
                if self.filter_flag[idx]['name_var'] is None and self.filter_flag[idx]['ac_processor'] is not None:
                    self.filter_flag[idx]['name_var'] = None ##not implemented, to get name_var from default ac_processor
                if self.filter_flag[idx]['flag_list'] is None and self.filter_flag[idx]['ac_processor'] is not None:
                    self.filter_flag[idx]['flag_list'] = None ##not implemented, to get flag_list from default ac_processor
                if not self.check_flag_list(self.filter_flag[idx]['name_var'],self.filter_flag[idx]['flag_list']):
                    check_qc = False
                if self.filter_flag[idx]['window_size'] == -1:
                    self.filter_flag[idx]['window_size'] = self.window_size
                n_rows_here = self.dataset.variables[self.filter_flag[idx]['name_var']][:].shape[1]
                n_cols_here = self.dataset.variables[self.filter_flag[idx]['name_var']][:].shape[2]
                check_w = self.check_window_size(self.filter_flag[idx]['window_size'], n_rows_here, n_cols_here)
                if not check_w:
                    check_qc = False

        if len(self.filter_spectral_th)>0:
            for idx in range(len(self.filter_spectral_th)):
                check_b = self.check_spectral_variable(self.filter_spectral_th[idx]['name_var'],self.filter_spectral_th[idx]['name_var_wl'])
                if check_b:
                    wl_min_abs = np.min(self.dataset.variables[self.filter_spectral_th[idx]['name_var_wl']][:])-self.max_diff_wl
                    wl_max_abs = np.max(self.dataset.variables[self.filter_spectral_th[idx]['name_var_wl']][:])+self.max_diff_wl
                    wl_min_here = self.filter_spectral_th[idx]['wl_min']
                    wl_max_here = self.filter_spectral_th[idx]['wl_min']
                    if not isinstance(wl_min_here, float):
                        print(f'[ERROR][QC_SAT] wl_min {wl_min_here} for filter_spectral_th_{idx} should be a float value')
                        check_qc = False
                    if not isinstance(wl_max_here, float):
                        print(f'[ERROR][QC_SAT] wl_max {wl_max_here} for filter_spectral_th_{idx} should be a float value')
                        check_qc = False
                    if isinstance(wl_min_here,float) and isinstance(wl_max_here,float):
                        if wl_min_here>wl_max_here:
                            print(f'[ERROR][QC_SAT] wl_max {wl_max_here} should be greater or equal to wl_min {wl_min_here} for filter_spectral_th_{idx}')
                            check_qc = False
                        else:
                            if wl_min_here<wl_min_abs or wl_min_here>wl_max_abs:
                                print(f'[ERROR][QC_SAT] wl_min {wl_min_here} should be in the spectral range {wl_min_abs} - {wl_max_abs} for filter_spectral_th_{idx}')
                                check_qc = False
                            if wl_max_here<wl_min_abs or wl_max_here>wl_max_abs:
                                print(f'[ERROR][QC_SAT] wl_max {wl_max_here} should be in the spectral range {wl_min_abs} - {wl_max_abs} for filter_spectral_th_{idx}')
                                check_qc = False

                    if self.filter_spectral_th[idx]['window_size']==-1:
                        self.filter_spectral_th[idx]['window_size']= self.window_size
                    n_rows_here = self.dataset.variables[self.filter_spectral_th[idx]['name_var']][:].shape[2]
                    n_cols_here = self.dataset.variables[self.filter_spectral_th[idx]['name_var']][:].shape[3]
                    check_w = self.check_window_size(self.filter_spectral_th[idx]['window_size'],n_rows_here,n_cols_here)
                    if not check_w:
                        check_qc = False
                else:
                    check_qc = False
                if not isinstance(self.filter_spectral_th[idx]['th_value'], float):
                    print(f'[ERROR][QC_SAT] Threshold th_value {self.filter_spectral_th[idx]['th_value']} for filter_spectral_th_{idx} should be a float')
                    check_qc = False
                if not self.filter_spectral_th[idx]['th_type'] in self.th_types:
                    print(f'[ERROR][QC_SAT] Threshold filter type th_type {self.filter_spectral_th[idx]['th_type']} for filter_spectral_th_{idx} should one of {self.th_types}')
                    check_qc = False

        if len(self.filter_var_th)>0:
            for idx in range(len(self.filter_var_th)):
                check_b = self.check_non_spectral_variable(self.filter_var_th[idx]['name_var'])
                if check_b:
                    n_rows_here = self.dataset.variables[self.filter_var_th[idx]['name_var']][:].shape[1]
                    n_cols_here = self.dataset.variables[self.filter_var_th[idx]['name_var']][:].shape[2]
                    if self.filter_var_th[idx]['window_size']==-1:
                        self.filter_var_th[idx]['window_size']= self.window_size
                    check_w = self.check_window_size(self.filter_var_th[idx]['window_size'],n_rows_here,n_cols_here)
                    if not check_w:
                        check_qc = False
                else:
                    check_qc = False
                if not isinstance(self.filter_var_th[idx]['th_value'], float):
                    print(f'[ERROR][QC_SAT] Threshold th_value {self.filter_var_th[idx]['th_value']} for filter_var_th_{idx} should be a float')
                    check_qc = False
                if not self.filter_var_th[idx]['th_type'] in self.th_types:
                    print(f'[ERROR][QC_SAT] Threshold filter type th_type {self.filter_var_th[idx]['th_type']} for filter_var_th_{idx} should one of {self.th_types}')
                    check_qc = False

        if len(self.filter_macropixel_spectral)>0:
            for idx in range(len(self.filter_macropixel_spectral)):
                check_b = self.check_spectral_variable(self.filter_macropixel_spectral[idx]['name_var'],self.filter_macropixel_spectral[idx]['name_var_wl'])
                if check_b:
                    wl_min_abs = np.min(self.dataset.variables[self.filter_macropixel_spectral[idx]['name_var_wl']][:]) - self.max_diff_wl
                    wl_max_abs = np.max(self.dataset.variables[self.filter_macropixel_spectral[idx]['name_var_wl']][:]) + self.max_diff_wl
                    wl_min_here = self.filter_macropixel_spectral[idx]['wl_min']
                    wl_max_here = self.filter_macropixel_spectral[idx]['wl_min']
                    if not isinstance(wl_min_here, float):
                        print(f'[ERROR][QC_SAT] wl_min {wl_min_here} for filter_macropixel_spectral_{idx} should be a float value')
                        check_qc = False
                    if not isinstance(wl_max_here, float):
                        print(f'[ERROR][QC_SAT] wl_max {wl_max_here} for filter_macropixel_spectral_{idx} should be a float value')
                        check_qc = False
                    if isinstance(wl_min_here,float) and isinstance(wl_max_here,float):
                        if wl_min_here > wl_max_here:
                            print(f'[ERROR][QC_SAT] wl_max {wl_max_here} should be greater or equal to wl_min {wl_min_here} for filter_macropixel_spectral_{idx}')
                            check_qc = False
                        else:
                            if wl_min_here < wl_min_abs or wl_min_here > wl_max_abs:
                                print(f'[ERROR][QC_SAT] wl_min {wl_min_here} should be in the spectral range {wl_min_abs} - {wl_max_abs} for filter_macropixel_spectral_{idx}')
                                check_qc = False
                            if wl_min_here < wl_min_abs or wl_min_here > wl_max_abs:
                                print(f'[ERROR][QC_SAT] wl_min {wl_min_here} should be in the spectral range {wl_min_abs} - {wl_max_abs} for filter_macropixels_spectral_th_{idx}')
                                check_qc = False
                    if self.filter_macropixel_spectral[idx]['window_size']==-1:
                        self.filter_macropixel_spectral[idx]['window_size']= self.window_size
                    n_rows_here = self.dataset.variables[self.filter_macropixel_spectral[idx]['name_var']][:].shape[2]
                    n_cols_here = self.dataset.variables[self.filter_macropixel_spectral[idx]['name_var']][:].shape[3]
                    check_w = self.check_window_size(self.filter_macropixel_spectral[idx]['window_size'],n_rows_here,n_cols_here)
                    if not check_w:
                        check_qc = False
                else:
                    check_qc = False
                if not isinstance(self.filter_macropixel_spectral[idx]['th_value'], float):
                    print(f'[ERROR][QC_SAT] Threshold th_value {self.filter_macropixel_spectral[idx]['th_value']} for filter_macropixel_spectral_{idx} should be a float')
                    check_qc = False
                if not self.filter_macropixel_spectral[idx]['th_type'] in self.th_types:
                    print(f'[ERROR][QC_SAT] Threshold filter type th_type {self.filter_macropixel_spectral[idx]['th_type']} for filter_macropixel_spectral_{idx} should one of {self.th_types}')
                    check_qc = False
                if not self.filter_macropixel_spectral[idx]['spatial_stat'] in self.macropixel_spatial_stats:
                    print(f'[ERROR][QC_SAT] spatial_stat {self.filter_macropixel_spectral[idx]['spatial_stat']} for filter_macropixel_spectral_{idx} should one of {self.macropixel_spatial_stats}')
                    check_qc = False
                if not self.filter_macropixel_spectral[idx]['spectral_stat'] in self.macropixel_spectral_stats:
                    print(f'[ERROR][QC_SAT] spectral_stat {self.filter_macropixel_spectral[idx]['spectral_stat']} for filter_macropixel_spectral_{idx} should one of {self.macropixel_spectral_stats}')
                    check_qc = False
                if self.filter_macropixel_spectral[idx]['use_outliers']:
                    if self.filter_macropixel_spectral[idx]['outliers_central_stat'] not in self.central_stat_values:
                        print(f'[ERROR][QC_SAT] outliers_central_stat {self.filter_macropixel_spectral[idx]['outliers_central_stat']} for filter_macropixel_spectral_{idx} should be one of {self.central_stat_values}')
                        check_qc = False
                    if self.filter_macropixel_spectral[idx]['outliers_dispersion_stat'] not in self.dispersion_stat_values:
                        print(f'[ERROR][QC_SAT] outliers_dispersion_stat {self.filter_macropixel_spectral[idx]['outliers_dispersion_stat']} for filter_macropixel_spectral_{idx} be one of {self.dispersion_stat_values}')
                        check_qc = False
                    factor_str = self.filter_macropixel_spectral[idx]['outliers_factor']
                    try:
                        self.filter_macropixel_spectral[idx]['factor'] = float(factor_str)
                    except Exception as ex:
                        print(f'[ERROR][QC_SAT] outliers_factor {factor_str} for filter_macropixel_spectral_{idx}  must be a float number. Exception: {ex}')
                        self.filter_macropixel_spectral[idx]['factor'] = None
                        check_qc = False

        if len(self.filter_macropixel_var)>0:
            for idx in range(len(self.filter_macropixel_var)):
                check_b = self.check_non_spectral_variable(self.filter_macropixel_var[idx]['name_var'])
                if check_b:
                    n_rows_here = self.dataset.variables[self.filter_macropixel_var[idx]['name_var']][:].shape[1]
                    n_cols_here = self.dataset.variables[self.filter_macropixel_var[idx]['name_var']][:].shape[2]
                    if self.filter_macropixel_var[idx]['window_size'] == -1:
                        self.filter_macropixel_var[idx]['window_size'] = self.window_size
                    check_w = self.check_window_size(self.filter_macropixel_var[idx]['window_size'], n_rows_here, n_cols_here)
                    if not check_w:
                        check_qc = False

                else:
                    check_qc = False

                if not isinstance(self.filter_macropixel_var[idx]['th_value'], float):
                    print(f'[ERROR][QC_SAT] Threshold th_value {self.filter_macropixel_var[idx]['th_value']} for filter_macropixel_var_{idx} should be a float')
                    check_qc = False

                if not self.filter_macropixel_var[idx]['th_type'] in self.th_types:
                    print(f'[ERROR][QC_SAT] Threshold filter type th_type {self.filter_macropixel_var[idx]['th_type']} for filter_macropixel_var_{idx} should one of {self.th_types}')
                    check_qc = False

                if not self.filter_macropixel_var[idx]['spatial_stat'] in self.macropixel_spatial_stats:
                    print(f'[ERROR][QC_SAT] spatial_stat {self.filter_macropixel_var[idx]['spatial_stat']} for filter_macropixel_var_{idx} should one of {self.macropixel_spatial_stats}')
                    check_qc = False

                if self.filter_macropixel_var[idx]['use_outliers']:
                    if self.filter_macropixel_var[idx]['outliers_central_stat'] not in self.central_stat_values:
                        print(f'[ERROR][QC_SAT] outliers_central_stat {self.filter_macropixel_var[idx]['outliers_central_stat']} for filter_macropixel_var_{idx} should be one of {self.central_stat_values}')
                        check_qc = False
                    if self.filter_macropixel_var[idx]['outliers_dispersion_stat'] not in self.dispersion_stat_values:
                        print(f'[ERROR][QC_SAT] outliers_dispersion_stat {self.filter_macropixel_var[idx]['outliers_dispersion_stat']} for filter_macropixel_var_{idx} be one of {self.dispersion_stat_values}')
                        check_qc = False
                    factor_str = self.filter_macropixel_var[idx]['outliers_factor']
                    try:
                        self.filter_macropixel_var[idx]['factor'] = float(factor_str)
                    except Exception as ex:
                        print(f'[ERROR][QC_SAT] outliers_factor {factor_str} for filter_macropixel_var_{idx}  must be a float number. Exception: {ex}')
                        self.filter_macropixel_var[idx]['factor'] = None
                        check_qc = False

        return check_qc

    def check_window_size(self,window_size,n_rows,n_cols):
        check_w = True
        if window_size % 2 == 0:
            print(f'[ERROR][QC_SAT] Window size ({window_size}) should be an uneven integer')
            check_w = False
        if n_rows >= 1 and n_cols >= 1:
            n_rc_max = min(n_rows, n_cols)
            if window_size > n_rc_max:
                print(f'[ERROR][QC_SAT] Window size ({self.window_size}) should be greater or equal than the maximum extrac size {n_rc_max}')
                check_w = False
        return check_w

    def check_spectral_variable(self,spectral_variable,bands_variable):
        check_qc = True
        if bands_variable is None or not bands_variable in self.dataset.variables:
            print(f'[ERROR][QC_SAT] Band wavelengths variable  {bands_variable} is not available in the dataset. Choose among:')
            print(f'[ERROR][QC_SAT] {list(self.dataset.variables)}')
            return False
        if spectral_variable is None or not spectral_variable in self.dataset.variables:
            print(f'[ERROR][QC_SAT] Spectral variable {spectral_variable} is not available in the dataset. Choose among:')
            print(f'[ERROR][QC_SAT] {list(self.dataset.variables)}')
            return False

        if len(self.dataset.variables[bands_variable].shape) != 1:
            print(f'[ERROR][QC_SAT] Band wavelengths variable {bands_variable} should be a 1D array')
            check_qc = False
            n_bands = -1
        else:
            n_bands = self.dataset.variables[bands_variable].shape[0]

        if len(self.dataset.variables[spectral_variable].shape) != 4:
            print(f'[ERROR][QC_SAT] Spectral variable {self.spectral_variable} should have 4 dimensions: satellite_id,satellite_bands,rows,columns')
            check_qc = False
        else:
            n_bands_spectral = self.dataset.variables[spectral_variable].shape[1]
            if n_bands_spectral != n_bands:
                print(f'[ERROR][QC_SAT] The number of bands in the spectral variable ({n_bands_spectral}) should be equal to the number of bands in the band wavelengths variable ({n_bands}) ')
                check_qc = False

        return check_qc

    def check_non_spectral_variable(self,name_var):
        check_qc = True
        if name_var is None or not name_var in self.dataset.variables:
            print(f'[ERROR][QC_SAT] Non-spectral variable {name_var} is not available in the dataset. Choose among: ')
            print(f'[ERROR][QC_SAT] {list(self.dataset.variables)}')
            return False

        if len(self.dataset.variables[name_var].shape) != 3:
            print(f'[ERROR][QC_SAT] Spectral variable {name_var} should have 3 dimensions: satellite_id, rows, columns')
            check_qc = False

        return check_qc

    def check_flag_list(self,name_variable,flag_list):
        if name_variable is None:
            print(f'[ERROR][QC_SAT] name_variable is required for flagging filters, it could not be None')
            return False
        if not name_variable in self.dataset.variables:
            print(f'[ERROR][QC_SAT] {name_variable} is not available in the dataset')
            return False
        if flag_list is None:
            print(f'[ERROR][QC_SAT] {flag_list} is required for flagging filters, it could not be None')
            return False
        # ##flag list could be given as: flag_meanings (string with space separated flag) or flag_list (comma separated list)
        # flag_list_var = None
        # if 'flag_meanings' in self.dataset.variables[name_variable].ncattrs():
        #     flag_list_var = [x.strip() for x in self.dataset.variables[name_variable].flag_meanings.split(' ')]
        # elif 'flag_list' in self.dataset.variables[name_variable].ncattrs():
        #     flag_list_var = [x.strip() for x in self.dataset.variables[name_variable].flag_list.split(',')]
        # if flag_list_var is None:
        #     print(f'[ERROR][QC_SAT] Flag list in not available for variable {name_variable}, attribute flag_meanings or flag_list is required')
        #     return False
        flag_list_var, flag_values_var = ffs.get_info_from_flag_variable(self.dataset.variables[name_variable],key_error='QC_SAT')
        if flag_values_var is None or flag_list_var is None:
            return False

        check = set(flag_list).issubset(set(flag_list_var))
        if not check:
            print(f'[ERROR][QC_SAT] {flag_list} flags are not available in the variable {name_variable} flag list: {flag_list_var}')
        return check

    def set_basic_dimensions(self):
        ##method to be called only if check parameters is True
        self.n_bands = len(self.wl_list)
        self.n_mu = self.dataset.variables[self.spectral_variable].shape[0]
        self.n_rows = self.dataset.variables[self.spectral_variable].shape[2]
        self.n_cols = self.dataset.variables[self.spectral_variable].shape[3]

    def compute_validity(self):
        spectral_data = self.dataset.variables[self.spectral_variable][:,self.indices_wl,:,:]
        spectral_data = np.ma.masked_invalid(spectral_data)##make sure that invalid values, NaN and so on are masked
        mask_invalid = self.compute_invalid_masks_array(spectral_data)
        if self.verbose:
            print(f'[INFO][QC_SAT] Number of pixels filtered using {self.filter_invalid} invalid filter: {np.sum(mask_invalid)}')
        mask_flag = self.compute_flag_mask_array()

    def get_window_dimensions(self,w_size=None,n_rows=None,n_cols=None):
        if w_size is None:
            w_size = self.window_size
        if n_rows is None:
            n_rows = self.n_rows
        if n_cols is None:
            n_cols = self.n_cols
        central_r = int(np.floor(n_rows / 2))
        central_c = int(np.floor(n_cols / 2))
        r_s = central_r - int(np.floor(w_size / 2))  # starting row
        r_e = central_r + int(np.floor(w_size / 2)) + 1  # ending row
        c_s = central_c - int(np.floor(w_size / 2))  # starting col
        c_e = central_c + int(np.floor(w_size / 2)) + 1  # ending col
        return central_r, central_c, r_s, r_e, c_s, c_e

    def compute_invalid_masks_array(self,spectral_data):
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_window_dimensions()
        spectral_data_window = np.moveaxis(spectral_data[:,:,r_s:r_e, c_s:c_e],1,-1)
        invalid_mask_all_bands = np.where(spectral_data_window.mask,1,0)
        invalid_mask = np.sum(invalid_mask_all_bands,axis=-1)
        if self.filter_invalid=='any':
            invalid_mask[invalid_mask>=1]=1
        elif self.filter_invalid=='all':
            invalid_mask[invalid_mask<self.n_bands]=0
            invalid_mask[invalid_mask==self.n_bands]=1
        return invalid_mask


    def compute_flag_mask_array(self):
        flag_mask = np.zeros((self.n_mu, self.window_size, self.window_size), dtype=np.uint64)
        for idx in range(len(self.filter_flag)):
            flag_mask_here = self.compute_flag_mask_array_impl(self.filter_flag[idx])
            if self.verbose:
                print(f'[INFO][QC_SAT] filter_flag_{idx}: Variable: {self.filter_flag[idx]["name_var"]}. Flagged pixels: {np.sum(flag_mask_here)}')
            if flag_mask_here is not None:
                flag_mask = flag_mask + flag_mask_here
            #self.info_flag[flag_band]['nflagged'] = np.sum(flag_mask.reshape((self.nmu, self.window_size * self.window_size)), axis=1)
        flag_mask[flag_mask > 0] = 1
        return flag_mask

    def compute_flag_mask_array_impl(self, f_flag):
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_window_dimensions(w_size=f_flag['window_size'])
        flag_array = self.dataset.variables[f_flag['name_var']][:,r_s:r_e,c_s:c_e]
        if np.issubdtype(flag_array.dtype, np.floating): ##floating points are not allowed
            flag_array = flag_array.astype('uint64')

        fw = ffs.start_flag_work_from_variable(self.dataset.variables[f_flag['name_var']])
        if fw is None:
            print(f'[ERROR][QC_SAT] Flag work object for variable {f_flag["name_var"]} could not be started')
            return np.zeros(flag_array.shape, dtype=np.int8)
        if self.verbose:
            print(f'[INFO][QC_SAT] Starting FlagWork object with variable {f_flag["name_var"]}. Data type: {fw.dType}')

        ##invalid
        mask_array = None
        if f_flag['flag_list'] is not None:
            mask_array = fw.mask(flag_array,f_flag['flag_list'])
            mask_array[mask_array >= 1] = 1
        if f_flag['flag_list_valid'] is not None:
            f_list = flag['flag_list_valid']
            if '$' in f_list:
                mask_array_v = np.where(flag_array==0,0,1)
                if len(f_list)>1:
                    flag_list.remove('$')
                    m_array = fw.mask(flag_array,f_list)
                    mask_array_v[m_array>=1] = 0
            else:
                m_array = fw.mask(flag_array, f_list)
                mask_array_v = np.where(m_array>=1,0,1)
            if mask_array is None:
                mask_array = mask_array_v
            else:
                mask_array = np.where(np.logical_and(mask_array_v==0,ask_array==0,mask_array),0,1)
        return mask_array
                
            



        # land = None
        # central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        # satellite_flag = self.info_flag[flag_band]['variable']
        # if satellite_flag is None:
        #     flag_mask = np.zeros((self.nmu, self.window_size, self.window_size), dtype=np.uint64)
        #     return flag_mask, land
        # # flag_meanings_ string separated by spaces or list
        # flag_meanings = satellite_flag.flag_meanings
        # if isinstance(flag_meanings, list):
        #     flag_meanings = ' '.join(flag_meanings)
        #
        # satellite_flag_band = satellite_flag[:, r_s:r_e, c_s:c_e]
        # # float32 is not allowed
        # if str(satellite_flag.dtype) == 'float32':
        #     satellite_flag_band = satellite_flag_band.astype('uint64')
        #
        # # flag list, it could be a list or a comma separated string
        # flag_list_tobe_applied = self.info_flag[flag_band]['flag_list']
        # if isinstance(flag_list_tobe_applied, str):
        #     flag_list_tobe_applied = [x.strip() for x in flag_list_tobe_applied.split(',')]
        #
        # if self.info_flag[flag_band]['ac_processor'] == 'POLYMER':
        #     flagging = flag.Class_Flags_Polymer(satellite_flag.flag_masks, flag_meanings)
        #     flag_mask = flagging.MaskGeneral(satellite_flag_band)
        #     flag_mask[np.where(flag_mask != 0)] = 1
        # elif self.info_flag[flag_band]['ac_processor'] == 'IDEPIX':
        #     flagging = flag.Class_Flags_Idepix(satellite_flag.flag_masks, flag_meanings)
        #     flag_mask = flagging.Mask(satellite_flag_band, flag_list_tobe_applied)
        #     flag_mask[np.where(flag_mask != 0)] = 1
        # else:
        #     ##we must be sure that flag_mask must be uint64
        #     satellite_flag_band = satellite_flag_band.astype('uint64')
        #     flag_masks = satellite_flag.flag_masks.astype('uint64')
        #     flagging = flag.Class_Flags_OLCI(flag_masks, flag_meanings)
        #     flag_mask = flagging.Mask(satellite_flag_band, flag_list_tobe_applied)
        #     flag_mask[np.where(flag_mask != 0)] = 1

        return flag_mask


    def check_rrs_variability(self):
        if len(self.wl_ref)<self.nbands:
            valid_bands = np.array([1 if wl_here in self.wl_ref else 0 for wl_here in self.sat_bands])
            print(f'[INFO][QC_SAT] Working with a subset of {np.sum(valid_bands)} bands (Total: {self.nbands})')
            if np.sum(valid_bands)<len(self.wl_ref):
                wl_to_check = self.sat_bands[valid_bands==0]
                print(f'[WARNING] Not all the bands given in the band list are available. The following satellite bands are missing: {wl_to_check}')
                print(f'[WARNING] Expected band list: {self.wl_ref}')


            self.indices_valid_bands = np.where(valid_bands==1)[0]

        flag_mask,land = self.compute_flag_mask_array()
        print(f'[INFO][QC_SAT]->Number of flagged pixels: {np.ma.sum(flag_mask)}/{np.ma.count(flag_mask)}')
        print(f'[INFO][QC_SAT]->Number of land pixels:  {np.ma.sum(land)}/{np.ma.count(land)}')
        mask_invalid = self.compute_invalid_masks_array()
        print(f'[INFO][QC_SAT]->Number of invalid rrs: {np.ma.sum(mask_invalid)}/{np.ma.count(mask_invalid)}')
        mask_th = self.compute_th_masks_array()
        print(f'[INFO][QC_SAT]->Number of pixels masked using user-defined thresholds: {np.ma.sum(mask_th)}/{np.ma.count(mask_th)}')

        final_mask = flag_mask + mask_invalid + mask_th
        final_mask[final_mask>0]=1

        print(f'[INFO][QC_SAT]->Number of masked pixels in the final mask: {np.ma.sum(final_mask)}/{np.ma.count(final_mask)}')


        ntotal_by_mu = self.window_size*self.window_size
        nmasked_by_mu = np.ma.sum(np.ma.reshape(final_mask,(self.nmu,self.window_size*self.window_size)),axis=1)
        nvalid_by_mu = ntotal_by_mu-nmasked_by_mu


        min_valid_pixels = self.min_valid_pixels

        if self.use_Bailey_Werdell:
            nland_by_mu = np.ma.sum(np.ma.reshape(land,(self.nmu,self.window_size*self.window_size)),axis=1)
            ntotalw_by_mu = ntotal_by_mu-nland_by_mu
            min_valid_pixels = np.floor(0.50 * ntotalw_by_mu) + 1
            min_valid_pixels[min_valid_pixels<self.min_valid_pixels]=self.min_valid_pixels

        min_pixel_condition = nvalid_by_mu >= min_valid_pixels
        print(f'[INFO][QC_SAT] Number of match-ups filtered because the number of valid pixels is lower than the required one: {np.count_nonzero(min_pixel_condition==False)}')
        masks_rrs = self.get_masks_rrs(final_mask)


        macropixel_filter = self.do_check_macropixel(final_mask,masks_rrs)

        macropixel_condition = macropixel_filter==0

        all_conditions = np.logical_and(min_pixel_condition,macropixel_condition)
        print(f'[INFO][QC_SAT] Final number of match-ups passing the satellite quality control: {np.sum(all_conditions)} / {self.nmu}')

        outliers_str = 'with' if self.apply_outliers else 'without'
        if outliers_str in masks_rrs:
            mask_rrs = masks_rrs[outliers_str]
        else:
            print(f'[WARNING] Mask with outliers is not available. Using mask without ouliers')
            mask_rrs = masks_rrs['without']

        std_rrs = self.get_stat_rrs('std', mask_rrs)
        cv_rrs = self.get_stat_rrs('CV',mask_rrs)
        nvalues_rrs = self.get_stat_rrs('n_values',mask_rrs)

        return std_rrs,cv_rrs,nvalues_rrs




    def check_validity_deprecated(self):


        if len(self.wl_ref)<self.nbands:
            valid_bands = np.array([1 if wl_here in self.wl_ref else 0 for wl_here in self.sat_bands])
            print(f'[INFO][QC_SAT] Working with a subset of {np.sum(valid_bands)} bands (Total: {self.nbands})')
            if np.sum(valid_bands)<len(self.wl_ref):
                wl_to_check = self.sat_bands[valid_bands==0]
                print(f'[WARNING] Not all the bands given in the band list are available. The following satellite bands are missing: {wl_to_check}')
                print(f'[WARNING] Expected band list: {self.wl_ref}')


            self.indices_valid_bands = np.where(valid_bands==1)[0]

        self.satellite_rrs = np.ma.masked_invalid(self.satellite_rrs)##make sure that nan,-inf,inf are masked

        flag_mask,land = self.compute_flag_mask_array()


        print(f'[INFO][QC_SAT]->Number of flagged pixels: {np.ma.sum(flag_mask)}/{np.ma.count(flag_mask)}')
        print(f'[INFO][QC_SAT]->Number of land pixels:  {np.ma.sum(land)}/{np.ma.count(land)}')
        mask_invalid = self.compute_invalid_masks_array()
        print(f'[INFO][QC_SAT]->Number of invalid rrs: {np.ma.sum(mask_invalid)}/{np.ma.count(mask_invalid)}')
        mask_th = self.compute_th_masks_array()
        print(f'[INFO][QC_SAT]->Number of pixels masked using user-defined thresholds: {np.ma.sum(mask_th)}/{np.ma.count(mask_th)}')

        final_mask = flag_mask + mask_invalid + mask_th
        final_mask[final_mask>0]=1
        print('mu31 final mask-->', np.sum(mask_invalid[31, :]))

        print(f'[INFO][QC_SAT]->Number of masked pixels in the final mask: {np.ma.sum(final_mask)}/{np.ma.count(final_mask)}')


        ntotal_by_mu = self.window_size*self.window_size
        nmasked_by_mu = np.ma.sum(np.ma.reshape(final_mask,(self.nmu,self.window_size*self.window_size)),axis=1)
        nvalid_by_mu = ntotal_by_mu-nmasked_by_mu


        min_valid_pixels = self.min_valid_pixels

        if self.use_Bailey_Werdell:
            nland_by_mu = np.ma.sum(np.ma.reshape(land,(self.nmu,self.window_size*self.window_size)),axis=1)
            ntotalw_by_mu = ntotal_by_mu-nland_by_mu
            min_valid_pixels = np.floor(0.50 * ntotalw_by_mu) + 1
            min_valid_pixels[min_valid_pixels<self.min_valid_pixels]=self.min_valid_pixels

        min_pixel_condition = nvalid_by_mu >= min_valid_pixels
        print(f'[INFO][QC_SAT] Number of match-ups filtered because the number of valid pixels is lower than the required one: {np.count_nonzero(min_pixel_condition==False)}')
        masks_rrs = self.get_masks_rrs(final_mask)


        macropixel_filter = self.do_check_macropixel(final_mask,masks_rrs)

        macropixel_condition = macropixel_filter==0

        all_conditions = np.logical_and(min_pixel_condition,macropixel_condition)
        print(f'[INFO][QC_SAT] Final number of match-ups passing the satellite quality control: {np.sum(all_conditions)} / {self.nmu}')

        outliers_str = 'with' if self.apply_outliers else 'without'
        if outliers_str in masks_rrs:
            mask_rrs = masks_rrs[outliers_str]
        else:
            print(f'[WARNING] Mask with outliers is not available. Using mask without ouliers')
            mask_rrs = masks_rrs['without']
        reported_rrs = self.get_stat_rrs(self.stat_value,mask_rrs)



        reported_rrs_unc = None
        if self.satellite_rrs_unc is not None:
            self.satellite_rrs_unc = np.ma.masked_invalid(self.satellite_rrs_unc)##make sure nan,inf,-inf are masked
            central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
            if self.indices_valid_bands is not None:
                rrs_unc = self.satellite_rrs_unc[:, self.indices_valid_bands, r_s:r_e, c_s:c_e]
            else:
                rrs_unc = self.satellite_rrs_unc[:, :, r_s:r_e, c_s:c_e]
            rrs_unc[mask_rrs] = np.ma.masked
            reported_rrs_unc = self.get_stat_spectral_impl(self.stat_value,rrs_unc)



        self.qc_sat_results = {
            'min_pixel_condition':min_pixel_condition,
            'macropixel_condition':macropixel_condition,
            'all_conditions': all_conditions,
            'reported_rrs':reported_rrs,
            'reported_rrs_unc':reported_rrs_unc
        }

        # print('----------------------------------------')
        # print(reported_rrs.shape)
        # print(reported_rrs[31,:])
        # indices = np.where(np.isnan(reported_rrs))
        # print(len(indices[0]))
        # print('--------------------')








    def compute_th_masks_array(self):
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        mask_thershold = np.zeros((self.nmu,self.window_size, self.window_size), dtype=np.uint64)

        for idx in range(len(self.th_masks)):
            th_mask = self.th_masks[idx]
            if th_mask['index_sat'] >= 0:
                band_here = self.satellite_rrs[:,th_mask['index_sat'], r_s:r_e, c_s:c_e]
            else:
                var_here = self.ncdataset.variables[th_mask['band_name']]
                band_here = var_here[:, r_s:r_e, c_s:c_e]

            #mask_thershold_here = np.zeros(band_here.shape, dtype=np.uint64)
            n_masked = 0
            if th_mask['type_th'] == 'greater':
                mask_thershold[band_here > th_mask['value_th']] = mask_thershold[band_here > th_mask['value_th']]+1
                n_masked = np.count_nonzero(band_here > th_mask['value_th'])
            elif th_mask['type_th'] == 'lower':
                mask_thershold[band_here < th_mask['value_th']] = mask_thershold[band_here < th_mask['value_th']]+1
                n_masked = np.count_nonzero(band_here < th_mask['value_th'])
            th_mask['n_masked'] = n_masked
            self.th_masks[idx] = th_mask

        mask_thershold[mask_thershold>0]=1

        return mask_thershold

    def get_masks_rrs(self,final_mask):

        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        if self.indices_valid_bands is not None:
            rrs = self.satellite_rrs[:,self.indices_valid_bands,r_s:r_e, c_s:c_e]
        else:
            rrs = self.satellite_rrs[:, :, r_s:r_e, c_s:c_e]
        print(f'[INFO][QC_SAT] Number of valid rrs values before masking: {np.ma.count(rrs)}')
        nbands_used = rrs.shape[1]
        for iband in range(nbands_used):
            rrs_band = np.ma.squeeze(rrs[:, iband, :, :])
            rrs_band[final_mask == 1] = np.ma.masked
            rrs[:,iband,:,:] = rrs_band[:,:,:]
        nvalid_total  = np.ma.count(rrs)
        nvalid_by_band = nvalid_total/nbands_used
        print(f'[INFO][QC_SAT] Number of valid rrs values after masking (without ouliers): {nvalid_total} By band: {nvalid_by_band}')

        masks = {'without':rrs.mask.copy()}

        if self.apply_outliers:
            central_stat = self.outliers_info['central_stat']
            rrs_central = None
            rrs_dispersion = None
            rrs_min_th = None
            rrs_max_th = None
            factor = self.outliers_info['factor']

            if central_stat=='avg':
                rrs_central = np.ma.mean(rrs.reshape((self.nmu,nbands_used,self.window_size*self.window_size)),axis=2)
            elif central_stat=='median':
                rrs_central = np.ma.median(rrs.reshape((self.nmu, nbands_used, self.window_size * self.window_size)),axis=2)

            dispersion_stat = self.outliers_info['dispersion_stat']
            if dispersion_stat=='std':
                rrs_dispersion = np.ma.std(rrs.reshape((self.nmu,nbands_used,self.window_size*self.window_size)),axis=2)
            elif dispersion_stat=='iqr':
                if factor<0:
                    ql = 25
                    qh = 75
                else:
                    ql = factor
                    qh = 100 - factor
                rrs_min_th = np.percentile(rrs.reshape((self.nmu,nbands_used,self.window_size*self.window_size)),ql,axis=2)
                rrs_max_th = np.percentile(rrs.reshape((self.nmu, nbands_used, self.window_size * self.window_size)),qh, axis=2)



            if rrs_central is not None and rrs_dispersion is not None and rrs_min_th is None and rrs_max_th is None:
                rrs_min_th = rrs_central  - (factor * rrs_dispersion)
                rrs_max_th = rrs_central  + (factor * rrs_dispersion)

            if rrs_min_th is not None and rrs_max_th is not None:
                rrs_min_th = np.repeat(rrs_min_th,self.window_size*self.window_size).reshape((self.nmu,nbands_used,self.window_size,self.window_size))
                rrs_max_th = np.repeat(rrs_max_th, self.window_size * self.window_size).reshape((self.nmu, nbands_used, self.window_size, self.window_size))
                rrs[rrs < rrs_min_th] = np.ma.masked
                rrs[rrs > rrs_max_th] = np.ma.masked

                masks['with'] = rrs.mask.copy()
                print(f'[INFO][QC_SAT] Number of valid rrs values after masking (with ouliers): {np.ma.count(rrs)}')
            else:
                print(f'[WARNING] Outliers could not be applied')


        return masks

    def get_stat_rrs(self,type_stat,mask_rrs):
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        if self.indices_valid_bands is not None:
            rrs = self.satellite_rrs[:, self.indices_valid_bands, r_s:r_e, c_s:c_e]
        else:
            rrs = self.satellite_rrs[:, :, r_s:r_e, c_s:c_e]
        rrs[mask_rrs] = np.ma.masked

        return self.get_stat_spectral_impl(type_stat,rrs)

    def get_stat_spectral_impl(self,type_stat,array):
        nbands_used = array.shape[1]
        if type_stat=='n_values':
            result = np.ma.count(array.reshape((self.nmu, nbands_used, self.window_size * self.window_size)), axis=2)
        elif type_stat=='avg':
            result = np.ma.mean(array.reshape((self.nmu,nbands_used,self.window_size*self.window_size)),axis=2)
        elif type_stat=='std':
            result = np.ma.std(array.reshape((self.nmu, nbands_used, self.window_size * self.window_size)), axis=2)
        elif type_stat=='median':
            result = np.ma.median(array.reshape((self.nmu, nbands_used, self.window_size * self.window_size)), axis=2)
        elif type_stat=='min':
            result = np.ma.min(array.reshape((self.nmu, nbands_used, self.window_size * self.window_size)), axis=2)
        elif type_stat=='max':
            result = np.ma.max(array.reshape((self.nmu, nbands_used, self.window_size * self.window_size)), axis=2)
        elif type_stat=='CV':
            avg = np.ma.mean(array.reshape((self.nmu, nbands_used, self.window_size * self.window_size)), axis=2)
            std = np.ma.std(array.reshape((self.nmu, nbands_used, self.window_size * self.window_size)), axis=2)
            result = (std / np.abs(avg)) * 100
        else:
            print(f'[WARNING] {type_stat} is not implemented in the computation of window stats')
            result = None

        return result

    def get_stats_non_spectral(self,array,type_stat,mask):
        if len(array.shape)==1:
            if type_stat=='n_values':
                result = np.array(array.mask).astype(np.int8)
            else:
                result = array[:]
            return result
        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        array = array[:, r_s:r_e, c_s:c_e]
        array[mask==1] = np.ma.masked

        if type_stat=='n_values':
            result = np.ma.count(array.reshape((self.nmu, self.window_size * self.window_size)), axis=1)
        elif type_stat=='avg':
            result = np.ma.mean(array.reshape((self.nmu,self.window_size*self.window_size)),axis=1)
        elif type_stat=='std':
            result = np.ma.std(array.reshape((self.nmu, self.window_size * self.window_size)), axis=1)
        elif type_stat=='median':
            result = np.ma.median(array.reshape((self.nmu, self.window_size * self.window_size)), axis=1)
        elif type_stat=='min':
            result = np.ma.min(array.reshape((self.nmu, self.window_size * self.window_size)), axis=1)
        elif type_stat=='max':
            result = np.ma.max(array.reshape((self.nmu, self.window_size * self.window_size)), axis=1)
        elif type_stat=='CV':
            avg = np.ma.mean(array.reshape((self.nmu, self.window_size * self.window_size)), axis=1)
            std = np.ma.std(array.reshape((self.nmu, self.window_size * self.window_size)), axis=1)
            result = (std / np.abs(avg)) * 100
        else:
            print(f'[WARNING] {type_stat} is not implemented in the computation of window stats')
            result = None

        return result



    def do_check_macropixel(self,final_mask,masks_rrs):
        macropixel_filter = np.zeros(self.nmu)
        ##filters based on rrs data
        required_rrs_stats = {}
        for check_stat in self.check_statistics:
            outliers_str = 'with' if check_stat['with_outliers'] else 'without'
            type_stat = check_stat['type_stat']
            ref = f'{type_stat}_{outliers_str}'
            if not ref in required_rrs_stats:
                if outliers_str in masks_rrs:
                    mask_rrs = masks_rrs[outliers_str]
                else:
                    print(f'[WARNING] Mask with outliers is not available. Using mask without ouliers')
                    mask_rrs = masks_rrs['without']
                required_rrs_stats[ref] = self.get_stat_rrs(type_stat,mask_rrs)

        for check_stat in self.check_statistics:
            outliers_str = 'with' if check_stat['with_outliers'] else 'without'
            type_stat = check_stat['type_stat']
            ref = f'{type_stat}_{outliers_str}'
            if not ref in required_rrs_stats:
                continue
            index_sat = check_stat['index_sat']
            if self.indices_valid_bands is not None:##recalculate index_sat if not all the bands are used
                wl_stat = self.sat_bands[index_sat]
                index_sat = int(np.argmin(np.abs(self.wl_ref-wl_stat)))
            stat_array = required_rrs_stats[ref]
            stat_array = stat_array[:,index_sat]
            if check_stat['type_th'] == 'greater':##false (+1) if stat>th
                print(f'[INFO][QC_SAT] RRS macro-pixel filter: {np.count_nonzero(stat_array>check_stat['value_th'])} match-ups filtered because {type_stat} at {self.sat_bands[index_sat]} nm > {check_stat["value_th"]}')
                macropixel_filter[stat_array>check_stat['value_th']] = macropixel_filter[stat_array>check_stat['value_th']]+1
            if check_stat['type_th'] == 'lower':##false (+1) if stat<th
                print(f'[INFO][QC_SAT] RRS macro-pixel filter: {np.count_nonzero(stat_array < check_stat['value_th'])} match-ups filtered because {type_stat} at {self.sat_bands[index_sat]} nm < {check_stat["value_th"]}')
                macropixel_filter[stat_array<check_stat['value_th']] = macropixel_filter[stat_array<check_stat['value_th']]+1

        ##filters based on non-rrs data
        for check_stat in self.check_statistics_norrs:
            name_band = check_stat['name_band']
            type_stat = check_stat['type_stat']
            array = check_stat['variable'][:]
            stat_array = self.get_stats_non_spectral(array,type_stat,final_mask)
            if check_stat['type_th'] == 'greater':##false (+1) if stat>th
                print(f'[INFO][QC_SAT] {name_band} macro-pixel filter: {np.count_nonzero(stat_array>check_stat['value_th'])} match-ups filtered because {type_stat} > {check_stat["value_th"]}')
                macropixel_filter[stat_array>check_stat['value_th']] = macropixel_filter[stat_array>check_stat['value_th']]+1
            if check_stat['type_th'] == 'lower':##false (+1) if stat<th
                print(f'[INFO][QC_SAT] {name_band} macro-pixel filter: {np.count_nonzero(stat_array < check_stat['value_th'])} match-ups filtered because {type_stat} < {check_stat["value_th"]}')
                macropixel_filter[stat_array<check_stat['value_th']] = macropixel_filter[stat_array<check_stat['value_th']]+1

        macropixel_filter[macropixel_filter>0]=1
        print(f'[INFO][QC_SAT] Total number of match-ups filtered based on macropixel filters: {np.ma.sum(macropixel_filter)}')

        return macropixel_filter




    ##type: 1: multiplied 2: divided
    def update_pi_correct(self, wl_list_rhow, type):
        for wl in wl_list_rhow:
            idx = self.get_index_sat_from_wlvalue(wl)
            if idx >= 0 and type == 1:
                self.pi_multiplied[idx] = True
            if idx >= 0 and type == 2:
                self.pi_divided[idx] = True

    def set_apply_invalid_mask_wl(self,wl,bvalue):
        for sat_index in range(self.nbands):
            wlhere = self.sat_bands[sat_index]
            sat_index_str = str(sat_index)
            if wlhere==wl:
                self.invalid_mask[sat_index_str]['apply_mask']=bvalue

    def update_invalid_mask(self):
        if self.wl_ref is None:
            return
        self.invalid_mask = {}
        for sat_index in range(self.nbands):
            apply_mask = False
            wlhere = self.sat_bands[sat_index]
            for wl in self.wl_ref:
                diffwl = abs(wlhere - wl)
                if diffwl < 1.0:
                    apply_mask = True
            sat_index_str = str(sat_index)
            self.invalid_mask[sat_index_str] = {
                'wavelength': self.sat_bands[sat_index],
                'apply_mask': apply_mask,
                'n_masked': 0,
                'ref': f'rrs_{self.sat_bands[sat_index]:.0f}_invalid'
            }

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
            if self.apply_outliers:
                cvalue = self.statistics[sat_index_str]['without_outliers'][self.outliers_info['central_stat']]
                dvalue = self.statistics[sat_index_str]['without_outliers'][self.outliers_info['dispersion_stat']]
                min_th = cvalue - (dvalue * self.outliers_info['factor'])
                max_th = cvalue + (dvalue * self.outliers_info['factor'])
                mask_outliers = np.zeros(rrs_valid.shape)
                mask_outliers[rrs_valid > max_th] = 1
                mask_outliers[rrs_valid < min_th] = 1
                n_outliers = np.sum(mask_outliers)
                if n_outliers > 0:
                    rrs_valid = rrs_valid[mask_outliers == 0]
                stats = self.compute_statistics_impl(self.statistics[sat_index_str]['with_outliers'], rrs_valid)
            self.statistics[sat_index_str]['with_outliers'] = stats

            ##uncertainties
            if self.satellite_rrs_unc is not None:
                rrs_here_unc = self.satellite_rrs_unc[index_mu, sat_index, r_s:r_e, c_s:c_e]
                rrs_valid_unc = rrs_here_unc[self.flag_mask == 0]
                stats_unc = self.compute_statistics_impl(self.statistics_unc[sat_index_str]['without_outliers'],rrs_valid_unc)
                self.statistics_unc[sat_index_str]['without_outliers'] = stats_unc
                if self.apply_outliers:
                    cvalue = self.statistics_unc[sat_index_str]['without_outliers'][self.outliers_info['central_stat']]
                    dvalue = self.statistics_unc[sat_index_str]['without_outliers'][self.outliers_info['dispersion_stat']]
                    min_th = cvalue - (dvalue * self.outliers_info['factor'])
                    max_th = cvalue + (dvalue * self.outliers_info['factor'])
                    mask_outliers_unc = np.zeros(rrs_valid_unc.shape)
                    mask_outliers_unc[rrs_valid_unc > max_th] = 1
                    mask_outliers_unc[rrs_valid_unc < min_th] = 1
                    n_outliers_unc = np.sum(mask_outliers_unc)
                    if n_outliers_unc > 0:
                        rrs_valid_unc = rrs_valid_unc[mask_outliers_unc == 0]
                    stats_unc = self.compute_statistics_impl(self.statistics_unc[sat_index_str]['with_outliers'], rrs_valid_unc)
                self.statistics_unc[sat_index_str]['with_outliers'] = stats_unc

        for check_stat_here in self.check_statistics_norrs:
            name_band = check_stat_here['name_band']
            var_here = check_stat_here['variable']
            if len(var_here.shape)==3:
                var_here_array = var_here[index_mu, r_s:r_e, c_s:c_e]
                var_here_valid = var_here_array[~var_here_array.mask]
            elif len(var_here.shape)==1:
                var_here_valid = var_here[index_mu]

            if not name_band in self.statistics_norrs:
                self.statistics_norrs[name_band] = {
                    'n_values': 0,
                    'avg': 0,
                    'std': 0,
                    'median': 0,
                    'min': 0,
                    'max': 0,
                    'CV': 0
                }
            stats = self.compute_statistics_impl(self.statistics_norrs[name_band], var_here_valid)
            self.statistics_norrs[name_band] = stats

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
            index_sat = str(check_stat['index_sat'])

            if index_sat in self.statistics:
                outliers_str = 'with_outliers'
                if not check_stat['with_outliers']:
                    outliers_str = 'without_outliers'
                val_here = self.statistics[index_sat][outliers_str][check_stat['type_stat']]
                if check_stat['type_th'] == 'greater' and val_here > check_stat['value_th']:
                    CHECK = False
                if check_stat['type_th'] == 'lower' and val_here < check_stat['value_th']:
                    CHECK = False
        for check_stat in self.check_statistics_norrs:
            name_band = check_stat['name_band']
            if name_band in self.statistics_norrs:
                val_here = self.statistics_norrs[name_band][check_stat['type_stat']]
                if check_stat['type_th'] == 'greater' and val_here > check_stat['value_th']:
                    CHECK = False
                if check_stat['type_th'] == 'lower' and val_here < check_stat['value_th']:
                    CHECK = False

        return CHECK

    def get_match_up_values_v2(self,index_mu):
        # self.qc_sat_results = {
        #     'min_pixel_condition': min_pixel_condition,
        #     'macropixel_condition': macropixel_condition,
        #     'all_conditions': all_conditions,
        #     'reported_rrs': reported_rrs,
        #     'reported_rrs_unc': reported_rrs_unc
        # }
        cond_min_pixels = self.qc_sat_results['min_pixel_condition'][index_mu]
        cond_stats = self.qc_sat_results['macropixel_condition'][index_mu]
        valid_mu = self.qc_sat_results['all_conditions'][index_mu]
        values = self.qc_sat_results['reported_rrs'][index_mu]
        values_unc = self.qc_sat_results['reported_rrs_unc'][index_mu] if  self.qc_sat_results['reported_rrs_unc'] is not None else None



        return cond_min_pixels, cond_stats, valid_mu, values, values_unc

    def get_match_up_values(self, index_mu):
        self.prepare_new_match_up()

        cond_min_pixels = self.compute_masks_and_check_roi(index_mu)


        cond_stats = False
        valid_mu = False

        wl_orig = []

        if self.wl_ref is None:
            print('ATTENTION: NON DEBERIA ARRIVARE QUI, SAREBBE UN ERRORRE', self.nbands)
            indexes_bands = range(self.nbands)
            wl_orig = self.sat_bands
        else:
            indexes_bands = []
            for wl in self.wl_ref:
                index = self.get_index_sat_from_wlvalue(wl)
                if index == -1:
                    print(f'[WARNING] No valid satellite band for wl: {wl}')
                wl_orig.append(self.sat_bands[index])
                indexes_bands.append(index)

        values = [0] * len(indexes_bands)
        values_unc = [0] * len(indexes_bands)

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
                if self.satellite_rrs_unc is not None:
                    values_unc[idx] = self.statistics_unc[sat_index_str][outliers_str][self.stat_value]
                else:
                    values_unc[idx] = -999
            if self.apply_band_shifting:
                values = bsc_qaa.bsc_qaa(values, wl_orig, self.wl_ref)

        return cond_min_pixels, cond_stats, valid_mu, values, values_unc

    def compute_masks_and_check_roi(self, index_mu):

        land = self.compute_flag_masks(index_mu)

        nv = self.NTP - np.sum(self.flag_mask)

        # if index_mu==362:
        #     print('After flag mask: ',nv)

        self.compute_invalid_masks(index_mu)
        nv = self.NTP - np.sum(self.flag_mask)

        # if index_mu == 362:
        #     print('After Invalid: ', nv)

        self.compute_th_masks(index_mu)
        self.NVP = self.NTP - np.sum(self.flag_mask)
        self.NTPW = self.NTP - np.sum(land, axis=(0, 1))
        # if index_mu == 362:
        #     print('After th: ',self.NVP)
        #     print(index_mu,'After th: ', self.NVP)
        #     print(f'[INFO][QC_SAT] Index mu: {index_mu}')
        #     print(f'[INFO][QC_SAT] Number total of pixels: {self.NTP}')
        #     print(f'[INFO][QC_SAT] Water pixels: {self.NTPW}')
        #     print(f'[INFO][QC_SAT] Valid (no-flag) pixels: {self.NVP}')

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

            #CHAMBELAc
            # if flag_mask_here is not None:
            #     print('aplica la chambella')
            #     central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions_inner(3)
            #     flag_mask_here[r_s:r_e, c_s:c_e] = 1
                # invalid = None
                # if self.ncdataset.platform=='A':
                #     invalid = [24]
                # if self.ncdataset.platform=='B':
                #     invalid = [19,22,24]
                # if invalid is not None:
                #     if index_mu in invalid:
                #         flag_mask_here[:,:] = 1

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

    def compute_invalid_masks(self, index_mu):

        central_r, central_c, r_s, r_e, c_s, c_e = self.get_dimensions()
        mask_invalid = np.zeros((self.window_size, self.window_size), dtype=np.uint64)
        for sat_index in range(self.nbands):
            sat_index_str = str(sat_index)
            if self.invalid_mask[sat_index_str]['apply_mask']:
                rrshere = self.satellite_rrs[index_mu, sat_index, r_s:r_e, c_s:c_e]
                mask_invalid_here = np.zeros(rrshere.shape, dtype=np.uint64)
                mask_invalid_here[rrshere.mask] = 1



                # if np.sum(mask_invalid_here)>0:
                #     print(index_mu, '->Index with invalid values: ',sat_index)
                n_masked = np.sum(mask_invalid_here)
                # if index_mu==6:
                #     print(sat_index,n_masked )
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

    def add_threhold_mask_range(self, wl_min, wl_max, value_th, type_th):
        for index_sat in range(self.nbands):
            if wl_min <= self.sat_bands[index_sat] <= wl_max:
                self.add_theshold_mask(index_sat, -1, value_th, type_th)

    def add_threshold_mask_norrs(self, band_name, value_th, type_th):
        th_mask = {
            'index_sat': -1,
            'band_name': band_name,
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

    def add_bands_norrs_statistics(self, name_band, type_stat, value_th, type_th):
        if self.ncdataset is None:
            print('[WARNING] Statistics for bands no rrs could not be added as dataset was not defined')
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
        self.check_statistics_norrs.append(check_val)

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

        # flag_meanings_ string separated by spaces or list
        flag_meanings = satellite_flag.flag_meanings
        if isinstance(flag_meanings, list):
            flag_meanings = ' '.join(flag_meanings)

        if satellite_flag is None:
            flag_mask = np.zeros((self.window_size, self.window_size), dtype=np.uint64)
            return flag_mask, land

        satellite_flag_band = satellite_flag[index_mu, r_s:r_e, c_s:c_e]
        # float32 is not allowed
        if str(satellite_flag.dtype) == 'float32':
            satellite_flag_band = satellite_flag_band.astype('uint64')

        # if index_mu==3:
        #     print('vistazo a los datos:')
        #     print(satellite_flag_band.shape)
        #     print(satellite_flag_band)

        # flag list, it coulb be a comma separated string
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
            # if index_mu==3:
            #     print(flag_list_tobe_applied)
            flag_mask = flagging.Mask(satellite_flag_band, flag_list_tobe_applied)
            # if index_mu==3:
            #     print(flag_list_tobe_applied)
            #     print(flag_mask)
            #     tal = [17536]
            #     tal = np.array(tal,dtype=np.uint64)
            #     print(flag_meanings)
            #     flag_meanings_l = [x.strip() for x in flag_meanings.split(' ')]
            #     for l in flag_meanings_l:
            #         ll = [l]
            #         ftal = flagging.Mask(tal,ll)
            #         print(l,ftal)
            flag_mask[np.where(flag_mask != 0)] = 1

            # if self.info_flag[flag_band][
            #     'ac_processor'] == 'C2RCC' and not flag_band == 'satellite_WQSF':  # C2RCC FLAGS
            #     valuePE = np.uint64(2147483648)
            #     flag_mask = np.ones(satellite_flag_band.shape, dtype=np.uint64)
            #     flag_mask[satellite_flag_band == valuePE] = 0
            # else:
            #     flag_mask = flagging.Mask(satellite_flag_band, flag_list_tobe_applied)
            #     flag_mask[np.where(flag_mask != 0)] = 1
            # for fl in self.info_flag[flag_band]['flag_list']:
            #     fltal = flagging.Mask(satellite_flag_band, ([fl]))
            #     ntal = np.count_nonzero(fltal)
            #     if ntal>0:
            #         print('----> ',fl,':',ntal)

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



    def get_dimensions_inner(self, wsize):
        # Dimensions
        nrows = self.satellite_rrs.shape[2]
        ncols = self.satellite_rrs.shape[3]
        central_r = int(np.floor(nrows / 2))
        central_c = int(np.floor(ncols / 2))
        r_s = central_r - int(np.floor(wsize / 2))  # starting row
        r_e = central_r + int(np.floor(wsize / 2)) + 1  # ending row
        c_s = central_c - int(np.floor(wsize / 2))  # starting col
        c_e = central_c + int(np.floor(wsize / 2)) + 1  # ending col
        return central_r, central_c, r_s, r_e, c_s, c_e

    def get_flag_defaults(self, ac_processor):
        flag_list = None
        flag_land = None
        flag_inlandwaters = None

        if ac_processor == 'STANDARD':
            flag_list = 'LAND,COASTLINE,CLOUD,CLOUD_AMBIGUOUS,CLOUD_MARGIN,INVALID,COSMETIC,SATURATED,SUSPECT,HISOLZEN,HIGHGLINT,SNOW_ICE,AC_FAIL,WHITECAPS,RWNEG_O2,RWNEG_O3,RWNEG_O4,RWNEG_O5,RWNEG_O6,RWNEG_O7,RWNEG_O8'
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

    def set_idepix_as_flag(self, satellite_idepix_flag):
        flag_list = ['IDEPIX_LAND', 'IDEPIX_COASTLINE', 'IDEPIX_INVALID', 'IDEPIX_CLOUD', 'IDEPIX_CLOUD_BUFFER',
                     'IDEPIX_CLOUD_SHADOW', 'IDEPIX_SNOW_ICE']
        flag_land = 'IDEPIX_LAND'
        flag_inlandwater = None
        self.info_flag = {}
        self.info_flag[satellite_idepix_flag.name] = {
            'variable': satellite_idepix_flag,
            'flag_list': flag_list,
            'flag_land': flag_land,
            'flag_inlandwater': flag_inlandwater,
            'ac_processor': 'IDEPIX',
            'nflagged': 0,
            'flag_stats': None
        }

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


def get_filter_list(options_config, prefix, key_values = None):
    if key_values is None:
        from QC_OPTIONS import QC_OPTIONS
        qc_options_general = QC_OPTIONS(None,False)
        retrieve_options,required = qc_options_general.gmanager.get_retrieve_options('QC_SAT')
        key_values = retrieve_options[prefix]['key_values']

    default_dict = {key:key_values[key]['default'] for key in key_values}
    index = 0
    exist_filter = True
    filter_list = []
    while exist_filter:
        key_filter = f'{prefix}{index}'
        exist_filter =  key_filter in options_config
        if exist_filter:
            dict_here = default_dict.copy()
            dict_here.update(options_config[key_filter])
            invalid_keys = options_config[key_filter].keys() - default_dict.keys()
            filter_list.append(dict_here)
            if len(invalid_keys)>0:
                print(f'[WARNING] The following {len(invalid_keys)}  keys in configuration file for {key_filter} are not valid:')
                print(f'[WARNING] --> {list(invalid_keys)}')
                print(f'[WARNING] --> valid expected keys: {list(default_dict.keys())}')

        index = index+1
    return filter_list