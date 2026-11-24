import numpy as np


class QC_COMMON_TOOLS:
    def __init__(self, dataset, max_diff_wl):
        self.dataset = dataset
        self.max_diff_wl = max_diff_wl
        self.th_types = ['greater', 'gt', 'greater_or_equal', 'gte', 'lower', 'lt', 'lower_or_equal', 'lte']

    def check_wl(self,name_var_wl,wl_min_here,wl_max_here,key_qc = '',key_log='-'):
        check_qc = True
        wl_min_abs = np.min(self.dataset.variables[name_var_wl][:]) - self.max_diff_wl
        wl_max_abs = np.max(self.dataset.variables[name_var_wl][:]) + self.max_diff_wl
        if not isinstance(wl_min_here, float):
            print(f'[ERROR][{key_qc}] wl_min {wl_min_here} for {key_log}  should be a float value')
            check_qc = False
        if not isinstance(wl_max_here, float):
            print(f'[ERROR][{key_qc}] wl_max {wl_max_here} for {key_log} should be a float value')
            check_qc = False
        if isinstance(wl_min_here, float) and isinstance(wl_max_here, float):
            if wl_min_here > wl_max_here:
                print(f'[ERROR][{key_qc}] wl_max {wl_max_here} should be greater or equal to wl_min {wl_min_here} for{key_log}')
                check_qc = False
            else:
                if wl_min_here < wl_min_abs or wl_min_here > wl_max_abs:
                    print(f'[ERROR][{key_qc}] wl_min {wl_min_here} should be in the spectral range {wl_min_abs} - {wl_max_abs} for {key_log}')
                    check_qc = False
                if wl_max_here < wl_min_abs or wl_max_here > wl_max_abs:
                    print(f'[ERROR][{key_qc}] wl_max {wl_max_here} should be in the spectral range {wl_min_abs} - {wl_max_abs} for {key_log}')
                    check_qc = False
        return check_qc

    def check_thresholds(self, th_value, th_min, th_max, th_types, key_qc = '',key_filter ='-',check_is_angle=False,is_angle=False):
        check_qc = True

        if th_value is not None:
            if not isinstance(th_value, float):
                print(f'[ERROR][{key_qc}] Threshold th_value {th_value} for {key_filter}  should be a float')
                check_qc = False
            if th_min is not None and th_max is not None:
                print(f'[WARNING][{key_qc}] Thresholds values th_min and th_max are ignored for {key_filter}, only th_value is considered')
            if th_types is None:
                print(f'[ERROR][{key_qc}] Threshold filter type th_type is required for {key_filter}. It should one of {self.th_types}')
                check_qc = False
            if not th_types[0] in self.th_types:
                print(f'[ERROR][{key_qc}] Threshold filter type th_type  {th_types[0]} for {key_filter} should one of {self.th_types}')
                check_qc = False

        if th_value is None:
            if th_min is not None and th_max is None:
                print(f'[ERROR][{key_qc}] Threshold th_max is required for {key_filter} as th_min is given.')
                check_qc = False
            elif th_min is None and th_max is not None:
                print(f'[ERROR][{key_qc}] Threshold th_min is required for {key_filter} as th_max is given.')
                check_qc = False
            elif th_min is not None and th_max is not None:
                if not isinstance(th_min, float):
                    print(f'[ERROR][{key_qc}] Threshold th_min {th_min} for {key_filter}  should be a float')
                    check_qc = False
                if not isinstance(th_max, float):
                    print(f'[ERROR][{key_qc}] Threshold th_max {th_max} for {key_filter}  should be a float')
                    check_qc = False
                if th_types is None:  ##set default to be used with th_min and th_max
                    th_types = ['gte', 'lt']

                if len(th_types) == 2:
                    if not th_types[0] in self.th_types[0:4]:
                        print(f'[ERROR][{key_qc}] Threshold filter type th_type (first element) {th_types[0]} for {key_filter}  should one of {self.th_types[0:4]}')
                        check_qc = False
                    elif not th_types[1] in self.th_types[4:]:
                        print(f'[ERROR][{key_qc}] Threshold filter type th_type (second element) {th_types[1]} for {key_filter}  should one of {self.th_types[4:]}')
                        check_qc = False
                    else:
                        if not check_is_angle or (check_is_angle and not is_angle):##th_max should be >=th_min except for angle variables
                            if th_max<th_min:
                                print(f'[ERROR] th_min {th_min }should be lower or equal to th_max {th_max} for {key_filter}')
                                check_qc = False
                else:
                    print(f'[ERROR] [{key_qc}] Threshold filter type th_type using th_min and th_max for {key_filter} should include two elements: {self.th_types[0:4]} , {self.th_types[4:]}')
                    check_qc = False
        return check_qc,th_types

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

 #return zero-ones array with same shape as data with 1 values meeting the condition (if invert, zero values meeting the condition)
def get_mask_threshold_impl(data,th_val,th_type,invert=False):
    if th_type=='greater' or th_type=='gt':
        mask_th = np.ma.where(data>th_val,1,0)
    elif th_type=='greater_or_equal' or th_type=='gte':
        mask_th = np.ma.where(data >= th_val, 1, 0)
    elif th_type=='lower' or th_type=='lt':
        mask_th = np.ma.where(data < th_val, 1, 0)
    elif th_type=='lower_or_equal' or th_type=='lte':
        mask_th = np.ma.where(data <= th_val, 1, 0)
    else:
        mask_th = None
    if invert and mask_th is not None:
        mask_th = np.ma.where(mask_th==1,0,1)
    return mask_th

def compute_stats_from_array(array,s_stat,axis=None):
    if s_stat == 'avg':
        result = np.ma.mean(array, axis=axis)
    elif s_stat == 'median':
        result = np.ma.median(array, axis=axis)
    elif s_stat == 'min':
        result = np.ma.min(array, axis=axis)
    elif s_stat == 'max':
        result = np.ma.max(array, axis=axis)
    elif s_stat == 'std':
        result = np.ma.std(array, axis=axis)
    elif s_stat == 'iqr':
        result = np.percentile(array, 75, axis=axis) - np.percentile(array, 25, axis=axis)
    elif s_stat == 'CV':
        result = (np.std(array, axis=axis) / np.abs(np.mean(array, axis=axis))) * 100
    return result

def get_thresholds_from_dict(dict_here):
    th_val = dict_here['th_value']
    th_min = dict_here['th_min']
    th_max = dict_here['th_max']
    th_types = dict_here['th_type']
    return th_val, th_min, th_max, th_types

def get_mask_threshold_from_dict(array,dict_here):
    th_val,th_min,th_max,th_types = get_thresholds_from_dict(dict_here)
    return get_mask_threshold(array, th_val, th_min, th_max, th_types)

def get_mask_threshold(array,th_val,th_min,th_max,th_types):
    if th_val is not None:
        mask_th = get_mask_threshold_impl(array, th_val, th_types[0])
        mask_th = mask_th * 2
    else:
        mask_th = get_mask_threshold_impl(array, th_min, th_types[0]) + get_mask_threshold_impl(array, th_max,th_types[1])
        mask_th[mask_th==1]=0


    return mask_th

def get_array_with_outliers_masked_from_dict(array,axis_window,dict_here,verbose=False):
    #print(dict_here,'-->',verbose)
    central_stat = dict_here['outliers_central_stat']
    dispersion_stat = dict_here['outliers_dispersion_stat']
    factor = dict_here['outliers_factor']
    return get_array_with_outliers_masked(array,central_stat,dispersion_stat,factor,axis_window,verbose=verbose)

##axis window should be (n_rows,n_cols)
def get_array_with_outliers_masked(array,central_stat,dispersion_stat,factor,axis_window,verbose=False):

    n_rows = array.shape[axis_window[0]]
    n_cols = array.shape[axis_window[1]]
    n_window = n_rows*n_cols

    array_c = None
    if central_stat=='avg':
        array_c = np.ma.mean(array, axis=axis_window)
    elif central_stat=='median':
        array_c = np.ma.median(array, axis=axis_window)

    array_d = None
    if dispersion_stat=='std':
        array_d = np.ma.std(array, axis=axis_window)
    elif dispersion_stat=='iqr':
        array_d = np.percentile(array, 75, axis=axis) - np.percentile(array, 25, axis=axis)

    if array_c is None or array_d is None:
        print(f'[WARNING][QC_SAT] Outliers could not be masked, using the original array')
        return array

    upper_array = array_c + (array_d*factor)
    lower_array = array_c - (array_d*factor)

    ##getting upper_array and lower_array with the same shape as array (match-ups x n_rows x n_cols x n_bands),
    repeat_shape = upper_array.shape + (n_rows,n_cols)
    upper_array = np.repeat(upper_array,n_window).reshape(repeat_shape)
    upper_array = np.moveaxis(upper_array,(-2, -1), axis_window)
    lower_array = np.repeat(lower_array,n_window).reshape(repeat_shape)
    lower_array = np.moveaxis(lower_array, (-2, -1), axis_window)

    ##masking array
    n_mask_before = np.ma.count_masked(array)
    array[array>upper_array] = np.ma.masked
    array[array<lower_array] = np.ma.masked
    n_masked = np.ma.count_masked(array) - n_mask_before

    if verbose:
        print(f'[INFO][QC_SAT] Number of pixels masked after applying outliers: {n_masked}')


    return array

def get_valid_match_ups_array(qc_sat,qc_ins):
    qc_sat_valid_match_ups = None
    qc_ins_valid_match_ups = None
    if len(qc_sat.validity_info) > 0:
        qc_sat_valid_match_ups = qc_sat.validity_info['qc_sat_valid_match_ups']
    else:
        print(f'[ERROR] Satellited quality control was not successfully completed to produce a valid match ups dataset')
        return [None]*2
    if len(qc_ins.validity_info) > 0:
        qc_ins_valid_match_ups = qc_ins.validity_info['qc_ins_valid_match_ups']
    else:
        print(f'[ERROR] In situ quality control was not successfully completed to produce a valid match ups dataset')
        return [None]*2

    valid_match_ups_array = np.ma.sum(np.vstack([qc_sat_valid_match_ups,qc_ins_valid_match_ups]),axis=0)
    valid_match_ups_array = np.ma.where(valid_match_ups_array==2,1,0)
    n_valid_match_ups = np.ma.sum(valid_match_ups_array)

    # for idx in range(len(qc_sat_valid_match_ups)):
    #     print('Idx: ',qc_sat_valid_match_ups[idx],qc_ins_valid_match_ups[idx])
    return valid_match_ups_array, n_valid_match_ups


