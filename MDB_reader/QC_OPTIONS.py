class QC_OPTIONS:

    def __init__(self, options):
        # config_file = ''
        # import configparser
        # options = configparser.ConfigParser()
        # options.read(config_file)
        self.options = options

    def get_wllist(self):
        section = 'QC_SAT'
        wllist = self.get_value_param(section, 'wllist', 'N/A', 'floatlist')
        if wllist != 'N/A':
            return wllist
        default = self.get_value_param(section, 'default', 'N/A', 'str')
        wllist = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.8, 681.3, 708.8, 753.8]
        if default == 'EUMETSAT':
            wllist = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.8, 681.3, 708.8, 753.8]

        return wllist

    def get_qcsat(self, qc_sat):
        section = 'QC_SAT'
        options_qcsat = {
            'wllist': {'valid': 0, 'value': None, 'type': 'floatlist'},
            'window_size': {'valid': 0, 'value': None, 'type': 'int'},
            'min_valid_pixels': {'valid': 0, 'value': None, 'type': 'int'},
            'use_Bailey_Werdell': {'valid': 0, 'value': None, 'type': 'boolean'},
            'stat_value': {'valid': 0, 'value': None, 'type': 'str'},
            'apply_outliers': {'valid': 0, 'value': None, 'type': 'boolean'},
            'outliers_info': {'valid': 0, 'value': None, 'type': 'dict',
                              'keys': {'central_stat': 'str', 'dispersion_stat': 'str', 'factor': 'float'}
                              },
            'info_flag_X': {'valid': 0, 'value': None, 'type': 'dict',
                            'keys': {'name': 'str', 'flag_list': 'str', 'flag_land': 'str', 'flag_inlandwater': 'str'}
                            },
            'rrs_th_X': {'valid': 0, 'value': None, 'type': 'dict',
                         'keys': {'wl_min': 'float', 'wl_max': 'float', 'th_value': 'float', 'th_type': 'str'}
                         },
            'macropixel_filter_rrs_X': {'valid': 0, 'value': None, 'type': 'dict',
                                        'keys': {'wl': 'float', 'stat': 'str', 'withoutliers': 'boolean',
                                                 'th_value': 'float', 'th_type': 'str'}
                                        },
            'macropixel_filter_band_X': {'valid': 0, 'value': None, 'type': 'dict',
                                         'keys': {'band': 'str', 'stat': 'str', 'th_value': 'float', 'th_type': 'str'}
                                         }

        }

        # defaults options---------------------------------------------------------------------------------------------
        default = self.get_value_param(section, 'default', 'N/A', 'str')
        if default == 'EUMETSAT':
            wllist = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.8, 681.3, 708.8, 753.8]
            qc_sat.wl_ref = wllist
            window_size = self.get_value_param(section, 'window_size', 3, 'int')
            if window_size is None:
                window_size = 3
            qc_sat.set_eumetsat_defaults(window_size)
        ##END OF DEFAULT OPTIONS----------------------------------------------------------------------------------------

        ##---REAMINING OPTIONS----------------
        for option in options_qcsat:
            if option.endswith('_X'):
                val = self.get_value_param_list(section, option, options_qcsat[option])
            else:
                if options_qcsat[option]['type'] == 'dict':
                    val = self.get_value_param_dict(section, option, 'N/A', options_qcsat[option]['keys'])
                else:
                    val = self.get_value_param(section, option, 'N/A', options_qcsat[option]['type'])
            if val is None:
                options_qcsat[option]['valid'] = -1  # invalid
            elif val != 'N/A':
                options_qcsat[option]['valid'] = 1  # valid
                options_qcsat[option]['value'] = val
                # print(option, '----------------')
                # print(val)

        ##CHECKING INVALID OPTIONS
        invalid_options = False
        for option in options_qcsat:
            if options_qcsat[option]['valid'] == -1:
                print(f'[ERROR] Option {option} is not valid. Please review configuration file')
                invalid_options = True
        if invalid_options:
            return None

        ##setting options
        for option in options_qcsat:
            if options_qcsat[option]['valid'] == 0:
                continue
            if option == 'wllist':
                qc_sat.wl_ref = options_qcsat[option]['value']
            elif option == 'window_size':
                qc_sat.window_size = options_qcsat[option]['value']
            elif option == 'min_valid_pixels':
                qc_sat.min_valid_pixels = options_qcsat[option]['value']
            elif option == 'use_Bailey_Werdell':
                qc_sat.use_Bailey_Werdell = options_qcsat[option]['value']
            elif option == 'stat_value':
                qc_sat.stat_value = options_qcsat[option]['value']
            elif option == 'apply_outliers':
                qc_sat.apply_outliers = options_qcsat[option]['value']
            elif option == 'outliers_info':
                qc_sat.outliers_info = options_qcsat[option]['value']
            elif option == 'info_flag_X':
                list_vals = options_qcsat[option]['value']
                for val in list_vals:
                    qc_sat.info_flag[val['name']] = {
                        'variable': None,
                        'flag_list': val['flag_list'],
                        'flag_land': val['flag_land'],
                        'flag_inlandwater': val['flag_inlandwater'],
                        'ac_processor': qc_sat.ac_processor,
                        'nflagged': 0,
                        'flag_stats': None
                    }
            elif option == 'rrs_th_X':
                list_vals = options_qcsat[option]['value']
                for val in list_vals:
                    qc_sat.add_threhold_mask_range(val['wl_min'], val['wl_max'], val['th_value'], val['th_type'])
            elif option == 'macropixel_filter_rrs_X':
                list_vals = options_qcsat[option]['value']
                for val in list_vals:
                    qc_sat.add_band_statistics(-1, val['wl'], val['stat'], val['withoutliers'], val['th_value'],
                                               val['th_type'])
            elif option == 'macropixel_filter_band_X':
                list_vals = options_qcsat[option]['value']

        return qc_sat

    def get_qc_insitu(self, qc_insitu):
        section = 'QC_INS'
        ##wl list
        wllist = self.get_wllist()
        qc_insitu.set_wllist_using_wlref(wllist)
        # other options
        options_qc_insitu = {
            'time_diff_max': {'valid': 0, 'value': None, 'type': 'float'},
            'wl_diff_max': {'valid': 0, 'value': None, 'type': 'int'},
            'apply_band_shift': {'valid': 0, 'value': None, 'type': 'boolean'},
            'filter_th_X': {'valid': 0, 'value': None, 'type': 'dict',
                            'keys': {'wlmin': 'float', 'wlmax': 'float', 'thmin': 'float', 'thmax': 'float'}
                            }
        }

        ##---REAMINING OPTIONS----------------
        for option in options_qc_insitu:
            if option.endswith('_X'):
                val = self.get_value_param_list(section, option, options_qc_insitu[option])
            else:
                if options_qc_insitu[option]['type'] == 'dict':
                    val = self.get_value_param_dict(section, option, 'N/A', options_qc_insitu[option]['keys'])
                else:
                    val = self.get_value_param(section, option, 'N/A', options_qc_insitu[option]['type'])
            if val is None:
                options_qc_insitu[option]['valid'] = -1  # invalid
            elif val != 'N/A':
                options_qc_insitu[option]['valid'] = 1  # valid
                options_qc_insitu[option]['value'] = val
                # print(option, '----------------')
                # print(val)

        ##CHECKING INVALID OPTIONS
        invalid_options = False
        for option in options_qc_insitu:
            if options_qc_insitu[option]['valid'] == -1:
                print(f'[ERROR] Option {option} is not valid. Please review configuration file')
                invalid_options = True
        if invalid_options:
            return None

        for option in options_qc_insitu:
            if options_qc_insitu[option]['valid'] == 0:
                continue
            if option == 'time_diff_max':
                qc_insitu.time_max = options_qc_insitu[option]['value'] * 60  ##in seconds
            elif option == 'wl_diff_max':
                qc_insitu.wl_diff_max = options_qc_insitu[option]['value']
            elif option == 'apply_band_shift':
                qc_insitu.apply_band_shift = options_qc_insitu[option]['value']
            elif option == 'filter_th_X':
                list_vals = options_qc_insitu[option]['value']
                for val in list_vals:
                    th_min = val['thmin']
                    th_max = val['thmax']
                    if th_min == -999:
                        th_min = None
                    if th_max == -999:
                        th_max = None
                    qc_insitu.set_thershold(th_min, th_max, val['wlmin'], val['wlmax'])

        return qc_insitu

    def get_value(self, section, key):
        value = None
        if self.options.has_option(section, key):
            value = self.options[section][key].strip()
        return value

    def get_value_param_dict(self, section, key, default, keys_dict):
        value = {}
        for kd in keys_dict:
            key_here = f'{key}.{kd}'
            val_here = self.get_value_param(section, key_here, 'N/A', keys_dict[kd])
            if val_here is None:
                return None
            if val_here == 'N/A':
                return default
            value[kd] = val_here
        return value

    def get_value_param_list(self, section, key, options):
        option_base = key[:-1]
        values = []
        index = 0
        while index >= 0:
            option_here = f'{option_base}{index}'
            if options['type'] == 'dict':
                val = self.get_value_param_dict(section, option_here, 'N/A', options['keys'])
            else:
                val = self.get_value_param(section, option_here, 'N/A', options['type'])
            if val is None:
                return None
            else:
                if val == 'N/A':
                    break
                else:
                    values.append(val)
                    index = index + 1
        return values

    def get_value_param(self, section, key, default, type):
        value = self.get_value(section, key)
        if value is None:
            return default
        if type == 'str':
            return value
        if type == 'int':
            try:
                return int(value)
            except:
                return None
        if type == 'float':
            try:
                return float(value)
            except:
                return None
        if type == 'boolean':
            if value == '1' or value.upper() == 'TRUE':
                return True
            elif value == '0' or value.upper() == 'FALSE':
                return False
            else:
                return None
        if type == 'rrslist':
            list_str = value.split(',')
            if len(list) == 0:
                return None
            list = []
            for vals in list_str:
                vals = vals.strip().replace('.', '_')
                list.append(f'RRS{vals}')
            return list
        if type == 'strlist':
            list_str = value.split(',')
            if len(list_str) == 0:
                return None
            list = []
            for vals in list_str:
                list.append(vals.strip())
            return list
        if type == 'floatlist':
            list_str = value.split(',')
            if len(list_str) == 0:
                return None
            list = []
            for vals in list_str:
                vals = vals.strip()
                try:
                    list.append(float(vals))
                except:
                    return None
            return list