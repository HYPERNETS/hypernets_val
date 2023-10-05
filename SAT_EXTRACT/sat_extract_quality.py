from netCDF4 import Dataset
import os
import numpy as np
import COMMON.Class_Flags_OLCI as flag


class QC_EXTRACT:

    def __init__(self, fconfig):
        self.final_mask = None
        self.info_flag = {}
        self.th_masks = []
        self.ncdataset = None
        self.options = None
        if fconfig is not None:
            self.start_quality_options_from_fconfig(fconfig)


    def start_quality_options_from_fconfig(self, fconfig):
        import configparser
        self.options = configparser.ConfigParser()
        self.options.read(fconfig)
        section = 'QC_EXTRACT'
        options_qcsat = {
            'info_flag_X': {'valid': 0, 'value': None, 'type': 'dict',
                            'keys': {'name': 'str', 'flag_list': 'str', 'ac_processor': 'str'}
                            },
            'band_th_X': {'valid': 0, 'value': None, 'type': 'dict',
                          'keys': {'band_name': 'str', 'th_value': 'float', 'th_type': 'str'}
                          }
        }
        for option in options_qcsat:
            if option.endswith('_X'):
                val = self.get_value_param_list(section, option, options_qcsat[option])
            if val is None:
                options_qcsat[option]['valid'] = -1  # invalid
            elif val != 'N/A':
                options_qcsat[option]['valid'] = 1  # valid
                options_qcsat[option]['value'] = val

        ##CHECKING INVALID OPTIONS
        invalid_options = False
        for option in options_qcsat:
            if options_qcsat[option]['valid'] == -1:
                print(f'[ERROR] Option {option} is not valid. Please review configuration file')
                invalid_options = True
            if invalid_options:
                return

        ##setting options
        for option in options_qcsat:
            if options_qcsat[option]['valid'] == 0:
                continue
            if option == 'info_flag_X':
                list_vals = options_qcsat[option]['value']
                if len(list_vals) > 0:
                    info_flag = {}
                for val in list_vals:
                    self.add_info_flag(val['name'], val['flag_list'], val['ac_processor'])
                for name_var in info_flag:
                    print(f'[INFO] Flag band: {name_var}')
                    flist = info_flag[name_var]['flag_list']
                    print(f'[INFO] Flag list: {flist}')
            elif option == 'band_th_X':
                list_vals = options_qcsat[option]['value']
                if len(list_vals) > 0:
                    self.th_masks = []
                for val in list_vals:
                    self.add_threshold_mask_norrs(val['band_name'], val['th_value'], val['th_type'])

    def compute_array_valid(self, fextract):
        self.reset()
        self.ncdataset = Dataset(fextract)
        self.compute_flag_masks()
        self.compute_th_masks()
        array = np.ones(self.final_mask.shape, dtype=np.int16)
        array[self.final_mask > 0] = 0  ##invalid values set to zero, valid values to ones
        self.ncdataset.close()
        return array

    def reset(self):
        self.ncdataset = None
        self.final_mask = None

    def compute_flag_masks(self):
        flag_mask = None  # np.zeros((self.window_size, self.window_size), dtype=np.uint64)

        for flag_band in self.info_flag.keys():
            flag_mask_here = self.compute_flag_mask_impl(flag_band)
            if flag_mask_here is not None:
                if flag_mask is None:
                    flag_mask = flag_mask_here
                else:
                    flag_mask = np.add(flag_mask, flag_mask_here)
                self.info_flag[flag_band]['nflagged'] = np.sum(flag_mask_here)

        if self.final_mask is None:
            self.final_mask = flag_mask

        self.final_mask[flag_mask > 0] = 1

    def compute_flag_mask_impl(self, flag_band):
        if self.ncdataset is None:
            return None
        if flag_band not in self.ncdataset.variables:
            return None
        flag_var = self.ncdataset.variables[flag_band]

        flag_meanings = flag_var.flag_meanings
        if isinstance(flag_meanings, list):
            flag_meanings = ' '.join(flag_meanings)

        flag_array = flag_var[0, :, :]
        # float32 is not allowed
        if str(flag_var.dtype) == 'float32':
            flag_array = flag_array.astype('uint64')

        # flag list, it coulb be a comma separated string
        flag_list_tobe_applied = self.info_flag[flag_band]['flag_list']
        if isinstance(flag_list_tobe_applied, str):
            flag_list_tobe_applied = [x.strip() for x in flag_list_tobe_applied.split(',')]

        if self.info_flag[flag_band]['ac_processor'] == 'POLYMER':
            flagging = flag.Class_Flags_Polymer(flag_var.flag_masks, flag_meanings)
            flag_mask = flagging.MaskGeneral(flag_array)
            flag_mask[np.where(flag_mask != 0)] = 1
        elif self.info_flag[flag_band]['ac_processor'] == 'IDEPIX':
            flagging = flag.Class_Flags_Idepix(flag_var.flag_masks, flag_meanings)
            flag_mask = flagging.Mask(flag_array, flag_list_tobe_applied)
            flag_mask[np.where(flag_mask != 0)] = 1
        else:
            ##we must be sure that flag_mask must be uint64
            flag_masks = flag_var.flag_masks.astype('uint64')
            flagging = flag.Class_Flags_OLCI(flag_masks, flag_meanings)
            flag_mask = flagging.Mask(flag_array, flag_list_tobe_applied)
            flag_mask[np.where(flag_mask != 0)] = 1

        return flag_mask

    def compute_th_masks(self):
        if self.ncdataset is None:
            return
        mask_thershold = None
        for idx in range(len(self.th_masks)):
            th_mask = self.th_masks[idx]
            if th_mask['index_sat'] >= 0:
                var_here = self.ncdataset.variables[th_mask['satellite_Rrs']]
                band_here = var_here[0, th_mask['index_sat'], :, :]
            else:
                var_here = self.ncdataset.variables[th_mask['band_name']]
                band_here = var_here[0, :, :]

            mask_thershold_here = np.zeros(band_here.shape, dtype=np.uint64)
            if th_mask['type_th'] == 'greater':
                mask_thershold_here[band_here > th_mask['value_th']] = 1
            elif th_mask['type_th'] == 'lower':
                mask_thershold_here[band_here < th_mask['value_th']] = 1
            n_masked = np.sum(mask_thershold_here)
            th_mask['n_masked'] = n_masked
            self.th_masks[idx] = th_mask
            if mask_thershold is None:
                mask_thershold = mask_thershold_here
            else:
                mask_thershold = mask_thershold + mask_thershold_here

        if self.final_mask is None:
            self.final_mask = mask_thershold

        self.final_mask[mask_thershold > 0] = 1

    def add_threshold_mask_norrs(self, band_name, value_th, type_th):
        th_mask = {
            'index_sat': -1,
            'band_name': band_name,
            'value_th': value_th,
            'type_th': type_th,
            'n_masked': 0
        }
        self.th_masks.append(th_mask)

    def add_info_flag(self, name, flag_list, ac_processor):
        self.info_flag[name] = {
            'flag_list': flag_list,
            'ac_processor': ac_processor,
            'nflagged': 0,
            'flag_stats': None
        }

    def get_value(self, section, key):
        value = None
        if self.options.has_option(section, key):
            value = self.options[section][key].strip()
        return value

    def get_value_param_dict(self, section, key, default, keys_dict):
        value = {}
        has_data = False
        for kd in keys_dict:
            key_here = f'{key}.{kd}'
            val_here = self.get_value_param(section, key_here, 'N/A', keys_dict[kd])
            if val_here is None:
                return None
            # if val_here == 'N/A':
            #     return default
            if val_here != 'N/A':
                has_data = True
            value[kd] = val_here

        if not has_data:
            return default

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
        if type == 'file':
            if os.path.exists(value):
                return value
            else:
                return default
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
