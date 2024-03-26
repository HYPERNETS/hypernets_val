import math
import os
import configparser


class OptionsManager():

    def __init__(self, config_file, options):
        self.options = None
        if options is None:
            if config_file is not None and os.path.exists(config_file):
                self.options = configparser.ConfigParser()
                self.options.read(config_file)
        else:
            self.options = options

        # config_file = 'a'
        # self.options = configparser.ConfigParser()
        # self.options.read(config_file)

    def get_virtual_flag_list(self):
        if self.options is None:
            return None
        slist = self.options.sections()
        sfinal = []
        for flag in slist:
            if self.options.has_option(flag, 'type'):
                if self.options[flag]['type'] == 'virtual_flag':
                    sfinal.append(flag)
        if len(sfinal) == 0:
            sfinal = None
        return sfinal

    def get_section_list(self, exclude_flags):
        if self.options is None:
            return None
        slist = self.options.sections()
        if exclude_flags is not None:
            sfinal = []
            for flag in slist:
                if flag not in exclude_flags:
                    sfinal.append(flag)
        else:
            sfinal = slist
        if len(sfinal) == 0:
            return None
        return sfinal

    ##when type option is selected, then the rest of options are only applied if type_group=type
    def read_options_as_dict(self, section, poptions):
        options_dict = {}
        type = None
        use_pow2_flags = False

        if self.options.has_option(section, 'type'):
            type = self.read_option(section, 'type', 'type', poptions, -1, use_pow2_flags)
        if type is None and self.options.has_option(section, 'typevirtual'):  ##typevirtual overwrites type
            type = self.read_option(section, 'typevirtual', 'typevirtual', poptions, -1, use_pow2_flags)
        if self.options.has_option(section, 'use_pow2_flags'):
            use_pow2_flags = self.options[section]['use_pow2_flags']
        options_dict = self.assign_options(options_dict, 'type', type)
        options_dict = self.assign_options(options_dict, 'use_pow2_flags', use_pow2_flags)

        for opt in poptions:
            if opt == 'type' or opt == 'typevirtual' or opt == 'use_pow2_flags':
                continue

            if type is not None and 'type_group' in poptions[opt]:

                if poptions[opt]['type_group'] != type:
                    continue


            if not opt.find('_index') >= 0:
                if self.options.has_option(section, opt):
                    value = self.read_option(section, opt, opt, poptions, -1, use_pow2_flags)
                    options_dict = self.assign_options(options_dict, opt, value)
                else:
                    if 'default' in poptions[opt].keys():
                        options_dict[opt] = poptions[opt]['default']
            else:
                if opt.find('_indexm')>=0:
                    idx = 0
                    has_option = True
                    while has_option:
                        inner_idx = 0
                        has_inner = True
                        while has_inner:
                            opt_here = opt.replace('_indexm', f'_{idx}_{inner_idx}')
                            #print(opt_here,idx, inner_idx)
                            if self.options.has_option(section, opt_here):
                                value = self.read_option(section, opt, opt_here, poptions, idx, use_pow2_flags)

                                options_dict = self.assign_options(options_dict, opt_here, value)
                                inner_idx = inner_idx + 1
                            else:
                                has_inner = False
                            if inner_idx==0:
                                has_option = False
                        idx = idx + 1
                else:
                    idx = 0
                    has_option = True
                    while has_option:
                        opt_here = opt.replace('_index', f'_{idx}')
                        if self.options.has_option(section, opt_here):
                            value = self.read_option(section, opt, opt_here, poptions, idx, use_pow2_flags)
                            options_dict = self.assign_options(options_dict, opt_here, value)
                        else:
                            has_option = False
                        idx = idx + 1


        return options_dict

    def assign_options(self, options_dict, key, value):
        keys = key.split('.')
        if len(keys) == 1:
            options_dict[key] = value
        elif len(keys) == 2:
            if keys[0] in options_dict.keys():
                options_dict[keys[0]][keys[1]] = value
            else:
                options_dict[keys[0]] = {keys[1]: value}
        return options_dict

    def read_option(self, section, opt, opt_here, poptions, idx, use_pow2_flags):
        default = None
        if 'default' in poptions[opt].keys():
            default = poptions[opt]['default']
        value = self.get_value_param(section, opt_here, default, poptions[opt]['type_param'])

        if 'list_values' in poptions[opt]:
            if value not in poptions[opt]['list_values']:
                value = None

        if poptions[opt]['type_param'] == 'strlist':
            # use_pow2_vflags = poptions[opt]['user_pow2_vflags']
            value = self.get_strlist_as_dict(value, opt, idx, use_pow2_flags)


        return value

    def get_strlist_as_dict(self, values, opt, idx, use_pow2_flags):
        if opt == 'flag_spatial_index':
            return self.get_dict_flag_spatial_index(values, idx, use_pow2_flags)
        if opt.startswith('flag_ranges_index'):
            return self.get_dict_flag_ranges_index(values, idx, use_pow2_flags)
        return values

    def get_dict_flag_ranges_index(self, values, idx, use_pow2_flags):
        value_dict = {}
        if use_pow2_flags:
            fvalue = int(math.pow(2, idx))
        else:
            fvalue = int(idx + 1)
        if len(values) == 3:
            value_dict = {
                'flag_var': values[0],
                'is_default': True,
                'flag_name': values[1],
                'flag_value': fvalue,
                'flag_condition': 'and'
            }
        if len(values) >= 4:
            try:
                minV = float(values[2])
            except:
                minV = None
            try:
                maxV = float(values[3])
            except:
                maxV = None
            value_dict = {
                'is_default': False,
                'flag_var': values[0],
                'flag_name': values[1],
                'flag_value': fvalue,
                'min_range': minV,
                'max_range': maxV,
                'flag_condition': 'and'
            }
            if len(values)==5:
                value_dict['flag_condition']=values[4]
        return value_dict

    def get_dict_flag_spatial_index(self, values, idx, use_pow2_flags):
        value_dict = {}
        if use_pow2_flags:
            fvalue = int(math.pow(2, idx))
        else:
            fvalue = int(idx + 1)

        if len(values) == 2:
            value_dict = {
                'is_default': True,
                'flag_name': values[1],
                'flag_value': fvalue
            }
        if len(values) == 5:
            value_dict = {
                'is_default': False,
                'lat_min': float(values[0]),
                'lat_max': float(values[1]),
                'lon_min': float(values[2]),
                'lon_max': float(values[3]),
                'flag_name': values[4],
                'flag_value': fvalue
            }
        return value_dict

    def get_value(self, section, key):
        value = None
        if self.options.has_option(section, key):
            value = self.options[section][key]
        return value

    def get_value_param(self, section, key, default, type):
        value = self.get_value(section, key)

        if value is None:
            return default
        if type == 'str':
            return value
        if type == 'file':
            if not os.path.exists(value.strip()):
                return default
            else:
                return value.strip()
        if type == 'directory':
            directory = value.strip()
            if not os.path.isdir(directory):
                try:
                    os.mkdir(directory)
                    return directory
                except:
                    return default
            else:
                return directory
        if type == 'int':
            return int(value)
        if type == 'float':
            return float(value)
        if type == 'boolean':
            if value == '1' or value.upper() == 'TRUE':
                return True
            elif value == '0' or value.upper() == 'FALSE':
                return False
            else:
                return True
        if type == 'rrslist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                vals = vals.strip().replace('.', '_')
                list.append(f'RRS{vals}')
            return list
        if type == 'strlist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                list.append(vals.strip())
            return list
        if type == 'floatlist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                vals = vals.strip()
                list.append(float(vals))
            return list
        if type == 'intlist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                vals = vals.strip()
                list.append(int(vals))
            return list
