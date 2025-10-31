import math,os,re,configparser



class OptionsManager():

    def __init__(self, config_file, options):
        self.options = None
        if options is None:
            if not os.path.exists(config_file):
                print(f'[ERROR] {config_file} does not exist')
            elif not os.path.isfile(config_file):
                print(f'[ERROR] {config_file} is not a valid file')
            else:
                try:
                    self.options = configparser.ConfigParser()
                    self.options.read(config_file)

                except Exception as ex:
                    print(f'[ERROR] Error parsing configuration file {config_file}: {ex}')
        else:
            self.options = options


    def add_section(self,new_section):
        if self.options is not None:
            self.options.add_section(new_section)

    def remove_option(self,section,option):
        if self.options is not None:
            if self.options.has_option(section,option):
                self.options.remove_option(section, option)

    def add_value(self,section,option,value):
        if self.options is not None:
            if not self.options.has_section(section):
                self.options.add_section(section)
            self.options.set(section,option,value)

    def save_copy_as_file(self,file_config):
        if self.options is not None:
            with open(file_config, 'w') as configw:
                self.options.write(configw)



    def is_valid(self):
        return False if self.options is None else True

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

    def get_section_list(self, exclude_sections):
        if self.options is None:
            return None
        slist = self.options.sections()
        if exclude_sections is not None:
            sfinal = []
            for section in slist:
                if section not in exclude_sections:
                    sfinal.append(section)
        else:
            sfinal = slist
        if len(sfinal) == 0:
            return None
        return sfinal

    def get_options_list(self,section):
        if self.options is None:
            return None
        if not self.options.has_section(section):
            print(f'[WARNING] Section {section} is not available')
            return None
        options = self.options.options(section)
        return options

    def get_retrieve_options(self,section):
        slist = self.get_options_list(section)
        if slist is None:
            print(f'[ERROR] Retrieve options were not found for section {section}')
            return [None] * 2
        soptions = {}
        required_list = None

        for op in slist:
            if op == 'required':
                required_list = self.get_value_param(section, op, None, 'strlist')
            else:
                list = self.get_value_param(section, op, None, 'strlist')

                type_param = list[0]
                default = None if list[1].upper() == 'NONE' else list[1]
                if default is not None:
                    default = self.get_value_param_impl(default,type_param,None)


                soptions[op] = {
                    'type_param': type_param,
                    'default': default
                }
                if len(list) == 3:
                    #soptions[op]['list_values'] = [s.strip() for s in list[2].split(',')]
                    if type_param.startswith('str'):##str or strlist:
                        soptions[op]['list_values'] = self.get_value_param_impl(list[2],'strlist',None)
                    elif type_param.startswith('int'):##int or intlist:
                        soptions[op]['list_values'] = self.get_value_param_impl(list[2],'intlist',None)
                    elif type_param.startswith('float'):##float or floatlist:
                        soptions[op]['list_values'] = self.get_value_param_impl(list[2],'floatlist',None)

        return soptions, required_list

    ##when type option is selected, then the rest of options are only applied if type_group=type
    ##it includes also overall options for virtual_flags
    def read_options_as_dict(self, section, poptions):
        options_dict = {}
        type = None
        use_pow2_flags = False

        if self.options.has_option(section, 'type'):
            type = self.read_option(section, 'type', 'type', poptions, -1, use_pow2_flags)

        if (type == 'virtual_flag' or type is None) and self.options.has_option(section,
                                                                                'typevirtual'):  ##typevirtual overwrites type
            type = self.read_option(section, 'typevirtual', 'typevirtual', poptions, -1, use_pow2_flags)
        if type is None:
            print(
                f'[ERROR] type (or typevirtual) are not define. Potential types: {poptions["typevirtual"]["list_values"]}')
            return None

        use_pow2_flags = self.get_value_param(section,'use_pow2_flags',False,'boolean')


        options_dict = self.assign_options(options_dict, 'type', type)
        options_dict = self.assign_options(options_dict, 'use_pow2_flags', use_pow2_flags)

        for opt in poptions:
            if opt == 'type' or opt == 'typevirtual' or opt == 'use_pow2_flags':
                continue

            if type is not None and 'type_group' in poptions[opt]:
                if type not in poptions[opt]['type_group']:
                    continue

            if not opt.find('_index') >= 0:
                if self.options.has_option(section, opt):
                    value = self.read_option(section, opt, opt, poptions, -1, use_pow2_flags)
                    options_dict = self.assign_options(options_dict, opt, value)
                else:
                    if 'default' in poptions[opt].keys():
                        options_dict[opt] = poptions[opt]['default']
            else:
                if opt.find('_indexm') >= 0:
                    idx = 0
                    has_option = True
                    while has_option:
                        inner_idx = 0
                        has_inner = True
                        while has_inner:
                            opt_here = opt.replace('_indexm', f'_{idx}_{inner_idx}')
                            # print(opt_here,idx, inner_idx)
                            if self.options.has_option(section, opt_here):
                                value = self.read_option(section, opt, opt_here, poptions, idx, use_pow2_flags)

                                options_dict = self.assign_options(options_dict, opt_here, value)
                                inner_idx = inner_idx + 1
                            else:
                                has_inner = False
                            if inner_idx == 0:
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
        if opt.startswith('flag_index'):
            return self.get_dict_flag_index(values,idx,use_pow2_flags)
        return values

    def get_dict_flag_ranges_index(self, values, idx, use_pow2_flags):
        if use_pow2_flags:
            fvalue = int(math.pow(2, idx))
        else:
            fvalue = int(idx + 1)
        value_dict = {
            'flag_name': values[0],
            'flag_value': fvalue,
            'condition_list': [],
        }
        for idx in range(1,len(values)):
            val_condition = values[idx]
            if val_condition.strip().startswith('[') and val_condition.strip().endswith(']'):
                conditions = val_condition.strip()[1:-1].split(';')

                if len(conditions)==2:
                    value_dict['condition_list'].append({
                        'name_var': conditions[0].strip(),
                        'name_flag': conditions[1].strip(),
                        'flag_or_range': 'flag'
                    })
            if val_condition.strip().startswith('(') and val_condition.strip().endswith(')'):
                val_condition = val_condition.strip()
                conditions = val_condition.strip()[1:-1].split(';')
                if len(conditions)==3:
                    value_dict['condition_list'].append({
                        'name_var': conditions[0].strip(),
                        'min_val': float(conditions[1].strip()) if conditions[1].strip().lower()!='none' else None,
                        'max_val': float(conditions[2].strip()) if conditions[2].strip().lower()!='none' else None,
                        'flag_or_range': 'range'
                    })

        ###DEPRECATED
        # if len(values) == 3:
        #     value_dict = {
        #         'flag_var': values[0],
        #         'flag_name': values[1],
        #         'flag_value': fvalue,
        #     }
        # if len(values) >= 4:
        #     try:
        #         minV = float(values[2])
        #     except:
        #         minV = None
        #     try:
        #         maxV = float(values[3])
        #     except:
        #         maxV = None
        #     value_dict = {
        #         'is_default': False,
        #         'flag_var': values[0],
        #         'flag_name': values[1],
        #         'flag_value': fvalue,
        #         'min_range': minV,
        #         'max_range': maxV,
        #         'flag_condition': 'and'
        #     }
        #     if len(values) == 5:
        #         value_dict['flag_condition'] = values[4]
        return value_dict

    def get_dict_flag_index(self, values, idx, use_pow2_flags):
        if use_pow2_flags:
            fvalue = int(math.pow(2, idx))
        else:
            fvalue = int(idx + 1)
        value_dict = {
            'flag_name': values[0],
            'flag_value': fvalue,
            'condition_list': []
        }
        nconditions = int((len(values)-1)/2)

        for idx in range(nconditions):
            il = (2*idx)+1

            value_dict['condition_list'].append(
                {
                    'flag_var': values[il].strip()[1:],
                    'flag_list': [x.strip() for x in values[il+1].strip()[:-1].split(';')]
                }
            )


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




    def get_options_as_dict(self,section,poptions,required):
        if not self.is_valid():
            return None
        if poptions is None and required is None:
            soptions, required = self.get_retrieve_options(section)
            if soptions is None:
                return None
        result = {}
        if not self.options.has_section(section):
            for option in poptions:
                result[option] = poptions[option]['default']
        else:
            for option in poptions:
                result[option] = self.get_option(section,option,poptions,None,None)

        if required is not None:
            for r in required:
                if not r in result:
                    print(f'[ERROR] Option {r} is required in section {section} of the configuration file.')
                    return None
                if result[r] is None:
                    print(f'[ERROR] Option {section}/{r}  of the configuration file is required.')
                    if poptions[r]['type_param']=='file' and self.options.has_option(section,r):
                        print(f'[ERROR] {section}/{r}: {self.options[section][r]} does not exist or is not a valid file')
                    return None

        return result

    def get_option(self,section,option,poptions,default,type_param):

        list_values = None
        if poptions is not None and option in poptions.keys():
            if default is None  and 'default' in poptions[option].keys():
                default = poptions[option]['default']
            if type_param is None and 'type_param' in poptions[option].keys():
                type_param = poptions[option]['type_param']
            if 'list_values' in poptions[option].keys():
                list_values = poptions[option]['list_values']

        if type_param is None:
            return None

        value = self.get_value_param(section, option, default, type_param)

        if list_values is not None:
            if not value in list_values:
                print(f'[WARNING] Section/option {section}/{option}: {value} is not a valid value. It should be in the list: {list_values}')
                value = None

        return value

    def get_value(self, section, key):
        value = None
        if self.options.has_option(section, key):
            value = self.options[section][key]
            value = value.strip()
        return value

    def get_value_param(self, section, key, default, type):
        value = self.get_value(section, key)
        if value is None:
            return default


        return self.get_value_param_impl(value,type,default)

    def get_value_param_impl(self,value,type,default):

        if type == 'str':
            return value.strip(f'"')

        if type == 'file':
            file = value.strip(f'"')
            if not os.path.exists(file):
                return default
            else:
                return file

        if type == 'directory' or type=='output_path':
            directory = value.strip(f'"')
            if not os.path.isdir(directory):
                try:
                    os.mkdir(directory)
                    return directory
                except:
                    return default
            else:
                return directory

        if type== 'input_path':
            input_path = value.strip(f'"')
            if not os.path.isdir(input_path):
                if default is not None and os.path.isdir(default):
                    return default
                else:
                    try:
                        os.mkdir(input_path)
                        return input_path
                    except:
                        return None
            else:
                return input_path

        if type == 'int':
            return int(value.strip(f'"'))

        if type == 'float':
            return float(value.strip(f'"'))

        if type == 'boolean':
            value = value.strip(f'"')
            if value == '1' or value.upper() == 'TRUE':
                return True
            elif value == '0' or value.upper() == 'FALSE':
                return False
            else:
                return True

        if type == 'rrslist':
            #list_str = value.split(',')
            list_str = [s.strip().strip('"') for s in re.split(r',(?=(?:[^"]*"[^"]*")*[^"]*$)', value)]
            list = []
            for vals in list_str:
                vals = vals.replace('.', '_')
                list.append(f'RRS{vals}')
            return list

        if type == 'strlist':
            list = [s.strip().strip('"') for s in re.split(r',(?=(?:[^"]*"[^"]*")*[^"]*$)', value)]

            return list

        if type == 'floatlist':
            list_str = [s.strip().strip('"') for s in re.split(r',(?=(?:[^"]*"[^"]*")*[^"]*$)', value)]
            list = []
            for vals in list_str:
                list.append(float(vals))
            return list

        if type == 'intlist':
            list_str = [s.strip().strip('"') for s in re.split(r',(?=(?:[^"]*"[^"]*")*[^"]*$)', value)]
            list = []
            for vals in list_str:
                list.append(int(vals))
            return list
