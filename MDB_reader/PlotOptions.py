import os.path
import MDBPlotDefaults as defaults


class PlotOptions:

    def __init__(self, options, config_file):

        if options is None and config_file is not None and os.path.exists(config_file):
            import configparser
            options = configparser.ConfigParser()
            options.read(config_file)

        self.options = options
        self.output_path = None



        ##GLOBAL OPTIONS
        self.global_options = {}#defaults.global_option
        self.mu_valid_variable = 'mu_valid'
        self.format_image = 'png'
        self.image_resolution = 300

        self.valid_stats = {
            'N': 0,
            'NMU': 0,
            'slope': 0.0,
            'intercept': 0.0,
            'PCC(r)': 0.0,
            'p_value': 0.0,
            'std_err': 0.0,
            'RMSD': 0.0,
            'RPD': 0.0,
            'APD': 0.0,
            'BIAS': 0.0,
            'DETER(r2)': 0.0,
            'SLOPE': 0.0,
            'OFFSET': 0.0,
            'XAVG': 0.0,
            'YAVG': 0.0,
            'CPRMSE': 0.0,
            'MAE': 0.0
        }


    def set_global_options(self):
        section = 'GLOBAL_OPTIONS'
        self.global_options = {}
        for goption in defaults.global_options:
            default = defaults.global_options[goption]['default']
            type = defaults.global_options[goption]['type']
            self.global_options[goption] = self.get_value_param(section,goption.strip(),default,type)
            if type=='str' and 'values' in defaults.global_options[goption].keys():
                values = defaults.global_options[goption]['values']
                if not self.global_options[goption] in values:
                    print(f'[ERROR] [{section}] {self.global_options[goption]} is not a valid  value for {goption}. Valid values: {values} ')

        #self.output_path = self.get_value_param(section, 'output_path', self.output_path, 'directory')
        # self.mu_valid_variable = self.get_value_param(section, 'mu_valid_variable', self.mu_valid_variable, 'str')


    def get_list_figures(self):
        sections = self.options.sections()
        list_figures = []
        for s in sections:
            apply = self.get_value_param(s, 'apply', False, 'boolean')
            if apply:
                list_figures.append(s)
        return list_figures


    def get_options(self, section):
        options_out = {'apply': self.get_value_param(section, 'apply', False, 'boolean')}
        if not options_out['apply']:
            return options_out
        options_out['type'] = self.get_value_param(section, 'type', None, 'str')
        if options_out['type'] is None:
            return options_out
        options_out['name'] = section


        options_out['multiple_plot'] = self.get_value_param(section, 'multiple_plot', None, 'str')
        # if options_out['type'] == 'csvtable':
        #     options_out = self.get_options_csv(section,options_out)
        if options_out['type'] == 'scatterplot':
            # options_out = self.get_group_options(section, options_out)
            # options_out = self.get_select_options(section, options_out)
            options_out = self.get_options_scatterplot(section, options_out)
        if options_out['type'].startswith('statstable'):
            options_out = self.get_options_csv_statstable(section, options_out)
        # if options_out['type'] == 'statswlplot':
        #     options_out = self.get_select_options(section, options_out)
        #     options_out = self.get_options_statswlplot(section, options_out)
        # if options_out['type'] == 'spectraplot':
        #     options_out['type_rrs'] = self.get_value_param(section, 'type_rrs', 'ins', 'str')
        #     if options_out['type_rrs'].startswith('flag'):
        #         options_out = self.get_group_options(section, options_out)
        #     options_out = self.get_options_spectraplot(section, options_out)
        #     options_out = self.get_select_options(section, options_out)
        # if options_out['type'] == 'flagplot':
        #     options_out = self.get_options_flag(section, options_out)

        return options_out

    def get_options_csv_statstable(self, section, options_out):
        options_out['xvar'] = self.get_value_param(section, 'xvar', 'mu_ins_rrs', 'str')
        options_out['yvar'] = self.get_value_param(section, 'yvar', 'mu_sat_rrs', 'str')
        options_out['params'] = self.get_value_param(section, 'params', self.valid_stats.keys(), 'strlist')
        options_out['log_scale'] = self.get_value_param(section, 'log_scale', False, 'boolean')
        options_out['use_rhow'] = self.get_value_param(section, 'use_rhow', False, 'boolean')
        options_out['flag'] = self.get_value_param(section, 'flag', 'GLOBAL', 'str')
        if self.output_path is not None:
            name_default = options_out['name'] + '.csv'
            file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(section, 'file_out', file_out_default, 'str')
        return options_out

    def get_options_scatterplot(self, section, options_out):

        options_out['type_scatterplot'] = self.get_value_param(section, 'type_scatterplot', 'rrs', 'str')

        options_out['xvar'] = self.get_value_param(section, 'xvar', 'mu_ins_rrs', 'str')
        options_out['yvar'] = self.get_value_param(section, 'yvar', 'mu_sat_rrs', 'str')

        options_out['legend'] = self.get_value_param(section, 'legend', True, 'boolean')
        options_out['legend_values'] = self.get_value_param(section, 'legend_values', None, 'strlist')
        options_out['include_stats'] = self.get_value_param(section, 'include_stats', False, 'boolean')

        options_out['regression_line_groups'] = self.get_value_param(section, 'regression_line_groups', False,
                                                                     'boolean')
        options_out['apply_density'] = self.get_value_param(section, 'apply_density', True, 'boolean')
        options_out['title'] = self.get_value_param(section, 'title', None, 'str')
        print(self.global_options)
        if self.global_options['output_path'] is not None:
            name_default = options_out['name'] + '.' + self.global_options['fig_extension']
            file_out_default = os.path.join(self.global_options['output_path'],name_default)
            #name_default = options_out['name'] + '.' + self.format_image
            #file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(section, 'file_out', file_out_default, 'str')
        options_out['log_scale'] = self.get_value_param(section, 'log_scale', False, 'boolean')
        options_out['use_rhow'] = self.get_value_param(section, 'use_rhow', False, 'boolean')
        options_out['min_xy'] = self.get_value_param(section, 'min_xy', None, 'float')
        options_out['max_xy'] = self.get_value_param(section, 'max_xy', None, 'float')
        options_out['ticks'] = self.get_value_param(section, 'ticks', None, 'floatlist')
        options_out['fontsizeaxis'] = self.get_value_param(section, 'fontsizeaxis', 12, 'float')
        options_out['fontsizelabels'] = self.get_value_param(section, 'fontsizelabels', 12, 'float')
        options_out['fontsizetitle'] = self.get_value_param(section, 'fontsizetitle', 12, 'float')
        options_out['fontsizestats'] = self.get_value_param(section, 'fontsizestats', 12, 'float')
        sfdefault = None
        unitsdefault = None
        if options_out['type_scatterplot'] == 'rrs':
            unitsdefault = r'sr$^-$$^1$'
            sfdefault = 1000
        if options_out['type_scatterplot'] == 'chla':
            unitsdefault = r'mg m$^-$$^3$'
            sfdefault = 1
            options_out['log_scale'] = True
        if options_out['type_scatterplot'] == 'kd':
            unitsdefault = r'm$^-$$^1$'
            sfdefault = 1
            options_out['log_scale'] = True
        xlabeldefault = defaults.xlabel_default
        ylabeldefault = defaults.ylabel_default
        options_out['scale_factor'] = self.get_value_param(section, 'scale_factor', sfdefault, 'float')
        options_out['xlabel'] = self.get_value_param(section, 'xlabel', xlabeldefault, 'str')
        options_out['ylabel'] = self.get_value_param(section, 'ylabel', ylabeldefault, 'str')
        options_out['units'] = self.get_value_param(section, 'units', unitsdefault, 'str')
        options_out['identity_line'] = self.get_value_param(section, 'identity_line', True, 'boolean')
        options_out['regression_line'] = self.get_value_param(section, 'regression_line', True, 'boolean')
        # marker, markersize, color, edgecolor, linewidth
        # o, 25, 'black', None, None
        options_out['marker'] = self.get_value_param(section, 'marker', ['o'], 'strlist')
        options_out['markersize'] = self.get_value_param(section, 'markersize', [25], 'intlist')
        options_out['color'] = self.get_value_param(section, 'color', ['black'], 'strlist')

        edgeColorDefault = None
        lineWidthDefault = None
        # if options_out['groupBy'] is not None:
        # if options_out['groupType'] == 'rrs':
        #     edgeColorDefault = 'gray'
        #     lineWidthDefault = 1.5
        # if options_out['groupType'] == 'flag':
        #     edgeColorDefault = 'black'
        #     lineWidthDefault = 0.25
        options_out['edgecolor'] = self.get_value_param(section, 'edgecolor', [edgeColorDefault], 'strlist')
        options_out['linewidth'] = self.get_value_param(section, 'linewidth', [lineWidthDefault], 'floatlist')

        options_out['xfigsize'] = self.get_value_param(section, 'xfigsize', 7, 'float')
        options_out['yfigsize'] = self.get_value_param(section, 'yfigsize', 7, 'float')
        options_out['widthspace'] = self.get_value_param(section, 'widthspace', 0.1, 'float')
        options_out['heightspace'] = self.get_value_param(section, 'heightspace', 0.1, 'float')
        options_out['stat_list'] = self.get_value_param(section, 'stat_list', None, 'strlist')
        options_out['stats_xpos'] = self.get_value_param(section, 'stats_xpos', 0.05, 'float')
        options_out['stats_ypos'] = self.get_value_param(section, 'stats_ypos', 0.70, 'float')
        options_out['individual_axis'] = self.get_value_param(section, 'individual_axis', False, 'boolean')
        return options_out

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
