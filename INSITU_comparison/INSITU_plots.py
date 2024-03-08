from INSITU_comparison import INSITUCOMPARISON
import numpy as np
import os
import INSITU_plots_defaults as defaults
from scipy.stats import gaussian_kde


class INSITU_plots():

    def __init__(self, insitucomparison):

        self.ic = insitucomparison
        self.output_path = None
        self.format_image = 'tif'

    def set_global_options(self, options):
        section = 'GLOBAL_OPTIONS'
        self.output_path = self.get_value_param(options, section, 'output_path', self.output_path, 'directory')

    def plot_from_options(self, options):
        self.set_global_options(options)
        plot_list = list(options.sections())
        for plot in plot_list:
            if plot == 'GLOBAL_OPTIONS':
                continue
            options_out = self.get_options(options, plot)

            # print(options_out)
            if options_out['apply']:
                print(f'[INFO] Starting plot: {plot}')
                self.plot_from_options_impl(options_out)

    def plot_from_options_impl(self, options_out):

        if options_out['type'] == 'scatterplot':
            self.plot_scatter_plot(options_out)

    def plot_scatter_plot(self, options_out):
        if options_out['type_scatterplot'] == 'single':
            self.ic.set_data_scatterplot_single(options_out['xvariable'], options_out['yvariable'],
                                                options_out['wlvariableref'])
            self.ic.compute_stats(False, False)
            self.plot_scatter_plot_impl(options_out,None)

        if options_out['type_scatterplot'] == 'spectral':
            if not options_out['selectByWavelength']:

                self.ic.set_data_scatterplot_spectral(None,options_out['xvariable'], options_out['yvariable'],
                                                options_out['wlvariablemu'])

                self.ic.compute_stats(False,False)
                self.plot_scatter_plot_impl(options_out,None)

            elif options_out['selectByWavelength']:

                wl_values = options_out['wl_values']
                if wl_values is None:
                    wl_values = self.ic.get_wl_list(options_out['wlvariablemu'])
                file_out_base = options_out['file_out'][:-4]
                for wl in wl_values:
                    wls = self.ic.get_wl_str_from_wl(wl)
                    self.ic.set_data_scatterplot_spectral(wl, options_out['xvariable'], options_out['yvariable'],
                                                          options_out['wlvariablemu'])

                    self.ic.compute_stats(False, False)
                    options_out['file_out'] = f'{file_out_base}{wls}.{self.format_image}'
                    self.plot_scatter_plot_impl(options_out,wl)


    def plot_scatter_plot_impl(self, options_out,wl):
        file_out = options_out['file_out']
        marker = options_out['marker']
        markersize = options_out['markersize']
        edgecolor = options_out['edgecolor']
        linewidth = options_out['linewidth']
        fontsizeaxis = options_out['fontsizeaxis']
        xlabel = options_out['xlabel']
        ylabel = options_out['ylabel']
        min_xy = options_out['min_xy']
        max_xy = options_out['max_xy']
        ticks = options_out['ticks']
        from MDB_reader.PlotScatter import PlotScatter
        plot = PlotScatter()
        plot.close_plot()
        plot.start_plot()

        xhere = np.asarray(self.ic.xdata, dtype=np.float)
        yhere = np.asarray(self.ic.ydata, dtype=np.float)
        # print('-----------------------------------')
        # print(xhere.shape)
        # print(yhere.shape)

        if options_out['groupBy'] is None:
            if options_out['apply_density']:
                try:
                    xy = np.vstack([xhere, yhere])
                    z = gaussian_kde(xy)(xy)
                    idx = z.argsort()
                    xhere, yhere, z = xhere[idx], yhere[idx], z[idx]
                    plot.set_cmap('jet')
                    plot.plot_data(xhere, yhere, marker[0], markersize[0], z, edgecolor[0], linewidth[0])
                except:
                    plot.plot_data(xhere, yhere, marker[0], markersize[0], 'k', edgecolor[0], linewidth[0])

        if options_out['title'] is not None:
            plot.set_title(options_out['title'])
        if min_xy is not None and max_xy is not None:
            plot.set_limits(min_xy, max_xy)
        if ticks is not None:
            plot.set_ticks(ticks, fontsizeaxis)

        plot.set_xaxis_title(xlabel)
        plot.set_yaxis_title(ylabel)

        plot.plot_regress_line(self.ic.xregress, self.ic.yregress, 'k')
        plot.plot_identity_line()

        if options_out['include_stats']:
            stat_list = options_out['stat_list']
            str_stats = defaults.get_str_stat_list(self.ic.valid_stats, stat_list, wl)

            plot.plot_text(options_out['stats_xpos'], options_out['stats_ypos'], str_stats)

        plot.set_equal_apect()

        plot.save_fig(file_out)
        plot.close_plot()

    ##OPTIONS
    def get_options(self, options, section):
        options_out = {'apply': self.get_value_param(options, section, 'apply', False, 'boolean')}
        if not options_out['apply']:
            return options_out
        options_out['type'] = self.get_value_param(options, section, 'type', None, 'str')
        if options_out['type'] is None:
            return options_out
        options_out['name'] = section
        options_out['multiple_plot'] = self.get_value_param(options, section, 'multiple_plot', None, 'str')

        # options_out = self.get_anot_options(options, section, options_out)
        if options_out['type'] == 'scatterplot':
            options_out = self.get_group_options(options, section, options_out)
            options_out = self.get_select_options(options, section, options_out)
            options_out = self.get_options_scatterplot(options, section, options_out)

        return options_out

    def get_group_options(self, options, section, options_out):
        options_out['groupBy'] = self.get_value_param(options, section, 'groupBy', None, 'str')
        options_out['groupValues'] = None
        # if options_out['groupBy'] is not None:
        #     options_out['groupType'] = 'float'
        #     var_name = options_out['groupBy']
        #     virtual_flag = False
        # if not var_name in self.mrfile.variables:
        #     print(
        #         f'[WARNING] {var_name} is not a variable defined in the MDB file. Checking {var_name} as virtual flag')
        #     virtual_flag = True
        #     flag_values, flag_meanings, flag_array = self.get_virtual_flag(options, var_name)
        #     print(flag_meanings)
        #     print(flag_values)
        #     print(flag_array.shape)
        #     options_out[var_name] = {
        #         'flag_values': flag_values,
        #         'flag_meanings': flag_meanings,
        #         'flag_array': flag_array
        #     }
        #     options_out['groupType'] = 'flag'
        #     print(f'[INFO] Set virtual flag: {var_name}')
        #
        # if var_name.startswith('flag') and not virtual_flag:
        #     flag_values = self.mrfile.nc.variables[var_name].flag_values
        #     flag_meanings_list = self.mrfile.nc.variables[var_name].flag_meanings.split(' ')
        #     flag_meanings = [x.strip() for x in flag_meanings_list]
        #     options_out[var_name] = {
        #         'flag_values': flag_values,
        #         'flag_meanings': flag_meanings
        #     }
        #     options_out['groupType'] = 'flag'
        #
        # if options_out['groupType'] == 'flag':
        #     flag_list_config = self.get_value_param(options, section, 'groupValues', None, 'strlist')
        #     if flag_list_config is None:
        #         options_out['groupValues'] = flag_values
        #     else:
        #         flag_values_config = []
        #         for flag_config in flag_list_config:
        #             if flag_config.strip() == 'GLOBAL':
        #                 flag_values_config.append(-1)
        #                 continue
        #             try:
        #                 iflag = flag_meanings.index(flag_config.strip())
        #                 flag_values_config.append(flag_values[iflag])
        #             except:
        #                 print(f'[WARNING] Flag {flag_config.strip()} is not in the list')
        #         options_out['groupValues'] = flag_values_config
        # if options_out['groupType'] == 'float':
        #     group_values = list(np.unique(np.array(self.mrfile.nc.variables[var_name])))
        #     options_out['groupValues'] = self.get_value_param(options, section, 'groupValues', group_values,
        #                                                       'floatlist')

        return options_out

    def get_select_options(self, options, section, options_out):

        options_out['selectByWavelength'] = self.get_value_param(options, section, 'selectByWavelength', False,
                                                                 'boolean')
        wl_values = self.get_value_param(options, section, 'wlvalues', None, 'floatlist')
        # if wl_values is None and options_out['selectByWavelength']:
        #     wl_values = list(np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])))
        options_out['wl_values'] = wl_values
        options_out['selectBy'] = self.get_value_param(options, section, 'selectBy', None, 'str')
        options_out['selectValues'] = None
        # if options_out['selectBy'] is not None:
        #     options_out['selectType'] = 'float'
        #     var_name = options_out['selectBy']
        #     virtual_flag = False
        #     if not var_name in self.mrfile.variables:
        #         print(
        #             f'[WARNING] {var_name} is not a variable defined in the MDB file. Checking {var_name} as virtual flag')
        #         virtual_flag = True
        #         flag_values, flag_meanings, flag_array = self.get_virtual_flag(options, var_name)
        #         options_out[var_name] = {
        #             'flag_values': flag_values,
        #             'flag_meanings': flag_meanings,
        #             'flag_array': flag_array
        #         }
        #         options_out['selectType'] = 'flag'
        #         print(f'[INFO] Set virtual flag: {var_name}')
        #
        #     if var_name.startswith('flag') and not virtual_flag:
        #         flag_values = self.mrfile.nc.variables[var_name].flag_values
        #         flag_meanings_list = self.mrfile.nc.variables[var_name].flag_meanings.split(' ')
        #         flag_meanings = [x.strip() for x in flag_meanings_list]
        #         options_out[var_name] = {
        #             'flag_values': flag_values,
        #             'flag_meanings': flag_meanings
        #         }
        #         options_out['selectType'] = 'flag'
        #
        #     if options_out['selectType'] == 'flag':
        #         flag_list_config = self.get_value_param(options, section, 'selectValues', None, 'strlist')
        #         if flag_list_config is None:
        #             options_out['selectValues'] = flag_values
        #         else:
        #             flag_values_config = []
        #             for flag_config in flag_list_config:
        #                 if flag_config.strip() == 'GLOBAL':
        #                     flag_values_config.append(-1)
        #                     continue
        #                 try:
        #                     iflag = flag_meanings.index(flag_config.strip())
        #                     flag_values_config.append(flag_values[iflag])
        #                 except:
        #                     print(f'[WARNING] Flag {flag_config.strip()} is not in the list')
        #             options_out['selectValues'] = flag_values_config
        #     if options_out['selectType'] == 'float':
        #         group_values = list(np.unique(np.array(self.mrfile.nc.variables[var_name])))
        #         options_out['selectValues'] = self.get_value_param(options, section, 'selectValues', group_values,
        #                                                            'floatlist')

        return options_out

    def get_options_scatterplot(self, options, section, options_out):
        options_out['type_scatterplot'] = self.get_value_param(options, section, 'type_scatterplot', 'spectral', 'str')
        options_out['xvariable'] = self.get_value_param(options, section, 'xvariable', 'mu_HYPSTAR_TO_AERONET_Lw',
                                                        'str')
        options_out['yvariable'] = self.get_value_param(options, section, 'yvariable', 'mu_AERONET_Lw', 'str')
        options_out['wlvariableref'] = self.get_value_param(options, section, 'wlvariableref',
                                                            'AERONET_nominal_wavelengths', 'str')

        options_out['wlvariablemu'] = self.get_value_param(options, section, 'wlvariablemu',
                                                            'mu_wavelength', 'str')

        options_out['legend'] = self.get_value_param(options, section, 'legend', True, 'boolean')
        options_out['legend_values'] = self.get_value_param(options, section, 'legend_values', None, 'strlist')
        options_out['include_stats'] = self.get_value_param(options, section, 'include_stats', False, 'boolean')

        options_out['regression_line_groups'] = self.get_value_param(options, section, 'regression_line_groups', False,
                                                                     'boolean')
        options_out['apply_density'] = self.get_value_param(options, section, 'apply_density', True, 'boolean')
        options_out['title'] = self.get_value_param(options, section, 'title', None, 'str')
        if self.output_path is not None:
            name_default = options_out['name'] + '.' + self.format_image
            file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(options, section, 'file_out', file_out_default, 'str')
        options_out['log_scale'] = self.get_value_param(options, section, 'log_scale', False, 'boolean')
        options_out['use_rhow'] = self.get_value_param(options, section, 'use_rhow', False, 'boolean')
        options_out['min_xy'] = self.get_value_param(options, section, 'min_xy', None, 'float')
        options_out['max_xy'] = self.get_value_param(options, section, 'max_xy', None, 'float')
        options_out['ticks'] = self.get_value_param(options, section, 'ticks', None, 'floatlist')
        options_out['fontsizeaxis'] = self.get_value_param(options, section, 'fontsizeaxis', 12, 'float')
        options_out['fontsizelabels'] = self.get_value_param(options, section, 'fontsizelabels', 12, 'float')
        options_out['fontsizetitle'] = self.get_value_param(options, section, 'fontsizetitle', 12, 'float')
        options_out['fontsizestats'] = self.get_value_param(options, section, 'fontsizestats', 12, 'float')
        sfdefault = None
        unitsdefault = None
        if options_out['type_scatterplot'] == 'rrs':
            unitsdefault = r'sr$^-$$^1$'
            if not options_out['use_rhow']:
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
        options_out['scale_factor'] = self.get_value_param(options, section, 'scale_factor', sfdefault, 'float')
        options_out['xlabel'] = self.get_value_param(options, section, 'xlabel', xlabeldefault, 'str')
        options_out['ylabel'] = self.get_value_param(options, section, 'ylabel', ylabeldefault, 'str')
        options_out['units'] = self.get_value_param(options, section, 'units', unitsdefault, 'str')
        options_out['identity_line'] = self.get_value_param(options, section, 'identity_line', True, 'boolean')
        options_out['regression_line'] = self.get_value_param(options, section, 'regression_line', True, 'boolean')
        # marker, markersize, color, edgecolor, linewidth
        # o, 25, 'black', None, None
        options_out['marker'] = self.get_value_param(options, section, 'marker', ['o'], 'strlist')
        options_out['markersize'] = self.get_value_param(options, section, 'markersize', [25], 'intlist')
        options_out['color'] = self.get_value_param(options, section, 'color', ['black'], 'strlist')

        edgeColorDefault = None
        lineWidthDefault = None
        if options_out['groupBy'] is not None:
            if options_out['groupType'] == 'rrs':
                edgeColorDefault = 'gray'
                lineWidthDefault = 1.5
            if options_out['groupType'] == 'flag':
                edgeColorDefault = 'black'
                lineWidthDefault = 0.25
        options_out['edgecolor'] = self.get_value_param(options, section, 'edgecolor', [edgeColorDefault], 'strlist')
        options_out['linewidth'] = self.get_value_param(options, section, 'linewidth', [lineWidthDefault], 'floatlist')

        options_out['xfigsize'] = self.get_value_param(options, section, 'xfigsize', 7, 'float')
        options_out['yfigsize'] = self.get_value_param(options, section, 'yfigsize', 7, 'float')
        options_out['widthspace'] = self.get_value_param(options, section, 'widthspace', 0.1, 'float')
        options_out['heightspace'] = self.get_value_param(options, section, 'heightspace', 0.1, 'float')
        options_out['stat_list'] = self.get_value_param(options, section, 'stat_list', None, 'strlist')
        options_out['stats_xpos'] = self.get_value_param(options, section, 'stats_xpos', 0.03, 'float')
        options_out['stats_ypos'] = self.get_value_param(options, section, 'stats_ypos', 0.70, 'float')
        options_out['individual_axis'] = self.get_value_param(options, section, 'individual_axis', False, 'boolean')
        return options_out

    def get_value(self, options, section, key):
        value = None
        if options.has_option(section, key):
            value = options[section][key]
        return value

    def get_value_param(self, options, section, key, default, type):

        value = self.get_value(options, section, key)
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
