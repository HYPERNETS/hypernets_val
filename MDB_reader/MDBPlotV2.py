import MDBPlotDefaults as defaults
from MDBFile import MDBFile
import numpy as np
import os


class MDBPlot:

    def __init__(self, path_mdbr_file):
        self.mrfile = MDBFile(path_mdbr_file)

        self.VALID = self.mrfile.VALID

        # plot general options
        self.title = ''
        self.file_name_base = ''
        self.format_image = 'jpg'

        self.output_path = None

        # validation scatterplot
        self.xdata = []
        self.ydata = []
        self.groupdata = []
        self.wldata = []
        self.xregress = []
        self.yregress = []

        # validation stats
        self.valid_stats = {
            'N': 0,
            'slope': 0.0,
            'intercept': 0.0,
            'r_value': 0.0,
            'p_value': 0.0,
            'std_err': 0.0,
            'rmse_val': 0.0,
            'mean_rel_diff': 0.0,
            'mean_abs_rel_diff': 0.0,
            'bias': 0.0,
            'r2': 0.0,
            'slope_typeII': 0.0,
            'offset_typeII': 0.0,
            'XAVG': 0.0,
            'YAVG': 0.0,
            'CPRMSE': 0.0,
            'MAE': 0.0
        }
        self.df_valid_stats = None

        self.units = r'sr$^-$$^1$'
        self.log_scale = False

        self.satellite = 'S3'
        self.platform = 'AB'

    def plot_from_options(self, options):
        plot_list = list(options.sections())
        for plot in plot_list:
            options_out = self.get_options(options, plot)
            print(options_out)
            if options_out['apply']:
                self.plot_from_options_impl(options_out)

    def plot_from_options_impl(self, options_out):
        if options_out['type'] == 'scatterplot':
            self.set_data_scatterplot(options_out['groupBy'], options_out['selectBy'], options_out['selectValues'])
            self.plot_scatter_plot(options_out)

    def set_data_scatterplot(self, groupBy, selectBy, valsSelect):
        rrs_ins = np.array(self.mrfile.nc.variables['mu_ins_rrs'])
        rrs_sat = np.array(self.mrfile.nc.variables['mu_sat_rrs'])
        id_all = np.array(self.mrfile.nc.variables['mu_satellite_id'])
        mu_valid = np.array(self.mrfile.nc.variables['mu_valid'])

        valid_all = self.get_array_all_from_arraymu(id_all, mu_valid)

        self.xdata = rrs_ins[valid_all == 1]
        self.ydata = rrs_sat[valid_all == 1]

    def get_array_all_from_arraymu(self, id_all_array, mu_array):
        array_out = np.zeros(id_all_array.shape, dtype=mu_array.dtype)
        for id in range(len(mu_array)):
            array_out[id_all_array == id] = mu_array[id]
        return array_out

    # MAIN FUNCTION TO PLOT SCATTERPLOT
    def plot_scatter_plot(self, options):

        if options['include_stats']:
            print('CALCULA ESTADISTICAS')
            # self.compute_statistics()

        ngroup = 1
        str_legend = []
        if len(self.groupdata) > 0:
            if options['groupValues'] is None:
                groupValues = np.unique(np.array(self.groupdata))
            else:
                groupValues = options['groupValues']
            ngroup = len(groupValues)
            if ngroup > 1 and options['legend'] and options['groupType'] == 'float':
                for g in ngroup:
                    str_legend.append(f'{g:.2f}')

        from PlotScatter import PlotScatter
        from scipy.stats import gaussian_kde
        plot = PlotScatter()
        plot.close_plot()
        plot.start_plot()

        if options['scale_factor'] is not None:
            self.xdata = self.xdata * options['scale_factor']
            self.ydata = self.ydata * options['scale_factor']
            if len(self.yregress) > 0 and len(self.xregress) > 0:
                self.yregress = np.array(self.yregress) * options['scale_factor']
                self.xregress = np.array(self.xregress) * options['scale_factor']

        if ngroup > 1:
            for g in ngroup:
                color = defaults.get_color_ref(g)
                xhere = self.xdata[self.groupdata == g]
                yhere = self.ydata[self.groupdata == g]
                # marker, markersize, color, edgecolor, linewidth
                # None, None, color, 'gray', 1.5
                plot.plot_data(xhere, yhere, options['marker'], options['markersize'], color, options['edgecolor'], options['linewidth'])
        else:  # density scatter plot

            xhere = np.asarray(self.xdata, dtype=np.float)
            yhere = np.asarray(self.ydata, dtype=np.float)

            # Density
            if options['apply_density']:
                xy = np.vstack([xhere, yhere])
                try:
                    z = gaussian_kde(xy)(xy)
                    idx = z.argsort()
                    xhere, yhere, z = xhere[idx], yhere[idx], z[idx]
                    plot.plot_data(xhere, yhere, options['marker'], options['markersize'], z, options['edgecolor'], options['linewidth'])
                    plot.set_cmap('jet')
                except:
                    plot.plot_data(xhere, yhere, options['marker'], options['markersize'], options['color'], options['edgecolor'], options['linewidth'])

            else:
                plot.plot_data(xhere, yhere, options['marker'], options['markersize'], options['color'], options['edgecolor'], options['linewidth'])

        if options['log_scale']:
            plot.set_log_scale()
            min_xy = options['min_xy'] = 0.1
            max_xy = options['max_xy'] = 100
        else:
            if options['min_xy'] is None and options['max_xy'] is None:
                min_x_data = np.min(self.xdata)
                min_y_data = np.min(self.ydata)
                min_xy = np.floor(np.min([min_x_data, min_y_data]))
                max_x_data = np.max(self.xdata)
                max_y_data = np.max(self.ydata)
                max_xy = np.ceil(np.max([max_x_data, max_y_data]))
            else:
                min_xy = options['min_xy']
                max_xy = options['max_xy']

        plot.set_limits(min_xy, max_xy)
        plot.set_xaxis_title(options['xlabel'])
        plot.set_yaxis_title(options['ylabel'])
        plot.set_equal_apect()
        if options['legend'] and len(str_legend) > 0:
            plot.set_legend(str_legend)
        if options['identity_line']:
            plot.plot_identity_line()

        # if include_stats:
        #     if self.units.startswith('mg'):  # chla
        #         str0 = 'N={:d}\nRMSD={:,.2f} UNITS\nAPD={:,.0f}%\nRPD={:,.0f}%\n$r^2$={:,.2f}\nbias={:,.2f} UNITS' \
        #             .format(self.valid_stats['N'],
        #                     self.valid_stats['rmse_val'],
        #                     self.valid_stats['mean_abs_rel_diff'],
        #                     self.valid_stats['mean_rel_diff'],
        #                     self.valid_stats['r_value'] ** 2,
        #                     self.valid_stats['bias'])
        #         str0 = str0.replace('UNITS', self.units)
        #         plot.plot_text(0.05, 0.68, str0)
        #
        #     if self.units.startswith('sr'):  # rrs
        #         str0 = 'N={:d}\nRMSD={:,.1e} UNITS\nAPD={:,.0f}%\nRPD={:,.0f}%\n$r^2$={:,.2f}\nbias={:,.1e} UNITS' \
        #             .format(self.valid_stats['N'],
        #                     self.valid_stats['rmse_val'],
        #                     self.valid_stats['mean_abs_rel_diff'],
        #                     self.valid_stats['mean_rel_diff'],
        #                     self.valid_stats['r_value'] ** 2,
        #                     self.valid_stats['bias'])
        #         str0 = str0.replace('UNITS', self.units)
        #         plot.plot_text(0.05, 0.70, str0)
        #
        #     if not self.log_scale:
        #         plot.plot_regress_line(self.xregress, self.yregress, 'black')

        if options['title'] is not None:
            title_here = options['title']
            # if nwl == 1 and not title_here.lower().find('in situ') >= 0:
            #     w = wl[0]
            #     strwl = f' λ = {w:0.2f} nm \n'
            #     title_here = self.title + strwl
            plot.set_title(title_here)

        if not options['file_out'] is None:
            plot.save_fig(options['file_out'])
            plot.close_plot()

    def get_options(self, options, section):
        options_out = {'apply': self.get_value_param(options, section, 'apply', False, 'boolean')}
        if not options_out['apply']:
            return options_out

        options_out['type'] = self.get_value_param(options, section, 'type', None, 'str')
        if options_out['type'] is None:
            return options_out
        if options_out['type'] == 'scatterplot':
            options_out['name'] = section
            options_out = self.get_group_options(options, section, options_out)
            options_out = self.get_select_options(options, section, options_out)
            options_out = self.get_options_scatterplot(options, section, options_out)

        return options_out

    def get_group_options(self, options, section, options_out):
        options_out['groupBy'] = self.get_value_param(options, section, 'groupBy', None, 'str')
        if options_out['groupBy'] is not None:
            options_out['groupValues'] = None
            options_out['groupType'] = 'float'
            if options_out['groupBy'].startswith('flag'):
                options_out['groupValues'] = self.get_value_param(options, section, 'groupValues', None, 'strlist')
                options_out['groupType'] = 'flag'
            if options_out['groupType'] == 'float':
                options_out['groupValues'] = self.get_value_param(options, section, 'groupValues', None, 'floatlist')
        return options_out

    def get_select_options(self, options, section, options_out):
        options_out['selectBy'] = self.get_value_param(options, section, 'selectBy', None, 'str')
        options_out['selectValues'] = None
        if options_out['selectBy'] is not None:
            options_out['selectType'] = 'float'
            if options_out['selectBy'].startswith('flag'):
                options_out['selectValues'] = self.get_value_param(options, section, 'selectValues', None, 'strlist')
                options_out['selectType'] = 'flag'
            if options_out['selectType'] == 'float':
                options_out['selectValues'] = self.get_value_param(options, section, 'selectValues', None, 'floatlist')
        return options_out

    def get_options_scatterplot(self, options, section, options_out):
        options_out['type_scatterplot'] = self.get_value_param(options, section, 'type_scatterplot', 'rrs', 'str')
        options_out['legend'] = self.get_value_param(options, section, 'legend', True, 'boolean')
        options_out['legend_values'] = self.get_value_param(options, section, 'legend_values', None, 'strlist')
        options_out['include_stats'] = self.get_value_param(options, section, 'include_stats', False, 'boolean')
        options_out['apply_density'] = self.get_value_param(options, section, 'apply_density', True, 'boolean')
        options_out['title'] = self.get_value_param(options, section, 'title', None, 'str')
        if self.output_path is not None:
            name_default = options_out['name'] + '.' + self.format_image
            file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(options, section, 'file_out', file_out_default, 'str')
        options_out['log_scale'] = self.get_value_param(options, section, 'log_scale', False, 'boolean')
        options_out['min_xy'] = self.get_value_param(options, section, 'min_xy', None, 'float')
        options_out['max_xy'] = self.get_value_param(options, section, 'max_xy', None, 'float')
        sfdefault = None
        if options_out['type_scatterplot'] == 'rrs':
            sfdefault = 1000
            xlabeldefault = defaults.xlabel_default
            ylabeldefault = defaults.ylabel_default
        options_out['scale_factor'] = self.get_value_param(options, section, 'scale_factor', sfdefault, 'float')
        options_out['xlabel'] = self.get_value_param(options, section, 'xlabel', xlabeldefault, 'float')
        options_out['ylabel'] = self.get_value_param(options, section, 'ylabel', ylabeldefault, 'float')
        options_out['identity_line'] = self.get_value_param(options,section,'identity_line',True,'boolean')

        #marker, markersize, color, edgecolor, linewidth
        #None, 25, 'black', None, None
        options_out['marker'] = self.get_value_param(options,section,'marker',None,'str')
        options_out['markersize'] = self.get_value_param(options, section, 'markersize', 25, 'int')
        options_out['color'] = self.get_value_param(options, section, 'color', 'black', 'str')
        if options_out['groupBy'] is not None:
            options_out['edgecolor'] = self.get_value_param(options, section, 'edgecolor', 'gray', 'str')
            options_out['linewidth'] = self.get_value_param(options, section, 'linewidth', 1.5, 'float')
        else:
            options_out['edgecolor'] = self.get_value_param(options, section, 'edgecolor', None, 'str')
            options_out['linewidth'] = self.get_value_param(options, section, 'linewidth', None, 'float')



        return options_out

    def get_value(self, options, section, key):
        value = None
        if options.has_option(section, key):
            value = options[section][key]
        return value

    def get_value_param(self, options, section, key, default, type):
        value = self.get_value(options, section, key)
        print(value)
        if value is None:
            return default
        if type == 'str':
            return value
        if type == 'file':
            if not os.path.exists(value.strip()):
                return default
            else:
                return value.strip()
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
