import MDBPlotDefaults as defaults
from MDBFile import MDBFile
import numpy as np
import pandas as pd
import os
from scipy import stats
import math
import COMMON.common_functions as cfs


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

        # params = {
        #     'N': 0,
        #     'SLOPE': 11,
        #     'OFFSET': 12,
        #     'XAVG': 13,
        #     'YAVG': 14,
        #     'DETER(r2)': 10,
        #     'RMSE': 6,
        #     'CPRMSE': 15,
        #     'BIAS': 9,
        #     'PCC(r)': 3,
        #     'RPD': 7,
        #     'APD': 8,
        #     'MAE': 16
        # }

        # self.df_valid_stats = None

        # self.units = r'sr$^-$$^1$'
        # self.log_scale = False

        # self.satellite = 'S3'
        # self.platform = 'AB'

    def plot_from_options(self, options):
        plot_list = list(options.sections())
        for plot in plot_list:
            options_out = self.get_options(options, plot)
            #print(options_out)
            if options_out['apply']:
                self.plot_from_options_impl(options_out)

    def plot_from_options_impl(self, options_out):
        if options_out['type'] == 'scatterplot':
            if not options_out['selectByWavelength'] and options_out['selectBy'] is None:
                self.set_data_scatterplot(options_out['groupBy'], None,None,None)
                self.plot_scatter_plot(options_out)

        if options_out['type'] == 'statstable_wl':
            params = options_out['params']#self.valid_stats.keys()
            flags = ['GLOBAL']
            flag_name = options_out['selectBy']
            if flag_name is not None and options_out['selectType']=='flag':
                flag_list = self.get_flag_list(options_out['selectValues'],options_out[flag_name]['flag_values'],options_out[flag_name]['flag_meanings'])
                flags = flags + flag_list
            table, indices = self.start_table_wl(flags,params,options_out['wl_values'])
            #global stats, it's always done
            self.set_data_scatterplot(None, None, None, None)
            self.compute_statistics()
            table = self.assign_stats_table_wl(table,indices,params,'GLOBAL',None)
            for wl in options_out['wl_values']:
                self.set_data_scatterplot(None,None,None,wl)
                self.compute_statistics()
                table = self.assign_stats_table_wl(table, indices, params, 'GLOBAL', wl)
            #results by flag
            if flag_name is not None:
                for idx in range(len(flag_list)):
                    flag = flag_list[idx]
                    flag_value = options_out['selectValues'][idx]
                    self.set_data_scatterplot(None, flag_name, flag_value, None)
                    self.compute_statistics()
                    table = self.assign_stats_table_wl(table, indices, params, flag, None)
                    for wl in options_out['wl_values']:
                        self.set_data_scatterplot(None, flag_name, flag_value, wl)
                        self.compute_statistics()
                        table = self.assign_stats_table_wl(table, indices, params,flag, wl)

            print(table)

    #statistics are computed in a previous step
    def assign_stats_table_wl(self,table,indices,params,flag,wl):
        if wl is None:
            col_name = 'ALL'
        else:
            col_name = f'{wl:.2f}'.replace('.','_')
        for param in params:
            if param in self.valid_stats:
                value = self.valid_stats[param]
                index_row = indices[flag][param]
                table.iloc[index_row].at[col_name] = value
        return table

    def start_table_wl(self,flags,params,wl_values):
        wl_list = [f'{x:.2f}'.replace('.','_') for x in wl_values]
        nrows = len(flags)*len(params)
        indices = {}
        index = 0
        col_names = ['FLAG','PARAM','ALL'] + wl_list
        table = pd.DataFrame(columns=col_names,index=range(nrows))
        for flag in flags:
            indices[flag] = {}
            for param in params:
                table.iloc[index].at['FLAG'] = flag
                table.iloc[index].at['PARAM'] = param
                indices[flag][param] = index
                index = index +1
        return table,indices

    def set_data_scatterplot(self, groupBy, selectBy, valSelect, wl_value):
        rrs_ins = np.array(self.mrfile.nc.variables['mu_ins_rrs'])
        rrs_sat = np.array(self.mrfile.nc.variables['mu_sat_rrs'])
        id_all = np.array(self.mrfile.nc.variables['mu_satellite_id'])
        mu_valid = np.array(self.mrfile.nc.variables['mu_valid'])

        valid_all = self.get_array_all_from_arraymu(id_all, mu_valid)

        if wl_value is not None:
            wl_array = np.array(self.mrfile.nc.variables['mu_wavelength'])
            valid_all[wl_array!=wl_value] = 0

        if selectBy is not None and valSelect is not None:
            select_array = np.array(self.mrfile.nc.variables[selectBy])
            if len(select_array)==len(mu_valid):
                select_array = self.get_array_all_from_arraymu(id_all, select_array)
            valid_all[select_array!=valSelect] = 0

        self.xdata = rrs_ins[valid_all == 1]
        self.ydata = rrs_sat[valid_all == 1]

        if groupBy is not None:
            group_array = np.array(self.mrfile.nc.variables[groupBy])
            if len(group_array) == len(mu_valid):
                group_array = self.get_array_all_from_arraymu(id_all, group_array)
            self.groupdata = group_array[valid_all == 1]

        print('setting group data', groupBy)

    def get_array_all_from_arraymu(self, id_all_array, mu_array):
        array_out = np.zeros(id_all_array.shape, dtype=mu_array.dtype)
        for id in range(len(mu_array)):
            array_out[id_all_array == id] = mu_array[id]
        return array_out

    # MAIN FUNCTION TO PLOT SCATTERPLOT
    def plot_scatter_plot(self, options):

        if options['include_stats']:
            self.compute_statistics()

        ngroup = 1
        str_legend = []
        groupValues = options['groupValues']
        if len(self.groupdata) > 0 and groupValues is not None:
            ngroup = len(groupValues)
            if ngroup > 1 and options['legend']:
                if options['groupType'] == 'float':
                    for g in groupValues:
                        str_legend.append(f'{g:.2f}')
                if options['groupType'] == 'flag':
                    flag_name = options['groupBy']
                    str_legend = self.get_flag_list(groupValues,options[flag_name]['flag_values'],options[flag_name]['flag_meanings'])

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
            for g in groupValues:
                if options['groupType'] == 'flag':
                    color = defaults.get_color_flag(g)
                else:
                    color = defaults.get_color_ref(g)
                xhere = self.xdata[self.groupdata == g]
                yhere = self.ydata[self.groupdata == g]
                # None, None, color, 'gray', 1.5
                plot.plot_data(xhere, yhere, options['marker'], options['markersize'], color, options['edgecolor'],
                               options['linewidth'])
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
                    plot.plot_data(xhere, yhere, options['marker'], options['markersize'], z, options['edgecolor'],
                                   options['linewidth'])
                    plot.set_cmap('jet')
                except:
                    plot.plot_data(xhere, yhere, options['marker'], options['markersize'], options['color'],
                                   options['edgecolor'], options['linewidth'])

            else:
                plot.plot_data(xhere, yhere, options['marker'], options['markersize'], options['color'],
                               options['edgecolor'], options['linewidth'])

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

        if options['include_stats']:
            if options['type_scatterplot'] == 'rrs':
                str0 = 'N={:d}\nRMSD={:,.1e} UNITS\nAPD={:,.0f}%\nRPD={:,.0f}%\n$r^2$={:,.2f}\nbias={:,.1e} UNITS' \
                    .format(self.valid_stats['N'],
                            self.valid_stats['RMSD'],
                            self.valid_stats['RPD'],
                            self.valid_stats['APD'],
                            self.valid_stats['DETER(r2)'],
                            self.valid_stats['BIAS'])
                str0 = str0.replace('UNITS', options['units'])
                plot.plot_text(0.05, 0.70, str0)
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
        if not options['log_scale'] and options['regression_line']:
            plot.plot_regress_line(self.xregress, self.yregress, 'black')

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
        if options_out['type'] == 'statstable_wl':
            options_out = self.get_select_options(options, section, options_out)
            options_out = self.get_options_statstable_wl(options,section,options_out)

        return options_out

    def get_group_options(self, options, section, options_out):
        options_out['groupBy'] = self.get_value_param(options, section, 'groupBy', None, 'str')
        if options_out['groupBy'] is not None:
            var_name = options_out['groupBy']
            if var_name in self.mrfile.nc.variables:
                options_out['groupValues'] = None
                options_out['groupType'] = 'float'
                if var_name.startswith('flag'):
                    flag_values = self.mrfile.nc.variables[var_name].flag_values
                    flag_meanings_list = self.mrfile.nc.variables[var_name].flag_meanings.split(' ')
                    flag_meanings = [x.strip() for x in flag_meanings_list]
                    options_out[var_name] = {
                        'flag_values': flag_values,
                        'flag_meanings': flag_meanings
                    }
                    options_out['groupValues'] = self.get_value_param(options, section, 'groupValues', flag_values,'intlist')
                    options_out['groupType'] = 'flag'
                if options_out['groupType'] == 'float':
                    group_values = list(np.unique(np.array(self.mrfile.nc.variables[var_name])))
                    options_out['groupValues'] = self.get_value_param(options, section, 'groupValues', group_values,
                                                                      'floatlist')
        return options_out

    def get_select_options(self, options, section, options_out):

        options_out['selectByWavelength'] = self.get_value_param(options,section,'selectByWavelength',False,'boolean')
        wl_values = self.get_value_param(options,section,'wlvalues',None,'floatlist')
        if wl_values is None and options_out['selectByWavelength']:
            wl_values = list(np.array(self.mrfile.nc.variables['satellite_bands']))
        options_out['wl_values'] = wl_values

        options_out['selectBy'] = self.get_value_param(options, section, 'selectBy', None, 'str')
        options_out['selectValues'] = None
        if options_out['selectBy'] is not None:
            var_name = options_out['selectBy']
            options_out['selectType'] = 'float'
            if var_name.startswith('flag'):
                flag_values = self.mrfile.nc.variables[var_name].flag_values
                flag_meanings_list = self.mrfile.nc.variables[var_name].flag_meanings.split(' ')
                flag_meanings = [x.strip() for x in flag_meanings_list]
                options_out[var_name] = {
                    'flag_values': flag_values,
                    'flag_meanings': flag_meanings
                }
                print('================================',flag_values)
                options_out['selectValues'] = self.get_value_param(options, section, 'selectValues', flag_values,'intlist')
                print('================================', options_out['selectValues'])
                options_out['selectType'] = 'flag'
            if options_out['selectType'] == 'float':
                group_values = list(np.unique(np.array(self.mrfile.nc.variables[var_name])))
                options_out['selectValues'] = self.get_value_param(options, section, 'selectValues', group_values,'floatlist')



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
        unitsdefault = None
        if options_out['type_scatterplot'] == 'rrs':
            unitsdefault = r'sr$^-$$^1$'
            sfdefault = 1000
        xlabeldefault = defaults.xlabel_default
        ylabeldefault = defaults.ylabel_default
        options_out['scale_factor'] = self.get_value_param(options, section, 'scale_factor', sfdefault, 'float')
        options_out['xlabel'] = self.get_value_param(options, section, 'xlabel', xlabeldefault, 'float')
        options_out['ylabel'] = self.get_value_param(options, section, 'ylabel', ylabeldefault, 'float')
        options_out['units'] = self.get_value_param(options, section, 'units', unitsdefault, 'str')
        options_out['identity_line'] = self.get_value_param(options, section, 'identity_line', True, 'boolean')
        options_out['regression_line'] = self.get_value_param(options, section, 'regression_line', True, 'boolean')
        # marker, markersize, color, edgecolor, linewidth
        # None, 25, 'black', None, None
        options_out['marker'] = self.get_value_param(options, section, 'marker', None, 'str')
        options_out['markersize'] = self.get_value_param(options, section, 'markersize', 25, 'int')
        options_out['color'] = self.get_value_param(options, section, 'color', 'black', 'str')

        edgeColorDefault = None
        lineWidthDefault = None
        if options_out['groupBy'] is not None:
            if options_out['groupType']=='rrs':
                edgeColorDefault = 'gray'
                lineWidthDefault = 1.5
            if options_out['groupType']=='flag':
                edgeColorDefault = 'black'
                lineWidthDefault = 0.25

        options_out['edgecolor'] = self.get_value_param(options, section, 'edgecolor', edgeColorDefault, 'str')
        options_out['linewidth'] = self.get_value_param(options, section, 'linewidth', lineWidthDefault, 'float')

        return options_out

    def get_options_statstable_wl(self,options,section,options_out):
        options_out['selectByWavelength'] = True ##option to be always true
        if options_out['wl_values'] is None:
            options_out['wl_values'] = list(np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])))

        options_out['params'] = self.get_value_param(options,section,'params',self.valid_stats.keys(),'strlist')

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

    def get_flag_list(self,values, allValues, allFlags):
        flag_list = []
        for val in values:
            indext = np.where(allValues==val)
            index = indext[0]
            if len(index)==1:
                indexf = index[0]
                flag_list.append(allFlags[indexf])
        return flag_list

    def compute_statistics(self):

        self.valid_stats['N'] = len(self.xdata)

        # Generated linear fit
        xdatal = []
        ydatal = []
        maxxy = 0
        for x, y in zip(self.xdata, self.ydata):
            xdatal.append(x)
            ydatal.append(y)
            if x > maxxy:
                maxxy = x
            if y > maxxy:
                maxxy = y

        slope, intercept, r_value, p_value, std_err = stats.linregress(xdatal, ydatal)

        self.xregress = []
        self.yregress = []
        self.xregress.append(0)
        self.yregress.append(intercept)
        for x in xdatal:
            yr = (x * slope) + intercept
            self.yregress.append(yr)
            self.xregress.append(x)
        yrmax = ((maxxy + 1) * slope) + intercept
        self.xregress.append(maxxy + 1)
        self.yregress.append(yrmax)

        self.valid_stats['slope'] = slope
        self.valid_stats['intercept'] = intercept
        self.valid_stats['PCC(r)'] = r_value
        self.valid_stats['p_value'] = p_value
        self.valid_stats['std_err'] = std_err

        ref_obs = np.asarray(self.xdata, dtype=np.float)
        sat_obs = np.asarray(self.ydata, dtype=np.float)

        # sat_avg = np.mean(sat_obs)
        # ref_avg = np.mean(ref_obs)
        # sat_minus_avg2 = 0
        # ref_minus_avg2 = 0
        # sat_ref = 0
        # cval = 0
        # for idx in range(len(sat_obs)):
        #     vsat = sat_obs[idx]
        #     vref = ref_obs[idx]
        #     val_sat = (vsat - sat_avg) * (vsat - sat_avg)
        #     sat_minus_avg2 = sat_minus_avg2 + val_sat
        #     val_ref = (vref - ref_avg) * (vref - ref_avg)
        #     ref_minus_avg2 = ref_minus_avg2 + val_ref
        #     val_here = (vsat - sat_avg) * (vref - ref_avg)
        #     sat_ref = sat_ref + val_here
        #     cval_here = math.pow(((vsat - sat_avg) - (vref - ref_avg)), 2)
        #     cval = cval + cval_here

        # num1 = math.pow((sat_minus_avg2 - ref_minus_avg2), 2)
        # num2 = math.pow(sat_ref, 2) * 4
        # num3 = math.pow((num1 + num2), 0.5)
        # num = sat_minus_avg2 - ref_minus_avg2 + num3
        # dem = 2 * sat_ref

        # results = regress2(ref_obs, sat_obs, _method_type_2="major axis")
        # self.valid_stats['slope_typeII'] = results['slope']
        # self.valid_stats['offset_typeII'] = results['intercept']

        self.valid_stats['RMSD'] = cfs.rmse(sat_obs, ref_obs)

        ref_mean = np.mean(ref_obs)
        sat_mean = np.mean(sat_obs)
        self.valid_stats['XAVG'] = ref_mean
        self.valid_stats['YAVG'] = sat_mean

        # CPRMSE
        xdiff = ref_obs - ref_mean
        ydiff = sat_obs - sat_mean
        cprmse = cfs.rmse(ydiff, xdiff)
        self.valid_stats['CPRMSE'] = cprmse

        # the mean of relative (signed) percent differences
        rel_diff = 100 * ((sat_obs - ref_obs) / ref_obs)
        self.valid_stats['RPD'] = np.mean(rel_diff)

        #  the mean of absolute (unsigned) percent differences
        self.valid_stats['APD'] = np.mean(np.abs(rel_diff))

        bias = np.mean(sat_obs - ref_obs)
        self.valid_stats['BIAS'] = bias

        mae = np.mean(np.abs(sat_obs - ref_obs))
        self.valid_stats['MAE'] = mae

        self.valid_stats['DETER(r2)'] = r_value * r_value
