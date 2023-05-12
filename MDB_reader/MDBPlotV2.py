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

        self.global_stats_ins_spectra = {}

    def plot_from_options(self, options):
        plot_list = list(options.sections())
        for plot in plot_list:
            options_out = self.get_options(options, plot)
            # print(options_out)
            if options_out['apply']:
                self.plot_from_options_impl(options_out)

    def plot_from_options_impl(self, options_out):
        if options_out['type'] == 'scatterplot':
            if options_out['selectByWavelength']:  # one scatterplot for wavelenght
                file_out_base = options_out['file_out']
                title_base = options_out['title']
                if options_out['selectBy'] is None:  # scatter plot global by wavelengh
                    if options_out['multiple_plot'] is None:
                        for wl in options_out['wl_values']:
                            self.set_data_scatterplot(options_out['groupBy'], None, None, wl)
                            options_out['file_out'] = self.get_file_out_name(file_out_base, wl, None)
                            options_out['title'] = self.get_title(title_base, wl, None, None)
                            self.plot_scatter_plot(options_out, None, -1, wl)
                    else:
                        rc = options_out['multiple_plot'].split(',')
                        nrow = int(rc[0].strip())
                        ncol = int(rc[1].strip())
                        ntot = nrow * ncol
                        index = 0
                        from PlotScatter import PlotScatter
                        plot_here = PlotScatter()
                        plot_here.xtitle_options['fontsize'] = options_out['fontsize']
                        plot_here.ytitle_options['fontsize'] = options_out['fontsize']
                        plot_here.plot_text_options['fontsize'] = options_out['fontsize']
                        plot_here.start_multiple_plot(nrow,ncol)
                        print(f'[INFO] Starting multiple plot with {nrow} rows and {ncol} cols')
                        print(options_out['wl_values'])
                        for wl in options_out['wl_values']:
                            self.set_data_scatterplot(options_out['groupBy'], None, None, wl)
                            options_out['title'] = self.get_title(title_base, wl, None, None)
                            options_out['file_out'] = None
                            self.plot_scatter_plot(options_out,plot_here,index,wl)
                            index = index +1
                        for index_blank in range(index,ntot):
                            plot_here.plot_blanck(index_blank)
                        if options_out['legend']:
                            str_legend = self.get_str_legend(options_out)
                            if len(str_legend) > 0:
                                plot_here.set_global_legend(str_legend)
                        file_out = self.get_file_out_name(file_out_base,None,None)
                        plot_here.save_fig(file_out)
                        plot_here.close_plot()


            if not options_out['selectByWavelength'] and options_out['selectBy'] is None:
                self.set_data_scatterplot(options_out['groupBy'], None, None, None)
                self.plot_scatter_plot(options_out,None,-1,-1)

        if options_out['type'] == 'statstable_wl':
            self.create_table_stats_wl(options_out)

        if options_out['type'] == 'statstable':
            self.create_table_stats(options_out, None)
            for wl in options_out['wl_values']:
                print('wl-->',wl)
                self.create_table_stats(options_out, wl)

        if options_out['type'] == 'statswlplot':
            self.plot_statistics_bywl(options_out)

        if options_out['type'] == 'spectraplot':
            self.plot_spectra_plot(options_out)

    def get_wl_str_from_wl(self, wl_value):
        wl_sat = np.array(self.mrfile.nc.variables['satellite_bands'])
        index_sat = np.argmin(np.abs(wl_sat - wl_value))
        wl_sat_value = wl_sat[index_sat]
        wl_sat_value_str = f'{wl_sat_value:.2f}'
        if wl_sat_value_str.endswith('.00'):
            return wl_sat_value_str[:-3]
        else:
            return wl_sat_value_str

    def get_file_out_name(self, file_out, wl, flag):
        if file_out is None:
            return None
        if wl is None and flag is None:
            return file_out
        if wl is not None:
            wls = self.get_wl_str_from_wl(wl)
            wls = wls.replace('.', '_')
        if wl is not None and flag is None:
            file_out = file_out[:-4] + '_' + wls + file_out[-4:]
        if wl is None and flag is not None:
            file_out = file_out[:-4] + '_' + flag + file_out[-4:]
        if wl is not None and flag is not None:
            file_out = file_out[:-4] + '_' + flag + '_' + wls + file_out[-4:]
        return file_out

    def get_file_out_flag_param(self, file_out, flag, param):
        if file_out is None:
            return None
        if param is None and flag is None:
            return file_out
        file_res = file_out[:-4]
        if flag is not None:
            file_res = file_res + '_' + flag
        if param is not None:
            file_res = file_res + '_' + param
        file_res = file_res + file_out[-4:]
        return file_res

    def get_title(self, title, wl, flag, param):
        if title is None:
            return None
        if wl is None and flag is None and param is None:
            return title
        if wl is not None:
            wls = self.get_wl_str_from_wl(wl)
            title = title.replace('$WL$', wls)
        if flag is not None:
            title = title.replace('$FLAG$', wls)

        if param is not None:
            title = title.replace('$PARAM$', param)

        return title

    def start_table_wl(self, flags, params, wl_values):
        # wl_list = [f'{x:.2f}'.replace('.', '_') for x in wl_values]
        wl_list = [self.get_wl_str_from_wl(x) for x in wl_values]
        nrows = len(flags) * len(params)
        indices = {}
        index = 0
        col_names = ['FLAG', 'PARAM', 'ALL'] + wl_list
        table = pd.DataFrame(columns=col_names, index=range(nrows))
        for flag in flags:
            indices[flag] = {}
            for param in params:
                table.iloc[index].at['FLAG'] = flag
                table.iloc[index].at['PARAM'] = param
                indices[flag][param] = index
                index = index + 1
        return table, indices

        # statistics are computed in a previous step

    def assign_stats_table_wl(self, table, indices, params, flag, wl):
        if wl is None:
            col_name = 'ALL'
        else:
            # col_name = f'{wl:.2f}'.replace('.', '_')
            col_name = self.get_wl_str_from_wl(wl)
        for param in params:
            if param in self.valid_stats:
                value = self.valid_stats[param]
                index_row = indices[flag][param]
                table.iloc[index_row].at[col_name] = value
        return table

    def start_table(self, flags, params):
        col_names = ['PARAM'] + flags
        nrows = len(params)
        indices = {}
        index = 0
        table = pd.DataFrame(columns=col_names, index=range(nrows))
        for param in params:
            table.iloc[index].at['PARAM'] = param
            indices[param] = index
            index = index + 1
        return table, indices

        # statistics are computed in a previous step

    def assign_table(self, table, indices, params, col_name):
        for param in params:
            if param in self.valid_stats:
                value = self.valid_stats[param]
                index_row = indices[param]
                table.iloc[index_row].at[col_name] = value
        return table

    def set_data_scatterplot(self, groupBy, selectBy, valSelect, wl_value):
        rrs_ins = np.array(self.mrfile.nc.variables['mu_ins_rrs'])
        rrs_sat = np.array(self.mrfile.nc.variables['mu_sat_rrs'])
        id_all = np.array(self.mrfile.nc.variables['mu_satellite_id'])
        mu_valid = np.array(self.mrfile.nc.variables['mu_valid'])
        #print('MVALID: ',np.sum(mu_valid))

        valid_all = self.get_array_all_from_arraymu(id_all, mu_valid)

        if wl_value is not None:
            wl_array = np.array(self.mrfile.nc.variables['mu_wavelength'])
            valid_all[wl_array != wl_value] = 0
        wl_array = np.array(self.mrfile.nc.variables['mu_wavelength'])
        valid_all[wl_array==832.8] = 0
        valid_all[wl_array==1613.7] =0
        valid_all[wl_array == 2202.4] = 0

        if selectBy is not None and valSelect is not None:
            select_array = np.array(self.mrfile.nc.variables[selectBy])
            if len(select_array) == len(mu_valid):
                select_array = self.get_array_all_from_arraymu(id_all, select_array)
            valid_all[select_array != valSelect] = 0

        self.xdata = rrs_ins[valid_all == 1]
        self.ydata = rrs_sat[valid_all == 1]

        if groupBy is not None:
            group_array = np.array(self.mrfile.nc.variables[groupBy])
            if len(group_array) == len(mu_valid):
                group_array = self.get_array_all_from_arraymu(id_all, group_array)
            self.groupdata = group_array[valid_all == 1]

        # print('setting group data', groupBy)

    def get_array_all_from_arraymu(self, id_all_array, mu_array):
        array_out = np.zeros(id_all_array.shape, dtype=mu_array.dtype)
        for id in range(len(mu_array)):
            array_out[id_all_array == id] = mu_array[id]
        return array_out

    def get_str_legend(self,options):
        str_legend = []
        ngroup = 1
        groupValues = options['groupValues']
        if len(self.groupdata) > 0 and groupValues is not None:
            ngroup = len(groupValues)
            if ngroup > 1:
                if options['groupType'] == 'float':
                    for g in groupValues:
                        str_legend.append(f'{g:.2f}')
                if options['groupType'] == 'flag':
                    flag_name = options['groupBy']
                    str_legend = self.get_flag_list(groupValues, options[flag_name]['flag_values'],
                                                    options[flag_name]['flag_meanings'])
        return str_legend

    # MAIN FUNCTION TO PLOT SCATTERPLOT
    def plot_scatter_plot(self, options, plot, index, wl):

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
                    str_legend = self.get_flag_list(groupValues, options[flag_name]['flag_values'],
                                                    options[flag_name]['flag_meanings'])


        from scipy.stats import gaussian_kde
        if plot is None and index==-1:
            from PlotScatter import PlotScatter
            plot = PlotScatter()
            plot.close_plot()
            plot.start_plot()
        if plot is not None and index>=0:
            plot.set_axhere_index(index)

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
                    plot.set_cmap('jet')
                    plot.plot_data(xhere, yhere, options['marker'], options['markersize'], z, options['edgecolor'],
                                   options['linewidth'])

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
        if index==-1 or plot.index_row==plot.nrow-1:
            plot.set_xaxis_title(options['xlabel'])
            if options['ticks'] is not None:
                plot.set_ticks(options['ticks'], options['fontsize'])
                added_ticks = True
        if index==-1 or plot.index_col==0:
            plot.set_yaxis_title(options['ylabel'])
            if options['ticks'] is not None:
                plot.set_ticks(options['ticks'], options['fontsize'])
        if index>=0 and plot.index_col>0:
            plot.set_yticks_labels_off(options['ticks'])
        if index>=0 and plot.index_row<(plot.nrow-1):
            plot.set_xticks_labels_off(options['ticks'])
        plot.set_equal_apect()

        if options['legend'] and len(str_legend) > 0 and index==-1:
            plot.set_legend(str_legend)
        if options['identity_line']:
            plot.plot_identity_line()

        if options['include_stats']:
            if options['type_scatterplot'] == 'rrs':
                if index==-1:
                    str0 = 'N={:d}\nRMSD={:,.1e} UNITS\nRPD={:,.0f}%\nAPD={:,.0f}%\n$r^2$={:,.2f}\nbias={:,.1e} UNITS' \
                        .format(self.valid_stats['N'],
                                self.valid_stats['RMSD'],
                                self.valid_stats['RPD'],
                                self.valid_stats['APD'],
                                self.valid_stats['DETER(r2)'],
                                self.valid_stats['BIAS'])
                    str0 = str0.replace('UNITS', options['units'])
                    plot.plot_text(0.05, 0.70, str0)
                if index>=0:
                    bias = self.valid_stats['BIAS']
                    r2 = self.valid_stats['DETER(r2)']
                    rmsd = self.valid_stats['RMSD']
                    str0 = f'{wl:.2f} nm\nbias={bias:.1e}\nr\u00b2={r2:.2f}'
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



    def plot_statistics_bywl(self, options_out):
        params = options_out['params']  # self.valid_stats.keys()
        flags = ['GLOBAL']
        flag_list = []
        flag_name = options_out['selectBy']
        if flag_name is not None and options_out['selectType'] == 'flag':
            flag_list = self.get_flag_list(options_out['selectValues'], options_out[flag_name]['flag_values'],
                                           options_out[flag_name]['flag_meanings'])
            flags = flags + flag_list
        wl_list = options_out['wl_values']
        params = options_out['params']
        xdata_plot = [float(self.get_wl_str_from_wl(x)) for x in wl_list]
        wl_col = [self.get_wl_str_from_wl(x) for x in wl_list]

        table, indices = self.start_table_wl(flags, params, options_out['wl_values'])
        legend = []
        if len(flag_list) == 0:
            for wl in options_out['wl_values']:
                self.set_data_scatterplot(None, None, None, wl)
                self.compute_statistics()
                table = self.assign_stats_table_wl(table, indices, params, 'GLOBAL', wl)
        else:
            for idx in range(len(flag_list)):
                flag = flag_list[idx]
                legend.append(flag)
                flag_value = options_out['selectValues'][idx]
                for wl in options_out['wl_values']:
                    self.set_data_scatterplot(None, flag_name, flag_value, wl)
                    self.compute_statistics()
                    table = self.assign_stats_table_wl(table, indices, params, flag, wl)

        from PlotSpectra import PlotSpectra

        # GLOBAL
        for param in params:
            plot = PlotSpectra()
            plot.close_plot()
            plot.start_plot()
            if len(flag_list) == 0:  ##solo global
                irow = indices['GLOBAL'][param]
                ydata_plot = np.array(table[wl_col].iloc[irow])
                self.plot_spectra_line_impl(plot, xdata_plot, ydata_plot, 0, options_out)
                file_out = self.get_file_out_flag_param(options_out['file_out'], 'GLOBAL', param)
            else:
                for idx in range(len(flag_list)):
                    flag = flag_list[idx]
                    irow = indices[flag][param]
                    ydata_plot = np.array(table[wl_col].iloc[irow])
                    self.plot_spectra_line_impl(plot, xdata_plot, ydata_plot, idx, options_out)
                file_out = self.get_file_out_flag_param(options_out['file_out'], flag_name, param)

            plot.set_xaxis_title(options_out['xlabel'])
            plot.set_yaxis_title(param)
            plot.set_xticks(xdata_plot, wl_col, 90, 8)
            plot.set_grid()

            plot.set_title(self.get_title(options_out['title'], None, None, param))
            if len(legend) > 0:
                plot.set_legend(legend)
            plot.set_tigth_layout()

            if file_out is not None:
                plot.save_fig(file_out)
            plot.close_plot()

    def plot_spectra_line_impl(self, plot, xdata_plot, ydata_plot, index, options_out):
        plot.xdata = xdata_plot
        line_color = options_out['line_color']
        lc = line_color[0]
        if 0 <= index < len(line_color):
            lc = line_color[index]

        line_type = options_out['line_type']
        lt = line_type[0]
        if 0 <= index < len(line_type):
            lt = line_type[index]

        line_width = options_out['line_width']
        lw = line_width[0]
        if 0 <= index < len(line_width):
            lw = line_width[index]

        marker = options_out['marker']
        m = marker[0]
        if 0 <= index < len(marker):
            m = marker[index]

        markersize = options_out['marker_size']
        ms = markersize[0]
        if 0 <= index < len(markersize):
            ms = markersize[index]

        plot.plot_single_line(ydata_plot, lc, lt, lw, m, ms)

    def plot_spectra_plot(self, options_out):

        type_rrs = options_out['type_rrs']

        if type_rrs.startswith('mu_'):  # mu_comparison, mu_sat, mu_ins
            print(options_out)
            self.plot_spectra_mu_comparison(options_out)
            return
        if type_rrs.startswith('flag_') and type_rrs!='flag_sat_insitu':
            self.plot_spectra_comparison_stats(options_out)
            return
        if type_rrs.startswith('flag_sat_insitu'):
            self.plot_comparison_stats_sat_insitu(options_out)
            return

        ##GETTING DATA
        spectra_data = options_out['plot_spectra']
        if type_rrs == 'ins':
            wavelength = self.mrfile.get_insitu_wl()
            if options_out['plot_stats'] or spectra_data[0].lower() == 'all':
                spectra, stats = self.mrfile.get_all_insitu_valid_spectra(options_out['scale_factor'])

        from PlotSpectra import PlotSpectra
        if not options_out['plot_stats']:
            stats = None

        make_by_mu = False
        muoptions = ['MU_all', 'MU_valid', 'MU_invalid', 'MU_selected']
        for m in muoptions:
            if m in options_out['plot_spectra']:
                make_by_mu = True
        if not make_by_mu:
            if options_out['plot_spectra'] is None or 'All' not in options_out['plot_spectra']:
                spectra = None
            pspectra = PlotSpectra()
            pspectra.plot_multiple_spectra(wavelength, spectra, stats, options_out['wl_min'], options_out['wl_max'])
            pspectra.set_xaxis_title(options_out['xlabel'])
            pspectra.set_yaxis_title(options_out['ylabel'])
            if options_out['title'] is not None:
                title_here = options_out['title']
                pspectra.set_title(title_here)
            pspectra.set_grid()
            pspectra.set_tigth_layout()
            if not options_out['file_out'] is None:
                pspectra.save_fig(options_out['file_out'])
            pspectra.close_plot()
        if make_by_mu:
            for index_mu in range(self.mrfile.n_mu_total):
                if index_mu == 0 or (index_mu % 100) == 0:
                    print(f'[INFO] Plotting spectra for MU: {index_mu}')
                if index_mu == 15:
                    self.plot_spectra_plot_mu(options_out, index_mu, wavelength, stats)

    def plot_spectra_plot_mu(self, options_out, index_mu, wavelength, stats):
        from PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        if stats is not None:
            imin, imax = pspectra.plot_multiple_spectra(wavelength, None, stats, options_out['wl_min'],
                                                        options_out['wl_max'])
        else:
            imin, imax = pspectra.get_imin_imax_from_wavelength(wavelength, options_out['wl_min'],
                                                                options_out['wl_max'])
            pspectra.xdata = wavelength[imin:imax]
        spectra_selected, spectra_valid, spectra_invalid = self.mrfile.get_mu_insitu_spectra(index_mu, options_out[
            'scale_factor'])

        if options_out['plot_spectra'][0] == 'MU_all':
            for spectrum in spectra_valid:
                hline = pspectra.plot_single_line(spectrum[imin:imax], 'black', 'solid', 1, None, 25)
            for spectrum in spectra_invalid:
                hline = pspectra.plot_single_line(spectrum[imin:imax], 'black', 'solid', 1, None, 25)
            for spectrum in spectra_selected:
                hline = pspectra.plot_single_line(spectrum[imin:imax], 'black', 'solid', 2, None, 25)

        str_legend = []
        hlines = []
        if 'MU_valid' in options_out['plot_spectra']:
            hline = None
            for spectrum in spectra_valid:
                hline = pspectra.plot_single_line(spectrum[imin:imax], 'green', 'solid', 1, None, 25)
            if hline is not None:
                hlines.append(hline)
                str_legend.append('Valid spectra')
        if 'MU_invalid' in options_out['plot_spectra']:
            hline = None
            for spectrum in spectra_invalid:
                hline = pspectra.plot_single_line(spectrum[imin:imax], 'red', 'solid', 1, None, 25)
            if hline is not None:
                hlines.append(hline)
                str_legend.append('Invalid spectra')
        if 'MU_selected' in options_out['plot_spectra']:
            hline = None
            for spectrum in spectra_selected:
                hline = pspectra.plot_single_line(spectrum[imin:imax], 'blue', 'solid', 2, None, 25)
            if hline is not None:
                hlines.append(hline)
                str_legend.append('Selected spectra')

        pspectra.set_xaxis_title(options_out['xlabel'])
        pspectra.set_yaxis_title(options_out['ylabel'])
        pspectra.set_y_range(-5,10)
        if options_out['title'] is not None:
            title_here = options_out['title'] + f' MU: {index_mu}'
            pspectra.set_title(title_here)
        pspectra.set_grid()
        if len(str_legend) > 0:
            pspectra.legend_options['bbox_to_anchor'] = (0.65, 1.0)
            pspectra.legend_options['framealpha'] = 1
            pspectra.set_legend_h(hlines, str_legend)
        pspectra.set_tigth_layout()
        if not options_out['file_out'] is None:
            file_out = options_out['file_out'][:-4]
            file_out = f'{file_out}_{index_mu}.{self.format_image}'
            if options_out['plot_spectra'][0] == 'MU_all':
                file_out = f'{file_out}_{index_mu}_MU_ALL.{self.format_image}'
            # print(file_out)
            pspectra.save_fig(file_out)
        pspectra.close_plot()

    def plot_spectra_mu_comparison(self, options_out):
        for index_mu in range(self.mrfile.n_mu_total):
            if index_mu == 0 or (index_mu % 100) == 0:
                print(f'[INFO] Plotting spectra for MU: {index_mu}')

            self.plot_spectra_mu_comparison_impl(index_mu, options_out)

    def plot_spectra_mu_comparison_impl(self, index_mu, options_out):
        wl, insitu_spectra, sat_spectra = self.mrfile.get_mu_spectra_insitu_and_sat(index_mu,
                                                                                    options_out['scale_factor'])
        if wl is None:
            return
        from PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        pspectra.xdata = wl
        wls = self.mrfile.get_sat_wl_as_strlist(wl)
        pspectra.set_xticks(wl, wls, 90, 8)
        if options_out['type_rrs'] == 'mu_comparison':
            hline1 = pspectra.plot_single_line(insitu_spectra, 'red', 'solid', 1, '.', 10)
            hline2 = pspectra.plot_single_line(sat_spectra, 'blue', 'solid', 1, '.', 10)
        if options_out['type_rrs'] == 'mu_sat':
            hline2 = pspectra.plot_single_line(sat_spectra, 'blue', 'solid', 1, '.', 10)
        if options_out['type_rrs'] == 'mu_ins':
            hline1 = pspectra.plot_single_line(insitu_spectra, 'red', 'solid', 1, '.', 10)

        pspectra.set_xaxis_title(options_out['xlabel'])
        pspectra.set_yaxis_title(options_out['ylabel'])
        if options_out['title'] is not None:
            title_here = options_out['title'] + f' MU: {index_mu}'
            pspectra.set_title(title_here)
        if options_out['type_rrs'] == 'mu_comparison':
            pspectra.set_legend_h([hline1, hline2], ['In situ Rrs', 'Satellite Rrs'])
        pspectra.set_grid()
        pspectra.set_tigth_layout()
        if not options_out['file_out'] is None:
            file_out = options_out['file_out'][:-4]
            file_out = f'{file_out}_{index_mu}.{self.format_image}'
            # print(file_out)
            pspectra.save_fig(file_out)
        pspectra.close_plot()

    def plot_comparison_stats_sat_insitu(self, options_out):
        from PlotSpectra import PlotSpectra
        wavelength, sat_stats, insitu_stats = self.mrfile.get_all_spectra_insitu_sat(options_out['scale_factor'])
        pspectra = PlotSpectra()
        imin, imax = pspectra.get_imin_imax_from_wavelength(wavelength, options_out['wl_min'], options_out['wl_max'])

        str_legend = ['insitu Rrs','satellite Rrs']
        ymin_sat, ymax_sat = pspectra.get_ymin_ymax_from_stats(sat_stats, imin, imax)
        ymin_insitu, ymax_insitu = pspectra.get_ymin_ymax_from_stats(insitu_stats, imin, imax)
        ymin = np.min([ymin_sat,ymin_insitu])
        ymax = np.max([ymax_sat,ymax_insitu])

        wl_list = wavelength[imin:imax]
        xdata_plot = [float(self.get_wl_str_from_wl(x)) for x in wl_list]
        wl_col = [self.get_wl_str_from_wl(x) for x in wl_list]
        pspectra.xdata = xdata_plot


        color = 'red'
        pspectra.stats_style['avg']['color'] = color
        pspectra.stats_style['avg']['marker'] = 'o'
        pspectra.stats_style['avg']['markersize'] = 5
        pspectra.stats_style['fill']['color'] = color
        pspectra.stats_style['fill']['framealpha'] = 0.5
        hlineinsitu = pspectra.plot_stats(insitu_stats, imin, imax)

        color = 'blue'
        pspectra.xdata = wavelength[imin:imax]
        pspectra.stats_style['avg']['color'] = color
        pspectra.stats_style['avg']['marker'] = 'o'
        pspectra.stats_style['avg']['markersize'] = 5
        pspectra.stats_style['fill']['color'] = color
        pspectra.stats_style['fill']['framealpha'] = 0.5
        hlinesat = pspectra.plot_stats(sat_stats, imin, imax)

        h_legend = [hlineinsitu,hlinesat]
        pspectra.set_xticks(xdata_plot,wl_col,90,8)
        pspectra.set_y_range(ymin, ymax)
        pspectra.set_xaxis_title(options_out['xlabel'])
        pspectra.set_yaxis_title(options_out['ylabel'])
        if options_out['title'] is not None:
            title_here = options_out['title']
            pspectra.set_title(title_here)
        pspectra.set_grid()
        pspectra.legend_options['bbox_to_anchor'] = (0.65, 1.0)
        pspectra.legend_options['framealpha'] = 1
        pspectra.set_legend_h(h_legend, str_legend)
        pspectra.set_tigth_layout()
        if not options_out['file_out'] is None:
            pspectra.save_fig(options_out['file_out'])
        pspectra.close_plot()

    def plot_spectra_comparison_stats(self, options_out):

        print(options_out)
        from PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        type_rrs = options_out['type_rrs']
        flag_values = options_out['groupValues']
        flag_name = options_out['groupBy']
        str_legend = self.get_flag_list(flag_values, options_out[flag_name]['flag_values'],
                                        options_out[flag_name]['flag_meanings'])
        h_legend = []
        if type_rrs == 'flag_ins_comparison':
            wavelength = self.mrfile.get_insitu_wl()
            imin, imax = pspectra.get_imin_imax_from_wavelength(wavelength, options_out['wl_min'],
                                                                options_out['wl_max'])

        ymin = None
        ymax = None
        for flag_value in flag_values:
            if type_rrs == 'flag_ins_comparison':
                wavelength = self.mrfile.get_insitu_wl()
                spectra_good_flag, stats_flag = self.mrfile.get_flag_insitu_valid_spectra(options_out['scale_factor'],
                                                                                          flag_name, flag_value)
                ymin_flag, ymax_flag = pspectra.get_ymin_ymax_from_stats(stats_flag, imin, imax)
                if ymin is None and ymax is None:
                    ymin = ymin_flag
                    ymax = ymax_flag
                else:
                    if ymin_flag < ymin:
                        ymin = ymin_flag
                    if ymax_flag > ymax:
                        ymax = ymax_flag
            print(ymin, ymax)

            # pspectra.plot_multiple_spectra(wavelength, None, stats_flag, options_out['wl_min'], options_out['wl_max'])
            color = defaults.get_color_flag(flag_value)
            pspectra.xdata = wavelength[imin:imax]
            pspectra.stats_style['avg']['color'] = color
            pspectra.stats_style['fill']['color'] = color
            pspectra.stats_style['fill']['framealpha'] = 0.5
            # print(pspectra.stats_style['fill']['color'])
            hline = pspectra.plot_stats(stats_flag, imin, imax)
            h_legend.append(hline)

        pspectra.set_y_range(ymin, ymax)
        pspectra.set_xaxis_title(options_out['xlabel'])
        pspectra.set_yaxis_title(options_out['ylabel'])
        if options_out['title'] is not None:
            title_here = options_out['title']
            pspectra.set_title(title_here)
        pspectra.set_grid()
        pspectra.legend_options['bbox_to_anchor'] = (0.65, 1.0)
        pspectra.legend_options['framealpha'] = 1
        pspectra.set_legend_h(h_legend, str_legend)
        pspectra.set_tigth_layout()
        if not options_out['file_out'] is None:
            pspectra.save_fig(options_out['file_out'])
        pspectra.close_plot()

    def plot_flag(self, info_flag, file_out):
        x = list(info_flag.keys())
        ylabel = 'nmacrow'
        ydata = []
        for flag in info_flag:
            ydata.append(info_flag[flag][ylabel])
        x = np.array(x)
        ydata = np.array(ydata)
        from matplotlib import pyplot as plt
        print(x)
        print(ydata)
        plt.figure()
        plt.barh(x, ydata)
        plt.grid(b=True, which='major', color='gray', linestyle='--')
        plt.xlabel('Total number of match-ups', fontsize=12)
        plt.gcf().tight_layout()
        if file_out is not None:
            plt.savefig(file_out, dpi=300)
        plt.close()

    def plot_flag_array(self, info_flag, flag, file_out):
        from matplotlib import pyplot as plt
        extent = [0, 25, 0, 25]
        cmap = plt.colormaps['jet']
        band = np.array(info_flag[flag]['parray'])

        plt.imshow(band, cmap=cmap, interpolation=None, extent=extent, vmin=30, vmax=50)

        plt.hlines(11, 11, 14, colors=['r'])
        plt.hlines(14, 11, 14, colors=['r'])
        plt.vlines(11, 11, 14, colors=['r'])
        plt.vlines(14, 11, 14, colors=['r'])
        plt.colorbar()
        plt.gcf().tight_layout()
        if file_out is not None:
            plt.savefig(file_out, dpi=300)
        plt.close()

    def create_table_stats_wl(self, options_out):
        params = options_out['params']  # self.valid_stats.keys()
        flags = ['GLOBAL']
        flag_name = options_out['selectBy']
        if flag_name is not None and options_out['selectType'] == 'flag':
            flag_list = self.get_flag_list(options_out['selectValues'], options_out[flag_name]['flag_values'],
                                           options_out[flag_name]['flag_meanings'])
            flags = flags + flag_list
        table, indices = self.start_table_wl(flags, params, options_out['wl_values'])
        # global stats, it's always done
        self.set_data_scatterplot(None, None, None, None)
        self.compute_statistics()
        table = self.assign_stats_table_wl(table, indices, params, 'GLOBAL', None)
        for wl in options_out['wl_values']:
            self.set_data_scatterplot(None, None, None, wl)
            self.compute_statistics()
            table = self.assign_stats_table_wl(table, indices, params, 'GLOBAL', wl)
        # results by flag
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
                    table = self.assign_stats_table_wl(table, indices, params, flag, wl)
        if not options_out['file_out'] is None:
            table.to_csv(options_out['file_out'], sep=';')

    def create_table_stats(self, options_out, wl):
        params = options_out['params']  # self.valid_stats.keys()
        flags = ['GLOBAL']
        flag_name = options_out['selectBy']
        flag_list = []
        if flag_name is not None and options_out['selectType'] == 'flag':
            # print(options_out['selectValues'])
            # print(options_out[flag_name]['flag_values'])
            # print(options_out[flag_name]['flag_meanings'])
            flag_list = self.get_flag_list(options_out['selectValues'], options_out[flag_name]['flag_values'],
                                           options_out[flag_name]['flag_meanings'])
            flags = flags + flag_list
        table, indices = self.start_table(flags, params)
        # global stats, it's always done
        self.set_data_scatterplot(None, None, None, wl)
        self.compute_statistics()
        table = self.assign_table(table, indices, params, 'GLOBAL')
        # stats by flag
        if len(flag_list) > 0:
            for idx in range(len(flag_list)):
                flag = flag_list[idx]
                flag_value = options_out['selectValues'][idx]
                self.set_data_scatterplot(None, flag_name, flag_value, wl)
                self.compute_statistics()
                table = self.assign_table(table, indices, params, flag)

        if not options_out['file_out'] is None:
            file_out = options_out['file_out']
            if wl is not None:
                wls = self.get_wl_str_from_wl(wl)
                wls = wls.replace('.', '_')
                file_out = file_out[:-4] + '_' + wls + '.csv'
            table.to_csv(file_out, sep=';')

    def get_options(self, options, section):
        # print(section)
        options_out = {'apply': self.get_value_param(options, section, 'apply', False, 'boolean')}
        if not options_out['apply']:
            return options_out
        options_out['type'] = self.get_value_param(options, section, 'type', None, 'str')
        if options_out['type'] is None:
            return options_out
        options_out['name'] = section
        options_out['multiple_plot'] = self.get_value_param(options,section,'multiple_plot',None,'str')
        if options_out['type'] == 'scatterplot':
            options_out = self.get_group_options(options, section, options_out)
            options_out = self.get_select_options(options, section, options_out)
            options_out = self.get_options_scatterplot(options, section, options_out)
        if options_out['type'].startswith('statstable'):
            options_out = self.get_select_options(options, section, options_out)
            options_out = self.get_options_statstable(options, section, options_out)
        if options_out['type'] == 'statswlplot':
            options_out = self.get_select_options(options, section, options_out)
            options_out = self.get_options_statswlplot(options, section, options_out)
        if options_out['type'] == 'spectraplot':
            options_out['type_rrs'] = self.get_value_param(options, section, 'type_rrs', 'ins', 'str')
            if options_out['type_rrs'].startswith('flag'):
                options_out = self.get_group_options(options, section, options_out)
            options_out = self.get_options_spectraplot(options, section, options_out)

        return options_out

    def get_group_options(self, options, section, options_out):
        options_out['groupBy'] = self.get_value_param(options, section, 'groupBy', None, 'str')
        options_out['groupValues'] = None
        if options_out['groupBy'] is not None:
            var_name = options_out['groupBy']
            if var_name in self.mrfile.nc.variables:
                options_out['groupType'] = 'float'
                if var_name.startswith('flag'):
                    flag_values = self.mrfile.nc.variables[var_name].flag_values
                    flag_meanings_list = self.mrfile.nc.variables[var_name].flag_meanings.split(' ')
                    flag_meanings = [x.strip() for x in flag_meanings_list]
                    options_out[var_name] = {
                        'flag_values': flag_values,
                        'flag_meanings': flag_meanings
                    }
                    #options_out['groupValues'] = self.get_value_param(options, section, 'groupValues', flag_values,'intlist')
                    flag_list_config = self.get_value_param(options, section, 'groupValues', None, 'strlist')
                    if flag_list_config is None:
                        print(flag_values)
                        options_out['groupValues'] = flag_values
                    else:
                        flag_values_config = []
                        for flag_config in flag_list_config:
                            try:
                                iflag = flag_meanings.index(flag_config.strip())
                                flag_values_config.append(flag_values[iflag])
                            except:
                                print(f'[WARNING] Flag {flag_config.strip()} is not in the list')
                        options_out['groupValues'] = flag_values_config
                    options_out['groupType'] = 'flag'
                if options_out['groupType'] == 'float':
                    group_values = list(np.unique(np.array(self.mrfile.nc.variables[var_name])))
                    options_out['groupValues'] = self.get_value_param(options, section, 'groupValues', group_values,
                                                                      'floatlist')
        return options_out

    def get_select_options(self, options, section, options_out):

        options_out['selectByWavelength'] = self.get_value_param(options, section, 'selectByWavelength', False,
                                                                 'boolean')
        wl_values = self.get_value_param(options, section, 'wlvalues', None, 'floatlist')
        if wl_values is None and options_out['selectByWavelength']:
            wl_values = list(np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])))
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
                flag_list_config = self.get_value_param(options, section, 'selectValues', None, 'strlist')
                if flag_list_config is None:
                    print('y aqui deberia llegar')
                    options_out['selectValues'] = flag_values
                else:
                    flag_values_config = []
                    for flag_config in flag_list_config:
                        try:
                            iflag = flag_meanings.index(flag_config.strip())
                            flag_values_config.append(flag_values[iflag])
                        except:
                            print(f'[WARNING] Flag {flag_config.strip()} is not in the list')
                    options_out['selectValues'] = flag_values_config
                #options_out['selectValues'] = self.get_value_param(options, section, 'selectValues', flag_values,'intlist')

                options_out['selectType'] = 'flag'
            if options_out['selectType'] == 'float':
                group_values = list(np.unique(np.array(self.mrfile.nc.variables[var_name])))
                options_out['selectValues'] = self.get_value_param(options, section, 'selectValues', group_values,
                                                                   'floatlist')

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
        options_out['ticks'] = self.get_value_param(options, section, 'ticks', None, 'floatlist')
        options_out['fontsize'] = self.get_value_param(options, section, 'fontsize', 12, 'float')
        sfdefault = None
        unitsdefault = None
        if options_out['type_scatterplot'] == 'rrs':
            unitsdefault = r'sr$^-$$^1$'
            sfdefault = 1000
        xlabeldefault = defaults.xlabel_default
        ylabeldefault = defaults.ylabel_default
        options_out['scale_factor'] = self.get_value_param(options, section, 'scale_factor', sfdefault, 'float')
        options_out['xlabel'] = self.get_value_param(options, section, 'xlabel', xlabeldefault, 'str')
        options_out['ylabel'] = self.get_value_param(options, section, 'ylabel', ylabeldefault, 'str')
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
            if options_out['groupType'] == 'rrs':
                edgeColorDefault = 'gray'
                lineWidthDefault = 1.5
            if options_out['groupType'] == 'flag':
                edgeColorDefault = 'black'
                lineWidthDefault = 0.25

        options_out['edgecolor'] = self.get_value_param(options, section, 'edgecolor', edgeColorDefault, 'str')
        options_out['linewidth'] = self.get_value_param(options, section, 'linewidth', lineWidthDefault, 'float')



        return options_out

    def get_options_statstable(self, options, section, options_out):
        options_out['selectByWavelength'] = True  ##option to be always true
        if options_out['wl_values'] is None:
            options_out['wl_values'] = list(np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])))

        options_out['params'] = self.get_value_param(options, section, 'params', self.valid_stats.keys(), 'strlist')

        if self.output_path is not None:
            name_default = options_out['name'] + '.csv'
            file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(options, section, 'file_out', file_out_default, 'str')

        return options_out

    def get_options_statswlplot(self, options, section, options_out):
        options_out['selectByWavelength'] = True  ##option to be always true
        if options_out['wl_values'] is None:
            options_out['wl_values'] = list(np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])))
        options_out['params'] = self.get_value_param(options, section, 'params', self.valid_stats.keys(), 'strlist')

        if self.output_path is not None:
            name_default = options_out['name'] + '.' + self.format_image
            file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(options, section, 'file_out', file_out_default, 'str')
        options_out['title'] = self.get_value_param(options, section, 'title', None, 'str')
        xlabeldefault = defaults.xlabel_wl_default
        ylabeldefault = defaults.ylabel_default
        options_out['xlabel'] = self.get_value_param(options, section, 'xlabel', xlabeldefault, 'str')
        options_out['ylabel'] = self.get_value_param(options, section, 'ylabel', ylabeldefault, 'str')
        line_color = ['black']
        if options_out['selectValues'] is not None:
            nvalues = len(options_out['selectValues'])
            line_color = defaults.get_color_list(nvalues)
        line_type = ['-']
        line_width = [1]
        marker = ['.']
        marker_size = [10]
        options_out['line_color'] = self.get_value_param(options, section, 'line_color', line_color, 'strlist')
        options_out['line_type'] = self.get_value_param(options, section, 'line_type', line_type, 'strlist')
        options_out['line_width'] = self.get_value_param(options, section, 'line_width', line_width, 'floatlist')
        options_out['marker'] = self.get_value_param(options, section, 'marker', marker, 'strlist')
        options_out['marker_size'] = self.get_value_param(options, section, 'marker_size', marker_size, 'floatlist')

        return options_out

    def get_options_spectraplot(self, options, section, options_out):
        options_out['type_rrs'] = self.get_value_param(options, section, 'type_rrs', 'ins', 'str')
        options_out['wl_min'] = self.get_value_param(options, section, 'wl_min', None, 'float')
        options_out['wl_max'] = self.get_value_param(options, section, 'wl_max', None, 'float')
        options_out['plot_stats'] = self.get_value_param(options, section, 'plot_stats', True, 'boolean')
        options_out['plot_spectra'] = self.get_value_param(options, section, 'plot_spectra', ['All'], 'strlist')
        if options_out['plot_spectra'][0].lower() == 'none':
            options_out['plot_spectra'] = None
        if self.output_path is not None:
            name_default = options_out['name'] + '.' + self.format_image
            file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(options, section, 'file_out', file_out_default, 'str')
        sfdefault = 1000
        options_out['scale_factor'] = self.get_value_param(options, section, 'scale_factor', sfdefault, 'float')
        options_out['title'] = self.get_value_param(options, section, 'title', None, 'str')
        xlabeldefault = defaults.xlabel_wl_default
        if options_out['type_rrs'] == 'ins':
            ylabeldefault = defaults.xlabel_default
        elif options_out['type_rrs'] == 'sat':
            ylabeldefault = defaults.ylabel_default
        elif options_out['type_rrs'].startswith('mu_') or options_out['type_rrs'].startswith('flag_'):
            ylabeldefault = defaults.ylabel_rrs_scaled
        options_out['xlabel'] = self.get_value_param(options, section, 'xlabel', xlabeldefault, 'str')
        options_out['ylabel'] = self.get_value_param(options, section, 'ylabel', ylabeldefault, 'str')

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

    def get_flag_list(self, values, allValues, allFlags):
        flag_list = []
        for val in values:
            indext = np.where(allValues == val)
            index = indext[0]
            if len(index) == 1:
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
        #rel_diff[rel_diff>500] = rel_diff[rel_diff>500]/10

        self.valid_stats['RPD'] = np.mean(rel_diff)

        #  the mean of absolute (unsigned) percent differences
        self.valid_stats['APD'] = np.mean(np.abs(rel_diff))

        bias = np.mean(sat_obs - ref_obs)
        self.valid_stats['BIAS'] = bias

        mae = np.mean(np.abs(sat_obs - ref_obs))
        self.valid_stats['MAE'] = mae

        self.valid_stats['DETER(r2)'] = r_value * r_value
