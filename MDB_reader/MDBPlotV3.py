from MDBFile import MDBFile
from PlotOptions import PlotOptions
import MDBPlotDefaults as defaults
import os
import numpy as np
import netCDF4
import math
from scipy import stats
import COMMON.common_functions as cfs


class MDBPlot:

    def __init__(self, path_mdbr_file):

        self.mrfile = None
        self.path_mdbr_file = path_mdbr_file
        self.VALID = False
        if path_mdbr_file is not None:
            self.mrfile = MDBFile(path_mdbr_file)
            self.VALID = self.mrfile.VALID

        self.global_options = None
        self.groupdata = []
        self.yregress = []
        self.xregress = []
        self.valid_stats = {}
        for s in defaults.valid_stats:
            self.valid_stats[defaults.valid_stats[s]['name']] = 0.0

    def compute_statistics(self, use_log_scale, use_rhow, type_regression):

        self.valid_stats['N'] = len(self.xdata)
        if self.valid_stats['N'] == 0:
            for key in self.valid_stats:
                self.valid_stats[key] = np.nan
            return

        self.valid_stats['NMU'] = self.valid_stats['N']
        self.valid_stats['NGROUP'] = self.valid_stats['N']

        # Generated linear fit
        xdatal = []
        ydatal = []
        maxxy = 0
        minxy = None
        for x, y in zip(self.xdata, self.ydata):
            if use_log_scale:
                if x > 0 and y > 0:
                    xdatal.append(math.log10(x))
                    ydatal.append(math.log10(y))
            else:
                # if np.isnan(x) or np.isnan(y):
                #     print(x, y)
                if use_rhow:
                    x = x * np.pi
                    y = y * np.pi
                xdatal.append(x)
                ydatal.append(y)
            if minxy is None:
                if x <= y:
                    minxy = x
                else:
                    minxy = y
            else:
                if x < minxy:
                    minxy = x
                if y < minxy:
                    minxy = y
            if x > maxxy:
                maxxy = x
            if y > maxxy:
                maxxy = y

        ##REGRESSION I
        slope, intercept, r_value, p_value, std_err = stats.linregress(xdatal, ydatal)
        self.valid_stats['slope_I'] = slope
        self.valid_stats['intercept_I'] = intercept
        self.valid_stats['PCC(r)'] = r_value
        self.valid_stats['p_value'] = p_value
        self.valid_stats['std_err_I'] = std_err

        ##REGRESSION II
        from pylr2 import regress2
        reg_2_valid = True
        try:
            results = regress2(np.array(xdatal, dtype=np.float64), np.array(ydatal, dtype=np.float64),
                               _method_type_2="reduced major axis")
            slope_II = results['slope']
            intercept_II = results['intercept']
            self.valid_stats['slope_II'] = slope_II
            self.valid_stats['intercept_II'] = intercept_II
            self.valid_stats['std_slope_II'] = results['std_slope']
            self.valid_stats['std_intercept_II'] = results['std_intercept']
        except:
            reg_2_valid = False

        if type_regression == 'I':
            self.xregress, self.yregress = self.get_regression_line(xdatal, ydatal, slope, intercept, minxy, maxxy)
        elif type_regression == 'II':
            if reg_2_valid:
                self.xregress, self.yregress = self.get_regression_line(xdatal, ydatal, slope_II, intercept_II, minxy,
                                                                        maxxy)

        ref_obs = np.asarray(self.xdata, dtype=np.float64)
        sat_obs = np.asarray(self.ydata, dtype=np.float64)
        if use_rhow:
            sat_obs = sat_obs * np.pi
            ref_obs = ref_obs * np.pi
        # the mean of relative (signed) percent differences
        rel_diff = 100 * ((sat_obs - ref_obs) / ref_obs)
        self.valid_stats['RPD'] = np.mean(rel_diff)
        #  the mean of absolute (unsigned) percent differences
        self.valid_stats['APD'] = np.mean(np.abs(rel_diff))
        if use_log_scale:
            sat_obs = np.log10(sat_obs)
            ref_obs = np.log10(ref_obs)
        self.valid_stats['RMSD'] = cfs.rmse(sat_obs, ref_obs)
        ref_mean = np.mean(ref_obs)
        sat_mean = np.mean(sat_obs)
        self.valid_stats['XAVG'] = ref_mean
        self.valid_stats['YAVG'] = sat_mean
        # CPRMSE
        xdiff = ref_obs - ref_mean
        ydiff = sat_obs - sat_mean
        cprmse = cfs.rmse(ydiff, xdiff)
        self.valid_stats['CRMSE'] = cprmse
        # bias
        bias = np.mean(sat_obs - ref_obs)
        self.valid_stats['BIAS'] = bias
        # mae
        mae = np.mean(np.abs(sat_obs - ref_obs))
        self.valid_stats['MAE'] = mae
        # deter(r2)
        self.valid_stats['DETER(r2)'] = r_value * r_value

        # print(self.valid_stats)

    def get_regression_line(self, xdatal, ydatal, slope, intercept, minxy, maxxy):
        if maxxy is None:
            maxxy = 0
            for x, y in zip(xdatal, ydatal):
                if x > maxxy:
                    maxxy = x
                if y > maxxy:
                    maxxy = y
        if slope is None and intercept is None:
            slope, intercept, r_value, p_value, std_err = stats.linregress(xdatal, ydatal)

        xregress = []
        yregress = []
        xregress.append(0)
        yregress.append(intercept)
        for x in xdatal:
            yr = (x * slope) + intercept
            yregress.append(yr)
            xregress.append(x)
        yrmax = ((maxxy + 1) * slope) + intercept
        xregress.append(maxxy + 1)
        yregress.append(yrmax)

        if minxy is not None:
            yrmin = ((minxy - 1) * slope) + intercept
            xregress.append(minxy - 1)
            yregress.append(yrmin)

        return xregress, yregress

    def plot_from_options_file(self, file_config):
        if not os.path.isfile(file_config):
            print(f'[ERROR] File config {file_config} is not a valid config file')
            return
        import configparser
        try:
            options = configparser.ConfigParser()
            options.read(file_config)
        except:
            print(f'[ERROR] Error reading file_config: {file_config}')
        self.plot_from_options(options)

    def plot_from_options(self, options):
        poptions = PlotOptions(options, None)
        poptions.set_global_options()
        self.global_options = poptions.global_options
        list_figures = poptions.get_list_figures()
        # print(list_figures)
        for figure in list_figures:
            print('------------------------------------------------------------------------------------------')
            print(f'[INFO] Starting figure: {figure}')
            options_figure = poptions.get_options(figure)
            if options_figure['selectBy'] is not None:
                options_figure = self.check_gs_options_impl(options_figure, 'selectBy', 'selectType', 'selectValues')
            self.plot_from_options_impl(options_figure)

    def plot_from_options_impl(self, options_figure):
        if options_figure['type'] == 'scatterplot':
            self.plot_scatterplot_from_options(options_figure)
        if options_figure['type'] == 'spectraplot':
            self.plot_spectraplot_from_options(options_figure)

    def plot_scatterplot_from_options(self, options_figure):

        ##WORKING WITH ALL THE DATA
        if options_figure['selectBy'] is None:
            if not options_figure['selectByWavelength']:  # GLOBAL SCATTERPLOT
                self.plot_global_scatterplot(options_figure)
            else:  # AN SCATTERPLOT BY WAVELENGTH
                if options_figure['multiple_plot'] is not None:  # single file
                    self.plot_multiple_wavelength_scatterplots_single_file(options_figure)
                else:  # multiple files
                    self.plot_multiple_wavelength_scatterplots_multiple_files(options_figure)

        ##WORKING WITH SELECTED OPTIONS
        if options_figure['selectBy'] is not None:
            selectValues = options_figure['selectValues']
            file_out_base = options_figure['file_out']
            title_base = options_figure['title']
            for svalue in selectValues:
                options_figure['selectValues'] = svalue
                flag = self.get_str_select_value(options_figure, svalue)
                options_figure['file_out'] = self.get_file_out_name(file_out_base, None, flag)
                options_figure['title'] = self.get_title(title_base, None, flag, None)
                if not options_figure['selectByWavelength']:  # GLOBAL SCATTERPLOT
                    self.plot_global_scatterplot(options_figure)
                else:  # AN SCATTERPLOT BY WAVELENGTH
                    if options_figure['multiple_plot'] is not None:  # single file
                        self.plot_multiple_wavelength_scatterplots_single_file(options_figure)
                    else:  # multiple files
                        self.plot_multiple_wavelength_scatterplots_multiple_files(options_figure)

    def plot_global_scatterplot(self, options_figure):
        if options_figure['apply_wavelength_color'] and options_figure['groupBy'] is None:
            options_figure['groupBy'] = 'mu_wavelength'
            options_figure['groupValues'] = options_figure['wlvalues']
            options_figure['groupType'] = 'wavelength'
        self.set_data_scatterplot(options_figure['groupBy'], options_figure['selectBy'], options_figure['selectValues'],
                                  None, options_figure)
        self.plot_scatter_plot(options_figure, None, -1, -1, -1)

    def plot_multiple_wavelength_scatterplots_single_file(self, options_figure):
        file_out_final = options_figure['file_out']
        title_base = options_figure['title']
        from PlotScatter import PlotScatter
        rc = options_figure['multiple_plot'].split(',')
        nrow = int(rc[0].strip())
        ncol = int(rc[1].strip())
        ntot = nrow * ncol

        index = 0
        plot_here = PlotScatter()
        plot_here.nrow = nrow
        plot_here.ncol = ncol
        plot_here.index_row = 0
        plot_here.index_col = 0
        plot_here.xtitle_options['fontsize'] = options_figure['fontsizelabels']
        plot_here.ytitle_options['fontsize'] = options_figure['fontsizelabels']
        plot_here.plot_text_options['fontsize'] = options_figure['fontsizestats']

        plot_here.start_multiple_plot_advanced(nrow, ncol, options_figure['xfigsize'],
                                               options_figure['yfigsize'], options_figure['widthspace'],
                                               options_figure['heightspace'])
        print(f'[INFO] Starting multiple plot with {nrow} rows and {ncol} cols')
        wl_values = options_figure['wlvalues']
        print(f'[INFO] Wavelenthts: {wl_values}')
        nblank = ntot - len(wl_values)
        if nblank >= nrow:
            print(
                f'[ERROR] Number of total axis {plot_here.nrow}x{plot_here.ncol} is higher than number of scatterplots:{len(wl_values)}')
            plot_here.close_plot()
            return
        print(f'[INFO] Number of total axis {plot_here.nrow}x{plot_here.ncol} Blank plots: {nblank}')
        if nblank > 0:
            index_col_adjust = plot_here.ncol - nblank
        for wl in wl_values:
            selectBy = None
            selectValue = None
            if options_figure['selectBy'] is not None:
                selectBy = options_figure['selectBy']
                selectValue = options_figure['selectValues']
            self.set_data_scatterplot(options_figure['groupBy'], selectBy, selectValue, wl, options_figure)
            if len(self.xdata) > 0 and len(self.ydata) > 0:
                options_figure['title'] = self.get_title(title_base, wl, None, None)
                print(f'[INFO] Plotting scatter plot for wavelength {wl} ({len(self.xdata)} points)')
                options_figure['file_out'] = None

                self.plot_scatter_plot(options_figure, plot_here, index, wl, index_col_adjust)

                plot_here.index_col = plot_here.index_col + 1
                if plot_here.index_col == plot_here.ncol:
                    plot_here.index_col = 0
                    plot_here.index_row = plot_here.index_row + 1
                index = index + 1
            else:
                print(f'[WARNING] No data for wavelength: {wl} nm')

        for index_blank in range(index, ntot):
            plot_here.plot_blanck(index_blank)

        if options_figure['legend']:
            str_legend = self.get_str_legend(options_figure)
            if len(str_legend) > 0:
                plot_here.set_global_legend(str_legend)

        plot_here.save_fig(file_out_final)
        plot_here.close_plot()

    def plot_multiple_wavelength_scatterplots_multiple_files(self, options_figure):
        file_out_base = options_figure['file_out']
        title_base = options_figure['title']
        wl_values = options_figure['wlvalues']
        print(f'[INFO] Wavelenthts: {wl_values}')
        for wl in wl_values:
            selectBy = None
            selectValue = None
            if options_figure['selectBy'] is not None:
                selectBy = options_figure['selectBy']
                selectValue = options_figure['selectValues']
            self.set_data_scatterplot(options_figure['groupBy'], selectBy, selectValue, wl, options_figure)
            if len(self.xdata) > 0 and len(self.ydata) > 0:
                print(f'[INFO] Plotting scatter plot for wavelength {wl} ({len(self.xdata)} points)')
                options_figure['title'] = self.get_title(title_base, wl, None, None)
                options_figure['file_out'] = self.get_file_out_name(file_out_base, wl, None)
                self.plot_scatter_plot(options_figure, None, -1, wl, -1)
            else:
                print(f'[WARNING] No data for wavelength: {wl} nm')

    # MAIN FUNCTION TO PLOT SCATTERPLOT
    def plot_scatter_plot(self, options, plot, index, wl, index_col_adjust):

        ##compute statistics if neeed
        use_rhow = options['use_rhow']

        if options['include_stats'] or options['regression_line']:
            use_log_scale = options['log_scale']
            self.compute_statistics(use_log_scale, use_rhow, options['type_regression'])

        # check groups and get legend if applicable
        ngroup = 1
        str_legend = []
        groupValues = None
        if 'groupValues' in options.keys():
            groupValues = options['groupValues']
        if groupValues is not None:
            ngroup = len(groupValues)
        if ngroup > 1 and len(self.groupdata) > 0 and options['legend']:
            str_legend = self.get_str_legend(options)
        if ngroup > 1:
            self.valid_stats['NGROUP'] = int(self.valid_stats['N'] / ngroup)
        if options['wlvalues'] is not None:
            nwl = len(options['wlvalues'])
            self.valid_stats['NMU'] = int(self.valid_stats['N'] / nwl)

        # start plot
        from scipy.stats import gaussian_kde
        if plot is None and index == -1:
            from PlotScatter import PlotScatter
            plot = PlotScatter()
            plot.close_plot()
            plot.start_plot()
        if plot is not None and index >= 0:
            plot.set_axhere_index(index)

        ##check x y data
        if options['scale_factor'] is not None:
            self.xdata = self.xdata * options['scale_factor']
            self.ydata = self.ydata * options['scale_factor']
            if len(self.yregress) > 0 and len(self.xregress) > 0:
                self.yregress = np.array(self.yregress) * options['scale_factor']
                self.xregress = np.array(self.xregress) * options['scale_factor']
        if use_rhow:
            self.xdata = self.xdata * np.pi
            self.ydata = self.ydata * np.pi

        # check scatter options
        colors = options['color']
        color = colors[0]
        markersizes = options['markersize']
        markersize = markersizes[0]
        markers = options['marker']
        marker = markers[0]
        edgecolors = options['edgecolor']
        edgecolor = edgecolors[0]
        linewidths = options['linewidth']
        linewidth = linewidths[0]

        ##plotting implementation
        ngroupReal = 0
        if ngroup > 1:  ##scatter plots with points coloured by group
            nmubygroup = [0] * ngroup
            for idx in range(ngroup):  # groupValues:
                g = groupValues[idx]
                if len(colors) == ngroup:
                    color = colors[idx]
                else:
                    if options['groupType'] == 'flag' or options['groupType'] == 'float':
                        if ngroup <= len(defaults.colors_default):
                            color = defaults.colors_default[idx]
                        else:
                            color = defaults.get_color_default(idx, 0, ngroup - 1)
                    else:
                        color = defaults.get_color_wavelength(g)
                xhere = self.xdata[self.groupdata == g]
                yhere = self.ydata[self.groupdata == g]

                if len(markers) == ngroup:
                    marker = markers[idx]
                else:
                    marker = markers[0]

                if len(markersizes) == ngroup:
                    markersize = markersizes[idx]
                else:
                    markersize = markersizes[0]

                if len(edgecolors) == ngroup:
                    edgecolor = edgecolors[idx]
                else:
                    edgecolor = edgecolors[0]

                if len(linewidths) == ngroup:
                    linewidth = linewidths[idx]
                else:
                    linewidth = linewidths[0]

                if len(xhere) > 0 and len(yhere) > 0:
                    nmubygroup[idx] = len(xhere)
                    ngroupReal = ngroupReal + 1

                plot.plot_data(xhere, yhere, marker, markersize, color, edgecolor, linewidth)
        else:  # density or normal scatter plot

            xhere = np.asarray(self.xdata, dtype=np.float64)
            yhere = np.asarray(self.ydata, dtype=np.float64)

            # Density
            if options['apply_density']:
                if options['log_scale']:
                    xherel = np.log10(xhere)
                    yherel = np.log10(yhere)
                    xy = np.vstack([xherel, yherel])
                else:
                    xy = np.vstack([xhere, yhere])

                try:
                    z = gaussian_kde(xy)(xy)
                    idx = z.argsort()
                    xhere, yhere, z = xhere[idx], yhere[idx], z[idx]
                    plot.set_cmap('jet')

                    plot.plot_data(xhere, yhere, marker, markersize, z, edgecolor, linewidth)

                except:

                    plot.plot_data(xhere, yhere, marker, markersize, color, edgecolor, linewidth)

            else:
                plot.plot_data(xhere, yhere, marker, markersize, color, edgecolor, linewidth)

        ##limitss
        if options['log_scale']:
            plot.set_log_scale()
        else:
            if options['min_xy'] is None or options['max_xy'] is None:
                min_xy, max_xy = self.get_min_max_xy()
        if options['min_xy'] is not None:
            min_xy = options['min_xy']
        if options['max_xy'] is not None:
            max_xy = options['max_xy']
        plot.set_limits(min_xy, max_xy)

        # ticks
        if options['ticks'] is None:
            ticks = self.get_ticks_from_min_max_xy(min_xy, max_xy)
        else:
            ticks = options['ticks']
        if ticks is not None:
            if options['log_scale']:
                tlabels = self.get_labels_for_log_ticks(ticks)
                plot.set_ticks_and_labels(ticks, tlabels, options['fontsizeaxis'])
            else:
                plot.set_ticks(ticks, options['fontsizeaxis'])

        ##x-y labels
        if options['individual_axis'] or index == -1:
            plot.set_xaxis_title(options['xlabel'])
            plot.set_yaxis_title(options['ylabel'])
        else:
            final_row = plot.nrow - 1
            prefinal_row = plot.nrow - 2
            if plot.index_row == final_row:
                plot.set_xaxis_title(options['xlabel'])
            if plot.index_col == 0:
                plot.set_yaxis_title(options['ylabel'])
            if plot.index_col > 0:
                plot.set_yticks_labels_off(ticks)
            if plot.index_row == prefinal_row and plot.index_col >= index_col_adjust >= 1:
                plot.set_xaxis_title(options['xlabel'])
            else:
                if plot.index_row < final_row:
                    plot.set_xticks_labels_off(ticks)

        plot.set_equal_apect()

        # legend
        if options['legend'] and len(str_legend) > 0 and index == -1:
            plot.set_legend(str_legend)

        # identity liine
        if options['identity_line']:
            plot.plot_identity_line()

        # stats
        if options['include_stats'] and options['stat_list'] is not None:
            str0 = self.get_str_stats(options, wl)
            xpos = options['stats_xpos']
            ypos = options['stats_ypos']
            plot.plot_text(xpos, ypos, str0)

        # regression lines
        if not options['log_scale']:
            if options['regression_line']:
                plot.plot_regress_line(self.xregress, self.yregress, 'black')
            if options['regression_line_groups']:
                if ngroup > 1:
                    for idx in range(len(groupValues)):
                        g = groupValues[idx]
                        if options['groupType'] == 'flag' or options['groupType'] == 'float':
                            # color = defaults.get_color_flag(g)
                            color = defaults.colors_default[idx]
                        else:
                            color = defaults.get_color_wavelength(g)
                        xhere = self.xdata[self.groupdata == g]
                        yhere = self.ydata[self.groupdata == g]
                        if len(xhere) > 0 and len(yhere) > 0:
                            type_regression = options['type_regression']
                            slope = self.valid_stats[f'slope_{type_regression}']
                            intercept = self.valid_stats[f'intercept_{type_regression}']
                            xregress, yregress = self.get_regression_line(xhere, yhere, slope, intercept, min_xy,
                                                                          max_xy)
                            plot.plot_regress_line(xregress, yregress, color)

        # title
        if options['title'] is not None:
            title_here = options['title']
            plot.set_title(title_here)
            plot.axhere.title.set_size(options['fontsizetitle'])

        ##saving to file
        if not options['file_out'] is None:
            plot.save_fig(options['file_out'])
            plot.close_plot()

        return plot

    def plot_spectraplot_from_options(self,options_figure):
        if options_figure['type_rrs']=='ins':
            self.plot_insitu_spectraplots(options_figure)

    def plot_insitu_spectraplots(self,options_figure):
        ##GETTING DATA
        wavelength = self.mrfile.get_insitu_wl()
        all_spectra, all_spectra_validity, spectra_stats = self.mrfile.get_all_insitu_spectra(options_figure['scale_factor'],options_figure['use_rhow'],options_figure['plot_stats'])

        # from PlotSpectra import PlotSpectra
        # if not options_out['plot_stats']:
        #     stats = None

    def check_select_group_options(self, options_figure):
        if options_figure['selectByWavelength'] or options_figure['apply_wavelength_color']:
            wl_values_ini = options_figure['wlvalues']
            if wl_values_ini is None:
                wl_values = np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])).tolist()
            else:
                wl_values_unique = np.unique(np.array(self.mrfile.nc.variables['mu_wavelength']))
                wl_values = []
                for wl in wl_values_ini:
                    imin = np.argmin(np.abs(wl - wl_values_unique))
                    if abs(wl - wl_values_unique[imin]) <= 1:
                        wl_values.append(wl_values_unique[imin])
                    else:
                        return None
            wl_values_str = []
            wl_values_sat = self.mrfile.satellite_bands
            for wl in wl_values:
                imin = np.argmin(np.abs(wl - wl_values_sat))
                wl_sat = wl_values_sat[imin]
                wl_sat_str = f'{wl_sat:.2f}'
                if wl_sat_str.endswith('.00'):
                    wl_sat_str = wl_sat_str[:-3]
                if wl_sat_str.endswith('0') and wl_sat_str.find('.') > 0:
                    wl_sat_str = wl_sat_str[:-1]
                wl_values_str.append(wl_sat_str)

            options_figure['wlvalues'] = wl_values
            options_figure['wlvalues_str'] = wl_values_str

        if options_figure['groupBy'] is not None:
            options_figure = self.check_gs_options_impl(options_figure, 'groupBy', 'groupType', 'groupValues')

        if options_figure is None:
            return None
        if options_figure['selectBy'] is not None:
            options_figure = self.check_gs_options_impl(options_figure, 'selectBy', 'selectType', 'selectValues')

        return options_figure

    def check_gs_options_impl(self, options_figure, by, type, values):
        var_group_name = options_figure[by]
        if options_figure[type] == 'flag':
            if var_group_name in self.mrfile.variables:
                flag_values = self.mrfile.nc.variables[var_group_name].flag_values
                flag_meanings_list = self.mrfile.nc.variables[var_group_name].flag_meanings.split(' ')
                flag_meanings = [x.strip() for x in flag_meanings_list]
                options_figure[var_group_name] = {
                    'flag_values': flag_values,
                    'flag_meanings': flag_meanings
                }
            else:  ##virtual flag
                print('virtual flag')
            if options_figure[values] is None:
                options_figure[values] = flag_values
            else:
                flag_list_config = options_figure[values]
                flag_values_config = []
                for flag_config in flag_list_config:
                    if flag_config.strip() == 'GLOBAL':
                        flag_values_config.append(-1)
                        continue
                    try:
                        iflag = flag_meanings.index(flag_config.strip())
                        flag_values_config.append(flag_values[iflag])
                    except:
                        print(f'[WARNING] Flag {flag_config.strip()} is not in the list')
                        return None
                options_figure[values] = flag_values_config

        if options_figure[type] == 'float':
            if var_group_name not in self.mrfile.variables:
                return None
            all_group_values = np.unique(np.array(self.mrfile.nc.variables[var_group_name]))
            if options_figure[values] is None:
                options_figure[values] = list(all_group_values)
            else:
                group_values_given = options_figure[values]
                group_values = []
                for val in group_values_given:
                    imin = np.argmin(np.abs(val - all_group_values))
                    if abs(val - all_group_values[imin]) < 0.1:
                        group_values.append((all_group_values[imin]))
                    else:
                        print(f'[WARNING] Value {val} is not in the variable {var_group_name}')
                        return None
                options_figure[values] = group_values
        return options_figure

    def set_data_scatterplot(self, groupBy, selectBy, valSelect, wl_value, options_out):
        rrs_ins = np.array(self.mrfile.nc.variables['mu_ins_rrs'])
        rrs_sat = np.array(self.mrfile.nc.variables['mu_sat_rrs'])
        id_mu = np.array(self.mrfile.nc.variables['mu_satellite_id'])

        mu_valid_variable = self.global_options['mu_valid_variable']
        mu_valid_satelliteid = np.array(self.mrfile.nc.variables[mu_valid_variable])

        valid_mu = self.get_array_muid_from_array_satelliteid(id_mu, mu_valid_satelliteid)

        valid_all = self.check_rrs_valid(valid_mu, rrs_ins, rrs_sat)

        if wl_value is not None:
            wl_array = np.array(self.mrfile.nc.variables['mu_wavelength'])
            valid_all[wl_array != wl_value] = 0

        if selectBy is not None and valSelect is not None:
            if selectBy in options_out.keys() and 'flag_array' in options_out[selectBy]:
                select_array = options_out[selectBy]['flag_array']
            else:
                select_array = np.array(self.mrfile.nc.variables[selectBy])
            if len(select_array) == len(mu_valid_satelliteid):
                select_array = self.get_array_muid_from_array_satelliteid(id_mu, select_array)
            valid_all[select_array != valSelect] = 0

        self.xdata = rrs_ins[valid_all == 1]
        self.ydata = rrs_sat[valid_all == 1]

        if groupBy is not None:
            if options_out['groupType'] == 'float' or options_out['groupType'] == 'wavelength':
                group_array = np.array(self.mrfile.nc.variables[groupBy])
            else:
                group_array, all_group_values, all_group_meanings = self.get_flag_array(options_out, 'groupBy')
            if len(group_array) == len(mu_valid_satelliteid):
                group_array = self.get_array_muid_from_array_satelliteid(id_mu, group_array)
            self.groupdata = group_array[valid_all == 1]

    # options_var: selectBy or groupBy
    def get_flag_array(self, options_out, option_var):
        var_flag = options_out[option_var]
        if var_flag in self.mrfile.variables:
            array_flag = np.array(self.mrfile.variables[var_flag][:])
            flag_values = self.mrfile.variables[var_flag].flag_values
            flag_meanings = self.mrfile.variables[var_flag].flag_meanings.split(' ')
        else:  ##previously virtual flag
            array_flag = options_out[var_flag]['flag_array']
            flag_values = options_out[var_flag]['flag_values']
            flag_meanings = options_out[var_flag]['flag_meanings']

        return array_flag, flag_values, flag_meanings

    def get_flag_list(self, values, allValues, allFlags):
        flag_list = []
        for val in values:
            if val == -1:
                flag_list.append('GLOBAL')
            indext = np.where(np.array(allValues) == val)
            index = indext[0]
            if len(index) == 1:
                indexf = index[0]
                flag_list.append(allFlags[indexf])
        return flag_list

    def get_flag_flag(self, val, allValues, allFlags):
        indext = np.where(allValues == val)
        index = indext[0]
        if len(index) == 1:
            indexf = index[0]
            return allFlags[indexf]
        return None

    def get_wl_str_from_wl(self, wl_value):
        wl_sat = np.array(self.mrfile.nc.variables['satellite_bands'])
        index_sat = np.argmin(np.abs(wl_sat - wl_value))
        wl_sat_value = wl_sat[index_sat]
        wl_sat_value_str = f'{wl_sat_value:.2f}'
        if wl_sat_value_str.endswith('.00'):
            return wl_sat_value_str[:-3]
        else:
            return wl_sat_value_str

    def get_str_stat(self, val, format_complete, units):
        format = format_complete.split('+')[0]
        if format == 'f0' or format == 'i' or format == 'e0':
            val_str = f'{val:.0f}'
        elif format == 'f1':
            val_str = f'{val:.1f}'
        elif format == 'f2':
            val_str = f'{val:.2f}'
        elif format == 'f3':
            val_str = f'{val:.3f}'
        elif format == 'f4':
            val_str = f'{val:.4f}'
        elif format == 'f5':
            val_str = f'{val:.5f}'
        elif format == 'f6':
            val_str = f'{val:.6f}'
        elif format == 'e1':
            val_str = f'{val:.1e}'
        elif format == 'e2':
            val_str = f'{val:.2e}'
        elif format == 'e3':
            val_str = f'{val:.3e}'
        elif format == 'e4':
            val_str = f'{val:.4e}'
        elif format == 'e5':
            val_str = f'{val:.5e}'
        elif format == 'e6':
            val_str = f'{val:.6e}'
        if len(format_complete.split('+')) == 2:
            if format_complete.split('+')[1].lower() == 'units':
                if len(units) > 0:
                    val_str = f'{val_str} {units}'
        return val_str

    def get_str_stats(self, options, wl):
        stat_list = options['stat_list']
        str0 = ''
        for stat in stat_list:
            if len(str0) > 0:
                str0 = f'{str0}\n'
            if stat.upper() == 'WL':
                wls = self.get_wl_str_from_wl(wl)
                str0 = f'{str0}{wls} nm'
            elif stat.upper() == 'EQUATION':
                typeRegression = options['type_regression']
                val_slope = self.valid_stats[defaults.valid_stats[f'SLOPE_{typeRegression}']['name']]
                val_slope = self.get_str_stat(val_slope, options[f'SLOPE_{typeRegression}_FORMAT'], '')
                val_offset = self.valid_stats[defaults.valid_stats[f'OFFSET_{typeRegression}']['name']]
                sign = '+'
                if val_offset < 0:
                    sign = '-'

                val_offset = self.get_str_stat(abs(val_offset), options[f'OFFSET_{typeRegression}_FORMAT'], '')

                str0 = f'{str0}y = {val_slope.strip()}x {sign} {val_offset.strip()}'
            else:
                val = self.valid_stats[defaults.valid_stats[stat.upper()]['name']]

                valstr = self.get_str_stat(val, options[f'{stat.upper()}_FORMAT'], options['units'])
                name_plot = options[f'{stat.upper()}_NAMEPLOT']
                str0 = f'{str0}{name_plot} = {valstr}'

        return str0

    def get_str_select_value(self, options, svalue):
        if options['selectType'] == 'float':
            str_out = f'{svalue:.2f}'
        if options['selectType'] == 'flag':
            flag_name = options['selectBy']
            str_out = self.get_flag_flag(svalue, options[flag_name]['flag_values'], options[flag_name]['flag_meanings'])
        return str_out

    def get_str_legend(self, options):
        if options['legend_values'] is not None:
            return options['legend_values']
        str_legend = []
        groupValues = options['groupValues']
        if len(self.groupdata) > 0 and groupValues is not None:
            ngroup = len(groupValues)
            if ngroup > 1:
                if options['groupType'] == 'float' or options['groupType'] == 'wavelength':
                    for g in groupValues:
                        str_legend.append(f'{g:.2f}')
                if options['groupType'] == 'flag':
                    flag_name = options['groupBy']
                    str_legend = self.get_flag_list(groupValues, options[flag_name]['flag_values'],
                                                    options[flag_name]['flag_meanings'])
                    if 'FUB' in str_legend:
                        index = str_legend.index('FUB')
                        if index >= 0:
                            str_legend[index] = 'S3 FUB-CSIRO'
                    if 'STANDARD' in str_legend:
                        index = str_legend.index('STANDARD')
                        if index >= 0:
                            str_legend[index] = 'WFR'
                    if 'POLYMER' in str_legend:
                        index = str_legend.index('POLYMER')
                        if index >= 0:
                            str_legend[index] = 'CMEMS-OLCI'
                    if 'CCIALL' in str_legend:
                        index = str_legend.index('CCIALL')
                        if index >= 0:
                            str_legend[index] = 'OC-CCI v.6 (complete time series)'
                    if 'CCI' in str_legend:
                        index = str_legend.index('CCI')
                        if index >= 0:
                            str_legend[index] = 'OC-CCI v.6 (OLCI period)'

        return str_legend

    def get_min_max_xy(self):
        min_xy = None
        max_xy = None
        if len(self.xdata) > 0:
            min_x_data = np.min(self.xdata)
            min_y_data = np.min(self.ydata)
            min_xy = np.floor(np.min([min_x_data, min_y_data]))
            max_x_data = np.max(self.xdata)
            max_y_data = np.max(self.ydata)
            max_xy = np.ceil(np.max([max_x_data, max_y_data]))
        return min_xy, max_xy

    def get_ticks_from_min_max_xy(self, min_xy, max_xy):
        ticks = []
        if min_xy is None or max_xy is None:
            return None
        dif = max_xy - min_xy
        increm = 1
        if dif >= 8:
            increm = 2

        for v in range(int(min_xy), int(max_xy) + 1, increm):
            if v <= max_xy:
                ticks.append(v)

        return ticks

    def get_labels_for_log_ticks(self, ticks):
        tlabels = []
        for t in ticks:
            tl = np.log10(t)
            if tl < 0:
                tls = str(t)
            else:
                tls = f'{t:.1f}'
            tlabels.append(tls)
        return tlabels

    def get_array_muid_from_array_satelliteid(self, id_mu, satellite_array):
        mu_array_out = np.zeros(id_mu.shape, dtype=satellite_array.dtype)
        for id in range(len(satellite_array)):
            mu_array_out[id_mu == id] = satellite_array[id]
        return mu_array_out

    def check_rrs_valid(self, valid_array, rrs1, rrs2):
        dfv = netCDF4.default_fillvals['f4']
        for id in range(len(valid_array)):
            if valid_array[id] == 1:
                if rrs1[id] == dfv or rrs2[id] == dfv:
                    valid_array[id] = 0
        return valid_array

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

    def get_title(self, title, wl, flag, param):
        if title is None:
            return None
        if wl is None and flag is None and param is None:
            return title
        if wl is not None:
            wls = self.get_wl_str_from_wl(wl)
            title = title.replace('$WL$', wls)
        if flag is not None:
            title = title.replace('$FLAG$', flag)
        if param is not None:
            title = title.replace('$PARAM$', param)

        return title
