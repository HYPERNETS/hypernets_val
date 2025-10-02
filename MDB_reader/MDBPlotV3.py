import pandas as pd
import pytz

from MDBFile import MDBFile
from PlotOptions import PlotOptions
import MDBPlotDefaults as defaults
import os
import numpy as np
import netCDF4
import math
from scipy import stats
import COMMON.common_functions as cfs
import COMMON.Class_Flags_OLCI as flag_class


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
        self.xdata = []
        self.ydata = []
        self.valid_stats = {}
        for s in defaults.valid_stats:
            self.valid_stats[defaults.valid_stats[s]['name']] = 0.0

        self.virtual_flags = {}
        self.mu_valid_variable = 'mu_valid'

    def close_mdb_file(self):
        self.mrfile.close()

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
        maxxy = None
        minxy = None
        for x, y in zip(self.xdata, self.ydata):
            if use_log_scale:
                if x > 0 and y > 0:
                    xdatal.append(math.log10(x))
                    ydatal.append(math.log10(y))
            else:
                if use_rhow:
                    x = x * np.pi
                    y = y * np.pi
                xdatal.append(x)
                ydatal.append(y)
            if minxy is None and maxxy is None:
                if x <= y:
                    minxy = x
                    maxxy = y
                else:
                    minxy = y
                    maxxy = x
            else:
                if x < minxy: minxy = x
                if y < minxy: minxy = y
                if x > maxxy: maxxy = x
                if y > maxxy: maxxy = y

        ##REGRESSION I
        slope, intercept, r_value, p_value, std_err = stats.linregress(xdatal, ydatal)
        self.valid_stats['slope_I'] = slope
        self.valid_stats['intercept_I'] = intercept
        self.valid_stats['PCC(r)'] = r_value
        self.valid_stats['p_value'] = p_value
        self.valid_stats['std_err_I'] = std_err

        ##REGRESSION II

        reg_2_valid = True
        try:
            from pylr2 import regress2
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

        if not reg_2_valid:
            type_regression = 'I'


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

        if use_log_scale:
            valid_array = np.logical_and(sat_obs > 0, ref_obs > 0)
            nvalid = np.count_nonzero(valid_array)
            self.valid_stats['N'] = nvalid
            self.valid_stats['NMU'] = nvalid
            self.valid_stats['NGROUP'] = nvalid
            sat_obs = sat_obs[valid_array]
            ref_obs = ref_obs[valid_array]
            self.xregress = np.pow(10,np.array(self.xregress))
            self.yregress = np.pow(10,np.array(self.yregress))

        # the mean of relative (signed) percent differences
        rel_diff = 100 * ((sat_obs - ref_obs) / ref_obs)
        self.valid_stats['RPD'] = np.mean(rel_diff)
        #  the mean of absolute (unsigned) percent differences
        self.valid_stats['APD'] = np.mean(np.abs(rel_diff))

        # the median of relative (signed) percent differences
        rel_diff = 100 * ((sat_obs - ref_obs) / ref_obs)
        self.valid_stats['MdRPD'] = np.median(rel_diff)
        #  the median of absolute (unsigned) percent differences
        self.valid_stats['MdAPD'] = np.median(np.abs(rel_diff))


        self.valid_stats['MIN_Y'] = np.ma.min(sat_obs)
        self.valid_stats['MAX_Y'] = np.ma.max(sat_obs)
        self.valid_stats['MIN_X'] = np.ma.min(ref_obs)
        self.valid_stats['MAX_X'] = np.ma.max(ref_obs)
        self.valid_stats['RANGE_Y'] = self.valid_stats['MAX_Y'] - self.valid_stats['MIN_Y']
        self.valid_stats['RANGE_X'] = self.valid_stats['MAX_X'] - self.valid_stats['MIN_X']

        if use_log_scale:
            sat_obs = np.log10(sat_obs)
            ref_obs = np.log10(ref_obs)

        self.valid_stats['RMSD'] = cfs.rmse(sat_obs, ref_obs)
        ref_mean = np.mean(ref_obs)
        sat_mean = np.mean(sat_obs)
        self.valid_stats['XAVG'] = ref_mean
        self.valid_stats['YAVG'] = sat_mean
        self.valid_stats['XMEDIAN'] = np.median(ref_obs)
        self.valid_stats['YMEDIAN'] = np.median(sat_obs)
        # CPRMSE
        xdiff = ref_obs - ref_mean
        ydiff = sat_obs - sat_mean
        cprmse = cfs.rmse(ydiff, xdiff)
        self.valid_stats['CRMSE'] = cprmse
        # bias (average)
        bias = np.mean(sat_obs - ref_obs)
        self.valid_stats['BIAS'] = bias
        # bias (median)
        meBias = np.median(sat_obs - ref_obs)
        self.valid_stats['MdBIAS'] = meBias

        # mad
        mad = np.mean(np.abs(sat_obs - ref_obs))
        self.valid_stats['MAD'] = mad

        # mdad
        mdad = np.median(np.abs(sat_obs - ref_obs))
        self.valid_stats['MdAD'] = mdad

        # deter(r2)
        self.valid_stats['DETER(r2)'] = r_value * r_value

        #joliff statistics
        std_x = np.std(ref_obs)
        std_y = np.std(sat_obs)
        sig_std = -1 if (std_y - std_x) > 0 else 1
        norm_std = std_y / std_x
        nuRMSD = np.ma.sqrt(1 + norm_std ** 2 - (2 * norm_std * r_value))
        nBias = np.ma.mean(sat_obs - ref_obs) / std_x
        suRMSD = sig_std * nuRMSD
        self.valid_stats['XSTD'] = std_x
        self.valid_stats['YSTD'] = std_y
        self.valid_stats['NORMSTD'] = norm_std
        self.valid_stats['nBIAS'] = nBias
        self.valid_stats['nuRMSD'] = nuRMSD
        self.valid_stats['suRMSD'] = suRMSD


        if use_log_scale:
            ##convert statistict to linear scale again
            stats_to_convert = ['RMSD', 'XAVG', 'YAVG', 'XMEDIAN','YMEDIAN','CRMSE', 'MAD','MdAD']
            for stat in stats_to_convert:
                self.valid_stats[stat] = np.power(10, self.valid_stats[stat])
            sign_stats_to_convert = ['BIAS','MdBIAS']
            for stat in sign_stats_to_convert:
                bias_neg = self.valid_stats[stat] < 0
                self.valid_stats[stat] = np.power(10, np.abs(self.valid_stats[stat]))
                if bias_neg:
                    self.valid_stats[stat] = self.valid_stats[stat] * (-1)

        # print(self.valid_stats)

    def compute_statistics_spectra(self,index_mu,options_figure):
        if not 'scale_factor' in options_figure:
            options_figure['scale_factor'] = 1000
        wl, insitu_spectra, sat_spectra, insitu_spectra_unc, sat_spectra_unc = self.mrfile.get_mu_spectra_insitu_and_sat(
            index_mu, options_figure['scale_factor'])
        if wl is None:
            return None

        indices_vis = np.logical_and(wl>=400.0,wl<=800.0)
        wl = wl[indices_vis]
        insitu_spectra = insitu_spectra[indices_vis]
        sat_spectra = sat_spectra[indices_vis]
        i560 = np.argmin(np.abs(wl-560.0))

        self.xdata = insitu_spectra
        self.ydata = sat_spectra
        self.compute_statistics(False,False,'II')
        self.valid_stats['INDEX_MU'] = index_mu
        self.valid_stats['SAM'] = np.rad2deg(np.acos(np.dot(insitu_spectra,sat_spectra)/(np.linalg.norm(insitu_spectra)*np.linalg.norm(sat_spectra))))
        Yinsitu = insitu_spectra/insitu_spectra[i560]
        Ysat = sat_spectra/sat_spectra[i560]
        self.valid_stats['CHI-SQUARE'] = np.sum(Yinsitu-Ysat/Yinsitu)

        return self.valid_stats

    def get_table_spectral_statistics(self,options_figure):
        wlvalues = options_figure['wlvalues']
        if wlvalues is None:
            wlvalues = list(np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])))

        col_names = ['GLOBAL']
        for wl in wlvalues:
            col_names.append(self.get_wl_str_from_wl(wl))
        df = pd.DataFrame(index=list(self.valid_stats.keys()), columns=col_names)

        #global
        self.set_data_scatterplot(None,None,None,None,options_figure)
        self.compute_statistics(False,False,'II')
        self.valid_stats['NMU'] = int(self.valid_stats['N'] / len(wlvalues))
        df = self.assign_stats_to_table(df, 'GLOBAL')

        for iwl,wl_value in enumerate(wlvalues):
            self.set_data_scatterplot(None,None,-1,wl_value,options_figure)
            self.compute_statistics(False, False, 'II')
            df = self.assign_stats_to_table(df, col_names[iwl+1])
        return df

    def get_table_match_up_statistics(self,options_figure):
        stat_list = list(self.valid_stats.keys())
        col_names = ['INDEX_MU']+stat_list+['SAM','CHI-SQUARE']
        df = pd.DataFrame(index=np.arange(self.mrfile.n_mu_total),columns=col_names)
        df['INDEX_MU'][:]=-1

        for index_mu in range(self.mrfile.n_mu_total):
            stats = self.compute_statistics_spectra(index_mu,options_figure)
            if stats is not None:
                #print(stats['INDEX_MU'],stats['N'],stats['SAM'],stats['RMSD'],stats['BIAS'],stats['CHI-SQUARE'])
                df = self.assign_stats_row_to_table(df,index_mu,stats)

        df = df[df['INDEX_MU']>=0][:]
        return df

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
    
    def get_option_from_file_config(self,file_config):
        if not os.path.isfile(file_config):
            print(f'[ERROR] File config {file_config} is not a valid config file')
            return None
        import configparser
        try:
            options = configparser.ConfigParser()
            options.read(file_config)
            return options
        except:
            print(f'[ERROR] Error reading file_config: {file_config}')
            return None

    def plot_from_options_file(self, file_config):
        options = self.get_option_from_file_config(file_config)
        if options is not None:
            self.plot_from_options(options)

    def plot_figure_from_options_file(self,file_config,figure):
        options = self.get_option_from_file_config(file_config)
        if options is not None:
            poptions = PlotOptions(options, None)
            poptions.set_global_options()
            self.global_options = poptions.global_options
            self.mu_valid_variable = self.global_options['mu_valid_variable']
            list_figures = poptions.get_list_figures()
            if figure not in list_figures:
                print(f'[ERROR] Figure {figure} is not available in the figure list')
                return
            list_virtual = poptions.get_list_virtual_flags()

            print('------------------------------------------------------------------------------------------')
            print(f'[INFO] Starting figure: {figure}')
            options_figure = poptions.get_options(figure)

            if options_figure is None:
                return
            if 'selectBy' in options_figure and options_figure['selectBy'] is not None:
                if options_figure['selectBy'] in list_virtual:
                    self.create_virtual_flag(poptions, options_figure['selectBy'])
                options_figure = self.check_gs_options_impl(options_figure, 'selectBy', 'selectType', 'selectValues')
            if 'groupBy' in options_figure and options_figure['groupBy'] is not None:
                if options_figure['groupBy'] in list_virtual:
                    self.create_virtual_flag(poptions, options_figure['groupBy'])
                options_figure = self.check_gs_options_impl(options_figure, 'groupBy', 'groupType', 'groupValues')
            self.plot_from_options_impl(options_figure)

    def plot_from_options(self, options):
        poptions = PlotOptions(options, None)
        poptions.set_global_options()
        self.global_options = poptions.global_options
        self.mu_valid_variable = self.global_options['mu_valid_variable']

        list_figures = poptions.get_list_figures()
        list_virtual = poptions.get_list_virtual_flags()

        for figure in list_figures:
            print('------------------------------------------------------------------------------------------')
            print(f'[INFO] Starting figure: {figure}')
            options_figure = poptions.get_options(figure)
            if options_figure is None:
                continue
            if 'selectBy' in options_figure and options_figure['selectBy'] is not None:
                if options_figure['selectBy'] in list_virtual:
                    self.create_virtual_flag(poptions, options_figure['selectBy'])
                options_figure = self.check_gs_options_impl(options_figure, 'selectBy', 'selectType', 'selectValues')
            if 'groupBy' in options_figure and options_figure['groupBy'] is not None:
                if options_figure['groupBy'] in list_virtual:
                    self.create_virtual_flag(poptions, options_figure['groupBy'])
                options_figure = self.check_gs_options_impl(options_figure, 'groupBy', 'groupType', 'groupValues')
            self.plot_from_options_impl(options_figure)

    def plot_from_options_impl(self, options_figure):
        if options_figure['type'] == 'scatterplot':
            self.plot_scatterplot_from_options(options_figure)
        if options_figure['type'] == 'spectraplot':
            self.plot_spectraplot_from_options(options_figure)
        if options_figure['type'] == 'timeseries':
            # self.plot_time_series(options_figure)
            if options_figure['groupBy'] is not None:
                self.plot_time_series_grouped(options_figure)
        if options_figure['type'] == 'mapplot':
            self.plot_map_plot(options_figure)
        if options_figure['type'] == 'imageplot':
            self.plot_image(options_figure)
        if options_figure['type'] == 'flagplot':
            self.plot_flag_plot_from_options(options_figure)
        if options_figure['type'] == 'singlestatstable':
            self.plot_single_stats_table(options_figure)
        if options_figure['type'] == 'spectralstatstable':
            self.plot_spectral_stats_table(options_figure)
        if options_figure['type'] == 'matchupsstatstable':
            self.plot_matchups_stats_table(options_figure)
        if options_figure['type'] == 'multiplestatsplot':
            self.plot_multiple_stats_plot(options_figure)
        if options_figure['type'] == 'multipleboundingbox':
            self.plot_multiple_bounding_box(options_figure)
        if options_figure['type'] == 'spectraparam':
            self.plot_spectra_params(options_figure)
        if options_figure['type'] == 'multipleplot':
            self.plot_multiple_plot(options_figure)

    def plot_multiple_plot(self,options_figure):
        print(options_figure)
        from PlotMultiple import PlotMultiple
        rc = options_figure['multiple_plot'].split(',')
        list_files = options_figure['multiple_files']
        file_out = options_figure['file_out']
        files_multiple = []
        for name in list_files:
            if os.path.exists(name):
                files_multiple.append(name)
            else:
                ext = os.path.basename(file_out)[os.path.basename(file_out).index('.'):]
                file_here = os.path.join(os.path.dirname(file_out),f'{name}{ext}')
                print(file_here)
                if os.path.exists(file_here):
                    files_multiple.append(file_here)

        nrow = int(rc[0].strip())
        ncol = int(rc[1].strip())
        ntot = nrow * ncol

        # print(ntot, len(files_multiple), file_out)
        # for idx, file in enumerate(files_multiple):
        #     print(idx, file)

        if ntot != len(files_multiple):
            return
        pm = PlotMultiple()
        xfigsize = options_figure['xfigsize']
        yfigsize = options_figure['yfigsize']
        wspace = options_figure['widthspace']
        hspace = options_figure['heightspace']
        pm.start_multiple_plot_advanced(nrow, ncol, xfigsize, yfigsize, wspace, hspace, True)
        index = 0
        for irow in range(nrow):
            for icol in range(ncol):
                pm.plot_image(files_multiple[index], irow, icol)
                index = index + 1

        #pm.plot_color_bar()
        ##annotations
        if not options_figure['anot_y_axis_'] is None:
            anots = options_figure['anot_y_axis_']
            style = options_figure['anot_default_style']
            fontsize = int(style[0])
            for anot in anots:
                anot_info = anots[anot]
                print(anot_info)
                xpos = float(anot_info[0])
                ypos_list = [float(x) for x in anot_info[1].split(';')]
                print(anot_info[2],type(anot_info[2]))
                text_list = [x.strip() for x in anot_info[2].split(';')]
                for ypos,text in zip(ypos_list,text_list):
                    print(xpos, ypos, text, fontsize)
                    pm.set_text_size(xpos, ypos, text, fontsize)

        if not options_figure['anot_'] is None:
            anots = options_figure['anot_']
            style = options_figure['anot_default_style']
            fontsize = int(style[0])
            for anot in anots:
                anot_info = anots[anot]
                ypos = float(anot_info[0])
                xpos = float(anot_info[1])
                text = anot_info[2]
                print(xpos, ypos, text, fontsize)
                pm.set_text_size(xpos, ypos, text, fontsize)

        pm.save_fig(file_out)

    def plot_spectra_params(self, options_figure):

        wlvalues = options_figure['wlvalues']
        if wlvalues is None:
            wlvalues = list(np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])))

        if options_figure['wl_min'] is not None:
            wlvalues = np.array(wlvalues)
            wlvalues = wlvalues[wlvalues>=options_figure['wl_min']]
        if options_figure['wl_max'] is not None:
            wlvalues = np.array(wlvalues)
            wlvalues = wlvalues[wlvalues<=options_figure['wl_max']]

        self.mrfile.var_mu_valid = self.mu_valid_variable
        stat_list = options_figure['stat_list']

        pspectra = None
        wl_stats = {}

        if options_figure['groupBy'] is None:
            for wl in wlvalues:
                self.set_data_scatterplot(None, None, None, wl, options_figure)
                self.compute_statistics(False, False, 'II')
                wl_stats[wl] = self.valid_stats.copy()

            for istat, stat in enumerate(stat_list):
                pspectra = self.start_spectra_params(wl_stats)
                self.add_line_to_spectra_params(options_figure, pspectra, stat, wl_stats, 0)
                self.close_spectra_params(options_figure, pspectra, stat, wl_stats)

        else:
            flag_name = options_figure['groupBy']
            for stat in stat_list:
                pspectra = None
                gflags = []
                for igroup,group in enumerate(options_figure['groupValues']):
                    wl_stats = {}.copy()
                    for wl in wlvalues:
                        self.set_data_scatterplot(None, options_figure['groupBy'], group, wl, options_figure)
                        self.compute_statistics(False, False, 'II')
                        wl_stats[wl] = self.valid_stats.copy()
                    if igroup==0:
                        pspectra = self.start_spectra_params(wl_stats)
                    gflags.append(self.get_flag_flag(group,options_figure[flag_name]['flag_values'], options_figure[flag_name]['flag_meanings']))
                    self.add_line_to_spectra_params(options_figure,pspectra,stat,wl_stats,igroup)
                if options_figure['legend'] and options_figure['legend_values'] is None:
                    options_figure['legend_values'] = gflags
                self.close_spectra_params(options_figure,pspectra,stat,wl_stats)




    def start_spectra_params(self,wl_stats):
        wl_values = list(wl_stats.keys())
        from PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        pspectra.xdata = wl_values
        return pspectra

    def add_line_to_spectra_params(self,options_figure, pspectra, stat, wl_stats, iseries):
        wl_values = list(wl_stats.keys())
        ydata = []
        for wl in wl_values:
            stat_ref = defaults.valid_stats[stat]['name']
            ydata.append(wl_stats[wl][stat_ref])

        colors = options_figure['color']
        color = colors[iseries] if iseries<len(colors) else colors[0]
        markers = options_figure['marker']
        marker  = markers[iseries] if iseries<len(markers) else markers[0]
        markersizes = options_figure['markersize']
        markersize = markersizes[iseries] if iseries < len(markersizes) else markersizes[0]
        linestyles = options_figure['linestyle']
        linestyle = linestyles[iseries] if iseries < len(linestyles) else linestyles[0]
        linewidths = options_figure['linewidth']
        linewidth = linewidths[iseries] if iseries < len(linewidths) else linewidths[0]
        markercolors = options_figure['markeredgecolor']
        markercolor = markercolors[iseries] if iseries< len(markercolors) else markercolors[0]

        #pspectra.plot_single_line(ydata, color, linestyle, linewidth, marker,markersize)
        pspectra.plot_single_linev2(ydata, color, linestyle, linewidth, marker,markersize,markercolor)

    def close_spectra_params(self,options_figure, pspectra, stat, wl_stats):
        wl_values = list(wl_stats.keys())
        wl_col = [self.get_wl_str_from_wl(x) for x in wl_values]

        xticks_size = 12
        if len(wl_col) >= 15:
            xticks_size = 10

        if len(wl_col) > 100:
            xdata_plot = np.array([350, 400, 450, 500, 550, 600, 650, 700])
            wl_col = [f'{x}' for x in xdata_plot]
            pspectra.set_xticks(xdata_plot, wl_col, 0, xticks_size)
        else:
            pspectra.set_xticks(wl_values, wl_col, 90, xticks_size)

        pspectra.set_grid()

        file_out = options_figure['file_out']
        if options_figure['ylabel'] is not None:
            pspectra.set_yaxis_title(options_figure['ylabel'])
        else:
            pspectra.set_yaxis_title(stat)
        pspectra.set_xaxis_title('Wavelength (nm)')

        if options_figure['legend']:
            pspectra.legend_options['markerscale']=2
            pspectra.legend_options['loc']='lower center'
            pspectra.legend_options['ncols'] = len(options_figure['legend_values'])
            pos = (0.50,-0.40)
            pspectra.legend_options['bbox_to_anchor'] = pos
            pspectra.set_legend(options_figure['legend_values'])

        #pspectra.prepare_poster()
        pspectra.set_tigth_layout()



        if file_out is not None:
            file_out = f'{file_out[:file_out.rindex(".")]}_{stat}{file_out[file_out.rindex("."):]}'
            pspectra.save_plot(file_out)


    def plot_spectra_params_impl(self, options_figure, pspectra, stat, wl_stats, iseries):
        wl_values = list(wl_stats.keys())
        if pspectra is None:
            from PlotSpectra import PlotSpectra
            pspectra = PlotSpectra()
            pspectra.xdata = wl_values

        ydata = []
        for wl in wl_values:
            stat_ref = defaults.valid_stats[stat]['name']
            ydata.append(wl_stats[wl][stat_ref])

        color = options_figure['color']
        marker = options_figure['marker']
        marker_size = options_figure['markersize']
        line_type = options_figure['linestyle']
        line_width = options_figure['linewidth']
        pspectra.plot_single_line(ydata, color[iseries], line_type[iseries], line_width[iseries], marker[iseries],
                                  marker_size[iseries])

        wl_col = [self.get_wl_str_from_wl(x) for x in wl_values]

        xticks_size = 12
        if len(wl_col) >= 15:
            xticks_size = 10

        if len(wl_col) > 100:
            xdata_plot = np.array([350, 400, 450, 500, 550, 600, 650, 700])
            wl_col = [f'{x}' for x in xdata_plot]
            pspectra.set_xticks(xdata_plot, wl_col, 0, xticks_size)
        else:
            pspectra.set_xticks(wl_values, wl_col, 90, xticks_size)

        pspectra.set_grid()

        file_out = options_figure['file_out']

        pspectra.set_yaxis_title(stat)
        pspectra.set_xaxis_title('Wavelength (nm)')
        pspectra.set_tigth_layout()
        if file_out is not None:
            file_out = f'{file_out[:file_out.rindex(".")]}_{stat}{file_out[file_out.rindex("."):]}'
            pspectra.save_plot(file_out)

        return pspectra

    def plot_multiple_bounding_box(self, options_figure):
        from matplotlib import pyplot as plt
        ##dimensions
        vars = options_figure['vars']
        groupsTicks = options_figure['groupsTicks']
        ##dimensions
        nseries = len(vars)
        ngroups = len(groupsTicks)
        width = 0.6
        intra = 0.1
        inter = 0.6

        all_pos_series = []

        width_series = round(nseries * (width + intra), 2)
        increm = round(width_series + inter, 2)
        for idx in range(nseries):
            start = inter + ((idx + 1) * (width + intra)) - (width / 2)
            start = round(start, 2)
            end = start + (increm * (ngroups - 1))
            end = round(end, 2)
            pos_series = np.arange(start, end + (increm / 2), increm)
            max_value = end + (width / 2) + inter + intra
            all_pos_series.append(pos_series)

        ticks = [0] * ngroups
        pos_ini = all_pos_series[0]
        pos_end = all_pos_series[-1]

        for idx in range(ngroups):
            ticks[idx] = pos_ini[idx] + ((pos_end[idx] - pos_ini[idx]) / 2)

        max_value = round(max_value, 2)

        ##getting data and plotting
        handles = []
        colors = options_figure['color']

        plt.figure(figsize=(8, 5.5))

        for idx, var in enumerate(vars):
            series_name = var
            array = self.mrfile.nc.variables[var][:]
            data_boxplot = []
            if len(colors) == 1:
                color = colors[0]
            elif len(colors) == nseries:
                color = colors[idx]

            for igroup, group in enumerate(options_figure['groups']):
                values = [x.strip() for x in options_figure['groupsValues'][igroup].split(';')]
                options_figure['selectBy'] = group
                options_figure['selectValues'] = values
                options_figure['selectType'] = 'flag'
                options_figure = self.check_gs_options_impl(options_figure, 'selectBy', 'selectType', 'selectValues')
                select_array, all_select_values, all_select_meanings = self.get_flag_array(options_figure, 'selectBy')
                for val_s in values:
                    val_n = all_select_values[all_select_meanings.index(val_s)]
                    data_here = array[select_array == val_n]
                    data_boxplot.append(data_here)

            positions = all_pos_series[idx]

            bbox = plt.boxplot(data_boxplot, positions=positions, widths=width, patch_artist=True, showmeans=True,
                               boxprops=dict(facecolor=color),
                               medianprops=dict(color='black'),
                               meanprops=dict(linewidth=0, marker='o', mec='black', mfc=color, markersize=4, mew=0.5),
                               flierprops=dict(markersize=0, markeredgecolor='none'))

            handles.append(bbox['boxes'][0])

            # hlegend = define_box_properties(bbox, series[series_name], series_name)
            # handles.append(hlegend)
            # if plot_average:
            #     for itime in range(n_time):
            #         plt.plot(positions[itime], data_average[itime], color=series[series_name], linewidth=0, marker='o',
            #                  markersize=8, mec='black', mew=1)

        plt.grid(which='major', color='lightgray', linestyle='--', axis='y')
        groupsTicks = options_figure['groupsTicks']
        plt.xticks(ticks, groupsTicks, fontsize=12)
        plt.yticks(fontsize=12)
        if options_figure['legend_values'] is not None:
            legend_values = options_figure['legend_values']
        else:
            legend_values = vars
        plt.legend(handles, legend_values, fontsize=12,loc='lower center', ncols=nseries, bbox_to_anchor=(0.5, -0.2))
        if options_figure['y_min'] is not None and options_figure['y_max'] is not None:
            plt.ylim([options_figure['y_min'], options_figure['y_max']])
        if options_figure['ylabel'] is not None:
            y_label = options_figure['ylabel']
            plt.ylabel(y_label, fontsize=12)

        plt.tight_layout()

        file_out = options_figure['file_out']
        if file_out is not None:
            plt.savefig(file_out, dpi=300)
        plt.close()

    def plot_spectral_stats_table(self,options_figure):
        df = self.get_table_spectral_statistics(options_figure)
        if options_figure['file_out'] is not None:
            file_out = options_figure['file_out']
            ext = file_out[file_out.rfind('.'):]
            file_out = file_out.replace(ext,'.csv')
            df.to_csv(file_out,sep=';')


    def plot_matchups_stats_table(self,options_figure):
        df = self.get_table_match_up_statistics(options_figure)
        if options_figure['file_out'] is not None:
            file_out = options_figure['file_out']
            ext = file_out[file_out.rfind('.'):]
            file_out = file_out.replace(ext,'.csv')
            df.to_csv(file_out,sep=';',index=None)

    def plot_single_stats_table(self, options_figure):
        yvar = options_figure['yvar']
        xvar = options_figure['xvar']
        if yvar is None or xvar is None:
            print(f'[ERROR] Please select xvar and yvar options for single stats table {options_figure["name"]}')
            return
        file_out = options_figure['file_out']
        if file_out is not None:
            name_out = os.path.basename(file_out)
            ext = name_out[name_out.find('.'):]
            name_out = name_out.replace(ext, '.csv')
            file_out = os.path.join(os.path.dirname(file_out), name_out)

        if options_figure['selectBy'] is None:
            self.set_data_scatterplot_general(None, None, None, options_figure)
            self.compute_statistics(False, False, 'II')
            col_names = ['GLOBAL']
            df = pd.DataFrame(index=list(self.valid_stats.keys()), columns=col_names)
            df = self.assign_stats_to_table(df, 'GLOBAL')
            df.to_csv(file_out, sep=';', index_label='METRIC')
        else:
            selectBy = options_figure['selectBy']
            selectValues = options_figure['selectValues']
            flagValues = self.get_flag_list(selectValues, options_figure[selectBy]['flag_values'],
                                            options_figure[selectBy]['flag_meanings'])
            col_names = ['GLOBAL'] + flagValues

            #print(self.valid_stats)
            df = pd.DataFrame(index=list(self.valid_stats.keys()), columns=col_names)

            self.set_data_scatterplot_general(None, None, None, options_figure)
            self.compute_statistics(True, False, 'II')
            df = self.assign_stats_to_table(df, 'GLOBAL')

            for idx, val in enumerate(selectValues):
                col_name = flagValues[idx]
                self.set_data_scatterplot_general(None, selectBy, val, options_figure)
                self.compute_statistics(True, False, 'II')
                df = self.assign_stats_to_table(df, col_name)

            df.to_csv(file_out, sep=';', index_label='METRIC')

    def assign_stats_to_table(self, df, col_name):
        for key in list(self.valid_stats.keys()):
            df.loc[key].at[col_name] = self.valid_stats[key]
        return df

    def assign_stats_row_to_table(self,df,index_mu,stats):
        for key in list(stats.keys()):
            df.loc[index_mu,key] = stats[key]
        return df

    def plot_flag_plot_from_options(self, options_figure):

        if options_figure['xlabel'] is None:
            options_figure['xlabel'] = '# of pixels'
        if options_figure['ylabel'] is None:
            options_figure['ylabel'] = 'Flag'

        if options_figure['type_flagplot'] == 'single':  ##flag plots for each match-up
            index_mu = options_figure['index_mu']
            if index_mu == -1:
                for imu in range(self.mrfile.n_mu_total):
                    print(f'[INFO] Plotting single flag plot for match-up {imu} / {self.mrfile.n_mu_total}')
                    self.plot_single_flag_plot(options_figure, imu)
            elif 0 <= index_mu < self.mrfile.n_mu_total:
                print(f'[INFO] Plotting single flag plot for match-up {index_mu}')
                self.plot_single_flag_plot(options_figure, index_mu)

    def plot_single_flag_plot(self, options_figure, index_mu):
        name_variable = options_figure['var_flag']
        if name_variable is None:
            return
        flag_list = options_figure['flag_list']
        array_flag, flag_values, flag_meanings = self.get_flag_array(options_figure, 'var_flag')
        if flag_list is None:
            flag_list = ' '.join(flag_meanings)
        window_sizes = options_figure['window_sizes']
        n_window = len(window_sizes)
        n_flag = len(flag_list)
        coordinates = np.zeros((n_window, 4), dtype=np.int32)
        series_names = []
        for iwindow in range(n_window):
            if window_sizes[iwindow] == -1:
                coordinates[iwindow, 0] = 0
                coordinates[iwindow, 1] = len(self.mrfile.dimensions['rows'])
                coordinates[iwindow, 2] = 0
                coordinates[iwindow, 3] = len(self.mrfile.dimensions['columns'])
                series_names.append('Complete')
            else:
                self.mrfile.window_size = window_sizes[iwindow]
                central_r, central_c, r_s, r_e, c_s, c_e = self.mrfile.get_dimensions()
                coordinates[iwindow, 0] = r_s
                coordinates[iwindow, 1] = r_e
                coordinates[iwindow, 2] = c_s
                coordinates[iwindow, 3] = c_e
                series_names.append(f'{window_sizes[iwindow]} x {window_sizes[iwindow]}')

        output_data = np.zeros((n_flag, n_window))
        array_flag = array_flag[index_mu]

        for idx, flag in enumerate(flag_list):
            flag_here = {
                '0': {
                    'is_default': False,
                    'flag_list': [flag],
                    'flag_value': 1,
                    'flag_meaning': 'FLAGGED'
                }
            }
            mask, flag_info = self.create_flag_mask(array_flag, options_figure['var_flag'], flag_here)
            for iwindow in range(n_window):
                mask_to_sum = mask[coordinates[iwindow, 0]:coordinates[iwindow, 1],
                              coordinates[iwindow, 2]:coordinates[iwindow, 3]]
                output_data[idx, iwindow] = np.sum(mask_to_sum)

        self.plot_flag_plot_impl(options_figure, output_data, series_names, flag_list, index_mu)

    def plot_flag_plot_impl(self, options_figure, output_data, series_names, flag_names, index_mu):
        nseries = len(series_names)
        nflag = output_data.shape[0]
        # colors = options_out['series_color']
        # if colors is None:
        #     colors = defaults.get_color_list(nseries)
        colors = ['blue', 'green']

        from matplotlib import pyplot as plt
        fig, ax = plt.subplots()
        height = 0.2
        xval = 0
        xticks_pos = []
        xticks_minor_pos = []
        handles = [None] * nseries
        heightbyseries = height / nseries
        for iflag in range(nflag):
            xini = xval
            xticks_minor_pos.append(xini - (heightbyseries / 2))
            for idx in range(nseries):
                hbar = plt.barh(xval, output_data[iflag, idx], height=heightbyseries, color=colors[idx])
                handles[idx] = hbar
                xval = xval + heightbyseries
            xfin = xval - heightbyseries
            xticks_minor_pos.append(xfin + (heightbyseries / 2))
            xoutput = (xini + xfin) / 2
            xticks_pos.append(xoutput)

        plt.ylabel(options_figure['ylabel'], fontsize=12)
        plt.xlabel(options_figure['xlabel'], fontsize=12)
        plt.yticks(xticks_pos, flag_names)
        ax.set_yticks(xticks_minor_pos, minor=True)
        plt.grid(which='minor', color='gray', linestyle='--', axis='y')
        plt.grid(which='major', color='gray', linestyle='--', axis='x')
        ax.tick_params(which='major', length=0, axis='y')
        ax.tick_params(which='minor', length=10, axis='y')
        if options_figure['legend']:
            legend_values = series_names if options_figure['legend_values'] is None else options_figure['legend_values']
            plt.legend(handles, legend_values, framealpha=1, loc='lower right')

        plt.tight_layout()

        if not options_figure['file_out'] is None:
            file_out = options_figure['file_out']
            if index_mu is not None:
                satellite_time = self.mrfile.sat_times[index_mu].strftime('%Y%m%d')
                file_out = f'{file_out[:-4]}_{satellite_time}_{index_mu}{file_out[-4:]}'
            plt.savefig(file_out, dpi=300)

        plt.close(fig)

    def plot_map_plot(self, options_figure):
        import cartopy
        import cartopy.crs as ccrs
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mticker

        if not self.VALID and not os.path.isfile(self.mrfile.file_path):
            print(f'[ERROR] {self.mrfile.file_path} shoud be a valid NetCDF file')
            return

        print(f'[INFO] [mapplot PLOT] Getting arrays...')
        latitude_var = self.mrfile.nc.variables[options_figure['latitude_variable']]
        longitude_var = self.mrfile.nc.variables[options_figure['longitude_variable']]
        if latitude_var.dimensions != longitude_var.dimensions:
            print(
                f'[ERROR] Latitude ({options_figure["latitude_variable"]}) and longitude ({options_figure["longitude_variable"]}) variable shoud have the same dimensions')
            return
        lat_array = self.mrfile.get_full_array_1D(options_figure['latitude_variable'],
                                                  options_figure['insitu_id_variable'], False)
        lon_array = self.mrfile.get_full_array_1D(options_figure['longitude_variable'],
                                                  options_figure['insitu_id_variable'], False)
        valid_array = np.ones(lat_array.shape)
        if options_figure['valid_variable_masked'] is not None:
            for mask_variable in options_figure['valid_variable_masked']:
                marray = self.mrfile.get_full_array_1D(mask_variable, options_figure['insitu_id_variable'], False)
                valid_array[marray.mask] = 0
        valid_lat_lon = np.ones(lat_array.shape)
        valid_lat_lon[lat_array.mask] = 0
        valid_lat_lon[lon_array.mask] = 0
        lat_array = lat_array[valid_lat_lon == 1]
        lon_array = lon_array[valid_lat_lon == 1]
        valid_array = valid_array[valid_lat_lon == 1]
        if options_figure['limit_to_valid']:
            lat_array = lat_array[valid_array == 1]
            lon_array = lon_array[valid_array == 1]
            valid_array = valid_array[valid_array == 1]

        print(f'[INFO] [mapplot PLOT] Getting geographical limits...')
        geo_limits = self.get_geo_limits(options_figure, lat_array, lon_array)
        extent = (geo_limits[2], geo_limits[3], geo_limits[0], geo_limits[1])

        print(
            f'[INFO] [mapplot PLOT] Maps limits: Latitude-> {geo_limits[0]} to {geo_limits[1]}; Longitude: {geo_limits[2]} to {geo_limits[3]}')

        print(f'[INFO] [mapplot PLOT] Plotting...')
        ax = plt.axes(projection=ccrs.PlateCarree(), extent=extent)

        # # ax.coastlines(linewidth=0.5)
        ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black', linewidth=0.5)

        gl = ax.gridlines(linewidth=0.5, linestyle='dotted', draw_labels=True)
        gl.xlabels_top = False
        gl.ylabels_right = False

        # # lon_labels = [-15,-10,-5,0,5,10,15,20,25]
        # # gl.xlocator = mticker.FixedLocator(lon_labels)
        #
        # # ax.set_yticks(ax.get_yticks())
        # # ax.set_xticks(ax.get_xticks())
        #
        if options_figure['groupBy'] is None:
            if options_figure['groupByValid']:
                lat_array_valid = lat_array[valid_array == 1]
                lon_array_valid = lon_array[valid_array == 1]
                lat_array_invalid = lat_array[valid_array == 0]
                lon_array_invalid = lon_array[valid_array == 0]
                if len(lat_array_valid) > 0 and len(lon_array_valid) > 0:
                    plt.plot(lon_array_valid.tolist(), lat_array_valid.tolist(),
                             color=options_figure['valid_style']['color'][0],
                             marker=options_figure['valid_style']['marker'][0],
                             markersize=options_figure['valid_style']['markersize'][0],
                             linestyle=options_figure['valid_style']['linestyle'][0],
                             linewidth=options_figure['valid_style']['linewidth'][0])
                if len(lat_array_invalid) > 0 and len(lon_array_invalid) > 0:
                    plt.plot(lon_array_invalid.tolist(), lat_array_invalid.tolist(),
                             color=options_figure['default_style']['color'][0],
                             marker=options_figure['default_style']['marker'][0],
                             markersize=options_figure['default_style']['markersize'][0],
                             linestyle=options_figure['default_style']['linestyle'][0],
                             linewidth=options_figure['default_style']['linewidth'][0])
            else:
                style = 'default_style'
                if options_figure['limit_to_valid']:
                    style = 'valid_style'
                plt.plot(lon_array.tolist(), lat_array.tolist(),
                         color=options_figure[style]['color'][0],
                         marker=options_figure[style]['marker'][0],
                         markersize=options_figure[style]['markersize'][0],
                         linestyle=options_figure[style]['linestyle'][0],
                         linewidth=options_figure[style]['linewidth'][0])
        else:
            pass
            # all_group_array, all_group_values, all_group_meanings = self.get_flag_array(options_out, 'groupBy')
            # group_values = options_out['groupValues']
        #
        #     print(options_out)
        #     if 'groupArraySelect' in options_out:
        #         groupArray = options_out['groupArraySelect']
        #     else:
        #         groupArray = all_group_array
        #
        #     ngroup = len(group_values)
        #     for idx in range(ngroup):
        #         gvalue = group_values[idx]
        #         ghere = self.get_flag_flag(gvalue, np.array(all_group_values), all_group_meanings)
        #         color = self.get_option_from_list(options_out['point_color'], idx, ngroup)
        #         if len(options_out['point_color']) == 1 and ngroup > 1:
        #             color = defaults.get_color_default(idx, 0, ngroup - 1)
        #         size = self.get_option_from_list(options_out['point_size'], idx, ngroup)
        #         print(f'[INFO] Plotting group: {ghere} with value: {gvalue} Color: {color}')
        #         array_lon_here = array_lon[groupArray == group_values[idx]]
        #         array_lat_here = array_lat[groupArray == group_values[idx]]
        #         plt.plot(array_lon_here.tolist(), array_lat_here.tolist(), color=color, marker='o', markersize=size,
        #                  linestyle='none')
        if options_figure['title'] is not None:
            plt.title(options_figure['title'])

        file_out = options_figure['file_out']
        if file_out is not None:
            if file_out.endswith('.tif'):
                plt.savefig(file_out, dpi=300, bbox_inches='tight', pil_kwargs={"compression": "tiff_lzw"})
            else:
                plt.savefig(file_out, dpi=300, bbox_inches='tight')

        # ax.close()
        plt.close()
        print(f'[INFO] [mapplot PLOT] Completed')

    def get_box(self, array_lat, array_lon, n):
        size = array_lat.shape[0]
        center = int(np.floor(size / 2))
        if n == -1:  # external box:
            increm = center
        elif n == 0:  # central pixel
            increm = 0
        else:  ##n shoud be even number
            increm = int(np.floor(n / 2))
        ini = center - increm
        end = center + increm
        if ini < 0: ini = 0
        if end >= size: end = size - 1



        lat_0_0 = array_lat[ini, ini]
        increm_lat_0_0 = abs((array_lat[ini + 1, ini + 1] - array_lat[ini, ini]) / 2)
        lat_0_n = array_lat[ini, end]
        increm_lat_0_n = abs((array_lat[ini + 1, end - 1] - array_lat[ini, end]) / 2)
        lat_n_0 = array_lat[end, ini]
        increm_lat_n_0 = abs((array_lat[end - 1, ini + 1] - array_lat[end, ini]) / 2)
        lat_n_n = array_lat[end, end]
        increm_lat_n_n = abs((array_lat[end - 1, end - 1] - array_lat[end, end]) / 2)

        lon_0_0 = array_lon[ini, ini]
        increm_lon_0_0 = abs((array_lon[ini + 1, ini + 1] - array_lon[ini, ini]) / 2)
        lon_0_n = array_lon[ini, end]
        increm_lon_0_n = abs((array_lon[ini + 1, end - 1] - array_lon[ini, end]) / 2)
        lon_n_0 = array_lon[end, ini]
        increm_lon_n_0 = abs((array_lon[end - 1, ini + 1] - array_lon[end, ini]) / 2)
        lon_n_n = array_lon[end, end]
        increm_lon_n_n = abs((array_lon[end - 1, end - 1] - array_lon[end, end]) / 2)

        ##from south to north
        if array_lat[0, center] < array_lat[size - 1, center]:
            lat_points = [lat_0_0 - increm_lat_0_0, lat_0_n - increm_lat_0_n, lat_n_n + increm_lat_n_n,
                          lat_n_0 + increm_lat_n_0, lat_0_0 - increm_lat_0_0]
        else:  # from north to south
            lat_points = [lat_0_0 + increm_lat_0_0, lat_0_n + increm_lat_0_n, lat_n_n - increm_lat_n_n,
                          lat_n_0 - increm_lat_n_0, lat_0_0 + increm_lat_0_0]

        ##from west to east
        if array_lon[center, 0] < array_lon[center, size - 1]:

            lon_points = [lon_0_0 - increm_lon_0_0, lon_0_n + increm_lon_0_n, lon_n_n + increm_lon_n_n,
                          lon_n_0 - increm_lon_n_0, lon_0_0 - increm_lon_0_0]
        else:  # from east to west

            lon_points = [lon_0_0 + increm_lon_0_0, lon_0_n - increm_lon_0_n, lon_n_n - increm_lon_n_n,
                          lon_n_0 + increm_lon_n_0, lon_0_0 + increm_lon_0_0]

        return lat_points, lon_points

    def get_geo_limits(self, options_figure, array_lat, array_lon):
        geo_limits = options_figure['geo_limits'] if 'geo_limits' in options_figure else None
        if geo_limits is None:
            min_lat = np.min(array_lat)
            max_lat = np.max(array_lat)
            min_lon = np.min(array_lon)
            max_lon = np.max(array_lon)
            # if options_out['plot_extracts'] is True:
            #     extract_list = options_out['plot_extract_list']
            #     if extract_list is None:
            #         extract_list = self.get_extract_list(dt.strptime('20230610', '%Y%m%d'))
            #         options_out['plot_extract_list'] = extract_list
            #
            #     for iextract in extract_list:
            #         lat_line_all, lon_line_all = self.get_polygon_extract(iextract, -1)
            #         min_lat_extract = np.min(np.array(lat_line_all))
            #         max_lat_extract = np.max(np.array(lat_line_all))
            #         min_lon_extract = np.min(np.array(lon_line_all))
            #         max_lon_extract = np.max(np.array(lon_line_all))
            #         if min_lat_extract < min_lat:
            #             min_lat = min_lat_extract
            #         if max_lat_extract > max_lat:
            #             max_lat = max_lat_extract
            #         if min_lon_extract < min_lon:
            #             min_lon = min_lon_extract
            #         if max_lon_extract > max_lon:
            #             max_lon = max_lon_extract

            # if abs(max_lat - min_lat) > 1:
            #     min_lat = np.floor(min_lat)
            #     max_lat = np.ceil(max_lat)
            # else:
            #     min_lat = min_lat - 0.01
            #     max_lat = max_lat + 0.01
            # if abs(max_lon - min_lon) > 1:
            #     min_lon = np.floor(min_lon)
            #     max_lon = np.ceil(max_lon)
            # else:
            #     min_lon = min_lon - 0.01
            #     max_lon = max_lon + 0.01
            # geo_limits = [float(min_lat), float(max_lat), float(min_lon), float(max_lon)]
            lat_box, lon_box = self.get_box(array_lat, array_lon, -1)
            min_lat = np.min(lat_box)
            max_lat = np.max(lat_box)
            min_lon = np.min(lon_box)
            max_lon = np.max(lon_box)
            geo_limits = [float(min_lat), float(max_lat), float(min_lon), float(max_lon)]

        return geo_limits

    def check_variable(self, config_ref_variable, name_variable):
        if name_variable is None:
            if config_ref_variable is not None:
                print(f'[ERROR] {config_ref_variable} should be defined in the configuration file')
            return False
        if name_variable not in self.mrfile.nc.variables:
            print(f'[ERROR] {name_variable} is not defined in {self.mrfile.file_path}')
            return False

    def plot_image(self, options_figure):
        if not self.VALID and not os.path.isfile(self.mrfile.file_path):
            print(f'[ERROR] {self.mrfile.file_path} shoud be a valid NetCDF file')
            return
        variables_to_check = ['latitude_variable', 'longitude_variable', 'plot_variable']
        for variable_to_check in variables_to_check:
            if self.check_variable(variables_to_check, options_figure[variable_to_check]):
                return
        index_mu = options_figure['index_mu']

        if index_mu == -1:
            for imu in range(self.mrfile.n_mu_total):
                if options_figure['apply_geo']:
                    self.plot_geo_image_impl(options_figure, imu)
                else:
                    self.plot_image_impl(options_figure, imu)
        else:
            if options_figure['apply_geo']:
                self.plot_geo_image_impl(options_figure, index_mu)
            else:
                self.plot_image_impl(options_figure, index_mu)

    def plot_image_impl(self, options_figure, index_mu):
        lat_array = self.mrfile.nc.variables[options_figure['latitude_variable']][index_mu]
        lon_array = self.mrfile.nc.variables[options_figure['longitude_variable']][index_mu]
        if self.mrfile.nc.variables[options_figure['plot_variable']].ndim == 4:
            index_band = options_figure['index_band']
            data_array = self.mrfile.nc.variables[options_figure['plot_variable']][index_mu, index_band, :, :]
        elif self.mrfile.nc.variables[options_figure['plot_variable']].ndim == 3:
            data_array = self.mrfile.nc.variables[options_figure['plot_variable']][index_mu, :, :]

        flag_info = None
        if options_figure['create_flag_mask']:
            data_array, flag_info = self.create_flag_mask(data_array, options_figure['plot_variable'],
                                                          options_figure['flag_'])

        if flag_info is not None:
            values_flag_info = list(flag_info.values())
            legend = list(flag_info.keys())
            values_flag_info.append(values_flag_info[-1] + 1)
            from matplotlib.colors import BoundaryNorm
            from matplotlib.colors import ListedColormap
            bnorm = BoundaryNorm(values_flag_info, ncolors=3)
            cmap = ListedColormap(['blue', 'beige', 'red'])

        import matplotlib.pyplot as plt
        plt.pcolormesh(data_array, norm=bnorm, cmap=cmap, edgecolors='dimgray', linewidths=0.01)

        ##box
        center = int(np.floor(data_array.shape[0] / 2))
        ybox = [center - 2, center - 2, center + 3, center + 3, center - 2]
        xbox = [center - 2, center + 3, center + 3, center - 2, center - 2]
        plt.plot(xbox, ybox, color='white', marker='o', markersize=0, linestyle='--', linewidth=1)
        # point
        plt.plot(center + 0.5, center + 0.5, color='white', marker='o', markersize=1)

        cbar = plt.colorbar()
        cbar.ax.get_yaxis().set_ticks([])

        cbar.ax.text(1.6, 1.5, legend[0], ha='center', va='center', rotation=90)
        cbar.ax.text(1.6, 3.0, legend[1], ha='center', va='center', rotation=90)
        cbar.ax.text(1.6, 4.5, legend[2], ha='center', va='center', rotation=90)

        if options_figure['title'] is not None:
            title = options_figure['title']
            title = title.replace('$INDEX_MU$', f'{index_mu}')
            if '$DATE$' in title:
                satellite_time = self.mrfile.sat_times[index_mu].strftime('%Y-%m-%d')
                title = title.replace('$DATE$', satellite_time)

            plt.title(title)

        file_out = options_figure['file_out']

        if file_out is not None:
            satellite_time = self.mrfile.sat_times[index_mu].strftime('%Y%m%d')
            file_out = f'{file_out[:-4]}_{satellite_time}_{index_mu}{file_out[-4:]}'
            if file_out.endswith('.tif'):
                plt.savefig(file_out, dpi=300, bbox_inches='tight', pil_kwargs={"compression": "tiff_lzw"})
            else:
                plt.savefig(file_out, dpi=300, bbox_inches='tight')

        plt.close()

    def plot_geo_image_impl(self, options_figure, index_mu):
        import cartopy
        import cartopy.crs as ccrs
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mticker

        lat_array = self.mrfile.nc.variables[options_figure['latitude_variable']][index_mu]
        lon_array = self.mrfile.nc.variables[options_figure['longitude_variable']][index_mu]
        if self.mrfile.nc.variables[options_figure['plot_variable']].ndim == 4:
            index_band = options_figure['index_band']
            data_array = self.mrfile.nc.variables[options_figure['plot_variable']][index_mu, index_band, :, :]
        elif self.mrfile.nc.variables[options_figure['plot_variable']].ndim == 3:
            data_array = self.mrfile.nc.variables[options_figure['plot_variable']][index_mu, :, :]

        geo_limits = self.get_geo_limits(options_figure, lat_array, lon_array)
        extent = (geo_limits[2], geo_limits[3], geo_limits[0], geo_limits[1])
        projection = ccrs.Mercator()
        ax = plt.axes(projection=projection, extent=extent)

        ##grid lines
        gl = ax.gridlines(linewidth=0.5, linestyle='dotted', draw_labels=True)
        gl.xlabels_top = False
        gl.ylabels_right = False
        # lon_ticks = [11.5,11.75,12,12.25,12.5,12.75,13,13.25]
        # gl.xlocator = mticker.FixedLocator(lon_ticks)

        ##external box
        lat_box, lon_box = self.get_box(lat_array, lon_array, -1)
        plt.plot(lon_box, lat_box, color='k', marker='o', markersize=0, linestyle='-', linewidth=0.5,
                 transform=ccrs.PlateCarree())

        ##5x5 box
        lat_box, lon_box = self.get_box(lat_array, lon_array, 5)
        plt.plot(lon_box, lat_box, color='k', marker='o', markersize=0, linestyle='--', linewidth=0.5,
                 transform=ccrs.PlateCarree())

        # in situ site location
        lat_point = 45.313900
        lon_point = 12.508300
        plt.plot(lon_point, lat_point, color='black', marker='.', markersize=1, transform=ccrs.PlateCarree())

        from matplotlib.colors import LogNorm
        from matplotlib.colors import Normalize
        plt.pcolormesh(lon_array, lat_array, data_array, transform=ccrs.PlateCarree(), cmap='jet',
                       norm=Normalize(vmin=400, vmax=700))
        plt.colorbar(shrink=0.6)

        import cartopy.feature as cfeature
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face',
                                                facecolor=cfeature.COLORS['land'])
        ax.add_feature(land_10m, zorder=0, edgecolor='black', linewidth=0.5)

        if options_figure['title'] is not None:
            title = options_figure['title']
            title = title.replace('$INDEX_MU$', f'{index_mu}')
            if '$DATE$' in title:
                satellite_time = self.mrfile.sat_times[index_mu].strftime('%Y-%m-%d')
                title = title.replace('$DATE$', satellite_time)

            plt.title(title)

        plt.tight_layout()

        file_out = options_figure['file_out']

        if file_out is not None:
            satellite_time = self.mrfile.sat_times[index_mu].strftime('%Y%m%d')
            file_out = f'{file_out[:-4]}_{satellite_time}_{index_mu}{file_out[-4:]}'
            if file_out.endswith('.tif'):
                plt.savefig(file_out, dpi=300, bbox_inches='tight', pil_kwargs={"compression": "tiff_lzw"})
            else:
                plt.savefig(file_out, dpi=300, bbox_inches='tight')

        plt.close()

    def create_flag_mask(self, array, name_variable, flags):
        variable = self.mrfile.nc.variables[name_variable]
        array_mask = array
        flag_info = None

        if 'flag_masks' in variable.ncattrs() and 'flag_meanings' in variable.ncattrs():

            all_flag_values = variable.flag_masks
            all_flag_meanings = variable.flag_meanings
            flag_obj = flag_class.Class_Flags_Polymer(all_flag_values, all_flag_meanings)
            array_mask = np.zeros(array.shape)

            flag_info = {}
            if flags is not None:
                default_value = 0
                for iflag in flags:
                    if flags[iflag]['is_default']:
                        default_value = flags[iflag]['flag_value']
                        meaning = flags[iflag]['flag_meaning']
                        flag_info[meaning] = default_value
                        array_mask[:] = default_value
                for iflag in flags:
                    if not flags[iflag]['is_default']:
                        flag_list = flags[iflag]['flag_list']
                        value = flags[iflag]['flag_value']
                        meaning = flags[iflag]['flag_meaning']
                        mask = flag_obj.Mask(array, flag_list)
                        array_mask[np.logical_and(mask > 0, array_mask == default_value)] = value
                        flag_info[meaning] = value

        return array_mask, flag_info

    def plot_time_series_grouped(self, options_figure):
        if not self.VALID and not os.path.isfile(self.mrfile.file_path):
            print(f'[ERROR] {self.mrfile.file_path} shoud be a valid NetCDF file')
            return

        groupVar = options_figure['groupBy']

        if groupVar in self.virtual_flags:
            time_array = self.virtual_flags[groupVar]['flag_array']

        selectVar = options_figure['selectBy']
        if selectVar is not None:
            if selectVar in self.virtual_flags:
                select_array = self.virtual_flags[selectVar]['flag_array']
        else:
            select_array = np.ones(time_array.shape)

        time_array[select_array == 0] = 0

        groupValues = options_figure['groupValues']
        groupValues.sort()

        dataVars = options_figure['data_var']
        col_name_ref = ['_num', '_sum', '_avg', '_median', '_std', '_min', '_max', '_p25', '_p75', '_range', '_iqr']
        col_names = []
        for dataVar in dataVars:
            for ref in col_name_ref:
                col_names.append(f'{dataVar}{ref}')

        dfTime = pd.DataFrame(index=groupValues, columns=col_names)

        # temp = self.mrfile.nc.variables['time'][:]
        # from datetime import datetime as dt
        # index_ini = -1
        # for idx, t in enumerate(temp):
        #     date_temp = dt.utcfromtimestamp(t)
        #     if date_temp.year >= 1998 and index_ini == -1:
        #         index_ini = idx
        #     if date_temp.year >= 2021:
        #         index_end = idx - 1
        #         break
        # print('INDEX INI->',index_ini)
        # print('INDEX_END->',index_end)

        # index_ini = 119 #first index 1998
        # index_end = 8519 # last index 2021
        # time_array = time_array[index_ini:index_end]

        for dataVar in dataVars:
            data_array = self.mrfile.nc.variables[dataVar]

            for val in groupValues:
                print(f'[INFO]-> Computing for {dataVar} / {val}')
                data_array_here = data_array[time_array == val]
                if np.ma.is_masked(data_array_here):
                    data_array_here = data_array_here[~data_array_here.mask]

                num = len(data_array_here)

                dfTime.loc[val, f'{dataVar}_num'] = num
                if num > 0:
                    dfTime.loc[val, f'{dataVar}_sum'] = np.sum(data_array_here)
                    dfTime.loc[val, f'{dataVar}_avg'] = np.mean(data_array_here)
                    dfTime.loc[val, f'{dataVar}_median'] = np.median(data_array_here)
                    dfTime.loc[val, f'{dataVar}_std'] = np.std(data_array_here)
                    dfTime.loc[val, f'{dataVar}_min'] = np.min(data_array_here)
                    dfTime.loc[val, f'{dataVar}_max'] = np.max(data_array_here)
                    dfTime.loc[val, f'{dataVar}_p25'] = np.percentile(data_array_here, 25)
                    dfTime.loc[val, f'{dataVar}_p75'] = np.percentile(data_array_here, 75)
            dfTime[f'{dataVar}_range'] = dfTime[f'{dataVar}_max'] - dfTime[f'{dataVar}_min']
            dfTime[f'{dataVar}_iqr'] = dfTime[f'{dataVar}_p75'] - dfTime[f'{dataVar}_p25']

        if 'file_out' in options_figure:
            file_out = options_figure['file_out']
            index_label = options_figure['index_label']
            file_csv = file_out[0:file_out.find('.')] + '.csv'
            print(f'[INFO] Saving grouped time data frame to: {file_csv}')
            dfTime.to_csv(file_csv, sep=';', index_label=index_label)

    def plot_time_series(self, options_figure):
        if not self.VALID and not os.path.isfile(self.mrfile.file_path):
            print(f'[ERROR] {self.mrfile.file_path} shoud be a valid NetCDF file')
            return
        time_var = options_figure['time_var']
        avg_vars = options_figure['avg_var']
        if time_var is None or avg_vars is None:
            print(f'[ERROR] time_var and avg_var should be defined in the configuration file')
            return
        if time_var not in self.mrfile.nc.variables:
            print(f'[ERROR] {time_var} is not defined in {self.mrfile.file_path}')
            return
        # for avg_var in avg_vars:
        #     if avg_var not in self.mrfile.nc.variables:
        #         print(f'[ERROR] {avg_var} is not defined in {self.mrfile.file_path}')
        #         return
        # print(avg_vars)
        time_array = self.mrfile.get_full_array_1D(options_figure['time_var'], options_figure['insitu_id_variable'],
                                                   False)
        # print(time_array.shape)
        valid_array = np.ones(time_array.shape)
        valid_array[time_array.mask] = 0
        time_array = time_array[valid_array == 1]
        # print(time_array.shape)

        if options_figure['type_time_axis'] == 'fix':
            instant_values, time_array_instants = self.get_fix_time_axis(time_array, options_figure['fix_axis_options'])

        if options_figure['type_time_axis'] == 'variable':
            ##temporal, for tara, getting stations
            # stations = np.arange(1,99,1)
            # stations[49:99] = stations[49:99]+1
            # stations = stations[valid_array==1]

            ##gettting stations
            ###bar plot
            ##types: WFR-chl,'WFR-spm', 'MSI-chl', 'MSI-spm' ,'GLOBAL', 'REGIONAL'
            type = 'WFR-spm'
            min_value = 0
            max_value = 25
            from PlotSpectra import PlotSpectra
            pspectra = PlotSpectra()
            valid_var = None
            refs = ['CMEMS_MULTI_GLOBAL', 'CMEMS_OLCI_GLOBAL', 'CMEMS_OLCI_REGIONAL']
            for avg_var in avg_vars:
                if avg_var.find('insitu') > 0:
                    continue
                if avg_var in refs:
                    var_array = self.get_avg_var_ref(avg_var)
                else:
                    var_array = self.mrfile.nc.variables[avg_var][:]
                if valid_var is None:
                    valid_var = np.zeros(var_array.shape)
                valid_var[~var_array.mask] = 1

            if type == 'WFR-spm':
                valid_var[30] = 0
                valid_var[54] = 0
                valid_var[12] = 0

            n_valid = np.sum(valid_var)
            instant_values = np.arange(0, n_valid, 1).astype(np.int32)
            pspectra.xdata = instant_values
            colors = ['green', 'blue', 'red', 'magenta']
            offset = 0.1
            if type == 'WFR-chl' or type == 'WFR-spm' or type == 'REGIONAL':
                nbars_total = 3  ##3: WFR, REGIONAL; 2: MSI; 4: GLOBAL
            elif type == 'MSI-chl' or type == 'MSI-spm':
                nbars_total = 2
            elif type == 'GLOBAL':
                nbars_total = 4
            for idx, avg_var in enumerate(avg_vars):
                if avg_var in refs:
                    var_array = self.get_avg_var_ref(avg_var)
                else:
                    var_array = self.mrfile.nc.variables[avg_var][:]
                var_array = var_array[valid_var == 1]
                pspectra.plot_single_bar_series(var_array, colors[idx], 0.8 / nbars_total, offset, 0)
                offset = offset + (0.8 / nbars_total)
            if options_figure['title'] is not None:
                pspectra.set_title(options_figure['title'])

            valid_time_array = time_array[valid_var == 1]
            # valid_stations = stations[valid_var==1]
            from datetime import datetime as dt
            xticks_values = [dt.utcfromtimestamp(x).strftime('%d-%b') for x in valid_time_array]
            if type == 'GLOBAL' or type == 'REGIONAL':
                for i in range(0, len(xticks_values), 2):
                    xticks_values[i] = ''

            xticks_minor = instant_values + 0.35
            pspectra.set_xticks_minor(xticks_minor, xticks_values, 90, 10)

            xticks_major = np.arange(0, n_valid + 1, 1).astype(np.int32)
            if type == 'WFR-chl' or type == 'REGIONAL':
                xticks_major = xticks_major - 0.1  ##WFR chl-a,REGIONAL
            else:
                xticks_major = xticks_major - 0.2  ##MSI,GLOBAL,WFR smp
            pspectra.set_xticks(xticks_major, [], 0, 0)

            pspectra.set_grid_bars(0.1)
            pspectra.set_grid_horizontal()
            pspectra.set_xaxis_title('Date')
            if type.endswith('spm'):
                pspectra.set_yaxis_title(r'SPM (g m$^-$$^3$)')
            else:
                pspectra.set_yaxis_title(r'chl (mg m$^-$$^3$)')

            if min_value is not None and max_value is not None:
                pspectra.set_y_range(min_value, max_value)

            ##WFR
            if type == 'WFR-chl':
                legend_str = ['in situ chl (HPLC)', 'satellite chl (NN)', 'satellite chl (OC4ME)']
                pspectra.legend_options['loc'] = 'lower center'
                pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.35)
                pspectra.legend_options['ncols'] = 3

            ##WFR-SPM
            if type == 'WFR-spm':
                legend_str = ['in situ SPM', 'satellite SPM (NN)']
                pspectra.legend_options['loc'] = 'lower center'
                pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.35)
                pspectra.legend_options['ncols'] = 2
                pspectra.set_yticks([0, 2, 5, 6, 8, 12], None, 0, None)

            ##MSI-chla
            if type == 'MSI-chl':
                legend_str = ['in situ chl (HPLC)', 'satellite chl']
                pspectra.legend_options['loc'] = 'lower center'
                pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.35)
                pspectra.legend_options['ncols'] = 2
                pspectra.set_yticks([0, 5, 10, 15, 20], None, 0, None)

            ##MSI-spm
            if type == 'MSI-spm':
                legend_str = ['in situ SPM', 'satellite SPM']
                pspectra.legend_options['loc'] = 'lower center'
                pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.35)
                pspectra.legend_options['ncols'] = 2
                pspectra.set_yticks([0, 5, 10, 15, 20, 25], None, 0, None)

            ##GLOBAL
            if type == 'GLOBAL':
                legend_str = ['in situ chl (HPLC)', 'OC-CCI chl', 'CMEMS-MULTI chl', 'CMEMS-OLCI chl']
                pspectra.legend_options['loc'] = 'lower center'
                pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.49)
                pspectra.legend_options['ncols'] = 2
                pspectra.set_yticks([0, 5, 10, 15, 20], None, 0, None)

            ##REGIONAL
            if type == 'REGIONAL':
                legend_str = ['in situ chl (HPLC)', 'CMEMS-MULTI chl', 'CMEMS-OLCI chl']
                pspectra.legend_options['loc'] = 'lower center'
                pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.46)
                pspectra.legend_options['ncols'] = 3
                pspectra.set_yticks([0, 5, 10, 15, 20], None, 0, None)

            pspectra.set_legend(legend_str)
            pspectra.set_tigth_layout()
            file_out = options_figure['file_out']
            if file_out is not None:
                pspectra.save_plot(file_out)

            return

        from PlotSpectra import PlotSpectra
        import MDBPlotDefaults as defaults
        pspectra = PlotSpectra()
        pspectra.xdata = instant_values

        style = pspectra.line_style_default.copy()
        style['linewidth'] = 0
        style['marker'] = 'o'
        style['markersize'] = 5
        colors = ['blue', 'red', 'green']

        ##temporal

        insitu_array = np.ma.masked_all((225,))
        global_array = np.ma.masked_all((225,))
        regional_array = np.ma.masked_all((225,))

        handles = []
        for index, avg_var in enumerate(avg_vars):
            # style['color'] = defaults.colors_default[index]
            style['color'] = colors[index]
            print(f'[INFO] [PLOT] Plotting variable: {avg_var}')
            # var_array =  self.mrfile.get_full_array_1D(avg_var,options_figure['insitu_id_variable'], False)
            var_array = self.mrfile.nc.variables[avg_var][:]

            if options_figure['type_time_axis'] == 'fix':
                if options_figure['method_fix_axis'] == 'all':
                    for ivalue in instant_values:
                        ydata = var_array[time_array_instants == ivalue]
                        # ydata = ydata[~ydata.mask]
                        # print(ivalue,len(ydata),ydata)
                        if len(ydata) == 0:
                            continue
                        # print(ydata,len(ydata))
                        ydata = ydata[0]

                        if index == 0:
                            # print(ivalue, '--->', ydata)
                            insitu_array[instant_values == ivalue] = ydata
                        if index == 1:
                            # print(ivalue, '--->', ydata)
                            global_array[instant_values == ivalue] = ydata
                        if index == 2:
                            regional_array[instant_values == ivalue] = ydata
                        # ydata = var_array[time_array_instants==ivalue]
                        # ydata = ydata[~ydata.mask]
                        # if len(ydata)==0:
                        #     continue
                        # xdata = np.zeros(len(ydata))
                        # xdata[:] = ivalue
                        # pspectra.plot_single_data(xdata,ydata,style)

        # print('??????????????????????????????????')
        # print(instant_values.shape)
        # print(pspectra.xdata.shape)
        # print(insitu_array.shape)
        # print(global_array.shape)
        # print(regional_array.shape)
        # print(style)
        style['linewidth'] = 1
        style['markersize'] = 0
        style['color'] = 'gray'
        pspectra.plot_data(insitu_array, style)
        style['color'] = 'blue'
        style['markersize'] = 5
        pspectra.plot_data(global_array, style)
        style['color'] = 'green'
        style['markersize'] = 5
        pspectra.plot_data(regional_array, style)
        pspectra.set_grid()

        #     var_array = np.array(self.mrfile.nc.variables[avg_var]).astype(np.float)
        #     var_array[var_array < -1] = np.nan
        #     if index == 1 and not options_figure['log_scale']:
        #         var_array[~np.isnan(var_array)] = var_array[~np.isnan(var_array)] / np.pi
        #     if options_figure['log_scale']:
        #         var_array[~np.isnan(var_array)] = np.log10(var_array[~np.isnan(var_array)])
        #     h = pspectra.plot_data(var_array, style)

        # time = np.array(self.mrfile.nc.variables[time_var])
        # from datetime import datetime as dt
        # time_ini_year = [
        #     dt.utcfromtimestamp(float(x)).replace(day=1, month=1, hour=0, minute=0, second=0, microsecond=0,
        #                                           tzinfo=pytz.UTC).timestamp() for x in time]
        # time_fin_year = [dt.utcfromtimestamp(float(x)).replace(tzinfo=pytz.UTC,
        #                                                        year=dt.utcfromtimestamp(float(x)).year + 1).timestamp()
        #                  for x in time_ini_year]
        # seconds_year = np.array(time_fin_year) - np.array(time_ini_year)
        # time_array = []
        # for there, ini_year, total_year in zip(time, time_ini_year, seconds_year):
        #     val = dt.utcfromtimestamp(there).year + ((there - ini_year) / total_year)
        #     time_array.append(val)
        # time_array = np.array(time_array)
        # width = (time_array[-1] - time_array[0]) / len(time_array)
        #
        # dispersion_min_var = options_figure['dispersion_min_var']
        # dispersion_max_var = options_figure['dispersion_max_var']

        # pspectra.xdata = time_array
        # style = pspectra.line_style_default.copy()
        # style['linewidth'] = 0
        # style['marker'] = 'o'
        # style['markersize'] = 1
        # handles = []
        # for index, avg_var in enumerate(avg_vars):
        #     style['color'] = defaults.colors_default[index]
        #     var_array = np.array(self.mrfile.nc.variables[avg_var]).astype(np.float)
        #     var_array[var_array < -1] = np.nan
        #     if index == 1 and not options_figure['log_scale']:
        #         var_array[~np.isnan(var_array)] = var_array[~np.isnan(var_array)] / np.pi
        #     if options_figure['log_scale']:
        #         var_array[~np.isnan(var_array)] = np.log10(var_array[~np.isnan(var_array)])
        #     h = pspectra.plot_data(var_array, style)
        #
        #     ##temporal, owt
        #     # h = pspectra.plot_single_bar_series(var_array,style['color'],width,0,0)
        #
        #     handles.append(h[0])
        #     if dispersion_min_var is not None and dispersion_max_var is not None:
        #         if len(dispersion_min_var) == len(avg_vars) and len(dispersion_max_var) == len(avg_vars):
        #             min_dispersion_array = np.array(self.mrfile.nc.variables[dispersion_min_var[index]])
        #             max_dispersion_array = np.array(self.mrfile.nc.variables[dispersion_max_var[index]])
        #             min_dispersion_array[min_dispersion_array < -1.0] = np.nan
        #             max_dispersion_array[max_dispersion_array < -1.0] = np.nan
        #             if index == 1 and not options_figure['log_scale']:
        #                 min_dispersion_array[~np.isnan(min_dispersion_array)] = min_dispersion_array[~np.isnan(
        #                     min_dispersion_array)] / np.pi
        #                 max_dispersion_array[~np.isnan(max_dispersion_array)] = max_dispersion_array[~np.isnan(
        #                     max_dispersion_array)] / np.pi
        #             if options_figure['log_scale']:
        #                 min_dispersion_array[~np.isnan(min_dispersion_array)] = np.log10(
        #                     min_dispersion_array[~np.isnan(min_dispersion_array)])
        #                 max_dispersion_array[~np.isnan(max_dispersion_array)] = np.log10(
        #                     max_dispersion_array[~np.isnan(max_dispersion_array)])
        #             pspectra.plot_iqr_basic(min_dispersion_array, max_dispersion_array, style['color'])
        #
        # ##temporal, for owt
        # # yticks = np.arange(1,19)
        # # pspectra.set_yticks(yticks,yticks,0,10)
        #
        # ##temporal
        # ymin = options_figure['y_min']
        # ymax = options_figure['y_max']
        # if ymin is not None or ymax is not None:
        #     pspectra.set_y_range(ymin, ymax)
        # if options_figure['ylabel'] is not None:
        #     pspectra.set_yaxis_title(options_figure['ylabel'])
        # if options_figure['xlabel'] is not None:
        #     pspectra.set_xaxis_title(options_figure['xlabel'])
        # pspectra.set_grid()
        # if options_figure['legend'] and options_figure['legend_values'] is not None:
        #     pspectra.legend_options['loc'] = 'lower center'
        #     pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.25)
        #     pspectra.legend_options['ncols'] = 2
        #     pspectra.legend_options['markerscale'] = 5
        #
        #     pspectra.set_legend_h(handles, options_figure['legend_values'])

        if options_figure['title'] is not None:
            pspectra.set_title(options_figure['title'])

        pspectra.set_tigth_layout()
        file_out = options_figure['file_out']
        if file_out is not None:
            pspectra.save_plot(file_out)

    def get_avg_var_ref(self, ref):
        dir_base = '/mnt/c/DATA_LUIS/TARA_TEST/station_match-ups/MDBs/MDBr_chla'
        from netCDF4 import Dataset
        if ref == 'CMEMS_MULTI_GLOBAL':
            file = os.path.join(dir_base, 'MDBr__CMEMS_MULTI_4KM_CMEMS-MULTI_20240101T000000_20240919T235959.nc')
            dataset = Dataset(file)
            var_array = dataset.variables['mu_satellite_CHL'][:]
            dataset.close()
            return var_array
        if ref == 'CMEMS_OLCI_GLOBAL':
            file = os.path.join(dir_base, 'MDBr__CMEMS_OLCI_300M_CMEMS-OLCI_20240101T000000_20240919T235959.nc')
            dataset = Dataset(file)
            var_array = dataset.variables['mu_satellite_CHL'][:]
            dataset.close()
            return var_array
        if ref == 'CMEMS_OLCI_REGIONAL':
            file = os.path.join(dir_base,
                                'MDBr__CMEMS_OLCI_300M_CMEMS-OLCI-REGIONAL_20240101T000000_20240919T235959.nc')
            dataset = Dataset(file)
            var_array = dataset.variables['mu_satellite_CHL'][:]
            dataset.close()
            return var_array

    def get_fix_time_axis(self, time_array, options):
        #print('aqui')
        #print(options)
        from datetime import datetime as dt
        time_array_instants = time_array.copy().astype(np.int32)
        format_abs = options['format_abs']
        min_time_abs = options['min_abs']
        max_time_abs = options['max_abs']
        check_min_time = False if min_time_abs is None else True
        check_max_time = False if max_time_abs is None else True

        for index, time in enumerate(time_array):
            time_here = dt.utcfromtimestamp(time)
            instant_here = int(time_here.strftime(format_abs))
            time_array_instants[index] = instant_here
            if check_min_time and instant_here < min_time_abs: time_array_instants[index] = -1
            if check_max_time and instant_here > max_time_abs: time_array_instants[index] = -1
            if not check_min_time and (
                    (min_time_abs is None) or (min_time_abs is not None and instant_here < min_time_abs)):
                min_time_abs = instant_here
            if not check_max_time and (
                    (max_time_abs is None) or (max_time_abs is not None and instant_here > max_time_abs)):
                max_time_abs = instant_here

        # print(min_time_abs,max_time_abs)
        instant_values = np.arange(min_time_abs, max_time_abs + 1, 1).astype(np.int32)
        # print(instant_values )

        return instant_values, time_array_instants

    def plot_scatterplot_from_options(self, options_figure):
        ##WORKING BY MATCH-UPS
        if options_figure['selectByMu']:
            index_mu = options_figure['index_mu']
            mu_valid = np.ones((self.mrfile.n_mu_total,))
            if self.mu_valid_variable in self.mrfile.variables:
                mu_valid = self.mrfile.variables[self.mu_valid_variable][:]
            file_out_base = options_figure['file_out']
            if index_mu == -1:
                for imu in range(self.mrfile.n_mu_total):
                    if mu_valid[imu] == 1:
                        print(f'[INFO] Plotting scatter plot for match-up {imu} / {self.mrfile.n_mu_total}')
                        self.plot_scatter_plot_mu(options_figure, imu, file_out_base)
            elif 0 <= index_mu < self.mrfile.n_mu_total and mu_valid[index_mu] == 1:
                print(f'[INFO] Plotting scatter plot for match-up {index_mu}')
                self.plot_scatter_plot_mu(options_figure, index_mu, file_out_base)

            return

        ##WORKING WITH ALL THE DATA
        if options_figure['selectBy'] is None or not options_figure['individual_plots']:
            if options_figure['type_scatterplot'] == 'rrs':
                if not options_figure['selectByWavelength']:  # GLOBAL SCATTERPLOT
                    self.plot_global_scatterplot(options_figure)
                else:  # AN SCATTERPLOT BY WAVELENGTH
                    if options_figure['multiple_plot'] is not None:  # single file
                        self.plot_multiple_wavelength_scatterplots_single_file(options_figure)
                    else:  # multiple files
                        self.plot_multiple_wavelength_scatterplots_multiple_files(options_figure)
            else:
                print(f'[INFO] Type scatterplot: {options_figure["type_scatterplot"]}')
                self.plot_general_scatterplot(options_figure)

        ##WORKING WITH SELECTED OPTIONS

        if options_figure['selectBy'] is not None and options_figure['individual_plots']:

            selectValues = options_figure['selectValues']

            file_out_base = options_figure['file_out']
            title_base = options_figure['title']
            for svalue in selectValues:
                #print('-->',svalue)
                options_figure['selectValues'] = svalue
                flag = self.get_str_select_value(options_figure, svalue)
                options_figure['file_out'] = self.get_file_out_name(file_out_base, None, flag)
                options_figure['title'] = self.get_title(title_base, None, flag, None)
                if options_figure['type_scatterplot'] != 'rrs':  ##GENERAL SCATTER PLOT
                    self.plot_general_scatterplot(options_figure)
                else:
                    if not options_figure['selectByWavelength']:  # GLOBAL SCATTERPLOT
                        self.plot_global_scatterplot(options_figure)
                    else:  # AN SCATTERPLOT BY WAVELENGTH
                        if options_figure['multiple_plot'] is not None:  # single file
                            self.plot_multiple_wavelength_scatterplots_single_file(options_figure)
                        else:  # multiple files
                            self.plot_multiple_wavelength_scatterplots_multiple_files(options_figure)

    def plot_general_scatterplot(self, options_figure):
        self.set_data_scatterplot_general(options_figure['groupBy'], options_figure['selectBy'],
                                          options_figure['selectValues'], options_figure)
        self.plot_scatter_plot(options_figure, None, -1, -1, -1)

    def plot_scatter_plot_mu(self, options_figure, index_mu, file_out_base):
        if options_figure['apply_wavelength_color'] and options_figure['groupBy'] is None:
            if options_figure['wlranges_min'] is not None and options_figure['wlranges_max'] is not None:
                options_figure = self.create_flag_array_wl_ranges(options_figure, 'wl_groups')
                options_figure['groupBy'] = 'wl_groups'
                options_figure['groupValues'] = self.virtual_flags['wl_groups']['flag_values']
                options_figure['groupType'] = 'flag'
            else:
                options_figure['groupBy'] = 'mu_wavelength'
                options_figure['groupValues'] = options_figure['wlvalues']
                options_figure['groupType'] = 'wavelength'

        if file_out_base is not None:
            satellite_time = self.mrfile.sat_times[index_mu].strftime('%Y%m%d')
            file_out = f'{file_out_base[:-4]}_{satellite_time}_{index_mu}{file_out_base[-4:]}'
            options_figure['file_out'] = file_out

        self.set_data_scatterplot(options_figure['groupBy'], 'mu_satellite_id', index_mu, None, options_figure)
        self.plot_scatter_plot(options_figure, None, -1, -1, -1)

    def plot_global_scatterplot(self, options_figure):
        if options_figure['apply_wavelength_color'] and options_figure['groupBy'] is None:
            if options_figure['wlranges_min'] is not None and options_figure['wlranges_max'] is not None:
                options_figure = self.create_flag_array_wl_ranges(options_figure, 'wl_groups')
                options_figure['groupBy'] = 'wl_groups'
                options_figure['groupValues'] = self.virtual_flags['wl_groups']['flag_values']
                options_figure['groupType'] = 'flag'
            else:
                if options_figure['wlvalues'] is None:
                    options_figure['wlvalues'] = np.unique(self.mrfile.nc.variables['mu_wavelength'][:])
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
        if options_figure['wlvalues'] is None:
            options_figure['wlvalues'] = np.unique(self.mrfile.nc.variables['mu_wavelength'][:])
        wl_values = options_figure['wlvalues']
        print(f'[INFO] Wavelenthts: {wl_values}')
        nblank = ntot - len(wl_values)
        if nblank >= nrow:
            print(
                f'[ERROR] Number of total axis {plot_here.nrow}x{plot_here.ncol} is higher than number of scatterplots:{len(wl_values)}')
            plot_here.close_plot()
            return
        print(f'[INFO] Number of total axis {plot_here.nrow}x{plot_here.ncol} Blank plots: {nblank}')

        index_col_adjust = plot_here.ncol - nblank if nblank > 0 else -1
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
        if options_figure['wlvalues'] is None:
            options_figure['wlvalues'] = np.unique(self.mrfile.nc.variables['mu_wavelength'][:])
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
            #print(self.valid_stats)
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
                print(f'[INFO] Number of data points for group {g}: {len(xhere)}')

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
                # if g==1:
                #     print(xhere.shape,yhere.shape,np.mean(xhere),np.mean(yhere))
                #     print(marker,markersize,color,edgecolor,linewidth)
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
                    # xhere = xhere[0:100000]
                    # yhere = yhere[0:100000]
                    xy = np.vstack([xhere, yhere])

                try:
                    print(f'[INFO] Computing density...')

                    z = gaussian_kde(xy)(xy)

                    # z = self.mrfile.variables['DISTANCE'][:]

                    print(f'[INFO]Sorting density...')
                    idx = z.argsort()
                    xhere, yhere, z = xhere[idx], yhere[idx], z[idx]
                    print(f'[INFO] Density values were sorted.')
                    plot.set_cmap('jet')


                    hscatter = plot.plot_data(xhere, yhere, marker, markersize, z, None, 0)

                    # file_kk = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/PLOTS/stal.csv'
                    # fw = open(file_kk,'w')
                    # fw.write('x;y')
                    # for idx in range(len(xhere)):
                    #     fw.write('\n')
                    #     fw.write(f'{xhere[idx]};{yhere[idx]}')
                    # fw.close()

                except:
                    print(f'[ERROR] Error creating density plot. Using default style')
                    plot.plot_data(xhere, yhere, marker, markersize, color, edgecolor, linewidth)

            else:
                plot.plot_data(xhere, yhere, marker, markersize, color, edgecolor, linewidth)


        # plot.colorbar(hscatter)

        ##limitss
        if options['log_scale']:
            plot.set_log_scale()
            if options['min_xy'] is None:
                options['min_xy'] = 0.1
                if options['type_scatterplot'] == 'kd':
                    options['min_xy'] = 0.01
            min_xy = options['min_xy']
            if options['max_xy'] is None:
                options['max_xy'] = 100
                if options['type_scatterplot'] == 'kd':
                    options['max_xy'] = 10
            max_xy = options['max_xy']
        else:
            if options['min_xy'] is None or options['max_xy'] is None:
                min_xy, max_xy = self.get_min_max_xy()
        if options['min_xy'] is not None:
            min_xy = options['min_xy']
        if options['max_xy'] is not None:
            max_xy = options['max_xy']

        # TEST FOR A PLOT
        # if wl==665:
        #     #max_xy = 0.008
        #     max_xy = 0.010

        min_x = min_xy
        max_x = max_xy
        min_y = min_xy
        max_y = max_xy
        ##differente limits for y and x
        if options['x_min'] is not None:
            min_x = options['x_min']
        if options['x_max'] is not None:
            max_x = options['x_max']
        if options['y_min'] is not None:
            min_y = options['y_min']
        if options['y_max'] is not None:
            max_y = options['y_max']

        plot.set_limits_X(min_x, max_x)
        plot.set_limits_Y(min_y, max_y)
        # plot.set_limits(min_xy, max_xy)

        # ticks
        if options['ticks'] is None and (options['x_ticks'] is None or options['y_ticks'] is None):
            if not options['log_scale']:
                # ticks = self.get_ticks_from_min_max_xy(min_xy, max_xy)
                x_ticks = self.get_ticks_from_min_max_xy(min_x, max_x)
                y_ticks = self.get_ticks_from_min_max_xy(min_y, max_y)

            else:
                min_tx = int(np.log10(min_x))
                max_tx = int(np.log10(max_x))
                x_ticks = [10 ** x for x in range(min_tx, max_tx + 1)]
                min_ty = int(np.log10(min_y))
                max_ty = int(np.log10(max_y))
                y_ticks = [10 ** x for x in range(min_ty, max_ty + 1)]
        else:
            if options['ticks'] is not None:
                x_ticks = options['ticks']
                y_ticks = options['ticks']
            elif options['x_ticks'] is not None and options['y_ticks'] is not None:
                x_ticks = options['x_ticks']
                y_ticks = options['y_ticks']

        ##TEST FOR A PLOT
        # if wl==665:
        #     # x_ticks = [0,0.002,0.004,0.006,0.008]
        #     # y_ticks = [0, 0.002, 0.004, 0.006, 0.008]
        #     x_ticks = [0, 0.002, 0.004, 0.006, 0.008,0.010]
        #     y_ticks = [0, 0.002, 0.004, 0.006, 0.008, 0.010]

        if x_ticks is not None and y_ticks is not None:
            if options['log_scale']:
                txlabels = self.get_labels_for_log_ticks(x_ticks)
                tylabels = self.get_labels_for_log_ticks(y_ticks)
                plot.set_ticks_and_labels_x(x_ticks, txlabels, options['fontsizeaxis'])
                plot.set_ticks_and_labels_y(y_ticks, tylabels, options['fontsizeaxis'])
            else:
                plot.set_ticks_x(x_ticks, options['fontsizeaxis'])
                plot.set_ticks_y(y_ticks, options['fontsizeaxis'])


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
                plot.set_yticks_labels_off(y_ticks)
            if plot.index_row == prefinal_row and plot.index_col >= index_col_adjust >= 1:
                plot.set_xaxis_title(options['xlabel'])
            else:
                if plot.index_row < final_row:
                    plot.set_xticks_labels_off(x_ticks)

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
            plot.plot_text_options['fontsize'] = options['fontsizestats']
            plot.plot_text(xpos, ypos, str0)

        # regression lines
        if options['log_scale']:
            if options['regression_line']:
                xr = np.array(self.xregress)
                yr = np.array(self.yregress)
                xr = xr[xr.argsort()]
                yr = yr[yr.argsort()]
                plot.plot_regress_line(xr, yr, 'black')
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
            plot.ax.title.set_size(options['fontsizetitle'])

        #print('PREPARE POSTER OPTION')
        #plot.prepare_poster()

        ##annotations
        if not options['anot_'] is None:
            anots = options['anot_']
            style = options['anot_default_style']
            fontsize = int(style[0])
            for anot in anots:
                anot_info = anots[anot]
                ypos = float(anot_info[0])
                xpos = float(anot_info[1])
                text = anot_info[2]
                print(xpos,ypos,text,fontsize)
                plot.set_text_size(xpos,ypos,text,fontsize)




        ##saving to file
        if not options['file_out'] is None:
            plot.save_fig(options['file_out'])
            plot.close_plot()

        return plot

    def plot_spectraplot_from_options(self, options_figure):

        if options_figure['xlabel'] is None:
            options_figure['xlabel'] = 'Wavelength (nm)'

        if options_figure['type_rrs'] == 'ins':
            self.plot_insitu_spectraplots(options_figure)

        if options_figure['type_rrs'] == 'mu_ins':
            index_mu = options_figure['index_mu']
            if index_mu == -1:
                for imu in range(self.mrfile.n_mu_total):
                    print(f'[INFO] Plotting ins spectra for match-up {imu} / {self.mrfile.n_mu_total}')
                    self.plot_insmu_spectraplot(options_figure, imu)
            elif 0 <= index_mu < self.mrfile.n_mu_total:
                print(f'[INFO] Plotting ins spectra for match-up {index_mu}')
                self.plot_insmu_spectraplot(options_figure, index_mu)

        if options_figure['type_rrs'] == 'mu_comparison' or options_figure['type_rrs'] == 'mu_sat':
            mu_valid = np.ones((self.mrfile.n_mu_total,))
            if self.mu_valid_variable in self.mrfile.variables:
                mu_valid = self.mrfile.variables[self.mu_valid_variable][:]
            index_mu = options_figure['index_mu']
            if index_mu == -1:
                for imu in range(self.mrfile.n_mu_total):
                    if mu_valid[imu] == 1:
                        self.plot_mu_spectraplot(options_figure, imu)
            elif 0 <= index_mu < self.mrfile.n_mu_total and mu_valid[index_mu] == 1:
                self.plot_mu_spectraplot(options_figure, index_mu)

        if options_figure['type_rrs'] == 'comparison_sat_insitu':
            print(f'[INFO] Plotting comparison sat-insitu')
            if options_figure['selectBy'] is None:
                if options_figure['groupBy'] is None:
                    self.plot_comparison_sat_insitu(options_figure,None,None,None)
                else:
                    self.plot_comparison_sat_insitu_grouped(options_figure)
            else:
                #var_select = options_figure['selectBy']
                select_values = options_figure['selectValues']
                file_out_base = options_figure['file_out']
                for svalue in select_values:
                    #options_figure['selectValues'] = svalue
                    flag = self.get_str_select_value(options_figure, svalue)
                    options_figure['file_out'] = self.get_file_out_name(file_out_base, None, flag)
                    #options_figure['title'] = self.get_title(title_base, None, flag, None)
                    self.plot_comparison_sat_insitu(options_figure, options_figure['selectBy'], svalue, None)

    def plot_comparison_sat_insitu_grouped(self, options_figure):
        wlvalues = options_figure['wlvalues']
        if wlvalues is None:
            wlvalues = list(np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])))
        self.mrfile.var_mu_valid = self.mu_valid_variable
        from PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        if options_figure['stat_plot_method'] is not None:
            if options_figure['stat_plot_method'] == 'iqr':
                pspectra.set_iqr_as_stats_plot()
            if options_figure['stat_plot_method'].startswith('std'):
                try:
                    factor = float(options_figure['stat_plot_method'].split(',')[1])
                except:
                    factor = 1.0
                pspectra.set_std_as_stats_plot(factor)
        imin, imax = pspectra.get_imin_imax_from_wavelength(np.array(wlvalues), options_figure['wl_min'],options_figure['wl_max'])

        flag_name = options_figure['groupBy']
        group_values = options_figure['groupValues']
        str_legend = ['insitu Rrs']
        sat_stats_array = []


        for gvalue in group_values:
            print(f'[INFO] Computing stats for group: {gvalue}')
            gvalue_flag = self.get_flag_flag(gvalue,options_figure[flag_name]['flag_values'],options_figure[flag_name]['flag_meanings'])
            str_legend.append(gvalue_flag)
            sat_stats, insitu_stats = self.mrfile.get_all_spectra_insitu_sat_with_wlvalues(options_figure['scale_factor'],wlvalues, options_figure['groupBy'], gvalue, None)
            sat_stats_array.append(sat_stats)

        if options_figure['y_min'] is None or options_figure['y_max'] is None:

            ymin_insitu, ymax_insitu = pspectra.get_ymin_ymax_from_stats(insitu_stats, imin, imax)
            ymin_list = [ymin_insitu]
            ymax_list = [ymax_insitu]
            for sat_stats in sat_stats_array:
                ymin_sat, ymax_sat = pspectra.get_ymin_ymax_from_stats(sat_stats, imin, imax)
                ymin_list.append(ymin_sat)
                ymax_list.append(ymax_sat)
            ymin = np.min(np.array(ymin_list))
            ymax = np.max(np.array(ymax_list))
        if options_figure['y_min'] is not None:
            ymin = options_figure['y_min']
        if options_figure['y_max'] is not None:
            ymax = options_figure['y_max']

        wl_list = wlvalues[imin:imax]
        xdata_plot = [float(self.get_wl_str_from_wl(x)) for x in wl_list]
        wl_col = [self.get_wl_str_from_wl(x) for x in wl_list]

        #CHANGE LYNE STYLE LINE FOR POSTER.
        #pspectra.stats_style['central']['linewidth']=2

        pspectra.xdata = xdata_plot
        h_legend = []
        ##insitu
        color = 'gray'
        pspectra.stats_style['central']['color'] = color
        pspectra.stats_style['central']['marker'] = 'o'
        pspectra.stats_style['central']['markersize'] = 5
        if len(wl_col) > 100:
            pspectra.stats_style['central']['markersize'] = 0.5
        pspectra.stats_style['fill']['color'] = color
        pspectra.stats_style['fill']['framealpha'] = 0.5
        hlineinsitu = pspectra.plot_stats(insitu_stats, imin, imax)
        h_legend.append(hlineinsitu[0])

        colors = ['blue','red']
        if len(sat_stats_array)==2:
            colors = ['blue','red']

        for igroup,sat_stats in enumerate(sat_stats_array):
            pspectra.stats_style['central']['color'] = colors[igroup]
            # pspectra.stats_style['dispersion']['color'] = colors[igroup]
            # pspectra.stats_style['dispersion']['linewidth'] = 1
            # pspectra.stats_style['dispersion']['linestyle'] = '--'
            pspectra.stats_style['fill']['color'] = colors[igroup]
            hlinesat = pspectra.plot_stats(sat_stats, imin, imax)
            h_legend.append(hlinesat[0])

        xticks_size = 12
        if len(xdata_plot) == 16 or len(xdata_plot) == 15:
            xticks_size = 10

        if options_figure['x_ticks'] is None:
            if len(wl_col) > 100:
                xdata_plot = np.array([350, 400, 450, 500, 550, 600, 650, 700])
                wl_col = [f'{x}' for x in xdata_plot]
                pspectra.set_xticks(xdata_plot, wl_col, 0, xticks_size)
            else:
                pspectra.set_xticks(xdata_plot, wl_col, 90, xticks_size)
        else:
            xdata_plot = options_figure['x_ticks']
            wl_col = [f'{x:.2f}'.replace('.00','') for x in xdata_plot]
            pspectra.set_xticks(xdata_plot, wl_col, 0, xticks_size)

        if len(xdata_plot) == 16 or len(xdata_plot) == 15:
            xticks, xlabels = pspectra.get_xticks()
            xlabels[8].set_visible(False)

        pspectra.set_yticks(None, None, None, 12)

        pspectra.set_y_range(ymin, ymax)

        pspectra.set_xaxis_title(options_figure['xlabel'])
        ylabel = options_figure['ylabel'] if options_figure['ylabel'] is not None else defaults.ylabel_rrs_scaled
        pspectra.set_yaxis_title(ylabel)
        if options_figure['title'] is not None:
            title_here = options_figure['title']
            pspectra.set_title(title_here)
        pspectra.set_grid()
        pspectra.legend_options['bbox_to_anchor'] = (0.65, 1.0)
        pspectra.legend_options['framealpha'] = 1
        if options_figure['legend']: ##FOR POSTER, ADD FONTSIZE=14 IN set_legend_h
            pspectra.set_legend_h(h_legend, str_legend)
        pspectra.set_tigth_layout()
        file_out = options_figure['file_out']

        ##POSTER STYLE
        #pspectra.prepare_poster()

        if not file_out is None:
            pspectra.save_fig(file_out)
        pspectra.close_plot()

    def plot_comparison_sat_insitu(self, options_figure,flag_name,flag_value,flag_array):
        wlvalues = options_figure['wlvalues']
        if wlvalues is None:
            wlvalues = list(np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])))
        self.mrfile.var_mu_valid = self.mu_valid_variable
        from PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        if options_figure['stat_plot_method'] is not None:
            if options_figure['stat_plot_method'] == 'iqr':
                pspectra.set_iqr_as_stats_plot()
            if options_figure['stat_plot_method'].startswith('std'):
                try:
                    factor = float(options_figure['stat_plot_method'].split(',')[1])
                except:
                    factor = 1.0
                pspectra.set_std_as_stats_plot(factor)


        sat_stats, insitu_stats = self.mrfile.get_all_spectra_insitu_sat_with_wlvalues(options_figure['scale_factor'],
                                                                                       wlvalues, flag_name, flag_value, flag_array)

        imin, imax = pspectra.get_imin_imax_from_wavelength(np.array(wlvalues), options_figure['wl_min'],
                                                            options_figure['wl_max'])
        str_legend = ['insitu Rrs', 'satellite Rrs']

        if options_figure['y_min'] is None or options_figure['y_max'] is None:
            ymin_sat, ymax_sat = pspectra.get_ymin_ymax_from_stats(sat_stats, imin, imax)
            ymin_insitu, ymax_insitu = pspectra.get_ymin_ymax_from_stats(insitu_stats, imin, imax)
            ymin = np.min([ymin_sat, ymin_insitu])
            ymax = np.max([ymax_sat, ymax_insitu])
        if options_figure['y_min'] is not None:
            ymin = options_figure['y_min']
        if options_figure['y_max'] is not None:
            ymax = options_figure['y_max']

        wl_list = wlvalues[imin:imax]
        xdata_plot = [float(self.get_wl_str_from_wl(x)) for x in wl_list]
        wl_col = [self.get_wl_str_from_wl(x) for x in wl_list]

        pspectra.xdata = xdata_plot

        color = 'red'
        pspectra.stats_style['central']['color'] = color
        pspectra.stats_style['central']['marker'] = 'o'
        pspectra.stats_style['central']['markersize'] = 5
        if len(wl_col) > 100:
            pspectra.stats_style['central']['markersize'] = 0.5
        pspectra.stats_style['fill']['color'] = color
        pspectra.stats_style['fill']['framealpha'] = 0.5
        hlineinsitu = pspectra.plot_stats(insitu_stats, imin, imax)

        color = 'blue'
        pspectra.xdata = wlvalues[imin:imax]
        pspectra.stats_style['central']['color'] = color
        pspectra.stats_style['central']['marker'] = 'o'
        pspectra.stats_style['central']['markersize'] = 5
        if len(wl_col) > 100:
            pspectra.stats_style['central']['markersize'] = 0.5
        pspectra.stats_style['fill']['color'] = color
        pspectra.stats_style['fill']['framealpha'] = 0.5
        hlinesat = pspectra.plot_stats(sat_stats, imin, imax)

        h_legend = [hlineinsitu[0], hlinesat[0]]

        xticks_size = 12

        if len(xdata_plot) == 16 or len(xdata_plot) == 15:
            xticks_size = 10
        if options_figure['x_ticks'] is None:
            if len(wl_col) > 100:
                xdata_plot = np.array([350, 400, 450, 500, 550, 600, 650, 700])
                wl_col = [f'{x}' for x in xdata_plot]
                pspectra.set_xticks(xdata_plot, wl_col, 0, xticks_size)
            else:
                pspectra.set_xticks(xdata_plot, wl_col, 90, xticks_size)
        else:
            xdata_plot = options_figure['x_ticks']
            wl_col = [f'{x}' for x in xdata_plot]
            pspectra.set_xticks(xdata_plot, wl_col, 0, xticks_size)

        if len(xdata_plot) == 16 or len(xdata_plot) == 15:
            xticks, xlabels = pspectra.get_xticks()
            xlabels[8].set_visible(False)
        pspectra.set_yticks(None, None, None, 12)

        pspectra.set_y_range(ymin, ymax)

        pspectra.set_xaxis_title(options_figure['xlabel'])
        ylabel = options_figure['ylabel'] if options_figure['ylabel'] is not None else defaults.ylabel_rrs_scaled
        pspectra.set_yaxis_title(ylabel)
        if options_figure['title'] is not None:
            title_here = options_figure['title']
            pspectra.set_title(title_here)
        pspectra.set_grid()
        pspectra.legend_options['bbox_to_anchor'] = (0.65, 1.0)
        pspectra.legend_options['framealpha'] = 1
        if options_figure['legend']:
            pspectra.set_legend_h(h_legend, str_legend)
        pspectra.set_tigth_layout()
        file_out = options_figure['file_out']
        if not file_out is None:
            pspectra.save_fig(file_out)
        pspectra.close_plot()

    def plot_mu_spectraplot(self, options_figure, index_mu):
        wl, insitu_spectra, sat_spectra, insitu_spectra_unc, sat_spectra_unc = self.mrfile.get_mu_spectra_insitu_and_sat(
            index_mu, options_figure['scale_factor'])

        if wl is None:
            return
        from PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        pspectra.xdata = wl
        if options_figure['wlticks'] is not None:
            wlt = options_figure['wlticks']
            wls = [f'{x:.2f}'.replace('.00','') for x in wlt]
        else:
            wlt = wl
            wls = self.mrfile.get_sat_wl_as_strlist(wl)

        pspectra.set_xticks(wlt, wls, 0, 12)

        if options_figure['type_rrs'] == 'mu_comparison':
            hline1 = pspectra.plot_single_line(insitu_spectra, 'red', 'solid', 2, '.', 0)
            hline2 = pspectra.plot_single_line(sat_spectra, 'blue', 'solid', 2, '.', 0)
            if insitu_spectra_unc is not None:
                insitu_min = insitu_spectra.copy()
                insitu_max = insitu_spectra.copy()
                insitu_min[~insitu_spectra_unc.mask] = insitu_spectra[~insitu_spectra_unc.mask] - insitu_spectra_unc[
                    ~insitu_spectra_unc.mask]
                insitu_max[~insitu_spectra_unc.mask] = insitu_spectra[~insitu_spectra_unc.mask] + insitu_spectra_unc[
                    ~insitu_spectra_unc.mask]
                pspectra.plot_iqr_basic(insitu_min, insitu_max, 'red')

            if sat_spectra_unc is not None:
                sat_min = sat_spectra.copy()
                sat_max = sat_spectra.copy()
                sat_min[~sat_spectra_unc.mask] = sat_spectra[~sat_spectra_unc.mask] - sat_spectra_unc[
                    ~sat_spectra_unc.mask]
                sat_max[~sat_spectra_unc.mask] = sat_spectra[~sat_spectra_unc.mask] + sat_spectra_unc[
                    ~sat_spectra_unc.mask]
                pspectra.plot_iqr_basic(sat_min, sat_max, 'blue')
        if options_figure['type_rrs'] == 'mu_sat':
            hline2 = pspectra.plot_single_line(sat_spectra, 'blue', 'solid', 1, '.', 10)
        if options_figure['type_rrs'] == 'mu_ins':
            hline1 = pspectra.plot_single_line(insitu_spectra, 'red', 'solid', 1, '.', 10)

        pspectra.set_xaxis_title(options_figure['xlabel'])
        pspectra.set_yaxis_title(options_figure['ylabel'])



        if options_figure['title'] is not None:
            title_here = options_figure['title'] + f' MU: {index_mu}'
            pspectra.set_title(title_here)
        if options_figure['type_rrs'] == 'mu_comparison':
            pspectra.legend_options['loc'] = 'lower center'
            pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.25)
            pspectra.legend_options['ncols'] = 2
            #pspectra.set_legend_h([hline1[0], hline2[0]], ['In situ Rrs', 'Satellite Rrs'])


        ### TEMPORAL FOR PRESENTATION VALIDATIION PACE (also comment lines 2521, 2530, 2557)
        # pspectra.set_vertical_line_impl(590,-1,4.2,'red','-')
        # pspectra.set_vertical_line_impl(610, -1,4.2, 'red', '-')
        # pspectra.kk()


        if options_figure['y_min'] is not None and options_figure['y_max'] is not None:
            pspectra.set_y_range(options_figure['y_min'], options_figure['y_max'])

        pspectra.set_grid()

        pspectra.prepare_poster()

        pspectra.set_tigth_layout()
        if not options_figure['file_out'] is None:
            file_out = options_figure['file_out']
            satellite_time = self.mrfile.sat_times[index_mu].strftime('%Y%m%d')
            file_out = f'{file_out[:-4]}_{satellite_time}_{index_mu}{file_out[-4:]}'
            pspectra.save_fig(file_out)
        pspectra.close_plot()

    def plot_insmu_spectraplot(self, options_figure, index_mu):
        ##GETTING DATA
        spectra_selected, spectra_valid, spectra_invalid, instrument_idx = self.mrfile.get_mu_insitu_spectra(index_mu, options_figure[
            'scale_factor'])
        n_selected = spectra_selected.shape[0]
        n_valid = spectra_valid.shape[0]
        wl = self.mrfile.insitu_bands
        if len(wl.shape)==2:
            if instrument_idx>=0:
                wl = np.squeeze(wl[instrument_idx,:])
            else:
                print(f'[ERROR] In situ spectra plot can not be plotted, ambigous identification of the wavelengths array')
                return

        from PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()

        ##Setting wavelenth data and ticks
        pspectra.xdata = wl
        if options_figure['wlticks'] is not None:
            wlt = options_figure['wlticks']
            wls = [f'{x}' for x in wlt]
        else:
            wlt = wl
            wls = []
            for wl_value in wl:
                wlstr = f'{wl_value:.2f}'
                wlstr = wlstr.replace('.00', '')
                wls.append(wlstr)
        pspectra.set_xticks(wlt, wls, 0, 12)

        ##Plotting data

        if n_valid >= 1:
            style = {'color': 'gray', 'linestyle': '-', 'linewidth': 1, 'marker': None, 'markersize': 10}
            for idx in range(n_valid):
                hvalid = pspectra.plot_data(spectra_valid[idx, :], style)
        if n_selected == 1:
            hselected = pspectra.plot_single_line(np.squeeze(spectra_selected), 'black', 'solid', 2, None, 10)

        if n_valid == 0 and n_selected == 0:
            pspectra.close_plot()
            return

        ##title, axis ranges and axis lables
        pspectra.set_xaxis_title(options_figure['xlabel'])
        pspectra.set_yaxis_title(options_figure['ylabel'])
        if options_figure['title'] is not None:
            title_here = options_figure['title'] + f' MU: {index_mu}'
            pspectra.set_title(title_here)
        if options_figure['y_min'] is not None and options_figure['y_max'] is not None:
            pspectra.set_y_range(options_figure['y_min'], options_figure['y_max'])

        ##legend
        pspectra.legend_options['loc'] = 'lower center'
        pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.25)
        pspectra.legend_options['ncols'] = 2
        if n_valid > 0:
            pspectra.set_legend_h([hselected[0], hvalid[0]], ['In situ selected spectrum', 'Other in situ spectra'])
        if n_valid == 0:
            pspectra.set_legend_h([hselected[0]], ['In situ selected spectrum'])

        ##compleing
        pspectra.set_grid()
        pspectra.set_tigth_layout()

        ##saving
        if not options_figure['file_out'] is None:
            file_out = options_figure['file_out']
            satellite_time = self.mrfile.sat_times[index_mu].strftime('%Y%m%d')
            file_out = f'{file_out[:-4]}_{satellite_time}_{index_mu}{file_out[-4:]}'
            pspectra.save_fig(file_out)
        pspectra.close_plot()

    def plot_insitu_spectraplots(self, options_figure):
        ##GETTING DATA
        wavelength = self.mrfile.get_insitu_wl()
        all_spectra, all_spectra_validity, spectra_stats = self.mrfile.get_all_insitu_spectra(
            options_figure['scale_factor'], options_figure['use_rhow'], options_figure['plot_stats'])

        from PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        pspectra.xdata = wavelength
        for ps in options_figure['plot_spectra']:
            if ps.lower() == 'none':
                continue
            if ps.lower() == 'valid':
                spectra_valid = all_spectra[all_spectra_validity == 1]
                for spectra in spectra_valid:
                    pspectra.plot_data(spectra, options_figure['valid_line_style'])

        pspectra.set_grid()

        if options_figure['xlabel'] is not None: pspectra.set_xaxis_title(options_figure['xlabel'])
        if options_figure['ylabel'] is not None: pspectra.set_yaxis_title(options_figure['ylabel'])
        pspectra.set_tigth_layout()
        if options_figure['file_out'] is not None: pspectra.save_plot(options_figure['file_out'])
        pspectra.close_plot()
        # if not options_out['plot_stats']:
        #     stats = None

    def plot_multiple_stats_plot(self, options_figure):
        if options_figure['type_point']=='wavelength':
            df = self.get_table_spectral_statistics(options_figure)
        elif options_figure['type_point']=='matchup':
            df = self.get_table_match_up_statistics(options_figure)

        if options_figure['type_plot']=='joliff':
            col_list = df.columns.tolist()
            ini = 1 if col_list[0]=='GLOBAL' else 0
            wldata = np.array([float(col_list[idx]) for idx in range(ini,len(col_list))])
            nBias = np.array(df.loc['nBIAS'][ini:]).astype(np.float64)
            suRMSD = np.array(df.loc['suRMSD'][ini:]).astype(np.float64)
            stat_array = None
            if options_figure['colorby'] == 'stat':
                stat_str = options_figure['stat_to_color']
                stat_list = df.index.tolist()
                if stat_str not in stat_str:
                    print(f'[ERROR] stat_to_color option: {stat_str} is not in the stat list ({stat_list})')
                    return
                stat_array = np.array(df.loc[stat_str][ini:]).astype(np.float64)

            self.plot_joliff_type(options_figure,suRMSD,nBias,wldata,stat_array)

        if options_figure['type_plot']=='nbias-sam':
            nBias = np.array(df['nBIAS']).astype(np.float64)
            sam = np.array(df['SAM']).astype(np.float64)
            stat_array = None
            if options_figure['colorby'] == 'stat':
                stat_str = options_figure['stat_to_color']
                stat_list = df.index.tolist()
                if stat_str not in stat_str:
                    print(f'[ERROR] stat_to_color option: {stat_str} is not in the stat list ({stat_list})')
                    return
                stat_array = np.array(df[stat_str]).astype(np.float64)
            nBiasAbs = np.abs(nBias)
            sam_signed = sam.copy()
            sam_signed[nBias<0] = sam_signed[nBias<0]*-1
            options_figure['theta_min']=-90
            options_figure['theta_max']=90
            options_figure['scale']='linear'
            self.plot_taylor_type(options_figure,nBiasAbs,sam_signed,stat_array)


    def plot_taylor_type(self,options_figure,rdata,angledata,stat_array):
        #print(options_figure)
        from PlotScatter import PlotScatter
        ps = PlotScatter()
        ps.start_plot_polar()

        ps.set_theta_zero_location(options_figure['theta_zero_location'])
        ps.set_theta_direction(options_figure['theta_direction'])

        ps.set_rscale(options_figure['scale'])
        colors = options_figure['color']
        markersizes = options_figure['markersize']
        markers = options_figure['marker']
        edgecolors = options_figure['edgecolor']
        linewidths = options_figure['linewidth']
        print(f'[INFO] Angle data: {np.min(angledata)} - {np.max(angledata)}')
        print(f'[INFO] Radius data: {np.min(rdata)} - {np.max(rdata)}')

        if options_figure['colorby']=='density':
            from scipy.stats import gaussian_kde
            ps.set_cmap('jet')
            xy = np.vstack([angledata, rdata])
            z = gaussian_kde(xy)(xy)
            idx = z.argsort()
            angledata, rdata, z = angledata[idx], rdata[idx], z[idx]
            hscatter = ps.plot_polar(angledata,rdata,True,markers[0],markersizes[0],z,edgecolors[0],linewidths[0])

        if options_figure['colorby'] == 'stat':
            ps.set_cmap('jet')
            hscatter = ps.plot_polar(angledata, rdata, True,markers[0], markersizes[0],stat_array, edgecolors[0], linewidths[0])

        if options_figure['colorby'] == 'none':
            ps.plot_polar(angledata, rdata, True,markers[0],markersizes[0],colors[0],edgecolors[0],linewidths[0])

        if options_figure['rlim'] is not None:
            ps.set_rlim(options_figure['rlim'])


        #rticks_pos = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]

        rticks_pos = np.array(ps.ax.get_yticks())
        rlabels_empty = [f'' for x in rticks_pos]
        ps.set_rticks_and_labels(rticks_pos, rlabels_empty, False)

        if options_figure['rlim'] is not None:
            rmin = options_figure['rlim'][1]*-1
            rmax = options_figure['rlim'][1]
            all_rticks = np.linspace(rmin,rmax,13)
            all_rticks_labels = ['{:g}'.format(x) for x in all_rticks]
            all_pos = np.linspace(0,1,13)
            ps.plot_text_options['fontsize'] = 10
            ps.plot_text_options['horizontalalignment'] = 'center'
            for label,pos in zip(all_rticks_labels,all_pos):
                ps.plot_text(pos, 0.21, label)


        ps.set_theta_range(options_figure['theta_min'], options_figure['theta_max'])
        if options_figure['colorbar']:
            ps.colorbar(hscatter)
            options_figure['legend']=False
        ps.set_equal_apect()
        ps.tight_layout()
        if options_figure['file_out'] is not None:
            ps.save_fig(options_figure['file_out'])
        ps.close_plot()


    def plot_joliff_type(self,options_figure,xdata,ydata,wldata,stat_array):
        #print(options_figure)
        #from matplotlib import pyplot as plt
        from PlotScatter import PlotScatter
        ps = PlotScatter()
        ps.start_plot()
        min_x_abs = np.min(xdata)
        min_y_abs = np.min(ydata)
        max_x_abs = np.max(xdata)
        max_y_abs = np.max(ydata)
        min_xy_abs = np.min([min_x_abs,min_y_abs])
        max_xy_abs = np.max([min_y_abs,max_y_abs])
        colors = options_figure['color']
        markersizes = options_figure['markersize']
        markers = options_figure['marker']
        edgecolors = options_figure['edgecolor']
        linewidths = options_figure['linewidth']
        handles = []

        if options_figure['colorby']=='density':
            from scipy.stats import gaussian_kde
            ps.set_cmap('jet')
            xy = np.vstack([xdata, ydata])
            z = gaussian_kde(xy)(xy)
            idx = z.argsort()
            xdata, ydata, z = xdata[idx], ydata[idx], z[idx]
            hscatter = ps.plot_data(xdata,ydata,markers[0],markersizes[0],z,edgecolors[0],linewidths[0])

        if options_figure['colorby']=='wavelength':
            ps.set_cmap('jet')
            hscatter = ps.plot_data(xdata, ydata, markers[0], markersizes[0], wldata, edgecolors[0], linewidths[0])

        if options_figure['colorby']=='wavelength_ranges':
            if options_figure['wlranges_min'] is None or options_figure['wlranges_max'] is None:
                ps.close_plot()
                print(f'[ERROR] wlranges_min and wlranges_max should be defined for colorby: wavelenth_ranges')
                return
            nranges,garray = self.get_basic_wl_ranges_info(options_figure,wldata)
            colors = options_figure['color']
            for irange in range(nranges):
                xdata_here = xdata[garray==irange]
                ydata_here = ydata[garray==irange]
                hscatter = ps.plot_data(xdata_here,ydata_here,markers[0],markersizes[0],colors[irange],edgecolors[0],linewidths[0])
                handles.append(hscatter)

        if options_figure['colorby'] == 'stat':
            ps.set_cmap('jet')
            hscatter = ps.plot_data(xdata, ydata, markers[0], markersizes[0],stat_array, edgecolors[0], linewidths[0])


        if options_figure['min_xy'] is not None or options_figure['max_xy'] is not None:
            min_xy = options_figure['min_xy'] if options_figure['min_xy'] is not None else min_xy_abs
            max_xy = options_figure['max_xy'] if options_figure['max_xy'] is not None else max_xy_abs
            ps.set_limits(min_xy,max_xy)

        if options_figure['x_min'] is not None or options_figure['x_max'] is not None:
            min_x = options_figure['x_min'] if options_figure['x_max'] is not None else min_x_abs
            max_x = options_figure['x_max'] if options_figure['x_max'] is not None else max_x_abs
            ps.set_limits_X(min_x,max_x)

        if options_figure['y_min'] is not None or options_figure['y_max'] is not None:
            min_y = options_figure['y_min'] if options_figure['y_min'] is not None else min_y_abs
            max_y = options_figure['y_max'] if options_figure['y_max'] is not None else max_y_abs
            ps.set_limits_Y(min_y,max_y)

        if options_figure['ticks'] is not None:
            x_ticks = options_figure['ticks']
            y_ticks = options_figure['ticks']
        else:
            if options_figure['x_ticks'] is None:
                x_ticks = self.get_ticks_from_min_max_xy(min_x_abs, max_x_abs)
            else:
                x_ticks = options_figure['x_ticks']
            if options_figure['y_ticks'] is None:
                y_ticks = self.get_ticks_from_min_max_xy(min_y_abs, max_y_abs)
            else:
                y_ticks = options_figure['y_ticks']

        ps.set_ticks_x(x_ticks, options_figure['fontsizeaxis'])
        ps.set_ticks_y(y_ticks, options_figure['fontsizeaxis'])

        # plt.grid(which='major', color='gray', linestyle=':', axis='both')
        # if extra_lines[0]:
        #     plt.plot([0, 0], [0, y_ranges[1]], color='k', marker=None, linewidth=1, linestyle='-')
        # if extra_lines[1]:
        #     plt.plot([x_ranges[0], x_ranges[1]], [1, 1], color='k', marker=None, linewidth=1, linestyle='-')

        if options_figure['xlabel'] is not None:
            xtitle = options_figure['xlabel']
        else:
            xtitle = 'suRMSD'
        if xtitle is not None:
            if options_figure['colorbar']:
                ps.plot_text(0.85, 0.45, xtitle)
            else:
                ps.plot_text(1.02,0.49,xtitle)

        if options_figure['ylabel'] is not None:
            ytitle = options_figure['ylabel']
        else:
            ytitle = 'nBIAS'
        if ytitle is not None:
            ps.plot_text(0.44, 1.02, 'nBIAS')


        # if xstat == 'DETER(r2)' or ystat == 'DETER(r2)':
        #     plt.legend(handles, col_names, ncol=3, frameon=True, framealpha=1,
        #                loc='upper center')  # ,bbox_to_anchor=(0.5,-0.4))
        # else:
        #     plt.legend(handles, col_names, loc='lower right', bbox_to_anchor=(1.9, 0))
        # if legend_info is not None:
        #     plt.legend(handles, col_names, loc=legend_info['loc'], bbox_to_anchor=legend_info['bbox_to_anchor'])

        c_ticks = x_ticks if options_figure['circular_spines'] else None
        ps.set_as_joliff(c_ticks)
        if options_figure['colorbar']:
            ps.colorbar(hscatter)
            options_figure['legend']=False
        if options_figure['legend']:
            ps.set_legend_h(handles,options_figure['legend_values'])


        ps.set_equal_apect()
        ps.tight_layout()
        if options_figure['file_out'] is not None:
            ps.save_fig(options_figure['file_out'])
        ps.close_plot()



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
                print(f'[INFO] Using virtual flag: {var_group_name}')
                flag_values = self.virtual_flags[var_group_name]['flag_values']
                flag_meanings = self.virtual_flags[var_group_name]['flag_meanings']
                options_figure[var_group_name] = {
                    'flag_values': flag_values,
                    'flag_meanings': flag_meanings
                }

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
                        print(f'[WARNING] Flag {flag_config.strip()} is not in the list.')
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

    def set_data_scatterplot_general(self, groupBy, selectBy, selectValues, options_out):
        self.xdata = []
        self.ydata = []
        xarray = self.mrfile.nc.variables[options_out['xvar']][:]
        yarray = self.mrfile.nc.variables[options_out['yvar']][:]
        if len(xarray.shape) > 1 or len(yarray.shape) > 1:
            return
        ndata = xarray.shape[0]
        if ndata != yarray.shape[0]:
            return

        valid_all = np.ones(xarray.shape)
        valid_all[xarray.mask] = 0
        valid_all[yarray.mask] = 0
        valid_all[np.isnan(xarray)] = 0
        valid_all[np.isnan(yarray)] = 0

        # for idx in range(len(xarray)):
        #     if valid_all[idx] == 1:
        #         perc = 100 * ((yarray[idx] - xarray[idx]) / xarray[idx])
        #         perc = abs(perc)
        #         if perc > 1000:
        #             valid_all[idx] = 0
        #         print(idx, ';', xarray[idx], ';', yarray[idx], ';', perc)

        # #smp wfr
        ##valid_all[27]=0

        if self.mu_valid_variable in self.mrfile.nc.variables:
            #print('===============================================================================>')
            mu_valid = self.mrfile.nc.variables[self.mu_valid_variable][:]
            valid_all = np.array(mu_valid)

        if selectBy is not None:
            select_array, all_select_values, all_select_meanings = self.get_flag_array(options_out, 'selectBy')
            # if len(select_array.shape) == 1 and select_array.shape[0] == ndata:
            #     select_data = select_array
            if len(select_array.shape) == 2 and select_array.shape[0] == ndata:
                select_array_1D = np.squeeze(select_array[:, 0])
                var_id = f'{options_out["xvar"]}_id'
                if var_id in self.mrfile.nc.variables:
                    var_id_array = np.array(self.mrfile.nc.variables[var_id])
                    max_id = np.max(var_id_array[:])
                    for id in range(1, max_id + 1):
                        select_array_1D[var_id_array == id] = select_array[var_id_array == id, id]
                select_array = select_array_1D

            valid_all_s = np.zeros(xarray.shape)
            try:
                for val in selectValues:
                    valid_all_s[np.logical_and(select_array == val, valid_all == 1)] = 1
            except:  ##selectValues is a single non-iterable value
                val = selectValues
                valid_all_s[np.logical_and(select_array == val, valid_all == 1)] = 1
            valid_all = valid_all_s

        self.xdata = np.array(xarray[valid_all == 1])
        self.ydata = np.array(yarray[valid_all == 1])

        if groupBy is not None:
            self.groupdata = []
            group_array, all_group_values, all_group_meanings = self.get_flag_array(options_out, 'groupBy')
            if len(group_array.shape) == 1 and group_array.shape[0] == ndata:
                self.groupdata = group_array[valid_all == 1]
            if len(group_array.shape) == 2 and group_array.shape[0] == ndata:
                group_array_1D = np.squeeze(group_array[:, 0])
                var_id = f'{options_out["xvar"]}_id'
                if var_id in self.mrfile.nc.variables:
                    var_id_array = np.array(self.mrfile.nc.variables[var_id])
                    max_id = np.max(var_id_array[:])
                    for id in range(1, max_id + 1):
                        group_array_1D[var_id_array == id] = group_array[var_id_array == id, id]
                self.groupdata = group_array_1D[valid_all == 1]

            # if len(group_array) == len(mu_valid_satelliteid):
            #     group_array = self.get_array_muid_from_array_satelliteid(id_mu, group_array)
            # self.groupdata = group_array[valid_all == 1]

    def set_data_scatterplot(self, groupBy, selectBy, valSelect, wl_value, options_out):
        rrs_ins = np.array(self.mrfile.nc.variables['mu_ins_rrs'])
        rrs_sat = np.array(self.mrfile.nc.variables['mu_sat_rrs'])
        id_mu = np.array(self.mrfile.nc.variables['mu_satellite_id'])

        mu_valid_variable = self.global_options['mu_valid_variable']
        mu_valid_satelliteid = np.array(self.mrfile.nc.variables[mu_valid_variable])

        valid_mu = self.get_array_muid_from_array_satelliteid(id_mu, mu_valid_satelliteid)

        valid_all = self.check_rrs_valid(valid_mu, rrs_ins, rrs_sat)

        if wl_value is not None:
            print(f'[INFO] Wavelength value: {wl_value}')
            wl_array = np.array(self.mrfile.nc.variables['mu_wavelength'])
            valid_all[wl_array != wl_value] = 0

        if wl_value is None and options_out['wl_min'] is not None:
            wl_array = np.array(self.mrfile.nc.variables['mu_wavelength'])
            valid_all[wl_array<options_out['wl_min']]=0

        if wl_value is None and options_out['wl_max'] is not None:
            wl_array = np.array(self.mrfile.nc.variables['mu_wavelength'])
            valid_all[wl_array>options_out['wl_max']]=0

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

        ##TEMPORAL POSTER
        # self.xdata = self.xdata * 1000
        # self.ydata = self.ydata * 1000

        if groupBy is not None:
            if options_out['groupType'] == 'float' or options_out['groupType'] == 'wavelength':
                group_array = np.array(self.mrfile.nc.variables[groupBy])
            else:
                group_array, all_group_values, all_group_meanings = self.get_flag_array(options_out, 'groupBy')
            if len(group_array) == len(mu_valid_satelliteid):
                group_array = self.get_array_muid_from_array_satelliteid(id_mu, group_array)
            self.groupdata = group_array[valid_all == 1]




    def get_basic_wl_ranges_info(self,options_figure,wldata):
        wlmin_values = options_figure['wlranges_min']
        wlmax_values = options_figure['wlranges_max']
        nranges = len(wlmin_values)
        ndata = len(wldata)
        garray = np.zeros((ndata,))
        colors_wl = []
        flag_meanings = []
        for idx in range(nranges):
            min_wl = wlmin_values[idx]
            max_wl = wlmax_values[idx]
            center_wl = min_wl + ((max_wl - min_wl) / 2)
            color_wl = defaults.get_color_wavelength(center_wl)
            colors_wl.append(color_wl)
            garray[np.logical_and(wldata >= min_wl, wldata < max_wl)] = idx
            flag_meaning = f'{min_wl}-{max_wl}'
            flag_meanings.append(flag_meaning)

        options_figure['legend_values'] = flag_meanings
        color_prev = options_figure['color']
        if len(color_prev) != len(wlmin_values):
            options_figure['color'] = colors_wl

        return nranges,garray

    def create_flag_array_wl_ranges(self, options_figure, name_fv):

        wlmin_values = options_figure['wlranges_min']
        wlmax_values = options_figure['wlranges_max']

        nranges = len(wlmin_values)
        mu_wavelength = self.mrfile.variables['mu_wavelength'][:]
        nmu = mu_wavelength.shape[0]
        array = np.zeros((nmu,))
        flag_values = []
        flag_meanings = []
        colors_wl = []
        for idx in range(nranges):
            min_wl = wlmin_values[idx]
            max_wl = wlmax_values[idx]
            center_wl = min_wl + ((max_wl - min_wl) / 2)
            color_wl = defaults.get_color_wavelength(center_wl)

            colors_wl.append(color_wl)
            flag_value = 2 ** idx
            flag_meaning = f'{min_wl}-{max_wl}'
            array[np.logical_and(mu_wavelength >= min_wl, mu_wavelength < max_wl)] = flag_value
            flag_values.append(flag_value)
            flag_meanings.append(flag_meaning)

        # options_figure[name_fv] = {
        #
        #     'flag_values':flag_values,
        #     'flag_meanings':flag_meanings
        # }
        self.virtual_flags[name_fv] = {
            'flag_array': array,
            'flag_values': flag_values,
            'flag_meanings': flag_meanings
        }
        options_figure['legend_values'] = flag_meanings
        color_prev = options_figure['color']
        if len(color_prev) != len(wlmin_values):
            options_figure['color'] = colors_wl

        return options_figure

    # options_var: selectBy or groupBy
    def get_flag_array(self, options_out, option_var):
        var_flag = options_out[option_var]
        if var_flag in self.mrfile.variables:
            array_flag = np.array(self.mrfile.variables[var_flag][:])
            if 'flag_values' in self.mrfile.variables[var_flag].ncattrs():
                flag_values = self.mrfile.variables[var_flag].flag_values
            elif 'flag_masks' in self.mrfile.variables[var_flag].ncattrs():
                flag_values = self.mrfile.variables[var_flag].flag_masks
            flag_meanings = self.mrfile.variables[var_flag].flag_meanings.split(' ')
        else:  ##previously virtual flag
            array_flag = self.virtual_flags[var_flag]['flag_array']
            flag_values = self.virtual_flags[var_flag]['flag_values']
            flag_meanings = self.virtual_flags[var_flag]['flag_meanings']

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
        allValues = np.array(allValues)
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
        elif wl_sat_value_str.endswith('0') and wl_sat_value_str.find('.')>0:
            return wl_sat_value_str[:-1]
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
                # print(stat.upper(), val)

                valstr = self.get_str_stat(val, options[f'{stat.upper()}_FORMAT'], options['units'])
                if stat.upper() == 'APD' or stat.upper() == 'RPD':
                    valstr = f'{valstr}%'
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
                        value_a = f'{g:.2f}'
                        if options['groupType'] == 'wavelength':
                            value_a = value_a.replace('.30', '.25')
                            value_a = value_a.replace('.80', '.75')

                        str_legend.append(value_a)
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

    def create_virtual_flag(self, poptions, vflag):
        if vflag in self.virtual_flags:
            return
        from FlagBuilder import FlagBuilder
        print(f'[INFO] Creating virtual flag: {vflag}')
        fbuilder = FlagBuilder(self.path_mdbr_file, None, poptions.options)
        flag_values, flag_names, array = fbuilder.create_flag_array(vflag, False)

        self.virtual_flags[vflag] = {
            'flag_array': array,
            'flag_values': flag_values,
            'flag_meanings': flag_names
        }
