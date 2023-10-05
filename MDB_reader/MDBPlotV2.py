import netCDF4

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

        self.mrfile = None
        self.path_mdbr_file = path_mdbr_file
        self.VALID = False
        if path_mdbr_file is not None:
            self.mrfile = MDBFile(path_mdbr_file)
            self.VALID = self.mrfile.VALID

        # plot general options
        self.title = ''
        self.file_name_base = ''
        self.format_image = 'png'

        self.output_path = None

        self.mu_valid_variable = 'mu_valid'

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

        self.global_stats_ins_spectra = {}

    def set_global_options(self, options):
        section = 'GLOBAL_OPTIONS'

        self.output_path = self.get_value_param(options, section, 'output_path', self.output_path, 'directory')
        self.mu_valid_variable = self.get_value_param(options, section, 'mu_valid_variable', self.mu_valid_variable,
                                                      'str')

    def plot_from_options(self, options):
        plot_list = list(options.sections())
        for plot in plot_list:
            if plot == 'GLOBAL_OPTIONS':
                continue
            options_out = self.get_options(options, plot)

            print(options_out)
            if options_out['apply']:
                self.plot_from_options_impl(options_out)

    def plot_from_options_impl(self, options_out):

        if options_out['type'] == 'scatterplot':
            file_out_base = options_out['file_out']
            title_base = options_out['title']
            if options_out['selectByWavelength']:  # one scatterplot for wavelenght
                if options_out['multiple_plot'] is not None:
                    self.plot_multiple_scatter_plots(options_out)

                if options_out['selectBy'] is None:  # one scatterplot global for wavelentht
                    if options_out['multiple_plot'] is None:

                        for wl in options_out['wl_values']:
                            self.set_data_scatterplot(options_out['groupBy'], None, None, wl,options_out)
                            options_out['file_out'] = self.get_file_out_name(file_out_base, wl, None)
                            options_out['title'] = self.get_title(title_base, wl, None, None)
                            if len(self.xdata) > 0 and len(self.ydata) > 0:
                                self.plot_scatter_plot(options_out, None, -1, wl)
                    else:
                        rc = options_out['multiple_plot'].split(',')
                        nrow = int(rc[0].strip())
                        ncol = int(rc[1].strip())
                        ntot = nrow * ncol
                        index = 0
                        from PlotScatter import PlotScatter
                        plot_here = PlotScatter()
                        plot_here.xtitle_options['fontsize'] = options_out['fontsizelabels']
                        plot_here.ytitle_options['fontsize'] = options_out['fontsizelabels']
                        plot_here.plot_text_options['fontsize'] = options_out['fontsizestats']

                        plot_here.start_multiple_plot_advanced(nrow, ncol, options_out['xfigsize'],
                                                               options_out['yfigsize'], options_out['widthspace'],
                                                               options_out['heightspace'])
                        print(f'[INFO] Starting multiple plot with {nrow} rows and {ncol} cols')
                        wl_values = options_out['wl_values']
                        print(f'[INFO] Wavelenthts: {wl_values}')
                        for wl in options_out['wl_values']:
                            self.set_data_scatterplot(options_out['groupBy'], None, None, wl,options_out)
                            if len(self.xdata) > 0 and len(self.ydata) > 0:
                                options_out['title'] = self.get_title(title_base, wl, None, None)
                                options_out['file_out'] = None

                                self.plot_scatter_plot(options_out, plot_here, index, wl)
                                index = index + 1
                            else:
                                print(f'[WARNING] No data for wavelength: {wl} nm')

                        for index_blank in range(index, ntot):
                            plot_here.plot_blanck(index_blank)
                        if options_out['legend']:
                            str_legend = self.get_str_legend(options_out)
                            if len(str_legend) > 0:
                                plot_here.set_global_legend(str_legend)
                        file_out = self.get_file_out_name(file_out_base, None, None)
                        plot_here.save_fig(file_out)
                        plot_here.close_plot()
                else:  # one scatterplot per flag and wavelenght
                    flag_name = options_out['selectBy']
                    if options_out['multiple_plot'] is None:
                        for wl in options_out['wl_values']:
                            for flag_value in options_out['selectValues']:

                                self.set_data_scatterplot(options_out['groupBy'], options_out['selectBy'], flag_value,
                                                          wl,options_out)
                                flag_here = self.get_flag_flag(flag_value, options_out[flag_name]['flag_values'],
                                                               options_out[flag_name]['flag_meanings'])
                                options_out['file_out'] = self.get_file_out_name(file_out_base, wl, flag_here)
                                # print(wl, flag_value, options_out['selectBy'], options_out['file_out'])
                                options_out['title'] = self.get_title(title_base, wl, flag_here, None)
                                if len(self.xdata) > 0 and len(self.ydata) > 0:
                                    self.plot_scatter_plot(options_out, None, -1, wl)

            if not options_out['selectByWavelength']:
                if options_out['selectBy'] is None:
                    self.set_data_scatterplot(options_out['groupBy'], None, None, None,options_out)
                    self.plot_scatter_plot(options_out, None, -1, -1)
                else:
                    flag_name = options_out['selectBy']
                    if options_out['multiple_plot'] is None:
                        for flag_value in options_out['selectValues']:
                            self.set_data_scatterplot(options_out['groupBy'], options_out['selectBy'], flag_value, None)
                            flag_here = self.get_flag_flag(flag_value, options_out[flag_name]['flag_values'],
                                                           options_out[flag_name]['flag_meanings'],options_out)
                            options_out['file_out'] = self.get_file_out_name(file_out_base, None, flag_here)

                            options_out['title'] = self.get_title(title_base, None, flag_here, None)
                            if len(self.xdata) > 0 and len(self.ydata) > 0:
                                self.plot_scatter_plot(options_out, None, -1, None)

        if options_out['type'] == 'statstable_wl':
            self.create_table_stats_wl(options_out)

        if options_out['type'] == 'statstable':
            self.create_table_stats(options_out, None)
            for wl in options_out['wl_values']:
                # print('wl-->', wl)
                self.create_table_stats(options_out, wl)

        if options_out['type'] == 'csvtable':
            self.create_csv_table(options_out)

        if options_out['type'] == 'statswlplot':
            self.plot_statistics_bywl(options_out)

        if options_out['type'] == 'spectraplot':
            self.plot_spectra_plot(options_out)

        if options_out['type'] == 'flagplot':
            self.plot_flag_fromoptions(options_out)

        if options_out['type'] == 'temporalheatmap':
            self.plot_temporal_heat_map_fromoptions(options_out)

        if options_out['type'] == 'validityplot':
            self.plot_validityplot_fromoptions(options_out)

        if options_out['type'] == 'multipleplot':
            self.plot_multiple_plot_from_options(options_out)

    def plot_validityplot_fromoptions(self, options_out):



        flag = options_out['flag']
        flag_list = options_out['flag_list']
        if 'mu_valid_common' in self.mrfile.variables:
            varvalid = 'mu_valid_common'
        else:
            varvalid = 'mu_valid'

        fcsv = None
        if options_out['file_csv'] is not None:
            fcsv = os.path.join(os.path.dirname(self.path_mdbr_file),options_out['file_csv'])
            if not os.path.isfile(fcsv):
                fcsv = None
        if fcsv is None:
            ntotal, nvalid, porc = self.mrfile.analyse_mu_flag(flag, flag_list, varvalid)
            self.plot_validityplot_impl(flag_list, ntotal, nvalid, porc, options_out)
        else:
            self.plot_validityplot_fromcsv(fcsv,flag,flag_list,options_out)

    def plot_validityplot_fromcsv(self,fcsv,flag_name,flag_list,options_out):
        file_out = options_out['file_out']
        import pandas as pd
        df = pd.read_csv(fcsv,sep=';')
        print(df)
        print(flag_list)
        output_type = options_out['output_type']
        print(output_type)
        from matplotlib import pyplot as plt
        plt.close()
        from matplotlib import pyplot as plt
        # hfig = plt.figure(layout='tight')
        hfig, ax = plt.subplots()
        height = 0.2
        xval = 0
        xticks_pos = []
        xticks_pos_minor = []
        nseriesplot = len(output_type)

        colors = options_out['series_color']
        if colors is None:
            colors = defaults.get_color_list(nseriesplot)
        legend = options_out['series_flag']
        heightbyseries = height / nseriesplot
        handles = []
        hbarvalidity = []
        hbarvalidity_values = []

        for fl in enumerate(flag_list):
            flag = fl[1]
            xvalini = xval - (heightbyseries / 2)
            xticks_pos_minor.append(xvalini)
            icolor = 0
            for ot in output_type:
                rslt_df = df.loc[df[flag_name] == flag]
                val = rslt_df.iloc[0].at[ot]
                print(flag,ot,val)
                hbar = plt.barh(xval, val, height=heightbyseries, color=colors[icolor])
                if fl[0] == 0:
                    handles.append(hbar)
                if ot=='Valid':
                    # invalid_insitu = rslt_df.iloc[0].at['INVALID IN SITU DATA']
                    # invalid_sat = rslt_df.iloc[0].at['INVALID SAT DATA']
                    # hbar = plt.barh(xval, invalid_insitu, left=val, height=heightbyseries, color='m')
                    # hbar = plt.barh(xval, invalid_sat, left=val+invalid_insitu, height=heightbyseries, color='r')
                    hbarvalidity.append(hbar)
                    porc = (val / rslt_df.iloc[0].at['Total_mu'])*100
                    hbarvalidity_values.append(porc)
                xval = xval + heightbyseries
                icolor = icolor + 1
            xvalfin = xval - (heightbyseries / 2)
            xticks_pos_minor.append(xvalfin)
            xticks_pos.append((xvalini + xvalfin) / 2)

        plt.yticks(xticks_pos, flag_list)
        ax.set_yticks(xticks_pos_minor, minor=True)
        plt.grid(b=True, which='minor', color='gray', linestyle='--', axis='y')
        plt.grid(b=True, which='major', color='gray', linestyle='--', axis='x')
        ax.tick_params(which='major', length=0, axis='y')
        ax.tick_params(which='minor', length=10, axis='y')

        for idx in range(len(hbarvalidity)):
            hbar = hbarvalidity[idx]
            val = hbarvalidity_values[idx]
            vals = [f'{val:.2f} %']

            plt.bar_label(hbar, vals)

        plt.xlabel('Number of match-ups', fontsize=12)

        plt.gcf().tight_layout()
        if legend is not None:
            handles.reverse()
            legend.reverse()
            plt.legend(handles, legend, framealpha=1)
        if file_out is not None:
            plt.savefig(file_out, dpi=300)
        plt.close(hfig)


    def plot_validityplot_impl(self, flag_list, ntotal, nvalid, porc, options_out):
        file_out = options_out['file_out']
        output_type = options_out['output_type']
        from matplotlib import pyplot as plt
        plt.close()
        from matplotlib import pyplot as plt
        # hfig = plt.figure(layout='tight')
        hfig, ax = plt.subplots()
        height = 0.2
        xval = 0
        xticks_pos = []
        xticks_pos_minor = []
        nseriesplot = len(output_type)
        if 'valid_ac' in output_type:
            nseriesplot = (nseriesplot - 1) + (nvalid.shape[1] - 1)
        colors = options_out['series_color']
        if colors is None:
            colors = defaults.get_color_list(nseriesplot)
        legend = options_out['series_flag']
        heightbyseries = height / nseriesplot
        handles = []
        hbarvalidity = []
        hbarvalidity_values = []

        for fl in enumerate(flag_list):
            xvalini = xval - (heightbyseries / 2)
            xticks_pos_minor.append(xvalini)
            icolor = 0
            for ot in enumerate(output_type):
                if ot[1] == 'total':
                    hbar = plt.barh(xval, ntotal[fl[0], 0], height=heightbyseries, color=colors[icolor])
                    if fl[0] == 0:
                        handles.append(hbar)
                    xval = xval + heightbyseries
                    icolor = icolor + 1
                if ot[1] == 'valid':
                    # print('valid->', colors[icolor], nvalid[fl[0], 0])
                    hbar = plt.barh(xval, nvalid[fl[0], 0], height=heightbyseries, color=colors[icolor])
                    if fl[0] == 0:
                        handles.append(hbar)
                    hbarvalidity.append(hbar)
                    hbarvalidity_values.append(porc[fl[0], 0])
                    xval = xval + heightbyseries
                    icolor = icolor + 1
                if ot[1] == 'valid_ac':
                    for iac in range(1, nvalid.shape[1]):
                        # print('valid_ac->', iac, colors[icolor], nvalid[fl[0], iac])
                        hbar = plt.barh(xval, nvalid[fl[0], iac], height=heightbyseries, color=colors[icolor])
                        if fl[0] == 0:
                            handles.append(hbar)
                        hbarvalidity.append(hbar)
                        hbarvalidity_values.append(porc[fl[0], iac])
                        xval = xval + heightbyseries
                        icolor = icolor + 1

            xvalfin = xval - (heightbyseries / 2)
            xticks_pos_minor.append(xvalfin)
            xticks_pos.append((xvalini + xvalfin) / 2)

            # xval = xval + heightbyseries/4

        plt.yticks(xticks_pos, flag_list)
        ax.set_yticks(xticks_pos_minor, minor=True)
        plt.grid(b=True, which='minor', color='gray', linestyle='--', axis='y')
        plt.grid(b=True, which='major', color='gray', linestyle='--', axis='x')
        ax.tick_params(which='major', length=0, axis='y')
        ax.tick_params(which='minor', length=10, axis='y')

        for idx in range(len(hbarvalidity)):
            hbar = hbarvalidity[idx]
            val = hbarvalidity_values[idx]
            vals = [f'{val:.2f} %']

            plt.bar_label(hbar, vals)

        plt.xlabel('Number of match-ups', fontsize=12)

        plt.gcf().tight_layout()
        if legend is not None:
            handles.reverse()
            legend.reverse()
            plt.legend(handles, legend, framealpha=1)
        if file_out is not None:
            plt.savefig(file_out, dpi=300)
        plt.close(hfig)

    def plot_temporal_heat_map_fromoptions(self, options_out):
        # print(options_out)
        file_out = options_out['file_out']
        files_multiple = []
        if options_out['flag'] is not None:
            flag = options_out['flag']
            flag_list = options_out['flag_list']
            ext = file_out[-4:]
            file_base = file_out[:-4]
            output_type = options_out['output_type']
            multiple_plot = False
            if len(output_type) > 1:
                multiple_plot = True
            index_fig = 1
            for otype in output_type:
                onlyvalid = False
                varvalid = None
                if otype == 'valid':
                    onlyvalid = True
                    varvalid = 'mu_valid'
                if otype == 'valid_common':
                    onlyvalid = True
                    varvalid = 'mu_valid_common'
                if otype == 'total_mu':
                    onlyvalid = False
                    varvalid = 'mu_valid'
                dfdata = self.mrfile.analyse_mu_temporal_flag(onlyvalid, varvalid, flag, flag_list)
                if multiple_plot:
                    file_out_here = f'{file_base}_{index_fig}{ext}'
                    files_multiple.append(file_out_here)
                else:
                    file_out_here = file_out
                self.plot_heat_map_impl(dfdata, file_out_here, options_out)
                index_fig = index_fig + 1

        if options_out['multiple_plot'] is not None:
            from PlotMultiple import PlotMultiple
            rc = options_out['multiple_plot'].split(',')
            nrow = int(rc[0].strip())
            ncol = int(rc[1].strip())
            ntot = nrow * ncol
            if ntot == len(files_multiple):
                pm = PlotMultiple()
                xfigsize = options_out['xfigsize']
                yfigsize = options_out['yfigsize']
                wspace = options_out['widthspace']
                hspace = options_out['heightspace']
                pm.start_multiple_plot_advanced(nrow, ncol, xfigsize, yfigsize, wspace, hspace, False)
                index = 0
                for irow in range(nrow):
                    for icol in range(ncol):
                        pm.plot_image(files_multiple[index], irow, icol)
                        index = index + 1

                pm.set_text(1800, -50, '(a)')
                pm.set_text(1800, 1400, '(b)')
                # file_out = self.get_file_out_flag_param(options_out['file_out'], None, None)
                pm.save_fig(file_out)

    def plot_heat_map_impl(self, dfdata, file_out, options_out):
        from matplotlib import cm
        from matplotlib import pyplot as plt
        import seaborn as sns
        ext = file_out[-4:]
        plt.close()
        h = plt.Figure()
        cmap = cm.get_cmap('RdYlBu_r')
        dfall_month_withnan = dfdata.copy()
        dfall_month_withnan[dfdata == 0] = np.nan
        vmin = options_out['vmin']
        vmax = options_out['vmax']
        sns.heatmap(dfall_month_withnan, vmin = vmin, vmax = vmax, annot=False, cmap=cmap, linewidths=1, linecolor='black')
        plt.xlabel('Year-Month')
        plt.ylabel('Site')
        plt.gcf().tight_layout()
        if file_out is not None:
            plt.savefig(file_out, dpi=300)
            file_out_csv = file_out.replace(ext, '.csv')
            dfall_month_withnan.to_csv(file_out_csv, sep=';')
        plt.close(h)

    def plot_multiple_plot_from_options(self, options_out):
        from PlotMultiple import PlotMultiple
        rc = options_out['multiple_plot'].split(',')
        files_multiple = options_out['multiple_files']
        file_out = options_out['file_out']
        nrow = int(rc[0].strip())
        ncol = int(rc[1].strip())
        ntot = nrow * ncol

        if ntot == len(files_multiple):
            pm = PlotMultiple()
            xfigsize = options_out['xfigsize']
            yfigsize = options_out['yfigsize']
            wspace = options_out['widthspace']
            hspace = options_out['heightspace']
            pm.start_multiple_plot_advanced(nrow, ncol, xfigsize, yfigsize, wspace, hspace, False)
            index = 0
            for irow in range(nrow):
                for icol in range(ncol):
                    pm.plot_image(files_multiple[index], irow, icol)
                    index = index + 1
            if ntot == 4:
                pm.set_text(-150, -50, '(a)')
                pm.set_text(1800, -50, '(b)')
                pm.set_text(-150, 1400, '(c)')
                pm.set_text(1800, 1400, '(d)')
            if ntot == 2:
                pm.set_text(1800, -50, '(a)')
                pm.set_text(1800, 1400, '(b)')
            # if len(legend) > 0 and len(legend) == len(handles):
            #     pm.set_global_legend(handles, legend)
            # file_out = self.get_file_out_flag_param(options_out['file_out'], None, None)
            pm.save_fig(file_out)

    def create_csv_table(self, options):
        potential_valiables = ['mu_id', 'mu_valid', 'mu_valid_common', 'flag_ac', 'flag_satellite', 'flag_sensor',
                               'flag_site',
                               'mu_sat_time', 'mu_ins_time', 'mu_time_diff']
        name_col = ['mu_id']
        norig = 1
        for varname in potential_valiables:
            if varname in self.mrfile.nc.variables:
                name_col.append(varname)
                norig = norig + 1

        # if 'mu_valid_common' in self.mrfile.nc.variables:
        #     name_col = ['mu_id', 'mu_valid', 'mu_valid_common', 'flag_ac', 'flag_satellite', 'flag_sensor', 'flag_site',
        #                 'mu_sat_time', 'mu_ins_time', 'mu_time_diff']
        #     norig = 10
        # else:
        #     name_col = ['mu_id', 'mu_valid', 'flag_ac', 'flag_satellite', 'flag_sensor', 'flag_site', 'mu_sat_time',
        #                 'mu_ins_time', 'mu_time_diff']
        #     norig = 9
        wl_data = np.array(self.mrfile.nc.variables['mu_wavelength'])
        wavelenghts = np.unique(wl_data)
        for w in wavelenghts:
            ws = f'{w:0.2f}'
            name_col.append(f'sat_{ws}')
            name_col.append(f'insitu_{ws}')
        nmu = self.mrfile.n_mu_total
        file_out = options['file_out']
        df_validation = pd.DataFrame(columns=name_col, index=list(range(nmu)))

        df_validation['mu_id'] = np.arange(0, nmu)
        for idx in range(1, norig):
            name_here = name_col[idx]
            df_validation[name_here] = np.array(self.mrfile.nc.variables[name_here])

        sat_data = np.array(self.mrfile.nc.variables['mu_sat_rrs'])
        ins_data = np.array(self.mrfile.nc.variables['mu_ins_rrs'])
        mu_id = np.array(self.mrfile.nc.variables['mu_satellite_id'])
        for idx in range(len(sat_data)):
            wl_here = wl_data[idx]
            sat_here = sat_data[idx]
            ins_here = ins_data[idx]
            index_here = mu_id[idx]
            ws = f'{wl_here:0.2f}'
            sat_col = f'sat_{ws}'
            ins_col = f'insitu_{ws}'
            # print(index_here,sat_col,ins_col)
            df_validation.loc[index_here, sat_col] = sat_here
            df_validation.loc[index_here, ins_col] = ins_here
        # print(df_validation)
        # print(file_out)
        df_validation.to_csv(file_out, sep=';')

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
            title = title.replace('$FLAG$', flag)

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

    def set_data_scatterplot(self, groupBy, selectBy, valSelect, wl_value,options_out):
        rrs_ins = np.array(self.mrfile.nc.variables['mu_ins_rrs'])
        rrs_sat = np.array(self.mrfile.nc.variables['mu_sat_rrs'])
        id_all = np.array(self.mrfile.nc.variables['mu_satellite_id'])
        mu_valid = np.array(self.mrfile.nc.variables[self.mu_valid_variable])
        # mu_valid = np.array(self.mrfile.nc.variables['mu_valid'])
        # print('MVALID: ',np.sum(mu_valid))

        valid_all = self.get_array_all_from_arraymu(id_all, mu_valid)
        valid_all = self.check_rrs_valid(valid_all, rrs_ins, rrs_sat)

        if wl_value is not None:
            wl_array = np.array(self.mrfile.nc.variables['mu_wavelength'])
            # wl_array_diff = np.abs(wl_array - wl_value)
            valid_all[wl_array != wl_value] = 0
            # valid_all[wl_array_diff > 1] = 0
        # wl_array = np.array(self.mrfile.nc.variables['mu_wavelength'])

        if selectBy is not None and valSelect is not None:
            if selectBy in options_out.keys() and 'flag_array' in options_out[selectBy]:
                select_array = options_out[selectBy]['flag_array']
            else:
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

    def check_rrs_valid(self, valid_array, rrs1, rrs2):
        dfv = netCDF4.default_fillvals['f4']
        for id in range(len(valid_array)):
            if valid_array[id] == 1:
                if rrs1[id] == dfv or rrs2[id] == dfv:
                    valid_array[id] = 0
        return valid_array

    def get_str_legend(self, options):
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

    # MAIN FUNCTION TO PLOT SCATTERPLOT
    def plot_scatter_plot(self, options, plot, index, wl):

        if options['include_stats'] or options['regression_line']:
            use_log_scale = options['log_scale']
            self.compute_statistics(use_log_scale)

        ngroup = 1
        str_legend = []
        groupValues = None
        if 'groupValues' in options.keys():
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
        if plot is None and index == -1:
            from PlotScatter import PlotScatter
            plot = PlotScatter()
            plot.close_plot()
            plot.start_plot()
        if plot is not None and index >= 0:
            plot.set_axhere_index(index)

        if options['scale_factor'] is not None:
            self.xdata = self.xdata * options['scale_factor']
            self.ydata = self.ydata * options['scale_factor']
            if len(self.yregress) > 0 and len(self.xregress) > 0:
                self.yregress = np.array(self.yregress) * options['scale_factor']
                self.xregress = np.array(self.xregress) * options['scale_factor']

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

        ngroupReal = 0

        if ngroup > 1:
            nmubygroup = [0] * ngroup
            for idx in range(ngroup):  # groupValues:
                g = groupValues[idx]
                if len(colors) == ngroup:
                    color = colors[idx]
                else:
                    if options['groupType'] == 'flag':

                        color = defaults.colors_default[idx]
                        # color = defaults.get_color_flag(g)
                    else:
                        color = defaults.get_color_ref(g)
                xhere = self.xdata[self.groupdata == g]
                yhere = self.ydata[self.groupdata == g]
                # None, None, color, 'gray', 1.5

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
        else:  # density scatter plot
            xhere = np.asarray(self.xdata, dtype=np.float)
            yhere = np.asarray(self.ydata, dtype=np.float)

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
                if len(self.xdata) > 0:
                    min_x_data = np.min(self.xdata)
                    min_y_data = np.min(self.ydata)
                    min_xy = np.floor(np.min([min_x_data, min_y_data]))
                    max_x_data = np.max(self.xdata)
                    max_y_data = np.max(self.ydata)
                    max_xy = np.ceil(np.max([max_x_data, max_y_data]))

            if options['min_xy'] is not None:
                min_xy = options['min_xy']
            if options['max_xy'] is not None:
                max_xy = options['max_xy']
        # print(min_xy, max_xy)
        # print(options['ticks'])

        if options['ticks'] is None:
            dif = max_xy - min_xy
            increm = 1
            if dif >= 8:
                increm = 2
            ticks = []
            for v in range(int(min_xy), int(max_xy) + 1, increm):
                if v <= max_xy:
                    ticks.append(v)
        else:
            ticks = options['ticks']

        plot.set_limits(min_xy, max_xy)
        if options['log_scale']:
            tlabels = []
            for t in ticks:
                tl = np.log10(t)
                if tl < 0:
                    tls = str(t)
                else:
                    tls = f'{t:.1f}'
                tlabels.append(tls)
            # print(ticks)
            # print(tlabels)
            plot.set_ticks_and_labels(ticks, tlabels, options['fontsizeaxis'])
        else:
            plot.set_ticks(ticks, options['fontsizeaxis'])

        if options['individual_axis'] or index == -1:
            plot.set_xaxis_title(options['xlabel'])
            plot.set_yaxis_title(options['ylabel'])
        else:
            if plot.index_row == plot.nrow - 1:
                plot.set_xaxis_title(options['xlabel'])
            if plot.index_col == 0:
                plot.set_yaxis_title(options['ylabel'])
            if plot.index_col > 0:
                plot.set_yticks_labels_off(ticks)
            if plot.index_row < (plot.nrow - 1):
                plot.set_xticks_labels_off(ticks)
        plot.set_equal_apect()

        if options['legend'] and len(str_legend) > 0 and index == -1:
            plot.set_legend(str_legend)
        if options['identity_line']:
            plot.plot_identity_line()

        if options['include_stats']:
            if options['type_scatterplot'] == 'rrs':
                if options['stat_list'] is not None:
                    stat_list = options['stat_list']
                    str0 = ''
                    for stat in stat_list:
                        if len(str0) > 0:
                            str0 = f'{str0}\n'
                        if stat == 'N':
                            val = self.valid_stats[stat]
                            str0 = f'{str0}N={val}'
                        if stat == 'NMATCH-UPS':
                            val = self.valid_stats['N']
                            if not options['selectByWavelength']:
                                nwl = np.unique(np.array(self.mrfile.variables['mu_wavelength'])).shape[0]
                                val = val / nwl
                            str0 = f'{str0}N={val:.0f}'

                        if stat == 'NMU':
                            if ngroupReal == 1:
                                val = self.valid_stats[stat]
                                str0 = f'{str0}N={val}'
                            else:
                                # iref = 1
                                valadded = []
                                for val in nmubygroup:
                                    if val in valadded:
                                        continue
                                    if val == 0:
                                        continue
                                    valadded.append(val)
                                if len(valadded) == 1:
                                    str0 = f'{str0}N={valadded[0]:.0f}'
                                else:
                                    str0 = f'{str0}N(1)={valadded[0]:.0f}'
                                    for idx in range(1, len(valadded)):
                                        iref = idx + 1
                                        str0 = f'{str0}\nN({iref})={val:.0f}'

                        if stat == 'r2':
                            val = self.valid_stats['DETER(r2)']
                            str0 = f'{str0}r\u00b2={val:.2f}'
                        if stat == 'RMSD' or stat == 'BIAS':
                            val = self.valid_stats[stat]
                            if stat == 'BIAS':
                                stat = stat.lower()
                            str0 = f'{str0}{stat}={val:.1e}'
                        if stat == 'RPD' or stat == 'APD':
                            val = self.valid_stats[stat]
                            str0 = f'{str0}{stat}={val:.0f}%'
                        if stat == 'WL':
                            wls = self.get_wl_str_from_wl(wl)
                            str0 = f'{str0}{wls} nm'
                            # str0 = f'{str0}{wl:.2f} nm'
                    xpos = options['stats_xpos']
                    ypos = options['stats_ypos']
                    plot.plot_text(xpos, ypos, str0)
                else:
                    if index == -1:
                        str0 = 'N={:d}\nRMSD={:,.1e} UNITS\nRPD={:,.0f}%\nAPD={:,.0f}%\n$r^2$={:,.2f}\nbias={:,.1e} UNITS' \
                            .format(self.valid_stats['N'],
                                    self.valid_stats['RMSD'],
                                    self.valid_stats['RPD'],
                                    self.valid_stats['APD'],
                                    self.valid_stats['DETER(r2)'],
                                    self.valid_stats['BIAS'])
                        str0 = str0.replace('UNITS', options['units'])
                        # print(str0)
                        plot.plot_text(0.05, 0.70, str0)
                    if index >= 0:
                        bias = self.valid_stats['BIAS']
                        r2 = self.valid_stats['DETER(r2)']
                        # rmsd = self.valid_stats['RMSD']
                        str0 = f'{wl:.2f} nm\nbias={bias:.1e}\nr\u00b2={r2:.2f}'
                        plot.plot_text(0.05, 0.70, str0)

            if options['type_scatterplot'] == 'chla' or options['type_scatterplot'] == 'kd':
                if options['stat_list'] is not None:
                    stat_list = options['stat_list']
                    str0 = ''
                    for stat in stat_list:
                        if len(str0) > 0:
                            str0 = f'{str0}\n'
                        if stat == 'N':
                            val = self.valid_stats[stat]
                            str0 = f'{str0}N={val}'
                        if stat == 'NMU':
                            val = self.valid_stats['N']
                            if ngroupReal > 1:
                                val = val / ngroupReal
                            str0 = f'{str0}N={val:.0f}'
                        if stat == 'r2':
                            val = self.valid_stats['DETER(r2)']
                            str0 = f'{str0}r\u00b2={val:.2f}'
                        if stat == 'RMSD' or stat == 'BIAS':
                            val = self.valid_stats[stat]
                            str0 = f'{str0}{stat}={val:.1e}'
                        if stat == 'RPD' or stat == 'APD':
                            val = self.valid_stats[stat]
                            str0 = f'{str0}{stat}={val:.0f}%'
                        if stat == 'WL':
                            str0 = f'{str0}{wl:.2f} nm'
                    xpos = options['stats_xpos']
                    ypos = options['stats_ypos']
                    plot.plot_text(xpos, ypos, str0)
                else:
                    str0 = 'N={:d}\nRPD={:,.0f}%\nAPD={:,.0f}%\n$r^2$={:,.2f}\nbias={:,.2f} UNITS' \
                        .format(self.valid_stats['N'],
                                self.valid_stats['RPD'],
                                self.valid_stats['APD'],
                                self.valid_stats['DETER(r2)'],
                                self.valid_stats['BIAS'])
                    # strunits = options['units']
                    # str0 = str0.replace('UNITS', f'log10({strunits})')
                    str0 = str0.replace('UNITS', '')
                    plot.plot_text(0.05, 0.74, str0)

        if not options['log_scale'] and options['regression_line']:
            if options['regression_line_groups']:
                if ngroup > 1:
                    for idx in range(len(groupValues)):
                        g = groupValues[idx]
                        if options['groupType'] == 'flag':
                            # color = defaults.get_color_flag(g)
                            color = defaults.colors_default[idx]
                        else:
                            color = defaults.get_color_ref(g)
                        xhere = self.xdata[self.groupdata == g]
                        yhere = self.ydata[self.groupdata == g]
                        if len(xhere) > 0 and len(yhere) > 0:
                            xregress, yregress = self.get_regression_line(xhere, yhere, None, None, min_xy, max_xy)
                            plot.plot_regress_line(xregress, yregress, color)
            else:
                plot.plot_regress_line(self.xregress, self.yregress, 'black')

        if options['title'] is not None:
            title_here = options['title']
            plot.set_title(title_here)
            plot.axhere.title.set_size(options['fontsizetitle'])

        if not options['file_out'] is None:
            plot.save_fig(options['file_out'])
            plot.close_plot()

        return plot

    def plot_statistics_bywl(self, options_out):
        flags = ['GLOBAL']
        flag_list = []
        flag_name = options_out['selectBy']
        if flag_name is not None and options_out['selectType'] == 'flag':
            flag_list = self.get_flag_list(options_out['selectValues'], options_out[flag_name]['flag_values'],
                                           options_out[flag_name]['flag_meanings'])
            flags = flag_list

        print(flag_list)

        wl_list = options_out['wl_values']
        params = list(options_out['params'])
        xdata_plot = [float(self.get_wl_str_from_wl(x)) for x in wl_list]
        wl_col = [self.get_wl_str_from_wl(x) for x in wl_list]

        table, indices = self.start_table_wl(flags, params, options_out['wl_values'])

        legend = []
        handles = []
        if len(flag_list) == 0 or 'GLOBAL' in flag_list:
            for wl in wl_list:
                self.set_data_scatterplot(None, None, None, wl,options_out)
                self.compute_statistics(False)
                table = self.assign_stats_table_wl(table, indices, params, 'GLOBAL', wl)
        if len(flag_list)>0:
            for idx in range(len(flag_list)):
                flag = flag_list[idx]
                flaga = flag
                if flag == 'FUB':
                    flaga = 'S3 FUB-CSIRO'
                if flag == 'STANDARD':
                    flaga = 'WFR'
                if flag == 'POLYMER':
                    flaga = 'CMEMS-OLCI'
                if flag == 'CCIALL':
                    flaga = 'OC-CCI v.6 (complete time series)'
                if flag == 'CCI':
                    flaga = 'OC-CCI v.6 (OLCI period)'

                legend.append(flaga)
                flag_value = options_out['selectValues'][idx]
                if flag_value>=0:
                    for wl in wl_list:
                        self.set_data_scatterplot(None, flag_name, flag_value, wl,options_out)
                        self.compute_statistics(False)
                        table = self.assign_stats_table_wl(table, indices, params, flag, wl)

        from PlotSpectra import PlotSpectra
        print(table)
        if options_out['multiple_plot'] is not None:
            from PlotMultiple import PlotMultiple
            rc = options_out['multiple_plot'].split(',')
            nrow = int(rc[0].strip())
            ncol = int(rc[1].strip())
            ntot = nrow * ncol
            files_multiple = ['']*ntot
            indices_files = list(range(ntot))
            if options_out['indices_files'] is not None:
                indices_files = options_out['indices_files']
            if options_out['multiple_files'] is not None:
                files_multiple = options_out['multiple_files']

        else:
            files_multiple = []
        multiple_ymin = options_out['multiple_ymin']
        multiple_ymax = options_out['multiple_ymax']

        # GLOBAL
        for iparam in range(len(params)):
            param = params[iparam]

            ymin = None
            if multiple_ymin is not None and len(multiple_ymin) == len(params):
                ymin = multiple_ymin[iparam]
            if ymin is None and options_out['ymin'] is not None:
                ymin = options_out['ymin']
            ymax = None
            if multiple_ymax is not None and len(multiple_ymax) == len(params):
                ymax = multiple_ymax[iparam]
            if ymax is None and options_out['ymax'] is not None:
                ymax = options_out['ymax']

            plot = PlotSpectra()
            plot.close_plot()
            plot.start_plot()
            if len(flag_list) == 0:  ##solo global
                irow = indices['GLOBAL'][param]
                ydata_plot = np.array(table[wl_col].iloc[irow])
                if param == 'RMSD' or param == 'BIAS':
                    ydata_plot = ydata_plot * 1000
                ydata_plot[ydata_plot > ymax] = ymax
                ydata_plot[ydata_plot < ymin] = ymin
                self.plot_spectra_line_impl(plot, xdata_plot, ydata_plot, 0, options_out)
                file_out = self.get_file_out_flag_param(options_out['file_out'], 'GLOBAL', param)
            else:
                for idx in range(len(flag_list)):
                    flag = flag_list[idx]
                    irow = indices[flag][param]
                    ydata_plot = np.array(table[wl_col].iloc[irow])
                    if param == 'RMSD' or param == 'BIAS':
                        ydata_plot = ydata_plot * 1000

                    if ymin is not None or ymax is not None:
                        options_out_oufr = {
                            'line_type': ':',
                            'line_color': 'k',
                            'line_width': options_out['line_width'],
                            'marker': options_out['marker'],
                            'marker_size': [0]
                        }

                        ydata_plot_oufr = []
                        xdata_plot_oufr = []
                        for ip in range(len(ydata_plot)):
                            if ymax >= ydata_plot[ip] >= ymin:
                                continue
                            if ip == 0:
                                xhere = np.array([xdata_plot[0], xdata_plot[1]])
                                yhere = np.array([ydata_plot[0], ydata_plot[1]])
                            elif ip == len(ydata_plot) - 1:
                                xhere = np.array([xdata_plot[ip - 1], xdata_plot[ip]])
                                yhere = np.array([ydata_plot[ip - 1], ydata_plot[ip]])
                            else:
                                xhere = np.array([xdata_plot[ip - 1], xdata_plot[ip], xdata_plot[ip + 1]])
                                yhere = np.array([ydata_plot[ip - 1], ydata_plot[ip], ydata_plot[ip + 1]])
                            ydata_plot_oufr.append(yhere)
                            xdata_plot_oufr.append(xhere)

                    h = self.plot_spectra_line_impl(plot, xdata_plot, ydata_plot, idx, options_out)

                    if ymin is not None or ymax is not None:
                        for ip in range(len(ydata_plot_oufr)):
                            xhere = xdata_plot_oufr[ip]
                            yhere = ydata_plot_oufr[ip]
                            # print(xhere, yhere)
                            # print(options_out_oufr)
                            self.plot_spectra_line_impl(plot, xhere, yhere, idx, options_out_oufr)

                    if param == params[0]:  # only for the first image
                        handles.append(h)
                file_out = self.get_file_out_flag_param(options_out['file_out'], flag_name, param)

            yminf, ymaxf = plot.get_y_range()
            if ymin is None:
                ymin = yminf
            if ymax is None:
                ymax = ymaxf
            plot.set_y_range(ymin, ymax)

            plot.set_xaxis_title(options_out['xlabel'])
            paramv = f'r$^2$'
            if param == 'RPD' or param == 'APD':
                paramv = f'{param}(%)'
            if param == 'RMSD' or param == 'BIAS':
                # paramv = f'{param}(sr$^-$$^1$)'
                paramv = f'{param} (10$^-$$^3$ sr$^-$$^1$)'
            plot.set_yaxis_title(paramv)
            plot.set_yticks(None, None, None, 12)
            xticks_size = 12
            if len(xdata_plot) == 16 or len(xdata_plot) == 15:
                xticks_size = 10

            plot.set_xticks(xdata_plot, wl_col, 90, xticks_size)
            if len(xdata_plot) == 16 or len(xdata_plot) == 15:
                xticks, xlabels = plot.get_xticks()
                xlabels[8].set_visible(False)
            plot.set_grid()

            plot.set_title(self.get_title(options_out['title'], None, None, param))
            if len(legend) > 0 and options_out['multiple_plot'] is None:
                plot.set_legend(legend)
            plot.set_tigth_layout()

            if file_out is not None:
                plot.save_fig(file_out)
                index_file = indices_files[iparam]
                #print(iparam,index_file,len(files_multiple))
                files_multiple[index_file] = file_out
                #files_multiple.append(file_out)
            plot.close_plot()

        if options_out['multiple_plot'] is not None:
            # from PlotMultiple import PlotMultiple
            # rc = options_out['multiple_plot'].split(',')
            # nrow = int(rc[0].strip())
            # ncol = int(rc[1].strip())
            # ntot = nrow * ncol
            if ntot == len(files_multiple):
                pm = PlotMultiple()
                xfigsize = options_out['xfigsize']
                yfigsize = options_out['yfigsize']
                wspace = options_out['widthspace']
                hspace = options_out['heightspace']
                pm.start_multiple_plot_advanced(nrow, ncol, xfigsize, yfigsize, wspace, hspace, True)
                index = 0
                for irow in range(nrow):
                    for icol in range(ncol):
                        pm.plot_image(files_multiple[index], irow, icol)
                        index = index + 1
                pm.set_text(-150, -50, '(a)')
                pm.set_text(1800, -50, '(b)')
                pm.set_text(-150, 1400, '(c)')
                pm.set_text(1800, 1400, '(d)')
                if len(legend) > 0 and len(legend) == len(handles):
                    pm.set_global_legend(handles, legend)
                file_out = self.get_file_out_flag_param(options_out['file_out'], None, None)

                pm.save_fig(file_out)

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

        h = plot.plot_single_line(ydata_plot, lc, lt, lw, m, ms)

        # ymin, ymax = plot.get_y_range()
        # if options_out['ymin'] is not None:
        #     ymin = options_out['ymin']
        # if options_out['ymax'] is not None:
        #     ymax = options_out['ymax']
        # plot.set_y_range(ymin,ymax)

        return h

    def plot_spectra_plot(self, options_out):

        type_rrs = options_out['type_rrs']

        if type_rrs.startswith('mu_'):  # mu_comparison, mu_sat, mu_ins
            self.plot_spectra_mu_comparison(options_out)
            return
        if type_rrs == 'flag_ins_comparison':
            self.plot_insitu_spectra_comparison_stats(options_out)
            return
        if type_rrs.startswith('flag_sat_insitu'):
            flag_ac = options_out['selectBy']
            flag_values = options_out['selectValues']
            if flag_ac is None:
                self.plot_comparison_stats_sat_insitu(options_out, None, None, None)
            else:

                self.plot_multiple_comparison_stats_sat_insitu(options_out, flag_ac, flag_values)
            return

        if type_rrs == 'flag_sat_comparison':
            flag_ac = options_out['selectBy']
            flag_values = options_out['selectValues']
            if flag_ac is None:
                self.plot_comparison_stats_sat(options_out, None, None, None)
            else:
                self.plot_multiple_comparison_stats_sat(options_out, flag_ac, flag_values)
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
        if options_out['plot_spectra'] is not None:
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
            pspectra.set_xticks(wavelength, wavelength, 90, 8)
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
        pspectra.set_y_range(-5, 10)
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

    def plot_multiple_comparison_stats_sat_insitu_deprecated(self, options_out, flag_name, flag_values):
        files_multiple = []
        handles = []
        # print(flag_name, flag_values)

        for flag_value in flag_values:
            flag = self.get_flag_flag(flag_value, options_out[flag_name]['flag_values'],
                                      options_out[flag_name]['flag_meanings'])
            if not options_out['file_out'] is None:
                file_out_base = options_out['file_out']
                file_out = file_out_base[:-4] + f'_{flag}.{self.format_image}'
                handles = self.plot_comparison_stats_sat_insitu(options_out, flag_name, flag_value, file_out)
                files_multiple.append(file_out)
                # print(file_out)

        if options_out['multiple_plot'] is not None:
            from PlotMultiple import PlotMultiple
            rc = options_out['multiple_plot'].split(',')
            nrow = int(rc[0].strip())
            ncol = int(rc[1].strip())
            ntot = nrow * ncol
            # print(nrow,ncol,ntot)
            if ntot == len(files_multiple):
                pm = PlotMultiple()
                xfigsize = options_out['xfigsize']
                yfigsize = options_out['yfigsize']
                wspace = options_out['widthspace']
                hspace = options_out['heightspace']
                pm.start_multiple_plot_advanced(nrow, ncol, xfigsize, yfigsize, wspace, hspace, False)
                index = 0
                for irow in range(nrow):
                    for icol in range(ncol):
                        pm.plot_image(files_multiple[index], irow, icol)
                        index = index + 1
                str_legend = ['insitu Rrs', 'satellite Rrs']
                if len(handles) == 2:
                    pm.set_global_legend(handles, str_legend)
                pm.save_fig(file_out_base)

    def plot_multiple_comparison_stats_sat_insitu(self, options_out, flag_name, flag_values):
        files_multiple = []
        handles = []
        multiple_ymin = options_out['multiple_ymin']
        multiple_ymax = options_out['multiple_ymax']
        # for flag_value in flag_values:
        for idx in range(len(flag_values)):
            flag_value = flag_values[idx]
            # print(flag_value,options_out[flag_name]['flag_values'],options_out[flag_name]['flag_meanings'])
            flag = self.get_flag_flag(flag_value, options_out[flag_name]['flag_values'],
                                      options_out[flag_name]['flag_meanings'])
            if not options_out['file_out'] is None:
                file_out_base = options_out['file_out']
                file_out = file_out_base[:-4] + f'_{flag}.{self.format_image}'
                if multiple_ymin is not None:
                    if len(flag_values) == len(multiple_ymin):
                        options_out['ymin'] = multiple_ymin[idx]
                if multiple_ymax is not None:
                    if len(flag_values) == len(multiple_ymax):
                        options_out['ymax'] = multiple_ymax[idx]
                handles = self.plot_comparison_stats_sat_insitu(options_out, flag_name, flag_value, file_out)
                files_multiple.append(file_out)
                # print(file_out)

        if options_out['multiple_plot'] is not None:
            from PlotMultiple import PlotMultiple
            rc = options_out['multiple_plot'].split(',')
            nrow = int(rc[0].strip())
            ncol = int(rc[1].strip())
            ntot = nrow * ncol
            # print(nrow,ncol,ntot)
            if ntot == len(files_multiple):
                pm = PlotMultiple()
                xfigsize = options_out['xfigsize']
                yfigsize = options_out['yfigsize']
                wspace = options_out['widthspace']
                hspace = options_out['heightspace']

                pm.start_multiple_plot_advanced(nrow, ncol, xfigsize, yfigsize, wspace, hspace, True)
                index = 0
                for irow in range(nrow):
                    for icol in range(ncol):
                        pm.plot_image(files_multiple[index], irow, icol)
                        index = index + 1
                str_legend = ['in situ Rrs', 'satellite Rrs']
                if len(handles) == 2:
                    pm.set_global_legend(handles, str_legend)

                anot = options_out['anot']
                if len(anot) > 0:
                    for an in anot:
                        pm.set_text(anot[an]['x'], anot[an]['y'], anot[an]['sval'])

                pm.save_fig(file_out_base)

    def plot_multiple_comparison_stats_sat(self, options_out, flag_name, flag_values):
        files_multiple = []
        h_legend = []
        str_legend = []
        multiple_ymin = options_out['multiple_ymin']
        multiple_ymax = options_out['multiple_ymax']

        for idx in range(len(flag_values)):
            flag_value = flag_values[idx]
            flag = self.get_flag_flag(flag_value, options_out[flag_name]['flag_values'],
                                      options_out[flag_name]['flag_meanings'])
            if not options_out['file_out'] is None:
                file_out_base = options_out['file_out']
                file_out = file_out_base[:-4] + f'_{flag}.{self.format_image}'
                if multiple_ymin is not None:
                    if len(flag_values) == len(multiple_ymin):
                        options_out['ymin'] = multiple_ymin[idx]
                if multiple_ymax is not None:
                    if len(flag_values) == len(multiple_ymax):
                        options_out['ymax'] = multiple_ymax[idx]
                str_legend, h_legend = self.plot_comparison_stats_sat(options_out, flag_name, flag_value, file_out)
                files_multiple.append(file_out)
                # print(file_out)

        if options_out['multiple_plot'] is not None:
            from PlotMultiple import PlotMultiple
            rc = options_out['multiple_plot'].split(',')
            nrow = int(rc[0].strip())
            ncol = int(rc[1].strip())
            ntot = nrow * ncol
            # print(nrow,ncol,ntot)
            if ntot == len(files_multiple):
                pm = PlotMultiple()
                xfigsize = options_out['xfigsize']
                yfigsize = options_out['yfigsize']
                wspace = options_out['widthspace']
                hspace = options_out['heightspace']

                pm.start_multiple_plot_advanced(nrow, ncol, xfigsize, yfigsize, wspace, hspace, True)
                index = 0
                for irow in range(nrow):
                    for icol in range(ncol):
                        pm.plot_image(files_multiple[index], irow, icol)
                        index = index + 1
                #str_legend = ['in situ Rrs', 'satellite Rrs']
                if len(str_legend) > 0:
                    pm.set_global_legend(h_legend, str_legend)

                anot = options_out['anot']
                if len(anot) > 0:
                    for an in anot:
                        pm.set_text(anot[an]['x'], anot[an]['y'], anot[an]['sval'])

                pm.save_fig(file_out_base)

    def plot_comparison_stats_sat_insitu(self, options_out, flag_name, flag_value, file_out):
        from PlotSpectra import PlotSpectra
        wlvalues = options_out['wl_values']
        allwlvalues = list(np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])))
        if wlvalues is None:
            wlvalues = allwlvalues

        self.mrfile.var_mu_valid = self.mu_valid_variable
        # print(f'[INFO] MU VALID VARIABLE: {self.mrfile.var_mu_valid}')
        flag_array = None
        if flag_name is not None and flag_name in options_out.keys() and 'flag_array' in options_out[flag_name].keys():
            flag_array = options_out[flag_name]['flag_array']
        sat_stats, insitu_stats = self.mrfile.get_all_spectra_insitu_sat_with_wlvalues(options_out['scale_factor'],
                                                                                       wlvalues, flag_name, flag_value,
                                                                                       flag_array)

        pspectra = PlotSpectra()
        if options_out['stat_plot_method'] is not None:
            if options_out['stat_plot_method'] == 'iqr':
                pspectra.set_iqr_as_stats_plot()
            if options_out['stat_plot_method'].startswith('std'):
                try:
                    factor = float(options_out['stat_plot_method'].split(',')[1])
                except:
                    factor = 1.0
                pspectra.set_std_as_stats_plot(factor)

        imin, imax = pspectra.get_imin_imax_from_wavelength(wlvalues, options_out['wl_min'], options_out['wl_max'])
        str_legend = ['insitu Rrs', 'satellite Rrs']

        if options_out['ymin'] is None or options_out['ymax'] is None:
            ymin_sat, ymax_sat = pspectra.get_ymin_ymax_from_stats(sat_stats, imin, imax)
            ymin_insitu, ymax_insitu = pspectra.get_ymin_ymax_from_stats(insitu_stats, imin, imax)
            ymin = np.min([ymin_sat, ymin_insitu])
            ymax = np.max([ymax_sat, ymax_insitu])
        if options_out['ymin'] is not None:
            ymin = options_out['ymin']
        if options_out['ymax'] is not None:
            ymax = options_out['ymax']

        wl_list = wlvalues[imin:imax]
        xdata_plot = [float(self.get_wl_str_from_wl(x)) for x in wl_list]
        wl_col = [self.get_wl_str_from_wl(x) for x in wl_list]

        pspectra.xdata = xdata_plot

        color = 'red'
        pspectra.stats_style['central']['color'] = color
        pspectra.stats_style['central']['marker'] = 'o'
        pspectra.stats_style['central']['markersize'] = 5
        pspectra.stats_style['fill']['color'] = color
        pspectra.stats_style['fill']['framealpha'] = 0.5
        hlineinsitu = pspectra.plot_stats(insitu_stats, imin, imax)

        color = 'blue'
        pspectra.xdata = wlvalues[imin:imax]
        pspectra.stats_style['central']['color'] = color
        pspectra.stats_style['central']['marker'] = 'o'
        pspectra.stats_style['central']['markersize'] = 5
        pspectra.stats_style['fill']['color'] = color
        pspectra.stats_style['fill']['framealpha'] = 0.5
        hlinesat = pspectra.plot_stats(sat_stats, imin, imax)

        h_legend = [hlineinsitu, hlinesat]

        xticks_size = 12

        if len(xdata_plot) == 16 or len(xdata_plot) == 15:
            xticks_size = 10
        pspectra.set_xticks(xdata_plot, wl_col, 90, xticks_size)
        if len(xdata_plot) == 16 or len(xdata_plot) == 15:
            xticks, xlabels = pspectra.get_xticks()
            xlabels[8].set_visible(False)
        pspectra.set_yticks(None, None, None, 12)
        pspectra.set_y_range(ymin, ymax)
        pspectra.set_xaxis_title(options_out['xlabel'])
        pspectra.set_yaxis_title(options_out['ylabel'])
        if options_out['title'] is not None:
            title_here = options_out['title']
            flag = self.get_flag_flag(flag_value, options_out[flag_name]['flag_values'],
                                      options_out[flag_name]['flag_meanings'])
            title_here = title_here.replace('$FLAG$', flag)
            if title_here == 'FUB':
                title_here = 'S3 FUB-CSIRO'
            if title_here == 'STANDARD':
                title_here = 'WFR'
            if title_here == 'POLYMER':
                title_here = 'CMEMS-OLCI'
            if title_here == 'CCIALL':
                title_here = 'OC-CCI v.6 (complete time series)'
            if title_here == 'CCI':
                title_here = 'OC-CCI v.6 (OLCI period)'
            pspectra.set_title(title_here)
        pspectra.set_grid()
        pspectra.legend_options['bbox_to_anchor'] = (0.65, 1.0)
        pspectra.legend_options['framealpha'] = 1
        # pspectra.set_legend_h(h_legend, str_legend)
        pspectra.set_tigth_layout()
        if file_out is None:
            file_out = options_out['file_out']
        if not file_out is None:
            pspectra.save_fig(file_out)
        pspectra.close_plot()

        return h_legend

    def plot_comparison_stats_sat(self, options_out, flag_name, flag_value, file_out):

        from PlotSpectra import PlotSpectra
        wlvalues = options_out['wl_values']
        allwlvalues = list(np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])))
        if wlvalues is None:
            wlvalues = allwlvalues

        self.mrfile.var_mu_valid = self.mu_valid_variable

        flag_array = None
        if flag_name is not None and flag_name in options_out.keys() and 'flag_array' in options_out[flag_name].keys():
            flag_array = options_out[flag_name]['flag_array']


        if options_out['groupBy'] is not None:
            group_array = np.array(self.mrfile.variables[options_out['groupBy']])
            group_values = options_out['groupValues']
        else:
            group_array = np.ones((self.mrfile.n_mu_total))
            group_values = [1]

        str_legend = []
        all_sat_stats = []
        insitu_stats = None
        for gvalue in group_values:
            sat_stats, insitu_stats_here = self.mrfile.get_all_spectra_insitu_sat_with_wlvalues_group(options_out['scale_factor'],
                                                                                       wlvalues, flag_name, flag_value,
                                                                                       flag_array,gvalue,group_array)
            if insitu_stats is None:
                insitu_stats = insitu_stats_here
            all_sat_stats.append(sat_stats)
            if options_out['groupBy'] is not None:
                all_meanings = options_out[options_out['groupBy']]['flag_meanings']
                all_values =  options_out[options_out['groupBy']]['flag_values'].tolist()
                str_legend.append(all_meanings[all_values.index(gvalue)])


        pspectra = PlotSpectra()
        if options_out['stat_plot_method'] is not None:
            if options_out['stat_plot_method'] == 'iqr':
                pspectra.set_iqr_as_stats_plot()
            if options_out['stat_plot_method'].startswith('std'):
                try:
                    factor = float(options_out['stat_plot_method'].split(',')[1])
                except:
                    factor = 1.0
                pspectra.set_std_as_stats_plot(factor)

        imin, imax = pspectra.get_imin_imax_from_wavelength(wlvalues, options_out['wl_min'], options_out['wl_max'])

        if options_out['ymin'] is None or options_out['ymax'] is None:
            ymin_sat, ymax_sat = pspectra.get_ymin_ymax_from_stats(sat_stats, imin, imax)
            ymin_insitu, ymax_insitu = pspectra.get_ymin_ymax_from_stats(insitu_stats, imin, imax)
            ymin = np.min([ymin_sat, ymin_insitu])
            ymax = np.max([ymax_sat, ymax_insitu])
        if options_out['ymin'] is not None:
            ymin = options_out['ymin']
        if options_out['ymax'] is not None:
            ymax = options_out['ymax']

        wl_list = wlvalues[imin:imax]
        xdata_plot = [float(self.get_wl_str_from_wl(x)) for x in wl_list]
        wl_col = [self.get_wl_str_from_wl(x) for x in wl_list]

        pspectra.xdata = xdata_plot
        nseries = len(group_values)
        colors = options_out['series_color']
        if colors is None:
            colors = defaults.get_color_list(nseries)

        h_legend = []
        for idx in range(nseries):
            sat_stats = all_sat_stats[idx]
            color = colors[idx]
            pspectra.stats_style['central']['color'] = color
            pspectra.stats_style['central']['marker'] = 'o'
            pspectra.stats_style['central']['markersize'] = 7
            pspectra.stats_style['central']['linewidth'] = 2
            pspectra.stats_style['dispersion']['color'] = color
            pspectra.stats_style['dispersion']['marker'] = 'o'
            pspectra.stats_style['dispersion']['markersize'] = 0
            pspectra.stats_style['dispersion']['linestyle'] = '--'
            pspectra.stats_style['dispersion']['linewidth'] = 1
            pspectra.stats_style['fill']['color'] = None
            pspectra.stats_style['fill']['framealpha'] = 0.5
            hline = pspectra.plot_stats(sat_stats, imin, imax)
            h_legend.append(hline)

        str_legend.append('in situ Rrs')
        pspectra.stats_style['central']['color'] = 'black'
        pspectra.stats_style['central']['marker'] = 'o'
        pspectra.stats_style['central']['markersize'] = 5
        pspectra.stats_style['central']['linewidth'] = 1
        pspectra.stats_style['dispersion']['markersize'] = 0
        pspectra.stats_style['dispersion']['linewidth'] = 0
        # pspectra.stats_style['dispersion']['color'] = 'black'
        # pspectra.stats_style['dispersion']['markersize'] = 0
        # pspectra.stats_style['dispersion']['linestyle'] = '--'
        # pspectra.stats_style['dispersion']['linewidth'] = 0
        pspectra.stats_style['fill']['color'] = 'black'
        pspectra.stats_style['fill']['framealpha'] = 0.5
        hline = pspectra.plot_stats(insitu_stats, imin, imax)
        h_legend.append(hline)
        # color = 'red'
        # pspectra.stats_style['central']['color'] = color
        # pspectra.stats_style['central']['marker'] = 'o'
        # pspectra.stats_style['central']['markersize'] = 5
        # pspectra.stats_style['fill']['color'] = color
        # pspectra.stats_style['fill']['framealpha'] = 0.5
        # hlineinsitu = pspectra.plot_stats(insitu_stats, imin, imax)
        #
        # color = 'blue'
        # pspectra.xdata = wlvalues[imin:imax]
        # pspectra.stats_style['central']['color'] = color
        # pspectra.stats_style['central']['marker'] = 'o'
        # pspectra.stats_style['central']['markersize'] = 5
        # pspectra.stats_style['fill']['color'] = color
        # pspectra.stats_style['fill']['framealpha'] = 0.5
        # hlinesat = pspectra.plot_stats(sat_stats, imin, imax)
        #
        # h_legend = [hlineinsitu, hlinesat]
        #
        xticks_size = 12
        #
        if len(xdata_plot) == 16 or len(xdata_plot) == 15:
            xticks_size = 10
        pspectra.set_xticks(xdata_plot, wl_col, 90, xticks_size)
        if len(xdata_plot) == 16 or len(xdata_plot) == 15:
            xticks, xlabels = pspectra.get_xticks()
            xlabels[8].set_visible(False)
        pspectra.set_yticks(None, None, None, 12)
        pspectra.set_y_range(ymin, ymax)
        pspectra.set_xaxis_title(options_out['xlabel'])
        pspectra.set_yaxis_title(options_out['ylabel'])
        if options_out['title'] is not None:
            title_here = options_out['title']
            flag = self.get_flag_flag(flag_value, options_out[flag_name]['flag_values'],
                                      options_out[flag_name]['flag_meanings'])
            title_here = title_here.replace('$FLAG$', flag)
            if title_here == 'FUB':
                title_here = 'S3 FUB-CSIRO'
            if title_here == 'STANDARD':
                title_here = 'WFR'
            if title_here == 'POLYMER':
                title_here = 'CMEMS-OLCI'
            if title_here == 'CCIALL':
                title_here = 'OC-CCI v.6 (complete time series)'
            if title_here == 'CCI':
                title_here = 'OC-CCI v.6 (OLCI period)'
            pspectra.set_title(title_here)
        pspectra.set_grid()
        pspectra.legend_options['bbox_to_anchor'] = (0.65, 1.0)
        pspectra.legend_options['framealpha'] = 1
        #pspectra.set_legend_h(h_legend, str_legend)
        pspectra.set_tigth_layout()
        if file_out is None:
            file_out = options_out['file_out']
        if not file_out is None:
            pspectra.save_fig(file_out)
        pspectra.close_plot()

        return str_legend,h_legend

    def plot_insitu_spectra_comparison_stats(self, options_out):
        # print(options_out)
        from PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        type_rrs = options_out['type_rrs']
        if type_rrs != 'flag_ins_comparison':
            return
        flag_values = options_out['groupValues']
        flag_name = options_out['groupBy']
        str_legend = self.get_flag_list(flag_values, options_out[flag_name]['flag_values'],
                                        options_out[flag_name]['flag_meanings'])
        h_legend = []

        wavelength = self.mrfile.get_insitu_wl()
        imin, imax = pspectra.get_imin_imax_from_wavelength(wavelength, options_out['wl_min'],options_out['wl_max'])

        ymin = None
        ymax = None
        for flag_value in flag_values:
            wavelength = self.mrfile.get_insitu_wl()
            spectra_good_flag, stats_flag = self.mrfile.get_flag_insitu_valid_spectra(options_out['scale_factor'],flag_name, flag_value)
            ymin_flag, ymax_flag = pspectra.get_ymin_ymax_from_stats(stats_flag, imin, imax)
            if ymin is None and ymax is None:
                ymin = ymin_flag
                ymax = ymax_flag
            else:
                if ymin_flag < ymin:
                    ymin = ymin_flag
                if ymax_flag > ymax:
                    ymax = ymax_flag

            color = defaults.get_color_flag(flag_value)
            pspectra.xdata = wavelength[imin:imax]
            pspectra.stats_style['central']['color'] = color
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

    def plot_flag_fromoptions(self, options_out):
        flag_options = options_out['flag_options']
        # for option in flag_options:
        #     print(option, flag_options[option])
        self.mrfile.window_size = options_out['window_size']
        limit_to_valid_ac = options_out['limit_to_valid_ac']
        valid_ac = None
        if limit_to_valid_ac is not None:
            valid_ac = self.mrfile.obtain_mu_valid_with_ac(options_out['limit_to_valid_ac'])


        flag_info = self.mrfile.analyze_sat_flags_advanced(flag_options, valid_ac)

        only_bigger_than_zero = options_out['only_bigger_than_zero']
        nseries = options_out['nseries']
        plot_param = options_out['plot_param']
        legend = options_out['series_names']

        # print(only_bigger_than_zero,nseries,plot_param)
        # print(options_out)
        file_out = options_out['file_out']
        xlabel = 'Total number of match-ups'
        if nseries == 0:
            self.plot_flag(flag_info, plot_param, only_bigger_than_zero, xlabel, file_out)
        else:
            self.plot_flag_series(flag_info, plot_param, xlabel, nseries, legend, file_out, options_out)

        # for finfo in flag_info:
        #     print(finfo, flag_info[finfo]['nmu'],flag_info[finfo]['nmacrow'])

    def plot_flag_series(self, flag_info, plot_param, xlabel, nseries, legend, file_out, options_out):
        xticks = []
        ylabel = plot_param
        data_figure = {}
        for flag in flag_info:
            name = flag_info[flag]['seriesoutput']
            if name not in xticks:
                xticks.append(name)
                data_figure[name] = {
                    'data': [0] * nseries,
                    'anot': [0] * nseries,
                    'nseriesplot': 0
                }
            for iserie in flag_info[flag]['seriesid']:
                ndata = flag_info[flag][ylabel]

                andata = flag_info[flag]['pmacrow']
                #andata = flag_info[flag]['ntotalw'] /(17*17*flag_info[flag]['nmu'])
                if ndata > 0:
                    data_here = data_figure[name]['data']
                    data_here[iserie - 1] = ndata
                    anot_here = data_figure[name]['anot']
                    anot_here[iserie - 1] = andata
                    data_figure[name]['data'] = data_here
                    data_figure[name]['anot'] = anot_here
                    data_figure[name]['nseriesplot'] = data_figure[name]['nseriesplot'] + 1

        only_bigger_than_zero = options_out['only_bigger_than_zero']

        colors = options_out['series_color']
        if colors is None:
            colors = defaults.get_color_list(nseries)

        from matplotlib import pyplot as plt
        hfig, ax = plt.subplots()
        height = 0.2
        xval = 0
        xticks_pos = []
        xticks_minor_pos = []
        handles = [None] * nseries
        heightbyseriesbase = height / nseries

        for name in data_figure:
            xini = xval
            nseriesplot = data_figure[name]['nseriesplot']
            if nseriesplot == 0:
                print(f'[WARNING] No flags for: {name}')
            data_here = data_figure[name]['data']
            if only_bigger_than_zero:
                heightbyseries = height / nseriesplot
            else:
                heightbyseries = heightbyseriesbase
            xticks_minor_pos.append(xini - (heightbyseries / 2))
            for idx in range(nseries - 1, -1, -1):
                do_plot = True
                if data_here[idx] == 0 and only_bigger_than_zero:
                    do_plot = False
                if do_plot:
                    hbar = plt.barh(xval, data_here[idx], height=heightbyseries, color=colors[idx])
                    handles[idx] = hbar
                    xval = xval + heightbyseries
            xfin = xval - heightbyseries
            xticks_minor_pos.append(xfin + (heightbyseries / 2))
            xoutput = (xini + xfin) / 2
            xticks_pos.append(xoutput)

        plt.xlabel(xlabel, fontsize=12)
        plt.yticks(xticks_pos, xticks)
        ax.set_yticks(xticks_minor_pos, minor=True)
        plt.grid(b=True, which='minor', color='gray', linestyle='--', axis='y')
        plt.grid(b=True, which='major', color='gray', linestyle='--', axis='x')
        ax.tick_params(which='major', length=0, axis='y')
        ax.tick_params(which='minor', length=10, axis='y')

        ##anotation
        xlim = ax.get_xlim()
        xrange = xlim[1] - xlim[0]
        xmin_new = (-1) * (xrange / 8)
        if xmin_new>-18:
            xmin_new = -18
        ax.set_xlim(left=xmin_new)
        xval = 0
        for name in data_figure:
            anot_here = data_figure[name]['anot']
            for idx in range(nseries - 1, -1, -1):
                val = anot_here[idx]
                vals = f'{val:.2f} %'
                plt.text(xmin_new + 1, xval - (heightbyseries / 2) + 0.005, vals, fontdict={'fontsize': 10})
                # plt.barh(xval,xmin_new, height=heightbyseries, color=colors[idx])
                xval = xval + heightbyseries

        if legend is not None:
            # plt.legend(legend)
            plt.legend(handles, legend, framealpha=1, loc='lower right')
        plt.tight_layout()

        if file_out is not None:
            plt.savefig(file_out, dpi=300)
        plt.close()

    def plot_flag(self, info_flag, plot_param, only_gtz, xlabel, file_out):
        x = list(info_flag.keys())
        ylabel = plot_param
        ydata = []
        for flag in info_flag:
            ydata.append(info_flag[flag][ylabel])
        x = np.array(x)
        ydata = np.array(ydata)
        if only_gtz:
            x = x[ydata > 0]
            ydata = ydata[ydata > 0]

        from matplotlib import pyplot as plt
        # print(x)
        # print(ydata)
        plt.figure()
        plt.barh(x, ydata)
        plt.grid(b=True, which='major', color='gray', linestyle='--')
        plt.xlabel(xlabel, fontsize=12)
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

    def plot_multiple_scatter_plots(self, options_out):
        file_out_base = options_out['file_out']
        title_base = options_out['title']
        rc = options_out['multiple_plot'].split(',')
        nrow = int(rc[0].strip())
        ncol = int(rc[1].strip())
        ntot = nrow * ncol
        index = 0
        from PlotScatter import PlotScatter
        plot_here = PlotScatter()
        plot_here.nrow = nrow
        plot_here.ncol = ncol
        plot_here.index_row = 0
        plot_here.index_col = 0
        plot_here.xtitle_options['fontsize'] = options_out['fontsizelabels']
        plot_here.ytitle_options['fontsize'] = options_out['fontsizelabels']
        plot_here.plot_text_options['fontsize'] = options_out['fontsizestats']

        plot_here.start_multiple_plot_advanced(nrow, ncol, options_out['xfigsize'],
                                               options_out['yfigsize'], options_out['widthspace'],
                                               options_out['heightspace'])
        print(f'[INFO] Starting multiple plot with {nrow} rows and {ncol} cols')
        wl_values = options_out['wl_values']
        print(f'[INFO] Wavelenthts: {wl_values}')
        for wl in options_out['wl_values']:
            selectBy = None
            selectValue = None
            if options_out['selectBy'] is not None:
                selectBy = options_out['selectBy']
                selectValue = options_out['selectValues']
            self.set_data_scatterplot(options_out['groupBy'], selectBy, selectValue, wl,options_out)
            if len(self.xdata) > 0 and len(self.ydata) > 0:
                options_out['title'] = self.get_title(title_base, wl, None, None)
                options_out['file_out'] = None

                self.plot_scatter_plot(options_out, plot_here, index, wl)
                plot_here.index_col = plot_here.index_col + 1
                if plot_here.index_col == plot_here.ncol:
                    plot_here.index_col = 0
                    plot_here.index_row = plot_here.index_row + 1
                index = index + 1
            else:
                print(f'[WARNING] No data for wavelength: {wl} nm')

        for index_blank in range(index, ntot):
            plot_here.plot_blanck(index_blank)
        if options_out['legend']:
            str_legend = self.get_str_legend(options_out)
            if len(str_legend) > 0:
                plot_here.set_global_legend(str_legend)
        file_out = self.get_file_out_name(file_out_base, None, None)
        plot_here.save_fig(file_out)
        plot_here.close_plot()

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
        self.set_data_scatterplot(None, None, None, None,options_out)
        self.compute_statistics(False)
        table = self.assign_stats_table_wl(table, indices, params, 'GLOBAL', None)
        for wl in options_out['wl_values']:
            self.set_data_scatterplot(None, None, None, wl,options_out)
            self.compute_statistics(False)
            table = self.assign_stats_table_wl(table, indices, params, 'GLOBAL', wl)
        # results by flag
        if flag_name is not None:
            for idx in range(len(flag_list)):
                flag = flag_list[idx]
                flag_value = options_out['selectValues'][idx]
                self.set_data_scatterplot(None, flag_name, flag_value, None,options_out)
                self.compute_statistics(False)
                table = self.assign_stats_table_wl(table, indices, params, flag, None)
                ng = len(flag_list) * len(options_out['wl_values'])
                for wl in options_out['wl_values']:
                    self.set_data_scatterplot(None, flag_name, flag_value, wl,options_out)
                    self.compute_statistics(False)
                    table = self.assign_stats_table_wl(table, indices, params, flag, wl)

        if options_out['formatted']:
            for index, row in table.iterrows():
                for col in table.columns:
                    if col == 'FLAG' or col == 'PARAM':
                        continue
                    val = table.at[index, col]
                    if row['PARAM'] == 'RMSD' or row['PARAM'] == 'BIAS':
                        table.at[index, col] = f'{val:.2e}'
                    if row['PARAM'] == 'APD' or row['PARAM'] == 'RPD':
                        table.at[index, col] = f'{val:.0f}'
                    if row['PARAM'] == 'DETER(r2)':
                        table.at[index, col] = f'{val:.2f}'

        # print(table)
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
        self.set_data_scatterplot(None, None, None, wl,options_out)
        use_log_scale = False
        if 'log_scale' in options_out:
            use_log_scale = options_out['log_scale']
        self.compute_statistics(use_log_scale)
        table = self.assign_table(table, indices, params, 'GLOBAL')
        # stats by flag
        if len(flag_list) > 0:
            for idx in range(len(flag_list)):
                flag = flag_list[idx]
                flag_value = options_out['selectValues'][idx]
                self.set_data_scatterplot(None, flag_name, flag_value, wl,options_out)
                self.compute_statistics(use_log_scale)
                table = self.assign_table(table, indices, params, flag)

        if not options_out['file_out'] is None:
            file_out = options_out['file_out']
            if wl is not None:
                wls = self.get_wl_str_from_wl(wl)
                wls = wls.replace('.', '_')
                file_out = file_out[:-4] + '_' + wls + '.csv'
            table.to_csv(file_out, sep=';')

    def get_options(self, options, section):
        options_out = {'apply': self.get_value_param(options, section, 'apply', False, 'boolean')}
        if not options_out['apply']:
            return options_out
        options_out['type'] = self.get_value_param(options, section, 'type', None, 'str')
        if options_out['type'] is None:
            return options_out
        options_out['name'] = section
        options_out['multiple_plot'] = self.get_value_param(options, section, 'multiple_plot', None, 'str')

        options_out = self.get_anot_options(options, section, options_out)

        if options_out['type'] == 'csvtable':
            options_out = self.get_options_csv(options, section, options_out)
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
            options_out = self.get_select_options(options, section, options_out)
        if options_out['type'] == 'flagplot':
            options_out = self.get_options_flag(options, section, options_out)
        if options_out['type'] == 'temporalheatmap':
            options_out = self.get_options_temporal_heatmap(options, section, options_out)
        if options_out['type'] == 'validityplot':
            options_out = self.get_options_validity_plot(options, section, options_out)
        if options_out['type'] == 'multipleplot':
            options_out = self.get_options_multiple_plot(options, section, options_out)
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
                    # options_out['groupValues'] = self.get_value_param(options, section, 'groupValues', flag_values,'intlist')
                    flag_list_config = self.get_value_param(options, section, 'groupValues', None, 'strlist')
                    if flag_list_config is None:
                        # print(flag_values)
                        options_out['groupValues'] = flag_values
                    else:
                        flag_values_config = []
                        for flag_config in flag_list_config:
                            if flag_config.strip()=='GLOBAL':
                                flag_values_config.append(-1)
                                continue
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
            options_out['selectType'] = 'float'
            var_name = options_out['selectBy']
            virtual_flag = False
            if not var_name in self.mrfile.variables:
                print(
                    f'[WARNING] {var_name} is not a variable defined in the MDB file. Checking {var_name} as virtual flag')
                virtual_flag = True
                flag_values, flag_meanings, flag_array = self.get_virtual_flag(options, var_name)
                options_out[var_name] = {
                    'flag_values': flag_values,
                    'flag_meanings': flag_meanings,
                    'flag_array': flag_array
                }
                options_out['selectType'] = 'flag'
                print(f'[INFO] Set virtual flag: {var_name}')

            if var_name.startswith('flag') and not virtual_flag:
                flag_values = self.mrfile.nc.variables[var_name].flag_values
                flag_meanings_list = self.mrfile.nc.variables[var_name].flag_meanings.split(' ')
                flag_meanings = [x.strip() for x in flag_meanings_list]
                options_out[var_name] = {
                    'flag_values': flag_values,
                    'flag_meanings': flag_meanings
                }
                options_out['selectType'] = 'flag'

            if options_out['selectType'] == 'flag':
                flag_list_config = self.get_value_param(options, section, 'selectValues', None, 'strlist')
                if flag_list_config is None:
                    options_out['selectValues'] = flag_values
                else:
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
                    options_out['selectValues'] = flag_values_config
            if options_out['selectType'] == 'float':
                group_values = list(np.unique(np.array(self.mrfile.nc.variables[var_name])))
                options_out['selectValues'] = self.get_value_param(options, section, 'selectValues', group_values,
                                                                   'floatlist')

        return options_out

    def get_anot_options(self, options, section, options_out):
        anot = {}
        idx = 0
        while idx >= 0:
            key = f'anot_{idx}'
            ostr = self.get_value_param(options, section, key, None, 'str')
            if ostr is None:
                break
            olist = ostr.split(',')
            try:
                xval = int(olist[0].strip())
                yval = int(olist[1].strip())
                sval = ','.join(olist[2:])
                anot[str(idx)] = {
                    'x': xval,
                    'y': yval,
                    'sval': sval
                }
                idx = idx + 1
            except:
                idx = idx + 1
                pass

        options_out['anot'] = anot
        return options_out

    def get_options_scatterplot(self, options, section, options_out):
        options_out['type_scatterplot'] = self.get_value_param(options, section, 'type_scatterplot', 'rrs', 'str')
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
        options_out['stats_xpos'] = self.get_value_param(options, section, 'stats_xpos', 0.05, 'float')
        options_out['stats_ypos'] = self.get_value_param(options, section, 'stats_ypos', 0.70, 'float')
        options_out['individual_axis'] = self.get_value_param(options, section, 'individual_axis', False, 'boolean')
        return options_out

    def get_options_statstable(self, options, section, options_out):
        options_out['selectByWavelength'] = True  ##option to be always true
        if options_out['wl_values'] is None:
            options_out['wl_values'] = list(np.unique(np.array(self.mrfile.nc.variables['mu_wavelength'])))

        options_out['params'] = self.get_value_param(options, section, 'params', self.valid_stats.keys(), 'strlist')
        options_out['formatted'] = self.get_value_param(options, section, 'params', False, 'boolean')

        if self.output_path is not None:
            ext = '.csv'
            if options_out['formatted']:
                ext = '_formatted.csv'
            name_default = options_out['name'] + ext
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

        options_out['xfigsize'] = self.get_value_param(options, section, 'xfigsize', 7, 'float')
        options_out['yfigsize'] = self.get_value_param(options, section, 'yfigsize', 7, 'float')
        options_out['widthspace'] = self.get_value_param(options, section, 'widthspace', 0.1, 'float')
        options_out['heightspace'] = self.get_value_param(options, section, 'heightspace', 0.1, 'float')

        options_out['widthspace'] = self.get_value_param(options, section, 'widthspace', 0.1, 'float')
        options_out['heightspace'] = self.get_value_param(options, section, 'heightspace', 0.1, 'float')

        options_out['ymin'] = self.get_value_param(options, section, 'ymin', None, 'float')
        options_out['ymax'] = self.get_value_param(options, section, 'ymax', None, 'float')
        options_out['multiple_ymin'] = self.get_value_param(options, section, 'multiple_ymin', None, 'floatlist')
        options_out['multiple_ymax'] = self.get_value_param(options, section, 'multiple_ymax', None, 'floatlist')

        list_files = self.get_value_param(options, section, 'multiple_files', None, 'strlist')
        multiple_files = None
        if list_files is not None:
            multiple_files = []
            for l in list_files:
                if os.path.exists(l):
                    multiple_files.append(l)
                else:
                    f = os.path.join(self.output_path, l)
                    if os.path.exists(f):
                        multiple_files.append(f)
                    else:
                        multiple_files.append('')
        options_out['multiple_files'] = multiple_files

        options_out['indices_files'] = self.get_value_param(options,section,'indices_files',None,'intlist')
        return options_out

    def get_options_spectraplot(self, options, section, options_out):
        options_out['type_rrs'] = self.get_value_param(options, section, 'type_rrs', 'ins', 'str')
        options_out['wl_min'] = self.get_value_param(options, section, 'wl_min', None, 'float')
        options_out['wl_max'] = self.get_value_param(options, section, 'wl_max', None, 'float')
        options_out['ymin'] = self.get_value_param(options, section, 'ymin', None, 'float')
        options_out['ymax'] = self.get_value_param(options, section, 'ymax', None, 'float')
        # options_out['xticks'] = self.get_value_param(options, section, 'xticks', None, 'floatlist')
        # options_out['yticks'] = self.get_value_param(options, section, 'yticks', None, 'floatlist')
        options_out['multiple_ymin'] = self.get_value_param(options, section, 'multiple_ymin', None, 'floatlist')
        options_out['multiple_ymax'] = self.get_value_param(options, section, 'multiple_ymax', None, 'floatlist')
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
        options_out['xfigsize'] = self.get_value_param(options, section, 'xfigsize', 7, 'float')
        options_out['yfigsize'] = self.get_value_param(options, section, 'yfigsize', 7, 'float')
        options_out['widthspace'] = self.get_value_param(options, section, 'widthspace', 0.1, 'float')
        options_out['heightspace'] = self.get_value_param(options, section, 'heightspace', 0.1, 'float')
        options_out['stat_plot_method'] = self.get_value_param(options, section, 'stat_plot_method', None, 'str')
        options_out['series_color'] = self.get_value_param(options, section, 'series_color', None, 'strlist')
        return options_out

    def get_options_flag(self, options, section, options_out):

        options_out['window_size'] = self.get_value_param(options, section, 'window_size', 3, 'int')
        options_out['limit_to_valid_ac'] = self.get_value_param(options, section, 'limit_to_valid_ac', None, 'str')
        options_out['nseries'] = self.get_value_param(options, section, 'nseries', 0, 'int')
        series_names = self.get_value_param(options, section, 'series_names', None, 'strlist')
        if options_out['nseries'] > 0:
            if series_names is None:
                series_names = []
            else:
                if len(series_names) != options_out['nseries']:
                    series_names = []
            if len(series_names) == 0:
                for idx in range(options_out['nseries']):
                    iserie = idx + 1
                    series_names.append(f'Series_{iserie}')
        options_out['series_names'] = series_names
        options_out['series_color'] = self.get_value_param(options, section, 'series_color', None, 'strlist')
        options_out['only_bigger_than_zero'] = self.get_value_param(options, section, 'only_bigger_than_zero', False,
                                                                    'boolean')
        options_out['plot_param'] = self.get_value_param(options, section, 'plot_param', 'nmacrow', 'str')

        flag_options = {}
        index_option = 0
        has_options = True
        while has_options:
            key_base = f'flag_option_{index_option}'
            flag_var_name = self.get_value_param(options, section, f'{key_base}.flag_var_name', None, 'str')
            if flag_var_name is None:
                break
            flag_list = self.get_value_param(options, section, f'{key_base}.flag_list', None, 'strlist')
            flag_ref = self.get_value_param(options, section, f'{key_base}.flag_ref', None, 'str')
            if flag_list is None and flag_ref is None:
                print(f'[ERROR] flag_list or flag_ref must be defined')
                return flag_options
            var_group_name = self.get_value_param(options, section, f'{key_base}.var_group_name', None, 'str')
            var_group_flag = self.get_value_param(options, section, f'{key_base}.var_group_flag', None, 'str')
            var_group_value = None
            if var_group_name is not None and var_group_flag is not None:
                var_group_value = self.mrfile.get_flag_value(var_group_name, var_group_flag)
            seriesid = self.get_value_param(options, section, f'{key_base}.seriesid', 0, 'intlist')
            seriesoutput = self.get_value_param(options, section, f'{key_base}.seriesoutput', None, 'str')

            # flag_ref = self.get_value_param(options, section, f'{key_base}.flag_ref', None, 'str')
            ##flag_list falg_method: eq
            if flag_list is not None:
                flags_output = self.get_value_param(options, section, f'{key_base}.flag_output', None, 'strlist')
                if flags_output is None:
                    flags_output = flag_list
                if len(flags_output) != len(flag_list):
                    print('[WARNING] Flags output list is diferent from flag list. Setting flags_output==flag_list')
                    flags_output = flag_list
                for idx in range(len(flag_list)):
                    flag_output = flags_output[idx]
                    seriesoutput = flag_output
                    if var_group_name is not None and var_group_value is None and len(seriesid) == len(series_names):
                        flag_output_base = flag_output
                        for iserie in seriesid:
                            serie_name = series_names[iserie - 1]
                            flag_output = f'{flag_output_base}_{serie_name}'
                            serie_value = self.mrfile.get_flag_value(var_group_name, serie_name)
                            # print(iserie,serie_name,flag_output)
                            flag_options[flag_output] = {
                                'flag_var_name': flag_var_name,
                                'flag_ref': flag_list[idx],
                                'flag_method': 'eq',
                                'var_group_name': var_group_name,
                                'var_group_value': serie_value,
                                'seriesid': [iserie],
                                'seriesoutput': seriesoutput
                            }
                    else:
                        flag_options[flag_output] = {
                            'flag_var_name': flag_var_name,
                            'flag_ref': flag_list[idx],
                            ' flag_method': 'eq',
                            'var_group_name': var_group_name,
                            'var_group_value': var_group_value,
                            'seriesid': seriesid,
                            'seriesoutput': seriesoutput
                        }
            if flag_ref is not None:
                flag_output = self.get_value_param(options, section, f'{key_base}.flag_output', None, 'str')
                if flag_output is None:
                    flag_output = flag_ref
                if seriesoutput is None:
                    seriesoutput = flag_output

                if var_group_name is not None and var_group_value is None and len(seriesid) == len(series_names):
                    flag_output_base = flag_output
                    for iserie in seriesid:
                        serie_name = series_names[iserie - 1]
                        flag_output = f'{flag_output_base}_{serie_name}'
                        serie_value = self.mrfile.get_flag_value(var_group_name, serie_name)
                        # print(iserie,serie_name,flag_output)
                        flag_options[flag_output] = {
                            'flag_var_name': flag_var_name,
                            'flag_ref': flag_ref,
                            'flag_method': 'list',
                            'var_group_name': var_group_name,
                            'var_group_value': serie_value,
                            'seriesid': [iserie],
                            'seriesoutput': seriesoutput
                        }
                else:
                    flag_options[flag_output] = {
                        'flag_var_name': flag_var_name,
                        'flag_ref': flag_ref,
                        'flag_method': 'list',
                        'var_group_name': var_group_name,
                        'var_group_value': var_group_value,
                        'seriesid': seriesid,
                        'seriesoutput': seriesoutput
                    }

            index_option = index_option + 1

        options_out['flag_options'] = flag_options
        if self.output_path is not None:
            name_default = options_out['name'] + '.' + self.format_image
            file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(options, section, 'file_out', file_out_default, 'str')

        return options_out

    def get_options_csv(self, options, section, options_out):
        if self.output_path is not None:
            name_default = options_out['name'] + '.csv'
            file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(options, section, 'file_out', file_out_default, 'str')

        return options_out

    def get_options_temporal_heatmap(self, options, section, options_out):
        if self.output_path is not None:
            name_default = options_out['name'] + '.' + self.format_image
            file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(options, section, 'file_out', file_out_default, 'str')
        flag_value = self.get_value_param(options, section, 'flag', None, 'str')

        if flag_value is not None and flag_value.lower() == 'none':
            flag_value = None
        options_out['flag'] = flag_value
        options_out['output_type'] = self.get_value_param(options, section, 'output_type', None, 'strlist')
        options_out['flag_list'] = self.get_value_param(options, section, 'flag_list', None, 'strlist')

        options_out['vmin'] = self.get_value_param(options, section, 'vmin', None, 'float')
        options_out['vmax'] = self.get_value_param(options, section, 'vmax', None, 'float')

        options_out['xfigsize'] = self.get_value_param(options, section, 'xfigsize', 7, 'float')
        options_out['yfigsize'] = self.get_value_param(options, section, 'yfigsize', 7, 'float')
        options_out['widthspace'] = self.get_value_param(options, section, 'widthspace', 0.1, 'float')
        options_out['heightspace'] = self.get_value_param(options, section, 'heightspace', 0.1, 'float')

        return options_out

    def get_options_validity_plot(self, options, section, options_out):
        if self.output_path is not None:
            name_default = options_out['name'] + '.' + self.format_image
            file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(options, section, 'file_out', file_out_default, 'str')
        flag_value = self.get_value_param(options, section, 'flag', None, 'str')

        if flag_value is not None and flag_value.lower() == 'none':
            flag_value = None
        options_out['flag'] = flag_value
        options_out['output_type'] = self.get_value_param(options, section, 'output_type', None, 'strlist')
        options_out['flag_list'] = self.get_value_param(options, section, 'flag_list', None, 'strlist')
        options_out['series_color'] = self.get_value_param(options, section, 'series_color', None, 'strlist')
        options_out['series_flag'] = self.get_value_param(options, section, 'series_flag', None, 'strlist')
        options_out['file_csv'] = self.get_value_param(options, section, 'file_csv', None, 'str')

        # options_out['xfigsize'] = self.get_value_param(options, section, 'xfigsize', 7, 'float')
        # options_out['yfigsize'] = self.get_value_param(options, section, 'yfigsize', 7, 'float')
        # options_out['widthspace'] = self.get_value_param(options, section, 'widthspace', 0.1, 'float')
        # options_out['heightspace'] = self.get_value_param(options, section, 'heightspace', 0.1, 'float')

        return options_out

    def get_options_multiple_plot(self, options, section, options_out):
        options_out['xfigsize'] = self.get_value_param(options, section, 'xfigsize', 7, 'float')
        options_out['yfigsize'] = self.get_value_param(options, section, 'yfigsize', 7, 'float')
        options_out['widthspace'] = self.get_value_param(options, section, 'widthspace', 0.1, 'float')
        options_out['heightspace'] = self.get_value_param(options, section, 'heightspace', 0.1, 'float')
        if self.output_path is not None:
            name_default = options_out['name'] + '.' + self.format_image
            file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(options, section, 'file_out', file_out_default, 'str')

        list_files = self.get_value_param(options, section, 'multiple_files', None, 'strlist')
        multiple_files = []
        if list_files is not None:
            for l in list_files:
                if os.path.exists(l):
                    multiple_files.append(l)
                else:
                    f = os.path.join(self.output_path, l)
                    if os.path.exists(f):
                        multiple_files.append(f)
        options_out['multiple_files'] = multiple_files

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
        # if type.startswith('multiple_strlist_'):
        #     idx = int(type.split('_')[2])
        #     list_str_general = value.split(';')
        #     if idx < len(list_str_general):
        #         value_n = list_str_general[idx]
        #         list_str = value_n.split(',')
        #         list = []
        #         for vals in list_str:
        #             list.append(vals.strip())
        #         return list
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

    def get_flag_list(self, values, allValues, allFlags):
        flag_list = []
        for val in values:
            if val==-1:
                flag_list.append('GLOBAL')
            indext = np.where(allValues == val)
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

    def get_virtual_flag(self, options, var_name):
        flag_values = None
        flag_meanings = None
        flag_array = None
        type_f = self.get_value_param(options, var_name, 'type', None, 'str')
        type_v = self.get_value_param(options, var_name, 'typevirtual', 'flags_c', 'str')
        if type_f != 'virtual_flag':
            return flag_values, flag_meanings, flag_array

        if type_v == 'flags_c':
            flag_names = self.get_value_param(options, var_name, 'flag_names', None, 'strlist')
            if flag_names is None:
                return flag_values, flag_meanings, flag_array

            if len(flag_names) == 2:
                flag_meanings_list_1 = self.mrfile.nc.variables[flag_names[0]].flag_meanings.split(' ')
                flag_meanings_list_1 = [x.strip() for x in flag_meanings_list_1]
                flag_meanings_list_2 = self.mrfile.nc.variables[flag_names[1]].flag_meanings.split(' ')
                flag_meanings_list_2 = [x.strip() for x in flag_meanings_list_2]
                flag_values_1 = self.mrfile.nc.variables[flag_names[0]].flag_values.tolist()
                flag_values_2 = self.mrfile.nc.variables[flag_names[1]].flag_values.tolist()
                flag_array = np.zeros(self.mrfile.nc.variables[flag_names[0]].shape, dtype=np.int32)
                flag_meanings = []
                flag_values = []
                flag_vm_dict = {}
                idx = 0
                for m1 in flag_meanings_list_1:
                    for m2 in flag_meanings_list_2:
                        meaning = f'{m1}_{m2}'
                        value = int(math.pow(2, idx))
                        flag_meanings.append(meaning)
                        flag_values.append(value)
                        flag_vm_dict[meaning] = value
                        idx = idx + 1

                for idx in range(flag_array.shape[0]):
                    v1 = self.mrfile.nc.variables[flag_names[0]][idx]
                    v2 = self.mrfile.nc.variables[flag_names[1]][idx]
                    m1 = flag_meanings_list_1[flag_values_1.index(v1)]
                    m2 = flag_meanings_list_2[flag_values_2.index(v2)]
                    m = f'{m1}_{m2}'
                    val = flag_vm_dict[m]
                    flag_array[idx] = val

                flag_values = np.array(flag_values)

                return flag_values, flag_meanings, flag_array

    def compute_statistics(self, use_log_scale):

        self.valid_stats['N'] = len(self.xdata)
        if self.valid_stats['N'] == 0:
            for key in self.valid_stats:
                self.valid_stats[key] = np.nan
            return

        # self.valid_stats['NMU'] = self.valid_stats['N']/ng

        # Generated linear fit
        xdatal = []
        ydatal = []
        maxxy = 0
        minxy = None
        for x, y in zip(self.xdata, self.ydata):
            if use_log_scale:
                xdatal.append(math.log10(x))
                ydatal.append(math.log10(y))
            else:
                if np.isnan(x) or np.isnan(y):
                    print(x, y)
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

        slope, intercept, r_value, p_value, std_err = stats.linregress(xdatal, ydatal)

        self.xregress, self.yregress = self.get_regression_line(xdatal, ydatal, slope, intercept, minxy, maxxy)
        # self.xregress = []
        # self.yregress = []
        # self.xregress.append(0)
        # self.yregress.append(intercept)
        # for x in xdatal:
        #     yr = (x * slope) + intercept
        #     self.yregress.append(yr)
        #     self.xregress.append(x)
        # yrmax = ((maxxy + 1) * slope) + intercept
        # self.xregress.append(maxxy + 1)
        # self.yregress.append(yrmax)

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
        self.valid_stats['CPRMSE'] = cprmse

        bias = np.mean(sat_obs - ref_obs)
        self.valid_stats['BIAS'] = bias

        mae = np.mean(np.abs(sat_obs - ref_obs))
        self.valid_stats['MAE'] = mae

        self.valid_stats['DETER(r2)'] = r_value * r_value
