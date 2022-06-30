import math

import pandas as pd

from MDBFile import MDBFile
import MDBPlotDefaults as defaults
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
from scipy import stats
import COMMON.common_functions as cfs
import os
from PlotSpectra import PlotSpectra
from PlotScatter import PlotScatter
from datetime import datetime as dt
import seaborn as sns
from scipy.stats import gaussian_kde
from pylr2 import regress2
from sklearn.metrics import r2_score


class MDBPlot:
    def __init__(self, mdb_file, df_param):
        self.mfile = mdb_file
        self.dfparam = df_param

        # plot general options
        self.title = ''
        self.file_name_base = ''
        self.format_image = 'jpg'

        # validation scatterplot
        self.xdata = []
        self.ydata = []
        self.wldata = []
        self.xregress = []
        self.yregress = []
        self.xlabel = defaults.xlabel_default
        self.ylabel = defaults.ylabel_default

        # spectra plot
        # self.ydata = []
        # self.xdata = []

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

    def get_df_val(self):
        dfval = None
        if self.mfile is not None and isinstance(self.mfile, MDBFile):
            dfval = self.mfile.df_validation_valid
            if len(dfval.index) == 0:
                self.mfile.prepare_df_validation()
                dfval = self.mfile.df_validation_valid
        if self.mfile is None and self.dfparam is not None:
            dfval = self.dfparam
        return dfval

    def plot_scatter_plot(self, title, legend, include_stats, file_out):
        if include_stats:
            self.compute_statistics()
        wl = self.wldata.unique()
        nwl = len(wl)
        str_legend = []
        if nwl > 0 and legend:
            for w in wl:
                str_legend.append(f'{w:.2f}')

        plot = PlotScatter()
        plot.close_plot()
        plot.start_plot()

        self.xdata = self.xdata * 1000
        self.ydata = self.ydata * 1000
        self.yregress = np.array(self.yregress) * 1000
        self.xregress = np.array(self.xregress) * 1000

        if nwl > 1:
            for w in wl:
                color = defaults.get_color_ref(w)
                xhere = self.xdata[self.wldata == w]
                yhere = self.ydata[self.wldata == w]
                plot.plot_data(xhere, yhere, None, None, color, 'gray', 1.5)
        else:  # density scatter plot
            xhere = np.asarray(self.xdata,dtype=np.float)
            yhere = np.asarray(self.ydata,dtype=np.float)
            xy = np.vstack([xhere, yhere])
            z = gaussian_kde(xy)(xy)
            idx = z.argsort()
            xhere, yhere, z = xhere[idx], yhere[idx], z[idx]
            plot.plot_data(xhere, yhere, None, 25, z, None, None)
            plot.set_cmap('jet')

        plot.set_equal_apect()
        max_x_data = np.min(self.xdata)
        max_y_data = np.max(self.ydata)
        max_xy = np.ceil(np.max([max_x_data, max_y_data]))
        # if max_xy > 10:
        #     max_xy = 10
        # if max_xy <= 4:
        #     max_xy = 5
        plot.set_limits(0, max_xy)
        plot.set_xaxis_title(self.xlabel)
        plot.set_yaxis_title(self.ylabel)
        if legend:
            plot.set_legend(str_legend)
        plot.plot_identity_line()

        if include_stats:
            str0 = 'N={:d}\nRMSD={:,.4f}\nAPD={:,.0f}%\nRPD={:,.0f}%\n$r^2$={:,.2f}\nbias={:,.4f}' \
                .format(self.valid_stats['N'],
                        self.valid_stats['rmse_val'],
                        self.valid_stats['mean_abs_rel_diff'],
                        self.valid_stats['mean_rel_diff'],
                        self.valid_stats['r_value'] ** 2,
                        self.valid_stats['bias'])
            # if nwl == 1:
            #     w = wl[0]
            #     strwl = f'λ = {w:0.2f} nm \n'
            #     str0 = strwl + str0
            plot.plot_text(0.05, 0.70, str0)

            plot.plot_regress_line(self.xregress, self.yregress, 'black')

        # data_plot = pd.concat([self.xdata, self.ydata], axis=1).astype(dtype=np.float)

        # sns.lmplot(data = data_plot,x='Ins_Rrs',y='Sat_Rrs',line_kws={'color': [0.3,0.3,0.3,1]})

        if title:
            title_here = self.title
            if nwl == 1:
                w = wl[0]
                strwl = f' λ = {w:0.2f} nm \n'
                title_here = self.title + strwl
            plot.set_title(title_here)

        if not file_out is None:
            plot.save_fig(file_out)
            plot.close_plot()

    def plot_scatter_plot_prev(self, title, legend, include_stats, file_out):
        wl = self.wldata.unique()
        nwl = len(wl)
        str_legend = []
        if nwl > 0 and legend:
            for w in wl:
                str_legend.append(f'{w:.2f}')

        plt.close()
        plt.figure()
        for xpoint, ypoint, wlpoint in zip(self.xdata, self.ydata, self.wldata):
            if nwl == 1:
                color_point = defaults.get_color_ref(wlpoint)
                plt.scatter(xpoint, ypoint, c=color_point, edgecolors='gray', linewidths=1.0)
            else:
                color_point = defaults.get_color_ref(wlpoint)
                # c=defaults.color_dict[f'{wlpoint:.2f}]'
                plt.scatter(xpoint, ypoint, c=color_point,
                            edgecolors='gray', linewidths=1.5)

        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel(self.xlabel, fontsize=12)
        plt.ylabel(self.ylabel, fontsize=12)
        if legend:
            plt.legend(str_legend, loc='upper left', bbox_to_anchor=(1.0, 1.0))
        xmin, xmax = plt.gca().get_xlim()
        ymin, ymax = plt.gca().get_ylim()
        xmin = np.min([xmin, ymin])
        xmax = np.max([xmax, ymax])
        plt.plot([xmin, xmax], [xmin, xmax], '--k')

        if include_stats:
            self.compute_statistics()
            str0 = 'N={:d}\nRMSD={:,.4f}\nAPD={:,.0f}%\nRPD={:,.0f}%\n$r^2$={:,.2f}\nbias={:,.4f}' \
                .format(self.valid_stats['N'],
                        self.valid_stats['rmse_val'],
                        self.valid_stats['mean_abs_rel_diff'],
                        self.valid_stats['mean_rel_diff'],
                        self.valid_stats['r_value'] ** 2,
                        self.valid_stats['bias'])
            if nwl == 1:
                w = wl[0]
                strwl = f'λ = {w:0.2f} nm \n'
                str0 = strwl + str0
            plt.text(0.05, 0.70, str0, horizontalalignment='left', fontsize=12, transform=plt.gca().transAxes)

        if title:
            plt.title(self.title)

        if not file_out is None:
            plt.savefig(file_out, dpi=300)

    def plot_spectra_plot(self, title, legend, file_out):
        plt.close()
        plt.figure()
        if not legend is None:
            df = pd.DataFrame(np.transpose(self.ydata), columns=legend, index=self.xdata)
            df.plot(lw=2, marker='.', markersize=10, xticks=self.xdata)
        else:
            df = pd.DataFrame(np.transpose(self.ydata), index=self.xdata)

            df.plot(lw=1, color='black', marker='.', markersize=10, legend=False, xticks=range(len(self.xdata)),
                    mec='gray')

        plt.xlabel(self.xlabel, fontsize=12)
        plt.ylabel(self.ylabel, fontsize=12)
        plt.xticks(rotation=90)
        plt.grid(b=True, which='major', color='gray', linestyle='--')
        if title:
            plt.title(self.title)
        plt.gcf().tight_layout()
        if not file_out is None:
            plt.savefig(file_out, dpi=300)

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
        self.valid_stats['r_value'] = r_value
        self.valid_stats['p_value'] = p_value
        self.valid_stats['std_err'] = std_err

        ref_obs = np.asarray(self.xdata,dtype=np.float)
        sat_obs = np.asarray(self.ydata,dtype=np.float)

        results = regress2(ref_obs, sat_obs, _method_type_2="reduced major axis")
        self.valid_stats['slope_typeII'] = results['slope']
        self.valid_stats['offset_typeII'] = results['intercept']

        self.valid_stats['rmse_val'] = cfs.rmse(sat_obs, ref_obs)

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
        rel_diff = 100 * (ref_obs - sat_obs) / ref_obs
        self.valid_stats['mean_rel_diff'] = np.mean(rel_diff)

        #  the mean of absolute (unsigned) percent differences
        self.valid_stats['mean_abs_rel_diff'] = np.mean(np.abs(rel_diff))

        bias = np.mean(sat_obs - ref_obs)
        self.valid_stats['bias'] = bias

        mae = np.mean(np.abs(sat_obs - ref_obs))
        self.valid_stats['MAE'] = mae

        self.valid_stats['r2'] = r_value * r_value



    def plot_all_scatter_plot(self, path_out):
        dfval = self.get_df_val()
        if dfval is None:
            return
        # ['Index','Index_MU','Index_Band','Sat_Time','Ins_Time','Time_Diff','Wavelenght','Ins_Rrs','Sat_Rrw','Valid']

        self.xdata = dfval[dfval['Valid']]['Ins_Rrs']
        self.ydata = dfval[dfval['Valid']]['Sat_Rrs']
        self.wldata = dfval[dfval['Valid']]['Wavelenght']

        show_title = False
        file_out = None
        if not path_out is None:
            show_title = True
            if self.mfile is not None:
                self.title = self.mfile.get_title()
            file_out = os.path.join(path_out, f'{self.get_file_name(None)}.{self.format_image}')

        self.plot_scatter_plot(show_title, True, True, file_out)

    def plot_all_scatter_plot_insitu(self, path_out):
        dfval = self.get_df_val()
        if dfval is None:
            return
        # ['Index','Index_MU','Index_Band','Sat_Time','Ins_Time','Time_Diff','Wavelenght','Ins_Rrs','Sat_Rrw','Valid']

        self.xdata = dfval[:]['PanthyrRRS']
        self.ydata = dfval[:]['HypstarRRS']
        self.wldata = dfval[:]['Wavelength']

        show_title = False
        file_out = None
        if not path_out is None:
            show_title = True
            if self.mfile is not None:
                self.title = self.mfile.get_title()
            file_out = os.path.join(path_out, f'{self.get_file_name(None)}.{self.format_image}')

        self.plot_scatter_plot(show_title, False, True, file_out)

    def plot_wavelength_scatter_plot(self, index_sat, wl, path_out):
        if wl is None and 0 <= index_sat < len(self.mfile.satellite_bands):
            wl = self.mfile.satellite_bands[index_sat]
        if wl is None:
            return
        wl = float(wl)
        wls = f'{wl:0.2f}'
        wls = wls.replace('.', '_')
        dfval = self.get_df_val()
        if dfval is None:
            return

        self.xdata = dfval[(dfval['Valid']) & (dfval['Wavelenght'] == wl)]['Ins_Rrs']
        self.ydata = dfval[(dfval['Valid']) & (dfval['Wavelenght'] == wl)]['Sat_Rrs']
        self.wldata = dfval[(dfval['Valid']) & (dfval['Wavelenght'] == wl)]['Wavelenght']
        # self.xdata = dfval[(dfval['Wavelength'] == wl)]['PanthyrRRS']
        # self.ydata = dfval[(dfval['Wavelength'] == wl)]['HypstarRRS']
        # self.wldata = dfval[(dfval['Wavelength'] == wl)]['Wavelength']

        print('NValues: ', len(self.xdata), 'Wavelength: ', wl)
        show_title = False
        file_out = None
        if not path_out is None:
            show_title = True
            if self.mfile is not None:
                self.title = self.mfile.get_title()
            file_out = os.path.join(path_out, f'{self.get_file_name(wls)}.{self.format_image}')

        self.plot_scatter_plot(show_title, False, True, file_out)

    def plot_wavelenght_scatter_plots(self, path_out, wllist):
        if wllist is None and self.mfile is not None:
            wllist = self.mfile.wlref
        for wl in wllist:
            self.plot_wavelength_scatter_plot(-1, wl, path_out)

    def compute_all_statistics(self, wllist):
        if wllist is None and self.mfile is not None:
            wllist = self.mfile.wlref
        col_names = ['Param', 'All']
        # for w in self.mfile.satellite_bands:
        for w in wllist:
            w_str = f'{w:0.2f}'
            col_names.append(w_str)
        self.df_valid_stats = pd.DataFrame(columns=col_names)
        dfval = self.get_df_val()
        if dfval is None:
            return

        self.xdata = dfval[dfval['Valid']]['Ins_Rrs']
        self.ydata = dfval[dfval['Valid']]['Sat_Rrs']
        self.wldata = dfval[dfval['Valid']]['Wavelenght']

        # self.xdata = dfval[:]['PanthyrRRS']
        # self.ydata = dfval[:]['HypstarRRS']
        # self.wldata = dfval[:]['Wavelength']

        self.compute_statistics()
        for param in self.valid_stats:
            row = {}
            for col in col_names:
                if col == 'Param':
                    row[col] = [param]
                elif col == 'All':
                    row[col] = [self.valid_stats[param]]
                else:
                    row[col] = 0
            # self.df_valid_stats = self.df_valid_stats.append(row, ignore_index=True)
            self.df_valid_stats = pd.concat([self.df_valid_stats, pd.DataFrame.from_dict(row)], ignore_index=True)

        nparam = len(self.df_valid_stats)

        for wl in wllist:
            print(f'Computing statistics for: {wl}')
            w_str = f'{wl:0.2f}'
            self.xdata = dfval[(dfval['Valid']) & (dfval['Wavelenght'] == wl)]['Ins_Rrs']
            self.ydata = dfval[(dfval['Valid']) & (dfval['Wavelenght'] == wl)]['Sat_Rrs']
            self.wldata = dfval[(dfval['Valid']) & (dfval['Wavelenght'] == wl)]['Wavelenght']

            # self.xdata = dfval[(dfval['Wavelength'] == wl)]['PanthyrRRS']
            # self.ydata = dfval[(dfval['Wavelength'] == wl)]['HypstarRRS']
            # self.wldata = dfval[(dfval['Wavelength'] == wl)]['Wavelength']

            self.compute_statistics()
            for iparam in range(nparam):
                param = self.df_valid_stats.at[iparam, 'Param']
                self.df_valid_stats.at[iparam, w_str] = self.valid_stats[param]

    def plot_spectra_param(self, param, path_out):
        if self.df_valid_stats is None:
            self.compute_all_statistics()

        file_out = None

        params = list(self.df_valid_stats.loc[:, 'Param'])

        if param in params or (param == 'r2' and 'r_value' in params):
            self.xdata = self.df_valid_stats.columns[2:len(self.df_valid_stats.columns)]
            if param == 'r2':
                self.ydata = self.df_valid_stats.loc[self.df_valid_stats['Param'] == 'r_value', self.xdata]
                self.ydata = self.ydata ** 2
                self.ylabel = f'r\N{SUPERSCRIPT TWO}'
            else:
                self.ydata = self.df_valid_stats.loc[self.df_valid_stats['Param'] == param, self.xdata]
                self.ylabel = param
            self.xlabel = defaults.xlabel_wl_default

            if self.mfile is not None:
                self.title = self.mfile.get_title()
            if not path_out is None:
                file_name = f'{self.get_file_name_param(param)}.{self.format_image}'
                file_out = os.path.join(path_out, file_name)
            self.plot_spectra_plot(True, None, file_out)

    def plot_all_spectra_param(self, path_out):
        if self.df_valid_stats is None:
            self.compute_all_statistics()
        params = list(self.df_valid_stats.loc[:, 'Param'])
        for param in params:
            if param == 'N' or param == 'p_value':
                continue
            if param == 'r_value':
                param = 'r2'
            self.plot_spectra_param(param, path_out)

    def plot_all_insitu_spectra(self, path_out, wlmin, wlmax):

        self.mfile.qc_insitu.compute_band_statistics(wlmin, wlmax)
        spectra, wldata = self.mfile.qc_insitu.get_all_valid_spectra(wlmin, wlmax)
        median_array = self.mfile.qc_insitu.get_stat_spectra(wlmin, wlmax, 'median')
        p25_array = self.mfile.qc_insitu.get_stat_spectra(wlmin, wlmax, 'p25')
        p75_array = self.mfile.qc_insitu.get_stat_spectra(wlmin, wlmax, 'p75')

        ps = PlotSpectra()
        ps.xdata = wldata
        ps.set_xaxis_title(defaults.xlabel_wl_default)
        ps.set_yaxis_title(defaults.label_insitu_default)
        nspectra = spectra.shape[0]
        for ispectra in range(nspectra):
            ps.plot_single_line(spectra[ispectra], 'gray', None, 1)
        ps.plot_single_line(median_array, 'black', None, 2)
        ps.plot_single_line(p25_array, 'black', '--', 1)
        ps.plot_single_line(p75_array, 'black', '--', 1)

        file_out = None
        if not path_out is None:
            file_out = os.path.join(path_out, f'{self.file_name_base}.{self.format_image}')
        if not file_out is None:
            ps.save_plot(file_out)

    def obtain_mu_info(self, dfall, dfvalid, path_out_base):
        if dfvalid is None:
            dfvalid = self.dfparam

        path_out = os.path.join(path_out_base, 'MUINFO')
        if not os.path.exists(path_out):
            os.mkdir(path_out)

        file_jpg = os.path.join(path_out, 'AllMUByMonth.jpg')
        fcsv = os.path.join(path_out, 'NMuAllMonth.csv')
        self.obtain_mu_info_impl(dfall, fcsv, file_jpg)

        file_jpg = os.path.join(path_out, 'ValidMUByMonth.jpg')
        fcsv = os.path.join(path_out, 'NMuValidMonth.csv')
        self.obtain_mu_info_impl(dfvalid, fcsv, file_jpg)

    def obtain_mu_info_impl(self, dfall, fcsv, fjpg):
        nall = len(dfall.index)
        first_date = dt.strptime(dfall.iloc[0].at['Sat_Time'], '%Y-%m-%d %H:%M')
        last_date = dt.strptime(dfall.iloc[nall - 1].at['Sat_Time'], '%Y-%m-%d %H:%M')
        year_min = first_date.year
        year_max = last_date.year + 1
        year = list(range(year_min, year_max))
        year.reverse()
        month = list(range(1, 13))
        dfall_month = pd.DataFrame(index=year, columns=month, dtype=np.float)
        dfall_month[:] = 0
        for index, row in dfall.iterrows():
            dif_wl = np.abs(np.float(row['Wavelenght']) - 412)
            if dif_wl < 5:
                date_here = dt.strptime(row['Sat_Time'], '%Y-%m-%d %H:%M')
                year_here = date_here.year
                month_here = date_here.month
                dfall_month.at[year_here, month_here] = dfall_month.at[year_here, month_here] + 1

        h = plt.Figure()
        cmap = cm.get_cmap('RdYlBu_r')
        dfall_month_withnan = dfall_month
        dfall_month_withnan[dfall_month == 0] = np.nan
        sns.heatmap(dfall_month_withnan, annot=False, cmap=cmap, linewidths=1, linecolor='black')
        plt.xlabel('Month')
        plt.ylabel('Year')

        plt.savefig(fjpg, dpi=300)
        plt.close(h)

        sum_by_month = np.zeros(12)
        for imonth in range(12):
            sum_by_month[imonth] = np.sum(dfall_month.loc[:, imonth + 1])
        sum_by_month = pd.DataFrame([sum_by_month], columns=month)
        dfall_month = pd.concat([dfall_month, sum_by_month])
        sum_by_year = np.zeros(len(dfall_month.index))
        ihere = 0
        for index, row in dfall_month.iterrows():
            sum_by_year[ihere] = np.sum(row)
            ihere = ihere + 1
        sum_by_year = pd.DataFrame(data=sum_by_year, columns=['All'], index=dfall_month.index)
        dfall_month = pd.concat([dfall_month, sum_by_year], axis=1)
        dfall_month.to_csv(fcsv, sep=';')

    def make_validation_mdbfile(self, path_out):
        file_data_valid = os.path.join(path_out, 'DataValid.csv')
        self.mfile.df_validation_valid.to_csv(file_data_valid, sep=';')
        file_data = os.path.join(path_out, 'Data.csv')
        self.mfile.df_validation.to_csv(file_data, sep=';')
        self.plot_all_scatter_plot(path_out)
        self.plot_wavelenght_scatter_plots(path_out, None)
        self.compute_all_statistics(None)
        file_results = os.path.join(path_out, 'Params.csv')
        self.df_valid_stats.to_csv(file_results, sep=';')
        self.plot_all_spectra_param(path_out)

    def save_validation_dfval_data(self, path_out):
        dfval = self.get_df_val()
        if dfval is None:
            return
        file_data_valid = os.path.join(path_out, 'DataValid.csv')
        dfval.to_csv(file_data_valid, sep=';')

    def make_validation_dfval(self, path_out, title, file_name_base, wllist):
        dfval = self.get_df_val()
        if dfval is None:
            return
        self.compute_all_statistics(wllist)
        self.title = title
        self.file_name_base = file_name_base
        self.plot_all_scatter_plot(path_out)
        self.plot_wavelenght_scatter_plots(path_out, wllist)

        self.compute_all_statistics(wllist)
        file_results = os.path.join(path_out, 'Params.csv')
        self.df_valid_stats.to_csv(file_results, sep=';')
        self.plot_all_spectra_param(path_out)

    def make_validation_dfval_insitu(self, path_out, title, file_name_base, wllist):
        dfval = self.get_df_val()
        if dfval is None:
            return
        self.title = title
        self.file_name_base = file_name_base
        # self.plot_all_scatter_plot_insitu(path_out)
        self.compute_all_statistics(wllist)
        file_results = os.path.join(path_out, 'Params.csv')
        self.df_valid_stats.to_csv(file_results, sep=';')
        #
        # wllist = [412.5, 442.5, 490, 510, 560, 665, 710, 780]
        # self.plot_wavelenght_scatter_plots(path_out,wllist)

        self.plot_all_spectra_param(path_out)

    def get_file_name(self, wl):
        if self.mfile is not None:
            return self.mfile.get_file_name(wl)
        else:
            if wl is None:
                file_name = self.file_name_base
            else:
                file_name = self.file_name_base + f'_{wl}'
            return file_name

    def get_file_name_param(self, param):
        if self.mfile is not None:
            return self.mfile.get_file_name_param(param)
        else:
            file_name = self.file_name_base + f'_{param}'
            return file_name
