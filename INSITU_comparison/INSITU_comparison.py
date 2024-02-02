from netCDF4 import Dataset
from datetime import datetime as dt
from datetime import timedelta
from scipy.stats import gaussian_kde
import pytz
import os
import numpy as np
import numpy.ma as ma
import math
from scipy import stats
from pylr2 import regress2


class INSITUCOMPARISON:

    def __init__(self, path_nc):
        self.path_nc = path_nc
        self.dataset_w = None

        self.options = {
            'n_sequences': 60,
            'sensors': {'AERONET': 11,
                        'HYPSTAR': 1538
                        }
        }

        # data for plotting

        ##scatterplots
        self.xdata = None
        self.ydata = None
        self.xregress = None
        self.yregress = None

        self.wl_data = None
        self.ydata_1 = None
        self.ydata_2 = None

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

    def start_file_w(self):
        self.dataset_w = Dataset(self.path_nc, 'w')
        self.dataset_w.createDimension('day_id', None)
        self.dataset_w.createDimension('sequence_id', self.options['n_sequences'])
        ##wavelength dimension
        for sensor in self.options['sensors']:
            dim_name = f'{sensor}_wavelength'
            dim_length = self.options['sensors'][sensor]
            self.dataset_w.createDimension(dim_name, dim_length)

    def create_aeronet_variables(self, areader, date_here):
        ##variables: time, wavelength
        dim_wavelenght = 'AERONET_wavelength'
        variable_time = 'AERONET_time'
        variable_wavelength = 'AERONET_nominal_wavelengths'
        spectral_variables = {'Rho', 'Exact_Wavelengths', 'Lw'}

        # self.dataset_w = Dataset(self.path_nc, 'w')
        var_time = self.dataset_w.createVariable(variable_time, 'f8', ('day_id', 'sequence_id'), zlib=True, complevel=6,
                                                 fill_value=-999.0)
        var_wl = self.dataset_w.createVariable(variable_wavelength, 'f4', (dim_wavelenght,), zlib=True, complevel=6,
                                               fill_value=-999.0)
        var_wl[:] = areader.get_all_nominal_wl()
        row_ini, row_fin = areader.get_indices_date(date_here.strftime('%Y-%m-%d'))

        if row_ini==-1 and row_fin==-1:
            return False

        time_list = areader.extract_time_sublist(row_ini, row_fin)
        idx = 0
        for t in time_list:
            timestamp = t.replace(tzinfo=pytz.utc).timestamp()
            var_time[0, idx] = timestamp
            idx = idx + 1

        for svariable in spectral_variables:
            name_variable = f'AERONET_{svariable}'
            var = self.dataset_w.createVariable(name_variable, 'f4', ('day_id', 'sequence_id', dim_wavelenght),
                                                zlib=True, complevel=6,
                                                fill_value=-999.0)

            array = areader.extract_spectral_data(svariable, False)
            var[0, 0:idx, :] = array[row_ini:row_fin + 1, :]

        return True

    def create_hypstar_variables(self, path_day, date_here):

        ##checking files and times
        time_list = []
        list_files = {}
        time_ini = date_here.replace(tzinfo=pytz.utc, hour=0, minute=0, second=0).timestamp()
        time_end = date_here.replace(tzinfo=pytz.utc, hour=23, minute=59, second=59).timestamp()
        nominal_wavelengths = None
        for name in os.listdir(path_day):
            if name.find('L2A_REF') < 0:
                continue
            file_here = os.path.join(path_day, name)
            dataset = Dataset(file_here, 'r')
            time_stamp = float(dataset.variables['acquisition_time'][0])
            time_here = dt.utcfromtimestamp(time_stamp)
            time_here_str = time_here.strftime('%Y%m%d%H%M')

            if time_ini <= time_stamp <= time_end:
                list_files[time_here_str] = file_here
                time_list.append(time_stamp)

            if nominal_wavelengths is None:
                nominal_wavelengths = np.array(dataset.variables['wavelength'][:])
            dataset.close()

        ##create hypstar variables
        dim_wavelenght = 'HYPSTAR_wavelength'
        variable_time = 'HYPSTAR_time'
        variable_wavelength = 'HYPSTAR_Nominal_Wavelengths'
        spectral_variables = {'water_leaving_radiance': 'Lw'}
        single_variables = {'rhof': 'Rho', 'quality_flag': 'quality_flag'}

        # self.dataset_w = Dataset(self.path_nc, 'w')
        var_time = self.dataset_w.createVariable(variable_time, 'f8', ('day_id', 'sequence_id'), zlib=True, complevel=6,
                                                 fill_value=-999.0)
        var_wl = self.dataset_w.createVariable(variable_wavelength, 'f4', (dim_wavelenght,), zlib=True, complevel=6,
                                               fill_value=-999.0)

        var_wl[:] = nominal_wavelengths[:]

        time_list.sort()
        idx = 0
        for t in time_list:
            var_time[0, idx] = t
            time_here_str = dt.utcfromtimestamp(t).strftime('%Y%m%d%H%M')
            file_here = list_files[time_here_str]
            dataset = Dataset(file_here, 'r')
            for var in spectral_variables:
                name_out = f'HYPSTAR_{spectral_variables[var]}'
                if name_out not in self.dataset_w.variables:
                    self.dataset_w.createVariable(name_out, 'f4', ('day_id', 'sequence_id', dim_wavelenght), zlib=True,
                                                  complevel=6,
                                                  fill_value=-999.0)
                self.dataset_w.variables[name_out][0, idx, :] = dataset.variables[var][:, 0]

            for var in single_variables:
                name_out = f'HYPSTAR_{single_variables[var]}'
                if name_out not in self.dataset_w.variables:
                    self.dataset_w.createVariable(name_out, 'f4', ('day_id', 'sequence_id'), zlib=True,
                                                  complevel=6,
                                                  fill_value=-999.0)
                self.dataset_w.variables[name_out][0, idx] = dataset.variables[var][0]
            dataset.close()
            idx = idx + 1

    def create_hypstar_to_aeronet_variables(self):
        aeronet_exact_wavelenghts = np.array(self.dataset_w.variables['AERONET_Exact_Wavelengths'][0, 0, :])
        hypspar_wavelengths = np.array(self.dataset_w.variables['HYPSTAR_Nominal_Wavelengths'][:])
        dim_wavelenght = 'AERONET_wavelength'

        ##NEAREST NEIGHBOUR METHOD
        indices_nearest = []
        for wl in aeronet_exact_wavelenghts:
            index = np.argmin(np.abs(hypspar_wavelengths - wl))
            indices_nearest.append(index)

        ##spectral variables
        spectral_variables = ['Lw']

        for variable in spectral_variables:
            var_hypstar = self.dataset_w.variables[f'HYPSTAR_{variable}']
            name_variable = f'HYPSTAR_TO_AERONET_{variable}'
            var_new = self.dataset_w.createVariable(name_variable, 'f4', ('day_id', 'sequence_id', dim_wavelenght),
                                                    zlib=True, complevel=6,
                                                    fill_value=-999.0)
            for idx in range(var_hypstar.shape[1]):
                array_hypstar = np.array(var_hypstar[0, idx, :])
                array_hypstar_to_aeronet = array_hypstar[indices_nearest] / 10
                var_new[0, idx, :] = array_hypstar_to_aeronet[:]

    def close_file_w(self):
        self.dataset_w.close()

    def add_match_ups_to_file(self, time_diff_minutes):
        now = dt.now().strftime('%Y%m%d%H%M%S')
        path_copy = os.path.join(os.path.dirname(self.path_nc), f'Temp_{now}.nc')
        self.dataset_w = self.copy_nc(self.path_nc, path_copy)

        self.dataset_w.createDimension('mu_id', None)

        new_variables = ['mu_wavelength', 'mu_day_id', 'mu_AERONET_sequence_id', 'mu_HYPSTAR_sequence_id',
                         'mu_time_diff', 'mu_AERONET_Lw', 'mu_HYPSTAR_TO_AERONET_Lw']

        for new_var_name in new_variables:
            data_type = 'f4'
            if new_var_name.endswith('id'):
                data_type = 'i2'
            self.dataset_w.createVariable(new_var_name, data_type, ('mu_id',), zlib=True, complevel=6)

        max_time_diff = time_diff_minutes * 60

        ndays = self.dataset_w.dimensions['day_id'].size
        nsequences = self.dataset_w.dimensions['sequence_id'].size
        nwl = self.dataset_w.dimensions['AERONET_wavelength'].size
        wl_list = np.array(self.dataset_w.variables['AERONET_nominal_wavelengths'][:])

        imu = 0
        for iday in range(ndays):
            # hypstar_time = np.array(self.dataset_w.variables['HYPSTAR_time'][iday,:])
            aeronet_time = np.array(self.dataset_w.variables['AERONET_time'][iday, :])
            for hsequence in range(nsequences):
                ht = self.dataset_w.variables['HYPSTAR_time'][iday, hsequence]
                qf = self.dataset_w.variables['HYPSTAR_quality_flag'][iday, hsequence]
                if ma.is_masked(ht):
                    continue
                if qf > 0:
                    continue
                asequence = np.argmin(np.abs(ht - aeronet_time))

                at = aeronet_time[asequence]
                time_diff = abs(aeronet_time[asequence] - ht)
                if time_diff <= max_time_diff:
                    print(dt.utcfromtimestamp(float(ht)), '<->', dt.utcfromtimestamp(float(at)), '->', time_diff)
                    imu_ini = imu
                    imu_fin = imu + nwl
                    imu = imu_fin
                    self.dataset_w.variables['mu_wavelength'][imu_ini:imu_fin] = wl_list[:]
                    self.dataset_w.variables['mu_day_id'][imu_ini:imu_fin] = iday
                    self.dataset_w.variables['mu_AERONET_sequence_id'][imu_ini:imu_fin] = asequence
                    self.dataset_w.variables['mu_HYPSTAR_sequence_id'][imu_ini:imu_fin] = hsequence
                    self.dataset_w.variables['mu_time_diff'][imu_ini:imu_fin] = time_diff
                    self.dataset_w.variables['mu_AERONET_Lw'][imu_ini:imu_fin] = self.dataset_w.variables['AERONET_Lw'][
                                                                                 iday, asequence, :]
                    self.dataset_w.variables['mu_HYPSTAR_TO_AERONET_Lw'][imu_ini:imu_fin] = self.dataset_w.variables[
                                                                                                'HYPSTAR_TO_AERONET_Lw'][
                                                                                            iday, asequence, :]

        self.close_file_w()
        os.rename(path_copy, self.path_nc)

    def copy_nc(self, ifile, ofile):
        with Dataset(ifile) as src:
            dst = Dataset(ofile, 'w', format='NETCDF4')
            # copy global attributes all at once via dictionary
            dst.setncatts(src.__dict__)
            # copy dimensions
            for name, dimension in src.dimensions.items():
                dst.createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else None))
            # copy all file data except for the excluded
            for name, variable in src.variables.items():
                dst.createVariable(name, variable.datatype, variable.dimensions)
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)

                dst[name][:] = src[name][:]
        return dst

    def compute_stats(self, use_log_scale, use_rhow):

        self.valid_stats['N'] = len(self.xdata)
        if self.valid_stats['N'] == 0:
            for key in self.valid_stats:
                self.valid_stats[key] = np.nan
            return

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

        # type I regression
        slope_I, intercept_I, r_value, p_value, std_err = stats.linregress(xdatal, ydatal)
        # self.xregress, self.yregress = self.get_regression_line(xdatal, ydatal, slope, intercept, minxy, maxxy)

        # type II regression
        results = regress2(np.array(xdatal, dtype=np.float64), np.array(ydatal, dtype=np.float64),
                           _method_type_2="reduced major axis")
        slope = results['slope']
        intercept = results['intercept']
        r_value = results['r']
        std_slope = results['std_slope']
        std_intercept = results['std_intercept']
        self.xregress, self.yregress = self.get_regression_line(xdatal, ydatal, slope, intercept, minxy, maxxy)

        self.valid_stats['slope_I'] = slope_I
        self.valid_stats['intercept_I'] = intercept_I
        self.valid_stats['slope'] = slope
        self.valid_stats['intercept'] = intercept
        self.valid_stats['PCC(r)'] = r_value
        self.valid_stats['p_value'] = p_value
        self.valid_stats['std_err'] = std_err
        self.valid_stats['std_slope'] = std_err
        self.valid_stats['std_intercept'] = std_intercept

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

        self.valid_stats['RMSD'] = self.rmse(sat_obs, ref_obs)
        ref_mean = np.mean(ref_obs)
        sat_mean = np.mean(sat_obs)
        self.valid_stats['XAVG'] = ref_mean
        self.valid_stats['YAVG'] = sat_mean

        # CPRMSE
        xdiff = ref_obs - ref_mean
        ydiff = sat_obs - sat_mean
        cprmse = self.rmse(ydiff, xdiff)
        self.valid_stats['CPRMSE'] = cprmse

        bias = np.mean(sat_obs - ref_obs)
        self.valid_stats['BIAS'] = bias

        mae = np.mean(np.abs(sat_obs - ref_obs))
        self.valid_stats['MAE'] = mae

        self.valid_stats['DETER(r2)'] = r_value * r_value

    def rmse(self, predictions, targets):
        return np.sqrt(((np.asarray(predictions) - np.asarray(targets)) ** 2).mean())

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

    def get_str_stat_list(self, stat_list, wl):
        str0 = ''
        for stat in stat_list:
            if len(str0) > 0:
                str0 = f'{str0}\n'
            if stat == 'N':
                val = self.valid_stats[stat]
                str0 = f'{str0}N={val}'
            if stat == 'NMATCH-UPS':
                val = self.valid_stats['N']
                self.dataset_w = Dataset(self.path_nc, 'r')
                nwl = np.unique(np.array(self.dataset_w.variables['mu_wavelength'])).shape[0]
                self.dataset_w.close()

                val = val / nwl
                str0 = f'{str0}N={val:.0f}'
            if stat == 'r2':
                val = self.valid_stats['DETER(r2)']
                str0 = f'{str0}R\u00b2={val:.2f}'
            if stat == 'RMSD' or stat == 'BIAS':
                val = self.valid_stats[stat]
                if stat == 'BIAS':
                    stat = stat.lower()
                # str0 = f'{str0}{stat}={val:.1e}'
                str0 = f'{str0}{stat}={val:.2f}'
            if stat == 'RPD' or stat == 'APD':
                val = self.valid_stats[stat]
                str0 = f'{str0}{stat}={val:.0f}%'
            if stat == 'WL':
                wls = self.get_wl_str_from_wl(wl)
                str0 = f'{str0}{wls} nm'
                # str0 = f'{str0}{wl:.2f} nm'

        return str0

    def get_wl_str_from_wl(self, wl_value):
        # wl_sat = np.array(self.mrfile.nc.variables['satellite_bands'])
        # index_sat = np.argmin(np.abs(wl_sat - wl_value))
        # wl_sat_value = wl_sat[index_sat]
        wl_sat_value_str = f'{wl_value:.2f}'
        if wl_sat_value_str.endswith('.00'):
            return wl_sat_value_str[:-3]
        else:
            return wl_sat_value_str

    def set_data_scatterplot(self, wl):
        if wl=='rho':
            return
        self.dataset_w = Dataset(self.path_nc, 'r')
        if wl is None:
            self.xdata = np.array(self.dataset_w.variables['mu_AERONET_Lw'])
            self.ydata = np.array(self.dataset_w.variables['mu_HYPSTAR_TO_AERONET_Lw'])
        else:
            xdata = np.array(self.dataset_w.variables['mu_AERONET_Lw'])
            ydata = np.array(self.dataset_w.variables['mu_HYPSTAR_TO_AERONET_Lw'])
            wl_array = np.array(self.dataset_w.variables['mu_wavelength'])
            self.xdata = xdata[wl_array == wl]
            self.ydata = ydata[wl_array == wl]
        self.close_file_w()


    def plot_all_spectra(self):
        path_plots = '/mnt/c/DATA_LUIS/ESA-POP_WORK/Technical_Meeting_2024Jan31'
        self.dataset_w = Dataset(self.path_nc, 'r')
        self.wl_data = np.array(self.dataset_w.variables['AERONET_nominal_wavelengths'][:])
        nwl = len(self.wl_data)

        hypstar_data = np.array(self.dataset_w.variables['mu_HYPSTAR_TO_AERONET_Lw'][:])
        aeronet_data = np.array(self.dataset_w.variables['mu_AERONET_Lw'][:])
        ndata = hypstar_data.shape[0]

        nspectra = ndata/nwl
        imu = 0
        for ispectra in range(int(nspectra)):
            file_out = os.path.join(path_plots,f'ComparisonSpectra_{ispectra+1}.tif')
            imu_ini = imu
            imu_end = imu+nwl
            self.ydata_1 = aeronet_data[imu_ini:imu_end]
            self.ydata_2 = hypstar_data[imu_ini:imu_end]
            iday = int(self.dataset_w.variables['mu_day_id'][imu_ini])
            aseq = int(self.dataset_w.variables['mu_AERONET_sequence_id'][imu_ini])
            hseq = int(self.dataset_w.variables['mu_HYPSTAR_sequence_id'][imu_ini])
            atime = dt.utcfromtimestamp(float(self.dataset_w.variables['AERONET_time'][iday,aseq]))
            htime = dt.utcfromtimestamp(float(self.dataset_w.variables['HYPSTAR_time'][iday, hseq]))
            atimestr = atime.strftime('%Y-%m-%d %H:%M')
            htimestr = htime.strftime('%Y-%m-%d %H:%M')
            legend = [f'AERONET-OC({atimestr})',f'HYPSTAR({htimestr})']
            title = f'SPECTRA #{ispectra+1}'
            self.plot_spectra_impl(file_out,legend,title)
            imu = imu_end

        self.close_file_w()

    def plot_spectra_impl(self,file_out,legend,title):
        xlabel = 'Wavelength (nm)'
        ylabel = 'Lw'
        from MDB_reader.PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        pspectra.xdata = self.wl_data
        hline = pspectra.plot_single_line(self.ydata_1, 'red', 'solid', 1, 'o', 8)
        hline = pspectra.plot_single_line(self.ydata_2, 'blue', 'solid', 1, 'o', 8)

        pspectra.set_xticks(self.wl_data,self.wl_data,90,8)
        pspectra.set_yaxis_title(ylabel)
        pspectra.set_xaxis_title(xlabel)
        pspectra.set_title(title)

        pspectra.set_grid()
        pspectra.legend_options['loc'] = 'lower center'
        pspectra.legend_options['bbox_to_anchor'] =  (0.1, -0.4, 0.80, 0.1)
        pspectra.legend_options['ncols'] = 2
        pspectra.set_legend(legend)
        pspectra.set_tigth_layout()
        pspectra.save_fig(file_out)




    def plot_all_scatterplots_wl(self):
        self.dataset_w = Dataset(self.path_nc,'r')
        wl_list = np.array(self.dataset_w.variables['AERONET_nominal_wavelengths'][:])
        self.close_file_w()
        for wl in wl_list:
            self.plot_scatterplot(wl)

    def plot_rho_scatterplot(self):
        self.dataset_w = Dataset(self.path_nc, 'r')
        # hypstar_Rho = np.array(self.dataset_w.variables['HYPSTAR_Rho'][:])
        # aeronet_Rho = np.array(self.dataset_w.variables['AERONET_Rho'][:,:,0])
        ndata = self.dataset_w.variables['mu_day_id'].shape[0]
        nwl = self.dataset_w.variables['AERONET_nominal_wavelengths'].shape[0]
        nspectra = int(ndata/nwl)

        #print(ndata)
        self.xdata=[0]*nspectra
        self.ydata=[0]*nspectra

        imu = 0
        for ispectra in range(int(nspectra)):

            iday = int(self.dataset_w.variables['mu_day_id'][imu])
            ihs = int(self.dataset_w.variables['mu_HYPSTAR_sequence_id'][imu])
            ias = int(self.dataset_w.variables['mu_AERONET_sequence_id'][imu])

            self.xdata[ispectra] = float(self.dataset_w.variables['AERONET_Rho'][iday,ias,0])
            self.ydata[ispectra] = float(self.dataset_w.variables['HYPSTAR_Rho'][iday, ihs])

            imu = imu + nwl

        self.close_file_w()

        self.plot_scatterplot('rho')

    def plot_scatterplot(self, wl):
        path_plots = '/mnt/c/DATA_LUIS/ESA-POP_WORK/Technical_Meeting_2024Jan31'
        marker = 'o'
        markersize = 20
        edgecolor = None
        linewidth = 0

        fontsizeaxis = 12

        if wl=='rho':
            xlabel = 'AERONET-OC Rho'
            ylabel = 'HYPSTAR Rho'
            min_xy = 0.025
            max_xy = 0.030
            ticks = [0.025,0.026,0.027,0.028,0.029,0.030]
        else:
            xlabel = r'AERONET-OC Lw (µW/cm$^2$/sr/nm)'
            ylabel = r'HYPSTAR Lw (µW/cm$^2$/sr/nm)'
            min_xy = None
            max_xy = None
            ticks  = None
            if wl is None or wl<600:
                min_xy = 0
                max_xy = 3
                ticks = [0, 0.5, 1, 1.5, 2, 2.5, 3]
            if wl>600: #and wl<700:
                min_xy = 0
                max_xy = 0.6
                ticks = [0, 0.2, 0.4, 0.6]
            # if wl>700 and wl<800:
            #     min_xy = 0
            #     max_xy = 0.05
            #     ticks = [0, 0.01, 0.02, 0.03,0.04,0.05]
            # if wl>800:
            #     min_xy = 0
            #     max_xy = 0.01
            #     ticks = [0, 0.002, 0.004, 0.006, 0.008, 0.01]
        stat_list = ['NMATCH-UPS', 'r2', 'RMSD', 'BIAS']
        stats_xpos = 0.05
        stats_ypos = 0.80
        self.set_data_scatterplot(wl)

        # print(self.xdata)
        # print(self.ydata)
        # print(min_xy,max_xy)

        self.compute_stats(False, False)

        from MDB_reader.PlotScatter import PlotScatter
        import MDB_reader.MDBPlotDefaults as defaults
        plot = PlotScatter()
        plot.close_plot()
        plot.start_plot()

        # global scatterplot
        if wl is None:
            file_out = os.path.join(path_plots, 'GlobalScatterplot.tif')
            self.dataset_w = Dataset(self.path_nc, 'r')
            groupData = np.array(self.dataset_w.variables['mu_wavelength'])
            self.dataset_w.close()
            groupValues = np.unique(groupData)
            ngroup = len(groupValues)
            str_legend = []
            for g in groupValues:
                str_legend.append(f'{g:.2f}')

            for idx in range(ngroup):  # groupValues:
                g = groupValues[idx]
                color = defaults.get_color_ref(g)
                xhere = self.xdata[groupData == g]
                yhere = self.ydata[groupData == g]
                plot.plot_data(xhere, yhere, marker, markersize, color, edgecolor, linewidth)

            plot.set_limits(min_xy, max_xy)
            plot.set_ticks(ticks, fontsizeaxis)
            plot.set_legend(str_legend)
            str_stats = self.get_str_stat_list(stat_list, None)
        else:
            if wl!='rho':
                wls = self.get_wl_str_from_wl(wl)
                file_out = os.path.join(path_plots, f'Wl_Scatterplot_{wls}.tif')
            else:
                file_out = os.path.join(path_plots, f'Rho_Scatterplot.tif')
            ##DENSITY
            xhere = np.asarray(self.xdata, dtype=np.float)
            yhere = np.asarray(self.ydata, dtype=np.float)
            xy = np.vstack([xhere, yhere])
            try:
                z = gaussian_kde(xy)(xy)
                idx = z.argsort()
                xhere, yhere, z = xhere[idx], yhere[idx], z[idx]
                plot.set_cmap('jet')
                plot.plot_data(xhere, yhere, marker, markersize, z, edgecolor, linewidth)
            except:
                plot.plot_data(xhere, yhere, marker, markersize, 'k', edgecolor, linewidth)
            stat_list = ['N', 'r2', 'RMSD', 'BIAS']
            str_stats = self.get_str_stat_list(stat_list, None)
            if wl!='rho':
                title = f'{wls} nm'
                plot.set_title(title)
            if min_xy is not None and max_xy is not None:
                plot.set_limits(min_xy, max_xy)
            if ticks is not None:
                plot.set_ticks(ticks, fontsizeaxis)

        plot.set_xaxis_title(xlabel)
        plot.set_yaxis_title(ylabel)

        #plot.plot_regress_line(self.xregress, self.yregress, 'k')
        plot.plot_identity_line()

        plot.plot_text(stats_xpos, stats_ypos, str_stats)
        plot.set_equal_apect()

        plot.save_fig(file_out)
        plot.close_plot()
