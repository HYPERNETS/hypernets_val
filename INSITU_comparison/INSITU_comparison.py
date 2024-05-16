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

from ComparisonOptions import ComparisonOptions


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

        self.coptions = ComparisonOptions()

        # data for plotting

        ##scatterplots
        self.xdata = None
        self.ydata = None
        self.xregress = None
        self.yregress = None

        self.wl_data = None
        self.ydata_1 = None
        self.ydata_2 = None
        self.spectra_stats = {}

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
        # spectral_variables = {'Rho', 'Exact_Wavelengths', 'Lw'}
        spectral_variables = self.coptions.aeronet_spectral_variables()
        single_variables = self.coptions.aeronet_single_variables()

        # self.dataset_w = Dataset(self.path_nc, 'w')
        var_time = self.dataset_w.createVariable(variable_time, 'f8', ('day_id', 'sequence_id'), zlib=True, complevel=6,
                                                 fill_value=-999.0)
        var_time.long_name = 'UNIX Time Instant'
        var_time.units = 'Seconds since 1970-01-01 00:00:00.000 UTC'
        var_wl = self.dataset_w.createVariable(variable_wavelength, 'f4', (dim_wavelenght,), zlib=True, complevel=6,
                                               fill_value=-999.0)
        var_wl.long_name = 'AERONET-OC Nominal Wavelengths'
        var_wl.units = 'nm'
        var_wl[:] = areader.get_all_nominal_wl()
        row_ini, row_fin = areader.get_indices_date(date_here.strftime('%Y-%m-%d'))

        if row_ini == -1 and row_fin == -1:
            print(f'[ERROR] Aeronet data were not available for day: {date_here}')
            return False

        time_list = areader.extract_time_sublist(row_ini, row_fin)
        idx = 0
        for t in time_list:
            timestamp = t.replace(tzinfo=pytz.utc).timestamp()
            var_time[0, idx] = timestamp
            idx = idx + 1

        for svariable in spectral_variables:
            name_variable = f'AERONET_{svariable}'
            # print(f'[INFO] Creating spectral variable:{name_variable}')
            var = self.dataset_w.createVariable(name_variable, 'f4', ('day_id', 'sequence_id', dim_wavelenght),
                                                zlib=True, complevel=6,
                                                fill_value=-999.0)
            var.long_name = spectral_variables[svariable]['long_name']
            if 'units' in spectral_variables[svariable]:
                var.units = spectral_variables[svariable]['units']

            array = areader.extract_spectral_data(svariable, False)
            var[0, 0:idx, :] = array[row_ini:row_fin + 1, :]

        ##rrs
        name_variable = 'AERONET_Rrs'
        var = self.dataset_w.createVariable(name_variable, 'f4', ('day_id', 'sequence_id', dim_wavelenght),
                                            zlib=True, complevel=6,
                                            fill_value=-999.0)
        var.long_name = 'Reflectance of the water column at the surface'
        var.units = 'sr-1'
        array = areader.extract_rrs(False)
        var[0, 0:idx, :] = array[row_ini:row_fin + 1, :]

        for svariable in single_variables:
            name_variable = f'AERONET_{svariable}'
            # print(f'[INFO] Creating single variable:{name_variable}')
            var = self.dataset_w.createVariable(name_variable, 'f4', ('day_id', 'sequence_id'),
                                                zlib=True, complevel=6,
                                                fill_value=-999.0)
            var.long_name = single_variables[svariable]['long_name']
            if 'units' in single_variables[svariable]:
                var.units = single_variables[svariable]['units']
            if svariable == 'Observing_Azimuth_Angle':
                var[0, 0:idx] = 270
            elif svariable == 'Observing_Zenith_Angle':
                var[0, 0:idx] = 40
            else:
                array = areader.extract_single_data(svariable)
                var[0, 0:idx] = array[row_ini:row_fin + 1]

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
        spectral_l2_variables = self.coptions.hypstar_l2_spectral_variables()
        spectral_l1_variables = self.coptions.hypstar_l1_spectral_variables()

        single_l2_variables = self.coptions.hypstar_l2_single_variables()
        single_l1_variables = self.coptions.hypstar_l1_single_variables()

        # self.dataset_w = Dataset(self.path_nc, 'w')
        var_time = self.dataset_w.createVariable(variable_time, 'f8', ('day_id', 'sequence_id'), zlib=True, complevel=6,
                                                 fill_value=-999.0)
        var_time.long_name = 'UNIX Time Instant'
        var_time.units = 'Seconds since 1970-01-01 00:00:00.000 UTC'
        var_wl = self.dataset_w.createVariable(variable_wavelength, 'f4', (dim_wavelenght,), zlib=True, complevel=6,
                                               fill_value=-999.0)
        var_wl.long_name = 'HYPSTAR Nominal Wavelengths'
        var_wl.units = 'nm'
        var_wl[:] = nominal_wavelengths[:]

        time_list.sort()
        idx = 0
        for t in time_list:
            var_time[0, idx] = t
            time_here_str = dt.utcfromtimestamp(t).strftime('%Y%m%d%H%M')
            file_here = list_files[time_here_str]
            dataset = Dataset(file_here, 'r')
            for var in spectral_l2_variables:
                name_out = f'HYPSTAR_{var}'
                if name_out not in self.dataset_w.variables:
                    var_new = self.dataset_w.createVariable(name_out, 'f4', ('day_id', 'sequence_id', dim_wavelenght),
                                                            zlib=True,
                                                            complevel=6,
                                                            fill_value=-999.0)

                    for at in dataset.variables[var].ncattrs():
                        if at == '_FillValue' or at == 'add_offset' or at == 'scale_factor':
                            continue
                        var_new.setncattr(at, dataset.variables[var].getncattr(at))
                var_out = dataset.variables[var][:, 0]
                if '_FillValue' in dataset.variables[var].ncattrs():
                    fill_value_orig = dataset.variables[var]._FillValue
                    var_out[var_out == fill_value_orig] = -999.0
                self.dataset_w.variables[name_out][0, idx, :] = var_out

            for var in single_l2_variables:
                name_out = f'HYPSTAR_{var}'
                if name_out not in self.dataset_w.variables:
                    var_new = self.dataset_w.createVariable(name_out, 'f4', ('day_id', 'sequence_id'), zlib=True,
                                                            complevel=6,
                                                            fill_value=-999.0)
                    for at in dataset.variables[var].ncattrs():
                        if at == '_FillValue' or at == 'add_offset' or at == 'scale_factor':
                            continue
                        var_new.setncattr(at, dataset.variables[var].getncattr(at))
                var_out = dataset.variables[var][0]
                if '_FillValue' in dataset.variables[var].ncattrs():
                    fill_value_orig = dataset.variables[var]._FillValue
                    if var_out == fill_value_orig:
                        var_out = -999.0
                self.dataset_w.variables[name_out][0, idx] = var_out
            dataset.close()

            file_l1 = file_here.replace('L2A_REF', 'L1C_ALL')
            if os.path.exists(file_l1):
                datasetL1 = Dataset(file_l1)
                for var in spectral_l1_variables:
                    name_out = f'HYPSTAR_{var}'
                    if name_out not in self.dataset_w.variables:
                        var_new = self.dataset_w.createVariable(name_out, 'f4',
                                                                ('day_id', 'sequence_id', dim_wavelenght),
                                                                zlib=True,
                                                                complevel=6,
                                                                fill_value=-999.0)
                        if var != 'std_upwelling_radiance':
                            for at in datasetL1.variables[var].ncattrs():
                                if at == '_FillValue' or at == 'add_offset' or at == 'scale_factor':
                                    continue
                                var_new.setncattr(at, datasetL1.variables[var].getncattr(at))
                        else:
                            var_new.long_name = 'standard deviation for upwelling radiance'
                            var_new.units = 'mW m^-2 nm^-1 sr^-1'
                    if var == 'std_upwelling_radiance':
                        var_out = np.ma.std(np.ma.array(datasetL1.variables['upwelling_radiance'][:, :]), axis=1)
                        var_out = var_out.filled(-999.0)
                        self.dataset_w.variables[name_out][0, idx, :] = var_out
                    elif datasetL1.variables[var].ndim == 2:
                        var_out = np.ma.mean(np.ma.array(datasetL1.variables[var][:, :]), axis=1)
                        var_out = var_out.filled(-999.0)
                        self.dataset_w.variables[name_out][0, idx, :] = var_out
                    elif datasetL1.variables[var].ndim == 1 and datasetL1.variables[var].dimensions[0] == 'wavelength':
                        var_out = datasetL1.variables[var][:]
                        if '_FillValue' in datasetL1.variables[var].ncattrs():
                            fill_value_orig = datasetL1.variables[var]._FillValue
                            var_out[var_out == fill_value_orig] = -999.0
                        self.dataset_w.variables[name_out][0, idx, :] = var_out

                for var in single_l1_variables:
                    name_out = f'HYPSTAR_{var}'
                    if name_out not in self.dataset_w.variables:
                        var_new = self.dataset_w.createVariable(name_out, 'f4', ('day_id', 'sequence_id'), zlib=True,
                                                                complevel=6,
                                                                fill_value=-999.0)
                        for at in datasetL1.variables[var].ncattrs():
                            if at == '_FillValue' or at == 'add_offset' or at == 'scale_factor':
                                continue
                            var_new.setncattr(at, datasetL1.variables[var].getncattr(at))
                    var_out = np.ma.mean(np.ma.array(datasetL1.variables[var][:]))
                    if ma.is_masked(var_out):
                        var_out = -999.0
                    self.dataset_w.variables[name_out][0, idx] = var_out
                datasetL1.close()

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
        spectral_variables = self.coptions.hypstar_aeronet_spectral_variables()

        for variable in spectral_variables:
            var_hypstar = self.dataset_w.variables[f'HYPSTAR_{variable}']
            var_name_out = spectral_variables[variable]['name']
            name_variable = f'HYPSTAR_TO_AERONET_{var_name_out}'
            var_new = self.dataset_w.createVariable(name_variable, 'f4', ('day_id', 'sequence_id', dim_wavelenght),
                                                    zlib=True, complevel=6,
                                                    fill_value=-999.0)
            scale_factor = 1.0
            if 'scale_factor' in spectral_variables[variable]:
                scale_factor = spectral_variables[variable]['scale_factor']
            var_new.setncatts(var_hypstar.__dict__)
            for at in spectral_variables[variable]:
                if at == 'name' or at == 'scale_factor':
                    continue
                var_new.setncattr(at, spectral_variables[variable][at])

            for idx in range(var_hypstar.shape[1]):
                array_hypstar = np.array(var_hypstar[0, idx, :])
                array_hypstar_to_aeronet = array_hypstar[indices_nearest]
                array_hypstar_to_aeronet[array_hypstar_to_aeronet != -999 - 0] = array_hypstar_to_aeronet[
                                                                                     array_hypstar_to_aeronet != -999 - 0] * scale_factor
                var_new[0, idx, :] = array_hypstar_to_aeronet[:]

    def close_file_w(self):
        self.dataset_w.close()

    def add_match_ups_to_file(self, time_diff_minutes):
        now = dt.now().strftime('%Y%m%d%H%M%S')
        path_copy = os.path.join(os.path.dirname(self.path_nc), f'Temp_{now}.nc')
        self.dataset_w = self.copy_nc(self.path_nc, path_copy)

        self.dataset_w.createDimension('mu_id', None)

        basic_mu_variables = self.coptions.basic_mu_variables()
        mu_variables = self.coptions.mu_variables()
        mu_variables_keys = self.coptions.mu_variables_keys()

        # new_variables = ['mu_wavelength', 'mu_day_id', 'mu_AERONET_sequence_id', 'mu_HYPSTAR_sequence_id',
        #                  'mu_time_diff', 'mu_AERONET_Lw', 'mu_HYPSTAR_TO_AERONET_Lw']

        for new_var_name in basic_mu_variables:
            data_type = 'f4'
            if new_var_name.endswith('id'):
                data_type = 'i2'
            var_new = self.dataset_w.createVariable(new_var_name, data_type, ('mu_id',), zlib=True, complevel=6)
            for at in basic_mu_variables[new_var_name]:
                var_new.setncattr(at, basic_mu_variables[new_var_name][at])

        for key in mu_variables_keys:
            for var in mu_variables:
                variable_in = f'{key}_{var}'
                if variable_in in self.dataset_w.variables:
                    new_var_name = f'mu_{variable_in}'
                    var_new = self.dataset_w.createVariable(new_var_name, 'f4', ('mu_id',), zlib=True, complevel=6)
                    for at in self.dataset_w.variables[variable_in].ncattrs():
                        if at == '_FillValue' or at == 'add_offset' or at == 'scale_factor':
                            continue
                        val = self.dataset_w.variables[variable_in].getncattr(at)
                        var_new.setncattr(at, val)

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

                    str = f'HYPSTAR time: {dt.utcfromtimestamp(float(ht))} AERONET-OC time: {dt.utcfromtimestamp(float(at))} Time diff: {time_diff} seconds'
                    print(f'[INFO] --> {str}')
                    # print(dt.utcfromtimestamp(float(ht)), '<->', dt.utcfromtimestamp(float(at)), '->', time_diff)
                    imu_ini = imu
                    imu_fin = imu + nwl
                    imu = imu_fin
                    self.dataset_w.variables['mu_wavelength'][imu_ini:imu_fin] = wl_list[:]
                    self.dataset_w.variables['mu_day_id'][imu_ini:imu_fin] = iday
                    self.dataset_w.variables['mu_AERONET_sequence_id'][imu_ini:imu_fin] = asequence
                    self.dataset_w.variables['mu_HYPSTAR_sequence_id'][imu_ini:imu_fin] = hsequence
                    self.dataset_w.variables['mu_time_diff'][imu_ini:imu_fin] = time_diff

                    for key in mu_variables_keys:
                        for var in mu_variables:
                            variable_in = f'{key}_{var}'
                            if variable_in in self.dataset_w.variables:
                                new_var_name = f'mu_{variable_in}'
                                self.dataset_w.variables[new_var_name][imu_ini:imu_fin] = self.dataset_w.variables[
                                                                                              variable_in][
                                                                                          iday, asequence, :]

                    # self.dataset_w.variables['mu_AERONET_Lw'][imu_ini:imu_fin] = self.dataset_w.variables['AERONET_Lw'][
                    #                                                              iday, asequence, :]
                    # self.dataset_w.variables['mu_HYPSTAR_TO_AERONET_Lw'][imu_ini:imu_fin] = self.dataset_w.variables[
                    #                                                                             'HYPSTAR_TO_AERONET_Lw'][
                    #                                                                         iday, asequence, :]

        self.close_file_w()
        os.rename(path_copy, self.path_nc)
        print(f'[INFO] --> Completed')

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
        self.xregress, self.yregress = self.get_regression_line(xdatal, ydatal, slope_I, intercept_I, minxy, maxxy)

        # type II regression
        from pylr2 import regress2
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
        # self.valid_stats['slope'] = slope
        # self.valid_stats['intercept'] = intercept
        self.valid_stats['PCC(r)'] = r_value
        self.valid_stats['p_value'] = p_value
        self.valid_stats['std_err'] = std_err
        self.valid_stats['std_slope'] = std_err
        # self.valid_stats['std_intercept'] = std_intercept

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
                str0 = f'{str0}{stat}={val:.4f}'
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

    def set_data_scatterplot_single(self, xvariable, yvariable, wlvariable):
        self.dataset_w = Dataset(self.path_nc, 'r')
        ndata = self.dataset_w.variables['mu_day_id'].shape[0]
        # nwl = self.dataset_w.variables['AERONET_nominal_wavelengths'].shape[0]
        nwl = self.dataset_w.variables[wlvariable].shape[0]
        nspectra = int(ndata / nwl)

        self.xdata = [0] * nspectra
        self.ydata = [0] * nspectra

        xvariable_ins = xvariable.split('_')[0]
        ixvariable = f'mu_{xvariable_ins}_sequence_id'
        yvariable_ins = yvariable.split('_')[0]
        iyvariable = f'mu_{yvariable_ins}_sequence_id'

        imu = 0

        for ispectra in range(int(nspectra)):
            iday = int(self.dataset_w.variables['mu_day_id'][imu])
            ix = int(self.dataset_w.variables[ixvariable][imu])
            iy = int(self.dataset_w.variables[iyvariable][imu])

            from datetime import datetime as dt

            time = dt.utcfromtimestamp(float(self.dataset_w.variables['HYPSTAR_time'][iday, ix]))
            date_str = time.strftime('%Y-%m-%d')
            time_str = time.strftime('%H:%M')
            time_diff = float(self.dataset_w.variables['mu_time_diff'][imu])
            # line = f'{date_str};{time_str};{time_diff};{float(self.dataset_w.variables[xvariable][iday, ix])};{float(self.dataset_w.variables[yvariable][iday, ix])}'
            #
            # print(line)

            if self.dataset_w.variables[xvariable].ndim == 3:
                self.xdata[ispectra] = float(self.dataset_w.variables[xvariable][iday, ix, 0])
            else:
                self.xdata[ispectra] = float(self.dataset_w.variables[xvariable][iday, ix])

            if self.dataset_w.variables[yvariable].ndim == 3:
                self.ydata[ispectra] = float(self.dataset_w.variables[yvariable][iday, iy, 0])
            else:
                self.ydata[ispectra] = float(self.dataset_w.variables[yvariable][iday, iy])

            imu = imu + nwl

        self.close_file_w()

    def set_data_scatterplot_spectral(self, wl, xvariable, yvariable, wlvariable):

        self.dataset_w = Dataset(self.path_nc, 'r')

        day_id = np.array(self.dataset_w['mu_day_id'][:])
        array_1020 = np.array(self.dataset_w['mu_HYPSTAR_TO_AERONET_Lw'][:])
        valid = np.ones(array_1020.shape)
        valid[array_1020>100000] = 0
        invalid_days = [21,30,32,35,63,66,80]
        for iday in invalid_days:
            valid[day_id==iday] = 0
        # valid[day_id==21] = 0
        # valid[day_id==35] = 0
        # valid[day_id==80] = 0


        if wl is None:
            self.xdata = np.array(self.dataset_w.variables[xvariable])
            self.ydata = np.array(self.dataset_w.variables[yvariable])
            self.xdata = self.xdata[valid==1]
            self.ydata = self.ydata[valid==1]
        else:
            xdata = np.array(self.dataset_w.variables[xvariable])
            ydata = np.array(self.dataset_w.variables[yvariable])
            wl_array = np.array(self.dataset_w.variables[wlvariable])  # 'mu_wavelength'

            xdata = xdata[valid == 1]
            ydata = ydata[valid == 1]
            wl_array = wl_array[valid==1]


            self.xdata = xdata[wl_array == wl]
            self.ydata = ydata[wl_array == wl]








        # print(self.xdata.shape, self.ydata.shape)
        # self.xdata = self.xdata[np.where(self.xdata < 100000)]
        # self.ydata = self.ydata[np.where(self.xdata < 100000)]
        # print(self.xdata.shape,self.ydata.shape)

        self.close_file_w()

        return valid

    def set_spectra_stats(self, variable1, variable2, wlvariable):
        self.dataset_w = Dataset(self.path_nc, 'r')
        wldata = np.array(self.dataset_w.variables[wlvariable])
        wlvalues = np.unique(wldata)
        self.wl_data = wlvalues
        nwl = len(wlvalues)
        self.spectra_stats = {
            variable1: {
                'mean': np.zeros((nwl,)),
                'std': np.zeros((nwl,)),
                'median': np.zeros((nwl,)),
                'p25': np.zeros((nwl,)),
                'p75': np.zeros((nwl,))
            },
            variable2: {
                'mean': np.zeros((nwl,)),
                'std': np.zeros((nwl,)),
                'median': np.zeros((nwl,)),
                'p25': np.zeros((nwl,)),
                'p75': np.zeros((nwl,))
            }
        }
        variables = [variable1, variable2]

        day_id = np.array(self.dataset_w['mu_day_id'][:])
        array_1020 = np.array(self.dataset_w['mu_HYPSTAR_TO_AERONET_Lw'][:])
        valid = np.ones(array_1020.shape)
        valid[array_1020 > 100000] = 0
        invalid_days = [21, 30, 32, 35, 63, 66, 80]
        for iday in invalid_days:
            valid[day_id == iday] = 0

        wldata = wldata[valid==1]

        for variable in variables:
            for idx in range(nwl):
                wl = wlvalues[idx]
                values_all = np.array(self.dataset_w.variables[variable])
                values_all = values_all[valid==1]
                values = values_all[wldata == wl]
                values = values[values < 100000]
                self.spectra_stats[variable]['mean'][idx] = np.mean(values)
                self.spectra_stats[variable]['std'][idx] = np.std(values)
                #print(idx, wl, np.median(values))
                self.spectra_stats[variable]['median'][idx] = np.median(values)
                self.spectra_stats[variable]['p25'][idx] = np.percentile(values, 25)
                self.spectra_stats[variable]['p75'][idx] = np.percentile(values, 75)

        self.close_file_w()

    def plot_spectra_stats(self,file_out,legend,ylabel,title):

        self.plot_spectra_impl(file_out, legend, ylabel,title)

    def plot_all_spectra(self):
        path_plots = '/mnt/c/DATA_LUIS/ESA-POP_WORK/Technical_Meeting_2024Jan31'
        self.dataset_w = Dataset(self.path_nc, 'r')
        self.wl_data = np.array(self.dataset_w.variables['AERONET_nominal_wavelengths'][:])
        nwl = len(self.wl_data)

        hypstar_data = np.array(self.dataset_w.variables['mu_HYPSTAR_TO_AERONET_Lw'][:])
        aeronet_data = np.array(self.dataset_w.variables['mu_AERONET_Lw'][:])
        ndata = hypstar_data.shape[0]

        nspectra = ndata / nwl
        imu = 0
        for ispectra in range(int(nspectra)):
            file_out = os.path.join(path_plots, f'ComparisonSpectra_{ispectra + 1}.tif')
            imu_ini = imu
            imu_end = imu + nwl
            self.ydata_1 = aeronet_data[imu_ini:imu_end]
            self.ydata_2 = hypstar_data[imu_ini:imu_end]
            iday = int(self.dataset_w.variables['mu_day_id'][imu_ini])
            aseq = int(self.dataset_w.variables['mu_AERONET_sequence_id'][imu_ini])
            hseq = int(self.dataset_w.variables['mu_HYPSTAR_sequence_id'][imu_ini])
            atime = dt.utcfromtimestamp(float(self.dataset_w.variables['AERONET_time'][iday, aseq]))
            htime = dt.utcfromtimestamp(float(self.dataset_w.variables['HYPSTAR_time'][iday, hseq]))
            atimestr = atime.strftime('%Y-%m-%d %H:%M')
            htimestr = htime.strftime('%Y-%m-%d %H:%M')
            legend = [f'AERONET-OC({atimestr})', f'HYPSTAR({htimestr})']
            title = f'SPECTRA #{ispectra + 1}'
            ylabel = ''
            self.plot_spectra_impl(file_out, legend, ylabel,title)
            imu = imu_end

        self.close_file_w()

    def plot_spectra_impl(self, file_out, legend, ylabel, title):

        variables = list(self.spectra_stats.keys())
        self.ydata_1 = self.spectra_stats[variables[0]]['median']
        self.ydata_2 = self.spectra_stats[variables[1]]['median']

        ydata_1_min = self.spectra_stats[variables[0]]['p25']
        ydata_1_max = self.spectra_stats[variables[0]]['p75']
        ydata_2_min = self.spectra_stats[variables[1]]['p25']
        ydata_2_max = self.spectra_stats[variables[1]]['p75']

        xlabel = 'Wavelength (nm)'
        #ylabel = 'Lt'
        from MDB_reader.PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        pspectra.xdata = self.wl_data
        hline1 = pspectra.plot_single_line(self.ydata_1, 'red', 'solid', 1, 'o', 5)

        pspectra.plot_iqr_basic(ydata_1_min, ydata_1_max, 'red')

        hline2 = pspectra.plot_single_line(self.ydata_2, 'blue', 'solid', 1, 'o', 5)
        pspectra.plot_iqr_basic(ydata_2_min, ydata_2_max, 'blue')

        handles = [hline1[0], hline2[0]]

        pspectra.set_xticks(self.wl_data, self.wl_data, 90, 8)
        pspectra.set_yaxis_title(ylabel)
        pspectra.set_xaxis_title(xlabel)
        pspectra.set_title(title)

        pspectra.set_grid()
        pspectra.legend_options['loc'] = 'lower center'
        pspectra.legend_options['bbox_to_anchor'] = (0.1, -0.4, 0.80, 0.1)
        pspectra.legend_options['ncols'] = 2
        pspectra.set_legend_h(handles, legend)
        pspectra.set_tigth_layout()
        pspectra.save_fig(file_out)

    def plot_all_scatterplots_wl(self):
        self.dataset_w = Dataset(self.path_nc, 'r')
        wl_list = np.array(self.dataset_w.variables['AERONET_nominal_wavelengths'][:])
        self.close_file_w()
        xvariables  = ['mu_HYPSTAR_TO_AERONET_Li_mean','mu_HYPSTAR_TO_AERONET_Lt_mean','mu_HYPSTAR_TO_AERONET_Lw']
        yvariables = ['mu_AERONET_Li_mean', 'mu_AERONET_Lt_mean', 'mu_AERONET_Lw']

        # self.dataset_w = Dataset(self.path_nc, 'r')
        # array_1020 = np.array(self.dataset_w['mu_HYPSTAR_TO_AERONET_Lw'][:])
        # array_865  = np.array(self.dataset_w['mu_AERONET_Lt_mean'][:])
        # array_400 = np.array(self.dataset_w['mu_AERONET_Li_mean'][:])
        # wavelength = np.array(self.dataset_w['mu_wavelength'][:])
        # day_id = np.array(self.dataset_w['mu_day_id'][:])
        # check_invalid = np.zeros(array_1020.shape)
        # check_invalid[np.logical_and(wavelength==1020.0,np.logical_and(array_1020 > 0.2,array_1020<100000))] = 1
        # check_invalid[np.logical_and(wavelength == 865.0, np.logical_and(array_865 > 0.5, array_865 < 100000))] = 1
        # check_invalid[np.logical_and(wavelength == 400.0, np.logical_and(array_400 > 20, array_400 < 100000))] = 1
        # print(np.sum(check_invalid))
        # day_id_invalid = day_id[check_invalid==1]
        # print(day_id_invalid)
        # self.close_file_w()

        # for idx in range(len(xvariables)):
        #     for wl in wl_list:
        #         self.plot_scatterplot(wl,xvariables[idx],yvariables[idx],'mu_wavelength')

        ##global scatterplots
        # for idx in range(len(xvariables)):
        #     self.plot_scatterplot(None,xvariables[idx],yvariables[idx],'mu_wavelength')
        #     self.plot_scatterplot('rrsdensity', xvariables[idx], yvariables[idx], 'mu_wavelength')



    def get_wl_list(self, wlvariable):
        self.dataset_w = Dataset(self.path_nc, 'r')
        wl_list = np.unique(np.array(self.dataset_w.variables[wlvariable][:]))
        self.close_file_w()
        return wl_list

    def plot_rho_scatterplot(self):
        self.dataset_w = Dataset(self.path_nc, 'r')
        # hypstar_Rho = np.array(self.dataset_w.variables['HYPSTAR_Rho'][:])
        # aeronet_Rho = np.array(self.dataset_w.variables['AERONET_Rho'][:,:,0])
        ndata = self.dataset_w.variables['mu_day_id'].shape[0]
        nwl = self.dataset_w.variables['AERONET_nominal_wavelengths'].shape[0]
        nspectra = int(ndata / nwl)

        # print(ndata)
        # self.xdata = [0] * nspectra
        # self.ydata = [0] * nspectra


        print('NData:',ndata,'Nwl:',nwl,'NSpectra:',nspectra)
        self.xdata = []
        self.ydata = []

        array_mu_1020 =np.array(self.dataset_w['mu_HYPSTAR_TO_AERONET_Lw'][:])
        invalid_days = [21, 30, 32, 35, 63, 66, 80]


        imu = 0
        for ispectra in range(int(nspectra)):

            iday = int(self.dataset_w.variables['mu_day_id'][imu])
            ihs = int(self.dataset_w.variables['mu_HYPSTAR_sequence_id'][imu])
            ias = int(self.dataset_w.variables['mu_AERONET_sequence_id'][imu])

            imu_ini = imu
            imu_fin = imu + nwl
            if iday in invalid_days:
                imu = imu + nwl
                continue
            if np.count_nonzero(array_mu_1020[imu_ini:imu_fin]>100000)>0:
                imu = imu + nwl
                continue


            # self.xdata[ispectra] = float(self.dataset_w.variables['AERONET_Rho'][iday, ias, 0])
            # self.ydata[ispectra] = float(self.dataset_w.variables['HYPSTAR_rhof'][iday, ihs])

            self.xdata.append(float(self.dataset_w.variables['AERONET_Rho'][iday, ias, 0]))
            self.ydata.append(float(self.dataset_w.variables['HYPSTAR_rhof'][iday, ihs]))

            imu = imu + nwl


        print(len(self.xdata))

        self.close_file_w()

        self.plot_scatterplot('rho',None,None,None)

    def plot_scatterplot(self, wl, xvariable, yvariable,wlvariable):

        path_plots = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/PLOTS'
        marker = 'o'
        markersize = 20
        edgecolor = None
        linewidth = 0

        fontsizeaxis = 12

        if wl == 'rho':
            xlabel = 'AERONET-OC Rho'
            ylabel = 'HYPSTAR Rho'
            min_xy = 0.025
            max_xy = 0.030
            ticks = [0.025, 0.026, 0.027, 0.028, 0.029, 0.030]
        else:
            if xvariable.find('Li')>0:
                ylabel = r'AERONET-OC Li [μW/(cm$^2$·sr·nm)]'
                xlabel = r'HYPSTAR Li [μW/(cm$^2$·sr·nm)]'
            elif xvariable.find('Lt')>0:
                ylabel = r'AERONET-OC Lt [μW/(cm$^2$·sr·nm)]'
                xlabel = r'HYPSTAR Lt [μW/(cm$^2$·sr·nm)]'
            elif xvariable.find('Lw')>0:
                ylabel = r'AERONET-OC [μW/(cm$^2$·sr·nm)]'
                xlabel = r'HYPSTAR Lw [μW/(cm$^2$·sr·nm)]'
            min_xy = None
            max_xy = None
            ticks = None
            # if wl is None or wl < 600:
            #     min_xy = 0
            #     max_xy = 3
            #     ticks = [0, 0.5, 1, 1.5, 2, 2.5, 3]
            # if wl > 600:  # and wl<700:
            #     min_xy = 0
            #     max_xy = 0.6
            #     ticks = [0, 0.2, 0.4, 0.6]

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
        #self.set_data_scatterplot(wl)

        if wl!='rho':
            if wl=='rrsdensity':
                valid_array = self.set_data_scatterplot_spectral(None, xvariable, yvariable, wlvariable)
            else:
                valid_array = self.set_data_scatterplot_spectral(wl,xvariable,yvariable,wlvariable)

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
            file_out = os.path.join(path_plots, f'GlobalScatterplot_{yvariable}.tif')
            self.dataset_w = Dataset(self.path_nc, 'r')
            groupData = np.array(self.dataset_w.variables['mu_wavelength'][:])
            if valid_array is not None:
                groupData = groupData[valid_array==1]
            self.dataset_w.close()
            groupValues = np.unique(groupData)

            ngroup = len(groupValues)
            str_legend = []
            for g in groupValues:
                str_legend.append(f'{g:.2f}')

            for idx in range(ngroup):  # groupValues:
                g = groupValues[idx]
                #color = defaults.get_color_ref(g)
                color = defaults.get_color_wavelength(g)
                xhere = self.xdata[groupData == g]
                yhere = self.ydata[groupData == g]
                plot.plot_data(xhere, yhere, marker, markersize, color, edgecolor, linewidth)
            if min_xy is not None and max_xy is not None:
                plot.set_limits(min_xy, max_xy)
            if ticks is not None:
                plot.set_ticks(ticks, fontsizeaxis)
            plot.set_legend(str_legend)
            str_stats = self.get_str_stat_list(stat_list, None)
        else:
            if wl != 'rho' and wl != 'rrsdensity':
                wls = self.get_wl_str_from_wl(wl)
                file_out = os.path.join(path_plots, f'Wl_Scatterplot_{yvariable}_{wls}.tif')
            else:
                if wl == 'rho':
                    file_out = os.path.join(path_plots, f'Rho_Scatterplot_{yvariable}.tif')
                elif wl == 'rrsdensity':
                    file_out = os.path.join(path_plots, f'Global_Scatterplot_{yvariable}_density.tif')
            ##DENSITY
            xhere = np.asarray(self.xdata, dtype=np.float)
            yhere = np.asarray(self.ydata, dtype=np.float)
            # print(xhere.shape,' but valid',valid_array.shape)
            # if valid_array is not None:
            #     xhere = xhere[valid_array==1]
            #     yhere = yhere[valid_array==1]

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
            if wl != 'rho' and wl != 'rrsdensity':
                title = f'{wls} nm'
                plot.set_title(title)
            if min_xy is not None and max_xy is not None:
                plot.set_limits(min_xy, max_xy)
            if ticks is not None:
                plot.set_ticks(ticks, fontsizeaxis)

        plot.set_xaxis_title(xlabel)
        plot.set_yaxis_title(ylabel)

        # plot.plot_regress_line(self.xregress, self.yregress, 'k')
        plot.plot_identity_line()

        plot.plot_text(stats_xpos, stats_ypos, str_stats)
        plot.set_equal_apect()

        plot.save_fig(file_out)
        plot.close_plot()
