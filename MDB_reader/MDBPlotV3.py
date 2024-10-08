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

        if use_log_scale:
            valid_array = np.logical_and(sat_obs > 0, ref_obs > 0)
            nvalid = np.count_nonzero(valid_array)
            self.valid_stats['N'] = nvalid
            self.valid_stats['NMU'] = nvalid
            self.valid_stats['NGROUP'] = nvalid
            sat_obs = sat_obs[valid_array]
            ref_obs = ref_obs[valid_array]

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

        if use_log_scale:
            ##convert statistict to linear scale again
            stats_to_convert = ['RMSD', 'XAVG', 'YAVG', 'CRMSE', 'MAE']
            for stat in stats_to_convert:
                self.valid_stats[stat] = np.power(10, self.valid_stats[stat])
            bias_neg = self.valid_stats['BIAS'] < 0
            self.valid_stats['BIAS'] = np.power(10, np.abs(self.valid_stats['BIAS']))
            if bias_neg:
                self.valid_stats['BIAS'] = self.valid_stats['BIAS'] * (-1)

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
        list_virtual = poptions.get_list_virtual_flags()
        # print(list_figures)
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
                if options_figure['selectBy'] in list_virtual:
                    self.create_virtual_flag(poptions, options_figure['groupBy'])
                options_figure = self.check_gs_options_impl(options_figure, 'groupBy', 'groupType', 'groupValues')
            self.plot_from_options_impl(options_figure)

    def plot_from_options_impl(self, options_figure):
        if options_figure['type'] == 'scatterplot':
            self.plot_scatterplot_from_options(options_figure)
        if options_figure['type'] == 'spectraplot':
            self.plot_spectraplot_from_options(options_figure)
        if options_figure['type'] == 'timeseries':
            self.plot_time_series(options_figure)
        if options_figure['type'] == 'mapplot':
            self.plot_map_plot(options_figure)
        if options_figure['type'] == 'imageplot':
            self.plot_image(options_figure)

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
            file_out = f'{file_out[:-4]}_{index_mu}{file_out[-4:]}'
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
            file_out = f'{file_out[:-4]}_{index_mu}{file_out[-4:]}'
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
        print('??????????????????????????????????')
        print(instant_values.shape)
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
                            print(ivalue, '--->', ydata)
                            insitu_array[instant_values == ivalue] = ydata
                        if index == 1:
                            print(ivalue, '--->', ydata)
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

        print('??????????????????????????????????')
        print(instant_values.shape)
        print(pspectra.xdata.shape)
        print(insitu_array.shape)
        print(global_array.shape)
        print(regional_array.shape)
        print(style)
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
        print('aqui')
        print(options)
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
            print(self.valid_stats)
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

                    plot.plot_data(xhere, yhere, marker, markersize, z, None, 0)

                except:
                    print(f'[ERROR] Error creating density plot. Using default style')
                    plot.plot_data(xhere, yhere, marker, markersize, color, edgecolor, linewidth)

            else:
                plot.plot_data(xhere, yhere, marker, markersize, color, edgecolor, linewidth)

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
        plot.set_limits(min_xy, max_xy)

        # ticks
        if options['ticks'] is None:
            if not options['log_scale']:
                ticks = self.get_ticks_from_min_max_xy(min_xy, max_xy)
            else:
                min_t = int(np.log10(min_xy))
                max_t = int(np.log10(max_xy))
                ticks = [10 ** x for x in range(min_t, max_t + 1)]
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
        if options['log_scale']:
            if options['regression_line']:
                xr = np.array([self.xregress[0], self.xregress[-1], self.xregress[-2]])
                yr = np.array([self.yregress[0], self.yregress[-1], self.yregress[-2]])
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
            plot.axhere.title.set_size(options['fontsizetitle'])

        ##saving to file
        if not options['file_out'] is None:
            plot.save_fig(options['file_out'])
            plot.close_plot()

        return plot

    def plot_spectraplot_from_options(self, options_figure):

        if options_figure['type_rrs'] == 'ins':
            self.plot_insitu_spectraplots(options_figure)

        if options_figure['type_rrs'] == 'mu_comparison':
            index_mu = options_figure['index_mu']
            if index_mu == -1:
                for imu in range(self.mrfile.n_mu_total):
                    self.plot_mu_spectraplot(options_figure, imu)
            elif index_mu >= 0 and index_mu < self.mrfile.n_mu_total:
                self.plot_mu_spectraplot(options_figure, index_mu)

    def plot_mu_spectraplot(self, options_figure, index_mu):
        wl, insitu_spectra, sat_spectra, insitu_spectra_unc, sat_spectra_unc = self.mrfile.get_mu_spectra_insitu_and_sat(
            index_mu,options_figure['scale_factor'])



        if wl is None:
            return
        from PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        pspectra.xdata = wl
        if options_figure['wlticks'] is not None:
            wlt = options_figure['wlticks']
            wls = [f'{x}' for x in wlt]
        else:
            wlt = wl
            wls = self.mrfile.get_sat_wl_as_strlist(wl)

        pspectra.set_xticks(wlt, wls, 0, 12)
        if options_figure['type_rrs'] == 'mu_comparison':
            hline1 = pspectra.plot_single_line(insitu_spectra, 'red', 'solid', 1, '.', 0)
            hline2 = pspectra.plot_single_line(sat_spectra, 'blue', 'solid', 1, '.', 0)
            if insitu_spectra_unc is not None:
                insitu_min =  insitu_spectra.copy()
                insitu_max  = insitu_spectra.copy()
                insitu_min[~insitu_spectra_unc.mask] = insitu_spectra[~insitu_spectra_unc.mask] - insitu_spectra_unc[
                    ~insitu_spectra_unc.mask]
                insitu_max[~insitu_spectra_unc.mask] = insitu_spectra[~insitu_spectra_unc.mask] + insitu_spectra_unc[
                    ~insitu_spectra_unc.mask]
                pspectra.plot_iqr_basic(insitu_min,insitu_max,'red')

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
            pspectra.set_legend_h([hline1[0], hline2[0]], ['In situ Rrs', 'Satellite Rrs'])

        if options_figure['y_min'] is not None and options_figure['y_max'] is not None:
            pspectra.set_y_range(options_figure['y_min'],options_figure['y_max'])
        pspectra.set_grid()
        pspectra.set_tigth_layout()
        if not options_figure['file_out'] is None:
            file_out = options_figure['file_out'][:-4]
            file_out = f'{file_out}_{index_mu}.tif'
            # print(file_out)
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

        for idx in range(len(xarray)):
            if valid_all[idx] == 1:
                perc = 100 * ((yarray[idx] - xarray[idx]) / xarray[idx])
                perc = abs(perc)
                if perc > 1000:
                    valid_all[idx] = 0
                print(idx, ';', xarray[idx], ';', yarray[idx], ';', perc)

        # #smp wfr
        ##valid_all[27]=0

        if 'mu_valid' in self.mrfile.nc.variables:
            print('===============================================================================>')
            mu_valid = self.mrfile.nc.variables['mu_valid'][:]
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
                print(stat.upper(), val)

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

    def create_virtual_flag(self, poptions, vflag):
        if vflag in self.virtual_flags:
            return
        from FlagBuilder import FlagBuilder
        fbuilder = FlagBuilder(self.path_mdbr_file, None, poptions.options)
        flag_values, flag_names, array = fbuilder.create_flag_array(vflag, False)

        self.virtual_flags[vflag] = {
            'flag_array': array,
            'flag_values': flag_values,
            'flag_meanings': flag_names
        }

        # options_flag = fbuilder.get_options_dict(vflag)
        # print(options_flag)
        # potential_types = fbuilder.flag_options['typevirtual']['list_values']
        # #poptions = PlotOptions()
        # type_v = poptions.get_value_param(vflag,'typevirtual',None,'str',potential_types)
        # print(type_v)

        # print('aqui')
