import os

import pandas as pd
import pytz
from netCDF4 import Dataset
from PlotSpectra import PlotSpectra
from datetime import datetime as dt
import numpy as np
from matplotlib import pyplot as plt
import cartopy
from matplotlib.colors import LogNorm
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap


def get_smooth_array(array, type, window_size):
    array_out = array.copy()
    leg = int(window_size / 2)
    if type == 'RunningAverage':
        array = np.ma.filled(array, 0)
        for idx in range(len(array)):
            index_ini = idx - leg
            index_end = idx + leg + 1
            if index_ini < 0: index_ini = 0
            if index_end > len(array): index_end = len(array)
            array_here = array[index_ini:index_end]
            val = np.mean(array_here)

            array_out[idx] = val

    return array_out


def define_box_properties(plot_name, color_code, label):
    for k, v in plot_name.items():
        plt.setp(plot_name.get(k), color=color_code)
    # use plot function to draw a small line to name the legend.
    h, = plt.plot([], c=color_code, label=label)
    return h


def plot_multiple_bounding_boxes(file_nc, var_time, vars_data, dims_bar, time_ranges, x_axis, y_axis_title, y_ticks,
                                 plot_average, series, file_out):
    print('PLOTTING...')
    from matplotlib import pyplot as plt

    # time and ticks
    dataset = Dataset(file_nc)
    time = dataset.variables[var_time][:]
    time_obj = [dt.utcfromtimestamp(t) for t in time]
    year_array = np.array([obj.year for obj in time_obj])
    month_array = np.array([obj.month for obj in time_obj])
    jjj_array = np.array([int(obj.strftime('%j')) for obj in time_obj])

    select_array = np.ones(time.shape)
    if 'year' in time_ranges.keys():
        year_ini = time_ranges['year'][0]
        year_end = time_ranges['year'][1]
        select_array[year_array < year_ini] = 0
        select_array[year_array > year_end] = 0
    if 'jday' in time_ranges.keys():
        j_ini = time_ranges['jday'][0]
        j_end = time_ranges['jday'][1]
        select_array[jjj_array < j_ini] = 0
        select_array[jjj_array > j_end] = 0

    ##subset array using time_ranges, defining ntime (number of data points)
    if x_axis == 'year':
        time_here = year_array[select_array == 1]
    if x_axis == 'month':
        time_here = month_array[select_array == 1]
    if x_axis == 'era':
        era_info = get_era_info()
        time_here = np.zeros(time.shape)
        index_era = 0
        era_ticks = ['S', 'SMM', 'MM', 'VM', 'VMOA', 'VMOAB', 'OAB']
        for era in era_info:
            if era == 'All': continue
            sdate = era_info[era]['start_time']
            edate = era_info[era]['end_time']
            time_here[np.logical_and(time >= sdate, time <= edate)] = index_era
            index_era = index_era + 1
        time_here = time_here[select_array == 1]

    time_here_ini = int(np.min(time_here))
    time_here_end = int(np.max(time_here))
    n_time = (time_here_end - time_here_ini) + 1

    ##dimensions
    nseries = len(vars_data)
    width = 0.6
    intra = 0.1
    inter = 0.6
    if dims_bar is not None:
        width = dims_bar[0]
        intra = dims_bar[1]
        inter = dims_bar[2]
    all_pos_series = []

    width_series = round(nseries * (width + intra), 2)
    increm = round(width_series + inter, 2)
    for idx in range(nseries):
        start = inter + ((idx + 1) * (width + intra)) - (width / 2)
        start = round(start, 2)
        end = start + (increm * (n_time - 1))
        end = round(end, 2)
        pos_series = np.arange(start, end + (increm / 2), increm)
        max_value = end + (width / 2) + inter + intra
        all_pos_series.append(pos_series)

    ticks = [0] * n_time
    pos_ini = all_pos_series[0]
    pos_end = all_pos_series[-1]

    for idx in range(n_time):
        ticks[idx] = pos_ini[idx] + ((pos_end[idx] - pos_ini[idx]) / 2)

    max_value = round(max_value, 2)

    ##getting data and plotting
    handles = []
    for idx, var in enumerate(vars_data):
        series_name = list(series.keys())[idx]
        array = dataset.variables[var][:]
        array_here = array[select_array == 1]
        data_boxplot = []
        data_average = []
        for there in range(time_here_ini, time_here_end + 1):
            data_here = array_here[time_here == there]
            data_here = data_here[data_here.mask == False]
            data_here = np.array(data_here).flatten()
            if plot_average:
                data_average.append(np.mean(data_here))
            data_boxplot.append(data_here)
        positions = all_pos_series[idx]

        bbox = plt.boxplot(data_boxplot, positions=positions, widths=width,
                           flierprops=dict(markersize=0, markeredgecolor='none'))
        hlegend = define_box_properties(bbox, series[series_name], series_name)
        handles.append(hlegend)
        if plot_average:
            for itime in range(n_time):
                plt.plot(positions[itime], data_average[itime], color=series[series_name], linewidth=0, marker='o',
                         markersize=8, mec='black', mew=1)

        if x_axis == 'era':
            plt.xticks(ticks, era_ticks)
        else:
            plt.xticks(ticks, np.arange(time_here_ini, time_here_end + 1), rotation=90)
    plt.xlim(0, max_value)
    plt.xlabel(x_axis.capitalize())

    plt.ylabel(y_axis_title)
    if y_ticks is not None:
        plt.yticks(y_ticks)
    plt.grid(which='major', color='lightgray', linestyle='--', axis='y')

    plt.legend(handles, list(series.keys()), loc='lower center', ncols=nseries, bbox_to_anchor=(0.5, -0.2))
    plt.tight_layout()
    plt.savefig(file_out, dpi=300)
    plt.close()

    dataset.close()


def plot_bounding_boxes(file_nc, var_time, var_data, time_ranges, x_axis, y_axis_title, y_ticks, plot_global,
                        plot_average, file_out):
    print('PLOTTING...')
    dataset = Dataset(file_nc)
    time = dataset.variables[var_time][:]
    array = dataset.variables[var_data][:]
    dataset.close()
    time_obj = [dt.utcfromtimestamp(t) for t in time]
    year_array = np.array([obj.year for obj in time_obj])
    month_array = np.array([obj.month for obj in time_obj])
    jjj_array = np.array([int(obj.strftime('%j')) for obj in time_obj])

    select_array = np.ones(time.shape)
    if 'year' in time_ranges.keys():
        year_ini = time_ranges['year'][0]
        year_end = time_ranges['year'][1]
        select_array[year_array < year_ini] = 0
        select_array[year_array > year_end] = 0
    if 'jday' in time_ranges.keys():
        j_ini = time_ranges['jday'][0]
        j_end = time_ranges['jday'][1]
        select_array[jjj_array < j_ini] = 0
        select_array[jjj_array > j_end] = 0

    ##subset array using time_ranges
    if x_axis == 'year':
        time_here = year_array[select_array == 1]
    if x_axis == 'month':
        time_here = month_array[select_array == 1]
    if x_axis == 'era':
        era_info = get_era_info()
        time_here = np.zeros(time.shape)
        index_era = 0
        era_ticks = ['S', 'SMM', 'MM', 'VM', 'VMOA', 'VMOAB', 'OAB']
        for era in era_info:
            if era == 'All': continue
            sdate = era_info[era]['start_time']
            edate = era_info[era]['end_time']
            time_here[np.logical_and(time >= sdate, time <= edate)] = index_era
            index_era = index_era + 1
        time_here = time_here[select_array == 1]

    array_here = array[select_array == 1]

    time_here_ini = int(np.min(time_here))
    time_here_end = int(np.max(time_here))
    n_time = (time_here_end - time_here_ini) + 1

    data_boxplot = []
    data_average = []
    for there in range(time_here_ini, time_here_end + 1):
        data_here = array_here[time_here == there]
        data_here = data_here[data_here.mask == False]
        data_here = np.array(data_here).flatten()
        if plot_average:
            data_average.append(np.mean(data_here))
        data_boxplot.append(data_here)
    from matplotlib import pyplot as plt

    plt.boxplot(data_boxplot)

    if plot_average:
        for x in range(1, n_time + 1):
            plt.plot(x, data_average[x - 1], color='white', linewidth=0, marker='o', markersize=10, mec='black', mew=1)

    if x_axis == 'era':
        plt.xticks(np.arange(1, n_time + 1), era_ticks)
    else:
        plt.xticks(np.arange(1, n_time + 1), np.arange(time_here_ini, time_here_end + 1), rotation=90)
    plt.xlabel(x_axis.capitalize())
    plt.ylabel(y_axis_title)
    if y_ticks is not None:
        plt.yticks(y_ticks)
    plt.grid(which='major', color='lightgray', linestyle='--', axis='y')

    if plot_global:
        array_here = array_here[array_here.mask == False]
        array_here = np.array(array_here).flatten()
        global_median = np.median(array_here)
        global_p25 = np.percentile(array_here, 25)
        global_p75 = np.percentile(array_here, 75)
        plt.plot([0.75, n_time + 0.25], [global_median, global_median], color='black', linestyle='-', linewidth=1,
                 marker=None)
        plt.plot([0.75, n_time + 0.25], [global_p25, global_p25], color='black', linestyle='--', linewidth=1,
                 marker=None)
        plt.plot([0.75, n_time + 0.25], [global_p75, global_p75], color='black', linestyle='--', linewidth=1,
                 marker=None)

    plt.tight_layout()
    plt.savefig(file_out, dpi=300)
    plt.close()


def plot_time_series_smoothed(file_nc, file_out, file_csv_year, file_csv_sensor, var_data, time_ranges, plot_global,
                              y_axis_title, y_ticks):
    print('PLOTTING...')
    dataset = Dataset(file_nc)
    time = dataset.variables['time'][:]
    array = dataset.variables[var_data][:]
    pspectra = PlotSpectra()
    pspectra.close_plot()

    time_obj = [dt.utcfromtimestamp(t) for t in time]
    year_array = np.array([obj.year for obj in time_obj])
    jjj_array = np.array([int(obj.strftime('%j')) for obj in time_obj])
    month_array = np.array([obj.month for obj in time_obj])

    select_array = np.ones(time.shape)
    if 'year' in time_ranges.keys():
        year_ini = time_ranges['year'][0]
        year_end = time_ranges['year'][1]
        select_array[year_array < year_ini] = 0
        select_array[year_array > year_end] = 0
    if 'jday' in time_ranges.keys():
        j_ini = time_ranges['jday'][0]
        j_end = time_ranges['jday'][1]
        select_array[jjj_array < j_ini] = 0
        select_array[jjj_array > j_end] = 0

    ##subset array using time_ranges and computing running average
    time_here = time[select_array == 1]
    array_here = array[select_array == 1]
    array_smooth = get_smooth_array(array, 'RunningAverage', 30)
    array_smooth_here = array_smooth[select_array == 1]

    ##getting year ticks
    year_ticks = []
    time_ticks = []
    prev = -1
    for t in time_here:
        tobj = dt.utcfromtimestamp(t)
        if tobj.year != prev:
            time_ticks.append(t)
            year_ticks.append(tobj.strftime('%Y'))
            prev = tobj.year
    time_ticks.append(time_here[-1])
    year_ticks.append('')

    ##plotting median by year
    if file_csv_year is not None:
        df = pd.read_csv(file_csv_year, sep=';')
        median_values = df[f'{var_data}_median']
        year_values = df['year']
        for idx in range(len(time_ticks) - 1):
            year_here = dt.utcfromtimestamp(time_ticks[idx]).year
            if year_values[idx] == year_here:
                yval = median_values[idx]
                xini = time_ticks[idx]
                xfin = time_ticks[idx + 1]
                pspectra.plot_single_data([xini, xfin], [yval, yval],
                                          {'color': 'black', 'linestyle': '-', 'linewidth': 1, 'marker': None,
                                           'markersize': 12})

    if file_csv_sensor is not None:
        df = pd.read_csv(file_csv_sensor, sep=';')
        median_values = df[f'{var_data}_median']
        eras_df = df['era']
        eras = get_era_info()
        style_era = {'color': 'red', 'linestyle': '-', 'linewidth': 1, 'marker': None, 'markersize': 12}
        for idx in range(len(eras_df)):
            era = eras_df[idx]
            print(eras[era])
            sdate = eras[era]['start_time']
            edate = eras[era]['end_time']
            x_ini = np.where(time_here == sdate)[0]
            if len(x_ini) == 0:
                x_ini = 0
            else:
                x_ini = x_ini[0]
            x_end = np.where(time_here == edate)[0]
            if len(x_end) == 0:
                x_end = len(time_here) - 1
            else:
                x_end = x_end[0]
            yval = median_values[idx]
            style_era['color'] = eras[era]['color']
            pspectra.plot_single_data([time_here[x_ini], time_here[x_end]], [yval, yval], style_era)

    ##plotting global medians
    if plot_global:
        array_here = np.ma.filled(array_here, 0)
        global_median = np.median(array_here[array_here != 0])
        global_p25 = np.percentile(array_here[array_here != 0], 25)
        global_p75 = np.percentile(array_here[array_here != 0], 75)

        pspectra.plot_single_data([time_here[0], time_here[-1]], [global_median, global_median],
                                  {'color': 'black', 'linestyle': '-', 'linewidth': 2, 'marker': None,
                                   'markersize': 12})
        pspectra.plot_single_data([time_here[0], time_here[-1]], [global_p25, global_p25],
                                  {'color': 'black', 'linestyle': '--', 'linewidth': 1, 'marker': None,
                                   'markersize': 12})
        pspectra.plot_single_data([time_here[0], time_here[-1]], [global_p75, global_p75],
                                  {'color': 'black', 'linestyle': '--', 'linewidth': 1, 'marker': None,
                                   'markersize': 12})

    ##plotting main data (running average)
    pspectra.xdata = time_here
    pspectra.plot_data(array_smooth_here,
                       {'color': 'blue', 'linestyle': '-', 'linewidth': 0.5, 'marker': None, 'markersize': 12})
    pspectra.set_xticks_time(time_ticks, year_ticks, 90, None)
    pspectra.set_grid_horizontal()

    ##ticks and titles
    if y_ticks is not None:
        print(y_ticks)
        pspectra.set_yticks(y_ticks, None, None, None)
    pspectra.set_yaxis_title_f(y_axis_title, 10)
    pspectra.set_xaxis_title_f('Year', 10)
    pspectra.set_tigth_layout()
    # saving
    pspectra.save_plot(file_out)
    pspectra.close_plot()


def start_full_figure():
    # start figure and axes

    fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
    fig.set_figwidth(15)
    fig.set_figheight(15)

    # coastlines
    ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black', linewidth=0.5)
    # ax.coastlines(resolution='10m')

    # grid lines
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=0.5, linestyle='dotted')
    gl.xlocator = mticker.FixedLocator([5, 10, 15, 20, 25, 30])
    gl.ylocator = mticker.FixedLocator([55, 57.5, 60, 62.5, 65])
    gl.right_labels = False
    gl.left_labels = True
    gl.bottom_labels = True
    gl.top_labels = False
    gl.xlabel_style = {'size': 15}
    gl.ylabel_style = {'size': 15}

    return fig, ax


def plot_mask_clara(file_mask, dir_out):
    print('PLOTTING MASK...')
    dataset = Dataset(file_mask)
    land_mask_clara = dataset.variables['Land_Mask_CFC'][:]
    lat_clara = dataset.variables['lat_cfc'][:]
    lon_clara = dataset.variables['lon_cfc'][:]
    land_mask_cci = dataset.variables['Land_Mask'][:]
    lat_cci = dataset.variables['lat'][:]
    lon_cci = dataset.variables['lon'][:]
    n_water_map = dataset.variables['NTotal_Water_Map_CFC'][:]
    cfc_mask = dataset.variables['CFC_Mask'][:]
    dataset.close()
    print(len(land_mask_cci[land_mask_cci==0]))
    print(len(cfc_mask[cfc_mask==0]))


    # CLARA MASK
    # lon_clara = lon_clara[32:105]
    # map = land_mask_clara[:, 32:105]
    # map = np.ma.masked_equal(map, 1)
    # fig, ax = start_full_figure()
    # ax.pcolormesh(lon_clara, lat_clara, map, cmap=mpl.colormaps[
    #     'bwr'])  # , vmin=vmin, vmax=vmax,cmap=mpl.colormaps['jet'])  # , norm=LogNorm(vmin=0.1, vmax=100))
    # file_out = os.path.join(dir_out, 'LandMaskClara.tif')
    # fig.savefig(file_out, dpi=300, bbox_inches='tight')
    # plt.close(fig)

    # CCI MASK
    # lon_cci = lon_cci[210:1186]
    # map = land_mask_cci[:, 210:1186]
    # map = np.ma.masked_equal(map, 1)
    # fig, ax = start_full_figure()
    # ax.pcolormesh(lon_cci, lat_cci, map, cmap=mpl.colormaps[
    #     'bwr'])  # , vmin=vmin, vmax=vmax,cmap=mpl.colormaps['jet'])  # , norm=LogNorm(vmin=0.1, vmax=100))
    # file_out = os.path.join(dir_out, 'LandMaskCCI.tif')
    # fig.savefig(file_out, dpi=300, bbox_inches='tight')
    # plt.close(fig)

    ##NWATER CCI CIN CLARA
    # lon_clara = lon_clara[32:105]
    # map = n_water_map[:, 32:105]
    # fig, ax = start_full_figure()
    # from matplotlib.colors import BoundaryNorm
    # from matplotlib.colors import ListedColormap
    # bnorm = BoundaryNorm([308,322,330,345,500], ncolors=4)
    # colors = ['cyan', 'blue','green','magenta']
    # cmap = ListedColormap(['cyan', 'blue','green','magenta'])
    # map[map==0] = np.ma.masked
    # ax.pcolormesh(lon_clara, lat_clara, map,norm = bnorm,cmap=cmap)
    # ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black', linewidth=0.5)
    # file_out = os.path.join(dir_out, 'NWaterCCIInClara.tif')
    # handles = []
    # from matplotlib.patches import Patch
    # for color in colors:
    #     handles.append(Patch(facecolor=color,edgecolor='black'))
    # plt.legend(handles,['308 OC-CCI pixels','322 OC-CCI pixels','330 OC-CCI pixels','345 OC-CCI pixels'],fontsize=15)
    # fig.savefig(file_out, dpi=300, bbox_inches='tight')
    # plt.close(fig)




# h = ax.pcolormesh(lon_array, lat_array, map, cmap=mpl.colormaps['jet'])
# cbar = fig.colorbar(h, cax=None, ax=ax, use_gridspec=True, fraction=0.03, format="$%.2f$")
# cbar.ax.tick_params(labelsize=15)

# cbar.set_label(label=f'Coverage ({units})', size=15)
# title = f'Chlorophyll a concentration - OLD ({units})'
# if dateherestr is not None:
#     title = f'{title} - {dateherestr}'
# ax.set_title(title, fontsize=20)
# file_chla_old = os.path.join(os.path.dirname(output_path), f'Img_Chla_OLD_{dateherestr}.png')


def plot_maps_clara(file_average, type, use_zeros_no_data):
    dir_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/PLOTS'
    if use_zeros_no_data:
        file_out = os.path.join(dir_out, f'Map_{type}_zerosinnodata.tif')
    else:
        file_out = os.path.join(dir_out, f'Map_{type}.tif')

    dataset = Dataset(file_average)
    array_map = dataset.variables['daily_cloud_free_map'][:]
    time = dataset.variables['time'][:]
    map = dataset.variables['Land_Mask_CFC'][:]
    lat_array = dataset.variables['lat_cfc'][:]
    lon_array = dataset.variables['lon_cfc'][:]
    dataset.close()

    ##DATA SELECTION
    time_obj = [dt.utcfromtimestamp(t) for t in time]
    year_array = np.array([obj.year for obj in time_obj])
    jjj_array = np.array([int(obj.strftime('%j')) for obj in time_obj])
    select_array = np.ones(time.shape)
    year_ini = 1998
    year_end = 2020
    select_array[year_array < year_ini] = 0
    select_array[year_array > year_end] = 0
    if type.endswith('_Summer'):
        j_ini = 161
        j_end = 270
        select_array[jjj_array < j_ini] = 0
        select_array[jjj_array > j_end] = 0

    array_map = array_map[select_array == 1, :]

    if type.startswith('PercentCoverageMap'):
        if use_zeros_no_data:
            array_map = np.ma.filled(array_map, 0)
        array_coverage = np.ma.mean(array_map, axis=0)
        if type.endswith('_Summer'):
            if use_zeros_no_data:
                vmin = 35
                vmax = 50
            else:
                vmin = 35
                vmax = 50
        else:
            if use_zeros_no_data:
                vmin = 17.5
                vmax = 32.5
            else:
                vmin = 35
                vmax = 45
        units = '%'

    if type.startswith('TotalCoverageMap'):
        array_coverage = np.ma.sum(array_map, axis=0) / 100
        vmin = None
        vmax = None
        units = 'pixels'

    ##CREATING MAP
    rows, columns = np.where(map == 0)
    map[rows, columns] = array_coverage[:]
    map = np.ma.masked_equal(map, 1)

    # SUBSETTING MAP
    lon_array = lon_array[32:105]
    map = map[:, 32:105]

    fig, ax = start_full_figure()

    h = ax.pcolormesh(lon_array, lat_array, map, vmin=vmin, vmax=vmax,
                      cmap=mpl.colormaps['jet'])  # , norm=LogNorm(vmin=0.1, vmax=100))
    # h = ax.pcolormesh(lon_array, lat_array, map, cmap=mpl.colormaps['jet'])
    cbar = fig.colorbar(h, cax=None, ax=ax, use_gridspec=True, fraction=0.03, format="$%.2f$")
    cbar.ax.tick_params(labelsize=15)

    cbar.set_label(label=f'Coverage ({units})', size=15)
    # title = f'Chlorophyll a concentration - OLD ({units})'
    # if dateherestr is not None:
    #     title = f'{title} - {dateherestr}'
    # ax.set_title(title, fontsize=20)
    # file_chla_old = os.path.join(os.path.dirname(output_path), f'Img_Chla_OLD_{dateherestr}.png')
    fig.savefig(file_out, dpi=300, bbox_inches='tight')


def plot_bars(xdata, ydata, x_title, y_title, range_y, file_out):
    pspectra = PlotSpectra()
    pspectra.xdata = xdata
    if range_y is not None:
        pspectra.set_y_range(range_y[0], range_y[1])
    pspectra.set_xticks(xdata, xdata, 90, None)
    pspectra.set_xaxis_title(x_title)
    if y_title is None:
        y_title = '# cloud-free pixels'
    pspectra.set_yaxis_title(y_title)
    pspectra.set_grid_horizontal()
    median = np.ma.median(ydata)
    p25 = np.percentile(ydata, 25)
    p75 = np.percentile(ydata, 75)
    pspectra.plot_single_bar_series(ydata, 'blue', 0.8, 0, 1)
    style_median = {'color': 'red', 'linestyle': '-', 'linewidth': 1, 'marker': None, 'markersize': 12}
    style_iqr = {'color': 'red', 'linestyle': '--', 'linewidth': 0.5, 'marker': None, 'markersize': 12}
    pspectra.plot_single_data([xdata[0], xdata[-1]], [median, median], style_median)
    pspectra.plot_single_data([xdata[0], xdata[-1]], [p25, p25], style_iqr)
    pspectra.plot_single_data([xdata[0], xdata[-1]], [p75, p75], style_iqr)
    pspectra.set_tigth_layout()
    pspectra.save_plot(file_out)
    pspectra.close_plot()


def plot_stacked_bars(xdata, ydata, x_title, y_title, range_y, ydata_secondary, y_title_secondary, file_out):
    # ydata series  in growing order
    pspectra = PlotSpectra()
    pspectra.xdata = xdata
    if range_y is not None:
        pspectra.set_y_range(range_y[0], range_y[1])
    pspectra.set_xticks(xdata, xdata, 90, None)
    pspectra.set_xaxis_title(x_title)
    pspectra.set_yaxis_title(y_title)
    pspectra.set_grid_horizontal()
    nseries = ydata.shape[1]
    import MDBPlotDefaults
    colors = MDBPlotDefaults.get_color_list(nseries)

    for iseries in range(nseries):

        if iseries == 0:
            ydata_here = ydata[:, iseries]
            pspectra.plot_single_bar_series(ydata_here, colors[iseries], 0.8, 0, 1)
        else:
            ydata_here = ydata[:, iseries]
            ydata_prev = ydata[:, iseries - 1]
            ydata_here = ydata_here - ydata_prev

            pspectra.plot_bar_series_bottom(ydata_here, colors[iseries], 0.8, ydata_prev, 1)

    if ydata_secondary is not None:
        if len(y_title_secondary) == 0:
            pspectra.plot_single_line(ydata_secondary, 'grey', '--', 0.5, 'o', 7)
        else:
            pspectra.plot_single_line_secondary_axis(ydata_secondary, y_title_secondary, 'black', '--', 0.5, 'o', 7,
                                                     'k', 1, 'white')

    pspectra.set_tigth_layout()

    pspectra.save_plot(file_out)
    pspectra.close_plot()


def plot_lines(file_csv, time_var, data_var, is_percent, y_axis_title, file_out):
    max_value = 80 if is_percent else 515

    df = pd.read_csv(file_csv, sep=';')

    xdata = df[time_var].to_numpy()
    ymedian = df[f'{data_var}_median'].to_numpy()
    ymedian[np.isnan(ymedian)] = 0
    if not is_percent:
        ymedian = ymedian / 100
    yp75 = df[f'{data_var}_p75'].to_numpy()
    yp75[np.isnan(yp75)] = 0
    if not is_percent:
        yp75 = yp75 / 100
    yp25 = df[f'{data_var}_p25'].to_numpy()
    yp25[np.isnan(yp25)] = 0
    if not is_percent:
        yp25 = yp25 / 100

    ymedian_smooth = get_smooth_array(ymedian, 'RunningAverage', 30)
    yp25_smooth = get_smooth_array(yp25, 'RunningAverage', 30)
    yp75_smooth = get_smooth_array(yp75, 'RunningAverage', 30)

    pspectra = PlotSpectra()
    pspectra.xdata = xdata

    pspectra.plot_data(ymedian_smooth,
                       {'color': 'blue', 'linestyle': '-', 'linewidth': 2, 'marker': None, 'markersize': 12})
    pspectra.plot_data(yp25_smooth,
                       {'color': 'blue', 'linestyle': '--', 'linewidth': 1, 'marker': None, 'markersize': 12})
    pspectra.plot_data(yp75_smooth,
                       {'color': 'blue', 'linestyle': '--', 'linewidth': 1, 'marker': None, 'markersize': 12})

    pspectra.set_tigth_layout()
    x_axis_title = 'Julian Day'
    pspectra.set_xaxis_title(x_axis_title)
    pspectra.set_yaxis_title(y_axis_title)
    pspectra.set_grid_horizontal()
    pspectra.plot_single_data([161, 161], [0, max_value],
                              {'color': 'black', 'linestyle': ':', 'linewidth': 1, 'marker': None, 'markersize': 12})
    pspectra.plot_single_data([270, 270], [0, max_value],
                              {'color': 'black', 'linestyle': ':', 'linewidth': 1, 'marker': None, 'markersize': 12})
    pspectra.set_tigth_layout()
    pspectra.save_plot(file_out)
    pspectra.close_plot()


def plot_scatter_plot(xdata, ydata, x_label, y_label, yx_range, groups, plot_lines, file_out):
    from PlotScatter import PlotScatter
    from pylr2 import regress2
    from scipy import stats
    from MDBPlotV3 import MDBPlot
    mplot = MDBPlot(None)
    pscatter = PlotScatter()

    if groups is not None:
        garray = groups['array']
        gvalues = groups['values']
        gmeanings = groups['meanings']
        colors = groups['colors']
        print(gvalues)
        print(gmeanings)
        print(colors)
        ngroups = len(gvalues)

        for idx in range(ngroups):
            gvalue = gvalues[idx]
            xdata_here = xdata[garray == gvalue]
            ydata_here = ydata[garray == gvalue]
            h = pscatter.plot_data(xdata_here, ydata_here, 'o', 20, colors[idx], 'grey', 0.5)
            # handles.append(h[0])
        # pscatter.set_legend(gmeanings)

    else:
        pscatter.plot_data(xdata, ydata, 'o', 12, 'blue', 'blue', 0)

    pscatter.set_xaxis_title(x_label)
    pscatter.set_yaxis_title(y_label)
    if plot_lines:
        pscatter.plot_identity_line()

    ##regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata, ydata)
    results = regress2(np.array(xdata, dtype=np.float64), np.array(ydata, dtype=np.float64),
                       _method_type_2="reduced major axis")
    slope_II = results['slope']
    intercept_II = results['intercept']
    if plot_lines:
        xregress, yregress = mplot.get_regression_line(xdata, ydata, slope_II, intercept_II, None, None)
        pscatter.plot_regress_line_options(xregress, yregress, 'black', '-', 1)

    str_equation = f'y={slope_II:.2f}x + {intercept_II:.2f}'
    str_r2 = r'R$^2$'
    r2 = r_value * r_value
    str_r2 = f'{str_r2}={r2:.2f}'
    pscatter.plot_text(0.03, 0.94, str_equation)
    pscatter.plot_text(0.03, 0.87, str_r2)
    if yx_range is not None:
        pscatter.set_limits_X(yx_range[0], yx_range[1])
        pscatter.set_limits_Y(yx_range[0], yx_range[1])

    # pscatter.set_equal_apect()
    pscatter.save_fig(file_out)
    pscatter.close_plot()


def plot_sensor_coverage_versus_clould_free(file_average, increm):
    dataset = Dataset(file_average)
    f_coverage = np.array(dataset.variables['f_coverage'][:])
    f_cloud_free = np.array(dataset.variables['f_cloud_free'][:])
    f_time = np.array(dataset.variables['f_time'][:])
    dataset.close()
    eras = {
        'All': {
            'start_time': dt(1997, 9, 4).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2024, 10, 31).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'black'
        },
        'SeaWiFS': {
            'start_time': dt(1997, 9, 4).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2002, 4, 28).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'blue'
        },
        'SeaWiFS_MERIS_MODIS': {
            'start_time': dt(2002, 4, 29).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2010, 12, 11).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'cyan'
        },
        'MERIS_MODIS': {
            'start_time': dt(2010, 12, 12).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2012, 4, 8).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'lime'
        },
        'VIIRS_MODIS': {
            'start_time': dt(2012, 4, 9).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2016, 4, 24).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'green'
        },
        'VIIRS_MODIS_OLCIA': {
            'start_time': dt(2016, 4, 25).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2018, 5, 13).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'salmon'
        },
        'VIIRS_MODIS_OLCIAB': {
            'start_time': dt(2018, 5, 14).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2019, 12, 31).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'red'
        },
        'OLCIAB': {
            'start_time': dt(2020, 1, 1).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2024, 10, 31).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'magenta'
        }
    }

    dir_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/PLOTS'
    for era in eras:
        print(f'Working with era: {era}')
        start_time = eras[era]['start_time']
        end_time = eras[era]['end_time']
        valid_era = np.logical_and(f_time >= start_time, f_time <= end_time)
        f_coverage_here = f_coverage[valid_era]
        f_cloud_free_here = f_cloud_free[valid_era]
        file_out = os.path.join(dir_out, f'SensorCoverage_vs_CloudFree_{era}.tif')
        median_f_coverage, popt = plot_sensor_coverage_versus_clould_free_impl(f_coverage_here, f_cloud_free_here,
                                                                               increm, None, file_out)
        eras[era]['median'] = median_f_coverage
        eras[era]['popt'] = popt

    file_out = os.path.join(dir_out, 'SensorCoverage_vs_CloudFree_Multiple.tif')
    plot_sensor_coverage_versus_clould_free_multiple_era(eras, increm, file_out)
    file_out = os.path.join(dir_out, 'SensorCoverage_vs_CloudFree_Multiple_Equations.tif')
    plot_sensor_coverage_versus_clould_free_multiple_era_eq(eras, file_out)


def plot_sensor_coverage_versus_clould_free_multiple_era_eq(eras, file_out):
    x_values = np.arange(0, 100.01, 0.01)

    pspectra = PlotSpectra()
    pspectra.xdata = x_values
    style = {'color': 'blue', 'linestyle': '-', 'linewidth': 1, 'marker': None, 'markersize': 3}

    legend_str = []
    legend_hanles = []
    for era in eras:
        if era == 'All':
            continue
        legend_str.append(era)
        popt = eras[era]['popt']
        values = sigmoid(x_values, *popt)

        style['color'] = eras[era]['color']

        hhere = pspectra.plot_data(values, style)
        legend_hanles.append(hhere[0])

    style = {'color': 'k', 'linestyle': '--', 'linewidth': 1, 'marker': None, 'markersize': 2}
    pspectra.plot_data(x_values, style)
    pspectra.legend_options['loc'] = 'lower center'
    pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.65)
    pspectra.legend_options['ncols'] = 2
    # pspectra.set_legend_h(legend_hanles,legend_str)
    pspectra.set_xaxis_title('Clould-free fraction (%)')
    pspectra.set_yaxis_title('Sensor coverage (%)')
    pspectra.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], None, None, None)
    pspectra.set_yticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], None, None, None)
    pspectra.set_grid()
    pspectra.set_tigth_layout()
    pspectra.save_plot(file_out)
    pspectra.close_plot()


def plot_sensor_coverage_versus_clould_free_multiple_era(eras, increm, file_out):
    porc_min = np.arange(0, 101, increm)
    porc_values = porc_min + (increm / 2)
    porc_values[-1] = 100
    pspectra = PlotSpectra()
    pspectra.xdata = porc_values
    style = {'color': 'blue', 'linestyle': '-', 'linewidth': 1, 'marker': 'o', 'markersize': 3}
    legend_str = []
    legend_hanles = []
    for era in eras:
        if era == 'All':
            continue
        legend_str.append(era)
        values = eras[era]['median']
        style['color'] = eras[era]['color']
        hhere = pspectra.plot_data(values, style)
        legend_hanles.append(hhere[0])
    style = {'color': 'k', 'linestyle': '--', 'linewidth': 1, 'marker': None, 'markersize': 2}
    pspectra.plot_data(porc_values, style)
    pspectra.legend_options['loc'] = 'lower center'
    pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.65)
    pspectra.legend_options['ncols'] = 2
    # pspectra.set_legend_h(legend_hanles,legend_str)
    pspectra.set_xaxis_title('Clould-free fraction (%)')
    pspectra.set_yaxis_title('Sensor coverage (%)')
    pspectra.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], None, None, None)
    pspectra.set_yticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], None, None, None)
    pspectra.set_grid()
    pspectra.set_tigth_layout()
    pspectra.save_plot(file_out)
    pspectra.close_plot()


def plot_sensor_coverage_versus_clould_free_impl(f_coverage, f_cloud_free, increm, y_label, file_out):
    porc_min = np.arange(0, 101, increm)
    nporc = len(porc_min)

    median_f_coverage = np.zeros((nporc,))
    avg_f_coverage = np.zeros((nporc,))
    p25_f_coverage = np.zeros((nporc,))
    p75_f_coverage = np.zeros((nporc,))
    for iporc in range(nporc):
        min_value = porc_min[iporc]
        max_value = porc_min[iporc] + increm
        values = f_coverage[np.logical_and(f_cloud_free >= min_value, f_cloud_free < max_value)]
        median_f_coverage[iporc] = np.median(values)
        avg_f_coverage[iporc] = np.mean(values)
        p25_f_coverage[iporc] = np.percentile(values, 25)
        p75_f_coverage[iporc] = np.percentile(values, 75)

    porc_values = porc_min + (increm / 2)
    porc_values[-1] = 100
    pspectra = PlotSpectra()
    pspectra.xdata = porc_values
    style = {'color': 'blue', 'linestyle': '-', 'linewidth': 1, 'marker': 'o', 'markersize': 3}
    pspectra.plot_data(median_f_coverage, style)
    style = {'color': 'red', 'linestyle': '-', 'linewidth': 0, 'marker': 'o', 'markersize': 3}
    pspectra.plot_data(avg_f_coverage, style)

    if y_label != 'Bloom coverage(%) - Image':
        style = {'color': 'k', 'linestyle': '--', 'linewidth': 1, 'marker': None, 'markersize': 2}
        pspectra.plot_data(porc_values, style)

    style = {'color': 'blue', 'linestyle': '--', 'linewidth': 0.5, 'marker': None, 'markersize': 3}
    pspectra.plot_data(p25_f_coverage, style)
    pspectra.plot_data(p75_f_coverage, style)
    # pspectra.plot_iqr_basic(p25_f_coverage,p75_f_coverage,'blue')

    pspectra.set_xaxis_title('Clould-free fraction (%)')
    if y_label is None:
        y_label = 'Sensor coverage (%)'

    pspectra.set_yaxis_title(y_label)
    pspectra.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], None, None, None)
    pspectra.set_yticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], None, None, None)
    if y_label == 'Bloom coverage(%) - Image':
        pspectra.set_y_range(0, 25)
        pspectra.set_yticks([0, 5, 10, 15, 20, 25], None, None, None)
    pspectra.set_grid()
    pspectra.set_tigth_layout()
    pspectra.save_plot(file_out)
    pspectra.close_plot()

    fit_curve = False
    popt = None
    if fit_curve:
        from scipy.optimize import curve_fit
        from scipy import stats
        p0 = [np.median(f_coverage), 1]
        popt, pcov = curve_fit(sigmoid, f_cloud_free, f_coverage, p0, method='trf', loss='soft_l1')
        print(f'N: {len(f_cloud_free)} x0: {popt[0]:.2f} k: {popt[1]:.4f}')
        ypredicted = sigmoid(f_cloud_free, *popt)
        slope, intercept, r_value, p_value, std_err = stats.linregress(f_coverage, ypredicted)

    return median_f_coverage, popt


def plot_linear_regression_lines(file_out):
    eras = get_era_info()
    from PlotScatter import PlotScatter
    pscatter = PlotScatter()
    # pscatter.close_plot()
    xdata = np.arange(0, 101)

    # style = {'color': 'blue', 'linestyle': '-', 'linewidth': 1.5, 'marker': None, 'markersize': 3}
    for era in eras:
        if era == 'All': continue
        print(era)
        regline = eras[era]['regline']
        ydata = (xdata * regline[0]) + regline[1]
        # style['color'] = eras[era]['color']
        pscatter.plot_regress_line(xdata, ydata, eras[era]['color'])

    # style = {'color': 'black', 'linestyle': '--', 'linewidth': 1, 'marker': None, 'markersize': 3}
    pscatter.plot_identity_line()

    pscatter.set_xaxis_title('Clould-free fraction (%)')
    pscatter.set_yaxis_title('Sensor coverage (%)')
    pscatter.set_ticks([0, 20, 40, 60, 80, 100], 0)

    pscatter.set_limits_Y(0, 100)
    pscatter.set_limits_X(0, 100)
    pscatter.set_equal_apect()

    pscatter.save_fig(file_out)
    pscatter.close_plot()


def plot_bloom_coverage(file_average):
    dir_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/CYANOBLOOM_EVOLUTION/COVERAGE_PLOTS'
    print('STARTED BLOOM COVERAGE')
    dataset = Dataset(file_average)
    daily_cloud_free_map = dataset.variables['daily_cloud_free_map'][:]
    n_coverage_map_corrected = dataset.variables['n_coverage_map_corrected'][:]
    n_bloom_map = dataset.variables['n_bloom_map'][:]
    f_cloud_free = dataset.variables['f_cloud_free'][:]
    f_time = dataset.variables['f_time'][:]
    time = dataset.variables['time'][:]
    n_coverage_corrected = dataset.variables['n_coverage_corrected'][:]
    daily_cloud_free = dataset.variables['daily_cloud_free_percent'][:]
    dataset.close()

    ##important this step
    n_bloom_map[daily_cloud_free_map.mask] = np.ma.masked

    ##by image
    file_out = os.path.join(dir_out, 'P_Bloom_vs_Cloud_Free_ByImage.tif')
    n_bloom = np.ma.sum(n_bloom_map, axis=1)
    n_bloom[np.where(n_coverage_corrected.mask)] = np.ma.masked
    p_bloom = (n_bloom / n_coverage_corrected) * 100
    p_bloom[np.where(np.logical_and(n_bloom == 0, n_coverage_corrected == 0))] = 0
    plot_sensor_coverage_versus_clould_free_impl(p_bloom, daily_cloud_free, 5, 'Bloom coverage(%) - Image', file_out)
    file_out = os.path.join(dir_out, 'P_Bloom_GreaterZero_vs_Cloud_Free_ByImage.tif')
    indices_g0 = np.where(p_bloom > 0)
    p_bloom = p_bloom[indices_g0]
    daily_cloud_free = daily_cloud_free[indices_g0]
    plot_sensor_coverage_versus_clould_free_impl(p_bloom, daily_cloud_free, 5, 'Bloom coverage(%) - Image', file_out)

    time = time[indices_g0]
    time_obj = [dt.utcfromtimestamp(t) for t in time]
    year_array = np.array([obj.year for obj in time_obj])
    groups = get_era_year_groups(year_array)
    file_out = os.path.join(dir_out, 'ScatterPlot_P_Bloom_GreaterZero_vs_Cloud_Free_ByImage.tif')
    plot_scatter_plot(daily_cloud_free, p_bloom, 'Cloud Free Fraction(%)', 'Bloom Coverage (%) - Image', [0, 100],
                      groups, True, file_out)
    for g in range(7):
        meaning = groups['meanings'][g]
        file_out = os.path.join(dir_out, f'ScatterPlot_P_Bloom_GreaterZero_vs_Cloud_Free_ByImage_{meaning}.tif')
        groups_here = {
            'array': groups['array'],
            'meanings': [meaning],
            'values': [groups['values'][g]],
            'colors': [groups['colors'][g]]
        }
        plot_scatter_plot(daily_cloud_free, p_bloom, 'Cloud Free Fraction(%)', 'Bloom Coverage (%) - Image', [0, 100],
                          groups_here, True, file_out)

    ##macropixel clara analysis
    # p_bloom_map = (n_bloom_map/n_coverage_map_corrected)*100
    # p_bloom_map[np.where(np.logical_and(n_bloom_map==0,n_coverage_map_corrected==0))]=0
    # p_bloom_map[p_bloom_map>100]=100.0
    # f_bloom = p_bloom_map[p_bloom_map.mask==False].flatten()
    # indices_g0 = np.where(f_bloom>0)
    # f_bloom = f_bloom[indices_g0]
    # f_cloud_free = f_cloud_free[indices_g0]
    # file_out = os.path.join(dir_out,'P_Bloom_WithoutZero_vs_Cloud_Free_ByMacropixelClara.tif')
    # plot_sensor_coverage_versus_clould_free_impl(f_bloom,f_cloud_free,5,'Bloom coverage(%)',file_out)


def prepare_df_bloom_coverage_by_year(file_average, file_mask, extra):
    dir_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/CYANOBLOOM_EVOLUTION/COVERAGE_PLOTS'
    dataset_m = Dataset(file_mask)
    # ntotal_byindex = dataset_m.variables[''][:]
    ntotal_water_cfc = dataset_m.variables['NTotal_Water_CFC'][:]
    dataset_m.close()
    dataset = Dataset(file_average)
    n_bloom_map = dataset.variables['n_bloom_map'][:]
    # daily_cloud_free_map = dataset.variables['daily_cloud_free_map'][:]
    n_coverage_map_corrected = dataset.variables['n_coverage_map_corrected'][:]
    n_expected_map = dataset.variables['n_expected_map'][:]
    n_coverage_corrected = dataset.variables['n_coverage_corrected'][:]
    n_expected = dataset.variables['n_expected_sum'][:]
    time = dataset.variables['time'][:]
    dataset.close()

    ##getting n_total
    n_time = len(time)
    n_total_by_day = np.zeros((n_time,))
    n_total_macro_by_day = np.zeros((n_time,))
    temp = n_coverage_map_corrected.copy()
    for iday in range(n_time):
        coverage_day = temp[iday, :]
        if extra == '_FullCoverage':
            expected_day = n_expected_map[iday, :]
            coverage_day[np.where(coverage_day < expected_day)] = np.ma.masked
            n_coverage_map_corrected[iday, :] = coverage_day
            expected_day[np.where(coverage_day < expected_day)] = np.ma.masked
            n_expected_map[iday, :] = expected_day
            n_bloom_map_day = n_bloom_map[iday, :]
            n_bloom_map_day[np.where(coverage_day < expected_day)] = np.ma.masked
            n_bloom_map[iday, :] = n_bloom_map_day
        if extra == '_FullCloudFree':
            expected_day = n_expected_map[iday, :]
            coverage_day[np.where(coverage_day < ntotal_water_cfc)] = np.ma.masked
            n_coverage_map_corrected[iday, :] = coverage_day
            expected_day[np.where(coverage_day < ntotal_water_cfc)] = np.ma.masked
            n_expected_map[iday, :] = expected_day
            n_bloom_map_day = n_bloom_map[iday, :]
            n_bloom_map_day[np.where(coverage_day < ntotal_water_cfc)] = np.ma.masked
            n_bloom_map[iday, :] = n_bloom_map_day

        indices_no_mask = np.where(coverage_day.mask == False)
        coverage_day[indices_no_mask] = ntotal_water_cfc[indices_no_mask]
        n_total_by_day[iday] = np.ma.sum(coverage_day)
        n_total_macro_by_day[iday] = len(indices_no_mask[0])

    if extra == '_FullCoverage' or extra == '_FullCloudFree':
        n_coverage_corrected = np.ma.sum(n_coverage_map_corrected, axis=1)
        n_expected = np.ma.sum(n_expected_map, axis=1)

    ##we obtain p_bloom_map, percentage of the real coverage classified as bloom
    n_bloom_map[n_coverage_map_corrected.mask] = np.ma.masked
    p_bloom_map = (n_bloom_map / n_coverage_map_corrected) * 100
    p_bloom_map[np.where(np.logical_and(n_bloom_map == 0, n_coverage_map_corrected == 0))] = 0
    p_bloom_map[np.where(np.logical_and(p_bloom_map.mask == False, p_bloom_map > 100))] = 100

    ##n_bloom_expected_map, number of bloom if we multiply the percentage of the real coverage by the expected number of pixels
    n_bloom_expected_map = np.ma.masked_all(p_bloom_map.shape)
    n_bloom_expected_map[p_bloom_map.mask == False] = (p_bloom_map[p_bloom_map.mask == False] * n_expected_map[
        p_bloom_map.mask == False]) / 100
    n_bloom_expected_map = np.ma.round(n_bloom_expected_map)
    # print('No adjusted: ',np.ma.count(n_bloom_expected_map),np.ma.sum(n_bloom_map),np.ma.sum(n_bloom_expected_map))
    # indices_equal = np.where(np.logical_and(n_bloom_map.mask==False,n_bloom_map==n_bloom_expected_map))
    # print('Equal',len(indices_equal[0]), ' of ',np.ma.count(n_bloom_expected_map))
    # indices_greater = np.where(np.logical_and(n_bloom_map.mask == False, n_bloom_map < n_bloom_expected_map))
    # print('Greater', len(indices_greater[0]), ' of ', np.ma.count(n_bloom_expected_map))
    indices_lower = np.where(np.logical_and(n_bloom_map.mask == False, n_bloom_map > n_bloom_expected_map))
    # print('Lower', len(indices_lower[0]), ' of ', np.ma.count(n_bloom_expected_map))
    n_bloom_expected_map_adjusted = n_bloom_expected_map.copy()
    n_bloom_expected_map_adjusted[indices_lower] = n_bloom_map[indices_lower]
    # indices_lower = np.where(np.logical_and(n_bloom_map.mask == False, n_bloom_map > n_bloom_expected_map_adjusted))
    # print('Lower Adjuted', len(indices_lower[0]), ' of ', np.ma.count(n_bloom_expected_map))
    # print('Adjusted',np.ma.count(n_bloom_expected_map), np.ma.sum(n_bloom_map), np.ma.sum(n_bloom_expected_map_adjusted))

    n_bloom_expected_map_max = np.ma.masked_all(p_bloom_map.shape)
    for iday in range(n_time):
        p_bloom_day = p_bloom_map[iday, :]
        n_bloom_expected_map_max_day = n_bloom_expected_map_max[iday, :]
        n_bloom_expected_map_max_day[p_bloom_day.mask == False] = (p_bloom_day[p_bloom_day.mask == False] *
                                                                   ntotal_water_cfc[p_bloom_day.mask == False]) / 100
        n_bloom_expected_map_max[iday, :] = n_bloom_expected_map_max_day[:]
    n_bloom_expected_map_max = np.round(n_bloom_expected_map_max)

    time_obj = [dt.utcfromtimestamp(t) for t in time]
    year_array = np.array([obj.year for obj in time_obj])
    jjj_array = np.array([int(obj.strftime('%j')) for obj in time_obj])
    select_array = np.ones(time.shape)
    select_array[year_array < 1998] = 0
    select_array[year_array > 2023] = 0
    select_array[jjj_array < 161] = 0
    select_array[jjj_array > 270] = 0

    # ##selection full sensor coverage
    # for itime in range(n_time):
    #     if n_coverage_corrected[itime]<n_expected[itime]:
    #         select_array[itime] = 0

    ##summer
    n_total_by_day_summer = n_total_by_day[select_array == 1]
    n_total_macro_by_day_summer = n_total_macro_by_day[select_array == 1]
    n_expected_summer = n_expected[select_array == 1]
    n_coverage_corrected_summer = n_coverage_corrected[select_array == 1]

    n_bloom_daily = np.ma.sum(n_bloom_map, axis=1)
    n_bloom_daily_summer = n_bloom_daily[select_array == 1]

    n_bloom_daily_expected = np.ma.sum(n_bloom_expected_map, axis=1)
    n_bloom_daily_expected_summer = n_bloom_daily_expected[select_array == 1]

    n_bloom_daily_expected_adjusted = np.ma.sum(n_bloom_expected_map_adjusted, axis=1)
    n_bloom_daily_expected_adjusted_summer = n_bloom_daily_expected_adjusted[select_array == 1]

    n_bloom_daily_expected_max = np.ma.sum(n_bloom_expected_map_max, axis=1)
    n_bloom_daily_expected_max_summer = n_bloom_daily_expected_max[select_array == 1]

    year_u = np.arange(1998, 2024)
    year_summer = year_array[select_array == 1]

    col_names = ['NTotalMacro', 'NTotal', 'NCoverageTotal', 'NExpectedTotal', 'NBloomTotal', 'NBloomExpected',
                 'NBloomExpectedAdjusted',
                 'NBloomExpectedMax']

    for idx in range(len(col_names)):
        col_names[idx] = f'{col_names[idx]}{extra}'

    df = pd.DataFrame(index=year_u, columns=col_names)
    for year in year_u:
        n_total_macro_year = n_total_macro_by_day_summer[year_summer == year]
        n_total_year = n_total_by_day_summer[year_summer == year]
        n_coverage_year = n_coverage_corrected_summer[year_summer == year]
        n_expected_year = n_expected_summer[year_summer == year]
        n_bloom_daily_year = n_bloom_daily_summer[year_summer == year]
        n_bloom_daily_expected_year = n_bloom_daily_expected_summer[year_summer == year]
        n_bloom_daily_expected_adjusted_year = n_bloom_daily_expected_adjusted_summer[year_summer == year]
        n_bloom_daily_expected_max_year = n_bloom_daily_expected_max_summer[year_summer == year]

        df.loc[year, f'NTotalMacro{extra}'] = np.nansum(n_total_macro_year)
        df.loc[year, f'NTotal{extra}'] = np.nansum(n_total_year)
        df.loc[year, f'NCoverageTotal{extra}'] = np.ma.sum(n_coverage_year)
        df.loc[year, f'NExpectedTotal{extra}'] = np.ma.sum(n_expected_year)
        df.loc[year, f'NBloomTotal{extra}'] = np.ma.sum(n_bloom_daily_year)
        df.loc[year, f'NBloomExpected{extra}'] = np.ma.sum(n_bloom_daily_expected_year)
        df.loc[year, f'NBloomExpectedAdjusted{extra}'] = np.ma.sum(n_bloom_daily_expected_adjusted_year)
        df.loc[year, f'NBloomExpectedMax{extra}'] = np.ma.sum(n_bloom_daily_expected_max_year)

    file_csv = os.path.join(dir_out, 'BloomCoverageSummer_ByYear.csv')
    if extra == '_FullCoverage':
        file_csv = os.path.join(dir_out, 'BloomCoverageSummer_ByYear_FullSensorCoverage.csv')
    if extra == '_FullCloudFree':
        file_csv = os.path.join(dir_out, 'BloomCoverageSummer_ByYear_FullCloudFree.csv')
    df.to_csv(file_csv, sep=';')


def plot_temporal_coverage_clara(dir_out):
    # temporal coverage clara
    file_out = os.path.join(dir_out, 'TemporalCoverageJDAY_CLARA.tif')
    from PlotSpectra import PlotSpectra
    pspectra = PlotSpectra()
    pspectra.start_plot_with_size(20, 1.5)
    xdata = np.arange(0, 365)
    data = np.ones((365,))
    data_group = np.zeros((365,))
    ## 1:partial corevareg,any year: 36-53, 280-310
    ## 2:partial coverage, all years: 54-67, 278-279
    ## 3: full coverage, some years: 68-89 242-277
    ## 4: full coverage, all years: 90-241
    data_group[35:53] = 1
    data_group[53:67] = 2
    data_group[67:89] = 3
    data_group[89:241] = 4
    data_group[241:277] = 3
    data_group[277:279] = 2
    data_group[279:310] = 1
    pspectra.xdata = xdata
    handles = pspectra.plot_single_bar_series_grouped(data, data_group, ['gray', 'blue', 'cyan', 'orange', 'red'], 1, 0,
                                                      0)
    xticks = []
    xticks_values = []
    xsep = []
    for imonth in range(1, 13, 1):
        date_here = dt(2020, imonth, 15)
        xticks.append(int(date_here.strftime('%j')))
        xticks_values.append(date_here.strftime('%b'))
        date_here = dt(2020, imonth, 1)
        xsep.append(int(date_here.strftime('%j')))
    xsep.append(365)
    xsep = np.array(xsep) - 1
    for x in xsep:
        pspectra.plot_single_data([x, x], [0, 1],
                                  {'color': 'k', 'linestyle': '-', 'linewidth': 1, 'marker': None, 'markersize': 0})
    xticks = np.array(xticks) - 1
    pspectra.set_xticks(xticks, xticks_values, 0, 18)
    pspectra.set_yticks([0], [''], 0, 0)
    str_legend = ['No Data', 'Partial coverage-some years', 'Partial coverage - all years',
                  'Full coverage - some years', 'Full coverage - all year']
    # pspectra.set_legend_h(handles,str_legend)

    # pspectra.set_grid()
    pspectra.set_tigth_layout()
    pspectra.save_plot(file_out)
    pspectra.close_plot()

def plot_temporal_coverage_cci(dir_out):
    file_average = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/CoverageAnalysis_BAL_MULTI_COMPLETED_19970904_20241031.nc'
    dataset = Dataset(file_average)
    n_coverage = dataset.variables['n_coverage'][:]
    n_coverage_adjusted = dataset.variables['n_coverage_corrected'][:]
    time = dataset.variables['time'][:]
    dataset.close()

    n_coverage_total = np.ma.sum(n_coverage)
    n_coverage_adjusted_total = np.ma.sum(n_coverage_adjusted)
    print(f'{n_coverage_total:1.2e}')
    print(f'{n_coverage_adjusted_total:1.2e}')
    increm = ((n_coverage_adjusted_total - n_coverage_total) / (n_coverage_total)) * 100
    print(f'{increm:0.2f}')

    time_obj = [dt.utcfromtimestamp(t) for t in time]
    year_array = np.array([obj.year for obj in time_obj])
    month_array = np.array([obj.month for obj in time_obj])
    jjj_array = np.array([int(obj.strftime('%j')) for obj in time_obj])
    jjj_array[year_array == 1997] = -1
    jjj_array[year_array == 2024] = -1
    data_group = np.zeros((365,))
    for jday in range(1, 366, 1):
        n_coverage_jday = n_coverage[jjj_array == jday]
        n_coverage_jday = np.ma.filled(n_coverage_jday, 0)
        n_coverage_jday[n_coverage_jday > 0] = 1
        ndays = np.sum(n_coverage_jday)
        if 1 <= ndays <= 8:
            data_group[jday - 1] = 1
        elif 9 <= ndays <= 16:
            data_group[jday - 1] = 2
        elif 17 <= ndays <= 25:
            data_group[jday - 1] = 3
        elif ndays == 26:
            data_group[jday - 1] = 4
    from matplotlib import colormaps as cm

    file_out = os.path.join(dir_out, 'TemporalCoverageJDAY_CCI_legend.tif')
    from PlotSpectra import PlotSpectra
    pspectra = PlotSpectra()
    # pspectra.start_plot_with_size(20, 1.5)
    xdata = np.arange(0, 365)
    data = np.ones((365,))
    pspectra.xdata = xdata
    handles = pspectra.plot_single_bar_series_grouped(data, data_group, ['gray', 'blue', 'cyan', 'orange', 'red'], 1, 0,
                                                      0)
    xticks = []
    xticks_values = []
    xsep = []
    for imonth in range(1, 13, 1):
        date_here = dt(2020, imonth, 15)
        xticks.append(int(date_here.strftime('%j')))
        xticks_values.append(date_here.strftime('%b'))
        date_here = dt(2020, imonth, 1)
        xsep.append(int(date_here.strftime('%j')))
    xsep.append(365)
    xsep = np.array(xsep) - 1
    for x in xsep:
        pspectra.plot_single_data([x, x], [0, 1],
                                  {'color': 'k', 'linestyle': '-', 'linewidth': 1, 'marker': None, 'markersize': 0})
    xticks = np.array(xticks) - 1
    pspectra.set_xticks(xticks, xticks_values, 0, 18)
    pspectra.set_yticks([0], [''], 0, 0)
    str_legend = ['No Data', '1-8 years', '9-16 years',
                  '17-25 years', 'All years (26)']
    pspectra.set_legend_h(handles, str_legend)

    # pspectra.set_grid()
    pspectra.set_tigth_layout()
    pspectra.save_plot(file_out)
    pspectra.close_plot()

def plot_maps_clara_cfc(dir_out):
    file_mask = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/BAL_Land_Mask_hr_CFC.nc'
    #file_cfc = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/CFCdm202307100000003UDAVPOSI1UD.nc'
    file_cfc = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/CFCdm201803020000003UDAVPOS01UD.nc'
    name_cfc = os.path.basename(file_cfc)
    file_out = os.path.join(dir_out,f'{name_cfc[:-3]}.tif')
    apply_mask = False
    use_cloud_free = False
    dataset_mask = Dataset(file_mask)
    land_mask = dataset_mask.variables['Land_Mask_CFC'][:]
    dataset_mask.close()
    dataset = Dataset(file_cfc)
    cfc_data = dataset.variables['cfc_day'][:]
    lat_array = dataset.variables['lat'][:]
    lon_array = dataset.variables['lon'][:]
    cfc_data = np.ma.squeeze(cfc_data)
    dataset.close()
    if apply_mask:
        cfc_data[land_mask==0]=np.ma.masked
    if use_cloud_free:
        cfc_data = 100 - cfc_data
    plot_baltic_map_impl(cfc_data,lat_array,lon_array,'cfc',file_out)

def plot_baltic_map_impl(data,lat_array,lon_array,type,file_out):
    start_full_figure()
    fig, ax = start_full_figure()

    if type.startswith('cfc'):
        h = ax.pcolormesh(lon_array, lat_array, data, vmin=0, vmax=100,cmap=mpl.colormaps['jet'])
    if type=='chl':
        h = ax.pcolormesh(lon_array, lat_array, data, norm=LogNorm(vmin=0.01, vmax=10))
    cbar = fig.colorbar(h, cax=None, ax=ax, use_gridspec=True, fraction=0.03, format="$%.2f$")
    cbar.ax.tick_params(labelsize=15)

    if type=='cfc':
        label = 'Cloud fraction area(%)'
    elif type=='cfc-free':
        label = 'Cloud-free fraction area(%)'
    elif type=='chl':
        label = 'chl'

    cbar.set_label(label=label, size=15)
    # title = f'Chlorophyll a concentration - OLD ({units})'
    # if dateherestr is not None:
    #     title = f'{title} - {dateherestr}'
    # ax.set_title(title, fontsize=20)
    # file_chla_old = os.path.join(os.path.dirname(output_path), f'Img_Chla_OLD_{dateherestr}.png')
    if type.startswith('cfc'):
        ax.coastlines(resolution='10m',color='white',linewidth=2)
    if type=='chl':
        ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black', linewidth=0.5)

    fig.savefig(file_out, dpi=300, bbox_inches='tight')

def plot_maps_cci(dir_out):
    file_mask = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/BAL_Land_Mask_hr_CFC.nc'
    # file_cfc = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/CFCdm202307100000003UDAVPOSI1UD.nc'
    file_chl = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/2009/220/C2009220-chl-bal-hr.nc'
    name_chl = os.path.basename(file_chl)
    file_out = os.path.join(dir_out, f'{name_chl[:-3]}.tif')
    apply_mask = False

    dataset = Dataset(file_chl)
    chl = dataset.variables['CHL'][:]
    chl = np.ma.squeeze(chl)
    lat = dataset.variables['lat'][:]
    lon = dataset.variables['lon'][:]
    dataset.close()

    plot_baltic_map_impl(chl,lat,lon,'chl',file_out)


def plot_era_sensors(dir_out):
    file_out = os.path.join(dir_out,'SensorEras.tif')
    start_date = dt(1997,9,4)
    end_date = dt(2024,11,30)

    sensors = {
        'SeaWiFS': {
            'start_time': dt(1997, 9, 4),
            'end_time': dt(2010, 12, 11),
            'yindex':5
        },
        'MERIS':{
            'start_time': dt(2002, 4, 29),
            'end_time': dt(2012, 4, 8),
            'yindex':4
        },
        'MODIS':{
            'start_time': dt(2002, 7, 4),
            'end_time': end_date,
            'yindex': 3
        },
        'VIIRS':{
            'start_time': dt(2012, 1, 2),
            'end_time': end_date,
            'yindex':2
        },
        'OLCI-A':{
            'start_time': dt(2016, 4, 25),
            'end_time': end_date,
            'yindex':1
        },
        'OLCI-B':{
            'start_time': dt(2018, 5, 14),
            'end_time': end_date,
            'yindex':0
        }
    }


    eras = {
        'SeaWiFS': {
            'yindex': [5],
            'color': 'blue',
            'start_time': dt(1997,9,4),
            'end_time': dt(2002,4,28)
        },
        'SeaWiFS_MERIS_MODIS': {
            'yindex': [5,4,3],
            'color': 'cyan',
            'start_time': dt(2002,4,29),
            'end_time': dt(2010,12,11)
        },
        'MERIS_MODIS':{
            'yindex': [4, 3],
            'color': 'lime',
            'start_time': dt(2010, 12, 12),
            'end_time': dt(2012, 4, 8)
        },
        'MERIS_MODIS_OVERWRITTEVIIRS': {
            'yindex': [2],
            'color': 'lime',
            'start_time': dt(2012, 1, 2),
            'end_time': dt(2012, 4, 8)
        },
        'VIIRS_MODIS':{
            'yindex': [3,2],
            'color': 'green',
            'start_time': dt(2012, 4, 9),
            'end_time': dt(2016, 4, 24)
        },

        'VIIRS_MODIS_OLCIA': {
            'yindex': [3,2,1],
            'color': 'salmon',
            'start_time': dt(2016, 4, 25),
            'end_time': dt(2018, 5, 13)
        },
        'VIIRS_MODIS_OLCIAB': {
            'yindex': [3, 2, 1,0],
            'color': 'red',
            'start_time': dt(2018, 5, 14),
            'end_time': dt(2019, 12, 31)
        },
        'OLCIAB': {
            'yindex': [1, 0],
            'color': 'magenta',
            'start_time': dt(2020, 1, 1),
            'end_time': end_date
        }
    }

    xticks = []
    xticks_minor = []
    labels_minor = []
    for year in range(1998,2026):
        date_here = dt(year,1,1)
        xticks.append((date_here-start_date).days)
        if year<=2024:
            date_here = dt(year, 6,30)
            xticks_minor.append((date_here - start_date).days)
            labels_minor.append(str(year))
    labels = ['']*len(xticks)




    plt.close('all')


    for sensor in sensors:
        stime = sensors[sensor]['start_time']
        etime = sensors[sensor]['end_time']
        yindex = sensors[sensor]['yindex']
        duration = (etime-stime).days+1
        sindex = (stime-start_date).days
        plt.barh(yindex,duration,left=sindex,align='center',color='white',edgecolor = 'black',linewidth=1)

    for era in eras:
        stime = eras[era]['start_time']
        etime = eras[era]['end_time']
        duration = (etime - stime).days + 1
        sindex = (stime - start_date).days
        indices = eras[era]['yindex']
        color = eras[era]['color']
        for yindex in indices:
            plt.barh(yindex, duration, left=sindex, align='center', color=color,edgecolor =color,linewidth=1)


    plt.xticks(xticks, labels)
    plt.grid(axis='x', linestyle=':', linewidth=0.5)

    plt.xticks(xticks_minor,labels_minor,rotation=90,minor=True)
    plt.tick_params(axis='x',which='minor',length=0)

    sensor_labels = list(sensors.keys())
    plt.yticks([5,4,3,2,1,0],sensor_labels)

    plt.savefig(file_out,dpi=300)


def plot_methods_plots():
    dir_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/METHODS_PLOTS'

    #plot_temporal_coverage_clara(dir_out)
    #plot_temporal_coverage_cci(dir_out)

    #plot_maps_clara_cfc(dir_out)
    #plot_maps_cci(dir_out)
    #plot_era_sensors(dir_out)
    #file_mask = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/BAL_Land_Mask_hr_CFC.nc'
    #plot_mask_clara(file_mask,dir_out)

    #plot_stats(dir_out)



def plot_bloom_coverage_by_year_from_df(file_csv):
    dir_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/CYANOBLOOM_EVOLUTION/COVERAGE_PLOTS'
    df = pd.read_csv(file_csv, sep=';')
    year = np.array(df['year'])
    n_bloom = df['NBloomTotal']
    n_bloom_expected = df['NBloomExpectedAdjusted']
    n_bloom_max = df['NBloomExpectedMax']

    yrange = [0, 3.1e06]

    cloud_free = df['cloud_free']
    sensor_expected = df['sensor_expected']
    sensor_coverage = df['sensor_coverage']

    # file_out = os.path.join(dir_out,'YearSeriesBlooms.tif')
    # plot_stacked_bars(year, np.stack([n_bloom],axis=1),'Year','# Bloom Pixels',yrange,None,file_out)
    #
    # file_out = os.path.join(dir_out, 'YearSeriesBloomsWithExpected.tif')
    # plot_stacked_bars(year, np.stack([n_bloom,n_bloom_expected],axis=1), 'Year', '# Bloom Pixels', yrange,None,None, file_out)
    #
    # file_out = os.path.join(dir_out, 'YearSeriesBloomsWithExpectedMax.tif')
    # plot_stacked_bars(year, np.stack([n_bloom, n_bloom_expected,n_bloom_max], axis=1), 'Year', '# Bloom Pixels', yrange,None, None,file_out)

    # file_out = os.path.join(dir_out,'YearPorcIncremWithExpected.tif')
    # plot_bars(year,np.array(df['increm_with_expected']),'Year','% Bloom increm',[0,100],file_out)
    # file_out = os.path.join(dir_out, 'YearPorcIncremWithMax.tif')
    # plot_bars(year, np.array(df['increm_with_max']), 'Year', '% Bloom increm', [0, 100], file_out)

    groups = get_era_year_groups(year)

    # file_out = os.path.join(dir_out,'ScatterPlotNBloomExpectedVSNBloom.tif')
    # plot_scatter_plot(n_bloom,n_bloom_expected,'# Blooms','# Expected Blooms',[0,3e06],groups,True,file_out)
    # file_out = os.path.join(dir_out, 'ScatterPlotNBloomMaxVSNBloom.tif')
    # plot_scatter_plot(n_bloom, n_bloom_max, '# Blooms', '# Maximum Blooms', [0,4e06],groups,True,file_out)

    # file_out = os.path.join(dir_out,'YearSeriesBloomsAndCloudFree.tif')
    # plot_stacked_bars(year, np.stack([n_bloom],axis=1),'Year','# Bloom Pixels',yrange,cloud_free,'Cloud Free Fraction(%)',file_out)
    # file_out = os.path.join(dir_out, 'YearSeriesBloomsAndSensorExpected.tif')
    # plot_stacked_bars(year, np.stack([n_bloom], axis=1), 'Year', '# Bloom Pixels', yrange, sensor_expected,'Observed CCI vs. Expected(%)', file_out)
    # file_out = os.path.join(dir_out, 'YearSeriesBloomsAndSensorCoverage.tif')
    # plot_stacked_bars(year, np.stack([n_bloom], axis=1), 'Year', '# Bloom Pixels', yrange, sensor_coverage, 'CCI Coverage(%)',file_out)

    # file_out = os.path.join(dir_out,'YearSeriesPBloomsAndCloudFree.tif')
    # plot_stacked_bars(year, np.stack([np.array(df['p_bloom'])],axis=1),'Year','% Bloom Pixels',[0,20],cloud_free,'Cloud Free Fraction(%)',file_out)
    # file_out = os.path.join(dir_out, 'YearSeriesPBloomsAndSensorExpected.tif')
    # plot_stacked_bars(year, np.stack([np.array(df['p_bloom'])], axis=1), 'Year', '% Bloom Pixels', [0,20], sensor_expected,'Observed CCI vs. Expected(%)', file_out)
    # file_out = os.path.join(dir_out, 'YearSeriesPBloomsAndSensorCoverage.tif')
    # plot_stacked_bars(year, np.stack([np.array(df['p_bloom'])], axis=1), 'Year', '% Bloom Pixels', [0,20], sensor_coverage, 'CCI Coverage(%)',file_out)

    # file_out = os.path.join(dir_out, 'ScatterPlot_NBloomsVSCloudFree.tif')
    # plot_scatter_plot(cloud_free, n_bloom, 'Cloud Free Fraction(%)', '#Blooms', None, groups, False, file_out)
    # file_out = os.path.join(dir_out, 'ScatterPlot_NBloomsVSSensorExpected.tif')
    # plot_scatter_plot(sensor_expected, n_bloom, 'Observed CCI vs. Expected(%)', '#Blooms', None, groups, False, file_out)
    # file_out = os.path.join(dir_out, 'ScatterPlot_NBloomsVSSensorCoverage.tif')
    # plot_scatter_plot(sensor_coverage, n_bloom, 'CCI Coverage(%)', '#Blooms', None, groups, False, file_out)
    #
    # file_out = os.path.join(dir_out, 'ScatterPlot_NBloomsExpectedVSCloudFree.tif')
    # plot_scatter_plot(cloud_free, n_bloom_expected, 'Cloud Free Fraction(%)', '#Expected Blooms', None, groups, False, file_out)
    # file_out = os.path.join(dir_out, 'ScatterPlot_NBloomsExpectedVSSensorExpected.tif')
    # plot_scatter_plot(sensor_expected, n_bloom_expected, 'Observed CCI vs. Expected(%)', '#Expected Blooms', None, groups, False,file_out)
    # file_out = os.path.join(dir_out, 'ScatterPlot_NBloomsExpectedVSSensorCoverage.tif')
    # plot_scatter_plot(sensor_coverage, n_bloom_expected, 'CCI Coverage(%)', '#Expected Blooms', None, groups, False, file_out)

    # file_out = os.path.join(dir_out, 'ScatterPlot_NBloomsMaxVSCloudFree.tif')
    # plot_scatter_plot(cloud_free, n_bloom_max, 'Cloud Free Fraction(%)', '#Max. Blooms', None, groups, False, file_out)
    # file_out = os.path.join(dir_out, 'ScatterPlot_NBloomsMaxVSSensorExpected.tif')
    # plot_scatter_plot(sensor_expected, n_bloom_max, 'Observed CCI vs. Expected(%)', '#Max. Blooms', None, groups, False,file_out)
    # file_out = os.path.join(dir_out, 'ScatterPlot_NBloomsMaxVSSensorCoverage.tif')
    # plot_scatter_plot(sensor_coverage, n_bloom_max, 'CCI Coverage(%)', '#Max. Blooms', None, groups, False, file_out)

    # file_out = os.path.join(dir_out, 'ScatterPlot_PBloomsVSCloudFree.tif')
    # plot_scatter_plot(cloud_free, df['p_bloom'], 'Cloud Free Fraction(%)', '% Blooms', [0,60], groups, True, file_out)
    # file_out = os.path.join(dir_out, 'ScatterPlot_PBloomsVSSensorExpected.tif')
    # plot_scatter_plot(sensor_expected, df['p_bloom'], 'Observed CCI vs. Expected(%)', '% Blooms', [0,125], groups, True, file_out)
    # file_out = os.path.join(dir_out, 'ScatterPlot_PBloomsVSSensorCoverage.tif')
    # plot_scatter_plot(sensor_coverage, df['p_bloom'], 'CCI Coverage(%)', '% Blooms', [0,60], groups, True, file_out)

    # file_out = os.path.join(dir_out, 'YearSeries_BloomAndPercentage.tif')
    # plot_stacked_bars(year, np.stack([n_bloom], axis=1), 'Year', '# Bloom Pixels', yrange, df['p_bloom'],'Bloom percent(%)', file_out)
    # file_out = os.path.join(dir_out, 'YearSeries_BloomExpectedAndPercentage.tif')
    # plot_stacked_bars(year, np.stack([n_bloom_expected], axis=1), 'Year', '# Expected Bloom Pixels', yrange, df['p_bloom_expected'],'Expected Bloom percent(%)', file_out)
    # file_out = os.path.join(dir_out, 'YearSeries_BloomMaxAndPercentage.tif')
    # plot_stacked_bars(year, np.stack([n_bloom_max], axis=1), 'Year', '# Max Bloom Pixels', yrange,df['p_bloom_max'], 'Max Bloom percent(%)', file_out)

    file_out = os.path.join(dir_out, 'YearSeriesBloomsFullCoverageAndNMacropixels.tif')
    plot_stacked_bars(year, np.stack([df['NBloomTotal_FullCoverage']], axis=1), 'Year', '# Bloom Pixels', [0, 2.2e6],
                      df['NTotalMacro_FullCoverage'], '#CLARA MACROPIXELS', file_out)

    file_out = os.path.join(dir_out, 'YearSeriesBloomsFullCloudFreeAndNMacropixels.tif')
    plot_stacked_bars(year, np.stack([df['NBloomTotal_FullCloudFree']], axis=1), 'Year', '# Bloom Pixels', [0, 1.5e6],
                      df['NTotalMacro_FullCloudFree'], '#CLARA MACROPIXELS', file_out)

    file_out = os.path.join(dir_out, 'YearSeriesBloomsFullCoverageAndNCCI.tif')
    plot_stacked_bars(year, np.stack([df['NBloomTotal_FullCoverage']], axis=1), 'Year', '# Bloom Pixels', [0, 2.2e6],
                      df['NCoverageTotal_FullCoverage'], '#CCI PIXELS', file_out)

    file_out = os.path.join(dir_out, 'YearSeriesBloomsFullCloudFreeAndNCCI.tif')
    plot_stacked_bars(year, np.stack([df['NBloomTotal_FullCloudFree']], axis=1), 'Year', '# Bloom Pixels', [0, 1.5e6],
                      df['NCoverageTotal_FullCloudFree'], '#CCIPIXELS', file_out)


def sigmoid(x, x0, k):
    y = 100 / (1 + np.exp(-k * (x - x0)))
    return (y)


def get_era_info():
    eras = {
        'All': {
            'start_time': dt(1997, 9, 4).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2024, 10, 31).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'black'
        },
        'SeaWiFS': {
            'start_time': dt(1997, 9, 4).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2002, 4, 28).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'blue',
            'regline': [0.67, -7.22]
        },
        'SeaWiFS_MERIS_MODIS': {
            'start_time': dt(2002, 4, 29).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2010, 12, 11).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'cyan',
            'regline': [1.04, -6.15]
        },
        'MERIS_MODIS': {
            'start_time': dt(2010, 12, 12).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2012, 4, 8).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'lime',
            'regline': [0.98, -8.03]
        },
        'VIIRS_MODIS': {
            'start_time': dt(2012, 4, 9).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2016, 4, 24).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'green',
            'regline': [1.07, -4.31]
        },
        'VIIRS_MODIS_OLCIA': {
            'start_time': dt(2016, 4, 25).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2018, 5, 13).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'salmon',
            'regline': [1.14, -2.38]
        },
        'VIIRS_MODIS_OLCIAB': {
            'start_time': dt(2018, 5, 14).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2019, 12, 31).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'red',
            'regline': [1.09, -1.80]
        },
        'OLCIAB': {
            'start_time': dt(2020, 1, 1).replace(tzinfo=pytz.utc).timestamp(),
            'end_time': dt(2024, 10, 31).replace(tzinfo=pytz.utc).timestamp(),
            'color': 'magenta',
            'regline': [1.00, -3.48]
        }
    }
    return eras


def get_era_year_groups(year):
    eras = {
        'SeaWiFS': [1997, 2001],
        'SeaWiFS-MERIS-MODIS': [2002, 2010],
        'MERIS-MODIS': [2011, 2011],
        'VIIRS-MODIS': [2012, 2015],
        'VIIRS-MODIS-OLCIA': [2016, 2017],
        'VIIRS-MODIS-OLCIAB': [2018, 2019],
        'OLCIAB': [2020, 2023]
    }
    colors = ['blue', 'cyan', 'lime', 'green', 'salmon', 'red', 'magenta']
    gmeanings = []
    gvalues = []
    garray = np.zeros(year.shape)
    index = 1
    for era in eras:
        year_range = eras[era]
        garray[np.where(np.logical_and(year >= year_range[0], year <= year_range[1]))] = index
        gmeanings.append(era)
        gvalues.append(index)
        index = index + 1
    groups = {
        'array': garray,
        'values': gvalues,
        'meanings': gmeanings,
        'colors': colors
    }

    return groups


def temp_doors():
    # file_coverage = '/mnt/c/DATA_LUIS/DOORS_WORK/COMPARISON_CMEMS_CERTO/Coverage_CMEMS_CERTO.nc'
    # from netCDF4 import Dataset
    # dataset = Dataset(file_coverage)
    # time = dataset.variables['time'][:]
    # chl_median = dataset.variables['CERTO_blended_chla_top_3_weighted_median'][:]
    # chl_N = dataset.variables['CERTO_blended_chla_top_3_weighted_N'][:]
    # owt_dominant = dataset.variables['CERTO_owt_dominant_OWT'][:]
    # dataset.close()
    # for idx in range(len(time)):
    #     time_here = dt.utcfromtimestamp(time[idx]).strftime('%Y-%m-%d')
    #     log_chl_median = np.log10(chl_median[idx])
    #     print(time_here, log_chl_median, owt_dominant[idx], chl_N[idx])

    # print('GETTING EXTRA DATES')
    # file_cmems = '/mnt/c/DATA_LUIS/DOORS_WORK/Extracts_2024/AERONET_OC/MDB_CMEMS_OLCI_300M_CMEMS_OBS-OC_BLK_BGC_20160401T000000_20231220T000000_AERONET_Galata_Platform.csv'
    # file_doors = '/mnt/c/DATA_LUIS/DOORS_WORK/MDBs/MDB_CSV/MDB_CERTO_OLCI_300M_CERTO-OLCI-L3_20160415T000000_20231110T000000_AERONET_Galata_Platform.csv'
    #
    # date_list = '/mnt/c/DATA_LUIS/DOORS_WORK/AERONET_DATE_CHECK/GalataDateList.csv'
    # file_out = '/mnt/c/DATA_LUIS/DOORS_WORK/AERONET_DATE_CHECK/GalataCERTO-OLCI.csv'
    # type = 'CERTO-OLCI'
    # file_ref = file_doors
    #
    # df = pd.read_csv(file_ref,sep=';')
    # time = df['SatelliteTimeStamp']
    # dates = []
    # for t in time:
    #     date_str = dt.utcfromtimestamp(t).strftime('%Y-%m-%d')
    #     if date_str not in dates:
    #         dates.append(date_str)
    #
    # fw = open(file_out,'w')
    # fw.write(f'DATE;{type}')
    # fr = open(date_list,'r')
    # for line in fr:
    #     date_here = line.strip()
    #     res = 0
    #     if date_here in dates:
    #         res = 1
    #     line_out = f'{date_here};{res}'
    #     fw.write('\n')
    #     fw.write(line_out)
    # fw.close()
    # fw.close()

    print('CHECKING SOURCES')
    file_dates = '/store3/DOORS/config_files/GalataDateList.csv'
    file_out = '/store3/DOORS/config_files/GalataDateList_SourcesCERTO_OLCI.csv'
    dir_sources = '/store/DOORS/CERTO_SOURCES'
    fw = open(file_out,'w')
    fw.write('DATE;SOURCE_CERTO_OLCI')
    fr = open(file_dates,'r')
    for line in fr:
        date_here_str = line.strip()
        date_here = dt.strptime(date_here_str,'%Y-%m-%d')
        dir_date = os.path.join(dir_sources,date_here.strftime('%Y'),date_here.strftime('%j'))
        name_out = f'CERTO_blk_{date_here.strftime("%Y%m%d")}_OLCI_RES300__final_l3_product.nc'
        file_nc = os.path.join(dir_date,name_out)
        res = 1 if os.path.exists(file_nc) else 0
        line_out = f'{date_here_str};{res}'
        fw.write('\n')
        fw.write(line_out)
    fr.close()
    fw.close()


def plot_stats(dir_out):
    print('here')
    plt.close('all')
    fig, ax_here = plt.subplots(subplot_kw=dict(projection='polar'))
    #ax_here.set_rscale('log')
    ax_here.set_rlim((0,100))
    ax_here.set_rticks([20, 40, 60, 80, 100])
    ax_here.set_theta_zero_location("N")
    ax_here.set_theta_direction(-1)

    APD = np.array([67,63,90,60])
    RPD = np.array([5,-5,42,-1])
    RPD = ((RPD-(-100))/200)*180
    angle = (RPD * np.pi) / 180
    file_out = os.path.join(dir_out,'test.tif')
    ax_here.scatter(angle, APD, marker='o',s=8, color='green')
    plt.savefig(file_out,dpi=300)


def main():
    file_mask = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/BAL_Land_Mask_hr_CFC.nc'
    # dir_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/MASK_PLOTS'
    # plot_mask_clara(file_mask,dir_out)
    #

    file_average = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/CoverageAnalysis_BAL_MULTI_COMPLETED_19970904_20241031.nc'
    dir_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/COVERAGE_ANALYSIS/PLOTS/'

    ##BAR TIME SERIES->csv must be done before using MDBReaderV2 -m PLOT
    # names_out = ['TotalCoverageByYear_Complete.tif','TotalCoverageByYear_Summer.tif','TotalCoverageByMonth.tif']
    # names_csv = ['time_series_year_complete.csv','time_series_year_summer.csv','time_series_month.csv']
    # time_vars = ['year','year','month']
    # data_var = 'daily_cloud_free_sum_sum'
    # for idx,name_out in enumerate(names_out):
    #     file_out = os.path.join(dir_base,name_out)
    #     file_csv = os.path.join(dir_base,names_csv[idx])
    #     time_var = time_vars[idx]
    #     df = pd.read_csv(file_csv, sep=';')
    #     time_array = df[time_var].to_numpy()
    #     data_array = df[data_var].to_numpy()
    #     data_array[np.isnan(data_array)] = 0
    #     data_array = data_array / 100
    #     range_y = None
    #     if idx==0: range_y = [50000,85000]
    #     if idx==1: range_y = [20000,45000]
    #     plot_bars(time_array, data_array, time_var.capitalize(), None,range_y,file_out)

    ##SMOOTH TIME SERIES
    # names_out = ['SmoothTimeSeries_DailyCloudFreeCoverage_YearComplete.tif',
    #              'SmoothTimeSeries_DailyCloudFreeCoverage_Summer.tif',
    #              'SmoothTimeSeries_TotalSensorCoverage_All.tif',
    #              'SmoothTimeSeries_TotalSensorCoverage_CLARA.tif',
    #              'SmoothTimeSeries_DailySensorCoverageVSExcpected.tif',
    #              'SmoothTimeSeries_DailySensorCoverage.tif'
    #              ]
    # names_csv_year = ['time_series_year_complete.csv']*6
    # names_csv_sensor = ['time_series_sensor.csv']*6
    # vars_data = ['daily_cloud_free_percent','daily_cloud_free_percent','n_coverage','n_coverage_corrected','p_coverage_cf','p_coverage']
    # y_axis_titles = ['Daily Cloud Free Fraction (%)','Daily Cloud Free Fraction (%)','Sensor coverage (#daily pixels)',
    #                  'Sensor coverage(#daily pixels)', 'Sensor coverage (%) vs. expected','Sensor coverage (%)']
    #
    # for idx,name_out in enumerate(names_out):
    #     var_data = vars_data[idx]
    #     y_ticks = None
    #     file_out = os.path.join(dir_base, name_out)
    #     plot_global = True
    #     if idx<=1:
    #         file_csv_year = os.path.join(dir_base, names_csv_year[idx])
    #         file_csv_sensor = None
    #         y_ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80]
    #     else:
    #         file_csv_year = None
    #         file_csv_sensor = os.path.join(dir_base,names_csv_sensor[idx])
    #     if idx >= 4:
    #         y_ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80,90,100]
    #         plot_global = False
    #     y_axis_title = y_axis_titles[idx]
    #     time_ranges = {
    #         'year': [1998, 2023]
    #     }
    #     if idx==1: time_ranges['jday'] = [161,270]
    #
    #     plot_time_series_smoothed(file_average,file_out,file_csv_year,file_csv_sensor,var_data,time_ranges,plot_global,y_axis_title,y_ticks)

    ##BOX PLOT COMPARISONS
    # names_out = ['BoxPlot_DailyCloudFreeByYear_Complete.tif',
    #              'BoxPlot_DailyCloudFreeByYear_Summer.tif',
    #              'BoxPlot_DailyCloudFreeByMonth.tif',
    #              'BoxPlot_DailyCloudFreeByEra.tif',
    #              'BoxPlot_DailyCoverageByEra.tif',
    #              'BoxPlot_DailyCoveragevsExpectedyEra.tif',
    #              'BoxPlot_MacropixelCloudFreeByEra.tif',
    #              'BoxPlot_MacropixelCoverageByEra.tif',
    #              'BoxPlot_MacropixelCoveragevsExpectedbyEra.tif']
    # vars_data = ['daily_cloud_free_percent']*4
    # vars_data = vars_data + ['p_coverage','p_coverage_cf','f_cloud_free','f_coverage','f_coverage_cf']
    # y_axis_titles = ['Daily Cloud Free Fraction (%)']*4
    # y_axis_titles = y_axis_titles + ['CCI Coverage','CCI Coverage VS. Expected','Cloud Free Fraction(%)','CCI Coverage(%)','CCI Coverage VS. Expected']
    #
    # x_axiss = ['year','year','month','era','era','era','era','era','era']
    # for idx,name_out in enumerate(names_out):
    #     if idx<=5:
    #         var_time='time'
    #         plot_global = True
    #         plot_average = False
    #     else:
    #         var_time='f_time'
    #         plot_global = False
    #         plot_average = True
    #
    #     file_out = os.path.join(dir_base,name_out)
    #     var_data = vars_data[idx]
    #     time_ranges = {
    #         'year': [1998, 2023]
    #     }
    #     if idx==1: time_ranges['jday'] = [161,270]
    #     y_axis_title = y_axis_titles[idx]
    #     y_ticks = [0,10,20,30,40,50,60,70,80,90,100]
    #     x_axis = x_axiss[idx]
    #     if x_axis=='era': time_ranges= {}
    #     plot_bounding_boxes(file_average,var_time,var_data,time_ranges,x_axis,y_axis_title,y_ticks,plot_global,plot_average, file_out)

    ##PLOT_BOX_MULTIPLE - COVERAGE BY IMAGE
    # file_out = os.path.join(dir_base,'MultipleBoxPlot_CoverageByImage.tif')
    # vars_data = ['daily_cloud_free_percent','p_coverage']
    # y_ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # series = {
    #     'CLARA': 'blue',
    #     'OC-CCI': 'green',
    # }
    # plot_multiple_bounding_boxes(file_average,'time',vars_data,None, {},'era','Coverage(%)',y_ticks,True,series,file_out)

    ##PLOT_BOX_MULTIPLE - COVERAGE BY IMAGE
    # file_out = os.path.join(dir_base,'MultipleBoxPlot_CoverageByMacropixel.tif')
    # vars_data = ['f_cloud_free','f_coverage']
    # y_ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # series = {
    #     'CLARA': 'blue',
    #     'OC-CCI': 'green',
    # }
    # plot_multiple_bounding_boxes(file_average,'f_time',vars_data,None, {},'era','Coverage(%)',y_ticks,True,series,file_out)

    # LINEAR REGRESSION LINE
    # file_out = os.path.join(dir_base, 'RegressionLines_DailyCoverage.tif')
    # plot_linear_regression_lines(file_out)

    ##SINGLE LINE TIME-SERIES
    # names_out = ['LineTimeSeries_TotalCloudFree_ByJulianDay',
    #              'LineTimeSeries_CloudFreeFraction_ByJulianDay']
    # names_csv = ['time_series_jday_complete.csv']*2
    # data_vars = ['daily_cloud_free_sum','daily_cloud_free_percent']
    # time_vars = ['jday']*2
    # y_axis_titles = ['Total cloud-free pixels','Cloud-free fraction']
    # for idx,name_out in enumerate(names_out):
    #     file_out = os.path.join(dir_base,name_out)
    #     file_csv = os.path.join(dir_base,names_csv[idx])
    #     data_var = data_vars[idx]
    #     is_percent = data_var.endswith('percent')
    #     plot_lines(file_csv,time_vars[idx],data_var,is_percent,y_axis_titles[idx],file_out)

    # PLOTTING SENSOR COVERAGE Versus CLOUD-FREE
    # plot_sensor_coverage_versus_clould_free(file_average,5)

    # PLOTTING BLOOM COVERAGE
    # plot_bloom_coverage(file_average)
    # prepare_df_bloom_coverage_by_year(file_average,file_mask,'_FullCloudFree')
    # file_csv = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/CYANOBLOOM_EVOLUTION/COVERAGE_PLOTS/BloomCoverageSummer_ByYear.csv'
    # plot_bloom_coverage_by_year_from_df(file_csv)

    # PLOTING MAPS CLARA RESOLUTION
    # refs = ['TotalCoverageMap_YearComplete','TotalCoverageMap_Summer','PercentCoverageMap_YearComplete','PercentCoverageMap_Summer']
    # for ref in refs:
    #     plot_maps_clara(file_average,ref,False)
    #     if ref.startswith('Percent'):
    #         plot_maps_clara(file_average,ref,True)

    ##METHODS PLOTS
    #plot_methods_plots()
    temp_doors()


    # ##until 2024-08-31: 0:9859
    # dataset_mask = Dataset(file_mask)
    # nwater = dataset_mask.variables['NTotal_Water_CFC'][:]
    # dataset_mask.close()
    # dataset = Dataset(file_average)
    # daily_cloud_free_map = dataset.variables['daily_cloud_free_map'][0:9859][:]
    # daily_cloud_free_percent = dataset.variables['daily_cloud_free_percent'][0:9859][:]
    # ncoverage_corrected = dataset.variables['p_coverage'][0:9859]
    # n_expected_map = dataset.variables['n_expected_map'][0:9859][:]
    # time = dataset.variables['time'][0:9859]
    # n_expected_sum = dataset.variables['n_expected_sum'][0:9859]
    # n_coverage_sum = dataset.variables['n_coverage'][0:9859]
    #
    # coverage_vs_expected = (n_coverage_sum/n_expected_sum)*100
    # for idx in range(len(coverage_vs_expected)):
    #     if coverage_vs_expected[idx]>500:
    #         print(n_expected_sum[idx],n_coverage_sum[idx],coverage_vs_expected[idx])
    # print(np.min(coverage_vs_expected))
    # print(np.max(coverage_vs_expected))
    # print(np.median(coverage_vs_expected))



    #
    # dataset.close()

    # print(np.ma.count(daily_cloud_free_percent))
    # print(np.ma.median(daily_cloud_free_percent))
    #
    # print(np.ma.count(ncoverage_corrected))
    # print(np.ma.median(ncoverage_corrected))
    # print(np.ma.mean(ncoverage_corrected))

    # print(np.ma.count(n_expected_map))
    # print(np.ma.sum(n_expected_map))

    # n_time = len(time)
    # n_water_daily = np.zeros((n_time,))
    # for itime in range(n_time):
    #     daily_cloud_free_map_here = np.ma.squeeze(daily_cloud_free_map[itime][:])
    #     n_water_here = nwater[daily_cloud_free_map_here.mask==False]
    #     n_water_daily[itime] = np.sum(n_water_here)
    #
    # n_water_total = np.sum(n_water_daily)
    # print(n_water_total)





if __name__ == '__main__':
    main()
