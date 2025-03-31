import os

import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import cartopy.feature as cfeature


class MapPlot:

    def __init__(self):
        self.path_extracts = None
        self.path_orig = None
        self.projection = ccrs.Mercator()
        self.ax = None

        self.box_style_default = {
            'color': 'k',
            'marker': 'o',
            'markersize':0,
            'linestyle': '--',
            'linewidth': 0.5
        }
        self.point_style_default = {
            'color': 'k',
            'marker': '.',
            'markersize': 1,
            'linestyle': None,
            'linewidth': 0
        }

    def start_map(self,geo_limits):
        extent = (geo_limits[2], geo_limits[3], geo_limits[0], geo_limits[1])
        self.ax = plt.axes(projection=self.projection, extent=extent)

    def plot_array(self,lat_array,lon_array,data_array,colormap,minv,maxv,use_log):
        if minv is None and maxv is None:
            plt.pcolormesh(lon_array, lat_array, data_array, transform=ccrs.PlateCarree(), cmap=colormap)
        else:
            if use_log:
                plt.pcolormesh(lon_array, lat_array, data_array, transform=ccrs.PlateCarree(), cmap=colormap,
                               norm=LogNorm(vmin=minv, vmax=maxv))
            else:
                plt.pcolormesh(lon_array, lat_array, data_array, transform=ccrs.PlateCarree(), cmap=colormap,
                               norm=Normalize(vmin=minv, vmax=maxv))

    def set_color_bar(self,shrink):
        if shrink is not None:
            shrink = 1
        plt.colorbar(shrink=shrink)

    def set_grid_lines(self):
        gl = self.ax.gridlines(linewidth=0.5, linestyle='dotted', draw_labels=True)
        return gl
    def top_labels(self,gl,b):
        gl.top_labels = b
        gl.xlabels_top = b
    def right_labels(self, gl, b):
        gl.right_labels = b
        gl.ylabels_right = b


    # def top_labels(self, gl, b):
    #     try:
    #             gl.top_labels = b
    #         except:
    #             pass
    #
    #     gl.top_labels = False
    #     gl.right_labels = False

    def plot_impl(self,lat,lon,style):
        style_here = self.box_style_default.copy()
        if style is not None:
            style_here.update(style)
        plt.plot(lon, lat, color=style_here['color'],
                 marker=style_here['marker'],
                 markersize=style_here['markersize'],
                 linestyle=style_here['linestyle'],
                 linewidth=style_here['linewidth'],
                 transform=ccrs.PlateCarree())

    def plot_box(self,lat_box,lon_box,style):
        style_here = self.box_style_default.copy()
        if style is not None:
            style_here.update(style)
        self.plot_impl(lat_box,lon_box,style_here)

    def plot_point(self,lat_point,lon_point,style):
        style_here = self.point_style_default.copy()
        if style is not None:
            style_here.update(style)
        self.plot_impl(lat_point,lon_point,style_here)

    def add_land(self):
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face',
                                                facecolor=cfeature.COLORS['land'])
        self.ax.add_feature(land_10m, zorder=0, edgecolor='black', linewidth=0.5)

    def add_title(self,title):
        plt.title(title)

    def set_tight_layout(self):
        plt.tight_layout()

    def save_map(self,file_out):
        if file_out.endswith('.tif'):
            plt.savefig(file_out, dpi=300, bbox_inches='tight', pil_kwargs={"compression": "tiff_lzw"})
        else:
            plt.savefig(file_out, dpi=300, bbox_inches='tight')

    def close_map(self):
        plt.close()



def plot_doors():
    from netCDF4 import Dataset
    from datetime import datetime as dt
    from MDBPlotV3 import MDBPlot
    mplot = MDBPlot(None)
    if os.path.isdir('/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/'):
        dir_extracts = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/DOORS_WORK/Extracts_2024/extracts_cmems_olci'
        dir_sources = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/DOORS_WORK/SOURCES'
        dir_out = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/DOORS_WORK/QL'
    else:
        dir_extracts = '/store3/DOORS/Extracts_2024/extracts_cmems_olci'
        dir_sources = '/dst04-data1/OC/OLCI/daily_v202311_bc'
        dir_out = '/store3/DOORS/quicklooks'



    extract_list = dict()
    for name in os.listdir(dir_extracts):
        print(name)
        date_here = dt.strptime(name.split('_')[4], '%Y%m%d')
        date_here_key = date_here.strftime('%Y%m%d')
        yyyy = date_here.strftime('%Y')
        jjj = date_here.strftime('%j')
        file_cmems = os.path.join(dir_sources,yyyy,jjj,f'O{yyyy}{jjj}-chl-bs-fr.nc')
        if not os.path.exists(file_cmems):
            continue
        print(f'[INFO] {name} --> {name.split('_')[4]}')
        file_extract = os.path.join(dir_extracts,name)
        dataset = Dataset(file_extract)
        lat_array = np.squeeze(dataset.variables['satellite_latitude'][:])
        lon_array = np.squeeze(dataset.variables['satellite_longitude'][:])
        dataset.close()
        lat_min = np.min(lat_array)
        lat_max = np.max(lat_array)
        lon_min = np.min(lon_array)
        lon_max = np.max(lon_array)
        lat_box,lon_box = mplot.get_box(lat_array,lon_array,-1)

        if date_here_key not in extract_list:
            extract_list[date_here_key]={
                'geo_limits': [lat_min,lat_max,lon_min,lon_max],
                'lat_centers': [lat_array[12,12]],
                'lon_centers':[lon_array[12,12]],
                'lat_boxes': [lat_box],
                'lon_boxes': [lon_box],
                'file_source': file_cmems,
                'title': f'CMEMS CHL-A {date_here.strftime("%Y-%m-%d")}',
                'file_out': os.path.join(dir_out,f'CMEMS_CHLA_{date_here_key}.png')
            }
        else:
            geo_limits = extract_list[date_here_key]['geo_limits']
            # lat_centers = extract_list[date_here.strftime('%Y%m%d')]['geo_limits']
            # lat_centers.append(date_here.strftime('%Y%m%d'),lat_array[12,12])
            # print(len(lat_centers))
            if lat_min < geo_limits[0]: geo_limits[0] = lat_min
            if lat_max > geo_limits[1]: geo_limits[1] = lat_max
            if lon_min < geo_limits[2]: geo_limits[2] = lon_min
            if lon_max > geo_limits[3]: geo_limits[3] = lon_max
            extract_list[date_here_key]['geo_limits'] = geo_limits
            extract_list[date_here_key]['lat_centers'].append(lat_array[12,12])
            extract_list[date_here_key]['lon_centers'].append(lon_array[12, 12])
            extract_list[date_here_key]['lat_boxes'].append(lat_box)
            extract_list[date_here_key]['lon_boxes'].append(lon_box)

    for name in extract_list:
        print(f'{name} ---> {extract_list[name]["geo_limits"]}, {len(extract_list[name]["lat_boxes"])}')
        plot_doors_impl(extract_list[name])

def plot_doors_impl(info):
    from netCDF4 import Dataset
    file_nc = info['file_source']
    dataset = Dataset(file_nc)
    lat = dataset.variables['lat'][:]
    lon = dataset.variables['lon'][:]
    chl = np.squeeze(dataset.variables['CHL'][:])
    dataset.close()
    lat_center = np.mean(info['lat_centers'])
    lon_center = np.mean(info['lon_centers'])
    ilat = np.argmin(np.abs(lat-lat_center))
    ilon = np.argmin(np.abs(lon-lon_center))

    nlat_map = int(len(lat)/4)
    nlon_map = int(len(lon)/4)
    ilat_min = int(ilat-(nlat_map/2))
    if ilat_min<0: ilat_min=0
    ilat_max = ilat_min+nlat_map
    ilon_min = int(ilon - (nlon_map / 2))
    if ilon_min < 0: ilon_min = 0
    ilon_max = ilon_min + nlon_map

    array_map = chl[ilat_min:ilat_max,ilon_min:ilon_max]
    lat_map = lat[ilat_min:ilat_max]
    lon_map = lon[ilon_min:ilon_max]
    geo_limits = [lat[ilat_min],lat[ilat_max],lon[ilon_min],lon[ilon_max]]

    geo_limits_extracts = info['geo_limits']
    lat_bbox_extracts = [geo_limits_extracts[0],geo_limits_extracts[1],geo_limits_extracts[1],geo_limits_extracts[0],geo_limits_extracts[0]]
    lon_bbox_extracts = [geo_limits_extracts[2],geo_limits_extracts[2],geo_limits_extracts[3],geo_limits_extracts[3],geo_limits_extracts[2]]

    lat_centers = info['lat_centers']
    lon_centers = info['lon_centers']

    mplot = MapPlot()
    mplot.start_map(geo_limits)
    gl = mplot.set_grid_lines()
    mplot.top_labels(gl,False)
    mplot.right_labels(gl,False)
    mplot.plot_array(lat_map,lon_map,array_map,'viridis',0.01,100,True)
    mplot.add_land()
    mplot.plot_box(lat_bbox_extracts,lon_bbox_extracts,None)
    for lat_c,lon_c in zip(lat_centers,lon_centers):
        mplot.plot_point(lat_c,lon_c,None)
    mplot.set_color_bar(0.5)
    mplot.add_title(info['title'])
    mplot.set_tight_layout()
    mplot.save_map(info['file_out'])
    mplot.close_map()



def main():
    print('[INFO] Started map plot')
    plot_doors()

if __name__ == '__main__':
    main()