import os
import shutil

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
    key = 'CERTO_TOP3'
    if os.path.isdir('/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/'):
        dir_extracts = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/DOORS_WORK/Extracts_2024/extracts_cmems_olci'
        dir_sources = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/DOORS_WORK/SOURCES'
        dir_out = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/DOORS_WORK/QL'
    else:
        dir_extracts = '/store3/DOORS/Extracts_2024/extracts_cmems_olci'
        if key=='CMEMS':
            dir_sources = '/dst04-data1/OC/OLCI/daily_v202311_bc'
        elif key.startswith('CERTO'):
            #dir_sources = '/store/DOORS/CERTO_SOURCES'
            dir_sources = '/store3/DOORS/CERTO_SOURCES'
        dir_out = '/store3/DOORS/quicklooks'



    extract_list = dict()
    for name in os.listdir(dir_extracts):
        print(f'[INFO] {name} --> {name.split("_")[4]}')
        date_here = dt.strptime(name.split('_')[4], '%Y%m%d')
        date_here_key = date_here.strftime('%Y%m%d')
        yyyy = date_here.strftime('%Y')
        jjj = date_here.strftime('%j')
        mm = date_here.strftime('%m')
        dd = date_here.strftime('%d')
        if key=='CMEMS':
            file_source = os.path.join(dir_sources,yyyy,jjj,f'O{yyyy}{jjj}-chl-bs-fr.nc')
        elif key.startswith('CERTO'):
            #file_source = os.path.join(dir_sources,yyyy,jjj,f'CERTO_blk_{yyyy}{mm}{dd}_OLCI_RES300__final_l3_product.nc')
            if date_here.year==2024:
                file_source = os.path.join(dir_sources, yyyy, jjj,
                                           f'CERTO_doors_cr3_olci_v4.16.1_{yyyy}{mm}{dd}_OLCI_RES300__final_l3_product.nc')
            else:
                file_source = os.path.join(dir_sources, yyyy, jjj,
                                       f'CERTO_blacksea_{yyyy}{mm}{dd}_OLCI_RES300__final_l3_product.nc')
        if not os.path.exists(file_source):
            continue
        try:
            dataset_s = Dataset(file_source)
            dataset_s.close()
        except:
            continue

        print(f'[INFO] Working...')
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
                'file_source': file_source,
                'title': f'{key} CHL-A {date_here.strftime("%Y-%m-%d")}',
                'file_out': os.path.join(dir_out,f'{key}_CHLA_{date_here_key}.png')
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
        plot_doors_impl(extract_list[name],key)

def plot_doors_impl(info: object, key: object) -> None:
    from netCDF4 import Dataset
    file_nc = info['file_source']
    dataset = Dataset(file_nc)
    lat = dataset.variables['lat'][:]
    lon = dataset.variables['lon'][:]
    if key=='CMEMS':
        chl = np.ma.squeeze(dataset.variables['CHL'][:])
    elif key=='CERTO_TOP3':
        chl = np.ma.squeeze(dataset.variables['blended_chla_top_3_weighted'][:])
    dataset.close()

    if 'geo_limits_image' in info:
        geo_limits_image = info['geo_limits_image']
        ilat1 = np.argmin(np.abs(geo_limits_image[0] - lat))
        ilat2 = np.argmin(np.abs(geo_limits_image[1] - lat))
        ilon1 = np.argmin(np.abs(geo_limits_image[2] - lon))
        ilon2 = np.argmin(np.abs(geo_limits_image[3] - lon))
        ilat_min = ilat1 if ilat1 < ilat2 else ilat2
        ilat_max = ilat2 if ilat2 > ilat1 else ilat1
        ilon_min = ilon1 if ilon1 < ilon2 else ilon2
        ilon_max = ilon2 if ilon2 > ilon1 else ilon1
        nlat_map = ilat_max - ilat_min
        nlon_map = ilon_max - ilon_min
        if ilat_min<0:
            ilat_min=0
            ilat_max = ilat_min+nlat_map
        if ilat_max>len(lat):
            ilat_max = len(lat)
            ilat_min = ilat_max-nlat_map
        if ilon_min < 0:
            ilon_min = 0
            ilon_max = ilon_min + nlon_map
        if ilon_max>len(lon):
            ilon_max = len(lon)
            ilon_min = ilon_max-nlon_map
    else:
        lat_center = np.mean(info['lat_centers'])
        lon_center = np.mean(info['lon_centers'])
        ilat = np.argmin(np.abs(lat-lat_center))
        ilon = np.argmin(np.abs(lon-lon_center))
        nlat_map = int(len(lat)/4)
        nlon_map = int(len(lon)/4)
        ilat_min = int(ilat-(nlat_map/2))
        if ilat_min<0:
            ilat_min=0
            ilat_max = ilat_min+nlat_map
        if ilat_max>len(lat):
            ilat_max = len(lat)
            ilat_min = ilat_max-nlat_map
        ilon_min = int(ilon - (nlon_map / 2))
        if ilon_min < 0:
            ilon_min = 0
            ilon_max = ilon_min + nlon_map
        if ilon_max>len(lon):
            ilon_max = len(lon)
            ilon_min = ilon_max-nlon_map




    array_map = chl[ilat_min:ilat_max,ilon_min:ilon_max]
    lat_map = lat[ilat_min:ilat_max]
    lon_map = lon[ilon_min:ilon_max]
    geo_limits = [lat[ilat_min],lat[ilat_max-1],lon[ilon_min],lon[ilon_max-1]]

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


def plot_multiple():
    from datetime import datetime as dt
    dir_plots = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/DOORS_WORK/QL'
    dir_out = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/DOORS_WORK/QL/MULTIPLE'
    for name in os.listdir(dir_plots):
        if name.startswith('CMEMS'):
            print(name)
            date_here = dt.strptime(name.split('_')[2][:-4], '%Y%m%d')
            if date_here.year<2023:
                continue
            file_cmems = os.path.join(dir_plots,name)
            name_certo = name.replace('CMEMS','CERTO_TOP3')
            file_certo = os.path.join(dir_plots,name_certo) if os.path.exists(os.path.join(dir_plots,name_certo)) else None
            file_out = os.path.join(dir_out,name.replace('CMEMS','CMEMS_CERTO_TOP3'))
            if file_certo is not None:
                from PlotMultiple import PlotMultiple
                pm = PlotMultiple()
                pm.start_multiple_plot(2,1)
                pm.plot_image(file_cmems,0,0)
                pm.plot_image(file_certo,1,0)
                pm.save_fig(file_out)
                pm.close_plot()
            else:
                shutil.copy(file_cmems,file_out)

def plot_figure_20190523():
    from netCDF4 import Dataset
    from datetime import datetime as dt
    dir_out = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/DOORS_WORK/QL/20190523'
    dir_sources = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/DOORS_WORK/SOURCES'
    dir_extracts = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/DOORS_WORK/Extracts_2024/extracts_cmems_olci'
    date_here = dt(2019,5,23)
    yyyy = date_here.strftime('%Y')
    jjj = date_here.strftime('%j')
    mm = date_here.strftime('%m')
    dd = date_here.strftime('%d')

    ##Getting extract info
    geo_extracts_limits = None
    lat_points = []
    lon_points = []
    for name in os.listdir(dir_extracts):
        if not name.find('20190523')>0:
            continue
        print(f'[INFO] Working with extract: {name}')
        file_extract = os.path.join(dir_extracts, name)
        dataset = Dataset(file_extract)
        lat_array = np.squeeze(dataset.variables['satellite_latitude'][:])
        lon_array = np.squeeze(dataset.variables['satellite_longitude'][:])
        dataset.close()
        geo_limits_here = [np.min(lat_array),np.max(lat_array),np.min(lon_array),np.max(lon_array)]
        if geo_extracts_limits is None:
            geo_extracts_limits = geo_limits_here
        else:
            if geo_limits_here[0] < geo_extracts_limits[0]: geo_extracts_limits[0] = geo_limits_here[0]
            if geo_limits_here[1] > geo_extracts_limits[1]: geo_extracts_limits[1] = geo_limits_here[1]
            if geo_limits_here[2] < geo_extracts_limits[2]: geo_extracts_limits[2] = geo_limits_here[2]
            if geo_limits_here[3] > geo_extracts_limits[3]: geo_extracts_limits[3] = geo_limits_here[3]
        lat_points.append(lat_array[12,12])
        lon_points.append(lon_array[12,12])

    ##cmems olci
    info = {
        'geo_limits': geo_extracts_limits,
        'lat_centers': lat_points,
        'lon_centers': lon_points,
        'file_source': os.path.join(dir_sources, yyyy, jjj, f'O{yyyy}{jjj}-chl-bs-fr.nc'),
        'title': f'CMEMS-OLCI CHL-A {date_here.strftime("%Y-%m-%d")}',
        'file_out': os.path.join(dir_out,f'CMEMS_OLCI_CHLA_{date_here.strftime("%Y%m%d")}.png'),
        'geo_limits_image': [43.6,45.6,28.6,32.6]
    }

    for lat,lon in zip(info['lat_centers'],info['lon_centers']):
        print(lat,';',lon)

    # file_cmems_olci = info['file_out']
    # #plot_doors_impl(info,'CMEMS')
    #
    # info['file_source']=os.path.join(dir_sources, yyyy, jjj, f'X{yyyy}{jjj}-chl-bs-hr.nc')
    # info['title']=f'CMEMS-MULTI CHL-A {date_here.strftime("%Y-%m-%d")}'
    # info['file_out'] = os.path.join(dir_out,f'CMEMS_MULTI_CHLA_{date_here.strftime("%Y%m%d")}.png')
    # file_cmems_multi = info['file_out']
    # # plot_doors_impl(info,'CMEMS')
    #
    # info['file_source'] =  os.path.join(dir_sources,yyyy,jjj,f'CERTO_blk_{yyyy}{mm}{dd}_OLCI_RES300__final_l3_product.nc')
    # info['title'] = f'CERTO CHL-A {date_here.strftime("%Y-%m-%d")}'
    # info['file_out'] = os.path.join(dir_out, f'CERTO_CHLA_{date_here.strftime("%Y%m%d")}.png')
    # file_certo = info['file_out']
    # # plot_doors_impl(info, 'CERTO_TOP3')
    #
    # file_out = os.path.join(dir_out,f'CHLA_{date_here.strftime("%Y%m%d")}.png')
    # from PlotMultiple import PlotMultiple
    # pm = PlotMultiple()
    # pm.start_multiple_plot(3,1)
    # pm.plot_image(file_cmems_olci,0,0)
    # pm.plot_image(file_certo, 1, 0)
    # pm.plot_image(file_cmems_multi, 2, 0)
    # pm.save_fig(file_out)
    # pm.close_plot()

def test():
    print('zona utm para coordenadas')
    from pyproj import CRS
    longitude = -8.85
    latitude = 42.59
    # Determine the UTM zone number
    zone_number = int((longitude + 180) / 6) + 1

    ZONE_LETTERS = "CDEFGHJKLMNPQRSTUVWXX"
    zone_letter = ZONE_LETTERS[int(latitude + 80) >> 3]

    # Determine the hemisphere
    hemisphere = 'north' if latitude >= 0 else 'south'

    # Create the UTM CRS
    utm_crs = CRS.from_dict({'proj': 'utm', 'zone': zone_number, 'south': hemisphere == 'south'})

    print(f'{zone_number}{zone_letter}')

def main():
    print('[INFO] Started map plot')
    #plot_doors()
    #plot_multiple()
    plot_figure_20190523()
    #test()

if __name__ == '__main__':
    main()