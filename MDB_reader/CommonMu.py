import os, pytz
import numpy as np
from netCDF4 import Dataset
from datetime import datetime as dt


class COMMON_MU:

    def __init__(self, options):
        print('[INFO] Common Match-ups class loaded')
        self.options = options

    def run(self):
        if self.options['type_cmu'] == 'multiple_resolution':
            self.run_multiple_resolution()

    def run_multiple_resolution(self):
        list_files = self.options['mdb_files']
        list_refs = self.options['mdb_refs']
        if list_refs is None:
            list_refs = [f'MDB_{idx}' for idx in range(len(list_files))]
        file_ref = list_files[0]
        sat_date_ref, lat_central, lon_central = self.get_insitu_lat_lon_nearest_extract_centre(file_ref)
        common_mu_arrays = None
        for ifile in range(1, len(list_files)):
            common_mu = self.find_common_mu_between_files(list_files[ifile],sat_date_ref,lat_central,lon_central)
            if ifile == 1:
                common_mu_arrays = common_mu.copy()
            else:
                common_mu_arrays = np.column_stack((common_mu_arrays, common_mu.copy()))
        common_mu_check = np.where(common_mu_arrays >= 0, 1, 0)
        common_mu_sum = np.sum(common_mu_check, axis=1)
        method = self.options['mdb_multiple_method']
        if method == 'any':
            valid = common_mu_sum > 0
        else:
            valid = common_mu_sum == len(list_files) - 1
        print(f'[INFO] Final common matchups using {method} method: {np.sum(valid)}')

    def get_insitu_lat_lon_nearest_extract_centre(self,file_ref):
        dref = Dataset(file_ref)
        lat_ref = dref.variables['satellite_latitude'][:]
        lon_ref = dref.variables['satellite_longitude'][:]
        sat_time_ref = dref.variables['satellite_time'][:]
        sat_date_ref = [dt.fromtimestamp(float(x)).astimezone(pytz.utc).replace(hour=0, minute=0, second=0,microsecond=0).timestamp() for x in sat_time_ref]
        sat_date_ref = np.array(sat_date_ref)
        lat_central, lon_central = get_lat_lon_central(lat_ref, lon_ref)
        insitu_lat = dref.variables['insitu_latitude'][:]
        insitu_lon = dref.variables['insitu_longitude'][:]
        insitu_spatial_index = dref.variables['insitu_spatial_index'][:]
        insitu_shape = insitu_lat.shape
        insitu_lat_f = insitu_lat.flatten()
        insitu_lon_f = insitu_lon.flatten()
        insitu_spatial_index = insitu_spatial_index.flatten()
        lat_central_r = np.tile(lat_central, (insitu_shape[1],))
        lon_central_r = np.tile(lon_central, (insitu_shape[1],))
        dif_sat_insitu = np.power((insitu_lat_f - lat_central_r), 2) + np.power((insitu_lon_f - lon_central_r), 2)
        dif_sat_insitu = np.ma.masked_where(insitu_spatial_index > 0, dif_sat_insitu)
        dif_sat_insitu = dif_sat_insitu.reshape(insitu_shape)
        indices_min = np.ma.argmin(dif_sat_insitu, axis=1)
        insitu_lat_ref = insitu_lat[np.arange(insitu_shape[0]), indices_min]
        insitu_lon_ref = insitu_lon[np.arange(insitu_shape[0]), indices_min]
        dref.close()

        return sat_date_ref,insitu_lat_ref,insitu_lon_ref


    def find_common_mu_between_files(self,file_work,sat_date_ref,lat_central,lon_central):

        dwork = Dataset(file_work)
        lat_work = dwork.variables['satellite_latitude'][:]
        lon_work = dwork.variables['satellite_longitude'][:]
        sat_time_work = dwork.variables['satellite_time'][:]
        sat_date_work = [dt.fromtimestamp(float(x)).astimezone(pytz.utc).replace(hour=0, minute=0, second=0,microsecond=0).timestamp() for x in sat_time_work]
        sat_date_work = np.array(sat_date_work)
        limits = get_geo_limits_central_pixel(lat_work, lon_work, 0.50)
        dwork.close()

        ndata = len(lat_central)
        common_mu = np.zeros(ndata)
        common_mu[:] = -1
        nmu = 0
        for ih in range(ndata):
            time_l = sat_date_ref[ih] == sat_date_work
            lat_l = np.logical_and(lat_central[ih] >= limits[:, 0], lat_central[ih] <= limits[:, 1])
            lon_l = np.logical_and(lon_central[ih] >= limits[:, 2], lon_central[ih] <= limits[:, 3])
            inside = np.logical_and(time_l, np.logical_and(lat_l, lon_l))
            ninside = np.count_nonzero(inside)
            if ninside == 1:
                nmu = nmu + 1
                common_mu[ih] = np.where(inside)[0][0]
                #print(ih,'->',common_mu[ih])
            elif ninside > 1:
                print('f[WARNING] There are more than one extract meeting the basic spatio/temporal condition')

        print(f'[INFO] Found {nmu} common match-ups from {ndata} extracts analysed')

        return common_mu


def get_lat_lon_central(lat_array, lon_array):
    nrows = lat_array.shape[1]
    ncols = lat_array.shape[2]
    cr = int(np.floor(nrows / 2))
    cl = int(np.floor(ncols / 2))
    return lat_array[:, cr, cl], lon_array[:, cr, cl]


#limits around the central pixel using a coverage value (0.5 dif_lat is the middle point between lat[r,c] and lat[r+1,c]
def get_geo_limits_central_pixel(lat_array, lon_array, coverage):
    ndata = lat_array.shape[0]
    nrows = lat_array.shape[1]
    ncols = lat_array.shape[2]
    cr = int(np.floor(nrows / 2))
    cl = int(np.floor(ncols / 2))
    limits = np.zeros((ndata, 4))
    dlat1 = np.abs(lat_array[:, cr, cl] - lat_array[:, cr - 1, cl]) * coverage
    dlat2 = np.abs(lat_array[:, cr, cl] - lat_array[:, cr + 1, cl]) * coverage
    dif_lat = (dlat1 + dlat2) / 2
    dlon1 = np.abs(lon_array[:, cr, cl] - lon_array[:, cr, cl - 1]) * coverage
    dlon2 = np.abs(lon_array[:, cr, cl] - lon_array[:, cr, cl + 1]) * coverage
    dif_lon = (dlon1 + dlon2) / 2
    limits[:, 0] = lat_array[:, cr, cl] - dif_lat
    limits[:, 1] = lat_array[:, cr, cl] + dif_lat
    limits[:, 2] = lon_array[:, cr, cl] - dif_lon
    limits[:, 3] = lon_array[:, cr, cl] + dif_lon
    return limits
