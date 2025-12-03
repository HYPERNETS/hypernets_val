import os, pytz, sys, __init__
import numpy as np
from netCDF4 import Dataset
from datetime import datetime as dt
from datetime import timezone
code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
import COMMON.common_functions as cfs

class COMMON_MU:

    def __init__(self, options):
        print('[INFO] Common Match-ups class loaded')
        self.options = options

    def run(self):
        if self.options['type_cmu'] == 'multiple_resolution':
            self.run_multiple_resolution()

        if self.options['type_cmu'] == 'single_mdb':
            self.run_single_mdb()

    def run_single_mdb(self):
        input_file = self.options['input_file']
        flag_among = self.options['flag_among']
        flag_among_values = self.options['flag_among_values']
        reference = self.options['reference']
        flag_array = None
        sat_time_array  = None
        sat_lat_array = None
        sat_lon_array = None
        rc_center = None

        dataset = Dataset(input_file,'r')
        if flag_among in dataset.variables:
            flag_array = dataset.variables[flag_among][:]
        if 'satellite_time' in dataset.variables:
            sat_time_array = dataset.variables['satellite_time'][:]
        if reference=='sat_time_lat_lon':
            if 'satellite_latitude' in dataset.variables:
                sat_lat_array = dataset.variables['satellite_latitude'][:]
                rc_center = int(np.floor(sat_lat_array.shape[1] / 2))

            if 'satellite_longitude' in dataset.variables:
                sat_lon_array = dataset.variables['satellite_longitude'][:]

        if flag_among_values is None:
            flag_among_values = np.unique(flag_array).tolist()

        dataset.close()

        if flag_array is None or sat_time_array is None:
            return

        if reference == 'sat_time_lat_lon' and (sat_lat_array is None or sat_lon_array is None):
            return


        check_mu = {}
        nflag_values = len(flag_array)

        for idx in range(nflag_values):
            refidx = dt.fromtimestamp(np.int64(sat_time_array[idx])).astimezone(timezone.utc).strftime('%Y%m%dT%H%M%S')
            if reference.startswith('sat_date'):
                refidx = dt.fromtimestamp(sat_time_array[idx]).astimezone(timezone.utc).strftime('%Y%m%d')
            flag_value = flag_array[idx]
            if flag_value in flag_among_values:
                if not refidx in check_mu:
                    check_mu[refidx] = {flag_value:[idx]}
                else:
                    if flag_value in check_mu[refidx]:
                        check_mu[refidx][flag_value].append(idx)
                    else:
                        check_mu[refidx][flag_value] = [idx]

        flag_value_ref = flag_among_values[0]
        flag_value_comp = flag_among_values[1:]
        index_used = [False]*nflag_values
        common_mu_array = np.zeros(nflag_values)
        index_common = 1
        for idx in range(nflag_values):
            flag_value = flag_array[idx]
            if flag_value!=flag_value_ref:
                continue
            refidx = dt.fromtimestamp(np.int64(sat_time_array[idx])).astimezone(timezone.utc).strftime('%Y%m%dT%H%M%S')
            if reference.startswith('sat_date'):
                refidx = dt.fromtimestamp(sat_time_array[idx]).astimezone(timezone.utc).strftime('%Y%m%d')
            ncomp = 0
            for fv in flag_value_comp:
                if fv in check_mu[refidx]:
                    indices_fv = check_mu[refidx][fv]
                    for index_fv in indices_fv:
                        if index_used[index_fv]:
                            continue
                        if cfs.is_central_pixel(sat_lat_array[index_fv,:,:],sat_lon_array[index_fv,:,:],sat_lat_array[idx,rc_center,rc_center],sat_lon_array[idx,rc_center,rc_center]):
                            common_mu_array[idx] = index_common
                            common_mu_array[index_fv] = index_common
                            index_used[idx] = True
                            index_used[index_fv] = True
                            ncomp = ncomp + 1
                            break
            if ncomp==len(flag_value_comp):
                index_common = index_common + 1


        print(f'[INFO] Number of common match-ups: {np.max(common_mu_array)} / {nflag_values}')


        return common_mu_array



    def get_mu_variable(self,common_var):
        if self.options['mu_valid_variable'] is None:
            print(f'[ERROR] mu_valid_variable is required if create_mu_variable is True')
            return None
        input_file = self.options['input_file']
        dataset = Dataset(input_file, 'r')
        if not common_var in dataset.variables:
            print(f'[ERROR] {common_var} variable is required for getting the common mu valid array but is not available in {os.path.basename(input_file)}')
            dataset.close()
            return None
        common_array = dataset.variables[common_var][:]
        mu_valid_var = self.options['mu_valid_variable']
        if not mu_valid_var in dataset.variables:
            print(f'[ERROR] {mu_valid_var} variable is required for getting the common mu valid array  but is not available in {os.path.basename(input_file)}')
            dataset.close()
            return None
        mu_valid_array = dataset.variables[mu_valid_var][:]
        dataset.close()

        output_array = np.zeros(mu_valid_array.shape).astype(np.int16)
        used = np.zeros(mu_valid_array.shape)
        output_array[:]=-1
        ndata = output_array.shape[0]
        nmu_x_index = -1
        nmu_common = 0
        for idx in range(ndata):
            if used[idx]==1:
                continue
            index_common = common_array[idx]
            if index_common<=0:
                used[idx]=1
                continue
            indices = np.where(common_array==index_common)
            nmu_x_index_here = len(indices[0])
            if nmu_x_index==-1:
                nmu_x_index = nmu_x_index_here
            else:
                if nmu_x_index!=nmu_x_index_here:
                    print(f'[WARNING] Index common {index_common} show a different number of match-ups')
                    continue
            used[indices]=1
            valid_indices = mu_valid_array[indices]
            if np.sum(valid_indices)==nmu_x_index:
                output_array[indices] = index_common
                nmu_common = nmu_common + 1


        print(f'[INFO] Number of valid match-ups: {nmu_common}')
        return output_array


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
