import os
import sys
from datetime import datetime
from datetime import timedelta

from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pandas as pd

code_home = os.path.abspath('../')
sys.path.append(code_home)
import COMMON.Class_Flags_OLCI as flag

DIMENSION_NAMES = ['satellite_id', 'rows', 'columns', 'satellite_bands', 'insitu_id', 'insitu_original_bands']
VAR_NAMES = ['satellite_time', 'satellite_PDU', 'satellite_latitude', 'satellite_longitude', 'satellite_bands',
             'satellite_Rrs', 'satellite_AOT_0865p50', 'satellite_WQSF', 'chl_oc4me', 'insitu_time', 'insitu_filename',
             'insitu_filepath', 'insitu_original_bands', 'insitu_rhow', 'time_difference']

ATRIB_NAMES = ['creation_time', 'satellite', 'platform', 'sensor', 'description', 'satellite_aco_processor',
               'satellite_proc_version', 'insitu_site_name', 'insitu_lat', 'insitu_lon', 'time_diff']


class MDBFile:

    def __init__(self, file_path):
        global DIMENSION_NAMES
        global VAR_NAMES
        global ATRIB_NAMES
        self.file_path = file_path
        self.VALID = True
        print(['Starting MDBFile...'])
        try:
            self.nc = Dataset(file_path)
            self.variables = self.nc.variables
            self.dimensions = self.nc.dimensions
            self.VALID = self.check_structure()
        # except:
        except Exception as e:
            self.VALID = False
            print(f'Exception: {e}')

        if self.VALID:
            self.n_mu_total = len(self.dimensions['satellite_id'])
            self.sat_times = []
            for st in self.variables['satellite_time']:
                self.sat_times.append(datetime(1970, 1, 1) + timedelta(seconds=int(st)))
            self.start_date = self.sat_times[0]
            self.end_date = self.sat_times[-1]
            self.insitu_bands = self.nc.variables['insitu_original_bands'][:]  # insitu_bands(insitu_bands)
            self.satellite_bands = self.nc.variables['satellite_bands'][:]
            self.n_insitu_bands = len(self.insitu_bands)
            self.n_satellite_bands = len(self.satellite_bands)

        # Sat filtering options
        self.flag_list = ''
        self.window_size = 3
        self.valid_min_pixels = 1
        self.delta_t = 7200
        self.set_default_sat_filtering_options()

        # Variables defining a specific MU. To load a MU, uses load_mu_data
        self.index_mu = 0
        self.mu_sat_time = []
        self.mu_insitu_time = []
        self.mu_curr_sat_rrs_mean = []
        self.mu_curr_ins_rrs = []

    # Checking atrib, var and dimensions names
    def check_structure(self):
        check_var = True
        check_dim = True
        check_atrib = True
        for var in VAR_NAMES:
            if var not in self.variables:
                check_var = False
        for dim in DIMENSION_NAMES:
            if dim not in self.dimensions:
                check_dim = False
        for atrib in ATRIB_NAMES:
            if atrib not in self.nc.ncattrs():
                check_atrib = False
        if check_var == False or check_dim == False or check_atrib == False:
            return False
        return True

    # Set default sat filtering options
    def set_default_sat_filtering_options(self):
        flag_list = 'CLOUD,CLOUD_AMBIGUOUS,CLOUD_MARGIN,INVALID,COSMETIC,SATURATED,SUSPECT,HISOLZEN,HIGHGLINT,SNOW_ICE,AC_FAIL,WHITECAPS,RWNEG_O2,RWNEG_O3,RWNEG_O4,RWNEG_O5,RWNEG_O6,RWNEG_O7,RWNEG_O8'
        flag_list = flag_list.replace(" ", "")
        self.flag_list = str.split(flag_list, ',')
        self.window_size = 3
        self.valid_min_pixels = 1
        self.delta_t = 7200

    # Funcion to load data from a specific MU
    def load_mu_data(self, index_mu):
        if not self.VALID:
            return False
        if index_mu < 0 or index_mu >= self.n_mu_total:
            print('Not valid index_mu')
            return False

        self.index_mu = index_mu
        self.mu_sat_time = self.sat_times[index_mu]
        time_difference = self.variables['time_difference'][index_mu]
        ins_time_index = np.argmin(np.abs(time_difference))
        ins_time = self.variables['insitu_time'][index_mu][ins_time_index]
        self.mu_insitu_time = datetime(1970, 1, 1) + timedelta(seconds=int(ins_time))

        # Dimensions
        central_r = int(np.floor(len(self.dimensions['rows']) / 2))
        central_c = int(np.floor(len(self.dimensions['columns']) / 2))
        r_s = central_r - int(np.floor(self.window_size / 2))  # starting row
        r_e = central_r + int(np.floor(self.window_size / 2)) + 1  # ending row
        c_s = central_c - int(np.floor(self.window_size / 2))  # starting col
        c_e = central_c + int(np.floor(self.window_size / 2)) + 1  # ending col

        # flags mask
        flagging = flag.Class_Flags_OLCI(self.variables['satellite_WQSF'].flag_masks,
                                         self.variables['satellite_WQSF'].flag_meanings)
        satellite_WQSF = self.variables['satellite_WQSF'][index_mu, r_s:r_e, c_s:c_e]
        if (str(self.flag_list[0])) != 'None':
            mask = flagging.Mask(satellite_WQSF, self.flag_list)
            mask[np.where(mask != 0)] = 1
        else:
            mask = np.full(satellite_WQSF.shape, 0, dtype=int)
        NGP = np.power(self.window_size, 2) - np.sum(mask)  # Number Good Pixels excluding flagged pixels
        land = flagging.Mask(satellite_WQSF, (['LAND']))
        inland_w = flagging.Mask(satellite_WQSF, (['INLAND_WATER']))
        land[np.where(inland_w > 0)] = 0
        NTP = np.power(self.window_size, 2) - np.sum(land, axis=(0, 1))
        # conditions minimum valid pixels
        if float(self.valid_min_pixels) == 1:  # 100% like Zibordi
            cond_min_valid_pxs = NGP == np.power(self.window_size, 2)
        else:
            cond_min_valid_pxs = NGP > (float(self.valid_min_pixels) * NTP + 1)

        if np.abs(time_difference[ins_time_index]) < self.delta_t and cond_min_valid_pxs:
            satellite_rrs = self.variables['satellite_Rrs'][index_mu]
            insitu_rrs = self.variables['insitu_rhow'][index_mu]
            insitu_rrs = insitu_rrs / np.pi  # transform from rhow to Rr
            self.mu_curr_ins_rrs = []
            self.mu_curr_sat_rrs_mean = []
            for sat_band_index in range(0, self.n_satellite_bands):
                wl = self.satellite_bands[sat_band_index]
                curr_sat_box = satellite_rrs[sat_band_index, r_s:r_e, c_s:c_e]
                ins_band_index = np.argmin(np.abs(wl - self.insitu_bands))
                # if debug:
                #     print(f'Closest in situ band to sat band {wl}: {insitu_bands[ins_band_index]} nm')
                curr_sat_box = np.ma.masked_where(curr_sat_box == -999, curr_sat_box)
                curr_sat_box = np.ma.masked_invalid(curr_sat_box)
                curr_sat_box = np.ma.masked_array(curr_sat_box, mask)
                numberValid = np.ma.count(curr_sat_box)
                if float(self.valid_min_pixels) == 1:  # 100% like Zibordi
                    cond_min_valid_pxs = numberValid == np.power(self.window_size, 2)
                else:
                    cond_min_valid_pxs = numberValid > (float(self.valid_min_pixels) * NTP + 1)

                if cond_min_valid_pxs:
                    curr_sat_box_mean = curr_sat_box.mean()
                    self.mu_curr_sat_rrs_mean.append(curr_sat_box_mean)
                    ins_value = insitu_rrs[ins_band_index, ins_time_index]
                    self.mu_curr_ins_rrs.append(ins_value)

            # print(self.mu_curr_sat_rrs_mean)
            # print(self.mu_curr_ins_rrs)

    ##PLOT FUNCTIONS
    def plot_spectra(self, hfig):
        print('Plot spectra...')
        legend = [f"S3{self.nc.platform}-{self.mu_sat_time.strftime('%d-%m-%Y %H:%M')}",
                  f"Hypstar-{self.mu_insitu_time.strftime('%d-%m-%Y %H:%M')}"]
        row_names = []
        for w in self.satellite_bands:
            wn = f'{w}'
            row_names.append(wn)
        ydata = np.array([self.mu_curr_sat_rrs_mean, self.mu_curr_ins_rrs])
        df = pd.DataFrame(np.transpose(ydata), columns=legend, index=row_names)
        title = f"Match-up #{self.index_mu} {self.mu_sat_time.strftime('%d-%m-%Y')}"

        if hfig is None:
            hfig = plt.figure()
        df.plot(lw=2, marker='.', markersize=10)
        rindex = list(range(16))
        plt.xticks(rindex, row_names, rotation=45, fontsize=12)
        plt.xlabel('Wavelength(nm)', fontsize=14)
        plt.ylabel('R$_{rs}$[1/sr]', fontsize=14)
        plt.title(title, fontsize=16)
        plt.gcf().subplots_adjust(bottom=0.18)
        plt.gcf().subplots_adjust(left=0.20)
        plt.gcf().canvas.draw()

        return hfig

        #
        # ofname = f"Spectra_S3{platform}_{sat_time.strftime('%Y%m%d')}.jpg"
        # ofname = os.path.join(path_out, ofname)
        # print('Name out:', ofname)
        #
        # plt.savefig(ofname, dpi=300)
        # plt.close()
        #
        # return ofname
