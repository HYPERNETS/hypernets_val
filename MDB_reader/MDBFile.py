from datetime import datetime
from datetime import timedelta

from netCDF4 import Dataset
import numpy as np

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
        print(['LOIS: Starting MDBFile...'])
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


    # Funcion to load data from a specific MU
    def load_mu_data(self, index_mu):
        if not self.VALID:
            return False
        if index_mu < 0 or index_mu >= self.n_mu_total:
            print('Not valid index_mu')
            return False

        self.index_mu = index_mu
        self.mu_sat_time = self.sat_times[index_mu]

        #Dimensions
        central_r = int(np.floor(len(self.dimensions['rows']) / 2))
        central_c = int(np.floor(len(self.dimensions['columns']) / 2))
        r_s = central_r - int(np.floor(self.window_size / 2))  # starting row
        r_e = central_r + int(np.floor(self.window_size / 2)) + 1  # ending row
        c_s = central_c - int(np.floor(self.window_size / 2))  # starting col
        c_e = central_c + int(np.floor(self.window_size / 2)) + 1  # ending col

        # flags mask
        flagging = flag.Class_Flags_OLCI(self.nc.satellite_WQSF_flag_masks, self.nc.satellite_WQSF_flag_meanings)
        satellite_WQSF = self.variables['satellite_WQSF'][r_s:r_e, c_s:c_e]
        if (str(self.flag_list[0])) != 'None':
            mask = flagging.Mask(satellite_WQSF, self.flag_list)
            mask[np.where(mask != 0)] = 1
        else:
            mask = np.full(satellite_WQSF.shape, 0, dtype=int)
        print(mask)


    ##PLOT FUNCTIONS
    def plot_spectra(self, index_mu, wavelenghts, ydata, platform, sat_time, insitu_time):
        print('plotea un espectro')
        # legend = [f"S3{platform}-{sat_time.strftime('%d-%m-%Y %H:%M')}",
        #           f"Hypstar-{insitu_time.strftime('%d-%m-%Y %H:%M')}"]
        # row_names = []
        # for w in wavelenghts:
        #     wn = f'{w}'
        #     row_names.append(wn)
        # df = pd.DataFrame(np.transpose(ydata), columns=legend, index=row_names)
        # title = f"Match-up #{index_mu} {sat_time.strftime('%d-%m-%Y')}"
        #
        # plt.figure()
        # df.plot(lw=2, marker='.', markersize=10)
        # rindex = list(range(16))
        # plt.xticks(rindex, row_names, rotation=45, fontsize=12)
        # plt.xlabel('Wavelength(nm)', fontsize=14)
        # plt.ylabel('R$_{rs}$[1/sr]', fontsize=14)
        # plt.title(title, fontsize=16)
        # plt.gcf().subplots_adjust(bottom=0.18)
        # plt.gcf().subplots_adjust(left=0.20)
        #
        # ofname = f"Spectra_S3{platform}_{sat_time.strftime('%Y%m%d')}.jpg"
        # ofname = os.path.join(path_out, ofname)
        # print('Name out:', ofname)
        #
        # plt.savefig(ofname, dpi=300)
        # plt.close()
        #
        # return ofname
