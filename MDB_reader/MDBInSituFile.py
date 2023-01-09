from netCDF4 import Dataset
from datetime import datetime
from datetime import timedelta
from QC_SAT import QC_SAT
import sys
import os
import numpy as np
import pandas as pd

DIMENSION_NAMES = ['satellite_id', 'rows', 'columns', 'satellite_bands']
VAR_NAMES = ['satellite_time', 'satellite_latitude', 'satellite_longitude', 'satellite_bands',
             'satellite_Rrs', 'insitu_time', 'time_difference']

ATRIB_NAMES = ['creation_time', 'satellite', 'platform', 'sensor', 'description', 'satellite_aco_processor',
               'satellite_proc_version', 'insitu_site_name', 'insitu_lat', 'insitu_lon']


class MDBInSituFile:

    def __init__(self, file_path):
        global DIMENSION_NAMES
        global VAR_NAMES
        global ATRIB_NAMES

        self.file_path = file_path
        self.VALID = True
        self.balmlp = None
        file_name = file_path.split('/')[-1]

        #self.delta_t = 7200 #2 hours
        self.delta_t = 10800 #3 hours


        # Variables defining a specific MU. To load a MU, uses load_mu_data
        self.index_mu = 0
        self.mu_sat_time = []
        self.mu_insitu_time = []
        self.satellite_rrs = []
        self.mu_curr_sat_rrs_mean = []

        self.df_validation = None
        self.df_validation_valid = None
        self.sat_retrievals = {}

        print(f'[INFO]Starting MDBInSituFile: {file_name}')
        try:
            self.nc = Dataset(file_path)
            self.variables = self.nc.variables
            self.dimensions = self.nc.dimensions
            self.flag_band_name = 'satellite_WQSF'
            self.VALID = self.check_structure()
            self.VALID = True
        # except:
        except Exception as e:
            self.VALID = False
            print(f'Exception starting MDB File: {e}')

        if not self.VALID:
            print(f'MDB File: {file_path} is not valid')

        if self.VALID:
            self.n_mu_total = len(self.dimensions['satellite_id'])
            print('[INFO ]Total mu: ', self.n_mu_total)
            self.sat_times = []
            for st in self.variables['satellite_time']:
                #self.sat_times.append(datetime(1970, 1, 1) + timedelta(seconds=int(st)))
                sat_time_here = datetime.fromtimestamp(float(st))
                self.sat_times.append(sat_time_here)
            # self.pdus = []
            # for st in self.variables['satellite_PDU']:
            #     self.pdus.append(st)
            self.insitu_times = []
            self.insitu_lats = []
            self.insitu_lons = []
            if 'insitu_time' in self.variables:
                for st in self.variables['insitu_time']:
                    self.insitu_times.append(datetime(1970, 1, 1) + timedelta(seconds=float(st)))
            if 'insitu_latitude' in self.variables:
                for st in self.variables['insitu_latitude']:
                    self.insitu_lats.append(float(st))
            if 'insitu_longitude' in self.variables:
                for st in self.variables['insitu_longitude']:
                    self.insitu_lons.append(float(st))


            self.start_date = self.sat_times[0]
            self.end_date = self.sat_times[-1]
            self.satellite_bands = self.nc.variables['satellite_bands'][:]
            self.n_satellite_bands = len(self.satellite_bands)
            self.info = {}
            atribs = self.nc.ncattrs()
            for at in atribs:
                if at == 'satellite_aco_processor':
                    self.info[at] = self.nc.getncattr(at).upper()
                else:
                    self.info[at] = self.nc.getncattr(at)
            if self.info['satellite_aco_processor'] == 'ATMOSPHERIC CORRECTION PROCESSOR: XXX':
                self.info['satellite_aco_processor'] = 'STANDARD'

            if self.info['satellite_aco_processor'] == 'CLIMATE CHANGE INITIATIVE - EUROPEAN SPACE AGENCY':
                self.info['satellite_aco_processor'] = 'CCI'

            self.wlref = self.satellite_bands
            self.wlref_sat_indices = list(range(len(self.satellite_bands)))
            # self.set_wl_ref_insitu_indices()

            # SATELLITE QUALITY CONTROL
            if self.nc.satellite_aco_processor == 'ACOLITE' or self.nc.satellite_aco_processor == 'Climate Change Initiative - European Space Agency':
                self.qc_sat = QC_SAT(self.variables['satellite_Rrs'], self.satellite_bands, None,
                                     self.info['satellite_aco_processor'])
            else:
                if self.flag_band_name in self.variables:
                    self.qc_sat = QC_SAT(self.variables['satellite_Rrs'], self.satellite_bands,self.variables[self.flag_band_name], self.info['satellite_aco_processor'])
                else:
                    self.qc_sat = QC_SAT(self.variables['satellite_Rrs'], self.satellite_bands,None, self.info['satellite_aco_processor'])

            # Variables to make validation...
            self.insitu_varnames = []
            self.col_names = ['Index_MU', 'Sat_Time', 'Ins_Time', 'Ins_Lat','Ins_Lon','Time_Diff', 'Valid']
            for var in self.nc.variables:
                if var.startswith('insitu_'):
                    name = var.split('_')[1]
                    if name.startswith('lat') or name.startswith('long') or name.startswith('time'):
                        continue
                    else:
                        self.insitu_varnames.append(var)

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

        if not self.flag_band_name in self.variables:
            if self.nc.satellite_aco_processor.upper() == 'POLYMER' and 'satellite_bitmask' in self.variables:
                self.flag_band_name = 'satellite_bitmask'

            if self.nc.satellite_aco_processor.upper() == 'C2RCC' and 'satellite_c2rcc_flags' in self.variables:
                self.flag_band_name = 'satellite_c2rcc_flags'

            if self.nc.satellite_aco_processor.upper() == 'FUB' and 'satellite_quality_flags':
                self.flag_band_name = 'satellite_quality_flags'

            if self.nc.satellite_aco_processor.upper() == 'ACOLITE':
                self.flag_band_name = 'NONE'

        return True

    def set_wl_ref(self, wllist):
        self.wlref = wllist
        self.wlref_sat_indices = []
        self.qc_sat.wl_ref = wllist
        for wl in wllist:
            index = np.argmin(np.abs(wl - self.satellite_bands))
            self.wlref_sat_indices.append(index)

    def set_wlsatlist_aswlref(self, wlsatlist):
        wllist = self.qc_sat.get_wl_sat_list_from_wlreflist(wlsatlist)
        print('[INFO] All sat available wavelenghts: ', self.satellite_bands)
        print('[INFO] Set wavelenghts for MU extraction to: ', wllist)
        if len(wllist) == len(wlsatlist):
            self.set_wl_ref(wllist)

    # start baltic chla retrievals, name_insitu_var must be the variable containing in situ chla
    def start_baltic_chla(self, path_code, name_insitu_var):
        sys.path.append(path_code)
        try:
            from balticmlp.balmlpensemble import BalMLP
            path_data = os.path.join(path_code, 'balticmlp')
            self.balmlp = BalMLP(path_data)
            print(f'[INFO] Loaded BalMLP module')
            if name_insitu_var is not None:
                self.sat_retrievals[name_insitu_var] = {
                    'function': 'baltic_chla',
                    'satvarnames': ['wgt3', 'wgt4', 'wgt5_412', 'wgt5_670', 'wgt6', 'chl3', 'chl4', 'chl5_412',
                                    'chl5_670', 'chl6', 'ens3', 'ens4']
                }

            return True
        except ModuleNotFoundError:
            print(f'[ERROR] BalMLP module is not found')
            return False

    def load_mu_datav2(self, index_mu):

        is_mu_valid = False
        status = 0

        if not self.VALID:
            status = -1  # 'NO VALID MDB FILE'
            return is_mu_valid, status

        if index_mu < 0 or index_mu >= self.n_mu_total:
            status = -2  # f'NO VALID MATCH-UP INDEX:{index_mu}'
            return is_mu_valid, status

        # Index match-up
        self.index_mu = index_mu

        # Sat rrs
        self.satellite_rrs = self.variables['satellite_Rrs'][index_mu]

        # Sat and instrument time
        self.mu_sat_time = self.sat_times[index_mu]
        if self.info['satellite_aco_processor'] == 'CCI':
            #self.mu_sat_time = self.mu_sat_time.replace(hour=11)
            self.mu_sat_time = self.mu_sat_time.replace(hour=13)

        self.mu_insitu_time = self.insitu_times[index_mu]

        #time_diff = float(self.nc.variables['time_difference'][index_mu])
        time_diff = abs((self.mu_sat_time-self.mu_insitu_time).total_seconds())
        time_condition = time_diff <= self.delta_t

        if not time_condition:
            status = -3  # f'IN SITU DATA OUT OF TIME WINDOW'
            return is_mu_valid, status

        cond_min_pixels, cond_stats, valid_mu, sat_values = self.qc_sat.get_match_up_values(index_mu)
        if not valid_mu:
            status = -6  # f'NO VALID SAT DATA'
            return is_mu_valid, status

        self.mu_curr_sat_rrs_mean = sat_values
        is_mu_valid = True
        return is_mu_valid, status

    def prepare_df_validation(self):
        print('[INFO] Preparing DF for validation...')
        ntot = self.n_mu_total

        for wl in self.wlref:
            wls =  '{0:3.0f}'.format(wl)
            self.col_names.append(f'{wls}')

        for varname in self.insitu_varnames:
            self.col_names.append(varname)
            if varname in self.sat_retrievals.keys():
                for satvarname in self.sat_retrievals[varname]['satvarnames']:
                    self.col_names.append(satvarname)
        print(self.col_names)
        self.df_validation = pd.DataFrame(columns=self.col_names, index=list(range(ntot)))
        nmu_valid = 0

        for index_mu in range(self.n_mu_total):
            if index_mu % 100 == 0:
                print(f'[INFO] MU: {index_mu} of {self.n_mu_total}')
            mu_valid, status = self.load_mu_datav2(index_mu)
            ##MANUAL
            if index_mu==112 or index_mu==288:
                mu_valid = False

            if mu_valid:
                nmu_valid = nmu_valid + 1

            # self.col_names = ['Index_MU', 'Sat_Time', 'Ins_Time', 'Time_Diff', 'Valid']
            time_diff = round(abs((self.mu_sat_time - self.mu_insitu_time).total_seconds() / 3600), 2)
            row = {
                'Index_MU': [index_mu],
                'Sat_Time': [self.mu_sat_time.strftime('%Y-%m-%d %H:%M')],
                'Ins_Time': [self.mu_insitu_time.strftime('%Y-%m-%d %H:%M')],
                'Time_Diff': [time_diff],
                'Valid': [mu_valid]  # [self.mu_valid_bands[sat_band_index]]
            }
            for idx in range(len(self.wlref)):
                wl = self.wlref[idx]
                wls = '{0:3.0f}'.format(wl)
                rrsvalue = self.mu_curr_sat_rrs_mean[idx]
                row[wls] = [rrsvalue]
            for varname in self.insitu_varnames:
                varvalue = self.nc.variables[varname][index_mu]
                row[varname] = [varvalue]
                if varname in self.sat_retrievals.keys():
                    satvarnames = self.sat_retrievals[varname]['satvarnames']
                    if mu_valid:
                        values = self.compute_sat_retrieval(self.sat_retrievals[varname]['function'])
                    else:
                        values = [-999] * len(satvarnames)
                    for idx in range(len(satvarnames)):
                        row[satvarnames[idx]] = [values[idx]]

            self.df_validation.iloc[index_mu] = pd.DataFrame.from_dict(row)

        self.df_validation_valid = self.df_validation[self.df_validation['Valid']][:]

        print(f'[INFO]# total match-ups: {self.n_mu_total} Valid: {nmu_valid})')
        return nmu_valid

    def prepare_df_extracts(self):
        print('[INFO] Preparing DF extracts...')
        ntot = self.n_mu_total
        for wl in self.wlref:
            wls = '{0:3.0f}'.format(wl)
            self.col_names.append(f'{wls}')
        self.df_validation = pd.DataFrame(columns=self.col_names, index=list(range(ntot)))
        nmu_valid = 0

        for index_mu in range(self.n_mu_total):
            if index_mu % 100 == 0:
                print(f'[INFO] MU: {index_mu} of {self.n_mu_total}')
            cond_min_pixels, cond_stats, valid_mu, values = self.qc_sat.get_match_up_values(index_mu)
            if valid_mu:
                nmu_valid = nmu_valid + 1
            self.mu_sat_time = self.sat_times[index_mu]
            self.mu_sat_time = self.mu_sat_time.replace(hour=12,minute=0,second=0)
            insitu_time = self.insitu_times[index_mu]
            time_dif = abs((self.mu_sat_time-insitu_time).total_seconds())/3600
            ins_lat = self.insitu_lats[index_mu]
            ins_lon = self.insitu_lons[index_mu]
            row = {
                'Index_MU': [index_mu],
                'Sat_Time': [self.mu_sat_time.strftime('%Y-%m-%d %H:%M')],
                'Ins_Time': [insitu_time.strftime('%Y-%m-%d %H:%M')],
                'Ins_Lat': [ins_lat],
                'Ins_Lon': [ins_lon],
                'Time_Diff': [time_dif],
                'Valid': [valid_mu],  # [self.mu_valid_bands[sat_band_index]]
            }
            for idx in range(len(self.wlref)):
                wl = self.wlref[idx]
                wls = '{0:3.0f}'.format(wl)
                rrsvalue = values[idx]
                row[wls] = [rrsvalue]
            self.df_validation.iloc[index_mu] = pd.DataFrame.from_dict(row)

        self.df_validation_valid = self.df_validation[self.df_validation['Valid']][:]
        print(f'[INFO]# total match-ups: {self.n_mu_total} Valid: {nmu_valid})')
        return nmu_valid

    def compute_sat_retrieval(self, function):
        if function == 'baltic_chla':
            return self.compute_baltic_chla_mu(-1)

        return None

    # index_mu==-1 (or <0) to no load data
    def compute_baltic_chla_mu(self, index_mu):
        if index_mu >= 0:
            is_mu_valid, load_info = self.load_mu_datav2(index_mu)
            if not is_mu_valid:
                return None

        rrs = np.zeros((1, len(self.mu_curr_sat_rrs_mean)))
        rrs[0, :] = np.array(self.mu_curr_sat_rrs_mean)

        self.balmlp.compute_chla_ensemble(rrs)
        values = self.balmlp.get_values()
        return values

    def add_baltic_chla(self,ofile,ncdest,closedest):
        if ncdest is None:
            ncdest = self.copy_nc(ofile)
        ncdest.createDimension('bal_chla_id', 12)
        chla_var = ncdest.createVariable('satellite_chla_bal_mlp', 'f4', ('satellite_id','bal_chla_id','rows','columns'),
                                                       fill_value=-999, zlib=True, complevel=6)
        if closedest:
            ncdest.close()

    def copy_nc(self, ofile):

        dst = Dataset(ofile, 'w', format='NETCDF4')

        # copy global attributes all at once via dictionary
        dst.setncatts(self.nc.__dict__)

        # copy dimensions
        for name, dimension in self.nc.dimensions.items():
            dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

        # copy all file data except for the excluded
        for name, variable in self.nc.variables.items():
            dst.createVariable(name, variable.datatype, variable.dimensions)
             # copy variable attributes all at once via dictionary
            dst[name].setncatts(self.nc[name].__dict__)
            dst[name][:] = self.nc[name][:]

        return dst

    def close_dataset(self):
        self.nc.close()
