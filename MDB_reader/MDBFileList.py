import datetime
import os.path

from MDBFile import MDBFile
from MDBPlot import MDBPlot
from QCBase import QCBase
import pandas as pd
from datetime import datetime as dt


class MDBFileList():
    def __init__(self):
        self.mdb_list = {}
        self.df_validation = None  # pd.DataFrame(columns=col_names)
        self.mu_dates = {}
        # self.date_list = {}

        self.wlref = None

        self.qc_base = QCBase()

    def set_wlsatlist_as_ref(self):
         for fmdb in self.mdb_list:
             if self.mdb_list[fmdb]['include']:
                mfile = MDBFile(self.mdb_list[fmdb]['path'])
                self.wlref = mfile.wlref
                
    def set_wlsatlist_from_wlreflist_asref(self,wlsatlist):
        for fmdb in self.mdb_list:
            if self.mdb_list[fmdb]['include']:
                mfile = MDBFile(self.mdb_list[fmdb]['path'])
                mfile.set_wlsatlist_aswlref(wlsatlist)
                self.wlref = mfile.wlref

    def set_wl_ref(self, wllist):
        self.wlref = wllist



    def add_mdb_file(self, path_mdb):
        mfile = MDBFile(path_mdb)
        if not mfile.VALID:
            return False
        if self.wlref is None:
            self.wlref = mfile.wlref

        name = path_mdb.split('/')[-1]
        self.mdb_list[name] = {
            'path': path_mdb,
            'include': True
        }
        for key in mfile.info:
            self.mdb_list[name][key] = mfile.info[key]

        if self.mdb_list[name]['satellite_aco_processor'] == 'Atmospheric Correction processor: xxx':
            self.mdb_list[name]['satellite_aco_processor'] = 'STANDARD'
        if self.mdb_list[name]['satellite_aco_processor'] == 'CLIMATE CHANGE INITIATIVE - EUROPEAN SPACE AGENCY':
            self.mdb_list[name]['satellite_aco_processor'] = 'CCI'

    def prepare_df_for_validation(self):
        for fmdb in self.mdb_list:
            if self.mdb_list[fmdb]['include']:
                mfile = MDBFile(self.mdb_list[fmdb]['path'])
                if not self.wlref is None:
                    mfile.set_wl_ref(self.wlref)

                if fmdb.find('CCI')>0:
                    mfile.set_hour_sat_time(11, 0) ##ONLY FOR CCI
                if fmdb.find('MULTI')>0:
                    mfile.set_hour_sat_time(11,0)
                if fmdb.find('OLCI-L3')>0:
                    mfile.set_hour_sat_time(11,0)

                #mfile.qc_sat.set_qc_from_qcbase(self.qc_base)
                mfile.qc_sat.set_eumetsat_defaults(3)
                if 'satellite_pixel_classif_flags' in mfile.nc.variables:
                    idepix_flag = mfile.nc.variables['satellite_pixel_classif_flags']
                    mfile.qc_sat.set_idepix_as_flag(idepix_flag)

                #mfile.qc_sat.add_band_statistics(-1, 400, 'avg', True, 0.003, 'greater')
                mfile.qc_sat.add_band_statistics(-1, 400, 'avg', True, 0, 'lower')

                mfile.qc_insitu.set_wllist_using_wlref(mfile.wlref)
                mfile.qc_insitu.check_indices_by_mu = True
                mfile.qc_insitu.apply_band_shift = True
                mfile.qc_insitu.set_thershold(0, None, 0, 600)
                mfile.qc_insitu.set_thershold(None, 0.005, 615, 625)
                mfile.qc_insitu.set_thershold(None, 0.01, 410, 415)

                if fmdb.find('CCI')>0:
                    mfile.qc_insitu.set_thershold(None, 0.003, 650, 670)  ##ONLY FOR CCI
                    mfile.qc_insitu.set_thershold(None, 0.004, 440, 450)  ##ONLY FOR CCI


                # mfile.qc_insitu.set_thershold(None,0.01,400,700)
                # mfile.qc_insitu.set_thershold(None, 0.001, 700, 800)
                #mfile.qc_sat.wl_ref = mfile.wlref
                #mfile.qc_sat.add_theshold_mask(-1,510,0.008,'greater')#STANDARD
                # mfile.qc_sat.add_theshold_mask(-1,412,0.006,'greater')#C2RCC
                # mfile.qc_sat.add_theshold_mask(-1, 620, 0.0055, 'greater')  # C2RCC
                mfile.prepare_df_validation()
                dfadd = mfile.df_validation
                dfadd['satellite'] = self.mdb_list[fmdb]['satellite'].upper()
                dfadd['platform'] = self.mdb_list[fmdb]['platform'].upper()
                dfadd['sensor'] = self.mdb_list[fmdb]['sensor'].upper()
                dfadd['site'] = self.mdb_list[fmdb]['insitu_site_name'].upper()
                ac_here = self.mdb_list[fmdb]['satellite_aco_processor'].upper()
                if ac_here == 'ATMOSPHERIC CORRECTION PROCESSOR: XXX':
                    ac_here = 'STANDARD'
                if ac_here == 'CLIMATE CHANGE INITIATIVE - EUROPEAN SPACE AGENCY':
                    ac_here = 'CCI'
                dfadd['ac'] = ac_here

                if self.df_validation is None:
                    self.df_validation = dfadd
                else:
                    self.df_validation = pd.concat([self.df_validation, dfadd], ignore_index=True)

                mu_dates_here = mfile.mu_dates
                if len(self.mu_dates) == 0:
                    for mu_date in mu_dates_here:
                        list_products_day = list([mu_dates_here[mu_date]])
                        self.mu_dates[mu_date] = list_products_day
                else:
                    for mu_date in mu_dates_here:
                        if mu_date in self.mu_dates.keys():
                            # keys_check = ['satellite', 'platform', 'sensor','site']
                            list_products_day = self.mu_dates[mu_date]
                            list_products_day.append(mu_dates_here[mu_date])
                            self.mu_dates[mu_date] = list_products_day

                        else:
                            list_products_day = list([mu_dates_here[mu_date]])
                            self.mu_dates[mu_date] = list_products_day

    def save_df_validation_to_file(self,path_out):
        if self.df_validation is not None:
            file_data = os.path.join(path_out, 'Data.csv')
            self.df_validation.to_csv(file_data,sep=';')

    def get_df_validation(self, params_include, param_agrup, strict_agrup):
        if self.df_validation is None:
            self.prepare_df_for_validation()

        dfvalid = pd.DataFrame(columns=list(self.df_validation.columns))

        if params_include is None:
            params_include = []

        if param_agrup is not None:
            groups = pd.Series(self.df_validation[param_agrup]).unique()
        else:
            groups = None

        index_valid = 0
        start_date = dt.now()
        end_date = dt(1970, 1, 1)
        for index, row in self.df_validation.iterrows():
            include = True
            for param in params_include:
                pinclude = False
                for val in params_include[param]:
                    if row[param] == val:
                        pinclude = True
                if not pinclude:
                    include = False
            if not include:
                continue
            if param_agrup is not None and strict_agrup:
                date_here = dt.strptime(row['Sat_Time'], '%Y-%m-%d %H:%M')
                if date_here < start_date:
                    start_date = date_here
                if date_here > end_date:
                    end_date = date_here
                sdate = date_here.strftime('%Y-%m-%d')
                list_products = self.mu_dates[sdate]
                valid_mu = True
                for g in groups:
                    if not self.check_validity(list_products, row, param_agrup, g):
                        valid_mu = False
                row['Valid'] = valid_mu

            if row['Valid']:
                #print(type(row), row)
                rowdf = pd.DataFrame.transpose(pd.DataFrame(row))
                if index_valid == 0:
                    dfvalid = rowdf
                else:
                    dfvalid = pd.concat([dfvalid, rowdf], ignore_index=True)
                index_valid = index_valid + 1
                #print(f'Index valid: {index_valid} of {len(dfvalid.index)}')

        print(f'[INFO] # Valid pooints: {index_valid}')
        if param_agrup is not None:
            print(len(dfvalid.index))
            for g in groups:
                dfg = dfvalid[dfvalid[param_agrup] == g]['Index_MU']
                imu = pd.Series(dfg).unique()
                print(f'Group: {g}  # match-ups: {len(imu)}')

        return dfvalid, groups, start_date, end_date



    def check_validity(self, list_products, row, param_agrup, g):
        for product in list_products:
            product_check = True
            for key in product.keys():
                if key != 'mu_valid':
                    value_check = row[key]
                    if key == param_agrup:
                        value_check = g
                    if product[key] != value_check:
                        product_check = False
                        break
            if product_check:
                return product['mu_valid']

        return False

    def make_validation_individual_mdb_files(self):
        for fmdb in self.mdb_list:
            mfile = MDBFile(self.mdb_list[fmdb]['path'])
            if not self.wlref is None:
                mfile.set_wl_ref(self.wlref)
            nvalid = mfile.prepare_df_validation()

            if nvalid == 0:
                print(f'[WARNING] Valid match-ups were not found. Skipping file: {fmdb}')
                return

            # path
            dir_name = os.path.dirname(self.mdb_list[fmdb]['path'])
            path_out = os.path.join(dir_name, fmdb[:-3])
            if not os.path.exists(path_out):
                os.mkdir(path_out)
            mplot = MDBPlot(mfile, None)
            mplot.make_validation_mdbfile(path_out)














