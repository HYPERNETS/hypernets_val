import datetime
import os.path

from MDBFile import MDBFile
from MDBPlot import MDBPlot
import pandas as pd
from datetime import datetime as dt



class MDBFileList():
    def __init__(self):
        self.mdb_list = {}
        self.df_validation = None  # pd.DataFrame(columns=col_names)
        self.mu_dates = {}
        # self.date_list = {}

        self.wlref = None

    def set_wl_ref(self, wllist):
        self.wlref = wllist

    def add_mdb_file(self, path_mdb):
        mfile = MDBFile(path_mdb)
        if not mfile.VALID:
            return False
        name = path_mdb.split('/')[-1]
        self.mdb_list[name] = {
            'path': path_mdb,
            'include': True
        }
        for key in mfile.info:
            self.mdb_list[name][key] = mfile.info[key]
        if self.mdb_list[name]['satellite_aco_processor'] == 'Atmospheric Correction processor: xxx':
            self.mdb_list[name]['satellite_aco_processor'] = 'STANDARD'

    def prepare_df_for_validation(self):
        for fmdb in self.mdb_list:
            if self.mdb_list[fmdb]['include']:
                mfile = MDBFile(self.mdb_list[fmdb]['path'])
                if not self.wlref is None:
                    mfile.set_wl_ref(self.wlref)
                mfile.prepare_df_validation()
                dfadd = mfile.df_validation
                dfadd['satellite'] = self.mdb_list[fmdb]['satellite'].upper()
                dfadd['platform'] = self.mdb_list[fmdb]['platform'].upper()
                dfadd['sensor'] = self.mdb_list[fmdb]['sensor'].upper()
                dfadd['site'] = self.mdb_list[fmdb]['insitu_site_name'].upper()
                ac_here = self.mdb_list[fmdb]['satellite_aco_processor'].upper()
                if ac_here == 'ATMOSPHERIC CORRECTION PROCESSOR: XXX':
                    ac_here = 'STANDARD'
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

    def get_df_validation(self, params_include, param_agrup, strict_agrup):
        if self.df_validation is None:
            self.prepare_df_for_validation()

        dfvalid = pd.DataFrame(columns=list(self.df_validation.columns))

        if param_agrup is not None:
            groups = pd.Series(self.df_validation[param_agrup]).unique()

        index_valid = 0
        start_date = dt.now()
        end_date = dt(1970,1,1)
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
                if date_here<start_date:
                    start_date = date_here
                if date_here>end_date:
                    end_date = date_here
                sdate = date_here.strftime('%Y-%m-%d')
                list_products = self.mu_dates[sdate]
                valid_mu = True
                for g in groups:
                    if not self.check_validity(list_products, row, param_agrup, g):
                        valid_mu = False
                row['valid'] = valid_mu

            if row['valid']:
                print(type(row), row)
                rowdf = pd.DataFrame.transpose(pd.DataFrame(row))
                if index_valid == 0:
                    dfvalid = rowdf
                else:
                    dfvalid = pd.concat([dfvalid, rowdf], ignore_index=True)
                index_valid = index_valid + 1
                print(f'Index valid: {index_valid} But {len(dfvalid.index)}')

        print('Nvalid: ', index_valid)
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
            mfile.prepare_df_validation()
            # path
            dir_name = os.path.dirname(self.mdb_list[fmdb]['path'])
            path_out = os.path.join(dir_name, fmdb[:-3])
            if not os.path.exists(path_out):
                os.mkdir(path_out)
            mplot = MDBPlot(mfile, None)
            mplot.make_validation_mdbfile(path_out)
