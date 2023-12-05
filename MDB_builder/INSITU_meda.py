from datetime import datetime as dt
from datetime import timedelta
import os

import pytz

from INSITU_base import INSITUBASE
from netCDF4 import Dataset
import numpy as np

class INSITU_MEDA(INSITUBASE):

    def __init__(self,mdb_options,verbose):
        if mdb_options is not None:
            self.path_meda = mdb_options.insitu_path_source
        else:
            self.path_meda = '/store2/data/meda/binary'
        self.file_name = 'meda_lam_opt_$DATE$_L1v2.nc'
        self.format_date = '%y%m%d'
        self.mdb_options = mdb_options
        self.verbose = verbose
        self.new_MDB = None
        self.wavelengths  = None

    def create_mdb_insitu_extract(self, extract_path, ofile,sat_time):
        file_meda  =self.check_file(sat_time)
        if file_meda is None:
            return False
        dataset = Dataset(file_meda)
        self.wavelengths = np.array(dataset.variables['bands'])
        self.mdb_options.insitu_options['n_insitu_bands'] = len(self.wavelengths)
        dataset.close()
        self.start_add_insitu(extract_path, ofile)

        return True

    def set_data(self,sat_time):

        file_meda = self.check_file(sat_time)
        time_max = self.mdb_options.insitu_options['time_window']  # time max in seconds
        dataset = Dataset(file_meda)
        time_array = np.array(dataset.variables['timetag'])
        rrs_array = np.array(dataset.variables['rrs'])
        if time_array.ndim==0:
            nobs = 1
            time_array = np.array([time_array])
        else:
            nobs = len(time_array)
        dataset.close()

        self.new_MDB.variables['insitu_original_bands'][:] = self.wavelengths[:]
        ihere = 0
        for index in range(nobs):
            ins_time = sat_time.replace(hour=0,minute=0,second=0,microsecond=0)
            ins_time = ins_time + timedelta(hours=float(time_array[index]))


            time_diff = abs((sat_time - ins_time).total_seconds())
            if time_diff < time_max:
                rrs_here = np.transpose(rrs_array[index, :])
                self.new_MDB.variables['insitu_Rrs'][0, :, ihere] = rrs_here
                self.new_MDB.variables['insitu_time'][0, ihere] = np.float64(ins_time.replace(tzinfo=pytz.utc).timestamp())
                ihere = ihere + 1

        if ihere>0:
            return True
        else:
            print(f'[WARNING] In situ data were not available for satellite time {sat_time} within the time window')
            return False

    def close_mdb(self):
        self.new_MDB.close()


    def check_file(self,date):

        path = os.path.join(self.path_meda,date.strftime('%Y'),date.strftime('%m'),date.strftime('%d'))
        if not os.path.isdir(path):
            return None
        date_name_file = date.strftime(self.format_date)
        name_file = self.file_name.replace('$DATE$',date_name_file)
        meda_file = os.path.join(path,name_file)
        if os.path.isfile(meda_file):
            return meda_file
        else:
            return None

    def get_datelist_file(self,output_file,start_date,end_date):
        f1 = open(output_file,'w')
        date_ref = start_date
        while date_ref<=end_date:
            meda_file = self.check_file(date_ref)
            if meda_file is not None:
                date_ref_str = date_ref.strftime('%Y-%m-%d')
                f1.write(date_ref_str)
                f1.write('\n')
            date_ref = date_ref + timedelta(hours=24)
        f1.close()
