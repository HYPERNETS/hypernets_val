import os.path

from INSITU_base import INSITUBASE
from netCDF4 import Dataset
import numpy as np
from datetime import datetime as dt

class INSITU_WISP3(INSITUBASE):

    def __init__(self, mdb_options, insitu_file, verbose):
        self.mdb_options = mdb_options
        self.insitu_file = insitu_file
        self.verbose = verbose
        self.new_MDB = None
        self.date_indices = {}
        self.wavelenghs = None
        self.start_date_indices_and_wavelengths()


    def start_date_indices_and_wavelengths(self):
        dataset = Dataset(self.insitu_file)
        time_array = np.array(dataset.variables['Time'])
        self.wavelenghs = np.array(dataset.variables['Nominal_Wavelenghts'])
        self.mdb_options.insitu_options['n_insitu_bands'] = len(self.wavelenghs)
        for itime in range(len(time_array)):
            time = time_array[itime]
            datestr = dt.utcfromtimestamp(float(time)).strftime('%Y-%m-%d')
            if datestr not in self.date_indices.keys():
                self.date_indices[datestr] = {
                    'index_min': itime,
                    'index_max': itime
                }
            else:
                if itime<self.date_indices[datestr]['index_min']:
                    self.date_indices[datestr]['index_min'] = itime
                if itime>self.date_indices[datestr]['index_max']:
                    self.date_indices[datestr]['index_max'] = itime
        dataset.close()

    def create_mdb_insitu_extract(self, extract_path, ofile):
        self.start_add_insitu(extract_path, ofile)

    def add_new_variables(self):
        for var_name in self.insitu_spectral_variables:
            type = self.insitu_spectral_variables[var_name]['type']
            var = self.new_MDB.createVariable(var_name, type, ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                              zlib=True, complevel=6)
            for at in self.insitu_spectral_variables[var_name]:
                if at == 'type' or at == 'name_orig':
                    continue
                var.setncattr(at, self.insitu_spectral_variables[var_name][at])
            var[:] = -999
        # self.new_MDB.close()
        if self.verbose:
            print('[INFO] Added new variables')


    def set_data(self,sat_time):
        time_max = self.mdb_options.insitu_options['time_window']  # time max in seconds
        #time_maxh = time_max / 3600

        date_ref_str = sat_time.strftime('%Y-%m-%d')
        if not date_ref_str in self.date_indices.keys():
            print(f'[WARNING] In situ data were not available for satellite date: {date_ref_str}')
            return False

        index_min = self.date_indices[date_ref_str]['index_min']
        index_max = self.date_indices[date_ref_str]['index_max']
        dataset = Dataset(self.insitu_file)
        time_array = np.array(dataset.variables['Time'])
        rrs_array = np.array(dataset.variables['RRS'])
        self.new_MDB.variables['insitu_original_bands'][:] = self.wavelenghs[:]
        ihere = 0
        for index in range(index_min,index_max+1):
            ins_time = dt.utcfromtimestamp(float(time_array[index]))
            time_diff = abs((sat_time-ins_time).total_seconds())
            if time_diff<time_max:
                rrs_here = np.transpose(rrs_array[index,:])
                #print(sat_time,ins_time,ihere,rrs_here.shape,time_diff,time_max)
                self.new_MDB.variables['insitu_Rrs'][0,:,ihere] = rrs_here
                self.new_MDB.variables['insitu_time'][0,ihere]= time_array[index]
                ihere = ihere + 1

        dataset.close()

        if ihere>0:
            return True
        else:
            print(f'[WARNING] In situ data were not available for satellite time {sat_time} within the time window')
            return False

    def close_mdb(self):
        self.new_MDB.close()