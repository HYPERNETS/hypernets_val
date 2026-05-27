import os
from INSITU_base import INSITUBASE
from netCDF4 import Dataset
import numpy as np

class INSITU_MEDA(INSITUBASE):

    def __init__(self,mdb_options,verbose):
        self.path_aeronet_nc = None
        if mdb_options is not None:
            self.path_aeronet_nc = mdb_options.insitu_path_source
        self.mdb_options = mdb_options
        self.verbose = verbose
        self.VALID = self.check_path_aeronet_nc()
        # else:
        #     self.path_meda = '/store2/data/meda/binary'
        # self.file_name = 'meda_lam_opt_$DATE$_L1v2.nc'
        # self.format_date = '%y%m%d'
        # self.mdb_options = mdb_options
        # self.verbose = verbose
        # self.new_MDB = None
        # self.wavelengths  = None