import __init__,os,sys
from datetime import datetime as dt
from datetime import timezone
import numpy as np
try:
    from INSITU_base import INSITUBASE
except:
    from MDB_builder.INSITU_base import INSITUBASE



code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
code_aeronet = os.path.join(os.path.dirname(code_home), 'aeronet')
sys.path.append(code_aeronet)

import COMMON.common_functions as cfs


class INSITU_AERONET:

    def __init__(self,insitu_options,verbose):

        self.path_aeronet_nc = insitu_options['path_aeronet']
        self.verbose = verbose
        #self.insitu_options = insitu_options
        self.site = insitu_options['site']
        self.file_list = {}
        self.date_list = []
        self.start_date = None
        self.end_date = None
        self.lat, self.lon = cfs.get_lat_lon_ins(self.site)
        self.fixed_site = True
        self.VALID = self.check_path_aeronet_nc()


    def check_path_aeronet_nc(self):
        try:
            from base.anet_nc_reader import AERONETReader
        except Exception as ex:
            print(f'[ERROR] aeronet code is not available. Exception: {ex}')
            if not os.path.isdir(code_aeronet):
                print(f'[ERROR] Folder with the code {code_aeronet} does not exist. You can download it by using GIT:')
                print(f'git clone https://github.com/luiscnr/aeronet.git')
            return False
        try:
            areader = AERONETReader(self.path_aeronet_nc)
            areader.dataset.close()
            return True
        except Exception as ex:
            print(f'[ERROR] AERONET-OC file {self.path_aeronet_nc} could not be started. Exception: {ex}')
            return False

    def get_ref_date(self,datehere):
        ref_date = f'AERONET_OC_DAY_{datehere.strftime("%Y%m%d")}'
        return ref_date

    # basic metadata for extracts for fixed sites
    def get_metadata_date_basic(self, datehere):
        insitu_time = np.array([datehere])
        insitu_lat = np.array([self.lat])
        insitu_lon = np.array([self.lon])
        insitu_indices = (np.array([0]),)
        return insitu_time, insitu_lat, insitu_lon, insitu_indices

    def prepare_data(self):
        print(f'[INFO] Launching prepare_data() method in INSITU_AERONET class (INSITU_aeronet.py)')
        if not self.VALID:
            return False
        from base.anet_nc_reader import AERONETReader
        areader = AERONETReader(self.path_aeronet_nc)
        time_list = areader.extract_time_list()
        areader.dataset.close()
        time_list = [t.replace(tzinfo=timezone.utc) for t in time_list]
        only_date_array = [t.strftime('%Y-%m-%d') for t in time_list]
        if len(only_date_array)==0:
            print(f'[ERROR] No valid dates retrieved from the AERONET-OC file.')
            return False
        only_date_array_unique = np.unique(only_date_array).tolist()
        for date_ts in only_date_array_unique:
            self.file_list[date_ts] = self.path_aeronet_nc
            self.date_list.append(dt.strptime(date_ts,'%Y-%m-%d'))

        self.date_list.sort()
        return True

    def prepare_csv_metadata(self,output_file):
        print(f'[INFO] Launching prepare_csv_metadata() method in INSITU_AERONET class (INSITU_aeronet.py)')
        if not self.prepare_data():
            return False

        lat_min = np.ma.min(self.lat) - 0.001
        lat_max = np.ma.max(self.lat) + 0.001
        lon_min = np.ma.min(self.lon) - 0.001
        lon_max = np.ma.max(self.lon) + 0.001

        fw = open(output_file, 'w')
        started = False

        for date_here in self.date_list:
            if self.start_date is not None and date_here<self.start_date:
                continue
            if self.end_date is not None and date_here>self.end_date:
                continue
            date_ts = date_here.strftime('%Y-%m-%d')
            if self.verbose:
                print(f'[INFO] Getting metadata for date: {date_ts}')
            line = f'{date_ts},{lat_min},{lat_max},{lon_min},{lon_max}'
            if started:
                fw.write('\n')
            if not started:
                started = True
            fw.write(line)

        fw.close()


        return True
