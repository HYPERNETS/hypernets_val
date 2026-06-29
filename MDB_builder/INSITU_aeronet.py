import pandas as pd

import __init__,os,sys
from datetime import datetime as dt
from datetime import timezone
import numpy as np
try:
    import INSITU_base as ISb
except:
    import MDB_builder.INSITU_base as ISb



code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
code_aeronet = os.path.join(os.path.dirname(code_home), 'aeronet')
sys.path.append(code_aeronet)


import COMMON.common_functions as cfs


class INSITU_AERONET:

    def __init__(self,insitu_options,verbose):
        self.insitu_type = 'AERONET_OC'
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

        ##for MDBBuilder
        self.wl_arrays = None

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

    def get_metadata_date(self,df,datehere):
        from base.anet_nc_reader import AERONETReader
        areader = AERONETReader(self.path_aeronet_nc)
        row_ini, row_fin = areader.get_indices_date(datehere.strftime('%Y-%m-%d'))
        areader.row_ini, areader.row_fin = row_ini,row_fin
        time_list = areader.extract_time_list()
        time_list_str = [x.strftime('%Y-%m-%dT%H:%M:%S') for x in time_list]
        indices = np.arange(row_ini,row_fin+1).astype(np.int32)
        extract_file = df['extract_file'][0]
        ntimes = len(time_list_str)
        data = {
            'extract_file': [extract_file]*ntimes,
            'insitu_index': indices,
            'insitu_time': time_list_str,
            'insitu_lat': [self.lat]*ntimes,
            'insitu_lon': [self.lon]*ntimes
        }
        df_new = pd.DataFrame(data=data)
        return df_new

    def get_insitu_time_list(self,row_ini,row_fin):
        from base.anet_nc_reader import AERONETReader
        areader = AERONETReader(self.path_aeronet_nc)
        areader.row_ini, areader.row_fin = row_ini,row_fin
        time_list = areader.extract_time_list()
        time_list_str = [x.strftime('%Y-%m-%dT%H:%M:%S') for x in time_list]
        areader.dataset.close()
        return time_list_str

    def get_insitu_instrument_id_array(self,indices,ninsitu_max):
        from base.anet_nc_reader import AERONETReader
        areader = AERONETReader(self.path_aeronet_nc)
        ewl = areader.dataset.variables['Exact_Wavelengths'][indices,:]
        areader.dataset.close()

        instrument_ids = np.zeros((1, ninsitu_max)).astype(np.int16)
        instrument_ids = np.ma.masked_equal(instrument_ids,0)
        for idx in range(ewl.shape[0]):
            diff = np.abs(np.sum(ewl[idx,:] - self.wl_arrays,axis=1))
            instrument_ids[0,idx] = int(np.argmin(diff))+1

        return instrument_ids

    def get_rrs_array(self,indices,ninsitu_max):
        from base.anet_nc_reader import AERONETReader
        areader = AERONETReader(self.path_aeronet_nc)
        areader.row_ini, areader.row_fin = indices[0],indices[-1]
        rrs_array = areader.extract_rrs(False)
        areader.dataset.close()
        rrs_array = np.transpose(rrs_array)
        return rrs_array



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
            date_here = dt.strptime(date_ts,'%Y-%m-%d')
            if self.start_date is not None and date_here<self.start_date:
                continue
            if self.end_date is not None and date_here>self.end_date:
                continue
            self.file_list[date_ts] = self.path_aeronet_nc
            self.date_list.append(date_here)

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

    def check_rrs_and_data_variables(self):
        if not self.VALID:
            return False
        try:
            from base.anet_nc_reader import AERONETReader
            areader = AERONETReader(self.path_aeronet_nc)
            ewl = areader.dataset.variables['Exact_Wavelengths'][:]
            wl_arrays = np.unique(ewl,axis=0)
            indices = [np.where(ewl==wl_arrays[idx])[0][0] for idx in range(wl_arrays.shape[0])]
            self.wl_arrays = wl_arrays[np.argsort(indices),:]
            areader.dataset.close()
            print(f'[INFO] Retrieving of wavelength arrays completed. Identified {self.wl_arrays.shape[0]} instruments.')
        except Exception as ex:
            print(f'[ERROR] Error retrieving wl_arrays from ExactWavelengths variable. Exception: {ex}')
        return True

    def create_mini_mdb_files(self, options, extract_dir, extracts, time_extracts, overwrite):
        from base.anet_nc_reader import AERONETReader
        areader = AERONETReader(self.path_aeronet_nc)
        options['n_insitu_bands'] = areader.nwl
        fcsv = os.path.join(extract_dir, f'MDBm_{time_extracts[0].strftime("%Y%m%d")}.csv')
        fw = open(fcsv, 'w')
        started = False
        dims = None
        for t in time_extracts:
            ref = t.strftime('%Y%m%dT%H%M%S')
            # print(ref,'-->',extracts[ref]['file'])
            file_out = ISb.get_mini_mdb_file_path(extract_dir, extracts[ref])
            write_line = False
            if os.path.isfile(file_out) and not overwrite:
                print(f'[INFO] Mini MDB file {os.path.basename(file_out)} already exists. Skipping...')
                write_line = True
            else:
                if self.verbose:
                    print(f'[INFO] Creating mini MDB file {os.path.basename(file_out)}')
                self.create_mini_mdb_file_impl(file_out, options, extracts[ref])
                if os.path.isfile(file_out):
                    write_line = True
                else:
                    print(f'[WARNING] File MDBm {os.path.basename(file_out)} could not be created. Skipping...')
            if write_line:
                fw, dims = ISb.add_line_csv_with_mdbm_info(fw, file_out, started)
                started = True

        fw.close()

        return dims

    def create_mini_mdb_file_impl(self,file_out, options,extract_info):
        # print(file_out)
        # print(options)
        # print(extract_info)
        ninsitu_real = self.check_ninsitu_real(extract_info)
        if ninsitu_real < 0:
            return
        if ninsitu_real == 0:
            print(f'[WARNING] No in situ data found for {extract_info["file"]}')
            return

        if ninsitu_real > options['ninsitu_max']:
            print(
                f'[WARNING] {ninsitu_real} is greater than the maximum number of in situ data points {options["ninsitu_max"]}. ')
            return
        if self.verbose:
            print(f'[INFO] Number of in situ data points for the extract: {ninsitu_real}')

        if self.wl_arrays is None:
            print(f'[ERROR] No wl_arrays found for {extract_info["file"]}]')
            return
        nids = self.wl_arrays.shape[0]
        options['instrument_ids'] = [f'AERONET_OC_{x+1}' for x in range(nids)]

        builder = ISb.Mini_MDB_Builder(options, self.verbose)
        builder.start_mini_mdb(extract_info['file'], file_out)
        ##in situ original bands
        if not builder.set_insitu_wavelengths_all_instruments(self.wl_arrays):
            builder.close_mini_mdb_file()
            os.remove(extract_info['file'])
            return
        ##instrument ids
        instrument_ids = self.get_insitu_instrument_id_array(extract_info['insitu_indices'],options['ninsitu_max'])
        if not builder.set_instrument_id(options['ninsitu_max'],instrument_ids):
            builder.close_mini_mdb_file()
            os.remove(extract_info['file'])
            return

        ##basic variables
        if not builder.set_insitu_basic_variables_from_dict(extract_info,is_fixed_site=True):
            builder.close_mini_mdb_file()
            os.remove(extract_info['file'])
            return

        ##Rrs
        rrs_array = self.get_rrs_array(extract_info['insitu_indices'],options['ninsitu_max'])
        if not builder.set_insitu_rrs(rrs_array):
            builder.close_mini_mdb_file()
            os.remove(extract_info['file'])
            return

        builder.close_mini_mdb_file()


    def check_ninsitu_real(self,extract_info):
        if not 'insitu_indices' in extract_info or extract_info['insitu_indices'] is None:
            print(f'[ERROR] insitu_indices is required in the extract info')
            return -1
        ninsitu_real = len(extract_info['insitu_indices'])
        if not 'insitu_lat' in extract_info or extract_info['insitu_lat'] is None:
            print(f'[ERROR] insitu_lat is required in the extract info')
            return -1
        if not 'insitu_lon' in extract_info or extract_info['insitu_lon'] is None:
            print(f'[ERROR] insitu_lon is required in the extract info')
            return -1
        if not 'insitu_time' in extract_info or extract_info['insitu_time'] is None:
            print(f'[ERROR] insitu_time is required in the extract info')
            return -1
        if not 'time_diff' in extract_info or extract_info['time_diff'] is None:
            print(f'[ERROR] time_diff is required in the extract info')
            return -1

        if len(extract_info['insitu_lat'])<ninsitu_real:
            print(f'[ERROR] Discrepancy in the number of in situ data points between insitu_lat {len(extract_info["insitu_lat"])} and insitu_indices {ninsitu_real}')
            return -1
        if len(extract_info['insitu_lon'])<ninsitu_real:
            print(f'[ERROR] Discrepancy in the number of in situ data points between insitu_lat {len(extract_info["insitu_lon"])} and insitu_indices {ninsitu_real}')
            return -1
        if len(extract_info['insitu_time'])<ninsitu_real:
            print(f'[ERROR] Discrepancy in the number of in situ data points between insitu_lat {len(extract_info["insitu_time"])} and insitu_indices {ninsitu_real}')
            return -1
        if len(extract_info['time_diff'])<ninsitu_real:
            print(f'[ERROR] Discrepancy in the number of in situ data points between time_diff {len(extract_info["time_diff"])} and insitu_indices {ninsitu_real}')
            return -1

        if 'insitu_spatial_index' in extract_info and extract_info['insitu_spatial_index'] is not None and len(extract_info['insitu_spatial_index'])<ninsitu_real:
            print(f'[ERROR] Discrepancy in the number of in situ data points between insitu_spatial_index {len(extract_info["insitu_spatial_index"])}and insitu_indices {ninsitu_real}')
            return -1



        return ninsitu_real