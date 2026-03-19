import os
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from datetime import datetime as dt
try:
    import INSITU_base as ISb
except:
    import MDB_builder.INSITU_base as ISb

class INSITU_SINGLE_CSV:

    def __init__(self, insitu_options, verbose):
        self.fixed_site = False
        self.insitu_type = 'SINGLE_CSV'
        self.verbose = verbose
        self.insitu_options = insitu_options
        self.site = self.insitu_options['site']
        self.file_list = {}
        self.date_list = []
        self.start_date = None
        self.end_date = None
        self.date_list_check = None
        #print(insitu_options)
        self.non_spectral_vars = None

    def get_csv_options(self):
        col_date = self.insitu_options['col_date']
        col_time = self.insitu_options['col_time']
        col_lat = self.insitu_options['col_lat']
        col_lon = self.insitu_options['col_lon']
        format_date = self.insitu_options['format_date']
        format_time = self.insitu_options['format_time']
        col_sep = self.insitu_options['col_sep']
        return col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep

    def check_csv_structure(self,file_csv):
        col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = self.get_csv_options()
        try:
            df = pd.read_csv(file_csv, sep=col_sep)
        except Exception as ex:
            print(f'[ERROR] File {file_csv} is not a valid csv separated by {col_sep}. Exception: {ex}')
            return False
        valid = True
        if not col_date in df.columns:
            print(f'[ERROR] col_date {col_date} is not available in the CSV file: {file_csv}')
            valid = False
        if col_time is not None and not col_time in df.columns:
            print(f'[ERROR] col_time {col_date} is not available in the CSV file: {file_csv}')
            valid = False
        if not col_lat in df.columns:
            print(f'[ERROR] col_lat {col_lat} is not available in the CSV file: {file_csv}')
            valid = False
        if not col_lon in df.columns:
            print(f'[ERROR] col_lon {col_lon} is not available in the CSV file: {file_csv}')
            valid = False
        return valid

    def check_data(self):
        file_csv = self.insitu_options['path_csv']
        if file_csv is None:
            print(f'[ERROR] path_csv option could not be None')
            return False
        if not os.path.isfile(file_csv):
            print(f'[ERROR] path_csv {file_csv} is not available or is not a valid file')
            return False
        if not self.check_csv_structure(file_csv):
            return False
        return True

    def get_date_list_from_dataframe(self,df, col_date, format_date, col_time, format_time):
        date_array_ts = df[col_date]
        time_array_ts = None
        if col_time is not None:
            format_datetime = f'{format_date}T{format_time}'
            time_array_ts = df[col_time]
        else:
            format_datetime = format_date
        only_date_array = []
        time_list = []
        for idx, xd in enumerate(date_array_ts):
            x = f'{xd}T{time_array_ts[idx]}' if time_array_ts is not None else xd
            try:
                time_here = dt.strptime(x, format_datetime)
                time_list.append(time_here)
                only_date_array.append(time_here.strftime('%Y-%m-%d'))
            except Exception as ex:
                print(f'[WARNING] {x} could not be parsed with format {format_datetime}. Exception: {ex}')
                continue

        return only_date_array, time_list

    def prepare_data(self):
        if not self.check_data():
            return False

        file_csv = self.insitu_options['path_csv']
        col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = self.get_csv_options()
        df = pd.read_csv(file_csv, sep=col_sep)
        only_date_array, time_list = self.get_date_list_from_dataframe(df, col_date, format_date, col_time, format_time)
        if len(only_date_array)==0:
            print(f'[ERROR] No valid dates retrieved from the CSV file.')
            return False
        only_date_array_unique = np.unique(only_date_array).tolist()
        if len(only_date_array_unique) == 0:
            print(f'[ERROR] No valid dates retrieved from the CSV file.')
            return False
        for date_ts in only_date_array_unique:
            self.file_list[date_ts] = file_csv
            self.date_list.append(dt.strptime(date_ts,'%Y-%m-%d'))
        if len(self.date_list)==0:
            print(f'[ERROR] No valid dates were retrieved from CSV files')
            return False
        self.date_list.sort()
        return True

    def get_metadata_date(self,datehere):
        date_ts = datehere.strftime('%Y-%m-%d')
        file_csv = self.file_list[date_ts]
        col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = self.get_csv_options()
        df = pd.read_csv(file_csv,sep=col_sep)
        only_date_array, insitu_time = self.get_date_list_from_dataframe(df, col_date, format_date, col_time, format_time)
        insitu_lat = np.ma.array(df[col_lat])
        insitu_lon = np.ma.array(df[col_lon])
        insitu_time = np.array(insitu_time)

        insitu_indices = np.where(np.array(only_date_array) == date_ts)

        if len(insitu_indices[0])>1:
            insitu_time_here = insitu_time[insitu_indices]
            sorted_indices_here = np.argsort(insitu_time_here)
            sorted_indices = insitu_indices[0][sorted_indices_here]
            insitu_indices = (sorted_indices,)

        insitu_time = insitu_time[insitu_indices]
        insitu_lat = insitu_lat[insitu_indices]
        insitu_lon = insitu_lon[insitu_indices]



        return insitu_time,insitu_lat,insitu_lon,insitu_indices

    def get_ref_date(self, datehere):
        path_s = self.insitu_options['path_csv']
        name = os.path.basename(path_s)
        name_ref = f'{name[0:name.rfind(".")]}_{datehere.strftime("%Y%m%d")}'
        return name_ref

    def check_rrs_and_data_variables(self):
        if not self.check_data():
            return False
        file_csv = self.insitu_options['path_csv']
        col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = self.get_csv_options()
        df = pd.read_csv(file_csv, sep=col_sep)
        col_names = df.columns.to_list()
        if self.insitu_options['non_spectral_variables'] is not None:
            self.non_spectral_vars = []
            for var_name in self.insitu_options['non_spectral_variables']:
                if var_name in col_names:
                    self.non_spectral_vars.append(var_name)
                else:
                    print(f'[ERROR] {var_name} is not included in the variable list')
                    return False
            if len(self.non_spectral_vars)==0:
                self.non_spectral_vars = None
        return True

    def create_mini_mdb_files(self,options,extract_dir,extracts,time_extracts,overwrite):
        options['n_insitu_bands'] = 0
        fcsv = os.path.join(extract_dir,f'MDBm_{time_extracts[0].strftime("%Y%m%d")}.csv')
        fw  = open(fcsv,'w')
        started = False
        dims = None
        for t in time_extracts:
            ref = t.strftime('%Y%m%dT%H%M%S')
            #print(ref,'-->',extracts[ref]['file'])
            file_out = ISb.get_mini_mdb_file_path(extract_dir,extracts[ref])
            write_line = False
            if os.path.isfile(file_out) and  not overwrite:
                print(f'[INFO] Mini MDB file {os.path.basename(file_out)} already exists. Skipping...')
                write_line = True
            else:
                if self.verbose:
                    print(f'[INFO] Creating mini MDB file {os.path.basename(file_out)}')
                self.create_mini_mdb_file_impl(file_out,options,extracts[ref])
                if os.path.isfile(file_out):
                   write_line = True
                else:
                    print(f'[WARNING] File MDBm {os.path.basename(file_out)} could not be created. Skipping...')
            if write_line:
                fw,dims = self.add_line_csv_with_MDBm_info(fw,file_out,started)
                started = True

        fw.close()

        return dims

    def add_line_csv_with_MDBm_info(self,fw,file_nc,started):
        if not started:
            first_line = 'name;satellite_id;insitu_id;instrument_id;satellite_bands;insitu_bands;rows;columns'
            fw.write(first_line)
        dset  = Dataset(file_nc)
        nsat = len(dset.dimensions['satellite_id'])
        ninsitu = len(dset.dimensions['insitu_id'])
        ninstrument = len(dset.dimensions['instrument_id'])
        nwlsat = len(dset.dimensions['satellite_bands']) if 'satellite_bands' in dset.dimensions else 0
        nwlinsitu = len(dset.dimensions['insitu_bands']) if 'insitu_bands' in dset.dimensions else 0
        rows = len(dset.dimensions['rows'])
        cols = len(dset.dimensions['columns'])
        line = f'{os.path.basename(file_nc)};{nsat};{ninsitu};{ninstrument};{nwlsat};{nwlinsitu};{rows};{cols}'
        dset.close()
        fw.write('\n')
        fw.write(line)
        dims = np.array([nsat,ninsitu,ninstrument,nwlsat,nwlinsitu,rows,cols])

        return fw,dims

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

    def create_mini_mdb_file_impl(self,file_out,options,extract_info):
        ninsitu_real = self.check_ninsitu_real(extract_info)
        if ninsitu_real<0:
            return
        if ninsitu_real==0:
            print(f'[WARNING] No in situ data found for {extract_info["file"]}')
            return

        if ninsitu_real>options['ninsitu_max']:
            print(f'[WARNING] {ninsitu_real} is greater than the maximum number of in situ data points {options["ninsitu_max"]}. ')
            return
        if self.verbose:
            print(f'[INFO] Number of in situ data points for the extract: {ninsitu_real}')

        file_csv = self.insitu_options['path_csv']
        col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = self.get_csv_options()
        df = pd.read_csv(file_csv, sep=col_sep)


        #insitu_indices = extract_info['insitu_indices']
        builder = ISb.Mini_MDB_Builder(options,self.verbose)
        builder.start_mini_mdb(extract_info['file'],file_out)
        builder.add_shipborne_variables()
        if options['n_insitu_bands']>0  and self.insitu_options['rrs_unc_format'] is not None:
            builder.add_rrs_uncentainty_variable()

        building_error = False
        if len(options['instrument_ids'])==1:
            if options['n_insitu_bands'] > 0:
                builder.set_insitu_wavelengths(0,self.wl_array)
            builder.set_instrument_id(ninsitu_real,1)
        #basic variables: insitu_lat,insitu_lon,insitu_time,time_diff,insitu_spatial_index

        if not builder.set_insitu_basic_variables_from_dict(extract_info):
            print(f'[ERROR] Error setting the basic variables. MDBm file could not be created')
            building_error = True

        indices = extract_info['insitu_indices']

        ##RRS and RRS_UNC NOT IMPLEMENTED FOR INSITU_singlecsv
        if options['n_insitu_bands'] > 0:
            if self.get_rrs_array():
                rrs_here = np.ma.transpose(self.rrs_array[indices,:])
                if not builder.set_spectral_variables('insitu_Rrs',rrs_here):
                    print(f'[ERROR] Error setting the in situ Rrs variable. MDBm file could not be created')
                    building_error = True
            else:
                print(f'[ERROR] Error retrieving the Rrs variable')
                building_error = True

            if self.insitu_options['rrs_unc_format'] is not None:
                if self.get_rrs_unc_array():
                    rrs_unc_here = np.ma.transpose(self.rrs_unc_array[indices, :])
                    if not builder.set_spectral_variables('insitu_Rrs_unc', rrs_unc_here):
                        print(f'[ERROR] Error setting the in situ Rrs unc variable. MDBm file could not be created')
                        building_error = True
                else:
                    print(f'[ERROR] Error retrieving the Rrs_unc variable from the SeaBass file')
                    building_error = True

        if self.non_spectral_vars is not None:
            for var_name in self.non_spectral_vars:
                var_name_nc = f'insitu_{var_name}' if not var_name.startswith('insitu_') else var_name
                #attrs = {'units': self.non_spectral_vars[var_name]}
                attrs = None##non implemented a method
                builder.add_non_spectral_variable(var_name_nc, attrs)
                array = np.array(df[var_name][:])
                array = array[indices]
                # array = np.array(self.sb.data[var_name])
                # array = array[indices]
                builder.set_non_spectral_variables(var_name_nc,array)

        builder.close_mini_mdb_file()
        if building_error:
            os.remove(file_out)

    def prepare_csv_metadata(self, output_file):

        if not self.check_data():
            return False

        file_csv = self.insitu_options['path_csv']
        col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = self.get_csv_options()
        df = pd.read_csv(file_csv, sep=col_sep)
        only_date_array, time_list = self.get_date_list_from_dataframe(df, col_date, format_date, col_time, format_time)
        lat_array = np.ma.array(df[col_lat][:])
        lon_array = np.ma.array(df[col_lon][:])
        fw = open(output_file, 'w')
        started = False
        only_date_array = np.array(only_date_array)
        date_array_unique = np.unique(only_date_array)


        for date_ts in date_array_unique:
            date_here = dt.strptime(date_ts,'%Y-%m-%d')
            if self.start_date is not None and self.end_date is not None:
                if date_here<self.start_date or date_here>self.end_date:
                    if self.verbose:
                        print(f'[INFO][singlecsv] Skipping date {date_here.strftime("%Y-%m-%d")} as it not in the range: {self.start_date.strftime("%Y-%m-%d")} to {self.end_date.strftime("%Y-%m-%d")}')
                    continue
            if self.verbose:
                print(f'[INFO][singlecsv] Getting metadata for date: {date_ts}')
            indices = np.where(only_date_array == date_ts)

            insitu_lat = lat_array[indices]
            insitu_lon = lon_array[indices]
            lat_min = np.ma.min(insitu_lat) - 0.001
            lat_max = np.ma.max(insitu_lat) + 0.001
            lon_min = np.ma.min(insitu_lon) - 0.001
            lon_max = np.ma.max(insitu_lon) + 0.001
            line = f'{date_ts},{lat_min},{lat_max},{lon_min},{lon_max}'
            if started:
                fw.write('\n')
            if not started:
                started = True
            fw.write(line)

        fw.close()
        if not started:
            print(f'[WARNING][singlecsv] No data were found for the given data range in the singlecsv dataset.')
            os.remove(output_file)
            return False
        else:
            return True