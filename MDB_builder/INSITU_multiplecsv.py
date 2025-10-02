import os
import shutil

import pandas as pd
from datetime import datetime as dt
import numpy as np

class INSITU_MULTIPLE_CSV():

    def __init__(self, insitu_options, verbose):
        self.fixed_site = False
        self.insitu_type = 'MULTIPLE_CSV'
        self.verbose = verbose
        self.insitu_options = insitu_options
        self.site = self.insitu_options['site']
        self.file_list = {}
        self.date_list = []
        self.start_date = None
        self.end_date = None

    def check_data(self):
        if not 'path_csv' in self.insitu_options:
            print(f'[ERROR] path_csv is not available in the options dictionary')
            return False
        path_csv = self.insitu_options['path_csv']
        if not os.path.isdir(path_csv):
            print(f'[ERROR] {path_csv} does not exist or is not a valid directory.')
            return False
        ##checking CSV files
        for name in os.listdir(path_csv):
            if not name.endswith('csv'):
                continue
            file_csv = os.path.join(path_csv, name)
            if self.verbose:
                print(f'[INFO] Checking in situ file: {file_csv}')
            if not self.check_csv_structure(file_csv):
                return False
        return True

    def check_csv_structure(self,file_csv):
        col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = self.get_csv_options()
        try:
            df = pd.read_csv(file_csv, sep=col_sep)
        except:
            print(f'[ERROR] File {file_csv} is not a valid csv separated by {col_sep}')
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

    def get_csv_options(self):
        col_date = self.insitu_options['col_date']
        col_time = self.insitu_options['col_time']
        col_lat = self.insitu_options['col_lat']
        col_lon = self.insitu_options['col_lon']
        format_date = self.insitu_options['format_date']
        format_time = self.insitu_options['format_time']
        col_sep = self.insitu_options['col_sep']
        return col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep

    #run only in check_data is true
    def prepare_data(self):
        if not self.check_data():
            return False
        path_csv = self.insitu_options['path_csv']
        col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = self.get_csv_options()

        for name in os.listdir(path_csv):
            if not name.endswith('csv'):
                continue
            file_csv = os.path.join(path_csv, name)
            df = pd.read_csv(file_csv, sep=col_sep)
            only_date_array, time_list = self.get_date_list_from_dataframe(df, col_date, format_date, col_time, format_time)
            only_date_array_unique = np.unique(only_date_array).tolist()
            if len(only_date_array_unique) == 0:
                print(f'[ERROR] No valid dates retrieved from the CSV file.')
                return False
            if len(only_date_array_unique) > 1:
                print(f'[ERROR] Each CSV file should also contain data for a single day. Check CSV parameters')
                return False
            date_ts = only_date_array_unique[0]
            self.file_list[date_ts] = file_csv
            self.date_list.append(dt.strptime(date_ts,'%Y-%m-%d'))

        if len(self.date_list)==0:
            print(f'[ERROR] No valid dates were retrieved from CSV files')
            return False
        self.date_list.sort()

        return True

    def get_ref_date(self, datehere):
        date_ts = datehere.strftime('%Y-%m-%d')
        file_csv = self.file_list[date_ts]
        name = os.path.basename(file_csv)
        return name[0:name.rfind('.')]




    def get_metadata_date(self,datehere):
        date_ts = datehere.strftime('%Y-%m-%d')
        file_csv = self.file_list[date_ts]
        col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = self.get_csv_options()
        df = pd.read_csv(file_csv,sep=col_sep)
        only_date_array, insitu_time = self.get_date_list_from_dataframe(df, col_date, format_date, col_time, format_time)
        insitu_lat = np.ma.array(df[col_lat])
        insitu_lon = np.ma.array(df[col_lon])
        #insitu_indices = np.where(np.array(only_date_array) == date_ts)
        ##sorting
        insitu_time = np.array(insitu_time)
        sorted_indices = np.argsort(insitu_time)
        insitu_time = insitu_time[sorted_indices]
        insitu_lat = insitu_lat[sorted_indices]
        insitu_lon = insitu_lon[sorted_indices]
        insitu_indices = (sorted_indices,)

        return insitu_time,insitu_lat,insitu_lon,insitu_indices

    def create_csv_metadata_for_date(self,insitu_date, file_csv_out):
        if self.file_list is None:
            self.prepare_data()
        if self.file_list is None:
            return False
        date_ts = insitu_date.strftime('%Y-%m-%d')
        if date_ts in self.file_list:
            file_csv_in = self.file_list[date_ts]
            try:
                shutil.copy(file_csv_in,file_csv_out)
            except Exception as ex:
                print(f'[ERROR] File {file_csv_in} could not be copied to {file_csv_out}. Please review permissions')
                return False
        else:
            print(f'[ERROR] CSV file was not found for date {date_ts}')
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
            except:
                print(f'[WARNING] {x} could not be parsed with format {format_datetime}')
                continue

        return only_date_array, time_list

    ##csv with metadata for downloading sources
    def prepare_csv_metadata(self,output_file):
        if not self.check_data():
            return False
        if self.verbose:
            print(f'[INFO] Preparing metatada file: {output_file}')
        path_csv = self.insitu_options['path_csv']
        col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = self.get_csv_options()

        fw = open(output_file,'w')
        started = False
        for name in os.listdir(path_csv):
            if not name.endswith('csv'):
                continue
            file_csv = os.path.join(path_csv, name)
            df = pd.read_csv(file_csv, sep=col_sep)
            only_date_array, time_list = self.get_date_list_from_dataframe(df, col_date, format_date, col_time, format_time)
            only_date_array_unique = np.unique(only_date_array).tolist()
            if len(only_date_array_unique) == 0:
                print(f'[ERROR] No valid dates retrieved from the CSV file.')
                fw.close()
                os.remove(output_file)
                return False
            if len(only_date_array_unique) > 1:
                print(f'[ERROR] Each CSV file should also contain data for a single day. Check CSV parameters')
                fw.close()
                os.remove(output_file)
                return False
            date_ts = only_date_array_unique[0]
            date_here = dt.strptime(date_ts,'%Y-%m-%d')
            if self.start_date is not None and self.end_date is not None:
                if date_here<self.start_date or date_here>self.end_date:
                    if self.verbose:
                        print(f'[INFO] Skipping date {date_here.strftime("%Y-%m-%d")} as it not in the range: {self.start_date.strftime("%Y-%m-%d")} to {self.end_date.strftime("%Y-%m-%d")}')
                    continue


            insitu_lat = np.ma.array(df[col_lat])
            insitu_lon = np.ma.array(df[col_lon])
            lat_min = np.ma.min(insitu_lat)-0.001
            lat_max = np.ma.max(insitu_lat)+0.001
            lon_min = np.ma.min(insitu_lon)-0.001
            lon_max = np.ma.max(insitu_lon)+0.001
            line = f'{date_ts},{lat_min},{lat_max},{lon_min},{lon_max}'
            if started:
                fw.write('\n')
            if not started:
                started = True
            fw.write(line)

        fw.close()
        return True