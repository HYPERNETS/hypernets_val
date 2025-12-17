import os
import pandas as pd
import numpy as np
from datetime import datetime as dt


class INSITU_SINGLE_CSV():

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
        print(insitu_options)

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
        #only_date_array = np.array(only_date_array)
        if len(only_date_array)==0:
            print(f'[ERROR] No valid dates retrieved from the CSV file.')
            return False
        print('98')
        only_date_array_unique = np.unique(only_date_array).tolist()
        print('100')
        if len(only_date_array_unique) == 0:
            print(f'[ERROR] No valid dates retrieved from the CSV file.')
            return False
        print('76')
        for date_ts in only_date_array_unique:
            self.file_list[date_ts] = file_csv
            self.date_list.append(dt.strptime(date_ts,'%Y-%m-%d'))

        if len(self.date_list)==0:
            print(f'[ERROR] No valid dates were retrieved from CSV files')
            return False
        self.date_list.sort()
        print('85')
        return True