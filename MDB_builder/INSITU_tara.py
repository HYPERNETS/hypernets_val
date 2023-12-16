import pytz
import numpy as np

from INSITU_base import INSITUBASE
from datetime import datetime as dt
import os
import pandas as pd
import COMMON.common_functions as cfs

class INSITU_TARA(INSITUBASE):

    def __init__(self, mdb_options, verbose):
        if mdb_options is not None:
            self.path_tara_metadata = mdb_options.insitu_path_metadata
            self.path_tara_data = mdb_options.insitu_path_source
        else:
            self.path_tara_metadata = '/mnt/c/DATA_LUIS/TARA_TEST/insitu_data/So-Rad_meta'

        self.mdb_options = mdb_options
        self.verbose = verbose
        self.new_MDB = None

        self.metadata_options = {
            'path_tara_metadata_org':None,
            'file_tara_metadata_name':'$DATE$_3C_metadata.csv',
            'file_tara_metadata_name_date_format': '%Y-%m-%d',
            'col_sep':',',
            'col_lat':'lat',
            'col_lon':'lon',
            'col_time':'timestamp',
            'format_time':'%Y-%m-%d %H:%M:%S.%f',
            'bands':{
                'gps_speed':{
                    'type':'f4'
                },
                'tilt_avg': {
                    'type': 'f4'
                },
                'tilt_std': {
                    'type': 'f4'
                },
                'rel_view_az':{
                    'type': 'f4'
                },
                'q_1':{
                    'type':'i1',
                },
                'q_2':{
                    'type':'i1'
                },
                'q_3':{
                    'type':'i1'
                }

            }
        }

        self.metadata = {}

    def check_match_up_conditions(self,sat_time,lat_array,lon_array,max_time_diff):
        window_size = lat_array.shape[0]
        ref = np.floor(window_size/2)
        for index in self.metadata:
            in_situ_lat = self.metadata[index]['lat']
            in_situ_lon = self.metadata[index]['lon']
            if cfs.contain_location(lat_array,lon_array,in_situ_lat,in_situ_lon):
                r, c = cfs.find_row_column_from_lat_lon(lat_array, lon_array, in_situ_lat, in_situ_lon)
                ipos = max([abs(r-ref),abs(c-ref)])
                print(index,in_situ_lat,in_situ_lon,r,c,ipos)
                self.metadata_options['index_spatial'] = ipos


    def retrieve_metadata_from_file(self,date_here):
        self.metadata = {}
        file_metadata = self.get_file_metadata(date_here)
        if file_metadata is None:
            return
        try:
            df = pd.read_csv(file_metadata,sep=self.metadata_options['col_sep'])
        except:
            return

        for index,row in df.iterrows():
            self.metadata[index] = {
                'lat': row[self.metadata_options['col_lat']],
                'lon': row[self.metadata_options['col_lon']],
                'time': dt.strptime(row[self.metadata_options['col_time']],self.metadata_options['format_time']).replace(tzinfo=pytz.utc),
                'index_spatial':-1,
                'index_temporal':-1,
                'time_diff':-1
            }
            for b in self.metadata_options['bands']:
                self.metadata[index][b] = float(row[b])

            print(self.metadata[index])


    def get_file_metadata(self,date_here):
        path_date = self.get_path_date(date_here)
        if path_date is None:
            return None

        name_file = self.metadata_options['file_tara_metadata_name']
        date_format = self.metadata_options['file_tara_metadata_name_date_format']
        name_file = name_file.replace('$DATE$',date_here.strftime(date_format))
        file_metadata = os.path.join(path_date,name_file)
        if not os.path.isfile(file_metadata):
            return None
        return file_metadata

    def get_path_date(self,date_here):
        org = self.metadata_options['path_tara_metadata_org']
        path_date = self.path_tara_metadata
        if org is not None:
            org_list = org.split('/')
            for ol in org_list:
                path_date = os.path.join(path_date,date_here.strftime(ol))
        if not os.path.isdir(path_date):
            return None
        return path_date


