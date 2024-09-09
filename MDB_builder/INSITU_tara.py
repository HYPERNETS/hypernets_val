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
            self.path_tara_metadata = '/mnt/c/DATA_LUIS/TARA_TEST/insitu_data/So-Rad_meta_V2'




        self.mdb_options = mdb_options



        self.verbose = verbose
        self.new_MDB = None

        self.data_options = {
            'path_tara_data_org': None,
            'file_tara_data_name': '$DATE$_3C_Rrs.csv',
            'file_tara_data_name_date_format': '%Y-%m-%d',
            'col_sep': ',',
        }

        self.metadata_options = {
            'path_tara_metadata_org': None,
            'file_tara_metadata_name': '$DATE$_3C_metadata_V2.csv',
            'file_tara_metadata_name_date_format': '%Y-%m-%d',
            'col_sep': ',',
            'col_lat': 'lat',
            'col_lon': 'lon',
            'col_time': 'timestamp',
            'format_time': '%Y-%m-%d %H:%M:%S.%f',
            'bands': {
                'gps_speed': {
                    'type': 'f4',
                    'long_name': 'GPS speed',
                    'units': 'knots',
                    'available': True
                },
                'tilt_avg': {
                    'type': 'f4',
                    'long_name': 'Average tilt angle',
                    'units': 'degrees',
                    'available': True
                },
                'tilt_std': {
                    'type': 'f4',
                    'long_name': 'Standard deviation of tilt angle',
                    'units': 'degrees',
                    'available': True
                },
                'rel_view_az': {
                    'type': 'f4',
                    'long_name': 'Relative view azimuth angle',
                    'units': 'degrees',
                    'available': True
                },
                'q_1': {
                    'type': 'i1',
                    'long_name': 'Quality Mask 1 (0: INVALID, 1: VALID)',
                    'available': True
                },
                'q_2': {
                    'type': 'i1',
                    'long_name': 'Quality Mask 2 (0: INVALID, 1: VALID)',
                    'available': True
                },
                'q_3': {
                    'type': 'i1',
                    'long_name': 'Quality Mask 3 (0: INVALID, 1: VALID)',
                    'available': True
                },
                'q_4': {
                    'type': 'i1',
                    'long_name': 'Quality Mask 4 (0: INVALID, 1: VALID)',
                    'available': True
                },
                'q_5': {
                    'type': 'i1',
                    'long_name': 'Quality Mask 5 (0: INVALID, 1: VALID)',
                    'available': True
                }

            }
        }


        self.update_options()

        #print(self.metadata_options)
        self.metadata = {}

        self.wavelengths = None

    def update_options(self):
        if not self.mdb_options.options.has_section('tara_options'):
            return
        for option in self.data_options:
            if self.mdb_options.options.has_option('tara_options',option):
                self.data_options[option] = self.mdb_options.options['tara_options'][option]
        for option in self.metadata_options:
            if self.mdb_options.options.has_option('tara_options', option):
                self.metadata_options[option] = self.mdb_options.options['tara_options'][option]

    def create_mdb_insitu_extract(self, extract_path, ofile, extract_time, max_time_diff):
        from netCDF4 import Dataset
        dextract = Dataset(extract_path)
        lat_array = dextract.variables['satellite_latitude'][0, :, :]
        lon_array = dextract.variables['satellite_longitude'][0, :, :]
        dextract.close()
        npoints = self.check_match_up_conditions(extract_time, lat_array, lon_array, max_time_diff)
        if self.verbose:
            print(
                f'[INFO] Total number of in situ data points in the area: {npoints[0]} Within the temporal window: {npoints[1]}')
            print(
                f'[INFO] Total number of in situ data points in the central pixel: {npoints[2]} Within the temporal window: {npoints[3]}')

        if npoints[3] == 0:
            print(
                f'[WARNING] No in situ data points in the central pixel and within the temporal window. Skipping extract...')
            return False

        file_data = self.get_file_data(extract_time)
        if file_data is None:
            print(f'[ERROR] No data file was available for time: {extract_time}')
            return False

        df = pd.read_csv(file_data, sep=self.data_options['col_sep'])
        if self.wavelengths is None:
            wl_values = []
            for col in df.columns:
                if col.startswith('#'):
                    val = float(col[1:].strip())
                else:
                    val = float(col)
                wl_values.append(val)
            self.wavelengths = np.array(wl_values)

        self.mdb_options.insitu_options['n_insitu_bands'] = len(self.wavelengths)

        self.start_add_insitu(extract_path, ofile)
        self.add_shipborne_variables()
        if 'bands' in self.metadata_options:
            for band in self.metadata_options['bands']:
                if not self.metadata_options['bands'][band]['available']:
                    continue
                ats = {}
                for at in self.metadata_options['bands'][band]:
                    if at == 'type' or at=='available':
                        continue

                    ats[at] = self.metadata_options['bands'][band][at]

                self.add_insitu_variable(band, self.metadata_options['bands'][band]['type'], ats)

        max_index_spatial = 12
        self.new_MDB.variables['insitu_original_bands'][:] = self.wavelengths[:]

        ihere = 0
        error_metadata_data = False
        ndata_added = 0
        for index in self.metadata:
            index_spatial = self.metadata[index]['index_spatial']
            index_temporal = self.metadata[index]['index_temporal']
            # print(index_temporal,index_temporal)
            if max_index_spatial >= index_spatial >= 0 and index_temporal == 1:
                try:
                    rrs_here = df.iloc[index].values
                except:
                    error_metadata_data = True
                    break
                ndata_added = ndata_added + 1
                self.new_MDB.variables['insitu_Rrs'][0, :, ihere] = rrs_here[:]
                self.new_MDB.variables['insitu_time'][0, ihere] = np.float64(self.metadata[index]['time'].timestamp())
                self.new_MDB.variables['insitu_latitude'][0, ihere] = self.metadata[index]['lat']
                self.new_MDB.variables['insitu_longitude'][0, ihere] = self.metadata[index]['lon']
                self.new_MDB.variables['time_difference'][0, ihere] = self.metadata[index]['time_diff']
                self.new_MDB.variables['insitu_spatial_index'][0, ihere] = self.metadata[index]['index_spatial']
                for b in self.metadata_options['bands']:
                    if not self.metadata_options['bands'][b]['available']:
                        continue
                    name_band = b
                    if not b.startswith('insitu_'):
                        name_band = f'insitu_{b}'
                    self.new_MDB.variables[name_band][0, ihere] = self.metadata[index][b]
                ihere = ihere + 1

        print(f'[INFO] Added spectra: {ndata_added}')
        if ndata_added!=npoints[1]:
            print(f'[ERROR] Number of spectra added should be equal to {npoints[1]}---------------------------------------------------------------------------***')
            return False




        if error_metadata_data:
            self.new_MDB.close()
            os.remove(ofile)
            print(f'[ERROR] Detected disagreement between data and metadata files. Skipped sat. extract')
            return False

        self.new_MDB.close()

        return True

    def check_match_up_conditions(self, sat_time, lat_array, lon_array, max_time_diff):
        npoints = [0, 0, 0, 0]  # nall,nalltime,ncentral,ncentraltime
        window_size = lat_array.shape[0]
        ref = np.floor(window_size / 2)
        for index in self.metadata:
            in_situ_lat = self.metadata[index]['lat']
            in_situ_lon = self.metadata[index]['lon']
            in_situ_time = self.metadata[index]['time']
            self.metadata[index]['index_spatial'] = -1
            self.metadata[index]['index_temporal'] = -1
            if cfs.contain_location(lat_array, lon_array, in_situ_lat, in_situ_lon):
                r, c = cfs.find_row_column_from_lat_lon(lat_array, lon_array, in_situ_lat, in_situ_lon)
                ipos = max([abs(r - ref), abs(c - ref)])
                self.metadata[index]['index_spatial'] = ipos

                sat_time = sat_time.replace(tzinfo=pytz.utc)
                in_situ_time = in_situ_time.replace(tzinfo=pytz.utc)
                # print(in_situ_time,type(in_situ_time))
                time_diff_here = abs((in_situ_time - sat_time).total_seconds())
                npoints[0] = npoints[0] + 1
                if ipos == 0:
                    npoints[2] = npoints[2] + 1
                if time_diff_here <= max_time_diff:
                    self.metadata[index]['index_temporal'] = 1
                    npoints[1] = npoints[1] + 1
                    if ipos == 0:
                        npoints[3] = npoints[3] + 1
                else:
                    self.metadata[index]['index_temporal'] = 0
                self.metadata[index]['time_diff'] = time_diff_here
                # print(index, in_situ_lat, in_situ_lon, r, c, ipos,sat_time,in_situ_time,time_diff_here,self.metadata[index]['index_temporal'])

            #print(index,self.metadata[index]['index_spatial'],self.metadata[index]['index_temporal'])

        return npoints

    def retrieve_metadata_from_file(self, date_here):
        self.metadata = {}
        file_metadata = self.get_file_metadata(date_here)
        if file_metadata is None:
            return
        try:
            df = pd.read_csv(file_metadata, sep=self.metadata_options['col_sep'])
        except:
            return
        col_names = df.columns.tolist()

        for index, row in df.iterrows():
            #print(row)
            self.metadata[index] = {
                'lat': row[self.metadata_options['col_lat']],
                'lon': row[self.metadata_options['col_lon']],
                'time': dt.strptime(row[self.metadata_options['col_time']],
                                    self.metadata_options['format_time']).replace(tzinfo=pytz.utc),
                'index_spatial': -1,
                'index_temporal': -1,
                'time_diff': -1
            }
            for b in self.metadata_options['bands']:
                if b in col_names:
                    self.metadata[index][b] = float(row[b])
                else:
                    self.metadata_options['bands'][b]['available'] = False

            # print(self.metadata[index])

    # def update_metadata(self):
    #     if len(self.metadata)==0:
    #         return
    #     for index in self.metadata:
    #         self.metadata[index]['index_spatial'] = -1
    #         self.metadata[index]['index_temporal'] = -1
    #         self.metadata[index]['time_diff'] = -1

    def get_time_distribution_from_metadata(self,date_here,date_ini,date_fin):
        self.retrieve_metadata_from_file(date_here)
        if len(self.metadata)==0:
            print('no data')
            return None
        date_ref = date_here.replace(hour=0,minute=0,second=0,tzinfo=pytz.UTC)
        # min_times = np.linspace(0,85500,96)
        # max_times = np.linspace(900,86400,96)
        distributions = np.zeros((96,4)) #0: total #1: q_3, #2:q4 #3 q5
        summary = np.zeros((4,))
        for index in self.metadata:
            withinrange = False
            itime = self.metadata[index]['time']
            if date_ini <= itime <= date_fin:
                withinrange = True
                summary[0] = summary[0]+1
            dif = (itime-date_ref).total_seconds()
            index_h = int(np.floor(dif/900.0))
            distributions[index_h,0] = distributions[index_h,0]+1

            if self.metadata[index]['q_3']==1:
                distributions[index_h, 1] = distributions[index_h, 1] + 1
                if withinrange:
                    summary[1] = summary[1] + 1
            if self.metadata[index]['q_4']==1:
                distributions[index_h, 2] = distributions[index_h, 2] + 1
                if withinrange:
                    summary[2] = summary[2] + 1
            if self.metadata[index]['q_5']==1:
                distributions[index_h, 3] = distributions[index_h, 3] + 1
                if withinrange:
                    summary[3] = summary[3] + 1

        return distributions,summary


    def get_file_metadata(self, date_here):
        org = self.metadata_options['path_tara_metadata_org']
        path_date = self.get_path_date(date_here, org, self.path_tara_metadata)
        if path_date is None:
            return None

        name_file = self.metadata_options['file_tara_metadata_name']
        date_format = self.metadata_options['file_tara_metadata_name_date_format']

        name_file = name_file.replace('$DATE$', date_here.strftime(date_format))
        file_metadata = os.path.join(path_date, name_file)
        if not os.path.isfile(file_metadata):
            return None
        return file_metadata

    def get_file_data(self, date_here):
        org = self.data_options['path_tara_data_org']
        path_date = self.get_path_date(date_here, org, self.path_tara_data)
        if path_date is None:
            return None

        name_file = self.data_options['file_tara_data_name']
        date_format = self.data_options['file_tara_data_name_date_format']
        name_file = name_file.replace('$DATE$', date_here.strftime(date_format))
        file_data = os.path.join(path_date, name_file)
        if not os.path.isfile(file_data):
            return None
        return file_data

    def get_path_date(self, date_here, org, path_date):
        if org is not None:
            org_list = org.split('/')
            for ol in org_list:
                path_date = os.path.join(path_date, date_here.strftime(ol))
        if not os.path.isdir(path_date):
            return None
        return path_date
