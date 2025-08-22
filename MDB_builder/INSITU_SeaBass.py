import os
from datetime import datetime as dt
import numpy as np
import pandas as pd
try:
    import SB_support as SBr
except:
    import MDB_builder.SB_support as SBr

class INSITU_SEABASS():

    def __init__(self, insitu_options, verbose):
        self.verbose = verbose
        self.insitu_options = insitu_options

        self.file_list = {}
        self.date_list = []



    def check_data(self):
        if not 'path_seabass' in self.insitu_options:
            print(f'[ERROR] path_seabass is not available in the options dictionary')
            return False
        path_s = self.insitu_options['path_seabass']
        if not os.path.isfile(path_s):
            print(f'[ERROR] {path_s} does not exist or is not a valid file.')
            return False
        var_date, var_time, var_lat, var_lon, format_date, format_time = self.get_seabass_options()
        try:
            sb = SBr.readSB(path_s)
        except Exception as ex:
            print(f'[ERROR] Error reading SeaBass file: {ex}')
            return False
        var_ok = SBr.check_seabass_variables(sb, var_date, var_lat, var_lon)
        if not var_ok:
            print(f'[ERROR] The following variables are avaialable: {[x for x in sb.variables]}')
            return False
        return True

    def prepare_data(self):
        if not self.check_data():
            return False
        path_s = self.insitu_options['path_seabass']
        var_date, var_time, var_lat, var_lon, format_date, format_time = self.get_seabass_options()
        sb = SBr.readSB(path_s)
        date_array, time_list,lat_array,lon_array,indices_array = self.get_matadata_arrays(sb, var_date, var_time, var_lat, var_lon, format_date, format_time)
        if date_array is None:
            return False
        name = os.path.basename(path_s)
        temp_folder = os.path.join(os.path.dirname(path_s),name[0:name.rfind('.')])
        if not os.path.isdir(temp_folder):
            try:
                os.mkdir(temp_folder)
            except:
                print(f'[ERROR] Temporary folder {temp_folder} could not be created. Please review permissions')
                return False

        date_array_unique = np.unique(date_array)
        for date_ts in date_array_unique:
            if self.verbose:
                print(f'[INFO] Getting metadata for date: {date_ts}')

            indices = np.where(date_array==date_ts)
            df = pd.DataFrame(index=np.arange(len(indices[0])),columns=['insitu_index','insitu_time','insitu_lat','insitu_lon'])
            df['insitu_index'] = indices_array[indices]
            df['insitu_time'] = time_list[indices]
            df['insitu_lat'] = lat_array[indices]
            df['insitu_lon'] = lon_array[indices]
            file_csv = os.path.join(temp_folder,f'Metadata_{date_ts}.csv')
            df.to_csv(file_csv,sep=';',index=None)
            self.file_list[date_ts] = file_csv
            self.date_list.append(dt.strptime(date_ts,'%Y-%m-%d'))

        if len(self.date_list) == 0:
            print(f'[ERROR] No valid dates were retrieved from CSV files')
            return False

        self.date_list.sort()

        return True

    def get_ref_date(self, datehere):
        path_s = self.insitu_options['path_seabass']
        name = os.path.basename(path_s)
        name_ref = f'{name[0:name.rfind(".")]}_{datehere.strftime("%Y%m%d")}'
        return name_ref

    def get_metadata_date(self, datehere):
        date_ts = datehere.strftime('%Y-%m-%d')
        file_csv = self.file_list[date_ts]
        df = pd.read_csv(file_csv, sep=';')

        insitu_time = np.ma.array([dt.strptime(x,'%Y-%m-%d %H:%M:%S') for x in df['insitu_time']])
        insitu_lat = np.ma.array(df['insitu_lat'])
        insitu_lon = np.ma.array(df['insitu_lon'])
        insitu_indices = np.ma.array(df['insitu_index'])

        sorted_indices = np.argsort(insitu_time)
        insitu_time = insitu_time[sorted_indices]
        insitu_lat = insitu_lat[sorted_indices]
        insitu_lon = insitu_lon[sorted_indices]
        insitu_indices = (insitu_indices[sorted_indices],)

        return insitu_time, insitu_lat, insitu_lon, insitu_indices

    ##csv with metadata for downloading sources
    def prepare_csv_metadata(self, output_file):
        if not self.check_data():
            return False

        path_s = self.insitu_options['path_seabass']
        var_date, var_time, var_lat, var_lon, format_date, format_time = self.get_seabass_options()
        sb = SBr.readSB(path_s)
        date_array, time_list, lat_array, lon_array, indices_array = self.get_matadata_arrays(sb, var_date, var_time,
                                                                                              var_lat, var_lon,
                                                                                              format_date, format_time)
        if date_array is None:
            return False

        fw = open(output_file, 'w')
        started = False
        date_array_unique = np.unique(date_array)
        for date_ts in date_array_unique:
            if self.verbose:
                print(f'[INFO] Getting metadata for date: {date_ts}')
            indices = np.where(date_array == date_ts)
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
        return True

    def get_matadata_arrays(self,sb,var_date,var_time,var_lat,var_lon,format_date,format_time):
        date_list_orig = sb.data[var_date]
        time_list_orig = sb.data[var_time] if var_time in sb.variables else None
        lat_array_orig = sb.data[var_lat]
        lon_array_orig = sb.data[var_lon]
        ndata = len(date_list_orig)
        indices_orig = np.arange(ndata)
        format_date_time = f'{format_date}T{format_time}'
        date_array = []
        time_list = []
        lat_array = []
        lon_array = []
        indices_array = []
        for idx, x in enumerate(date_list_orig):

            try:
                date_here = dt.strptime(str(x), format_date)
            except:
                print(
                    f'[ERROR] Error parsing dates in SeaBass file: {x} could not be parsed using {format_date} format.')
                print(
                    f'[ERROR] Plase review SEABASS_SELECTION/format_date in the config. file. Expected format: {sb.variables[var_date][1]}')
                return [None] * 2
            # if start_date is not None and end_date is not None:
            #     if date_here < start_date or date_here > end_date:
            #         continue

            lat_array.append(lat_array_orig[idx])
            lon_array.append(lon_array_orig[idx])
            indices_array.append(indices_orig[idx])

            date_array.append(date_here.strftime('%Y-%m-%d'))
            if time_list_orig is None:
                time_list.append(dt.strptime(str(x), format_date))
            else:
                try:
                    val_s = f'{str(x)}T{str(time_list_orig[idx])}'
                    time_list.append(dt.strptime(val_s, format_date_time))
                except:
                    print(
                        f'[ERROR] Error parsing dates in SeaBass file: {x} could not be parsed using {format_date_time} format.')
                    print(
                        f'[ERROR] Plase review SEABASS_SELECTION/format_date in the config. file. Expected format: {sb.variables[var_date][1]}T{sb.variables[var_time][1]}')
                    return [None] * 2

        if len(date_array) == 0:
            print(f'[ERROR] No data was found for the given temporal range')
            return [None] * 2
        date_array = np.array(date_array)
        time_list = np.array(time_list)
        lat_array = np.array(lat_array)
        lon_array = np.array(lon_array)
        indices_array = np.array(indices_array)

        return date_array, time_list, lat_array, lon_array, indices_array

    def get_seabass_options(self):
        var_date = self.insitu_options['var_date']
        var_time = self.insitu_options['var_time']
        var_lat = self.insitu_options['var_lat']
        var_lon = self.insitu_options['var_lon']
        format_date = self.insitu_options['format_date']
        format_time = self.insitu_options['format_time']
        return var_date, var_time, var_lat, var_lon, format_date, format_time
