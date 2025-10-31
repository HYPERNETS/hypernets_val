import os
from datetime import datetime as dt
from netCDF4 import Dataset
import numpy as np
import pandas as pd
try:
    import SB_support as SBr
    import INSITU_base as ISb
except:
    import MDB_builder.SB_support as SBr
    import MDB_builder.INSITU_base as ISb




class INSITU_SEABASS():

    def __init__(self, insitu_options, verbose):

        self.sb = None

        self.fixed_site = False
        self.insitu_type = 'SINGLE_SEABASS'
        self.verbose = verbose
        self.insitu_options = insitu_options
        self.site = self.insitu_options['site']
        self.file_list = {}
        self.date_list = []
        self.start_date = None
        self.end_date = None

        self.basic_arrays = None

        ##RRS & RRS UNC
        self.wl_array = None
        self.nwl = 0
        self.rrs_band_list = None
        self.rrs_unc_band_list = None
        self.rrs_array = None
        self.rrs_unc_array = None
        self.non_spectral_vars = None

    def open_seabass(self):
        if not 'path_seabass' in self.insitu_options:
            print(f'[ERROR] path_seabass is not available in the options dictionary')
            return False
        path_s = self.insitu_options['path_seabass']
        if not os.path.isfile(path_s):
            print(f'[ERROR] {path_s} does not exist or is not a valid file.')
            return False
        if self.sb is None:
            path_s = self.insitu_options['path_seabass']
            try:
                self.sb = SBr.readSB(path_s, no_warn=True)
            except Exception as ex:
                print(f'[ERROR] Error reading SeaBass file: {ex}')
                return False
        return True

    def check_data(self):
        if self.sb is not None:
            return True
        if not self.open_seabass():
            return False
        var_date, var_time, var_lat, var_lon, format_date, format_time = self.get_seabass_options()
        var_ok = SBr.check_seabass_variables(self.sb, var_date, var_lat, var_lon)
        if not var_ok:
            self.sb = None
            return False
        return True

    def get_temp_folder(self):
        path_s = self.insitu_options['path_seabass']
        name = os.path.basename(path_s)
        temp_folder = os.path.join(os.path.dirname(path_s), name[0:name.rfind('.')])
        if not os.path.isdir(temp_folder):
            try:
                os.mkdir(temp_folder)
            except:
                print(f'[ERROR] Temporary folder {temp_folder} could not be created. Please review permissions')
                return None
        return temp_folder

    def prepare_data(self):

        date_array, time_list,lat_array,lon_array,indices_array = self.get_metadata_arrays()
        if date_array is None:
            return False
        temp_folder = self.get_temp_folder()

        if temp_folder is None:
            return False

        date_array_unique = np.unique(date_array)
        for date_ts in date_array_unique:
            date_ts_obj = dt.strptime(date_ts, '%Y-%m-%d')
            if self.start_date is not None and self.end_date is not None:
                if date_ts_obj<self.start_date or date_ts_obj>self.end_date:
                    continue
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
            self.date_list.append(date_ts_obj)

        if len(self.date_list) == 0:
            print(f'[ERROR] No valid dates were retrieved from CSV files')
            return False

        self.date_list.sort()

        return True

    def get_rrs_array(self):
        if self.rrs_array is not None:
            return True
        if not self.check_data():
            return False
        if self.nwl == 0:
            if not self.check_rrs_and_data_variables():
                return False

        self.rrs_array = np.zeros((self.sb.length, self.nwl))
        for iwl, band in enumerate(self.rrs_band_list):
            self.rrs_array[:, iwl] = np.array(self.sb.data[band])
        if self.verbose:
            print(f'[INFO] Rrs data were extracted')
        return True

    def get_rrs_unc_array(self):
        if self.insitu_options['rrs_unc_format'] is None:
            return False
        if self.rrs_unc_array is not None:
            return True
        if not self.check_data():
            return False
        if self.nwl == 0:
            if not self.check_rrs_and_data_variables():
                return False

        self.rrs_unc_array = np.zeros((self.sb.length, self.nwl))
        for iwl, band in enumerate(self.rrs_unc_band_list):
            self.rrs_unc_array[:, iwl] = np.array(self.sb.data[band])
        if self.verbose:
            print(f'[INFO] Rrs uncentainty data were extracted')
        return True


    def check_rrs_and_data_variables(self):
        if not self.check_data():
            return False

        col_names = list(self.sb.variables.keys())


        ##RRS
        if self.insitu_options['rrs_format'] is not None: ##RRS info is available
            self.rrs_band_list, self.wl_array = self.check_rrs_list(self.insitu_options['rrs_format'],self.insitu_options['rrs_bandlist'],self.insitu_options['rrs_list'],col_names)
            if self.rrs_band_list is None or self.wl_array is None:
                return False
            self.nwl = len(self.wl_array)

            if self.verbose:
                print(f'[INFO] {self.nwl} in situ Rrs bands identified in the SeaBass file')


        ##RRS UNCERTAINTY
        if self.insitu_options['rrs_unc_format'] is not None:
            self.rrs_unc_band_list, wl_array_unc = self.check_rrs_list(self.insitu_options['rrs_unc_format'],self.insitu_options['rrs_unc_bandlist'],self.insitu_options['rrs_list'],col_names)
            if self.rrs_unc_band_list is None or wl_array_unc is None:
                return False
            if self.nwl!=len(wl_array_unc):
                print(f'[ERROR] Rrs uncentainties are available for {len(wl_array_unc)} band but Rrs are available for {self.nwl} bands')
                return False
            if not (self.wl_array==wl_array_unc).all():
                print(f'[ERROR] Rrs wavelengths are different for Rrs uncentainty wavelengths')
                return False
            if self.verbose:
                print(f'[INFO] {self.nwl} in situ Rrs uncentainty bands identified in the SeaBass file')

        if self.insitu_options['non_spectral_variables'] is not None:
            self.non_spectral_vars = {}
            for var_name in self.insitu_options['non_spectral_variables']:
                var_name = var_name.lower()
                if var_name in col_names:
                    try:
                        self.non_spectral_vars[var_name] = self.sb.variables[var_name][1]
                    except:
                        self.non_spectral_vars[var_name] = None
            if len(self.non_spectral_vars)==0:
                self.non_spectral_vars = None


        return True

    def check_rrs_list(self,rrs_format,rrs_band_list,rrs_list,col_names):
        wl = []
        if rrs_band_list is not None:
            if rrs_list is not None:
                if len(rrs_list) != len(rrs_band_list):
                    print(f'[ERROR] Inconsistency between SeaBass options for rrs_list and rrs_bandlist (or rrs_unc_band_list). Both lists have a different number of elements')
                    return [None]*2
                for rrs, band in zip(rrs_list,rrs_band_list):
                    try:
                        wl.append(float(rrs))
                    except:
                        print(f'[ERROR] {rrs} value in rrs_list is not a valid number')
                        return [None]*2
                    expected_band = rrs_format.replace('$BAND$', rrs)
                    if expected_band != band:
                        expected_band = expected_band.replace('.', '_')
                    if expected_band != band:
                        print(f'[ERROR] Inconsistency between SeaBass options rrs_list and rrs_bandlist (or rrs_unc_band_list). Expected band should be {rrs_format.replace("$BAND$", rrs)}, but {band} was found')
                        return [None]*2
            else:
                prefix = rrs_format[0:rrs_format.index('$BAND$')]
                suffix = rrs_format[rrs_format.index('$BAND$') + 6:]
                for band in rrs_band_list:
                    rrs = band.replace(prefix, '')
                    rrs = rrs.replace(suffix, '')
                    try:
                       wl.append(float(rrs))
                    except:
                        print(f'[ERROR] {rrs} value in rrs_band_list is not a valid number')
                        return [None]*2

        if rrs_band_list is None:
            rrs_band_list = []
            if rrs_list is not None:
                for rrs in rrs_list:
                    try:
                        rrs_band_list.append(rrs_format.replace('$BAND$',rrs))
                        wl.append(float(rrs))
                    except:
                        print(f'[ERROR] {rrs} value in rrs_list is not a valid number')
                        return [None]*2
            else:

                prefix = rrs_format[0:rrs_format.index('$BAND$')]
                suffix = rrs_format[rrs_format.index('$BAND$') + 6:]

                for col in col_names:
                    if col.startswith(prefix) and col.endswith(suffix):
                        iini = len(prefix)
                        ifin = len(col)-len(suffix)
                        val = col[iini:ifin].replace('_','.')
                        try:
                            wl.append(float(val))
                            rrs_band_list.append(col)
                        except:
                            continue

        if len(rrs_band_list)!=len(wl):
            print(f'[ERROR] Inconsistency in the lenghts of rrs_band_list and wl')
            return [None]*2
        if len(rrs_band_list)==0:
            print(f'[ERROR] No bands found using format {rrs_format}')
            return [None]*2
        if not self.check_exist_variables(rrs_band_list,col_names):
            return [None]*2
        wl = np.array(wl)
        return rrs_band_list,wl


    def check_exist_variables(self,var_to_check,col_names):
        if not self.check_data():
            return False
        if col_names is None:
            col_names = list(self.sb.variables.keys())
        for var in var_to_check:
            if var not in col_names:
                print(f'[ERROR] Variable {var} is not valid variables in the SeaBass file')
                return False

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

        date_array, time_list, lat_array, lon_array, indices_array = self.get_metadata_arrays()
        if date_array is None:
            return False

        fw = open(output_file, 'w')
        started = False
        date_array_unique = np.unique(date_array)
        for date_ts in date_array_unique:
            date_here = dt.strptime(date_ts,'%Y-%m-%d')
            if self.start_date is not None and self.end_date is not None:
                if date_here<self.start_date or date_here>self.end_date:
                    if self.verbose:
                        print(f'[INFO][SeaBass] Skipping date {date_here.strftime("%Y-%m-%d")} as it not in the range: {self.start_date.strftime("%Y-%m-%d")} to {self.end_date.strftime("%Y-%m-%d")}')
                    continue
            if self.verbose:
                print(f'[INFO][SeaBass] Getting metadata for date: {date_ts}')
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
        if not started:
            print(f'[WARNING][SeaBass] No data were found for the given data range in the SeaBass dataset.')
            os.remove(output_file)
            return False
        else:
            return True

    def create_csv_metadata_for_date(self,insitu_date, file_csv_out):
        if self.basic_arrays is None:
            date_array, time_list, lat_array, lon_array, indices_array = self.get_metadata_arrays()
            self.basic_arrays = {
                'date':date_array,
                'time':time_list,
                'lat':lat_array,
                'lon':lon_array
            }
        else:
            date_array = self.basic_arrays['date']
            time_list = self.basic_arrays['time']
            lat_array = self.basic_arrays['lat']
            lon_array = self.basic_arrays['lon']

        if date_array is None:
            return False

        date_ts = insitu_date.strftime('%Y-%m-%d')
        indices = np.where(date_array == date_ts)
        ndata = len(indices[0])
        if ndata==0:
            return False
        insitu_lat = lat_array[indices]
        insitu_lon = lon_array[indices]
        insitu_time = [x.strftime('%Y-%m-%d %H:%M:%S.%f') for x in time_list[indices]]
        df = pd.DataFrame(index=np.arange(ndata),columns=['date','lat','lon'])
        df.loc[:,'lat'] = insitu_lat[:]
        df.loc[:,'lon'] = insitu_lon[:]
        df.loc[:,'date'] = insitu_time[:]
        df.to_csv(file_csv_out,sep=';',index=None)

    def get_metadata_arrays(self):

        if not self.check_data():
            return [None]*5
        var_date, var_time, var_lat, var_lon, format_date, format_time = self.get_seabass_options()
        date_list_orig = self.sb.data[var_date]
        time_list_orig = self.sb.data[var_time] if var_time in self.sb.variables else None
        lat_array_orig = self.sb.data[var_lat]
        lon_array_orig = self.sb.data[var_lon]
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
                    f'[ERROR] Plase review SEABASS_SELECTION/format_date in the config. file. Expected format: {self.sb.variables[var_date][1]}')
                return [None] * 5
            if self.start_date is not None and self.end_date is not None:
                if date_here<self.start_date or date_here>self.end_date:
                    continue
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
                        f'[ERROR] Plase review SEABASS_SELECTION/format_date in the config. file. Expected format: {self.sb.variables[var_date][1]}T{self.sb.variables[var_time][1]}')
                    return [None] * 5

        if len(date_array) == 0:
            print(f'[ERROR] No data was found for the given temporal range')
            return [None] * 5
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

    def create_mini_mdb_files(self,options,extract_dir,extracts,time_extracts,overwrite):
        options['n_insitu_bands'] = self.nwl
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
        nwlsat = len(dset.dimensions['satellite_bands'])
        nwlinsitu = len(dset.dimensions['insitu_bands'])
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


        #insitu_indices = extract_info['insitu_indices']
        builder = ISb.Mini_MDB_Builder(options,self.verbose)
        builder.start_mini_mdb(extract_info['file'],file_out)
        builder.add_shipborne_variables()
        if self.insitu_options['rrs_unc_format'] is not None:
            builder.add_rrs_uncentainty_variable()

        building_error = False
        if len(options['instrument_ids'])==1:
            builder.set_insitu_wavelengths(0,self.wl_array)
            builder.set_instrument_id(ninsitu_real,1)
        #basic variables: insitu_lat,insitu_lon,insitu_time,time_diff,insitu_spatial_index

        if not builder.set_insitu_basic_variables_from_dict(extract_info):
            print(f'[ERROR] Error setting the basic variables. MDBm file could not be created')
            building_error = True

        indices = extract_info['insitu_indices']
        if self.get_rrs_array():
            rrs_here = np.ma.transpose(self.rrs_array[indices,:])
            if not builder.set_spectral_variables('insitu_Rrs',rrs_here):
                print(f'[ERROR] Error setting the in situ Rrs variable. MDBm file could not be created')
                building_error = True
        else:
            print(f'[ERROR] Error retrieving the Rrs variable from the SeaBass file')
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
                attrs = {'units': self.non_spectral_vars[var_name]}
                builder.add_non_spectral_variable(var_name_nc, attrs)
                array = np.array(self.sb.data[var_name])
                array = array[indices]
                builder.set_non_spectral_variables(var_name_nc,array)

        builder.close_mini_mdb_file()
        if building_error:
            os.remove(file_out)
