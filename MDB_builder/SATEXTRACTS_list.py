import os,pytz,sys,__init__
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from datetime import datetime as dt
code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
import COMMON.common_functions as cfs


class SAT_EXTRACTS_LIST:
    def __init__(self, boptions, verbose):

        self.boptions = boptions
        self.verbose = verbose

        self.set_time_default = False
        self.hour_default = 11
        self.minute_default = 0

        if self.boptions.insitu_options['time_sat_default'] is not None:

            str_time = self.boptions.insitu_options['time_sat_default']
            try:
                time = dt.strptime(str_time,'%H:%M')
                self.set_time_default = True
                self.hour_default = int(time.hour)
                self.minute_default  = int(time.minute)
                print(f'[INFO] Time sat default set to: {str_time}')
            except:
                pass

        self.file_name_format = self.boptions.param_sat['file_name_format']
        self.file_name_date_format = self.boptions.param_sat['file_name_date_format']

        self.path_org = self.boptions.param_sat['path_org']

    def creating_copy_correction_sattime(self,input_file,output_file,sat_time_new):
        from datetime import timezone
        input_dataset = Dataset(input_file)
        ncout = Dataset(output_file, 'w', format='NETCDF4')

        # copy global attributes all at once via dictionary
        ncout.setncatts(input_dataset.__dict__)

        # copy dimensions
        for name, dimension in input_dataset.dimensions.items():
            ncout.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))

        for name, variable in input_dataset.variables.items():
            fill_value = None
            if '_FillValue' in list(variable.ncattrs()):
                fill_value = variable._FillValue

            ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                                 shuffle=True, complevel=6)

            # copy variable attributes all at once via dictionary
            ncout[name].setncatts(input_dataset[name].__dict__)

            # copy data
            if name == 'satellite_time':
                ftime = sat_time_new.replace(tzinfo=timezone.utc).timestamp()
                ncout[name][:] = ftime

            else:
                ncout[name][:] = input_dataset[name][:]

        ncout.close()

        input_dataset.close()

    def get_list_as_dict(self):
        start_date = self.boptions.start_date
        end_date = self.boptions.end_date
        site = self.boptions.param_insitu

        prefix = self.boptions.param_sat['prefix']

        sat_extract_dir = self.boptions.satellite_path_source

        sat_list = {}

        # print(self.boptions.param_sat)
        nadded = 0

        for name in os.listdir(sat_extract_dir):
            if not name.endswith('.nc'):
                continue
            if prefix is not None and not name.startswith(prefix):
                continue

            fextract = os.path.join(sat_extract_dir, name)
            if self.verbose:
                print(f'[INFO] Checking extract file: {name}')

            try:
                dataset = Dataset(fextract)
            except:
                print(f'[WARNING] Extract {name} is not a valid NetCDF file. Skipping...')
                continue
            time_here,correct_time = self.check_time(name, dataset, start_date, end_date)

            if time_here is None:
                dataset.close()
                continue

            site_here = self.check_site(name, dataset, site['station_name'])
            if site_here is None:
                dataset.close()
                continue
            platform_name = self.boptions.param_sat['satellite'] + self.boptions.param_sat['platform']
            platform = self.check_platform(dataset, platform_name)
            if platform is None:
                dataset.close()
                continue
            ac = self.check_atm_correction(dataset, self.boptions.param_sat['ac'])
            if ac is None:
                dataset.close()
                continue
            sensor = self.check_sensor(dataset, self.boptions.param_sat['sensor'])
            if sensor is None:
                dataset.close()
                continue
            resolution = self.check_resolution(name, dataset, self.boptions.param_sat['resolution'])
            if resolution is None:
                dataset.close()
                continue
            dataset.close()

            if correct_time:
                output_file = os.path.join(sat_extract_dir,name[:-3]+'_COPY.NC')
                self.creating_copy_correction_sattime(fextract,output_file,time_here)
                os.rename(output_file,fextract)
                dataset = Dataset(fextract)
                ftime = dataset.variables['satellite_time']
                time_check = dt.utcfromtimestamp(float(dataset.variables['satellite_time'][0]))
                dataset.close()
                print(f'[WARNING] Time was corrected and set as UTC to: {time_check}')


            time_ref = time_here.strftime('%Y%m%d')
            if site_here=='SHIPBORNE': #more than one extract

                time_ref_base = time_ref


                if time_ref_base not in sat_list:
                    time_ref = f'{time_ref_base}_1'
                    sat_list[time_ref_base] = {
                        time_ref: {
                            'path': fextract,
                            'time': time_here,  # .strftime('%Y%m%dT%H%M%S'),
                            'site': site_here,
                            'sensor': sensor,
                            'platform': platform,
                            'ac': ac,
                            'resolution': resolution
                        }
                    }
                else:
                    index = 1
                    time_ref = f'{time_ref_base}_{index}'
                    while time_ref in sat_list[time_ref_base].keys():
                        index = index + 1
                        time_ref = f'{time_ref_base}_{index}'
                    sat_list[time_ref_base][time_ref] = {
                        'path': fextract,
                        'time': time_here,  # .strftime('%Y%m%dT%H%M%S'),
                        'site': site_here,
                        'sensor': sensor,
                        'platform': platform,
                        'ac': ac,
                        'resolution': resolution
                    }


            else:
                sat_list[time_ref] = {
                    'path': fextract,
                    'time': time_here,#.strftime('%Y%m%dT%H%M%S'),
                    'site': site_here,
                    'sensor': sensor,
                    'platform': platform,
                    'ac': ac,
                    'resolution': resolution
                }
            nadded = nadded + 1
        if self.verbose:
            print(f'[INFO] Number of extract files added to the list: {nadded} ')
        return sat_list

    def check_time(self,fname,dataset,start_date,end_date):
        correct = False
        datetime_here = dt.utcfromtimestamp(float(dataset.variables['satellite_time'][0]))
        if self.file_name_format is not None and self.file_name_date_format is not None:
            istart = self.file_name_format.index('$DATE$')
            iend = istart + len(dt.now().strftime(self.file_name_date_format))
            if istart>0 and iend>0:
                datetime_here_name_str = fname[istart:iend]
                try:
                    datetime_here_name = dt.strptime(datetime_here_name_str,self.file_name_date_format).astimezone(pytz.utc)
                    tdif = abs((datetime_here-datetime_here_name).total_seconds())
                    if tdif>120:
                        correct = True
                        datetime_here = datetime_here_name
                except:
                    pass

        #if datetime_here.hour == 0 and datetime_here.minute == 0 and datetime_here.second == 0:
        if self.set_time_default:
            datetime_here = datetime_here.replace(hour=self.hour_default, minute=self.minute_default)
        if start_date <= datetime_here <= end_date:
            return datetime_here,correct
        else:
            return None,correct

    def check_time_deprecated(self, fname, dataset, start_date, end_date):
        datetime_here = dt.utcfromtimestamp(float(dataset.variables['satellite_time'][0]))

        try:
            datetime_here_name = dt.strptime(fname.split('_')[7], '%Y%m%dT%H%M%S')
            if datetime_here_name > datetime_here:
                datetime_here = datetime_here_name
        except:
            try:
                datetime_here_name = dt.strptime(fname[1:8], '%Y%j')
                if datetime_here_name > datetime_here:
                    datetime_here = datetime_here_name
            except:
                pass

        if datetime_here.hour == 0 and datetime_here.minute == 0 and datetime_here.second == 0:
            datetime_here = datetime_here.replace(hour=self.hour_default, minute=self.minute_default)
        if start_date <= datetime_here <= end_date:
            return datetime_here
        else:
            return None

    def check_resolution(self, fname, dataset, res_name):
        res_here = None
        if 'resolution' in dataset.ncattrs():
            res_here = dataset.resolution
        if res_here is None:
            if fname.find(res_name) > 0:
                return res_name
            else:
                print(f'[WARNING] Resolution set to {res_name} despite of not being defined in the extract file')
                return res_name
        else:
            if res_here.upper() == res_name.upper():
                return res_name
            else:
                print(
                    f'[WARNING] Extract resolution {res_here} was not selected in the config file. Skipping extract...')
                return None

    def check_sensor(self, dataset, sensor_name):
        sensor_here = None
        if 'sensor' in dataset.ncattrs():
            sensor_here = dataset.sensor
            if sensor_here.startswith('MODIS Moderate Resolution Imaging Spectroradiometer,'):
                sensor_here = 'MULTI'
            if sensor_here.startswith('SeaWiFS,'):
                sensor_here = 'MULTI'
            if sensor_here.startswith('MSI'):
                sensor_here = 'MSI'

        if sensor_here is None:
            print(f'[WARNING] Sensor set to {sensor_name} despite of not being defined in the extract file')
            return sensor_name
        else:
            if sensor_here.upper() == sensor_name.upper():
                return sensor_name
            else:
                print(
                    f'[WARNING] Extract sensor {sensor_here} was not selected in the config file. Skipping extract...')
                return None

    def check_atm_correction(self, dataset, ac_name):
        ac_here = None
        if 'satellite_aco_processor' in dataset.ncattrs():
            ac_here = dataset.satellite_aco_processor
            if ac_here.startswith('Atmospheric Correction processor:'):
                if ac_here.upper().find(ac_name.upper()) > 0:
                    ac_here = ac_name.upper()
                else:
                    ac_here = 'STANDARD'
            if ac_here == 'Climate Change Initiative - European Space Agency':
                ac_here = 'CCI'

        if len(ac_here) == 0:
            ac_here = None
        if ac_here is None:
            print(
                f'[WARNING] Atmospheric correction set to {ac_name.upper()} despite of not being defined in the extract file')
            return ac_name.upper()
        else:
            if ac_here.upper() == ac_name.upper():
                return ac_name
            else:
                print(
                    f'[WARNING] Extract atmospheric processor {ac_here} was not selected in the config file. Skipping extract...')
                return None

    def check_platform(self, dataset, platform_name):
        platform_here = None
        if 'satellite' in dataset.ncattrs() and 'platform' in dataset.ncattrs():
            platform_here = f'{dataset.satellite}{dataset.platform}'
        if len(platform_here) == 0:
            platform_here = None
        if platform_here is None:
            print(f'[WARNING] Platform set to {platform_name} despite of not being defined in the extract file')
            return platform_name
        else:
            if platform_here.upper() == platform_name.upper():
                return platform_name
            else:
                print(
                    f'[WARNING] Extract platform {platform_here} was not selected in the config file {platform_name}. Skipping extract...')
                return None

    def check_site(self, fname, dataset, site_name):
        site_here = None
        if 'insitu_site_name' in dataset.ncattrs():
            site_here = dataset.insitu_site_name
        if site_here is None:
            if fname.upper().find(site_name.upper()) > 0:
                return site_name
            else:
                if site_name!='SHIPBORNE':
                    print(f'[WARNING] Site name set to {site_name} despite of not being defined in the extract file')
                return site_name
        else:
            if site_here.upper() == site_name.upper():
                return site_name
            else:
                print(f'[WARNING] Extract site {site_here} was not selected in the config file. Skipping extract...')
                return None

        # self.param_sat = {
        #     'satellite': sat_satellite.upper(),
        #     'sensor': sat_sensor.upper(),
        #     'platform': sat_platform.upper(),
        #     'resolution': sat_res.upper(),
        #     'ac': atm_corr.upper(),
        #     'prefix': prefix
        # }


class EXTRACT_LIST:

    def __init__(self, mdb_options,insitu_type, sat_type, verbose):
        self.mo = mdb_options
        self.insitu_type = insitu_type
        self.sat_type =sat_type
        self.verbose = verbose
        self.extract_files_by_date = None
        self.csv_files_by_date = None

    def prepare_extract_list(self,insituBase,allow_partial_mdb):
        self.csv_files_by_date = self.check_csv_extract_files(insituBase,allow_partial_mdb)
        if self.csv_files_by_date is None:
            pass ##self.extract_file_by_date should be prepare

    def get_valid_extracts_date(self,insituBase,date_here):
        mdb_options = self.mo.get_mdb_options()
        sat_options = self.mo.get_general_options('satellite_options')
        extract_path = mdb_options['extract_dir']
        if extract_path is None:
            return [None]*3
        time_diff_mu = mdb_options['time_diff_match_up']
        time_diff_mu = time_diff_mu * 60 ##from minutes to seconds
        ninsitu_max = mdb_options['ninsitu_max']
        time_diff_tv = mdb_options['time_diff_temporal_variability']
        insitu_time_day, insitu_lat_day, insitu_lon_day, insitu_indices_day =  [None]*4
        if time_diff_tv>0: ##if greater than zero, we need all the metadata for the day to complete the in situ metadata
                           ##with in situ data points outside the central pixel
            time_diff_tv = time_diff_tv * 60 ##from minutes to seconds
            insitu_time_day, insitu_lat_day, insitu_lon_day, insitu_indices_day = insituBase.get_metadata_date(date_here)

        info_extracts = {}
        time_extracts = []
        key = date_here.strftime('%Y%m%d')

        nall = 0


        if self.csv_files_by_date is not None and key in self.csv_files_by_date:
            file_csv = self.csv_files_by_date[key]

            df = pd.read_csv(file_csv,sep=';')
            extract_names_array = df['extract_file'][:]
            insitu_times_array = df['insitu_time'][:]
            insitu_times_array_ts = np.array([dt.strptime(x,'%Y-%m-%dT%H:%M:%S').replace(tzinfo=pytz.utc).timestamp() for x in insitu_times_array])
            extract_names = np.unique(extract_names_array)
            nall = len(extract_names)



            for extract_name in extract_names:

                file_extract = os.path.join(extract_path,extract_name)
                if os.path.exists(file_extract):

                    insitu_times_extract = insitu_times_array_ts[extract_names_array==extract_name]
                    insitu_indices_extract = df['insitu_index'][extract_names_array == extract_name].to_numpy()
                    insitu_lat_extract = df['insitu_lat'][extract_names_array == extract_name].to_numpy()
                    insitu_lon_extract = df['insitu_lon'][extract_names_array == extract_name].to_numpy()
                    info = self.check_extract_file(file_extract,insitu_indices_extract,insitu_times_extract,insitu_lat_extract,insitu_lon_extract,time_diff_mu,ninsitu_max)

                    if info is None:
                        return [None]*3
                    if len(info)==0:
                        continue

                    info = self.check_attributes(info, sat_options, insituBase)
                    if info is None:
                        return [None]*3

                    nvalid = len(info['insitu_indices'])

                    if time_diff_tv>0 and nvalid<ninsitu_max and nvalid<len(insitu_indices_day[0]):

                        info = self.check_insitu_variability_extract(info,insitu_time_day, insitu_lat_day, insitu_lon_day, insitu_indices_day[0],time_diff_tv,ninsitu_max)
                        if info is None:
                            return [None]*3

                    time_min_diff = dt.fromtimestamp(info['time_min_diff']).astimezone(pytz.utc)
                    ref = time_min_diff.strftime('%Y%m%dT%H%M%S')

                    info_extracts[ref] = info
                    time_extracts.append(time_min_diff)

        if len(time_extracts)>0:
            time_extracts.sort()



        return info_extracts,time_extracts,nall


    def check_insitu_variability_extract(self,info,insitu_time_day, insitu_lat_day, insitu_lon_day, insitu_indices_day,time_diff_tv,ninsitu_max):

        time_diff_prev = info['time_diff'][:]
        pos_min_time_diff_prev = np.argmin(time_diff_prev)
        index_min_time_diff = info['insitu_indices'][pos_min_time_diff_prev]
        pos_ref = int(np.where(insitu_indices_day==index_min_time_diff)[0][0])
        pos_min = int(pos_ref - np.floor(ninsitu_max / 2)) if (pos_ref - np.floor(ninsitu_max / 2)) > 0 else 0
        pos_max = pos_min + ninsitu_max
        if pos_max>=len(insitu_indices_day):
            pos_max = len(insitu_indices_day)
            pos_min = pos_max - ninsitu_max
            if pos_min<0:
                pos_min=0


        nvalid_new = pos_max-pos_min
        insitu_indices_new = insitu_indices_day[pos_min:pos_max]
        insitu_time_new = np.array([x.replace(tzinfo=pytz.utc).timestamp() for x in insitu_time_day[pos_min:pos_max]]).astype(np.float64)
        insitu_lat_new = insitu_lat_day[pos_min:pos_max]
        insitu_lon_new = insitu_lon_day[pos_min:pos_max]
        satellite_ts = info['satellite_time']

        time_diff_new = np.abs(satellite_ts-insitu_time_new)
        valid_new = time_diff_new<time_diff_tv

        print('-------')
        print(valid_new,np.sum(valid_new),nvalid_new)

        insitu_spatial_index_new = np.zeros(nvalid_new)
        nc_sat  = Dataset(info['file'])
        lat_array = np.squeeze(nc_sat.variables['satellite_latitude'][:])
        lon_array = np.squeeze(nc_sat.variables['satellite_longitude'][:])
        rc_center = int(np.floor(lat_array.shape[0]/2))
        nc_sat.close()

        for idx in range(nvalid_new):
            if not valid_new[idx]:
                insitu_spatial_index_new[idx] = -1
                continue
            r, c = cfs.find_row_column_from_lat_lon(lat_array.astype(np.float64), lon_array.astype(np.float64), insitu_lat_new[idx], insitu_lon_new[idx])
            print(idx,r,c)
            if np.isnan(r) and np.isnan(c):
                valid_new[idx] = False
                insitu_spatial_index_new[idx]=-1
            else:
                insitu_spatial_index_new[idx] = max(abs(r-rc_center),abs(c-rc_center))

        print(insitu_spatial_index_new)

        info['insitu_time'] = insitu_time_new[valid_new]
        info['insitu_lat'] = insitu_lat_new[valid_new]
        info['insitu_lon'] = insitu_lon_new[valid_new]
        info['insitu_indices'] = insitu_indices_new[valid_new]
        info['insitu_spatial_index'] = insitu_spatial_index_new[valid_new]
        info['time_diff'] = time_diff_new[valid_new]

        tf = info['time_diff'].copy()
        isi = info['insitu_spatial_index']


        tf[isi>0] = np.finfo(np.float32).max
        pos_min_tf = np.argmin(tf)

        time_min_diff = info['insitu_time'][pos_min_tf]

        if time_min_diff!=info['time_min_diff']:
            print(f'[INFO] Inconsistency in the in situ point with the minimal time difference.')
            return None

        return info

    def check_extract_file(self,file_extract,insitu_indices_extract,insitu_time_extract,insitu_lat_extract,insitu_lon_extract,time_diff_mu,ninsitu_max):
        ##check the extract file. insitu_times is the timestamp of measurements inside the central pixel of the extract
        info = {}
        if self.verbose:
            print(f'[INFO] -> Checking extract file using time difference with central pixel for extract {os.path.basename(file_extract)}')
        try:
            dataset = Dataset(file_extract)
        except:
            print(f'[WARNING] Extract {file_extract} is not a valid NetCDF file. Skipping...')
            return None
        info_atrs = {}
        atrs = ['site','satellite','platform','sensor','satellite_aco_processor','satellite_proc_version','res']
        for at in atrs:
            info_atrs[at] = dataset.getncattr(at) if at in dataset.ncattrs() else None

        satellite_time = float(dataset.variables['satellite_time'][0])
        time_diff= np.abs(satellite_time-insitu_time_extract)
        valid_ref = time_diff<time_diff_mu
        dataset.close()

        nvalid = np.count_nonzero(valid_ref)
        if nvalid==0:
            print(f'[WARNING] No insitu data points were found in the central pixel considering a time difference of {time_diff_mu/60:.0f} minutes')
            return info


        time_diff_valid = time_diff[valid_ref]
        insitu_time_valid = insitu_time_extract[valid_ref]
        insitu_indices_valid = insitu_indices_extract[valid_ref]
        insitu_lat_valid = insitu_lat_extract[valid_ref]
        insitu_lon_valid = insitu_lon_extract[valid_ref]

        min_ref = np.argmin(time_diff_valid)
        time_min_diff = insitu_time_valid[min_ref]

        if nvalid > ninsitu_max:
            if self.verbose:
                print(f'[INFO] Number of valid in situ data points ({nvalid}) is greater than the maximum number of in situ points {ninsitu_max}. Limiting to {ninsitu_max} points')
            pos_min = min_ref-np.floor(ninsitu_max/2) if (min_ref-np.floor(ninsitu_max/2))>0 else 0
            pos_max = pos_min + ninsitu_max
            if pos_max>len(insitu_time_valid):
                pos_max = len(insitu_time_valid)
                pos_min = pos_max - ninsitu_max
            pos_min = int(pos_min)
            pos_max = int(pos_max)
            insitu_time_valid = insitu_time_valid[pos_min:pos_max]
            time_diff_valid = time_diff_valid[pos_min:pos_max]
            insitu_indices_valid = insitu_indices_valid[pos_min:pos_max]
            insitu_lat_valid = insitu_lat_valid[pos_min:pos_max]
            insitu_lon_valid = insitu_lon_valid[pos_min:pos_max]
            nvalid = ninsitu_max


        info = {
            'file': file_extract,
            'insitu_indices': insitu_indices_valid,
            'insitu_time': insitu_time_valid,
            'insitu_lat': insitu_lat_valid,
            'insitu_lon': insitu_lon_valid,
            'insitu_spatial_index': np.zeros(nvalid),
            'satellite_time':satellite_time,
            'time_diff': time_diff_valid,
            'time_min_diff': time_min_diff,
            'site': info_atrs['site'],
            'satellite': info_atrs['satellite'],
            'platform': info_atrs['platform'],
            'sensor': info_atrs['sensor'],
            'proc_version': info_atrs['satellite_proc_version'],
            'aco_processor': info_atrs['satellite_aco_processor'],
            'res': info_atrs['res']
        }

        return info

    def check_attributes(self,info, sat_options, insituBase):
        if insituBase.fixed_site:
            expected_site = insituBase.site
            if expected_site!=info['site']:
                print(f'[ERROR] Attribute site in the extract {info["site"]} is different from the site option {expected_site}')
                return None
        else:
            info['site'] = 'SHIPBORNE' if insituBase.site is None else insituBase.site
        sat_attrs = ['satellite', 'platform', 'sensor', 'aco_processor', 'proc_version', 'res']
        for at in sat_attrs:
            if info[at] is None:
                info[at] = sat_options[at] if sat_options[at] is not None else ''
            else:
                if sat_options[at] is not None and sat_options[at].upper()!=info[at].upper():
                    print(f'[ERROR] Satellite attribute inconsistency. Found {info[at]} in the satellite extract file, but expected {sat_options[at]}')
                    return None
        return info




    def check_csv_extract_files(self,insituBase,allow_partial_mdb):

        extract_path = self.mo.get_extract_path()
        if extract_path is None:
            return None
        date_list = insituBase.date_list
        files_csv = {}
        nfiles_dates = 0
        for date_h in date_list:
            file_csv_here = os.path.join(extract_path, f'{insituBase.get_ref_date(date_h)}_extracts.csv')
            if os.path.exists(file_csv_here):
                files_csv[date_h.strftime('%Y%m%d')] = file_csv_here
                nfiles_dates = nfiles_dates+1
            # else:
            #     print(f'[WARNING] CSV extract file for date {date_h.strftime("%Y-%m-%d")} is not available')

        if nfiles_dates==0:
            print(f'[ERROR] CSV extract files are not available for the given data range.')
            return files_csv
        if nfiles_dates<len(date_list):
            print(f'[WARNING] CSV extract files were only retrieved for {nfiles_dates} of {len(date_list)} dates')
            if allow_partial_mdb:
                print(f'[WARNING] MDB build will continue using only {nfiles_dates} dates')
                return files_csv
            else:
                print(f'[ERROR] MDB builder was halted because CSV extract reference files were not available for all the potential dates.  Options:')
                print(f' --> To obtain a MDB file with only the available dates, run the same script with the option mdb_options/allow_partial_mdb: True')
                print(f' --> To create the CSV extract reference file for the missing dates, run the script with mdb_options/force_reference_csv: True')
                return None

        return files_csv