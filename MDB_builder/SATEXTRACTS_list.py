import os
from netCDF4 import Dataset
from datetime import datetime as dt
import pytz


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
                print(f'[WARNING] Resulution set to {res_name} despite of not being defined in the extract file')
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
