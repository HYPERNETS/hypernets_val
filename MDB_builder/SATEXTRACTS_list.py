import os
from netCDF4 import Dataset
from datetime import datetime as dt

class SAT_EXTRACTS_LIST:
    def __init__(self, boptions, verbose):
        self.boptions = boptions
        self.verbose = verbose

    def get_list_as_dict(self):
        start_date = self.boptions.start_date
        end_date = self.boptions.end_date
        site = self.boptions.param_insitu


        prefix = self.boptions.param_sat['prefix']
        sat_extract_dir = self.boptions.satellite_path_source

        sat_list = {}

        for name in os.listdir(sat_extract_dir):
            if not name.endswith('.nc'):
                continue
            if prefix is not None and not name.startswith(prefix):
                continue


            fextract = os.path.join(sat_extract_dir, name)
            if self.verbose:
                print(f'[INFO] Checking extract file: {name}')
            dataset = Dataset(fextract)
            time_here = self.check_time(dataset,start_date,end_date)
            if time_here is None:
                dataset.close()
                continue
            site_here = self.check_site(name, dataset, site['station_name'])
            if site_here is None:
                dataset.close()
                continue
            platform = self.check_platform(dataset,self.boptions.param_sat['sat_platform'])
            if platform is None:
                dataset.close()
                continue
            ac = self.check_atm_correction(dataset,self.boptions.param_sat['ac'])
            if ac is None:
                dataset.close()
                continue
            sensor = self.check_sensor(dataset,self.boptions['sensor'])
            if sensor is None:
                dataset.close()
                continue
            resolution = self.check_resolution(name,dataset,self.boptions['resolution'])
            if resolution is None:
                dataset.close()
                continue

            sat_list[name] = {
                'path':fextract,
                'time': time_here.strftime('%Y%m%dT%H%M%S'),
                'site': site_here,
                'sensor':sensor,
                'platform':platform,
                'ac':ac,
                'resolution':resolution
            }

        return sat_list

    def check_time(self,dataset,start_date,end_date):
        datetime_here = dt.fromtimestamp(float(dataset.variables['satellite_time'][0]))
        if start_date <= datetime_here <= end_date:
            return datetime_here
        else:
            return None

    def check_resolution(self,fname,dataset,res_name):
        res_here = None
        if 'resolution' in dataset.ncattrs:
            res_here = dataset.resolution
        if res_here is None:
            if fname.find(res_name)>0:
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

    def check_sensor(self,dataset,sensor_name):
        sensor_here = None
        if 'sensor' in dataset.ncattrs:
            sensor_here = dataset.sensor
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

    def check_atm_correction(self,dataset,ac_name):
        ac_here = None
        if 'satellite_aco_processor' in dataset.ncattrs():
            ac_here = dataset.satellite_aco_processor
            if ac_here.startswith('Atmospheric Correction'):
                ac_here = 'STANDARD'
        if ac_here is None:
            print(f'[WARNING] Atmospheric correction set to {ac_name} despite of not being defined in the extract file')
            return ac_name
        else:
            if ac_here.upper() == ac_name.upper():
                return ac_name
            else:
                print(
                    f'[WARNING] Extract atmospheric processor {ac_here} was not selected in the config file. Skipping extract...')
                return None

    def check_platform(self,dataset,platform_name):
        platform_here = None
        if 'satellite' in dataset.ncattrs() and 'platform' in dataset.ncattrs():
            platform_here = f'{dataset.satellite}{dataset.platform}'

        if platform_here is None:
            print(f'[WARNING] Platform set to {platform_name} despite of not being defined in the extract file')
            return platform_name
        else:
            if platform_here.upper()==platform_name.upper():
                return platform_name
            else:
                print(f'[WARNING] Extract platform {platform_here} was not selected in the config file. Skipping extract...')
                return None


    def check_site(self, fname, dataset, site_name):
        site_here = None
        if 'site' in dataset.ncattrs():
            site_here = dataset.site
        if site_here is None:
            if fname.upper().find(site_name.upper()) > 0:
                return site_name
            else:
                print(f'[WARNING] Site name set to {site_name} despite of not being defined in the extract file')
                return site_name
        else:
            if site_here.upper()==site_name.upper():
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
