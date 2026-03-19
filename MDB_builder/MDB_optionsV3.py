import os,sys,__init__
code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
from OPTIONS.OptionsManager import OptionsManager



class MDBBuilderOptions:

    def __init__(self,config_file,verbose):
        self.verbose = verbose
        self.omanager = None
        self.gmanager = None

        general_options_file = os.path.join(code_home,'OPTIONS','general_options.ini')
        self.gmanager = OptionsManager(general_options_file,None)
        self.omanager = OptionsManager(config_file, None)

        self.is_valid = self.gmanager.is_valid() and self.omanager.is_valid()

    def get_general_options(self,section):
        if not self.is_valid:
            return None
        soptions, required = self.gmanager.get_retrieve_options(section)

        if soptions is not None:
            options = self.omanager.get_options_as_dict(section, soptions,required)
            return options
        else:
            return None

    def get_insitu_options(self,insitu_type):
        if insitu_type is None:
            return None
        insitu_options_file = os.path.join(code_home, 'OPTIONS', 'insitu_options.ini')
        om = OptionsManager(insitu_options_file, None)
        if not om.is_valid():
            return None
        insitu_type_list = om.get_option('GLOBAL','type_list',None,None,'strlist')
        if insitu_type_list is None:
            print(f'[ERROR] In situ type list (option GLOBAL/type_list) in file {insitu_options_file} is unavailable or not valid')
            return None
        if insitu_type not in insitu_type_list:
            print(f'[ERROR] In situ type {insitu_type} is not defined. Please choose among: {insitu_type_list}')
            return None
        sopt, required = self.get_retrive_options_insitu(insitu_type)
        if sopt is None:
            print(f'[ERROR] Section {insitu_type} defining in situ type {insitu_type} is not available in {insitu_options_file}.')
            return None
        insitu_options = self.omanager.get_options_as_dict(insitu_type, sopt, required)
        if insitu_options is None:
            print(f'[ERROR] In situ options for in situ type {insitu_type} could not be retrieved')
            return None
        return insitu_options

    def get_satellite_options(self,sat_type):
        if sat_type is None:
            return None
        sat_options_file = os.path.join(code_home, 'OPTIONS', 'satellite_options.ini')
        om = OptionsManager(sat_options_file, None)
        if not om.is_valid():
            return None
        sat_type_list = om.get_option('GLOBAL','type_list',None,None,'strlist')
        if sat_type_list is None:
            print(f'[ERROR] Satellite type list (option GLOBAL/type_list) in file {sat_options_file} is unavailable or not valid')
            return None
        if sat_type not in sat_type_list:
            print(f'[ERROR] Satellite type {sat_type} is not defined. Please choose among: {sat_type_list}')
            return None
        sopt, required = self.get_retrive_options_sat(sat_type)
        if sopt is None:
            print(f'[ERROR] Section {sat_type} defining in satellite type {sat_type} is not available in {sat_options_file}.')
            return None
        sat_options = self.omanager.get_options_as_dict(sat_type, sopt, required)
        if sat_options is None:
            print(f'[ERROR] Satellite options for satellite type {sat_type} could not be retrieved')
            return None
        return sat_options

    def get_retrive_options_insitu(self,section):
        insitu_options_file = os.path.join(code_home, 'OPTIONS', 'insitu_options.ini')
        om = OptionsManager(insitu_options_file, None)
        if om.is_valid():
            return om.get_retrieve_options(section)
        else:
            return None

    def get_retrive_options_sat(self,section):
        insitu_options_file = os.path.join(code_home, 'OPTIONS', 'satellite_options.ini')
        om = OptionsManager(insitu_options_file, None)
        if om.is_valid():
            return om.get_retrieve_options(section)
        else:
            return None

    def get_mdb_options(self):
        return self.get_general_options('mdb_options')

    def get_file_date_list_check(self):
        options = self.get_mdb_options()
        return options['file_date_list_check'] if options is not None else None

    def get_extract_path(self):
        options = self.get_mdb_options()
        return options['extract_dir'] if options is not None else None

    def allow_partial_mdb(self):
        options = self.get_mdb_options()
        return options['allow_partial_mdb']

    def overwrite(self):
        options = self.get_general_options('file_path')
        return options['overwrite']

    def overwrite_mini(self):
        options = self.get_mdb_options()
        return options['overwrite_mini']

    def delete_mini(self):
        options = self.get_mdb_options()
        return options['delete_mini']

    def get_mdb_name(self,insitu_type,start_date,end_date):
        insitu_options = self.get_insitu_options(insitu_type)
        site = insitu_options['site'].replace(' ', '_')
        info = self.get_general_options('satellite_options')
        start = start_date.strftime("%Y%m%d")
        end = end_date.strftime("%Y%m%d")
        platform = '' if info['platform'] is None else info['platform']
        name = f'MDB_{info["satellite"]}{platform}_{info["sensor"]}_{site}_{start}_{end}.nc'
        return name
