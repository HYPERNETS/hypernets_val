import numpy as np
from netCDF4 import Dataset
from datetime import datetime as dt
import numpy.ma as ma
#import pandas as pd
import argparse,os,subprocess,pytz,__init__,sys,shutil

from statsmodels.datasets.macrodata.data import variable_names

code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
import COMMON.common_functions as cfs
import MDB_builder.INSITU_base as ibase
from OPTIONS.OptionsManager import OptionsManager
from datetime import timedelta
from multiprocessing import Pool
from sat_extract import SatExtractOptions
import sat_extract as sextract

parser = argparse.ArgumentParser(description="Tool for downloading satellite sources starting from in situ data defined in the configurations file")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-c', "--config_file", help="Config File.", required=True)
parser.add_argument('-sd', "--startdate", help="The Start Date - format YYYY-MM-DD ")
parser.add_argument('-ed', "--enddate", help="The End Date - format YYYY-MM-DD ")
parser.add_argument('-sat_type', "--satellite_type", help="Satellite data type to be donwloaded",choices=['PACE'],default='PACE',required=True)
args = parser.parse_args()

def main():
    print('[INFO] Creating satellite extracts')
    if not args.config_file:
        return

    options = SatExtractOptions(args.config_file, args.verbose)
    if not options.is_valid:
        if not options.gmanager.is_valid():
            print('[ERROR] Problem processing global options')
        if not options.omanager.is_valid():
            print(f'[ERROR] Problem processing the configuration file {args.config_file}')
        return
    insitu_type, insitu_options = options.check_insitu_options()
    if insitu_type is None:
        return
    sat_type = args.satellite_type
    if args.verbose:
        print(f'[INFO] Satellite type: {sat_type}')
        print(f'[INFO] In sity type: {insitu_type}')
    insituBase = ibase.get_insitu_object(insitu_type,insitu_options,args.verbose)

    if insituBase is None:
        return
    input_path_info = options.get_input_path_info()
    if input_path_info is None:
        return

    now = dt.now().timestamp()
    file_metadata = os.path.join(input_path_info['path_source'],f'SourceDonwloadMetadata_{sat_type}_{insitu_type}_{now}.csv')
    if not insituBase.prepare_csv_metadata(file_metadata):
        return

    if sat_type=='PACE':
        download_pace_data(file_metadata,input_path_info)

    os.remove(file_metadata)

def download_pace_data(file_metadata,input_path_info):
    eistools_folder = os.path.join(os.path.dirname(code_home),'eistools')
    if not os.path.isdir(eistools_folder):
        print(f'[ERROR] Donwload PACE source required the package eistools')
        print(f'[ERROR] {eistools_folder} is not avaiable')
        return

    sys.path.append(eistools_folder)
    try:
        from nasa_download import NASA_DOWNLOAD
        ndownload = NASA_DOWNLOAD()
    except Exception as ex:
        print(f'[ERROR] NASA_DOWNLOAD class could not be loaded: {ex}')
        return


    fr = open(file_metadata)
    list_aop = []
    list_dates = []
    for line in fr:
        line_s = line.strip().split(',')
        date_here = dt.strptime(line_s[0],"%Y-%m-%d")
        region = [float(line_s[1]),float(line_s[2]),float(line_s[3]),float(line_s[4])]
        if args.verbose:
            print(f'[INFO] Getting granule list for {line_s[0]} in the area {region}')
        list_aop_here = ndownload.get_list_date('PACE_AOP',None, region,None,None, date_here, False)
        for l in list_aop_here:
            if l not in list_aop:
                list_aop.append(l)
                list_dates.append(date_here)

    ngranules = len(list_aop)
    if ngranules>0:
        for idx in range(ngranules):
            granule_aop = list_aop[idx]
            path_out = sextract.get_path_date( input_path_info['path_source'], input_path_info['org'],list_dates[idx],True)
            ndownload.download_granule(granule_aop,path_out,False)

            granule_geo = granule_aop.replace('OC_AOP', 'OC_GEO')
            file_geo = os.path.join(path_out, granule_geo)

            if not os.path.exists(file_geo):
                granule_l1b = granule_aop.replace('L2.OC_AOP','L1B')
                version = granule_aop.split('.')[4]
                granule_l1b = granule_l1b.replace(version,'V3')
                granule_l1b = granule_l1b.replace('.NRT','')
                ndownload.download_granule(granule_l1b, path_out, False)

                file_l1b = os.path.join(path_out,granule_l1b)
                if os.path.exists(file_l1b):
                    create_file_geo(file_l1b,file_geo)
                if os.path.exists(file_geo):
                    os.remove(file_l1b)



    fr.close()

def create_file_geo(file_l1b,file_geo):
    variable_names = []
    try:
        dataset_in = Dataset(file_l1b)
        geo_dataset = dataset_in['geolocation_data']

        for name in geo_dataset.variables:
            variable_names.append(name)

        nscans = len(dataset_in.dimensions['scans'])
        npixels =len(dataset_in.dimensions['pixels'])
    except Exception as ex:
        print(f'[ERROR] Error opening file level 1B: {file_l1b} Exception: {ex}')
        return
    if len(variable_names)>0 and nscans>0 and npixels>0:

        ncout = Dataset(file_geo,'w')
        ncout.createDimension('scans',nscans)
        ncout.createDimension('pixels',npixels)
        for name_var in variable_names:
            type = 'f4'
            fill_value = -999.0
            if name_var == 'height':
                type = 'i2'
                fill_value = -32767
            elif name_var == 'quality_flag' or name_var == 'watermask':
                type = 'u1'
                fill_value = 255

            var = ncout.createVariable(name_var,type,('scans','pixels'),fill_value=fill_value, zlib=True, complevel=6)
            var[:] = geo_dataset.variables[name_var][:]
            attributes = ['long_name','description','coordinates','units']
            for at in attributes:
                if at in geo_dataset.variables[name_var].ncattrs():
                    var.setncattr(at,geo_dataset.variables[name_var].getncattr(at))
            if 'valid_min' in geo_dataset.variables[name_var].ncattrs():
                valid_min = geo_dataset.variables[name_var].valid_min
                valid_max = geo_dataset.variables[name_var].valid_max
                if name_var.startswith('solar') or name_var.startswith('sensor'):
                    valid_min = float(valid_min * 0.01)
                    valid_max = float(valid_max * 0.01)
                var.valid_min = valid_min
                var.valid_max = valid_max
        ncout.close()

    dataset_in.close()


if __name__ == '__main__':
    main()
