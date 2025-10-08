
from netCDF4 import Dataset
from datetime import datetime as dt

import argparse,os,subprocess,pytz,__init__,sys,shutil

code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
import MDB_builder.INSITU_base as ibase
import COMMON.args_functions as arf
from OPTIONS.OptionsManager import OptionsManager
from sat_extract import SatExtractOptions
import sat_extract as sextract

parser = argparse.ArgumentParser(description="Tool for downloading satellite sources starting from in situ data defined in the configurations file")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-c', "--config_file", help="Config File.", required=True)
parser.add_argument('-sd', "--start_date", help="The Start Date - format YYYY-MM-DD ")
parser.add_argument('-ed', "--end_date", help="The End Date - format YYYY-MM-DD ")
parser.add_argument('-insitu',"--insitu_type",help="In situ type")
parser.add_argument('-sat', "--sat_type", help="Satellite data type to be donwloaded",choices=['PACE','OLCI','CMEMS','CCI'],default='PACE',required=True)
args = parser.parse_args()

class DownloadOptions:

    def __init__(self,config_file):
        download_options_file = os.path.join(code_home, 'OPTIONS', 'download_options.ini')
        self.do = OptionsManager(download_options_file, None)
        self.om = OptionsManager(config_file,None)
        self.is_valid = self.do.is_valid() and self.om.is_valid()

    def get_eumetsat_download_options(self):
        return self.get_options('eumetsat_download')

    def get_cmems_download_options(self):
        return self.get_options('cmems_download')

    def get_options(self, section):
        if not self.is_valid:
            return None
        soptions, required = self.do.get_retrieve_options(section)


        if soptions is not None:
            options = self.om.get_options_as_dict(section, soptions, required)
            return options
        else:
            return None




def main():
    print('[INFO] Source download for extracts preparation')
    if not args.config_file:
        return

    options = SatExtractOptions(args.config_file, args.verbose)
    if not options.is_valid:
        if not options.gmanager.is_valid():
            print('[ERROR] Problem processing global options')
        if not options.omanager.is_valid():
            print(f'[ERROR] Problem processing the configuration file {args.config_file}')
        return
    insitu_type, insitu_options = options.check_insitu_options(args.insitu_type)
    if insitu_type is None:
        return
    sat_type = args.sat_type
    if args.verbose:
        print(f'[INFO] Satellite type: {sat_type}')
        print(f'[INFO] In situ type: {insitu_type}')
    insituBase = ibase.get_insitu_object(insitu_type,insitu_options,args.verbose)

    if insituBase is None:
        return
    input_path_info = options.get_input_path_info()
    if input_path_info is None:
        return


    now = dt.now().timestamp()
    file_metadata = os.path.join(input_path_info['path_source'],f'SourceDonwloadMetadata_{sat_type}_{insitu_type}_{now}.csv')
    insituBase.start_date, insituBase.end_date = arf.get_start_end_date_from_args(args)

    if not insituBase.prepare_csv_metadata(file_metadata):
        return

    doptions = DownloadOptions(args.config_file)


    if sat_type=='PACE':
        download_pace_data(file_metadata,input_path_info)

    if sat_type=='OLCI':
        eumetsat_options = doptions.get_eumetsat_download_options()
        download_olci_data(file_metadata,input_path_info,eumetsat_options)

    if sat_type=='CMEMS':
        cmems_options = doptions.get_cmems_download_options()
        if cmems_options is not None:
            download_cmems_data(file_metadata,input_path_info,cmems_options)

    if sat_type=='CCI':
        download_cci_data(file_metadata,input_path_info)




    #os.remove(file_metadata)

def download_cci_data(file_metadata,input_path_info):
    eistools_folder = os.path.join(os.path.dirname(code_home), 'eistools')
    if not os.path.isdir(eistools_folder):
        print(f'[ERROR] Donwload OLCI sources requires the package eistools')
        print(f'[ERROR] {eistools_folder} is not avaiable')
        return
    sys.path.append(eistools_folder)
    try:
        from html_download import OC_CCI_V6_Download
        cciDownload = OC_CCI_V6_Download()
    except Exception as ex:
        print(f'[ERROR] cciDownload class could not be loaded: {ex}')
        return
    fr = open(file_metadata)
    already_download = []  ##refs to files already download, equal to date_utm or only data
    for line in fr:
        line_s = line.strip().split(',')
        datehere = dt.strptime(line_s[0], "%Y-%m-%d")
        ref = f'{datehere.strftime("%Y%m%d")}'
        if ref not in already_download:
            status = cciDownload.download_date(datehere,input_path_info['path_source'])
            if status==0 and cciDownload.check_file_date(input_path_info['path_source'],datehere):
                already_download.append(ref)
            else:
                print(f'[ERROR] Error downloading OC-CCI file for date: {ref}')

def download_cmems_data(file_metadata,input_path_info,options):
    eistools_folder = os.path.join(os.path.dirname(code_home), 'eistools')
    if not os.path.isdir(eistools_folder):
        print(f'[ERROR] Donwload OLCI sources requires the package eistools')
        print(f'[ERROR] {eistools_folder} is not avaiable')
        return
    sys.path.append(eistools_folder)
    try:
        from cmems_lois import CMEMS_LOIS
        cdownload = CMEMS_LOIS(args.verbose)
    except Exception as ex:
        print(f'[ERROR] EUMDAC_LOIS class could not be loaded: {ex}')
        return

    dataset_name_file = options['dataset_name_file']
    if dataset_name_file is None:
        dataset_name_file = f'$DATE$_{options["dataset"]}.nc'
    if args.verbose:
        print(f'[INFO] Dataset name file: {dataset_name_file}')


    fr = open(file_metadata)
    already_download = [] ##refs to files already download, equal to date_utm or only data
    for line in fr:
        line_s = line.strip().split(',')
        datehere = dt.strptime(line_s[0], "%Y-%m-%d")
        datefile = datehere.strftime(options['dataset_name_format_date'])
        namefile = dataset_name_file.replace('$DATE$', datefile)
        options['start_date'] = datehere
        options['end_date'] = datehere
        options['date_list'] = None

        if namefile.find('$UTM$')>0:
            namefile_base = namefile
            region = [float(line_s[1]), float(line_s[2]), float(line_s[3]), float(line_s[4])]
            lat_points = [region[0],region[0], region[1],region[1]]
            lon_points = [region[2], region[3], region[2], region[3]]
            for lat_p,lon_p in zip(lat_points,lon_points):
                utm_zone = sextract.get_utm_zone(lat_p,lon_p)
                ref = f'{datehere.strftime("%Y%m%d")}_{utm_zone}'
                if ref not in already_download:
                    namefile = namefile_base.replace('$UTM$',utm_zone)
                    options['remote_name'] = namefile
                    launch_download_cmems(cdownload, options, input_path_info, datefile, utm_zone)
                    already_download.append(ref)
        else:
            ref = f'{datehere.strftime("%Y%m%d")}'
            if ref not in already_download:
                options['remote_name'] = namefile
                launch_download_cmems(cdownload,options,input_path_info, datefile,None)
                already_download.append(ref)



def launch_download_cmems(cdownload,options,input_path_info,datefile,utm_zone):
    cdownload.make_cmems_download(options, True, input_path_info['path_source'], input_path_info['org'], False)
    if options['extra_dataset'] is not None:


        for index in range(len(options['extra_dataset'])):
            options_extra = options.copy()
            edataset = options['extra_dataset'][index]
            if options['extra_dataset_name_file'] is not None:
                dataset_name_file_e = options['extra_dataset_name_file'][index]
            else:
                dataset_name_file_e = f'$DATE$_{edataset}.nc'

            if options['extra_bucket'] is not None:
                ebucket = options['extra_bucket'][index].strip()
                if ebucket.strip()!='*':
                    options_extra['bucket'] = ebucket

            if options['extra_tag'] is not None:
                etag = options['extra_tag'][index].strip()
                if etag!='*':
                    options_extra['tag'] = etag
            #print(options_extra)
            namefile_e = dataset_name_file_e.replace('$DATE$', datefile)
            if utm_zone is not None:
                namefile_e = namefile_e.replace('$UTM$', utm_zone)
            options_extra['dataset'] = edataset
            options_extra['remote_name'] = namefile_e

            cdownload.make_cmems_download(options_extra, True, input_path_info['path_source'], input_path_info['org'],False)

# def get_utm_zone(latitude,longitude):
#     #from pyproj import CRS
#     #longitude = -8.85
#     #latitude = 42.59
#     # Determine the UTM zone number
#     zone_number = int((longitude + 180) / 6) + 1
#     ZONE_LETTERS = "CDEFGHJKLMNPQRSTUVWXX"
#     zone_letter = ZONE_LETTERS[int(latitude + 80) >> 3]
#     # Determine the hemisphere
#     #hemisphere = 'north' if latitude >= 0 else 'south'
#     #utm_crs = CRS.from_dict({'proj': 'utm', 'zone': zone_number, 'south': hemisphere == 'south'})
#     return f'{zone_number}{zone_letter}'



def download_olci_data(file_metadata,input_path_info,options):
    eistools_folder = os.path.join(os.path.dirname(code_home), 'eistools')
    if not os.path.isdir(eistools_folder):
        print(f'[ERROR] Donwload OLCI sources requires the package eistools')
        print(f'[ERROR] {eistools_folder} is not available')
        return
    sys.path.append(eistools_folder)
    try:
        from eumdac_lois import EUMDAC_LOIS
        edownload = EUMDAC_LOIS(args.verbose,options['user'])
    except Exception as ex:
        print(f'[ERROR] EUMDAC_LOIS class could not be loaded: {ex}')
        return

    fr = open(file_metadata)
    list_granules = []
    list_dates = []
    list_products = []
    for line in fr:
        line_s = line.strip().split(',')
        date_here = dt.strptime(line_s[0], "%Y-%m-%d")
        region = [float(line_s[1]), float(line_s[2]), float(line_s[3]), float(line_s[4])]
        if args.verbose:
            print(f'[INFO] Getting granule list for {line_s[0]} in the area {region}')
        products, product_names, collection_id = edownload.search_olci_by_bbox(date_here, options['resolution'], options['level'],
                                                                          region, -1, -1,options['timeliness'])

        if products is None or len(products)==0:
            # fr.close()
            # os.remove(file_metadata)
            continue
        for product,name in zip(products,product_names):
            if name not in list_granules:
                list_granules.append(name)
                list_products.append(product)
                list_dates.append(date_here)

    ngranules = len(list_granules)
    if ngranules == 0:
        print(f'[WARNING] No granules for download were found. ')
        return

    for idx in range(ngranules):
        product = list_products[idx]
        path_out = sextract.get_path_date(input_path_info['path_source'], input_path_info['org'], list_dates[idx], True)
        if path_out is None:
            print(f'[ERROR] Download requires a valid output directory. Skipping...')
            continue
        edownload.download_product(product, path_out, False)

def download_pace_data(file_metadata,input_path_info):
    eistools_folder = os.path.join(os.path.dirname(code_home),'eistools')
    if not os.path.isdir(eistools_folder):
        print(f'[ERROR] Donwload PACE sources requires the package eistools')
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
        if list_aop_here is None:
            fr.close()
            os.remove(file_metadata)
            return
        for l in list_aop_here:
            if l not in list_aop:
                list_aop.append(l)
                list_dates.append(date_here)



    ngranules = len(list_aop)
    if ngranules==0:
        print(f'[WARNING] No granules for download were found. ')
        return


    for idx in range(ngranules):
        granule_aop = list_aop[idx]
        path_out = sextract.get_path_date(input_path_info['path_source'], input_path_info['org'],list_dates[idx],True)
        if path_out is None:
            print(f'[ERROR] Download requires a valid output directory. Skipping...')
            continue
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
        else:
            print(f'[INFO] Geo file {granule_geo} already exists. Skipping...')

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
