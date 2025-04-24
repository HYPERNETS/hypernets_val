import argparse,os,configparser,sys,sat_extract
import pandas as pd
import numpy as np
import sys
from datetime import datetime as dt
from datetime import timedelta
from netCDF4 import Dataset
import pytz


code_home = os.path.dirname(os.path.dirname(sat_extract.__file__))
import COMMON.common_functions as cfs

parser = argparse.ArgumentParser(description="Create a spectra list in csv format")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-c', "--config_file", help="Config File.",required=True)
parser.add_argument('-o', "--output_file", help="Output CSV File.",required=True)
parser.add_argument('-wce',"--wce_expression",help="Wild card expression to limit input metadata files",default=None)
parser.add_argument('-stm',"--single_to_multiple",help="Mode to create multiple csv metadata file (one for day) starting from a single csv. Output should be a directory",action="store_true")

args = parser.parse_args()

def run_multiple_csv(options,output_file):
    path_csv = options['MULTIPLE_CSV_SELECTION']['path_csv']
    if not os.path.isdir(path_csv):
        print(f'[ERROR] Path to csv files {path_csv} was not found or is not a valid directory')
        return
    arg_wce = args.wce_expression


    col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = get_csv_options_from_file_config(
        options, 'MULTIPLE_CSV_SELECTION')
    extract_options = get_cmems_extract_options(options, 'MULTIPLE_CSV_SELECTION')
    path_source, org, wce, time_start, time_stop = get_find_product_info(options)

    if not os.path.isdir(path_source):
        try:
            os.mkdir(path_source)
        except:
            print(f'[ERROR] Satellite source path {path_source} does not exist and could not be created')
            return

    cmems_download_options = get_cmems_download_options(options)
    line_list = {}
    first_line = None
    all_rrs_data = {}

    for name in os.listdir(path_csv):
        if not name.endswith('csv'):
            continue

        if arg_wce is not None and not name.find(arg_wce)>=0:
            continue
        namefile = name[:-4]

        try:
            file_csv = os.path.join(path_csv, name)
            df = pd.read_csv(file_csv, sep=col_sep)
        except:
            print(f'[ERROR] File {path_csv} is not a valid csv separated by {col_sep}')
            return

        if args.verbose:
            print(f'[INFO] --------------------------------------------------------------------------------------')
            print(f'[INFO] Working with csv {namefile} Obtaining product(s)...')

        date_array_ts = df[col_date]
        ndata = len(date_array_ts)
        date_array = [None]*ndata
        date_array_str = [None]*ndata
        lat_array = np.array(df[col_lat])
        lon_array = np.array(df[col_lon])

        only_date_array_unique = []
        for idx in range(ndata):
            try:
                date_here = dt.strptime(date_array_ts[idx], format_date)
                date_str_here = date_here.strftime('%Y-%m-%d')
                if date_str_here not in only_date_array_unique:
                    only_date_array_unique.append(date_str_here)
                date_array_str[idx] = date_str_here
                date_array[idx] = date_here
            except:
                continue

        product_list = {}

        lat_map = None
        lon_map = None
        for date_str in only_date_array_unique:
            if args.verbose:
                print(f'[INFO] Checking available products for day {date_str} with single file option set to {extract_options["use_single_file"]}')
            list_files = None
            if extract_options['use_single_file']:
                fproduct = get_cmems_product_day_strict(path_source, org, dt.strptime(date_str, '%Y-%m-%d'),extract_options['dataset_name_file'],extract_options['dataset_name_format_date'],cmems_download_options)
            else:
                list_files = get_cmems_multiple_product_day(path_source, org, dt.strptime(date_str, '%Y-%m-%d'),extract_options['dataset_name_file'],extract_options['dataset_name_format_date'],extract_options['dataset_var_list'])
                fproduct = list_files[0] if list_files is not None else None

            if fproduct is None:
                print(f'[WARNING] No product(s) is(are) available for date: {date_str}')
                continue
            product_list[date_str]={
                'fproduct': fproduct,
                'list_files': list_files,
                'key_list': [],
                'row_list': [],
                'col_list': [],
                'rrs_data': None
            }
            if lat_map is None or lon_map is None:
                lat_map, lon_map = get_lat_lon_arrays(options, fproduct)

        if len(product_list)==0:
            print(f'[WARNING] No satellite products were available for the csv file {name}. Skipping...')
            continue
        elif args.verbose:
            print(f'[INFO] {len(product_list)} were obtained for csv file {name}')

        if args.verbose:
            print(f'[INFO] Getting row/col for each lat/long location to idenfity unique sites...')
        for idx in range(ndata):
            if args.verbose and (idx%100)==0:
                print(f'[INFO] Row: {idx}')
            datehere_str = date_array_str[idx]
            datehere  = date_array[idx]
            if datehere_str not in product_list.keys():
                continue
            lathere = lat_array[idx]
            lonhere = lon_array[idx]
            if np.isnan(lathere) or np.isnan(lonhere):
                continue
            fproduct = product_list[datehere_str]['fproduct']
            rc = get_rc(options, fproduct, lathere, lonhere, lat_map, lon_map)
            if rc is None:
                print(f'[WARNING] In situ location out of the limits of the satellite product. Skipping...')
                continue
            rint = int(rc[0])
            cint = int(rc[1])
            key = f'{datehere_str}_{rint}_{cint}'
            if key in line_list.keys():
                line = line_list[key]
                if datehere<line[3]:
                    line[3] = datehere
                if datehere>line[4]:
                    line[4] = datehere
                line[5] = line[5]+1
                continue

            if first_line is None:
                if extract_options['use_single_file']:
                    rrs_var_list = extract_options['rrs_var_list']
                else:
                    rrs_list = extract_options['rrs_list']
                    rrs_var_list = [f'RRS{get_wls(wl)}' for wl in rrs_list]
                first_line = 'Date;SatLat;SatLong;InSituFirst;InSituLast;NInSitu;Ref;' + ";".join(rrs_var_list)

            sat_lat, sat_long = get_lat_long_values(lat_map,lon_map,rint,cint)
            key_date_list = product_list[datehere_str]['key_list']
            if key not in key_date_list:
                product_list[datehere_str]['key_list'].append(key)
                product_list[datehere_str]['row_list'].append(rint)
                product_list[datehere_str]['col_list'].append(cint)
            # if extract_options['use_single_file']:

            line_list[key]= [datehere_str,f'{sat_lat}',f'{sat_long}',datehere,datehere,1,f'{key}']#+[f'{x}' for x in rrs_data]



        for date_str in product_list:
            print(f'[INFO] Extracting {len(product_list[date_str]["key_list"])} sites (pixels) for date {date_str}...')
            if extract_options['use_single_file']:
                rrs_var_list = extract_options['rrs_var_list']
                rrs_data = get_spectral_data_from_product(fproduct, rrs_var_list, product_list[date_str]['row_list'],product_list[date_str]['col_list'])
            else:
                rrs_list = extract_options['rrs_list']
                list_files = product_list[datehere_str]['list_files']
                rrs_data = get_spectral_data_from_list_files(list_files, rrs_list, product_list[date_str]['row_list'],product_list[date_str]['col_list'])
            rrs_data = np.ma.filled(rrs_data,-999.0)
            key_date_list = product_list[datehere_str]['key_list']
            for idx in range(rrs_data.shape[0]):
                all_rrs_data[key_date_list[idx]]=rrs_data[idx,:]


    if first_line is None or len(line_list)==0:
        return

    fw = open(output_file, 'w')
    fw.write(first_line)
    for line in line_list:
        line_list_here = line_list[line]
        line_list_here[3] = line_list_here[3].strftime('%H:%M:%S')
        line_list_here[4] = line_list_here[4].strftime('%H:%M:%S')
        line_list_here[5] = f'{line_list_here[5]:.0f}'
        key_ref = line_list_here[6]
        rrs_data = [f'{x}' for x in all_rrs_data[key_ref]]
        line_list_here = line_list_here + rrs_data
        line_str = ";".join(line_list_here)
        fw.write('\n')
        fw.write(line_str)
    fw.close()


def run_single_to_multiple(options,output_dir):
    section = 'SINGLE_TO_MULTIPLE'
    if not options.has_option(section,'path_csv'):
        print(f'[ERROR] path_csv option in required is section {section} for {section} mode.')
        return
    path_csv = options[section]['path_csv']
    if not os.path.isfile(path_csv):
        print(f'[ERROR] Path to csv file {path_csv} was not found or is not a valid file')
        return

    if args.verbose:
        print(f'[INFO] Reading file: {path_csv}')
    col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = get_csv_options_from_file_config(
        options, section)
    format_date_out = '%Y-%m-%d %H:%M:%S.%f'
    df = pd.read_csv(path_csv,sep=col_sep)
    dateobj_array = df[col_date]
    lat_array = df[col_lat]
    lon_array = df[col_lon]
    date_in_array = np.empty(dateobj_array.shape).astype(str)
    date_out_array = np.empty(dateobj_array.shape).astype(str)
    for idx,x in enumerate(dateobj_array):
        dobj = dt.strptime(x,format_date)
        date_in_array[idx] = dobj.strftime('%Y-%m-%d')
        date_out_array[idx] = dobj.strftime(format_date_out)

    date_list = np.unique(date_in_array)
    for date in date_list:
        if args.verbose:
            print(f'[INFO] Preparing metadata file for date {date}')

        date_here = date_out_array[date_in_array==date]
        lat_here = lat_array[date_in_array==date]
        lon_here = lon_array[date_in_array==date]
        nhere = date_here.shape[0]
        df_out = pd.DataFrame(index=np.arange(nhere),columns=['timestamp','lat','lon'])
        df_out['timestamp'] = np.array(date_here)
        df_out['lat'] = np.array(lat_here)
        df_out['lon'] = np.array(lon_here)


        # if date=='2023-04-04':
        #     print(df_out)
        #     print(lat_here)
        #     # valid = np.logical_and(np.isnan(lat_here) == False, np.isnan(lon_here) == False)
        #     # nvalid = np.sum(valid)
        #     # valid_lat = pd.isna(lat_here)
        #     # print(valid_lat.loc[0])
        #     # print(lat_here.loc[0])
        #     # lat_valid = np.isnan(np.array(lat_here))
        #     # print(date, nhere, nvalid,type(lat_here))
        #     # print(lat_valid,lat_here[0])
        #     break
        file_out = os.path.join(output_dir,f'{date}_metadata.csv')
        df_out.to_csv(file_out,sep=",",index=False)
        del df_out


def get_cmems_product_day_strict(path_source, org, datehere, dataset_name_file, dataset_name_format_date,
                                 cmems_download_options):
    product_day = get_cmems_product_day(path_source, org, datehere, dataset_name_file, dataset_name_format_date,
                                        cmems_download_options, False)

    if product_day is None:
        product_day = get_cmems_product_day(path_source, org, datehere, dataset_name_file, dataset_name_format_date,
                                            cmems_download_options, True)

    return product_day


def get_cmems_product_day(path_source, org, datehere, dataset_name_file, dataset_name_format_date,
                          cmems_download_options, use_myint):
    path_day = path_source
    if org is not None:
        if org == 'YYYYjjj':
            yearstr = datehere.strftime('%Y')
            jjjstr = datehere.strftime('%j')
            path_year = os.path.join(path_source, yearstr)
            path_day = os.path.join(path_source, yearstr, jjjstr)

    formats_date = dataset_name_format_date.split(',')

    if len(formats_date) == 1:
        datefile = datehere.strftime(dataset_name_format_date)
        namefile = dataset_name_file.replace('$DATE$', datefile)
    elif len(formats_date) == 2:
        datefile = datehere.strftime(formats_date[0].strip())
        namefile = dataset_name_file.replace('$DATE1$', datefile)
        datefile2 = datehere.strftime(formats_date[1].strip())
        namefile = namefile.replace('$DATE2$', datefile2)
    if use_myint:
        namefile = namefile.replace('my', 'myint')

    file = os.path.join(path_day, f'{namefile}')

    if not os.path.exists(file):
        ods = None
        if org == 'YYYYjjj':
            ods = '%Y/%j'
            path_year = os.path.join(path_source, yearstr)
            if not os.path.isdir(path_year):
                os.mkdir(path_year)
            if not os.path.isdir(path_day):
                os.mkdir(path_day)


        ##DOWNLOAD CMEMS SOURCES
        if cmems_download_options is not None:
            code_eistools = os.path.join(os.path.dirname(code_home), 'eistools')
            if os.path.exists(code_eistools):
                sys.path.append(code_eistools)
                from cmems_lois import CMEMS_LOIS
                clois = CMEMS_LOIS(args.verbose)
                cmems_download_options['start_date'] = datehere
                cmems_download_options['end_date'] = datehere
                cmems_download_options['date_list'] = None
                cmems_download_options['remote_name'] = namefile
                clois.make_cmems_download(cmems_download_options, True, path_source, ods, True)
            else:
                print(f'[WARNING] Code {code_eistools} for downloading is not available')

    if os.path.exists(file) and os.stat(file).st_size == 0:
        os.remove(file)
        if org == 'YYYYjjj':
            if len(os.listdir(path_day)) == 0:
                os.rmdir(path_day)
            if len(os.listdir(path_year)) == 0:
                os.rmdir(path_year)

    if not os.path.exists(file):
        print(f'[WARNING] File: {file} does not exist for date: {datefile}. Skiping...')
        return None

    return file

def main():
    print('[INFO] Creating spectra list')

    if not os.path.exists(args.config_file):
        print(f'[ERROR] File {args.config_file} does not exist')
        return
    options = configparser.ConfigParser()
    options.read(args.config_file)

    if args.single_to_multiple:
        output_dir = args.output_file
    else:
        output_file = args.output_file
        if not output_file.endswith('.csv'):
            print(f'[ERROR] Output file {output_file} should be a csv file')
            return
        output_dir = os.path.dirname(output_file)
    if not os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir)
        except:
            print(f'[ERROR] {output_dir} does not exist and could not be created. Please review permissions.')
            return
    if not os.access(output_dir,os.W_OK):
        print(f'[ERROR] You do not have write permissions in {output_dir} . Please review.')
        return

    if args.single_to_multiple and options.has_section('SINGLE_TO_MULTIPLE'):
        run_single_to_multiple(options,output_dir)

    if options.has_section('MULTIPLE_CSV_SELECTION') and options.has_option('MULTIPLE_CSV_SELECTION', 'path_csv'):
        run_multiple_csv(options,output_file)




def get_csv_options_from_file_config(options, section):
    if section == 'CSV_SELECTION':
        col_date = 'date'
        col_lat = 'lat'
        col_lon = 'lon'
        col_sep = ';'
        format_date = '%Y-%m-%dT%H:%M'
    elif section == 'MULTIPLE_CSV_SELECTION':
        col_date = 'timestamp'
        col_lat = 'lat'
        col_lon = 'lon'
        format_date = '%Y-%m-%d %H:%M:%S.%f'
        col_sep = ','
    else:
        col_date = 'dt'
        col_lat = 'lat'
        col_lon = 'lon'
        format_date = '%d-%b-%Y %H:%M:%S'
        col_sep = ','

    if options.has_option(section, 'col_date'):
        col_date = options[section]['col_date']
    if options.has_option(section, 'col_lat'):
        col_lat = options[section]['col_lat']
    if options.has_option(section, 'col_lon'):
        col_lon = options[section]['col_lon']
    if options.has_option(section, 'format_date'):
        format_date = options[section]['format_date']
    if options.has_option(section, 'col_sep'):
        col_sep = options[section]['col_sep']

    col_time = None
    format_time = '%H:%M:%S'
    if options.has_option(section, 'col_time'):
        col_time = options[section]['col_time']
    if options.has_option(section, 'format_time'):
        format_time = options[section]['format_time']

    return col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep

def get_cmems_extract_options(options, section):
    use_single_file = True
    n_bands = 0
    if options.has_option(section, 'use_single_file'):  ##SINGLE FILE SELECTION
        usf = options[section]['use_single_file']
        if usf.strip().lower() == 'true' or usf.strip() == '1':
            use_single_file = True
        else:
            use_single_file = False

    dataset_name_file = None
    dataset_name_format_date = '%Y%m%d'
    if options.has_option(section,'dataset_name_file'):
        dataset_name_file = options[section]['dataset_name_file']
    if options.has_option(section,'dataset_name_format_date'):
        dataset_name_format_date = options[section]['dataset_name_format_date']

    dataset_var_list = None
    dataset_var_list_out = None
    if options.has_option(section, 'dataset_var_list'):
        s = options[section]['dataset_var_list']
        if s.strip() != '':
            dataset_var_list = [x.strip() for x in s.split(',')]
            dataset_var_list_out = dataset_var_list
            if options.has_option(section, 'dataset_var_list_out'):
                s = options[section]['dataset_var_list_out']
                dataset_var_list_out = [x.strip() for x in s.split(',')]

    rrs_list = []
    rrs_var_list = []
    if use_single_file:
        is_reflectance = False
        if options.has_option(section, 'var_rrs_list') and options.has_option(section, 'var_rrs_format'):

            var_rrs_list = options[section]['var_rrs_list'].strip()
            var_rrs_format = options[section]['var_rrs_format'].strip()
            is_reflectance = True
            for r in var_rrs_list.split(','):
                try:
                    rs = r.strip().replace('.', '_')
                    rrs_here = float(r.strip())
                    rrs_list.append(rrs_here)
                    var_here = var_rrs_format.replace('$BAND$', rs)
                    rrs_var_list.append(var_here)
                except:
                    is_reflectance = False
                    break
            if is_reflectance:
                n_bands = len(rrs_list)
    else:
        is_reflectance = True
        for var in dataset_var_list:
            try:
                value = float(var)
                rrs_list.append(value)
            except:
                is_reflectance = False
        if is_reflectance:
            n_bands = len(dataset_var_list)

    options = {
        'use_single_file': use_single_file,
        'is_reflectance': is_reflectance,
        'n_bands': n_bands,
        'dataset_name_file': dataset_name_file,
        'dataset_name_format_date': dataset_name_format_date,
        'dataset_var_list': dataset_var_list,
        'dataset_var_list_out': dataset_var_list_out,
        'rrs_list': rrs_list,
        'rrs_var_list': rrs_var_list
    }

    return options

def get_find_product_info(options):
    path_source = options['file_path']['sat_source_dir']
    org = None
    if options.has_option('file_path', 'sat_source_dir_organization') and options['file_path'][
        'sat_source_dir_organization']:
        org = options['file_path']['sat_source_dir_organization']
    wce = '*'
    if options.has_option('file_path', 'wce') and options['file_path']['wce']:
        wce = options['file_path']['wce']
    time_start = None
    time_stop = None
    section = 'Time_and_sites_selection'
    if options.has_section(section):
        if options.has_option(section, 'time_start'):
            time_start = dt.strptime(options['Time_and_sites_selection']['time_start'], '%Y-%m-%d')
        if options.has_option(section, 'time_stop'):
            time_stop = dt.strptime(options['Time_and_sites_selection']['time_stop'], '%Y-%m-%d') + timedelta(hours=24)
    # print('temp')
    return path_source, org, wce, time_start, time_stop

def get_cmems_download_options(options):
    section = 'satellite_options'
    options_download = ['download_product', 'download_dataset', 'download_endpoint', 'download_bucket', 'download_tag']
    keys_download = ['product', 'dataset', 'endpoint', 'bucket', 'tag']
    keys_present = [False, False, True, False, False]  ##endpoint is not required, always True
    cmems_donwload_options = {'start_date': None, 'end_date': None}
    for idx, op in enumerate(options_download):
        if options.has_option(section, op):
            cmems_donwload_options[keys_download[idx]] = options[section][op]
            keys_present[idx] = True
        else:
            cmems_donwload_options[keys_download[idx]] = None

    if keys_present.count(True) == 5:
        return cmems_donwload_options
    else:
        return None

def get_lat_lon_arrays(options, file_nc):
    nc_sat = Dataset(file_nc, 'r')
    var_lat, var_lon = get_lat_long_var_names(options)
    lat, lon = get_lat_long_arrays(nc_sat, var_lat, var_lon)
    nc_sat.close()
    return lat, lon

def get_lat_long_var_names(options):
    var_lat = 'lat'
    var_lon = 'lon'
    if options.has_option('satellite_options', 'lat_variable'):
        var_lat = options['satellite_options']['lat_variable']
    if options.has_option('satellite_options', 'lon_variable'):
        var_lon = options['satellite_options']['lon_variable']
    return var_lat, var_lon

def get_lat_long_arrays(nc_sat, var_lat, var_lon):
    vlat = nc_sat.variables[var_lat]
    vlon = nc_sat.variables[var_lon]
    # if vlat.ndim == 2 and vlon.ndim == 2:
    #     lat = nc_sat.variables[var_lat][:, :]
    #     lon = nc_sat.variables[var_lon][:, :]
    # if vlat.ndim == 1 and vlon.ndim == 1:
    lat = nc_sat.variables[var_lat][:]
    lon = nc_sat.variables[var_lon][:]
    return lat, lon

def get_lat_long_values(lat_array,lon_array,row,col):
    if len(lat_array.shape)==1 and len(lon_array.shape)==1:
        lat_v = lat_array[row]
        lon_v = lon_array[col]
    else:
        lat_v = lat_array[row,col]
        lon_v = lon_array[row,col]
    return lat_v,lon_v

def get_spectral_data_from_product(fproduct,rrs_var_list,rint,cint):
    if np.isscalar(rint) and np.isscalar(cint):
        rint = [rint]
        cint = [cint]
    ndata = len(rint)
    nbands = len(rrs_var_list)
    data = np.zeros((ndata,nbands))
    nc_sat = Dataset(fproduct)
    for iband in range(nbands):
        var =  rrs_var_list[iband]
        array = np.ma.squeeze(nc_sat.variables[var][:])
        vals = array[rint,cint]
        data[:, iband] = vals
    nc_sat.close()
    return data

def get_spectral_data_from_list_files(list_files,rrs_list,rint,cint):
    if np.isscalar(rint) and np.isscalar(cint):
        rint = [rint]
        cint = [cint]
    ndata = len(rint)
    nbands = len(list_files)
    data = np.zeros((ndata,nbands))


    for iband in range(nbands):
        file_in = list_files[iband]
        wl = rrs_list[iband]
        nc_sat = Dataset(file_in)
        wls = get_wls(wl)
        var = f'RRS{wls}'
        if not var in nc_sat.variables:
            print('aqui no llega')
            var = None
            for name in nc_sat.variables:
                if name.find(wls)>0:
                    var = name
                    break


        array = np.ma.squeeze(nc_sat.variables[var][:])
        vals = array[rint,cint]
        data[:,iband] = vals

        nc_sat.close()
    return data

def get_wls(wl):
    wls = f'{wl:.2f}'
    wls = wls.replace('.', '_')
    if wls.endswith('_00'):
        wls = wls[:-3]
    if wls.find('_') > 0 and wls.endswith('0'):
        wls = wls[:-1]
    return wls

def get_cmems_multiple_product_day(path_source, org, datehere, dataset_name_file, dataset_name_format_date,
                                   dataset_var_list):
    path_day = path_source
    if org is not None:
        if org == 'YYYYjjj':
            yearstr = datehere.strftime('%Y')
            jjjstr = datehere.strftime('%j')
            path_day = os.path.join(path_source, yearstr, jjjstr)
    strdate = datehere.strftime(dataset_name_format_date)
    ncfiles = []

    for var in dataset_var_list:
        name = dataset_name_file
        var = var.replace('.', '_')
        name = name.replace('$DATE$', strdate)
        name = name.replace('$BAND$', var)
        fname = os.path.join(path_day, name)
        if os.path.exists(fname):
            ncfiles.append(fname)
        else:
            print(f'[WARNING] File: {fname} was not found')

    if len(ncfiles) == len(dataset_var_list):
        return ncfiles
    else:
        return None

def get_info_from_row(row, col_date, col_time, format_date, format_time, col_lat, col_lon):
    try:
        datehere = None
        if col_time is not None:
            try:
                datetimerow = f'{row[col_date].strip()}T{row[col_time].strip()}'
                if np.isreal(datetimerow):
                    datetimerow = f'{datetimerow:.0f}'
                format_datetime = f'{format_date}T{format_time}'
                datehere = dt.strptime(datetimerow, format_datetime).replace(tzinfo=pytz.utc)
            except:
                pass
        if datehere is None:
            datetimerow = row[col_date]
            if np.isreal(datetimerow):
                datetimerow = f'{datetimerow:.0f}'
            format_datetime = format_date

            datehere = dt.strptime(datetimerow, format_datetime).replace(tzinfo=pytz.utc)
            #datehere = datehere.replace(hour=12, minute=0, second=0,microsecond=0).replace(tzinfo=pytz.utc)
        lathere = float(row[col_lat])
        lonhere = float(row[col_lon])
        if np.isnan(lathere):
            return [datehere,None,lonhere]
        if np.isnan(lonhere):
            return [datehere,lonhere,None]

        return datehere, lathere, lonhere
    except:
        return [None]*3

def get_rc(options, file_nc, insitu_lat, insitu_lon, lat, lon):
    if lat is None or lon is None:
        lat, lon = get_lat_lon_arrays(options, file_nc)
    rc = None
    if cfs.contain_location(lat, lon, insitu_lat, insitu_lon) == 1:
        if lat.ndim == 1 and lon.ndim == 1:
            r = np.argmin(np.abs(lat - insitu_lat))
            c = np.argmin(np.abs(lon - insitu_lon))
        else:
            r, c = cfs.find_row_column_from_lat_lon(lat, lon, insitu_lat, insitu_lon)

        if (0 <= r < lat.shape[0]) and (0 <= c < lon.shape[0]):
            rc = [r, c]

    return  rc

if __name__ == '__main__':
    main()