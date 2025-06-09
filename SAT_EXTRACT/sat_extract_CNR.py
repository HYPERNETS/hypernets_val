import argparse, os, sys
import shutil

import pandas as pd
from netCDF4 import Dataset
from datetime import datetime as dt
import numpy.ma as ma
import numpy as np
from multiprocessing import Pool


# import user defined functions from other .py
import sat_extract

code_home = os.path.dirname(os.path.dirname(sat_extract.__file__))
sys.path.append(code_home)
from SAT_EXTRACT.sat_extract import SatExtract

parser = argparse.ArgumentParser(description="Satellite extracts from multiple files (one file for variable) available in the CNR server.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-c', "--config_file", help="Config File.", required=True)
parser.add_argument('-sd', "--startdate", help="The Start Date - format YYYY-MM-DD ")
parser.add_argument('-ed', "--enddate", help="The End Date - format YYYY-MM-DD ")
parser.add_argument('-no_concat', "--no_concatenate", help="Use internaly for sbatch mode",action="store_true")
parser.add_argument('-make_concat', "--make_concatenate", help="Use internaly for sbatch mode to make final concatenation",action="store_true")
parser.add_argument('-p', "--product_file", help="Image file.")
args = parser.parse_args()


def add_variable_multiple(newEXTRACT, extract, variable_list, variable_list_out):
    list_files = extract['list_files']
    limits = extract['limits']
    start_idx_y = limits[0]
    stop_idx_y = limits[1]
    start_idx_x = limits[2]
    stop_idx_x = limits[3]
    nvar = len(list_files)
    for idx in range(nvar):
        file_in = list_files[idx]
        variable_in = variable_list[idx]
        variable_out = f'satellite_{variable_list_out[idx]}'
        nc_in = Dataset(file_in, 'r')
        if variable_in in nc_in.variables:
            var_in = nc_in.variables[variable_in]
        elif variable_in.upper() in nc_in.variables:
            var_in = nc_in.variables[variable_in.upper()]
        else:
            print(f'[ERROR] Variable {variable_in} is not available in the file. Extract not created')
            newEXTRACT.close_file()
            return None
        var_array = ma.array(var_in[:])
        var_array = np.array(var_array.filled(-999.0))

        if variable_out not in newEXTRACT.EXTRACT.variables:
            variable = newEXTRACT.create_2D_variable_general(variable_out, var_array, limits)
        else:
            variable = newEXTRACT.EXTRACT.variables[variable_out]
            variable[0, :, :] = var_array[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]

        for at in var_in.ncattrs():
            if at == '_FillValue' or at == 'add_offset' or at == 'scale_factor':
                continue
            variable.setncattr(at, var_in.getncattr(at))
        nc_in.close()

    return newEXTRACT


def add_reflectance_multiple(newEXTRACT, extract, wl_list):
    if not 'satellite_bands' in newEXTRACT.EXTRACT.dimensions:
        print(f'[ERROR] Dimension satellite bands is not defined')
        return
    if 'global_at' in extract:
        global_at = extract['global_at']
    if '1' in extract:
        global_at = extract['1']['global_at']
    list_files = extract['list_files']
    limits = extract['limits']
    start_idx_y = limits[0]
    stop_idx_y = limits[1]
    start_idx_x = limits[2]
    stop_idx_x = limits[3]

    nwl = len(list_files)

    if 'satellite_Rrs' in newEXTRACT.EXTRACT.variables:
        satellite_Rrs = newEXTRACT.EXTRACT.variables['satellite_Rrs']
    else:
        satellite_Rrs = newEXTRACT.create_rrs_variable(global_at['sensor'])
    wavelengths = []
    for iwl in range(nwl):
        wl = float(wl_list[iwl])
        wavelengths.append(float(wl))
        input_dataset = Dataset(list_files[iwl])
        for name, variable in input_dataset.variables.items():
            # wls = str(wl).replace('.', '_')
            wls = f'{wl:.2f}'
            wls = wls.replace('.', '_')
            if wls.endswith('_00'):
                wls = wls[:-3]
            if wls.find('_') > 0 and wls.endswith('0'):
                wls = wls[:-1]
            # print(name,
            #       '-------------------------------------------------------------------------------------------------------------------------------->',
            #       wls)
            ifind = name.find(wls)
            if ifind >= 0:
                if variable.ndim == 3:
                    bandarray = ma.array(variable[:, :, :])
                    satellite_Rrs[0, iwl, :, :] = bandarray[0, start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
                elif variable.ndim == 2:
                    bandarray = ma.array(variable[:, :])
                    satellite_Rrs[0, iwl, :, :] = bandarray[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]
        input_dataset.close()

    if 'satellite_bands' in newEXTRACT.EXTRACT.variables:
        satellite_bands = newEXTRACT.EXTRACT.variables['satellite_bands']
        satellite_bands[:] = wavelengths
    else:
        newEXTRACT.create_satellite_bands_variable(wavelengths)

    return newEXTRACT

##start the extract without using original file
def start_extract(ofname, extract_info):
    newEXTRACT = SatExtract(ofname)
    if not newEXTRACT.FILE_CREATED:
        print(f'[ERROR] File {ofname} could not be created')
        return False

    if args.verbose:
        print(f'[INFO] Starting file: {ofname}')

    window = extract_info['limits']
    newEXTRACT.set_global_attributes(extract_info['global_at'])

    newEXTRACT.create_dimensions(extract_info['size_box'], extract_info['n_bands'])
    newEXTRACT.create_lat_long_variables(extract_info['lat_array'], extract_info['lon_array'], window)
    newEXTRACT.create_satellite_time_variable(extract_info['satellite_time'])

    return newEXTRACT


def get_cmems_multiple_product_day(path_source, org, datehere, dataset_name_file, dataset_name_format_date,
                                   dataset_var_list):
    # path_day = path_source
    # if org is not None:
    #     if org == 'YYYYjjj':
    #         yearstr = datehere.strftime('%Y')
    #         jjjstr = datehere.strftime('%j')
    #         path_day = os.path.join(path_source, yearstr, jjjstr)
    path_day = sat_extract.get_path_date(path_source,org,datehere,False)
    if path_day is None:
        return None
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

def create_multiple_csv_sbatch(options,output_path,mp_options):
    import COMMON.sbatch_scripter as sbs

    temp_path = sat_extract.create_dir(os.path.join(output_path,f'Temp_{str(dt.now().timestamp()).replace(".","_")}'))
    if temp_path is None:
        print(f'[ERROR] Temporary path for sbatch files could not be created. Please review permissions.')
        return
    if args.verbose:
        print(f'[INFO] Temporary path: {temp_path}')
    path_csv = options['MULTIPLE_CSV_SELECTION']['path_csv']
    if not os.path.isdir(path_csv):
        print(f'[ERROR] Path to csv files {path_csv} was not found or is not a valid directory')
        return

    col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = sat_extract.get_csv_options(
        options, 'MULTIPLE_CSV_SELECTION')

    npool = os.cpu_count() if mp_options['ncores']<0 else 1 if mp_options['ncores']==0 else mp_options['ncores']
    nincluded = 0
    index_folder = 1
    folder_csv = os.path.join(temp_path,f'CSV_{index_folder}')
    sat_extract.create_dir(folder_csv)
    sbatch_files = []
    sbatch_log_files = []
    dir_code = os.path.dirname(sat_extract.__file__)
    line_py_base = f'python {dir_code}/sat_extract_CNR.py -c'

    for name in os.listdir(path_csv):
        if not name.endswith('csv'):
            continue
        try:
            file_csv = os.path.join(path_csv, name)
            pd.read_csv(file_csv, sep=col_sep)
        except:
            print(f'[ERROR] File {path_csv} is not a valid csv separated by {col_sep}')
            return
        if nincluded==npool:
            if args.verbose:
                print(f'[NFO] Creating sbatch file...')
            file_config = os.path.join(temp_path,f'config_file_{index_folder}.ini')
            options.set('MULTIPLE_CSV_SELECTION','path_csv',folder_csv)
            options.set('multiprocessing','use_sbatch','False')
            with open(file_config,'w') as configw:
                options.write(configw)
            file_sbatch = os.path.join(temp_path,f'sbatch_script_{index_folder}.slurm')
            scripter = sbs.SBATCH_SCRIPTER(file_sbatch)
            scripter.start_script(mp_options,True)
            scripter.add_blank_lines(2)
            scripter.add_line(f'{line_py_base} {file_config} -v')
            scripter.close_script()
            sbatch_files.append(file_sbatch)
            sbatch_log_files.append(os.path.join(temp_path,f'sbatch_script_log_{index_folder}.log'))
            # preparing next---
            nincluded = 0
            index_folder = index_folder + 1
            folder_csv = os.path.join(temp_path, f'CSV_{index_folder}')
            sat_extract.create_dir(folder_csv)
            ##-------------------------------

        if args.verbose:
            print(f'[INFO] Adding CSV file {name}')
        file_csv_copy = os.path.join(folder_csv,name)
        shutil.copy(file_csv,file_csv_copy)
        nincluded = nincluded +1

    #final sbatch file
    if args.verbose:
        print(f'[INFO] Creating final sbatch file...')
    file_config = os.path.join(temp_path, f'config_file_{index_folder}.ini')
    options.set('MULTIPLE_CSV_SELECTION', 'path_csv', folder_csv)
    print(options['MULTIPLE_CSV_SELECTION']['path_csv'])
    with open(file_config, 'w') as configw:
        options.write(configw)
    file_sbatch = os.path.join(temp_path, f'sbatch_script_{index_folder}.slurm')
    scripter = sbs.SBATCH_SCRIPTER(file_sbatch)
    scripter.start_script(mp_options, True)
    scripter.add_line(f'{line_py_base} {file_config} -v')
    scripter.close_script()
    sbatch_files.append(file_sbatch)
    sbatch_log_files.append(os.path.join(temp_path, f'sbatch_script_log_{index_folder}.log'))

    #file out sh
    file_out_sh = os.path.join(temp_path,f'launch_multiple_sbatch.sh')
    sbs.prepare_sh_script_with_multiple_sbatch(file_out_sh,sbatch_files,sbatch_log_files,mp_options['sbatch_max_cores'])
    print(f'[INFO] SH file: {file_out_sh} has been created.')
    if mp_options['sbatch_launch']:
        import subprocess
        cmd = f'sh {file_out_sh}'
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]Error lunching script: {err}')

def create_single_seabass_sbatch(options,output_path,mp_options):
    import COMMON.sbatch_scripter as sbs

    temp_path = sat_extract.create_dir(os.path.join(output_path, f'Temp_{str(dt.now().timestamp()).replace(".", "_")}'))
    if temp_path is None:
        print(f'[ERROR] Temporary path for sbatch files could not be created. Please review permissions.')
        return
    if args.verbose:
        print(f'[INFO] Temporary path: {temp_path}')

    path_seabass = options['SEABASS_SELECTION']['path_seabass']
    if not os.path.isfile(path_seabass):
        print(f'[ERROR] Path to SeaBass file {path_seabass} was not found or is not a valid file')
        return
    import MDB_builder.SB_support as SBr
    sb = SBr.readSB(path_seabass)
    var_date, var_time, var_lat, var_lon, format_date, format_time = sat_extract.get_seabass_options(
        options, 'SEABASS_SELECTION')
    var_ok = SBr.check_seabass_variables(sb, var_date, var_lat, var_lon)
    if not var_ok:
        print(f'[ERROR] The following variables are avaialble: {[x for x in sb.variables]}')
        return

    ##get list dates
    start_date, end_date = sat_extract.get_start_end_date_from_args(args)
    date_array, time_list = sat_extract.get_date_list_from_seabass(sb, var_date, var_time, format_date, format_time,start_date, end_date)
    if date_array is None:  ##it could be none by a parsing error or because there is not data in the temoporal range between start_date and end_date
        return
    date_array_unique = np.unique(date_array)
    date_array_unique = np.array([dt.strptime(x,'%Y-%m-%d').timestamp() for x in date_array_unique])
    date_array_unique = np.sort(date_array_unique)

    npool = os.cpu_count() if mp_options['ncores'] < 0 else 1 if mp_options['ncores'] == 0 else mp_options['ncores']
    min_ts_abs = None
    max_ts_abs = None
    date_ts_folder = []
    index_folder = 1
    sbatch_files = []
    sbatch_log_files = []
    dir_code = os.path.dirname(sat_extract.__file__)
    line_py_base = f'python {dir_code}/sat_extract_CNR.py -c'


    for date_ts in date_array_unique:
        if len(date_ts_folder)==npool:
            if args.verbose:
                print(f'[NFO] Creating sbatch file...')
            file_config = os.path.join(temp_path,f'config_file_{index_folder}.ini')
            file_sb = os.path.join(temp_path, f'SB_TEMP_{index_folder}.sb')
            shutil.copy(path_seabass, file_sb)
            options.set('SEABASS_SELECTION','path_seabass',file_sb)
            options.set('multiprocessing','use_sbatch','False')
            with open(file_config,'w') as configw:
                options.write(configw)
            file_sbatch = os.path.join(temp_path,f'sbatch_script_{index_folder}.slurm')
            scripter = sbs.SBATCH_SCRIPTER(file_sbatch)
            scripter.start_script(mp_options,True)
            scripter.add_blank_lines(2)
            min_ts = np.min(np.array(date_ts_folder))
            max_ts = np.max(np.array(date_ts_folder))
            if min_ts_abs is None:
                min_ts_abs = min_ts
                max_ts_abs = max_ts
            else:
                min_ts_abs = min_ts if min_ts < min_ts_abs else min_ts_abs
                max_ts_abs = max_ts if max_ts > max_ts_abs else max_ts_abs

            scripter.add_line(f'{line_py_base} {file_config} -sd {dt.fromtimestamp(min_ts).strftime("%Y-%m-%d")} -sd {dt.fromtimestamp(max_ts).strftime("%Y-%m-%d")} -v')
            scripter.close_script()
            sbatch_files.append(file_sbatch)
            sbatch_log_files.append(os.path.join(temp_path,f'sbatch_script_log_{index_folder}.log'))
            # preparing next---
            date_ts_folder = []
            index_folder = index_folder + 1
            ##-------------------------------

        if args.verbose:
            print(f'[INFO] Adding date: {dt.fromtimestamp(date_ts)}')
        date_ts_folder.append(date_ts)

    # final sbatch file
    if args.verbose:
        print(f'[INFO] Creating final sbatch file...')
    file_config = os.path.join(temp_path, f'config_file_{index_folder}.ini')
    file_sb = os.path.join(temp_path, f'SB_TEMP_{index_folder}.sb')
    shutil.copy(path_seabass, file_sb)
    options.set('SEABASS_SELECTION', 'path_seabass', file_sb)
    options.set('multiprocessing', 'use_sbatch', 'False')
    with open(file_config, 'w') as configw:
        options.write(configw)
    file_sbatch = os.path.join(temp_path, f'sbatch_script_{index_folder}.slurm')
    scripter = sbs.SBATCH_SCRIPTER(file_sbatch)
    scripter.start_script(mp_options, True)
    scripter.add_blank_lines(2)
    min_ts = np.min(np.array(date_ts_folder))
    max_ts = np.max(np.array(date_ts_folder))
    min_ts_abs = min_ts if min_ts < min_ts_abs else min_ts_abs
    max_ts_abs = max_ts if max_ts > max_ts_abs else max_ts_abs
    scripter.add_line(
        f'{line_py_base} {file_config} -sd {dt.fromtimestamp(min_ts).strftime("%Y-%m-%d")} -ed {dt.fromtimestamp(max_ts).strftime("%Y-%m-%d")} -no_concat -v')
    scripter.close_script()
    sbatch_files.append(file_sbatch)
    sbatch_log_files.append(os.path.join(temp_path, f'sbatch_script_log_{index_folder}.log'))

    #sbatch file for contanenation
    file_sbatch = os.path.join(temp_path, f'sbatch_script_contact.slurm')
    scripter = sbs.SBATCH_SCRIPTER(file_sbatch)
    scripter.start_script(mp_options, True)
    scripter.add_blank_lines(2)
    scripter.add_line(
        f'{line_py_base} {args.config_file} -sd {dt.fromtimestamp(min_ts_abs).strftime("%Y-%m-%d")} -ed {dt.fromtimestamp(max_ts_abs).strftime("%Y-%m-%d")} -make_concat -v')
    scripter.close_script()
    file_log = os.path.join(temp_path, f'sbatch_script_log_concat.log')

    # file out sh
    file_out_sh = os.path.join(temp_path, f'launch_multiple_sbatch.sh')
    sbs.prepare_sh_script_with_multiple_sbatch(file_out_sh, sbatch_files, sbatch_log_files,
                                               mp_options['sbatch_max_cores'])
    ##add contatenation sbatch
    fw = open(file_out_sh,'a')
    fw.write('\n')
    fw.write('\n')
    fw.write(f'sbatch --dependency=afterany:$jobid0:$jobid{len(sbatch_files)-1} --output={file_log} {file_sbatch})')
    fw.close()

    print(f'[INFO] SH file: {file_out_sh} has been created.')
    if mp_options['sbatch_launch']:
        import subprocess
        cmd = f'sh {file_out_sh}'
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]Error lunching script: {err}')

def run_multiple_csv(options,output_path, overwrite, ncores):
    path_csv = options['MULTIPLE_CSV_SELECTION']['path_csv']
    if not os.path.isdir(path_csv):
        print(f'[ERROR] Path to csv files {path_csv} was not found or is not a valid directory')
        return

    col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = sat_extract.get_csv_options(
        options, 'MULTIPLE_CSV_SELECTION')
    extract_options = sat_extract.get_cnr_extract_options(options)
    input_path_info = sat_extract.get_input_path_info(options)
    if input_path_info is None:
        return

    ##common for all days
    lat_array = None
    lon_array = None

    params_list = []

    for name in os.listdir(path_csv):
        if not name.endswith('csv'):
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

        only_date_array_unique, time_list = sat_extract.get_date_list_from_dataframe(df,col_date,format_date,col_time,format_time)
        if len(only_date_array_unique) == 0:
            print(f'[ERROR] No valid dates retrieved from the CSV file.')
            return
        if len(only_date_array_unique)>1:
            print(f'[ERROR] Each CSV file should also contain data for a single day. Check CSV parameters')
            return

        date_str = only_date_array_unique[0]
        if args.verbose:
                print(f'[INFO] Checking available products for day {date_str} in {input_path_info["path_source"]} with org. {input_path_info["org"]}')



        list_files = get_cmems_multiple_product_day(input_path_info['path_source'], input_path_info['org'],
                                                        dt.strptime(date_str, '%Y-%m-%d'),
                                                        extract_options['dataset_name_file'],
                                                        extract_options['dataset_name_format_date'],
                                                        extract_options['dataset_var_list'])
        if list_files is not None:
           satellite_time = sat_extract.get_satellite_time_l3(options,list_files[0],dt.strptime(date_str, '%Y-%m-%d'))
           if satellite_time is None:
               print(f'[WARNING] Satellite time could not be defined. Skipping...')
               continue
           extract_info= {
               'global_at': sat_extract.get_satellite_global_atrib_from_options(options),
               'list_files': list_files,
               'insitu_time': time_list,
               'insitu_lat': np.ma.array(df[col_lat]),
               'insitu_lon': np.ma.array(df[col_lon]),
               'satellite_time': satellite_time,
           }
           if lat_array is None and lon_array is None:
               lat_array, lon_array = sat_extract.get_lat_lon_arrays(options, list_files[0])

           if ncores==0:
               run_extract_day(extract_options,extract_info,lat_array,lon_array,output_path, overwrite)
           else:
               if ncores > 0 or ncores == -1:
                   params_list.append([extract_options,extract_info,lat_array,lon_array,output_path, overwrite])


    if ncores>0 or ncores<0:
        if args.verbose:
            print(f'[INFO] Starting parallel processing. Number of dates: {len(params_list)}')
            print(f'[INFO] CPUs: {os.cpu_count()}')
            print(f'[INFO] Parallel processes: {ncores}')
        npool = os.cpu_count() if ncores<0 else ncores
        poolhere = Pool(npool)
        poolhere.map(run_parallel_extract_day, params_list)

def run_single_seabass(options,output_path,overwrite, ncores, concatenate_daily_csv_extract):
    path_seabass = options['SEABASS_SELECTION']['path_seabass']
    if not os.path.isfile(path_seabass):
        print(f'[ERROR] Path to SeaBass file {path_seabass} was not found or is not a valid file')
        return

    path_extract_output = os.path.join(output_path,f'{path_seabass[0:path_seabass.rfind('.')]}_extracts.csv')
    if os.path.isfile(path_extract_output):
        try:
            os.remove(path_extract_output)
        except:
            print(f'[ERROR] Path extract output {path_extract_output} could not be remove. Please review permission or if the file is open')
            return
    if args.verbose:
        print(f'[INFO] Path extract output: {path_extract_output}')

    #from MDB_builder.SB_support import readSB
    import MDB_builder.SB_support as SBr
    sb = SBr.readSB(path_seabass)
    var_date, var_time, var_lat, var_lon, format_date, format_time = sat_extract.get_seabass_options(
        options, 'SEABASS_SELECTION')
    extract_options = sat_extract.get_cnr_extract_options(options)
    input_path_info = sat_extract.get_input_path_info(options)
    if input_path_info is None:
        return

    var_ok = SBr.check_seabass_variables(sb,var_date,var_lat,var_lon)
    if not var_ok:
        print(f'[ERROR] The following variables are avaialable: {[x for x in sb.variables]}')
        return

    start_date, end_date = sat_extract.get_start_end_date_from_args(args)
    date_array, time_list = sat_extract.get_date_list_from_seabass(sb,var_date,var_time,format_date,format_time,start_date,end_date)
    if date_array is None:##it could be none by a parsing error or because there is not data in the temoporal range between start_date and end_date
        return

    date_array_unique = np.unique(date_array)

    lat_array = np.array(sb.data[var_lat])
    lon_array = np.array(sb.data[var_lon])
    lat_image = None
    lon_image = None

    params_list = []
    path_extract_output_list = []
    path_extract_output_date = os.path.join(output_path,
                                            f'{path_seabass[0:path_seabass.rfind('.')]}_$DATE$_extracts.csv')
    ##remove files by date if they exist
    for date_str in date_array_unique:
        path_extract_output_here = path_extract_output_date.replace('$DATE$', date_str)
        if os.path.isfile(path_extract_output_here):
            try:
                os.remove(path_extract_output_here)
            except:
                print(
                    f'[ERROR] Path extract output date {path_extract_output_here} could not be remove. Please review permission or if the file is open')
                return

    for date_str in date_array_unique:
        indices_here = np.where(date_array==date_str)
        lat_here = lat_array[indices_here]
        lon_here = lon_array[indices_here]
        time_list_here = time_list[indices_here]
        if args.verbose:
            print(f'[INFO] Checking available products for day {date_str} in {input_path_info["path_source"]} with org. {input_path_info["org"]}')

        list_files = get_cmems_multiple_product_day(input_path_info['path_source'], input_path_info['org'],
                                                        dt.strptime(date_str, '%Y-%m-%d'),
                                                        extract_options['dataset_name_file'],
                                                        extract_options['dataset_name_format_date'],
                                                        extract_options['dataset_var_list'])

        if list_files is not None:
            satellite_time = sat_extract.get_satellite_time_l3(options, list_files[0],
                                                               dt.strptime(date_str, '%Y-%m-%d'))
            if satellite_time is None:
                print(f'[WARNING] Satellite time could not be defined. Skipping...')
                continue
            path_extract_output_here = path_extract_output_date.replace('$DATE$',date_str)
            path_extract_output_list.append(path_extract_output_here)
            extract_info = {
                'global_at': sat_extract.get_satellite_global_atrib_from_options(options),
                'list_files': list_files,
                'insitu_time': time_list_here,
                'insitu_lat': lat_here,
                'insitu_lon': lon_here,
                'insitu_indices': indices_here,
                'path_extract_output': path_extract_output_here,
                'satellite_time': satellite_time
            }
            if lat_image is None and lon_image is None:
                lat_image, lon_image = sat_extract.get_lat_lon_arrays(options, list_files[0])
            if ncores == 0:
                run_extract_day(extract_options, extract_info, lat_image, lon_image, output_path, overwrite)
            else:
                if ncores > 0 or ncores == -1:
                    params_list.append([extract_options, extract_info, lat_image, lon_image, output_path, overwrite])


    if ncores>0 or ncores<0:
        if args.verbose:
            print(f'[INFO] Starting parallel processing. Number of dates: {len(params_list)}')
            print(f'[INFO] CPUs: {os.cpu_count()}')
            print(f'[INFO] Parallel processes: {ncores}')
        npool = os.cpu_count() if ncores<0 else ncores
        poolhere = Pool(npool)
        poolhere.map(run_parallel_extract_day, params_list)


    if concatenate_daily_csv_extract:
        if args.verbose:
            print(f'[INFO] Concatenating date files...')
        sat_extract.concatenate_csv(path_extract_output_list,path_extract_output,True)
        if args.verbose:
            print(f'[INFO] Completed')

def run_parallel_extract_day(params):
    run_extract_day(params[0],params[1],params[2],params[3],params[4],params[5])

def run_extract_day(extract_options,extract_info,lat_array,lon_array,output_path, overwrite):
    insitu_lat = extract_info['insitu_lat']
    insitu_lon = extract_info['insitu_lon']
    insitu_time = extract_info['insitu_time']
    insitu_indices = extract_info['insitu_indices']
    path_extract_output = extract_info['path_extract_output']
    if not os.path.exists(path_extract_output):
        few = open(path_extract_output,'w')
        few.write('extract_file;insitu_index;insitu_time;insitu_lat;insitu_lon')
    else:
        few = open(path_extract_output,'a')

    ntimes = len(insitu_time)


    site_list = []
    format_datetime = '%Y-%m-%dT%H:%M:%S'

    for itime in range(ntimes):


        limits, rc = sat_extract.get_geo_info(extract_options, None, insitu_lat[itime], insitu_lon[itime], lat_array, lon_array)
        if limits is None:
            print(
                f'[WARNING] In situ location out of the limits of the satellite product. Skipping...')
            continue

        global_at = extract_info['global_at'].copy()
        datehere_str = extract_info['satellite_time'].strftime('%Y%m%d')
        site = f'{sat_extract.get_satellite_ref(global_at)}_{datehere_str}_{rc[0]}_{rc[1]}'
        ofname = os.path.join(output_path, f'extract_{site}.nc')

        if os.path.exists(ofname) and not overwrite:
            few.write('\n')
            few.write(
                f'extract_{site}.nc;{insitu_indices[0][itime]};{insitu_time[itime].strftime(format_datetime)};{insitu_lat[itime]};{insitu_lon[itime]}')
            print(f'[WARNING] [{itime + 1}/{ntimes}] Satellite extract extract_{site}.nc already exists. {itime+1}/{ntimes} Skiping...')
            continue

        if overwrite and site in site_list:
            few.write('\n')
            few.write(
                f'extract_{site}.nc;{insitu_indices[0][itime]};{insitu_time[itime].strftime(format_datetime)};{insitu_lat[itime]};{insitu_lon[itime]}')
            print(f'[WARNING] [{itime + 1}/{ntimes}] Satellite extract for {site}  has been already done.  Skiping...')
            continue

        if args.verbose:
            print(f'[INFO] [{itime + 1}/{ntimes}] Preparing extract for site {site}...')

        global_at_here = sat_extract.add_insitu_global_atrib(global_at, site, insitu_lat[itime], insitu_lon[itime], None)
        extract_info_here ={
            'global_at': global_at_here,
            'limits': limits,
            'size_box': extract_options['size_box'],
            'lat_array': lat_array,
            'lon_array': lon_array,
            'satellite_time': extract_info['satellite_time'],
            'n_bands': extract_options['n_bands'],
            'list_files': extract_info['list_files']
        }
        newExtract = start_extract(ofname, extract_info_here)
        if extract_options['is_reflectance']:
            newExtract = add_reflectance_multiple(newExtract, extract_info_here, extract_options['rrs_list'])
        else:
            newExtract = add_variable_multiple(newExtract, extract_info_here,
                                               extract_options['dataset_var_list'],
                                               extract_options['dataset_var_list_out'])
        if newExtract is None:
            print(f'[ERROR] Error creating extract for site {site}')
            os.remove(ofname)
            few.write('\n')
            few.write(
                f'noextract;{insitu_indices[0][itime]};{insitu_time[itime].strftime(format_datetime)};{insitu_lat[itime]};{insitu_lon[itime]}')
            continue
        else:
            print(f'[INFO] Extract {ofname} completed.')
            site_list.append(site)
            few.write('\n')
            few.write(f'extract_{site}.nc;{insitu_indices[0][itime]};{insitu_time[itime].strftime(format_datetime)};{insitu_lat[itime]};{insitu_lon[itime]}')
    few.close()

def run_time_and_sites_selection(options,output_path):
    extract_options = sat_extract.get_cnr_extract_options(options)
    if args.verbose:
        print('------------------------------')
    # product_list = {}
    # for date_here in date_list:
    #     if args.verbose:
    #         print(f'[INFO] Checking available products for day: {date_str}')
    #
    #     list_files = get_cmems_multiple_product_day(path_source, org, date_here,
    #                                                 extract_options['dataset_name_file'],
    #                                                 extract_options['dataset_name_format_date'],
    #                                                 extract_options['dataset_var_list'])

    # path_source, org, wce, time_start, time_stop = get_find_product_info(options)
    # path_to_list = get_path_list_products(options)
    # if path_to_list is not None:
    #     product_path_list = get_list_products(path_to_list, path_source, org, wce, time_start, time_stop)
    #     if len(product_path_list) == 0:
    #         if args.verbose:
    #             print('-----------------')
    #         print(f'[WARNING] No valid datasets found on:  {path_source}')
    #         print('------------------------------')
    #         print('[INFO] COMPLETED. 0 sat extract files were created')
    #         return
    #     ncreated = 0
    #     for filepath in product_path_list:
    #         if args.verbose:
    #             print('------------------------------')
    #             print(f'[INFO]Extracting sat extract for product: {filepath}')
    #         nhere = launch_create_extract(filepath, options)
    #         ncreated = ncreated + nhere
    #     print('------------------------------')
    #     print(f'COMPLETED. {ncreated} sat extract files were created')


def make_seabass_concatenation(options,output_path):
    start_date, end_date = sat_extract.get_start_end_date_from_args(args)
    if start_date is None or end_date is None:
        print(f'[ERROR] Start date and end date are required.')
        return
    path_seabass = options['SEABASS_SELECTION']['path_seabass']
    if not os.path.isfile(path_seabass):
        print(f'[ERROR] Path to SeaBass file {path_seabass} was not found or is not a valid file')
        return

    path_extract_output = os.path.join(output_path, f'{path_seabass[0:path_seabass.rfind('.')]}_extracts.csv')
    path_extract_output_date = os.path.join(output_path,
                                            f'{path_seabass[0:path_seabass.rfind('.')]}_$DATE$_extracts.csv')

    file_list = []
    work_date = start_date
    from datetime import timedelta
    while work_date<=end_date:
        file_date = path_extract_output_date.replace('$DATE$',work_date.strftime("%Y-%m-%d"))
        if os.path.isfile(file_date):
            file_list.append(file_date)
        work_date = work_date + timedelta(hours=24)
    sat_extract.concatenate_csv(file_list,path_extract_output,True)


def main():

    print('[INFO] Creating satellite extracts')
    if not args.config_file:
        return
    if not os.path.exists(args.config_file):
        print(f'[ERROR] File {args.config_file} does not exist')
        return

    options = sat_extract.config_reader(args.config_file)
    path_output = sat_extract.get_output_path(options,args.verbose)
    if path_output is None:
        print(f'ERROR: {path_output} is not valid')
        return
    overwrite = sat_extract.overwrite(options)
    mp_options = sat_extract.get_multiprocessing_options(options)
    if mp_options is None:
        return

    ##MULTIPLE CSV FILE SELECTION (e.g., FOR TARA METADATA FILES)
    if options.has_section('MULTIPLE_CSV_SELECTION') and options.has_option('MULTIPLE_CSV_SELECTION', 'path_csv'):
        if mp_options['use_sbatch']:
            create_multiple_csv_sbatch(options,path_output,mp_options)
        else:
            run_multiple_csv(options,path_output,overwrite,mp_options['ncores'])

    ##SINGLE SEABASS FILE SELECTION (e.g., FOR HYPERCP)
    if options.has_section('SEABASS_SELECTION'):
        if not options.has_option('SEABASS_SELECTION', 'path_seabass'):
            print(f'[ERROR] Option path_seabass is required for SEABASSS_SELECTION')
            return
        if args.make_concatenate and args.startdate and args.enddate:
            make_seabass_concatenation(options,path_output)
        else:
            if mp_options['use_sbatch']:
                create_single_seabass_sbatch(options,path_output,mp_options)
            else:
                run_single_seabass(options,path_output,overwrite,mp_options['ncores'],args.no_concatenate)



    if options.has_section('CSV_SELECTION') and options.has_option('CSV_SELECTION', 'path_csv'):
        # path_csv = options['CSV_SELECTION']['path_csv']
        # if not os.path.exists(path_csv):
        #     print(f'[ERROR] Path csv: {path_csv} does not exist')
        #     return
        #
        # with open(path_csv) as f:
        #     first_line = f.readline().strip()
        #
        # path_csv_out = f'{path_csv[:-4]}_out.csv'
        # if options.has_option('CSV_SELECTION', 'path_csv_out'):
        #     path_csv_out = options['CSV_SELECTION']['path_csv_out']
        #
        # fcsv_out = open(path_csv_out, 'w')
        # fcsv_out.write(f'{first_line};Extract;Index')
        #
        # use_single_file = False
        # n_bands = 0
        # if options.has_option('CSV_SELECTION', 'use_single_file'):  ##SINGLE FILE SELECTION
        #     usf = options['CSV_SELECTION']['use_single_file']
        #     if usf.strip().lower() == 'true' or usf.strip() == '1':
        #         use_single_file = True
        #
        # dataset_name_file = options['CSV_SELECTION']['dataset_name_file']
        # dataset_name_format_date = options['CSV_SELECTION']['dataset_name_format_date']
        # s = options['CSV_SELECTION']['dataset_var_list']
        # dataset_var_list = [x.strip() for x in s.split(',')]
        # dataset_var_list_out = dataset_var_list
        # if options.has_option('CSV_SELECTION', 'dataset_var_list_out'):
        #     s = options['CSV_SELECTION']['dataset_var_list_out']
        #     dataset_var_list_out = [x.strip() for x in s.split(',')]
        #
        # if use_single_file:
        #     rrs_list = []
        #     rrs_var_list = []
        #     is_reflectance = False
        #     if options.has_option('CSV_SELECTION', 'var_rrs_list') and options.has_option('CSV_SELECTION',
        #                                                                                   'var_rrs_format'):
        #
        #         var_rrs_list = options['CSV_SELECTION']['var_rrs_list'].strip()
        #         var_rrs_format = options['CSV_SELECTION']['var_rrs_format'].strip()
        #         is_reflectance = True
        #         for r in var_rrs_list.split(','):
        #             try:
        #                 rs = r.strip().replace('.', '_')
        #      # path_csv = options['MULTIPLE_CSV_SELECTION']['path_csv']
        # if not os.path.isdir(path_csv):
        #     print(f'[ERROR] Path to csv files {path_csv} was not found or is not a valid directory')
        #     return
        #
        # col_date, col_time, col_lat, col_lon, format_date, format_time, col_sep = get_csv_options_from_file_config(
        #     options, 'MULTIPLE_CSV_SELECTION')
        #
        # extract_options = get_cmems_extract_options(options, 'MULTIPLE_CSV_SELECTION')
        # path_source, org, wce, time_start, time_stop = get_find_product_info(options)
        # size_box = get_box_size(options)
        # # var_lat, var_lon = get_lat_long_var_names(options)
        # cmems_download_options = get_cmems_download_options(options)
        #
        # lat_array = None  ##CMEMS EXTRACTS WORKS ALWAYS WITH THE SAME LAT/LON ARRAYS FOR LOCATION
        # lon_array = None
        # ncreated = 0
        #
        # for name in os.listdir(path_csv):
        #     if not name.endswith('csv'):
        #         continue
        #     namefile = name[:-4]
        #
        #     try:
        #         file_csv = os.path.join(path_csv, name)
        #         df = pd.read_csv(file_csv, sep=col_sep)
        #     except:
        #         print(f'[ERROR] File {path_csv} is not a valid csv separated by {col_sep}')
        #         return
        #
        #     if args.verbose:
        #         print(f'[INFO] --------------------------------------------------------------------------------------')
        #         print(f'[INFO] Working with csv {namefile} Obtaining product(s)...')
        #
        #     ##1: STEP 1: Obtain the products for the date list
        #     date_array_ts = df[col_date]
        #     only_date_array = []
        #     for x in date_array_ts:
        #         try:
        #             only_date_array.append(dt.strptime(x, format_date).strftime('%Y-%m-%d'))
        #         except:
        #             continue
        #     only_date_array_unique = np.unique(only_date_array).tolist()
        #     product_list = {}
        #     for date_str in only_date_array_unique:
        #         if args.verbose:
        #             print(f'[INFO] Checking available products for day: {date_str}')
        #
        #         list_files = get_cmems_multiple_product_day(path_source, org, dt.strptime(date_str, '%Y-%m-%d'),
        #                                                         extract_options['dataset_name_file'],
        #                                                         extract_options['dataset_name_format_date'],
        #                                                         extract_options['dataset_var_list'])
        #         product_list[date_str] = list_files
        #
        #     ##2. STEP 2: Create extracts for each date
        #     for date_str in product_list:
        #
        #
        #         list_files = product_list[date_str]
        #         fproduct = list_files[0] if list_files is not None else None
        #
        #         if fproduct is None:
        #             print(f'[WARNING] No product(s) is(are) available for date: {date_str}')
        #             continue
        #         else:
        #             if args.verbose:
        #                 print(f'[INFO] Working with {len(list_files)} files for date {date_str}')
        #
        #         if lat_array is None or lon_array is None:
        #             lat_array, lon_array = get_lat_lon_arrays(options, fproduct)
        #
        #         ##CHECKING EXTRACT INFO OPTION (INCLUDING SATELLITE TIME)
        #         cmems_time = '11:00'
        #         if options.has_option('satellite_options', 'satellite_time'):
        #             cmems_time = options['satellite_options']['satellite_time'].strip()
        #         satellite_time = get_satellite_time_from_global_attributes(fproduct)
        #         if satellite_time is None:
        #             try:
        #                 datehere_str = dt.strptime(date_str, '%Y-%m-%d').strftime('%Y%m%d')
        #                 satellite_time = dt.strptime(f'{datehere_str}T{cmems_time}',
        #                                              '%Y%m%dT%H:%M').replace(tzinfo=pytz.utc)
        #             except:
        #                 print(f'{cmems_time} is not a valid satellite time option. Skipping')
        #                 continue
        #         extract_info = {
        #             'global_at': None,
        #             'satellite_time': satellite_time,
        #             'size_box': size_box,
        #             'n_bands': extract_options['n_bands'],
        #             'limits': None,
        #             'lat_array': lat_array,
        #             'lon_array': lon_array,
        #             'file': fproduct,
        #             'list_files': list_files
        #         }
        #
        #         ##CHECKING EXTRACT FOR EACH ROW
        #         for idx, row in df.iterrows():
        #             datehere, lathere, lonhere = get_info_from_row(row, col_date, col_time, format_date, format_time,
        #                                                            col_lat, col_lon)
        #             if datehere is None or lathere is None or lonhere is None:
        #                 print(f'[WARNING] Row {idx} is not valid. Date, latitute and/or longitude could not be parsed')
        #                 continue
        #             if datehere.strftime('%Y-%m-%d') != date_str:
        #                 continue
        #             limits, rc = get_geo_info(options, fproduct, lathere, lonhere, lat_array, lon_array)
        #             if limits is None:
        #                 print(
        #                     f'[WARNING] In situ location out of the limits of the satellite product. Skipping...')
        #                 continue
        #
        #             global_at = get_satellite_global_atrib_from_options(options)
        #             datehere_str = datehere.strftime('%Y%m%d')
        #             site = f'{get_satellite_ref(global_at)}_{datehere_str}_{rc[0]}_{rc[1]}'
        #             ofname = os.path.join(path_output, f'extract_{site}.nc')
        #             if os.path.exists(ofname):
        #                 print(f'[WARNING] Satellite extract extract_{site}.nc already exists. Skiping...')
        #                 continue
        #             global_at = add_insitu_global_atrib(global_at, site, lathere, lonhere, None)
        #             extract_info['global_at'] = global_at
        #             extract_info['limits'] = limits
        #
        #             newExtract = start_extract_2(ofname, extract_info)
        #
        #
        #             if extract_options['is_reflectance']:
        #                 newExtract = add_reflectance_multiple(newExtract, extract_info, extract_options['rrs_list'])
        #             else:
        #                 newExtract = add_variable_multiple(newExtract, extract_info,extract_options['dataset_var_list'],extract_options['dataset_var_list_out'])
        #
        #             if newExtract is None:
        #                 os.remove(ofname)
        #                 continue
        #
        #             newExtract.close_file()
        #             ncreated = ncreated + 1
        #
        # print(f'[INFO] Extract generation for MULTIPLE_CSV_FILE was completed. ')
        # print(f'[INFO] {ncreated} extracts were created')
        # return           rrs_here = float(r.strip())
        #                 rrs_list.append(rrs_here)
        #                 var_here = var_rrs_format.replace('$BAND$', rs)
        #                 rrs_var_list.append(var_here)
        #             except:
        #                 is_reflectance = False
        #                 break
        #         if is_reflectance:
        #             n_bands = len(rrs_list)
        # else:
        #     is_reflectance = True
        #     for var in dataset_var_list:
        #         try:
        #             float(var)
        #         except:
        #             is_reflectance = False
        #     if is_reflectance:
        #         n_bands = len(dataset_var_list)
        #
        # try:
        #     df = pd.read_csv(path_csv, sep=';')
        # except:
        #     print(f'[ERROR] File {path_csv} is not a valid csv separated by ;')
        #     return
        # col_date = 'date'
        # col_lat = 'lat'
        # col_lon = 'lon'
        # format_date = '%Y-%m-%dT%H:%M'
        # if options.has_option('CSV_SELECTION', 'col_date'):
        #     col_date = options['CSV_SELECTION']['col_date']
        # if options.has_option('CSV_SELECTION', 'col_lat'):
        #     col_lat = options['CSV_SELECTION']['col_lat']
        # if options.has_option('CSV_SELECTION', 'col_lon'):
        #     col_lon = options['CSV_SELECTION']['col_lon']
        # if options.has_option('CSV_SELECTION', 'format_date'):
        #     format_date = options['CSV_SELECTION']['format_date']
        # col_time = None
        # if options.has_option('CSV_SELECTION', 'col_time'):
        #     col_time = options['CSV_SELECTION']['col_time']
        #     format_time = '%H:%M:%S'
        # if options.has_option('CSV_SELECTION', 'format_time'):
        #     format_time = options['CSV_SELECTION']['format_time']
        # csv_flags = None
        # csv_flags_meanings = None
        # if options.has_option('CSV_SELECTION', 'csv_flags'):
        #     csv_flags_str = options['CSV_SELECTION']['csv_flags']
        #     csv_flags = [x.strip() for x in csv_flags_str.split(',')]
        #
        # path_source, org, wce, time_start, time_stop = get_find_product_info(options)
        # size_box = get_box_size(options)
        # var_lat, var_lon = get_lat_long_var_names(options)
        #
        # extract_list = {}
        # if csv_flags is not None:
        #     csv_flags_meanings = {}
        #
        # for idx, row in df.iterrows():
        #     list = [str(x).strip() for x in row.to_list()]
        #     line_orig = ";".join(list)
        #     try:
        #         datehere = None
        #         if col_time is not None:
        #             try:
        #                 datetimerow = f'{row[col_date].strip()}T{row[col_time].strip()}'
        #                 if np.isreal(datetimerow):
        #                     datetimerow = f'{datetimerow:.0f}'
        #                 format_datetime = f'{format_date}T{format_time}'
        #                 datehere = dt.strptime(datetimerow, format_datetime).replace(tzinfo=pytz.utc)
        #             except:
        #                 pass
        #         if datehere is None:
        #             datetimerow = row[col_date]
        #             if np.isreal(datetimerow):
        #                 datetimerow = f'{datetimerow:.0f}'
        #             format_datetime = format_date
        #             datehere = dt.strptime(datetimerow, format_datetime).replace(tzinfo=pytz.utc)
        #             datehere = datehere.replace(hour=12, minute=0, second=0).replace(tzinfo=pytz.utc)
        #     except:
        #         print(f'[WARNING] Row {idx} is not valid. Date/Time could not be parsed. Skipping...')
        #         fcsv_out.write('\n')
        #         fcsv_out.write(f'{line_orig};NaN;-1')
        #         continue
        #     lathere = float(row[col_lat])
        #     lonhere = float(row[col_lon])
        #     if np.isnan(lathere) or np.isnan(lonhere):
        #         print(f'[WARNING] Row {idx} is not valid. Latitude or longitude could not be parsed. Skipping...')
        #         fcsv_out.write('\n')
        #         fcsv_out.write(f'{line_orig};NaN;-1')
        #         continue
        #     if csv_flags is not None:
        #         for f in csv_flags:
        #             val = row[f].strip()
        #             val = val.replace(' ', '_')
        #             if f not in csv_flags_meanings.keys():
        #                 csv_flags_meanings[f] = [val]
        #             else:
        #                 lflags = csv_flags_meanings[f]
        #                 if val not in lflags:
        #                     lflags.append(val)
        #                     csv_flags_meanings[f] = lflags
        #
        #
        #
        #
        #     list_files = get_cmems_multiple_product_day(path_source, org, datehere, dataset_name_file,
        #                                                     dataset_name_format_date, dataset_var_list)
        #
        #     if list_files is not None:
        #
        #         limits, rc = get_geo_info(options, list_files[0], lathere, lonhere, None, None)
        #         if limits is not None:
        #             global_at = get_satellite_global_atrib_from_options(options)
        #             datehere_str = datehere.strftime('%Y%m%d')
        #             site = f'{get_satellite_ref(global_at)}_{datehere_str}_{rc[0]}_{rc[1]}'
        #             # print(f'[INFO] Site: {site}')
        #             other = None
        #             if csv_flags is not None:
        #                 other = {}
        #                 for f in csv_flags:
        #                     val = row[f].strip()
        #                     val = val.replace(' ', '_')
        #                     other[f] = val
        #
        #             global_at = add_insitu_global_atrib(global_at, site, lathere, lonhere, other)
        #
        #             ofname = os.path.join(path_output, f'extract_{site}.nc')
        #             cmems_time = '11:00'
        #             if options.has_option('satellite_options', 'satellite_time'):
        #                 cmems_time = options['satellite_options']['satellite_time'].strip()
        #             try:
        #                 satellite_time = dt.strptime(f'{datehere_str}T{cmems_time}', '%Y%m%dT%H:%M').astimezone(
        #                         pytz.utc)
        #                 print(satellite_time)
        #             except:
        #                 print(f'{cmems_time} is not a valid satellite time option. Skipping')
        #                 fcsv_out.write('\n')
        #                 fcsv_out.write(f'{line_orig};NaN;-1')
        #                 continue
        #
        #             if site not in extract_list.keys():
        #                 extract_list[site] = {
        #                     'ninsitu': 1,
        #                     'satellite_time': satellite_time,
        #                     'ofname': ofname,
        #                     'limits': limits,
        #                     'list_files': list_files,
        #                     'size_box': size_box,
        #                     'n_bands': n_bands,
        #                     'var_lat': var_lat,
        #                     'var_lon': var_lon,
        #                     '1': {
        #                         'insitu_time': datehere,
        #                         'global_at': global_at,
        #                     }
        #                 }
        #                 fcsv_out.write('\n')
        #                 fcsv_out.write(f'{line_orig};{site}.nc;0')
        #             else:
        #                 idx = extract_list[site]['ninsitu'] + 1
        #                 fcsv_out.write('\n')
        #                 fcsv_out.write(f'{line_orig};{site}.nc;{idx - 1}')
        #                 idxs = str(idx)
        #                 extract_list[site]['ninsitu'] = idx
        #                 extract_list[site][idxs] = {
        #                     'insitu_time': datehere,
        #                     'global_at': global_at,
        #                 }
        #         else:
        #             fcsv_out.write('\n')
        #             fcsv_out.write(f'{line_orig};NaN;-1')
        #     else:
        #         print(f'[WARNING] Files not found for date: {datehere.strftime("%Y-%m-%d")}. Skipping...')
        #         fcsv_out.write('\n')
        #         fcsv_out.write(f'{line_orig};NaN;-1')
        #
        # fcsv_out.close()
        #
        # for site in extract_list:
        #     extract = extract_list[site]
        #     ofname = extract['ofname']
        #     ofname_temp = f'{ofname[:-3]}_temp.nc'
        #     if os.path.exists(ofname):
        #         newExtract = copy_extract(ofname, ofname_temp)
        #     else:
        #         newExtract = start_extract(extract, ofname)
        #
        #     if is_reflectance:
        #         newExtract = add_reflectance_multiple(newExtract, extract, dataset_var_list)
        #     else:
        #         newExtract = add_variable_multiple(newExtract, extract, dataset_var_list, dataset_var_list_out)
        #     if newExtract is None:
        #             os.remove(ofname)
        #             continue
        #     nidx = 50
        #
        #     newExtract = add_insitu_basic_info(newExtract, extract, 1, nidx, csv_flags_meanings)
        #
        #     nhere = extract['ninsitu']
        #     if nhere > 1:
        #         for idx in range(2, nhere + 1):
        #             newExtract = add_insitu_basic_info(newExtract, extract, idx, nidx, csv_flags_meanings)
        #     newExtract.close_file()
        #
        #     if os.path.exists(ofname_temp):
        #         os.rename(ofname_temp, ofname)
        return



    if options.has_section('TIME_AND_SITES_SELECTION'):
        run_time_and_sites_selection(options,path_output)
        return


    #if options.has_option('file_path', 'path_insitu') and options.has_option('file_path', 'path_insitu_code'):
    #    run_insitu_option(options)
    #    return




    #use_single_file = extract_options['use_single_file']





    # work only with the specified product file
    # if args.product_file:
    #     print('------------------------------')
    #     filepath = args.product_file
    #     if args.verbose:
    #         print(f'[INFO]Extracting sat extract for product: {filepath}')
    #     ncreated = launch_create_extract(filepath, options)
    #     print('------------------------------')
    #     print(f'[INFO] COMPLETED. {ncreated} sat extract files were created')
    # else:
    #     print('------------------------------')
    #     path_source, org, wce, time_start, time_stop = get_find_product_info(options)
    #     path_to_list = get_path_list_products(options)
    #     if path_to_list is not None:
    #         product_path_list = get_list_products(path_to_list, path_source, org, wce, time_start, time_stop)
    #         if len(product_path_list) == 0:
    #             if args.verbose:
    #                 print('-----------------')
    #             print(f'[WARNING] No valid datasets found on:  {path_source}')
    #             print('------------------------------')
    #             print('[INFO] COMPLETED. 0 sat extract files were created')
    #             return
    #         ncreated = 0
    #         for filepath in product_path_list:
    #             if args.verbose:
    #                 print('------------------------------')
    #                 print(f'[INFO]Extracting sat extract for product: {filepath}')
    #             nhere = launch_create_extract(filepath, options)
    #             ncreated = ncreated + nhere
    #         print('------------------------------')
    #         print(f'COMPLETED. {ncreated} sat extract files were created')


if __name__ == '__main__':
    main()


 # if options.has_option('file_path', 'path_skie') and options.has_option('file_path', 'path_skie_code'):
    #     path_skie = options['file_path']['path_skie']
    #     path_skie_code = options['file_path']['path_skie_code']
    #     if not os.path.exists(path_skie):
    #         print(f'[ERROR] Skie file: {path_skie} is not avialable')
    #         return
    #     if not os.path.exists(path_skie_code) or not os.path.isdir(path_skie_code):
    #         print(f'[ERROR] Skie code folder: {path_skie_code} is not avialable')
    #         return
    #     if args.verbose:
    #         print(f'[INFO] Path skie: {path_skie}')
    #         print(f'[INFO] Path skie: {path_skie_code}')
    #     sys.path.append(path_skie_code)
    #     try:
    #         from skie_csv import SKIE_CSV
    #     except ModuleNotFoundError:
    #         print(f'[ERROR] Skie module is not found')
    #         return
    #     skie_file = SKIE_CSV(path_skie)
    #     skie_file.start_list_dates()
    #     skie_file.extract_wl_colnames()
    #
    #     print('------------------------------')
    #     path_source, org, wce, time_start, time_stop = get_find_product_info(options)
    #     path_to_list = get_path_list_products(options)
    #     if path_to_list is not None:
    #         product_path_list = get_list_products(path_to_list, path_source, org, wce, time_start, time_stop)
    #         if len(product_path_list) == 0:
    #             if args.verbose:
    #                 print('-----------------')
    #             print(f'[WARNING] No valid datasets found on:  {path_source}')
    #             print('------------------------------')
    #             print('[INFO] COMPLETED. 0 sat extract files were created')
    #             return
    #         ncreated = 0
    #         for filepath in product_path_list:
    #             if args.verbose:
    #                 print('------------------------------')
    #                 pname = filepath.split('/')[-1]
    #                 print(f'[INFO]Extracting sat extract for product: {pname}')
    #             nhere = launch_create_extract_skie(filepath, skie_file, options)
    #             ncreated = ncreated + nhere
    #         print('------------------------------')
    #         print(f'COMPLETED. {ncreated} sat extract files were created')
    #     return

#     for idx, row in df.iterrows():
#         datehere, lathere, lonhere = get_info_from_row(row, col_date, col_time, format_date, format_time,
#                                                        col_lat, col_lon)
#         if datehere is None or lathere is None or lonhere is None:
#             print(f'[WARNING] Row {idx} is not valid. Date, latitute and/or longitude could not be parsed')
#             continue
#         if extract_options['use_single_file']:
#
#             fproduct = get_cmems_product_day(path_source, org, datehere, extract_options['dataset_name_file'],
#                                              extract_options['dataset_name_format_date'],
#                                              cmems_download_options,False)
#
#             if fproduct is None:##check also myint
#                 fproduct = get_cmems_product_day(path_source, org, datehere, extract_options['dataset_name_file'],
#                                              extract_options['dataset_name_format_date'],
#                                              cmems_download_options,True)
#
#             if fproduct is not None:
#
#                 limits, rc = get_geo_info(options, fproduct, lathere, lonhere, None, None)
#
#                 if limits is not None:
#                     global_at = get_satellite_global_atrib_from_options(options)
#                     datehere_str = datehere.strftime('%Y%m%d')
#                     site = f'{get_satellite_ref(global_at)}_{datehere_str}_{rc[0]}_{rc[1]}'
#                     print(f'[INFO] Site: {site}')
#                     other = None
#                     #         if csv_flags is not None:
#                     #             other = {}
#                     #             for f in csv_flags:
#                     #                 val = row[f].strip()
#                     #                 val = val.replace(' ', '_')
#                     #                 other[f] = val
#                     #
#                     global_at = add_insitu_global_atrib(global_at, site, lathere, lonhere, other)
#
#                     ofname = os.path.join(path_output, f'extract_{site}.nc')
#                     cmems_time = '11:00'
#                     if options.has_option('satellite_options', 'satellite_time'):
#                         cmems_time = options['satellite_options']['satellite_time'].strip()
#
#                     satellite_time = get_satellite_time_from_global_attributes(fproduct)
#                     if satellite_time is None:
#                         try:
#                             satellite_time = dt.strptime(f'{datehere_str}T{cmems_time}',
#                                                          '%Y%m%dT%H:%M').replace(tzinfo=pytz.utc)
#                         except:
#                             print(f'{cmems_time} is not a valid satellite time option. Skipping')
#                             continue
#
#                     if site not in extract_list.keys():
#                         extract_list[site] = {
#                             'ninsitu': 1,
#                             'satellite_time': satellite_time,
#                             'ofname': ofname,
#                             'limits': limits,
#                             'file': fproduct,
#                             'size_box': size_box,
#                             'n_bands': extract_options['n_bands'],
#                             'var_lat': var_lat,
#                             'var_lon': var_lon,
#                             '1': {
#                                 'insitu_time': datehere,
#                                 'global_at': global_at,
#                             }
#                         }
#                     else:
#                         idx = extract_list[site]['ninsitu'] + 1
#                         idxs = str(idx)
#                         extract_list[site]['ninsitu'] = idx
#                         extract_list[site][idxs] = {
#                             'insitu_time': datehere,
#                             'global_at': global_at,
#                         }
#                 else:
#                     print(
#                         f'[WARNING] In situ location out of the limites of the satellite product. Skipping...')
#             else:
#                 print(f'[WARNING] Files not found for date: {datehere.strftime("%Y-%m-%d")}. Skipping...')
#
# if extract_options['use_single_file']:
#     for site in extract_list:
#         extract = extract_list[site]
#         nhere = extract['ninsitu']
#         ofname = extract['ofname']
#         if os.path.exists(ofname):
#             print(f'[WARNING] {ofname} already exist. Skipping...')
#             continue
#         else:
#             newExtract = start_extract(extract, ofname)
#
#         print(f'[INFO] Site: {site} Number of spectra: {nhere}')
#
#         if extract_options['is_reflectance']:
#             newExtract = add_reflectance_single(newExtract, extract, extract_options['rrs_list'],
#                                                 extract_options['rrs_var_list'])
#
#         newExtract = add_variable_single(newExtract, extract, extract_options['dataset_var_list'],
#                                          extract_options['dataset_var_list_out'],
#                                          extract_options['rrs_var_list'])
#         # nidx = 50
#         # newExtract = add_insitu_basic_info(newExtract, extract, 1, nidx, None)
#         # if nhere > 1:
#         #     for idx in range(2, nhere + 1):
#         #         newExtract = add_insitu_basic_info(newExtract, extract, idx, nidx, None)
#         newExtract.close_file()