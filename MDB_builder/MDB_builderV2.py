import os
import argparse
import configparser

import pandas as pd
import pytz

from MDB_builder_options import MDBBuilderOptions
from SATEXTRACTS_list import SAT_EXTRACTS_LIST
from INSITU_hypernets import INSITU_HYPERNETS_DAY
from datetime import datetime as dt

parser = argparse.ArgumentParser(
    description="Create Match-up DataBase files (MDB) files from satellite extracts and in situ L2 HYPERNETS files.")

parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-o', "--output",
                    help="Output file. Required with --listdates or single concatenation")
parser.add_argument('-edir', "--sat_extract_dir",
                    help="Input sat. extract dir. Optional for --listdates, required for single concatenation")
parser.add_argument('-site', "--sitename", help="Site name. Only required with --listdates")
parser.add_argument('-ld', "--listdates",
                    help="Option to obtain a date list for a specific HYPERNETS site (-site option).",
                    action="store_true")
parser.add_argument('-sd', "--start_date", help="Start date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-ed', "--end_date", help="End date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-nd', "--nodelfiles", help="Do not delete temp files.", action="store_true")
parser.add_argument('-chvar',"--check_variables_extract", help="Check variables of extracts in single CSV files", action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")

args = parser.parse_args()

import MDB_builder_options
import sys

code_home = os.path.dirname(os.path.dirname(MDB_builder_options.__file__))
# code_home = os.path.abspath('../')
sys.path.append(code_home)

# import BRDF.brdf_olci as brdf
from COMMON import common_functions as cfs


def test():
    print('test')
    from INSITU_tara import INSITU_TARA
    im = INSITU_TARA(None, True)
    from datetime import datetime as dt
    date_here = dt(2023, 4, 3)
    im.retrieve_metadata_from_file(date_here)
    file_extract = '/mnt/c/DATA_LUIS/TARA_TEST/extracts/S3B_OL_2_WFR____20230403T111113_20230403T111413_20230404T221050_0179_078_037_2160_MAR_O_NT_003_SEN3_extract_549_4638.nc'
    from netCDF4 import Dataset
    dataset = Dataset(file_extract)
    lat_array = dataset.variables['satellite_latitude'][0, :, :]
    lon_array = dataset.variables['satellite_longitude'][0, :, :]
    import numpy as np
    index_time = float(np.array(dataset.variables['satellite_time'][0]))
    sat_time = dt.fromtimestamp(index_time).replace(tzinfo=pytz.utc)
    dataset.close()
    max_time_diff = 180 * 60
    im.check_match_up_conditions(sat_time, lat_array, lon_array, max_time_diff)

    return True


def main():
    print('[INFO] Creating MDB files!')

    # b = test()
    # if b:
    #     return

    # Option to create a date list (HYPSTAR)
    if args.sitename and args.output and args.listdates:

        sat_extract_dir = None
        if args.sat_extract_dir:
            sat_extract_dir = args.sat_extract_dir
        start_date = None
        end_date = None
        if args.start_date and args.end_date:
            try:
                start_date = dt.strptime(args.start_date, '%Y-%m-%d')
                end_date = dt.strptime(args.end_date, '%Y-%m-%d')
            except:
                print(
                    f'[ERROR] Parameters start date: {args.start_date} and/or end date: {args.end_date} are not in the correct format (YYYY-mm-dd)')
                return
        if args.sitename == 'LAIT':
            from INSITU_meda import INSITU_MEDA
            idm = INSITU_MEDA(None, args.verbose)
            idm.get_datelist_file(args.output, start_date, end_date)
        else:
            ihd = INSITU_HYPERNETS_DAY(None, None, args.verbose)
            ihd.save_list_dates_to_file(args.output, args.sitename, start_date, end_date, sat_extract_dir)
        return

    # option to make single extract concatenation from files in a input path
    if args.sat_extract_dir and args.output:
        output_file = args.output
        if os.path.isdir(output_file):
            print(f'[ERROR] Output should be a NetCDF file name')
            return
        if not output_file.endswith('.nc'):
            output_file = f'{output_file}.nc'
        if os.path.exists(output_file):
            print(f'[WARNING] Output file {output_file} alreaday exist. Skipping')
            return
        output_path = os.path.dirname(output_file)
        if not os.path.isdir(output_path):
            print(f'[ERROR] Output path {output_path} does not exist.')
            return
        input_path = args.sat_extract_dir
        if not os.path.isdir(input_path):
            print(f'[ERROR] Input path {input_path} should be a directory')

        input_files = []
        for name in os.listdir(input_path):
            if name.endswith('.nc'):
                input_file = os.path.join(input_path, name)
                input_files.append(input_file)
        if len(input_files) == 0:
            print(f'[WARNING] No input NC files found in {input_path}')
            return
        concatenate_nc_impl(input_files, output_path, output_file)
        return

    if not args.config_file:
        print('[ERROR] Configuration file (-c or --config_file) option is required')
        return
    if os.path.isfile(args.config_file):
        options = configparser.ConfigParser()
        options.read(args.config_file)
        mo = MDBBuilderOptions(options, args.verbose)
    else:
        print('[ERROR] Configuration file does not exist')
        return

    valid_options = mo.set_compulsory_options(cfs)
    if not valid_options:
        return

    # SYKE or INSITU MODE, in situ are already included in the minifiles, builder makes a simple concatenation
    if mo.insitu_type == 'SYKE' or mo.insitu_type == 'INSITU':
        print(f'[INFO] Entering simple concatenation mode...')
        ##make_simple_builder(options, satellite_path_source, path_out)

    ##sat extracts options
    mo.get_param_sat_extracts()

    ##in situ options
    mo.get_insitu_options()

    ##dates
    if not mo.insitu_type.startswith('SINGLE_CSV'):
        if args.verbose:
            print(f'[INFO] Checking available dates--------------------------------------------------------------START')
        mo.get_dates()
        if args.verbose:
            print(f'[INFO] Start date for MDB_builder: {mo.start_date}')
            print(f'[INFO] End date for MDB_builder:  {mo.end_date}')
            print(f'[INFO] Checking available dates--------------------------------------------------------------STOP')

    ##retrieving sat extract list (except for SINGLE_CSV mode, with extract name included in the csv file)
    extract_list = None
    if not mo.insitu_type.startswith('SINGLE_CSV'):
        if args.verbose:
            print(f'[INFO] Obtaining extract list----------------------------------------------------------------START')
        slist = SAT_EXTRACTS_LIST(mo, args.verbose)
        extract_list = slist.get_list_as_dict()
        if args.verbose:
            print(f'[INFO] Obtaining extract list----------------------------------------------------------------STOP')

    ##checking in situ files
    if args.verbose:
        print(f'[INFO] Generating MDB extract files----------------------------------------------------------START')

    if mo.insitu_type == 'WISP3':
        insitu_file = mo.get_single_insitu_file()
        if insitu_file is not None:
            create_mdb_wisp_file(mo, extract_list, insitu_file)
        return

    if mo.insitu_type == 'MEDA':
        create_mdb_meda_file(mo, extract_list)
        return

    if mo.insitu_type == 'TARA':
        create_mdb_tara_file(mo, extract_list)
        return

    if mo.insitu_type == 'SINGLE_CSV_RRS':
        insitu_file = mo.get_single_insitu_file()
        create_mdb_single_csv_rrs(mo, extract_list, insitu_file)
        return

    if mo.insitu_type == 'SINGLE_CSV_VAR':
        insitu_file = mo.get_single_insitu_file()
        mo.get_dates_from_insitu_file()
        #print(mo.get_mdb_path())
        create_mdb_single_csv_var(mo, insitu_file)
        return

    ihd = INSITU_HYPERNETS_DAY(mo, None, args.verbose)
    if args.verbose:
        if mo.insitu_options['apply_rsync']:
            print(f'[INFO] Checking SSH access: {ihd.CHECK_SSH}')
        time_maxh = ihd.mdb_options.insitu_options['time_window'] / 3600
        print(f'[INFO] Maximum time window: {time_maxh:0.2f} hours')
    ins_sensor = 'HYPSTAR'
    mdb_extract_files = []
    for extract in extract_list:
        if args.verbose:
            print(f'[INFO] Working with extract: {extract} *******************')
        bad_spectra_times = {}
        extract_name = os.path.basename(extract_list[extract]['path'])
        ofile = mo.get_mdb_extract_path(extract_name, ins_sensor)
        if os.path.exists(ofile):
            from netCDF4 import Dataset
            import numpy as np
            dtal = Dataset(ofile)
            if 'insitu_original_bands' not in dtal.variables:
                dtal.close()
                print(f'[ERROR] Error in MDB extract file: {ofile}. Skipping...')
                continue
            insitu_bands = np.array(dtal.variables['insitu_original_bands'])
            vtal = insitu_bands[0]
            dtal.close()
            if vtal == -999.0:
                print(f'[ERROR] Error in MDB extract file. Skipping...')
                continue
            mdb_extract_files.append(ofile)
            print(f'[WARNING] MDB extract file already exits. Skipping...')
            continue
        # date_here_str = extract_list[extract]['time']
        # date_here = dt.strptime(date_here_str, '%Y%m%dT%H%M%S')
        date_here = extract_list[extract]['time']
        insitu_files = ihd.get_insitu_files(date_here)
        if insitu_files is None and mo.insitu_options['apply_rsync'] and ihd.check_ssh():
            ihd.get_files_day_ssh(date_here, True)
            insitu_files = ihd.get_insitu_files(date_here)

        # print(mo.insitu_options)
        if mo.insitu_options['bad_spectra_file_list'] is not None:
            prefix = mo.insitu_options['bad_spectra_prefix']
            time_format = mo.insitu_options['bad_spectra_format_time']
            f1 = open(mo.insitu_options['bad_spectra_file_list'])
            for line in f1:
                if prefix is not None:
                    if not line.strip().startswith(prefix):
                        continue
                    datestr = line.replace(prefix, '').strip()
                else:
                    datestr = line.strip()
                date_py = dt.strptime(datestr, time_format)
                date_py_str = date_py.strftime('%Y%m%d%H%M')
                bad_spectra_times[date_py_str] = 1
            f1.close()

        if len(bad_spectra_times) > 0 and args.verbose:
            for bad_time in bad_spectra_times:
                print(bad_time)
                print(f'[INFO] Spectrum at {bad_time} is invalid')

        if not insitu_files is None:
            ninsitu = len(insitu_files)
            if args.verbose:
                print(f'[INFO] Number of in situ files for the extract: {ninsitu}')
                print(f'[INFO] Creating MDB in situ extract: {ofile}')
            ihd.create_mdb_insitu_extract(extract_list[extract]['path'], ofile)
            idx = 0
            for insitu_file in insitu_files:
                b = ihd.set_data(insitu_file, idx, date_here, mo.get_sat_extracts_info())
                if b:
                    idx = idx + 1

            ihd.close_mdb()

            if os.path.exists(ofile):
                mdb_extract_files.append(ofile)
                if args.verbose:
                    print(f'[INFO] MDB extract file was created')

    nextract_files = len(mdb_extract_files)
    if args.verbose:
        print(f'[INFO] {nextract_files} were created/added')
        print(f'[INFO] Generating MDB extract files----------------------------------------------------------STOP')
    if nextract_files == 0:
        print(f'[INFO] Completed. No MDB file was created')
        return

    if args.verbose:
        print(f'[INFO] Generating final MDB file-------------------------------------------------------------START')
    path_mdb = mo.get_mdb_path()
    concatenate_nc_impl(mdb_extract_files, mo.path_out, path_mdb)
    if args.verbose:
        print(f'[INFO] Generating final MDB file-------------------------------------------------------------STOP')

    if len(bad_spectra_times) > 0 and args.verbose:
        for bad_time in bad_spectra_times:
            print(bad_time)
            print(f'[INFO] Spectrum at {bad_time} is invalid')


def create_mdb_tara_file(mo, extract_list):
    from INSITU_tara import INSITU_TARA
    it = INSITU_TARA(mo, args.verbose)

    date_list = list(extract_list.keys())
    date_list.sort()
    mdb_extract_files = []
    for date in date_list:
        print(f'[INFO] Working for date: {date} ****************************************************************')
        date_here = dt.strptime(date, '%Y%m%d')
        ##checking of the file dates
        file_metadata = it.get_file_metadata(date_here)
        if file_metadata is None:
            print(f'[WARNING] Metadata file for date: {date} is not available. Skipping extracts...')
            continue
        file_data = it.get_file_data(date_here)
        if file_data is None:
            print(f'[WARNING] Data file for date: {date} is not available. Skipping extracts...')
            continue
        extract_list_here = extract_list[date]
        mdb_extract_files = create_mdb_tara_file_impl(mo, extract_list_here, mdb_extract_files)

    nextract_files = len(mdb_extract_files)
    if args.verbose:
        print(f'[INFO] *******************')
        print(f'[INFO] {nextract_files} MDB extract files were created/added')
        print(f'[INFO] Generating MDB extract files----------------------------------------------------------STOP')
    if nextract_files == 0:
        print(f'[INFO] Completed. No MDB file was created')
        return

    path_mdb = mo.get_mdb_path()
    if args.verbose:
        print(f'[INFO] Generating final MDB file-------------------------------------------------------------START')
    concatenate_nc_impl(mdb_extract_files, mo.path_out, path_mdb)
    if args.verbose:
        print(f'[INFO] Generating final MDB file-------------------------------------------------------------STOP')


##all these extracts are associated to the same date, hence to the same metadata and data files.
def create_mdb_tara_file_impl(mo, extract_list, mdb_extract_files):
    from INSITU_tara import INSITU_TARA
    from netCDF4 import Dataset
    import numpy as np
    ins_sensor = 'HyperBOOST'

    it = INSITU_TARA(mo, args.verbose)
    metadata_started = False

    max_time_diff = mo.insitu_options['time_window']

    for extract in extract_list:
        if args.verbose:
            print(f'[INFO] Working with extract: {extract} *******************')
        extract_name = os.path.basename(extract_list[extract]['path'])
        ofile = mo.get_mdb_extract_path(extract_name, ins_sensor)
        if os.path.exists(ofile):
            dtal = Dataset(ofile)
            if 'insitu_original_bands' not in dtal.variables:
                dtal.close()
                print(f'[ERROR] Error in MDB extract file: {ofile}. Skipping...')
                continue
            insitu_bands = np.array(dtal.variables['insitu_original_bands'])
            vtal = insitu_bands[0]
            dtal.close()
            if vtal == -999.0:
                print(f'[ERROR] Error in MDB extract file. Skipping...')
                continue
            mdb_extract_files.append(ofile)
            print(f'[WARNING] MDB extract file already exits and added to the list')
            continue

        extract_time = extract_list[extract]['time']
        if not metadata_started:
            it.retrieve_metadata_from_file(extract_time)
            metadata_started = True

        create = it.create_mdb_insitu_extract(extract_list[extract]['path'], ofile, extract_list[extract]['time'],
                                              max_time_diff)

        if os.path.exists(ofile) and create:
            mdb_extract_files.append(ofile)
            if args.verbose:
                print(f'[INFO] MDB extract file was created and added to the list')

    return mdb_extract_files


def create_mdb_meda_file(mo, extract_list):
    ins_sensor = 'MEDA'
    mdb_extract_files = []
    from INSITU_meda import INSITU_MEDA
    im = INSITU_MEDA(mo, args.verbose)

    for extract in extract_list:
        if args.verbose:
            print(f'[INFO] Working with extract: {extract} *******************')
        # bad_spectra_times = {}
        extract_name = os.path.basename(extract_list[extract]['path'])
        ofile = mo.get_mdb_extract_path(extract_name, ins_sensor)
        if os.path.exists(ofile):
            from netCDF4 import Dataset
            import numpy as np
            dtal = Dataset(ofile)
            if 'insitu_original_bands' not in dtal.variables:
                dtal.close()
                print(f'[ERROR] Error in MDB extract file: {ofile}. Skipping...')
                continue
            insitu_bands = np.array(dtal.variables['insitu_original_bands'])
            vtal = insitu_bands[0]
            dtal.close()
            if vtal == -999.0:
                print(f'[ERROR] Error in MDB extract file. Skipping...')
                continue
            mdb_extract_files.append(ofile)
            print(f'[WARNING] MDB extract file already exits. Skipping...')
            continue
        create = im.create_mdb_insitu_extract(extract_list[extract]['path'], ofile, extract_list[extract]['time'])
        if not create:
            print(f'[WARNING] In situ file is not available for {extract}')
            continue
        b = im.set_data(extract_list[extract]['time'])
        im.close_mdb()

        if os.path.exists(ofile):
            if not b:
                os.remove(ofile)
            else:
                mdb_extract_files.append(ofile)
                if args.verbose:
                    print(f'[INFO] MDB extract file was created')

    nextract_files = len(mdb_extract_files)
    if args.verbose:
        print(f'[INFO] {nextract_files} were created/added')
        print(f'[INFO] Generating MDB extract files----------------------------------------------------------STOP')
    if nextract_files == 0:
        print(f'[INFO] Completed. No MDB file was created')
        return

    path_mdb = mo.get_mdb_path()
    if args.verbose:
        print(f'[INFO] Generating final MDB file-------------------------------------------------------------START')
    concatenate_nc_impl(mdb_extract_files, mo.path_out, path_mdb)
    if args.verbose:
        print(f'[INFO] Generating final MDB file-------------------------------------------------------------STOP')


def create_mdb_wisp_file(mo, extract_list, insitu_file):
    ins_sensor = 'WISP3'
    time_maxh = mo.insitu_options['time_window'] / 3600
    mdb_extract_files = []
    from INSITU_wisp3 import INSITU_WISP3
    iw = INSITU_WISP3(mo, insitu_file, args.verbose)

    for extract in extract_list:
        if args.verbose:
            print(f'[INFO] Working with extract: {extract} *******************')
        # bad_spectra_times = {}
        extract_name = os.path.basename(extract_list[extract]['path'])
        ofile = mo.get_mdb_extract_path(extract_name, ins_sensor)
        if os.path.exists(ofile):
            from netCDF4 import Dataset
            import numpy as np
            dtal = Dataset(ofile)
            if 'insitu_original_bands' not in dtal.variables:
                dtal.close()
                print(f'[ERROR] Error in MDB extract file: {ofile}. Skipping...')
                continue
            insitu_bands = np.array(dtal.variables['insitu_original_bands'])
            vtal = insitu_bands[0]
            dtal.close()
            if vtal == -999.0:
                print(f'[ERROR] Error in MDB extract file. Skipping...')
                continue
            mdb_extract_files.append(ofile)
            print(f'[WARNING] MDB extract file already exits. Skipping...')
            continue
        iw.create_mdb_insitu_extract(extract_list[extract]['path'], ofile)
        b = iw.set_data(extract_list[extract]['time'])
        iw.close_mdb()

        if os.path.exists(ofile):
            if not b:
                os.remove(ofile)
            else:
                mdb_extract_files.append(ofile)
                if args.verbose:
                    print(f'[INFO] MDB extract file was created')

    nextract_files = len(mdb_extract_files)
    if args.verbose:
        print(f'[INFO] {nextract_files} were created/added')
        print(f'[INFO] Generating MDB extract files----------------------------------------------------------STOP')
    if nextract_files == 0:
        print(f'[INFO] Completed. No MDB file was created')
        return

    path_mdb = mo.get_mdb_path()
    if args.verbose:
        print(f'[INFO] Generating final MDB file-------------------------------------------------------------START')
    concatenate_nc_impl(mdb_extract_files, mo.path_out, path_mdb)
    if args.verbose:
        print(f'[INFO] Generating final MDB file-------------------------------------------------------------STOP')


def create_mdb_single_csv_rrs(mo, insitu_file):
    import numpy as np
    first_rrs = mo.insitu_options['first_rrs']
    if first_rrs is None:
        print('Option first_rrs is not defined in [insitu_options] section in the configuration file')
        return
    last_rrs = mo.insitu_options['last_rrs']
    if last_rrs is None:
        print('Option last_rrs is not defined in [insitu_options] section in the configuration file')
        return
    try:
        df = pd.read_csv(insitu_file, sep=';')
    except:
        print(f'[ERROR] File {insitu_file} is not a valid csv file')
        return
    col_names = df.columns
    if first_rrs not in col_names:
        print(f'[ERROR] {first_rrs} is not a valid column name in the CSV file {insitu_file}')
        return
    if last_rrs not in col_names:
        print(f'[ERROR] {last_rrs} is not a valid column name in the CSV file {insitu_file}')
        return
    wl_values = []
    wl_cols = []
    be_added = False
    prefix = mo.insitu_options['rrs_prefix']
    suffix = mo.insitu_options['rrs_suffix']
    iini = 0 if prefix is None else len(prefix)
    ns = 0 if suffix is None else len(suffix)
    for col in col_names:
        if col == first_rrs:
            be_added = True
        if be_added:
            rrs_s = col[iini:len(col) - ns]
            rrs_s = rrs_s.replace('_', '.')
            try:
                wl_values.append(float(rrs_s))
                wl_cols.append(col)
            except:
                print(f'[ERROR] {rrs_s} is not a valid col name containing a rrs value')
                return
        if col == last_rrs:
            be_added = False
    print(wl_values)
    print(len(wl_values))
    if mo.insitu_options['n_insitu_bands'] != len(wl_values):
        mo.insitu_options['n_insitu_bands'] = len(wl_values)
        print(f'[WARNING] Number of bands in the configuration file was wrong, set to:{len(wl_values)}')

    for index, row in df.iterrows():
        rrs_array = np.array(df.loc[index, wl_cols]).astype(np.float64)
        rrs_array[np.isnan(rrs_array)] = -999.0
        if mo.insitu_options['fill_value'] is not None:
            rrs_array[rrs_array == mo.insitu_options['fill_value']] = -999.0
        print('----')


def create_mdb_single_csv_var(mo, insitu_file):
    import numpy as np
    from netCDF4 import Dataset
    try:
        df = pd.read_csv(insitu_file, sep=';')
    except:
        print(f'[ERROR] File {insitu_file} is not a valid csv file')
        return

    cols_to_check = ['col_extract', 'col_lat', 'col_lon', 'col_date']
    for col in cols_to_check:
        if not mo.insitu_options[col] in df.columns.tolist():
            print(
                f'[ERROR]{mo.insitu_options[col]} is not a valid column in the input csv file. Please review {col} in the [insitu_options] of the configuration file')
            return
    col_vars = mo.insitu_options['col_vars']
    if col_vars is None:
        print(
            f'[ERROR] No valid variable names were given. Please set col_vars key in the [insitu_options] of the configuration file')
        return
    for var in col_vars:
        if var not in df.columns.tolist():
            print(
                f'[ERROR] Variable {var} is not a valid variable in the csv file {insitu_file}. Please review col_vars key in the [insitu_options] of the configuration file ')
            return

    col_flags = mo.insitu_options['col_flags']
    flags_info = {}
    if col_flags is not None:
        for var in col_flags:
            if var not in df.columns.tolist():
                print(f'[ERROR] Variable {var} is not a valid variable in the csv file {insitu_file}. Please review col_flags key in the [insitu_options] of the configuration file ')
                return
            flag_meanings_list = np.unique(df[var]).tolist()
            flag_values = [pow(2,x) for x in range(len(flag_meanings_list))]
            flag_meanings = " ".join(flag_meanings_list)
            flags_info[var]={
                'flag_meanings_list': flag_meanings_list,
                'flag_values': flag_values,
                'attrs': {
                    'flag_values': flag_values,
                    'flag_meanings': flag_meanings
                }
            }


    # print(mo.insitu_options)
    col_extract = mo.insitu_options['col_extract']

    ###CHECKIN VARIABLES IN EXTRACTS
    variables_to_exclude = []
    if args.check_variables_extract:
        if args.verbose:
            print(f'[INFO] Checking extract variables...')
        all_variables = {}
        nfiles = 0
        for index, row in df.iterrows():
            name_extract = f'{row[col_extract]}'
            file_extract = name_extract
            if not os.path.exists(file_extract):
                file_extract = os.path.join(mo.satellite_path_source, name_extract)
            if not os.path.exists(file_extract):
                name_extract = f'extract_{row[col_extract]}'
                file_extract = os.path.join(mo.satellite_path_source, name_extract)
            if not os.path.exists(file_extract):
                continue
            dataset = Dataset(file_extract,'r')
            for var_name in dataset.variables:
                if var_name not in all_variables.keys():
                    all_variables[var_name]=1
                else:
                    all_variables[var_name] = all_variables[var_name] + 1
            dataset.close()
            nfiles = nfiles + 1
        for var_name in all_variables:
            if all_variables[var_name]<nfiles:
                variables_to_exclude.append(var_name)
        if len(variables_to_exclude)>0:
            print(f'[WARNING] The following {len(variables_to_exclude)} satellite variables will be excluded from the MDB:')
            for var_name in variables_to_exclude:
                print(f'-->{var_name}')
    ##END CHECKING VARIABLES


    extra_extracts = mo.insitu_options['col_extract_extra']
    if extra_extracts is None:
        extra_extracts = {}

    extra_extracts_info = {}
    for key in extra_extracts:
        vals = extra_extracts[key]
        path = vals[0]
        variables = [x.strip() for x in vals[1:]]
        if os.path.exists(path) and key in df.columns.tolist():
            extra_extracts_info[key] = {
                'path': vals[0],
                'variables': variables
            }
            extracts = np.array(df[key])
            for extract in extracts:
                file_ref = os.path.join(path, f'extract_{extract}')

                if os.path.exists(file_ref):
                    dref = Dataset(file_ref)
                    for variable in variables:
                        if variable in dref.variables:
                            extra_extracts_info[key][f'{variable}.dtype'] = dref.variables[variable].datatype
                            extra_extracts_info[key][f'{variable}.dims'] = dref.variables[variable].dimensions
                            extra_extracts_info[key][f'{variable}.atts'] = dref.variables[variable].__dict__
                    dref.close()
                    break
            ##checking
            all_var = True
            for variable in variables:
                if not f'{variable}.dtype' in extra_extracts_info[key]: all_var = False
                if not f'{variable}.dims' in extra_extracts_info[key]: all_var = False
                if not f'{variable}.atts' in extra_extracts_info[key]: all_var = False

            if not all_var:
                del extra_extracts_info[key]
                print(
                    f'[ERROR] col_extract_extra {key} is not correctly defined. Please review if the following variable names: ')
                print(f'[ERROR] {variables}')
                print(f'[ERROR] are available in extract files in {path}')
                return
        else:
            if not os.path.exists(path): print(
                f'[ERROR] Extract path: {path} for col_extract_extra {key} does not exist')
            if not key in df.columns.tolist(): print(
                f'[ERROR] col_extract_extra {key} is not available in the CSV file {insitu_file}')
            return

    extract_list = {}
    nrows = len(df.index)
    from INSITU_base import INSITUBASE
    ibase = INSITUBASE(mo)
    for index, row in df.iterrows():
        name_extract = f'{row[col_extract]}'
        file_extract = name_extract
        if not os.path.exists(file_extract):
            file_extract = os.path.join(mo.satellite_path_source,name_extract)
        if not os.path.exists(file_extract):
            name_extract = f'extract_{row[col_extract]}'
            file_extract = os.path.join(mo.satellite_path_source, name_extract)
        if not os.path.exists(file_extract):
            print(f'[WARNING] Extract file {row[col_extract]} is not available. Skipping row {index}')
            continue
        mdb_extract = mo.get_mdb_extract_path(row[col_extract], None)
        insitu_lat = get_float_value_from_row(index, row, mo.insitu_options['col_lat'])
        insitu_lon = get_float_value_from_row(index, row, mo.insitu_options['col_lon'])
        insitu_time = get_datetime_from_row(index, row, mo.insitu_options['col_date'], mo.insitu_options['format_date'],
                                            mo.insitu_options['col_time'], mo.insitu_options['format_time'],
                                            mo.insitu_options['insitu_time'])
        if insitu_time is None:
            print(f'[WARNING] Please review date and time formats in config file and/or data in the CSV file')
            break
        if (index % 100) == 0: print(
            f'[INFO] Row: {index}/{nrows} In situ time: {insitu_time.strftime("%Y-%m-%dT%H:%M:%S")} Latitude: {insitu_lat} Longitude: {insitu_lon}')

        if name_extract not in extract_list:
            if os.path.exists(mdb_extract):
                extract_list[name_extract] = {
                    'path': mdb_extract,
                    'lat': [insitu_lat],
                    'lon': [insitu_lon],
                    'time': [insitu_time.timestamp()]
                }
                for col in col_vars:
                    extract_list[name_extract][col] = [float(row[col])]
                if col_flags is not None:
                    for col in col_flags:
                        val = str(row[col])
                        index_val = flags_info[col]['flag_meanings_list'].index(val)
                        extract_list[name_extract][col] = int(flags_info[col]['flag_values'][index_val])
                if (index % 100) == 0: print(f'[INFO] MDB extract file already exists. Skipping...')
                continue
                #os.remove(mdb_extract)
            ibase.start_add_insitu_no_rrs(file_extract, mdb_extract,variables_to_exclude)
            ibase.add_shipborne_variables()
            ##adding variables and data for extra_extrats
            if len(extra_extracts_info) > 0:
                for extra in extra_extracts_info:
                    variables = extra_extracts_info[extra]['variables']
                    for variable in variables:
                        ibase.add_general_variable(variable, extra_extracts_info[extra][f'{variable}.dims'],
                                                   extra_extracts_info[extra][f'{variable}.dtype'],
                                                   extra_extracts_info[extra][f'{variable}.atts'])
                    path = extra_extracts_info[extra]['path']
                    file_extract = os.path.join(path,f'extract_{row[extra]}')
                    if os.path.exists(file_extract):
                        for variable in variables:
                            ibase.add_data_variable(file_extract,variable,ibase.get_name_satellite_variable(variable))
            for col in col_vars:
                name_var = f'insitu_{col}' if not col.startswith('insitu_') else col
                ats = mo.insitu_options[col] if col in mo.insitu_options else None
                ibase.add_insitu_variable(name_var, 'f4', ats)
            if col_flags is not None:
                for col in col_flags:
                    name_var = col
                    if not col.startswith('insitu_flag_'):
                        if col.startswith('insitu_'):
                            name_var = f'insitu_flag_{col.split("_")[1]}'
                        else:
                            name_var = f'insitu_flag_{col}'
                    ibase.add_insitu_variable(name_var, 'i4', flags_info[col]['attrs'])

            ibase.close()
            extract_list[name_extract] = {
                'path': mdb_extract,
                'lat': [insitu_lat],
                'lon': [insitu_lon],
                'time': [insitu_time.timestamp()]
            }
            for col in col_vars:
                extract_list[name_extract][col] = [float(row[col])]
            if col_flags is not None:
                for col in col_flags:
                    val = str(row[col])
                    index_val = flags_info[col]['flag_meanings_list'].index(val)

                    extract_list[name_extract][col] = [int(flags_info[col]['flag_values'][index_val])]

        else:
            extract_list[name_extract]['lat'].append(insitu_lat)
            extract_list[name_extract]['lon'].append(insitu_lon)
            extract_list[name_extract]['time'].append(insitu_time.timestamp())
            for col in col_vars:
                extract_list[name_extract][col].append(float(row[col]))
            if col_flags is not None:
                for col in col_flags:
                    val = str(row[col])
                    index_val = flags_info[col]['flag_meanings_list'].index(val)
                    print('AQUI DEBERIAMOS ESTAR O QUE NARICES, APPENDING...')
                    extract_list[name_extract][col].append(int(flags_info[col]['flag_values'][index_val]))



    mdb_paths_to_concatenate = []
    for extract in extract_list:
        print(f'[INFO] Setting data for extract {extract}')
        mdb_path = extract_list[extract]['path']
        time_array = np.array(extract_list[extract]['time'])
        indices_sort = np.argsort(time_array)
        time_array = time_array[indices_sort]
        lat_array = np.array(extract_list[extract]['lat'])[indices_sort]
        lon_array = np.array(extract_list[extract]['lon'])[indices_sort]

        max_time_diff = mo.insitu_options['time_window']
        dataset_w = Dataset(mdb_path, 'a')


        if max_time_diff > 0:
            time_diff = np.abs(time_array - dataset_w.variables['satellite_time'][0])
        elif max_time_diff <= 0:
            ##redundant check to be sure, satellite and in situ time on the same day
            max_time_diff = 0
            sat_time = dt.utcfromtimestamp(float(dataset_w.variables['satellite_time'][0])).strftime('%Y%m%d')
            time_diff = np.zeros(time_array.shape)
            for idx in range(len(time_array)):
                ins_time = dt.utcfromtimestamp(time_array[idx]).strftime('%Y%m%d')
                if sat_time != ins_time:
                    time_diff[idx] = 1

        time_array = time_array[time_diff <= max_time_diff]
        lat_array = lat_array[time_diff <= max_time_diff]
        lon_array = lon_array[time_diff <= max_time_diff]
        n_valid = len(time_array)

        if n_valid == 0:
            print(f'[WARNING] No valid data were found for extract {extract} Removing MDB path and skipping...')
            dataset_w.close()
            os.remove(mdb_path)
            continue

        dataset_w.variables['insitu_latitude'][0, 0:n_valid] = lat_array[:]
        dataset_w.variables['insitu_longitude'][0, 0:n_valid] = lon_array[:]
        dataset_w.variables['insitu_time'][0, 0:n_valid] = time_array[:]

        dataset_w.variables['time_difference'][0, 0:n_valid] = time_diff[:]


        for col in col_vars:
            name_var = f'insitu_{col}' if not col.startswith('insitu_') else col
            array = np.array(extract_list[extract][col])[indices_sort]
            if max_time_diff > (-1):
                array = array[time_diff <= max_time_diff]
            dataset_w.variables[name_var][0, 0:n_valid] = array[:]

        if col_flags is not None:
            for col in col_flags:
                name_var = col
                if not col.startswith('insitu_flag_'):
                    if col.startswith('insitu_'):
                        name_var = f'insitu_flag_{col.split("_")[1]}'
                    else:
                        name_var = f'insitu_flag_{col}'


                array = np.array(extract_list[extract][col])[indices_sort]
                if max_time_diff > (-1):
                    array = array[time_diff <= max_time_diff]
                dataset_w.variables[name_var][0, 0:n_valid] = array[:]

        dataset_w.close()
        if os.path.exists(mdb_path):
            mdb_paths_to_concatenate.append(mdb_path)

    if len(mdb_paths_to_concatenate) == 0:
        print(f'[WARNING] No MDB mini-files found for concatenation, no MDB file was created')
        return

    path_mdb = mo.get_mdb_path()
    if args.verbose:
        print(f'[INFO] Generating final MDB file-------------------------------------------------------------START')
    concatenate_nc_impl(mdb_paths_to_concatenate, mo.path_out, path_mdb)
    if args.verbose:
        print(f'[INFO] Generating final MDB file-------------------------------------------------------------STOP')


def get_datetime_from_row(index, row, col_date, format_date, col_time, format_time, insitu_time):
    date_row = str(row[col_date])
    format_row = format_date
    if col_time is not None:
        date_row = f'{date_row}{row[col_time]}'
        format_row = f'{format_row}{format_time}'
    try:
        date_here = dt.strptime(date_row, format_row).replace(tzinfo=pytz.UTC)
    except:
        print(f'[WARNING] Datetime {date_row} for row {index} could not be parsed using the format {format_row}')
        return None

    if insitu_time is not None:
        date_here_str = f'{date_here.strftime("%Y-%m-%d")}T{insitu_time}'
        try:
            date_here = dt.strptime(date_here_str, '%Y-%m-%dT%H:%M')
        except:
            print(
                f'[WARNING] insitu_time {insitu_time} for row {index} could not be parsed using the format %H:%M. Please review config file')
            return None
    date_here = date_here.replace(tzinfo=pytz.utc)
    return date_here


def get_float_value_from_row(index, row, col_name):
    try:
        value = float(row[col_name])
    except:
        value = -999.0
        print(f'[WARNING] {col_name} value {row[col_name]} for row {index} is a not a valid float value')
    return value


def concatenate_nc_impl(list_files, path_out, ncout_file):
    if len(list_files) == 0:
        print(f'[WARNING] No MDB sat extract files were found. Please review')
        return
    import subprocess
    nfiles_ref = 100
    if len(list_files) > nfiles_ref:
        if args.verbose:
            print(f'[INFO] Preparing contatenation of {len(list_files)} files...')
        list_files_tmp = []
        for icent in range(0, len(list_files), nfiles_ref):
            if args.verbose:
                print(f'[INFO] Concatening: {icent} / {len(list_files)}')
            indextmp = int(icent / nfiles_ref)
            list_files_here = list_files[icent:icent + nfiles_ref]
            if icent == 74:
                print(list_files_here)
            ncout_file_tmp = os.path.join(path_out, f'Temp_{indextmp}.nc')
            list_files_tmp.append(ncout_file_tmp)
            list_files_here.append(ncout_file_tmp)
            # concatenation
            cmd = [f"ncrcat -O -h"] + list_files_here
            cmd = " ".join(cmd)
            prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            out, err = prog.communicate()
            if err:
                print(f'[ERROR]{err}')

        list_files_tmp.append(ncout_file)
        cmd = [f"ncrcat -O -h"] + list_files_tmp
        cmd = " ".join(cmd)
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]{err}')

        if not args.nodelfiles:
            [os.remove(f) for f in list_files]

        for f in list_files_tmp:
            fname = f.split('/')[-1]
            if fname.startswith('Temp_'):
                os.remove(f)


    else:
        list_files.append(ncout_file)
        # concatenation
        cmd = [f"ncrcat -O -h"] + list_files
        cmd = " ".join(cmd)

        # print(f'CMD="{cmd}"')
        # os.system(cmd)
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]{err}')

        if not args.nodelfiles:
            [os.remove(f) for f in list_files[:-1]]
    if args.verbose:
        print(f'[INFO] Concatenated MDB file created: {ncout_file}')


# %%
if __name__ == '__main__':
    main()
