import os
import argparse
import configparser

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
    date_here = dt(2023,4,3)
    im.retrieve_metadata_from_file(date_here)
    file_extract = '/mnt/c/DATA_LUIS/TARA_TEST/extracts/S3B_OL_2_WFR____20230403T111113_20230403T111413_20230404T221050_0179_078_037_2160_MAR_O_NT_003_SEN3_extract_549_4638.nc'
    from netCDF4 import Dataset
    dataset = Dataset(file_extract)
    lat_array = dataset.variables['satellite_latitude'][0,:,:]
    lon_array = dataset.variables['satellite_longitude'][0,:,:]
    import numpy as np
    index_time = float(np.array(dataset.variables['satellite_time'][0]))
    sat_time = dt.fromtimestamp(index_time).replace(tzinfo=pytz.utc)
    dataset.close()
    max_time_diff = 180 * 60
    im.check_match_up_conditions(sat_time,lat_array,lon_array,max_time_diff)


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
            idm = INSITU_MEDA(None,args.verbose)
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
    if args.verbose:
        print(f'[INFO] Checking available dates--------------------------------------------------------------START')
    mo.get_dates()
    if args.verbose:
        print(f'[INFO] Start date for MDB_builder:{mo.start_date}')
        print(f'[INFO] End date for MDB_builder: {mo.end_date}')
        print(f'[INFO] Checking available dates--------------------------------------------------------------STOP')

    ##retrieving sat extract list
    if args.verbose:
        print(f'[INFO] Obtaining extract list----------------------------------------------------------------START')
    slist = SAT_EXTRACTS_LIST(mo, args.verbose)
    extract_list = slist.get_list_as_dict()
    # for e in extract_list:
    #     print(e, extract_list[e])
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
        create_mdb_meda_file(mo,extract_list)
        return

    if mo.insitu_type == 'TARA':
        create_mdb_tara_file(mo,extract_list)
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
    from netCDF4 import Dataset
    import numpy as np
    ins_sensor = 'HyperBOOST'
    mdb_extract_files = []
    it = INSITU_TARA(mo,args.verbose)
    extract_time_prev = None
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
        if extract_time!=extract_time_prev:
            it.retrieve_metadata_from_file(extract_time)
            extract_time_prev = extract_time
        create = it.create_mdb_insitu_extract(extract_list[extract]['path'], ofile, extract_list[extract]['time'],max_time_diff)

        if os.path.exists(ofile) and create:
            mdb_extract_files.append(ofile)
            if args.verbose:
                print(f'[INFO] MDB extract file was created and added to the list')

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

def create_mdb_meda_file(mo, extract_list):
    ins_sensor = 'MEDA'
    mdb_extract_files = []
    from INSITU_meda import INSITU_MEDA
    im = INSITU_MEDA(mo,args.verbose)

    for extract in extract_list:
        if args.verbose:
            print(f'[INFO] Working with extract: {extract} *******************')
        #bad_spectra_times = {}
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
        create = im.create_mdb_insitu_extract(extract_list[extract]['path'],ofile,extract_list[extract]['time'])
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
    iw = INSITU_WISP3(mo,insitu_file,args.verbose)

    for extract in extract_list:
        if args.verbose:
            print(f'[INFO] Working with extract: {extract} *******************')
        #bad_spectra_times = {}
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
        iw.create_mdb_insitu_extract(extract_list[extract]['path'],ofile)
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
        # print(cmd)
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
