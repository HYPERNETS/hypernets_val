import os
import argparse
import configparser
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
parser.add_argument('-sd',"--start_date", help="Start date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-ed',"--end_date", help="End date. Optional with --listdates (YYYY-mm-dd)")
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


def main():
    print('[INFO] Creating MDB files!')

    # Option to a list dates
    if args.sitename and args.output and args.listdates:
        ihd = INSITU_HYPERNETS_DAY(None, None, args.verbose)
        sat_extract_dir = None
        if args.sat_extract_dir:
            sat_extract_dir = args.sat_extract_dir
        start_date = None
        end_date = None
        if args.start_date and args.end_date:
            try:
                start_date = dt.strptime(args.start_date,'%Y-%m-%d')
                end_date = dt.strptime(args.end_date,'%Y-%m-%d')
            except:
                print(f'[ERROR] Parameters start date: {args.start_date} and/or end date: {args.end_date} are not in the correct format (YYYY-mm-dd)')
                return

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
    #     print(e,extract_list[e]['time'])
    if args.verbose:
        print(f'[INFO] Obtaining extract list----------------------------------------------------------------STOP')

    ##checking in situ files
    if args.verbose:
        print(f'[INFO] Generating MDB extract files----------------------------------------------------------START')
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
            print(f'[INFO] Working with extract file: {extract} *******************')
        bad_spectra_times = {}
        ofile = mo.get_mdb_extract_path(extract, ins_sensor)
        if os.path.exists(ofile):
            from netCDF4 import Dataset
            import numpy as np
            dtal = Dataset(ofile)
            insitu_bands = np.array(dtal.variables['insitu_original_bands'])
            vtal = insitu_bands[0]
            dtal.close()
            if vtal==-999.0:
                print(f'[ERROR] Error in MDB extract file. Skipping...')
                continue
            mdb_extract_files.append(ofile)
            print(f'[WARNING] MDB extract file already exits. Skipping...')
            continue
        #date_here_str = extract_list[extract]['time']
        #date_here = dt.strptime(date_here_str, '%Y%m%dT%H%M%S')
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
            if icent == 290:
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
