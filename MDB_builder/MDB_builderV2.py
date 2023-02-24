import os
import argparse
import configparser
from MDB_builder_options import MDBBuilderOptions
from SATEXTRACTS_list import SAT_EXTRACTS_LIST
from INSITU_hypernets import INSITU_HYPERNETS_DAY
from datetime import datetime as dt

parser = argparse.ArgumentParser(
    description="Create Match-up DataBase files (MDB) files from satellite extracts and in situ L2 files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
# parser.add_argument('-sd', "--startdate", help="The Start Date - format YYYY-MM-DD ")
# parser.add_argument('-ed', "--enddate", help="The End Date - format YYYY-MM-DD ")
parser.add_argument('-site', "--sitename", help="Site name. Only with --listdates")
# parser.add_argument('-ins', "--insitu", help="Satellite sensor name.", choices=['PANTHYR', 'HYPERNETS'])  # ,'HYPSTAR'])
# parser.add_argument('-pi', "--path_to_ins", help="Path to in situ sources.")
# parser.add_argument('-sat', "--satellite", help="Satellite sensor name.", choices=['OLCI', 'MSI'])
parser.add_argument('-c', "--config_file", help="Config File.", required=True)
# parser.add_argument('-ps', "--path_to_sat", help="Path to satellite extracts.")
parser.add_argument('-o', "--output", help="Output file. Only with --listdates")
# parser.add_argument('-res', "--resolution", help="Resolution OL_2: WRR or WFR (for OLCI)")
parser.add_argument('-ld', "--listdates", help="Option to obtain a date list for a specific HYPERNETS site (-site option).", action="store_true")
parser.add_argument('-nd', "--nodelfiles", help="Do not delete temp files.", action="store_true")

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

    #Option to a list dates
    if args.sitename and args.output and args.listdates:
        ihd = INSITU_HYPERNETS_DAY(None, args.verbose)
        ihd.save_list_dates_to_file(args.output,args.sitename,None,None)
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
    if args.verbose:
        print(f'[INFO] Obtaining extract list----------------------------------------------------------------STOP')

    ##checking in situ files
    if args.verbose:
        print(f'[INFO] Generating MDB extract files----------------------------------------------------------START')
    ihd = INSITU_HYPERNETS_DAY(mo, args.verbose)
    ins_sensor = 'HYPSTAR'
    mdb_extract_files = []
    for extract in extract_list:
        if args.verbose:
            print(f'[INFO] Working with extract file: {extract} *******************')
        ofile = mo.get_mdb_extract_path(extract, ins_sensor)
        if os.path.exists(ofile):
            mdb_extract_files.append(ofile)
            print(f'[WARNING] MDB extract file already exits. Skipping...')
            continue
        date_here_str = extract_list[extract]['time']
        date_here = dt.strptime(date_here_str, '%Y%m%dT%H%M%S')
        insitu_files = ihd.get_insitu_files(date_here)
        if insitu_files is None and mo.insitu_options['apply_rsync'] and ihd.check_ssh():
            ihd.get_files_day_ssh(date_here, True)
            insitu_files = ihd.get_insitu_files(date_here)

        bad_spectra_times = {}
        print(mo.insitu_options)
        if mo.insitu_options['insitu_bad_spectra_file_list'] is not None:
            prefix = mo.insitu_options['bad_spectra_prefix']
            time_format = mo.insitu_options['bad_spectra_format_time']
            f1 = open(mo.insitu_options['insitu_bad_spectra_file_list'])
            for line in f1:
                if prefix is not None:
                    if not line.strip().startswith(prefix):
                        continue
                    datestr = line.replace(prefix,'')
                else:
                    datestr = line.strip()

                date_py = dt.strptime(datestr,time_format)
                date_py_str = date_py.strftime('%Y%m%d%H%M')
                bad_spectra_times[date_py_str] = 1
            f1.close()

        if len(bad_spectra_times)>0 and args.verbose:
            for bad_time in bad_spectra_times:
                print(f'[INFO] Spectrum at {bad_time} is invalid')

        if not insitu_files is None:
            ninsitu = len(insitu_files)
            if args.verbose:
                print(f'[INFO] Number of in situ files for the extract: {ninsitu}')
            ihd.create_mdb_insitu_extract(extract_list[extract]['path'], ofile)
            for idx in range(ninsitu):
                insitu_file = insitu_files[idx]
                # print(insitu_file)
                ihd.set_data(insitu_file, idx, date_here, mo.get_sat_extracts_info())

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
