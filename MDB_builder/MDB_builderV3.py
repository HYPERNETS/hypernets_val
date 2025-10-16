import argparse,__init__,os,sys

import numpy as np
import pandas as pd

from MDB_optionsV3 import MDBBuilderOptions
from SATEXTRACTS_list import EXTRACT_LIST
import INSITU_base as ibase
code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
import COMMON.args_functions as arf

def main():
    print('[INFO] Creating MDB files!')

    if args.verbose:
        print(f'[INFO] Opening configuration file: {args.config_file}')
    mo = MDBBuilderOptions(args.config_file,args.verbose)
    if not mo.is_valid:
        return
    insitu_options = mo.get_insitu_options(args.insitu_type)
    if insitu_options is None:
        return
    insituBase = ibase.get_insitu_object(args.insitu_type, insitu_options, args.verbose)
    if insituBase is None:
        return
    start_date, end_date = arf.get_start_end_date_from_args(args)
    insituBase.start_date, insituBase.end_date = start_date,end_date
    if args.verbose:
        print(f'[INFO] Preparing in situ data...')
    if not insituBase.prepare_data():
        return

    if args.verbose:
        print(f'[INFO] Checking in situ Rrs and data variables...')
    if not insituBase.check_rrs_and_data_variables():
        return

    ##retrieving sat extract list
    if args.verbose:
        print(f'[INFO] Obtaining extract list----------------------------------------------------------------START')
    slist = EXTRACT_LIST(mo, args.insitu_type,args.sat_type, args.verbose)
    slist.prepare_extract_list(insituBase,mo.allow_partial_mdb())
    if slist.csv_files_by_date is None or len(slist.csv_files_by_date)==0:
        print(f'[ERROR] Extract list could not be retrieved. Exiting without completing the MDB file.')
        return
    else:
        if args.verbose:
            print(f'[INFO] Obtaining extract list----------------------------------------------------------------STOP')

    options = mo.get_mdb_options()

    extract_dir = os.path.join(options['output_dir'], 'MDBm')
    if not os.path.isdir(extract_dir):
        try:
            os.mkdir(extract_dir)
        except:
            print(
                f'[ERROR] MDBm folder could not be created in {options["output_dir"]}. Mini files could not be created. Please review permissions.')
            return
    ow = mo.overwrite_mini()
    dims_mdb = None

    for date_here in insituBase.date_list:
        if args.verbose:
            print(f'[INFO] Checking valid extracts for date: {date_here.strftime("%Y-%m-%d")}****************************')
        extracts,time_extracts,nall = slist.get_valid_extracts_date(insituBase,date_here)
        if extracts is None or time_extracts is None:
            print(f'[ERROR] Error retrieving extracts for date: {date_here.strftime("%Y-%m-%d")}, building is stopped. Launch again the extract only for that date before retrying the builder could solve the problem')
            return
        if len(time_extracts)>0:
            dims = insituBase.create_mini_mdb_files(options,extract_dir,extracts,time_extracts,ow)
            if dims is not None:
                dims_mdb = dims
    files_concat = []
    min_date  = None
    max_date = None
    for date_here in  insituBase.date_list:
        file_csv_mdbm = os.path.join(extract_dir,f'MDBm_{date_here.strftime("%Y%m%d")}.csv')
        if os.path.exists(file_csv_mdbm):
            df = pd.read_csv(file_csv_mdbm,sep=';')
            dims_min = df[['satellite_id','insitu_id','instrument_id','satellite_bands','insitu_bands','rows','columns']].min().to_numpy()
            dims_max = df[['satellite_id', 'insitu_id', 'instrument_id', 'satellite_bands', 'insitu_bands', 'rows','columns']].max().to_numpy()
            if np.equal(dims_mdb.all(),dims_min.all()) and np.equal(dims_mdb.all(),dims_max.all()):
                for name in df['name']:
                    file_name = os.path.join(extract_dir,name)
                    if os.path.exists(file_name):
                        files_concat.append(f'"{file_name}"')
                        if min_date is None and max_date is None:
                            min_date = date_here
                            max_date = date_here
                        else:
                            if date_here<min_date:
                                min_date = date_here
                            if date_here>max_date:
                                max_date = date_here
                    else:
                        print(f'[WARNING] {file_name} is not available for the concatenation...')
            os.remove(file_csv_mdbm)

    if len(files_concat)==0:
        print(f'[ERROR] No MDBm files were obtained for the concatenation. MDB file will no be created.')
        return

    file_mdb = os.path.join(options['output_dir'],mo.get_mdb_name(args.insitu_type,min_date,max_date))
    if args.verbose:
        print(f'[INFO] Number of MDBm files to be concatenated: {len(files_concat)}')
        print(f'[INFO] Output MDB file: {file_mdb}')

    concatenate_nc_impl(files_concat,file_mdb,mo.delete_mini())




def concatenate_nc_impl(list_files, ncout_file, delete_files):
    if len(list_files) == 0:
        print(f'[WARNING] No MDB sat extract files were found. Please review')
        return
    import subprocess
    nfiles_ref = 100
    path_out = os.path.dirname(ncout_file)
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
            ncout_file_tmp = f'"{ncout_file_tmp}"'
            list_files_tmp.append(ncout_file_tmp)
            list_files_here.append(ncout_file_tmp)
            # concatenation
            cmd = [f"ncrcat -O -h"] + list_files_here
            cmd = " ".join(cmd)
            prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            out, err = prog.communicate()
            if err:
                print(f'[ERROR]{err}')
        list_files_tmp.append(f'"{ncout_file}"')
        cmd = [f"ncrcat -O -h"] + list_files_tmp
        cmd = " ".join(cmd)
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]{err}')

        if delete_files:
            for f in list_files[:-1]:
                fhere = f.replace(f'"', '').strip()
                os.remove(fhere)

        for f in list_files_tmp:
            f = f.replace(f'"','').strip()
            fname = f.split('/')[-1]
            if fname.startswith('Temp_'):
                os.remove(f)


    else:
        list_files.append(f'"{ncout_file}"')
        # concatenation
        cmd = [f"ncrcat -O -h"] + list_files
        cmd = " ".join(cmd)

        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]{err}')

        if delete_files:
            for f in list_files[:-1]:
                fhere = f.replace(f'"','').strip()
                print('deleting,',fhere)
                os.remove(fhere)
    if args.verbose:
        print(f'[INFO] Concatenated MDB file created: {ncout_file}')



# %%
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Create Match-up DataBase files (MDB) files from satellite extracts and in situ L2 HYPERNETS files.")

    parser.add_argument('-c', "--config_file", help="Config File.",required=True)
    parser.add_argument('-o', "--output",
                        help="Output file. Required with --listdates or single concatenation")
    parser.add_argument('-edir', "--sat_extract_dir",
                        help="Input sat. extract dir. Optional for --listdates, required for single concatenation or metadata")
    parser.add_argument('-ifolder', "--insitu_folder", help="In situ data folder. Optional with --listdates")
    parser.add_argument('-site', "--sitename", help="Site name. Only required with --listdates")
    parser.add_argument('-ld', "--listdates",
                        help="Option to obtain a date list for a specific HYPERNETS site (-site option).",
                        action="store_true")
    parser.add_argument('-sd', "--start_date", help="Start date. Optional with --listdates (YYYY-mm-dd)")
    parser.add_argument('-ed', "--end_date", help="End date. Optional with --listdates (YYYY-mm-dd)")
    parser.add_argument('-nd', "--nodelfiles", help="Do not delete temp files.", action="store_true")
    parser.add_argument('-chvar', "--check_variables_extract", help="Check variables of extracts in single CSV files",
                        action="store_true")
    parser.add_argument('-insitu', "--insitu_type", help="In situ type",required=True)
    parser.add_argument('-sat', "--sat_type", help="Satellite type",required=True)
    parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")

    args = parser.parse_args()
    main()