import argparse,__init__,os,sys
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
    if not insituBase.prepare_data():
        return

    insituBase.check_rrs_and_data_variables()

    ##retrieving sat extract list
    if args.verbose:
        print(f'[INFO] Obtaining extract list----------------------------------------------------------------START')
    slist = EXTRACT_LIST(mo, args.insitu_type,args.sat_type, args.verbose)
    slist.prepare_extract_list(insituBase)
    if args.verbose:
        print(f'[INFO] Obtaining extract list----------------------------------------------------------------STOP')

    for date_here in insituBase.date_list:
        extracts,time_extracts = slist.get_valid_extracts_date(insituBase,date_here)
        if len(time_extracts)>0:
            insituBase.create_mini_mdb_files(mo,extracts,time_extracts)





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