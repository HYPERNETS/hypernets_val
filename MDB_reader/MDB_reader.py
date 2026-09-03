import __init__,argparse, warnings
import COMMON.args_functions as afs
import COMMON.common_functions as cfs
from MDBFile import MDBFile
from QC_OPTIONS import QC_OPTIONS

# import datetime
# from datetime import datetime as dt
# from datetime import timedelta
# from datetime import timezone
# import matplotlib.pyplot as plt
# import numpy as np
# import pytz
# import os.path
# import argparse
# import warnings
#
# import pandas as pd

warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', RuntimeWarning)



# try:
#     from MDB_reader.MDBFile import MDBFile
# except:
#     from MDBFile import MDBFile
# from MDB_builder.INSITU_base import INSITUBASE


parser = argparse.ArgumentParser(description="Match-ups extraction from MDB files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument("-m", "--mode", help="Mode",choices=["GENERATEMU", "GENERATEMU_S"],required=True)
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-i', "--input_path", help="Input MDB path")
parser.add_argument('-arep', "--allow_repeated", help="Set to allow multiple match-ups on the same date, e.g., for transects or multiple AC applied to the same satellite extracts",action="store_true")
#parser.add_argument('-rmdb', "--reduce_mdbr", help="MDBr should be reduced to only one insitu_id",action="store_true")
#parser.add_argument('-version',"--version_plot", help="Plot version", default='V3', choices=['V2','V3'])
args = parser.parse_args()


def main():
    print(f'[INFO] Started MDBReader with mode: {args.mode}')
    if args.mode == 'GENERATEMU':
        required_args = {'input_path':{'type':'input_file_nc'},'config_file':{'type':'input_file_ini'}}
        args_d = afs.get_args_as_dict(args,required_args,False)
        if args_d is None:
            return
        output_path = cfs.get_mdb_output_path(args_d['input_path'])
        if output_path is None:
            return
        if args.verbose:
            print(f'[INFO] MDB input path: {args_d["input_path"]}')
            print(f'[INFO] Configuration file: {args_d["config_file"]}')
            print(f'[INFO] MDB output path: {output_path}')
        mfile = MDBFile(args_d['input_path'])

        if not mfile.VALID:
            return

        if not mfile.check_repeated() and not args.allow_repeated:
            print(f'To allow repeated ids (e.g. if your are working with shipborne data) use --allow_repeated option with GENERATE_MU')
            print(f'To remove repeated satellite ids before run GENERATE_MU,  you could run:')
            print(f'python MDB_readerV2.py -m REMOVEREP -i {args.input_path} -v')
            return
        qc_options = QC_OPTIONS(args_d['config_file'],args.verbose)
        if not qc_options.is_valid:
            return
        qc_sat = qc_options.get_qc_sat(mfile.nc)


if __name__ == '__main__':
    main()
