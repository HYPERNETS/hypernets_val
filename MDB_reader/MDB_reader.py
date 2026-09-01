import argparse
import os.path

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

def main():
    mode = args.mode
    print(f'[INFO] Started MDBReader with mode: {args.mode}')
    if args.mode == 'GENERATEMU':

        input_path = args.input_path
        if not os.path.isfile(input_path):
            print(f'[ERROR] Input path {input_path} does not exist')
            return
        output_folder = os.path.dirname(input_path)
        if args.output and os.path.isdir(args.output):
            output_folder = args.output
        output_path = get_mdb_output_path(input_path, output_folder)
        reader = MDB_READER(input_path, True)
        if not reader.mfile.VALID:
            return
        if not reader.mfile.check_repeated() and not allow_repeated:
            print(f'To allow repeated ids (e.g. if your are working with shipborne data) use --allow_repeated option with GENERATE_MU')
            print(f'To remove repeated satellite ids before run GENERATE_MU,  you could run:')
            print(f'python MDB_readerV2.py -m REMOVEREP -i {input_path} -v')
            return
        if args.config_file and not os.path.exists(args.config_file):
            print(f'[ERROR] {args.config_file} is not available')
            return
if __name__ == '__main__':
    main()
