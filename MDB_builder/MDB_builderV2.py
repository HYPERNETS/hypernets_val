import os,sys
import argparse
import configparser
from MDB_builder_options import MDBBuilderOptions

parser = argparse.ArgumentParser(
    description="Create Match-up DataBase files (MDB) files from satellite extracts and in situ L2 files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
# parser.add_argument('-sd', "--startdate", help="The Start Date - format YYYY-MM-DD ")
# parser.add_argument('-ed', "--enddate", help="The End Date - format YYYY-MM-DD ")
# parser.add_argument('-site', "--sitename", help="Site name.", choices=['VEIT', 'BEFR', 'BSBE'])
# parser.add_argument('-ins', "--insitu", help="Satellite sensor name.", choices=['PANTHYR', 'HYPERNETS'])  # ,'HYPSTAR'])
# parser.add_argument('-pi', "--path_to_ins", help="Path to in situ sources.")
# parser.add_argument('-sat', "--satellite", help="Satellite sensor name.", choices=['OLCI', 'MSI'])
parser.add_argument('-c', "--config_file", help="Config File.", required=True)
# parser.add_argument('-ps', "--path_to_sat", help="Path to satellite extracts.")
# parser.add_argument('-o', "--output", help="Path to output")
# parser.add_argument('-res', "--resolution", help="Resolution OL_2: WRR or WFR (for OLCI)")
parser.add_argument('-nl', "--nolist", help="Do not create satellite and in situ lists.", action="store_true")
parser.add_argument('-nd', "--nodelfiles", help="Do not delete temp files.", action="store_true")

args = parser.parse_args()

code_home = os.path.abspath('../')
sys.path.append(code_home)

# import BRDF.brdf_olci as brdf
import COMMON.common_functions as cfs




def main():
    print('[INFO] Creating MDB files!')
    if os.path.isfile(args.config_file):
        options = configparser.ConfigParser()
        options.read(args.config_file)
        mo = MDBBuilderOptions(options,args.verbose)
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
    mo.get_param_sat_extracts(options)


# %%
if __name__ == '__main__':
    main()
