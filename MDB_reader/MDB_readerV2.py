import os.path
import sys
from MDBFile import MDBFile

import argparse
from MDB_builder.INSITU_base import INSITUBASE

parser = argparse.ArgumentParser(
    description="Match-ups extraction from MDB files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-i', "--input_path", help="Input MDB path")
parser.add_argument('-o', "--output", help="Path to output")

# parser.add_argument('-sd', "--startdate", help="The Start Date - format YYYY-MM-DD ")
# parser.add_argument('-ed', "--enddate", help="The End Date - format YYYY-MM-DD ")
# parser.add_argument('-site', "--sitename", help="Site name.", choices=['VEIT', 'BEFR', 'BSBE'])
# parser.add_argument('-ins', "--insitu", help="Satellite sensor name.", choices=['PANTHYR', 'HYPERNETS'])  # ,'HYPSTAR'])
# parser.add_argument('-pi', "--path_to_ins", help="Path to in situ sources.")
# parser.add_argument('-sat', "--satellite", help="Satellite sensor name.", choices=['OLCI', 'MSI'])

# parser.add_argument('-ps', "--path_to_sat", help="Path to satellite extracts.")
# parser.add_argument('-o', "--output", help="Path to output")
# parser.add_argument('-res', "--resolution", help="Resolution OL_2: WRR or WFR (for OLCI)")
# parser.add_argument('-nl', "--nolist", help="Do not create satellite and in situ lists.", action="store_true")
# parser.add_argument('-nd', "--nodelfiles", help="Do not delete temp files.", action="store_true")

args = parser.parse_args()


class MDB_READER():

    def __init__(self, path_mdb, start_mdb):
        if path_mdb is not None:
            self.path_mdb = path_mdb
            if start_mdb:
                self.mfile = MDBFile(path_mdb)

    def create_mdb_results_file(self, fout):
        if not self.mfile.VALID:
            return

        if self.mfile.df_validation is None:
            self.mfile.prepare_df_validation()

        ndata = len(self.mfile.df_validation.index)
        ibase = INSITUBASE(None)
        new_MDB = ibase.copy_nc(self.mfile.file_path, fout)
        new_MDB.createDimension('mu_id',ndata)

        import numpy as np
        import numpy.ma as ma
        new_variables = {
            'mu_wavelength':{
                'namedf':'Wavelenght',
                'fillvalue': None,
                'type': 'f4'
            },
            'mu_satellite_id':{
                'namedf': 'Index_MU',
                'fillvalue': None,
                'type': 'i4'
            },
            'mu_ins_rrs':{
                'namedf': 'Ins_Rrs',
                'fillvalue': -999,
                'type': 'f4'
            },
            'mu_sat_rrs':{
                'namedf': 'Sat_Rrs',
                'fillvalue': -999,
                'type': 'f4'
            }
        }
        for new_var_name in new_variables:
            #type = new_variables[new_var_name]['type']
            print(new_var_name)
            new_var = new_MDB.createVariable(new_var_name,new_variables[new_var_name]['type'],('mu_id',), zlib=True,complevel=6)
            print('obtener la array...')
            array = np.array(self.mfile.df_validation.loc[:, new_variables[new_var_name]['namedf']])

            fillValue = new_variables[new_var_name]['fillvalue']
            if  fillValue is not None:
                array = ma.masked_array(array, mask= array == fillValue)
            print(array.shape,type(array))
            new_var[:] = [array[:]]

        new_variables_sat_mu = {
            'mu_sat_time': {
                'namedf': 'Sat_Time',
                'fillvalue': None,
                'type': 'f4'
            },
            'mu_ins_time': {
                'namedf': 'Ins_Time',
                'fillvalue': None,
                'type': 'f4'
            },
            'mu_time_diff': {
                'namedf': 'Time_Diff',
                'fillvalue': None,
                'type': 'f4'
            },
            'mu_valid':{
                'namedf': 'Valid',
                'fillvalue': None,
                'type': 'i8'
            }
        }

        if self.mfile.df_mu is None:
            self.mfile.prepare_df_mu()

        new_MDB.close()

        # mu_wavelength = new_MDB.createVariable('mu_wavelength', 'f4', ('mu_id',), zlib=True,complevel=6)
        # array = np.array(self.mfile.df_validation.loc[:,'Wavelenght'])
        # array = ma.masked_array(array,mask=array==-999)
        # mu_wavelength[:] = [[array]]
        #
        # mu_satellite_id = new_MDB.createVariable('mu_satellite_id','u4', ('mu_id',), zlib=True,complevel=6)
        # mu_ins_rrs = new_MDB.createVariable('mu_ins_rrs', 'f4', ('mu_id',), zlib=True, complevel=6)
        # mu_sat_rrs = new_MDB.createVariable('mu_sat_rrs', 'f4', ('mu_id',), zlib=True, complevel=6)


    def set_defaults_olci_wfr_hypstar(self):
        wllist = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.8, 681.3, 708.8, 753.8]
        self.mfile.set_wl_ref(wllist)

        ##Satellite quality control
        self.mfile.qc_sat.wl_ref = wllist
        self.mfile.qc_sat.set_eumetsat_defaults(3)

        # In situ quality control
        self.mfile.qc_insitu.set_wllist_using_wlref(wllist)
        self.mfile.qc_insitu.set_thershold(0, None, 100, 1000)


def get_mdb_output_path(input_path, output_folder):
    input_name = os.path.basename(input_path)
    if input_name.startswith('MDB_'):
        output_name = f'MDBr_{input_name[3:]}'
    output_path = os.path.join(output_folder, output_name)
    return output_path


def main():
    print('Started MDBReader!')

    ##WORKING WITH DEFAULT OPTIONS FROM A INPUT MDB FILE
    if args.input_path:
        input_path = args.input_path
        if not os.path.exists(input_path):
            print(f'[ERROR] Input path {input_path} does not exist')
            return
        output_folder = os.path.dirname(input_path)
        if args.output and os.path.isdir(args.output):
            output_folder = args.output
        output_path = get_mdb_output_path(input_path, output_folder)
        reader = MDB_READER(input_path, True)
        reader.set_defaults_olci_wfr_hypstar()
        reader.create_mdb_results_file(output_path)


if __name__ == '__main__':
    main()
