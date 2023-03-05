import datetime
import os.path
import sys
import argparse
import warnings

warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', RuntimeWarning)

from MDBFile import MDBFile
from MDB_builder.INSITU_base import INSITUBASE

parser = argparse.ArgumentParser(
    description="Match-ups extraction from MDB files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument("-m", "--mode", help="Mode", choices=["GENERATEMU", "CONCATENATE", "REMOVEREP", "PLOT", "TEST"],
                    required=True)
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

        ibase = INSITUBASE(None)
        new_MDB = ibase.copy_nc(self.mfile.file_path, fout)
        new_MDB.createDimension('mu_id', None)

        import numpy as np
        import numpy.ma as ma
        new_variables = {
            'mu_wavelength': {
                'namedf': 'Wavelenght',
                'fillvalue': None,
                'type': 'f4'
            },
            'mu_satellite_id': {
                'namedf': 'Index_MU',
                'fillvalue': None,
                'type': 'i4'
            },
            'mu_ins_rrs': {
                'namedf': 'Ins_Rrs',
                'fillvalue': -999,
                'type': 'f4'
            },
            'mu_sat_rrs': {
                'namedf': 'Sat_Rrs',
                'fillvalue': -999,
                'type': 'f4'
            }
        }
        for new_var_name in new_variables:
            new_var = new_MDB.createVariable(new_var_name, new_variables[new_var_name]['type'], ('mu_id',), zlib=True,
                                             complevel=6)
            array = np.array(self.mfile.df_validation.loc[:, new_variables[new_var_name]['namedf']])
            fillValue = new_variables[new_var_name]['fillvalue']
            if fillValue is not None:
                array = ma.masked_array(array, mask=array == fillValue)
            # unlimited dimension,
            new_var[:] = array[:]

        new_variables_sat_mu = {
            'mu_sat_time': {
                'namedf': 'Sat_Time',
                'fillvalue': -999.0,
                'type': 'f4'
            },
            'mu_ins_time': {
                'namedf': 'Ins_Time',
                'fillvalue': -999.0,
                'type': 'f4'
            },
            'mu_time_diff': {
                'namedf': 'Time_Diff',
                'fillvalue': -999.0,
                'type': 'f4'
            },
            'mu_valid': {
                'namedf': 'mu_valid',
                'fillvalue': -999,
                'type': 'i8'
            }
        }

        if self.mfile.df_mu is None:
            self.mfile.prepare_df_mu()

        # print(self.mfile.df_mu)

        for new_var_name in new_variables_sat_mu:
            new_var = new_MDB.createVariable(new_var_name, new_variables_sat_mu[new_var_name]['type'],
                                             ('satellite_id',),
                                             zlib=True,
                                             complevel=6)
            array = np.array(self.mfile.df_mu.loc[:, new_variables_sat_mu[new_var_name]['namedf']])
            if new_var_name == 'mu_sat_time' or new_var_name == 'mu_ins_time':
                array_t = []
                for val in array:
                    try:
                        valdt = datetime.datetime.strptime(val, '%Y-%m-%d %H:%M')
                        array_t.append(valdt.timestamp())
                    except:
                        array_t.append(-999.0)
                array = np.array(array_t)
            # print(new_var_name, '->', array.shape)
            fillValue = new_variables_sat_mu[new_var_name]['fillvalue']
            if fillValue is not None:
                array = ma.masked_array(array, mask=array == fillValue)
            for idx in range(self.mfile.n_mu_total):
                new_var[idx] = [array[idx]]

        new_MDB.close()

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


def get_flag_lists(input_path, ats_in, flag_bands):
    flag_lists = {}
    for name in os.listdir(input_path):
        if not name.endswith('.nc'):
            continue
        file_in = os.path.join(input_path, name)
        reader = MDB_READER(file_in, True)
        for idx in range(len(ats_in)):
            name_band = flag_bands[idx]
            at_names = ats_in[idx]
            at_value = get_at_value(reader, at_names)
            if name_band not in flag_lists:
                flag_lists[name_band] = {
                    'input_ats': at_names,
                    'flag_values': [at_value]
                }
            else:
                if at_value not in flag_lists[name_band]['flag_values']:
                    flag_lists[name_band]['flag_values'].append(at_value)

        reader.mfile.close()
    return flag_lists


def get_at_value(reader, at_names):
    if type(at_names) == list:
        at_value = []
        for at in at_names:
            if at in reader.mfile.info:  # list(reader.mfile.nc.ncattrs()):
                at_value.append(reader.mfile.info[at])
        at_value = ''.join(at_value)
    else:
        if at_names in reader.mfile.info:  # list(reader.mfile.nc.ncattrs()):
            at_value = reader.mfile.info[at_names]  # nc.getncattr(at_names)

    return at_value


def creating_copy_with_flag_bands(reader, file_out, flag_lists, satellite_id_ref):
    # reader = MDB_READER('', True)
    from netCDF4 import Dataset
    import numpy as np
    ncout = Dataset(file_out, 'w', format='NETCDF4')

    # copy global attributes all at once via dictionary
    ncout.setncatts(reader.mfile.nc.__dict__)

    # copy dimensions
    for name, dimension in reader.mfile.nc.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    # copy variables
    for name, variable in reader.mfile.nc.variables.items():
        fill_value = None
        if '_FillValue' in list(reader.mfile.nc.ncattrs()):
            fill_value = variable._FillValue
        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(reader.mfile.nc[name].__dict__)
        # copy variable data
        if name == 'mu_satellite_id' and satellite_id_ref > 0:
            # array = reader.mfile.nc[name][:] + satellite_id_ref
            ncout[name][:] = reader.mfile.nc[name][:] + satellite_id_ref
        else:
            ncout[name][:] = reader.mfile.nc[name][:]

    ##creating flag variables
    for name_flag in flag_lists:
        var = ncout.createVariable(name_flag, 'i4', ('satellite_id',), zlib=True, shuffle=True, complevel=6)
        flag_list = flag_lists[name_flag]['flag_values']
        flag_meanings = ' '.join(flag_list)
        flag_values = np.power(2, np.arange(len(flag_list)))
        var.flag_meanings = flag_meanings
        var.flag_values = flag_values
        at_in = flag_lists[name_flag]['input_ats']
        value = get_at_value(reader, at_in)
        id = flag_list.index(value)
        var[:] = flag_values[id]

    ncout.close()
    reader.mfile.close()
    return True


def concatenate_nc_impl(list_files, path_out, ncout_file):
    if len(list_files) == 0:
        print(f'[WARNING] No files were found. Please review')
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
            ncout_file_tmp = os.path.join(path_out, f'TempList_{indextmp}.nc')
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

        [os.remove(f) for f in list_files[:-1]]
    if args.verbose:
        print(f'[INFO] Concatenated MDB file created: {ncout_file}')


def main():
    mode = args.mode
    print(f'Started MDBReader with mode: {mode}')

    if args.mode == 'TEST':
        fconfig = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/config_qc.ini'
        import configparser
        from QC_OPTIONS import QC_OPTIONS
        # from QC_SAT import  QC_SAT
        options = configparser.ConfigParser()
        options.read(fconfig)
        qco = QC_OPTIONS(options)
        qco.getting_qcsat(None)

    ##REMOVE REPETEADED SATELLITE ID FROM SINGLE MDF FILES
    if args.input_path and mode == 'REMOVEREP':
        print('Running option: Remove repeated')
        input_path = args.input_path
        if not os.path.exists(input_path):
            print(f'[ERROR] Input path {input_path} does not exist')
            return
        reader = MDB_READER(input_path, True)
        if not reader.mfile.check_repeated():
            reader.mfile.remove_repeated()
        else:
            print(f'[INFO] No repeated ids were found')
        return

    ##CREATING MDB READER FILES WORKING WITH DEFAULT OPTIONS FROM A INPUT MDB FILE
    if args.input_path and mode == 'GENERATEMU':
        input_path = args.input_path
        if not os.path.exists(input_path):
            print(f'[ERROR] Input path {input_path} does not exist')
            return
        output_folder = os.path.dirname(input_path)
        if args.output and os.path.isdir(args.output):
            output_folder = args.output
        output_path = get_mdb_output_path(input_path, output_folder)
        reader = MDB_READER(input_path, True)
        if not reader.mfile.VALID:
            return
        if not reader.mfile.check_repeated():
            print('[ERROR] To remove repeated satellite ids you could run:')
            print(f'python MDB_readerV2.py -m REMOVEREP -i {input_path} -v')
            return
        if args.config_file and os.path.exists(args.config_file):
            if args.verbose:
                print(f'[INFO] Using file: {args.config_file} to set quality control options...')
            import configparser
            from QC_OPTIONS import QC_OPTIONS
            options = configparser.ConfigParser()
            options.read(args.config_file)
            qco = QC_OPTIONS(options)
            wllist = qco.get_wllist()
            reader.mfile.set_wl_ref(wllist)
            reader.mfile.qc_sat.ncdataset = reader.mfile.nc
            reader.mfile.qc_sat = qco.get_qcsat(reader.mfile.qc_sat,reader.mfile.nc)
            reader.mfile.qc_insitu.ncdataset = reader.mfile.nc
            reader.mfile.qc_insitu = qco.get_qc_insitu(reader.mfile.qc_insitu)

        else:
            reader.set_defaults_olci_wfr_hypstar()
        reader.create_mdb_results_file(output_path)
        return

    ##CONCATENATING MDB READER FILES INCLUDING FLAGS BANDS FOR SATELLITE+PLATFORM, SENSOR, AC, SITE
    if args.input_path and mode == 'CONCATENATE':
        input_path = args.input_path
        if not os.path.exists(input_path):
            print(f'[ERROR] {input_path} does not exist')
        if not os.path.isdir(input_path):
            print(f'[ERROR] {input_path} should be a directory')
            return
        if args.output:
            if os.path.isdir(args.output):
                output_folder = args.output
                ncout_file = os.path.join(output_folder, 'MDBrc.nc')
            else:
                ncout_file = args.output
        else:
            ncout_file = os.path.join(input_path, 'MDBrc.nc')

        ats_in = [['satellite', 'platform'], 'sensor', 'satellite_aco_processor', 'insitu_site_name']
        flag_bands = ['flag_satellite', 'flag_sensor', 'flag_ac', 'flag_site']
        flag_lists = get_flag_lists(input_path, ats_in, flag_bands)
        idfile = 1
        satellite_id_ref = 0
        list_files = []
        for name in os.listdir(input_path):
            if not name.endswith('.nc'):
                continue
            file_in = os.path.join(input_path, name)
            reader = MDB_READER(file_in, True)
            file_out = os.path.join(input_path, f'Temp_{idfile}.nc')
            creating_copy_with_flag_bands(reader, file_out, flag_lists, satellite_id_ref)
            satellite_id_ref = satellite_id_ref + reader.mfile.n_mu_total
            list_files.append(file_out)
            idfile = idfile + 1
        concatenate_nc_impl(list_files, input_path, ncout_file)

    ##PLOTTING
    if args.input_path and args.config_file and mode == 'PLOT':
        input_path = args.input_path
        if not os.path.exists(input_path):
            print(f'[ERROR] {input_path} does not exist')
            return
        config_file = args.config_file
        if not os.path.exists(config_file):
            print(f'[ERROR] {config_file} does not exist')
            return
        if args.output:
            output_path = args.output
        else:
            path_mdb = os.path.dirname(input_path)
            output_path = os.path.join(path_mdb, 'PLOTS')
        if not os.path.exists(output_path):
            try:
                os.mkdir(output_path)
            except:
                pass
        if not os.path.isdir(output_path):
            print(f'[ERROR] Ouput path: {output_path} does not exist or is not a directory')

        from MDBPlotV2 import MDBPlot
        mplot = MDBPlot(input_path)
        mplot.output_path = output_path

        import configparser
        options = configparser.ConfigParser()
        options.read(config_file)
        print(mplot.VALID)
        mplot.plot_from_options(options)


if __name__ == '__main__':
    main()
