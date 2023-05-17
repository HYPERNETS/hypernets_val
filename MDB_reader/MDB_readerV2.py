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
                'type': 'i2'
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
                'fillvalue': -1,
                'type': 'i1'
            },
            'mu_insitu_id': {
                'namedf': 'mu_insitu_id',
                'fillvalue': -1,
                'type': 'i1'
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

        ##validitiy of spectrums
        if args.verbose:
            print('[INFO] Check validity spectra...')
        new_var = new_MDB.createVariable('insitu_valid', 'i1', ('satellite_id', 'insitu_id'), zlib=True, complevel=6,
                                         fill_value=-1)
        for index_mu in range(self.mfile.n_mu_total):
            if (index_mu%100)==0 and args.verbose:
                print(f'[INFO] MU: {index_mu} of {self.mfile.n_mu_total}')
            validity_spectra = self.mfile.qc_insitu.check_validity_spectra_mu(index_mu)
            new_var[index_mu] = validity_spectra[:]

        new_MDB.close()
        if args.verbose:
            print(f'[INFO] Completed')

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
        # print(name)
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


def creating_copy_correcting_band(reader, file_out):
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
        # print(name)
        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(reader.mfile.nc[name].__dict__)

        # copy variable data
        ncout[name][:] = reader.mfile.nc[name][:]
        if name == 'satellite_Rrs':
            ncout[name][:, 1, :, :] = reader.mfile.nc[name][:, 1, :, :] * np.pi

    ncout.close()
    reader.mfile.close()
    return True


def creating_copy_with_eight_bands(reader, file_out):
    # reader = MDB_READER('', True)
    from netCDF4 import Dataset
    import numpy as np
    ncout = Dataset(file_out, 'w', format='NETCDF4')

    # copy global attributes all at once via dictionary
    ncout.setncatts(reader.mfile.nc.__dict__)

    # copy dimensions
    for name, dimension in reader.mfile.nc.dimensions.items():
        if name == 'satellite_bands':
            ncout.createDimension(name, 8)
        else:
            ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

    indices_good = [0, 1, 2, 3, 4, 5, 6, 8]

    # copy variables
    for name, variable in reader.mfile.nc.variables.items():
        fill_value = None
        if '_FillValue' in list(reader.mfile.nc.ncattrs()):
            fill_value = variable._FillValue
        # print(name)
        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(reader.mfile.nc[name].__dict__)

        # copy variable data
        if name == 'satellite_bands':
            ncout[name][:] = reader.mfile.nc[name][indices_good]
        elif name == 'satellite_Rrs':
            ncout[name][:, :, :, :] = reader.mfile.nc[name][:, indices_good, :, :]
        else:
            ncout[name][:] = reader.mfile.nc[name][:]

    ncout.close()
    reader.mfile.close()
    return True

def creating_copy_with_valid_indices(reader,file_out,indices_sat_valid):
    from netCDF4 import Dataset
    import numpy as np
    ncout = Dataset(file_out, 'w', format='NETCDF4')

    # copy global attributes all at once via dictionary
    ncout.setncatts(reader.mfile.nc.__dict__)

    # copy dimensions
    for name, dimension in reader.mfile.nc.dimensions.items():
        ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

        # copy variables
    for name, variable in reader.mfile.nc.variables.items():
        fill_value = None
        if '_FillValue' in list(reader.mfile.nc.ncattrs()):
            fill_value = variable._FillValue
        # print(name)
        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(reader.mfile.nc[name].__dict__)

        if name=='mu_valid':
            mu_valid_array = np.array(reader.mfile.nc[name][:])
            mu_valid_array[indices_sat_valid==False] = 0
            ncout[name][:] = mu_valid_array[:]
        else:
            ncout[name][:] = reader.mfile.nc[name][:]

    ncout.close()
    reader.mfile.close()
    return True

def getting_common_matchups():
    # getting common match-ups
    dir_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/CONTATENATED_WITHOUT_COMMON_MATCHUPS'
    file_in = os.path.join(dir_base, 'MDBrc_ALL.nc')
    mu = {}
    from netCDF4 import Dataset
    import numpy as np
    from datetime import datetime as dt
    from datetime import timedelta
    dataset = Dataset(file_in)
    sat_time = np.array(dataset.variables['satellite_time'][:])
    satellite = np.array(dataset.variables['flag_satellite'][:])
    processor = np.array(dataset.variables['flag_ac'][:])
    valid = np.array(dataset.variables['mu_valid'])
    nsat = sat_time.shape[0]
    for idx in range(nsat):
        if not valid[idx]:
            continue
        time = sat_time[idx]
        time_here = dt(1970, 1, 1) + timedelta(seconds=int(time))
        time_here_str = time_here.strftime('%Y%m%d')
        sat = 'A'
        if satellite[idx] == 2:
            sat = 'B'
        ac = 'ACOLITE'
        if processor[idx] == 2:
            ac = 'C2RCC'
        if time_here_str not in mu:
            mu[time_here_str] = {
                'A': {
                    'ACOLITE': False,
                    'C2RCC': False
                },
                'B': {
                    'ACOLITE': False,
                    'C2RCC': False
                }

            }

        if ac == 'ACOLITE':
            mu[time_here_str][sat]['ACOLITE'] = True
        elif ac == 'C2RCC':
            mu[time_here_str][sat]['C2RCC'] = True
    list_dates = []
    list_sat = []
    cm = {}
    for datestr in mu:
        add_a = False
        add_b = False
        if mu[datestr]['A']['ACOLITE'] and mu[datestr]['A']['C2RCC']:
            list_dates.append(datestr)
            list_sat.append('A')
            add_a = True
        if mu[datestr]['B']['ACOLITE'] and mu[datestr]['B']['C2RCC']:
            list_dates.append(datestr)
            list_sat.append('B')
            add_b = True
        if sat is not None:
            if datestr not in cm.keys():
                # print('llega aqui')
                cm[datestr] = {
                    'A': False,
                    'B': False
                }
            if add_a:
                cm[datestr]['A'] = True
            if add_b:
                cm[datestr]['B'] = True
            # print(cm[datestr]['A'])
            # print('-------------------')

    return list_dates, list_sat, cm


def getting_common_satellite_id(file_in, cm):
    from netCDF4 import Dataset
    import numpy as np
    from datetime import datetime as dt
    from datetime import timedelta
    dataset = Dataset(file_in)
    sat_time = np.array(dataset.variables['satellite_time'][:])
    satellite = np.array(dataset.variables['flag_satellite'][:])
    valid = np.array(dataset.variables['mu_valid'][:])
    nsat = sat_time.shape[0]
    indices = []
    indices_valid = [False] * nsat
    used = {}


    for idx in range(nsat):
        if not valid[idx]:
            continue
        time = sat_time[idx]
        time_here = dt(1970, 1, 1) + timedelta(seconds=int(time))
        time_here_str = time_here.strftime('%Y%m%d')
        if time_here_str not in cm.keys():
            continue
        if time_here_str not in used.keys():
            used[time_here_str] = {
                'A': 0,
                'B': 0
            }

        sat = 'A'
        if satellite[idx] == 2:
            sat = 'B'

        if cm[time_here_str][sat] and used[time_here_str][sat]<2:
            indices.append(idx)
            indices_valid[idx] =  True
            used[time_here_str][sat] = used[time_here_str][sat]+1

    print(len(indices))

    indices_mu = np.array(dataset.variables['mu_satellite_id'][:])
    nmu = indices_mu.shape[0]
    indices_mu_common = []
    indices_mu_common_valid = [False] * nmu
    for idx in range(nmu):
        index_here = indices_mu[idx]
        if indices_valid[index_here]:
            indices_mu_common.append(idx)
            indices_mu_common_valid[idx] = True


    return np.array(indices_valid,dtype=bool),np.array(indices_mu_common_valid,dtype=bool)


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
        # fconfig = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/config_qc.ini'
        # import configparser
        # from QC_OPTIONS import QC_OPTIONS
        # # from QC_SAT import  QC_SAT
        # options = configparser.ConfigParser()
        # options.read(fconfig)
        # qco = QC_OPTIONS(options)
        # qco.getting_qcsat(None)

        # dir_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/LPAR'
        # dir_deprecated = os.path.join(dir_base,'DEPRECATED')
        # for name in os.listdir(dir_deprecated):
        #     if name.startswith('MDB_'):
        #         file_in = os.path.join(dir_deprecated,name)
        #         file_out = os.path.join(dir_base,name)
        #         print(file_in,'-->',file_out)
        #         reader = MDB_READER(file_in, True)
        #         creating_copy_correcting_band(reader, file_out)

        # dir_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/ACOLITE_11BANDS'
        # dir_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/ALLr'
        # for name in os.listdir(dir_base):
        #     file_in = os.path.join(dir_base, name)
        #     file_out = os.path.join(dir_out,name)
        #     print(file_in,'-->',file_out)
        #     reader = MDB_READER(file_in, True)
        #     creating_copy_with_eight_bands(reader,file_out)

        # name_a = 'MDB_S2A_MSI_20M_ACOLITE_20220101T000000_20230313T235959_HYPSTAR_VEIT.nc'
        # file_a = os.path.join(dir_base, 'DEPRECATED', name_a)
        # file_ao = os.path.join(dir_base,name_a)
        # reader = MDB_READER(file_a, True)
        # creating_copy_correcting_band(reader,file_ao)
        #
        # name_b = 'MDB_S2B_MSI_20M_ACOLITE_20220101T000000_20230313T235959_HYPSTAR_VEIT.nc'
        # file_b = os.path.join(dir_base, 'DEPRECATED',name_b)
        # file_bo = os.path.join(dir_base, name_b)
        # reader = MDB_READER(file_b, True)
        # creating_copy_correcting_band(reader, file_bo)

        list_dates, list_sat, cm = getting_common_matchups()
        na = 0
        nb = 0
        for datestr in cm:
            if cm[datestr]['A']:
                na = na + 1
            if cm[datestr]['B']:
                nb = nb + 1
        ncommon = na + nb
        print(na, nb, ncommon)

        dir_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI'
        file_in = os.path.join(dir_base, 'CONTATENATED_WITHOUT_COMMON_MATCHUPS','MDBrc_ALL.nc')
        file_out = os.path.join(dir_base,'MRDrc_ALL.nc')
        indices, indices_mu = getting_common_satellite_id(file_in, cm)


        import numpy as np
        print(np.sum(indices))
        # print(np.sum(indices_mu))
        reader = MDB_READER(file_in, True)
        creating_copy_with_valid_indices(reader,file_out,indices)



        return

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
            reader.mfile.qc_sat = qco.get_qcsat(reader.mfile.qc_sat, reader.mfile.nc)
            reader.mfile.qc_insitu.ncdataset = reader.mfile.nc
            reader.mfile.qc_insitu = qco.get_qc_insitu(reader.mfile.qc_insitu)
            if reader.mfile.qc_insitu is None:
                return
            if not reader.mfile.qc_insitu.apply_nir_correction:
                reader.mfile.qc_insitu.insitu_rrs = reader.mfile.variables['insitu_Rrs_nosc']
                if args.verbose:
                    print('[INFO] NIR Correction is not applied. Using insitu_Rrs_nosc as input variable')

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
        # print(mplot.VALID)
        #mplot.plot_from_options(options)

        # path_img = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/PLOTS'
        # from PlotMultiple import PlotMultiple
        # file_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/PLOTS/'

        # path_img = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/PLOTS'
        # from PlotMultiple import PlotMultiple
        # file_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/PLOTS/'

        path_img = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/CMEMSMULTI/VEIT/PLOTS'
        from PlotMultiple import PlotMultiple
        file_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/CMEMSMULTI/VEIT/PLOTS'

        # mplot.mrfile.save_wl_images(path_img)
        # mplot.mrfile.save_rgb_images(path_img)

        ##FLAG ANALYSIS
        # mplot.mrfile.window_size = 3
        # flag_l = 'LAND,COASTLINE,CLOUD,CLOUD_AMBIGUOUS,CLOUD_MARGIN,INVALID,COSMETIC,SATURATED,SUSPECT,HISOLZEN,HIGHGLINT,SNOW_ICE,AC_FAIL,WHITECAPS,RWNEG_O2,RWNEG_O3,RWNEG_O4,RWNEG_O5,RWNEG_O6,RWNEG_O7,RWNEG_O8'
        # #flag_l = 'IDEPIX_INVALID, IDEPIX_CLOUD, IDEPIX_CLOUD_AMBIGUOUS, IDEPIX_CLOUD_SURE, IDEPIX_CLOUD_BUFFER, IDEPIX_CLOUD_SHADOW, IDEPIX_SNOW_ICE, IDEPIX_BRIGHT, IDEPIX_WHITE, IDEPIX_COASTLINE, IDEPIX_LAND, IDEPIX_CIRRUS_SURE, IDEPIX_CIRRUS_AMBIGUOUS, IDEPIX_CLEAR_LAND, IDEPIX_BRIGHTWHITE, IDEPIX_VEG_RISK, IDEPIX_MOUNTAIN_SHADOW, IDEPIX_POTENTIAL_SHADOW, IDEPIX_CLUSTERED_CLOUD_SHADOW'
        # flag_list = [x.strip() for x in flag_l.split(',')]
        # #flag_list = ['CLOUD']
        # flag_info = mplot.mrfile.analyze_sat_flags('satellite_WQSF',flag_list)
        # for flag in flag_info:
        #    print(flag,flag_info[flag]['nmacrow'],flag_info[flag]['pmacrow'])
        # file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/PLOTS/FlagTotal.jpg'
        # #file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/LOTS/FlagTotal.jpg'
        # mplot.plot_flag(flag_info,file_out)
        # # file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/PLOTS/Water.jpg'
        # # mplot.plot_flag_array(flag_info,'CLOUD',file_out)

        # MU ANALYSIS
        # file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/PLOTS/MUALL_SITE.jpg'
        # mplot.mrfile.analyse_mu_temporal_flag(False,'flag_site',file_out)
        # file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/PLOTS/MUVALID_SITE.jpg'
        # mplot.mrfile.analyse_mu_temporal_flag(True, 'flag_site',file_out)
        # file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/PLOTS/MUVALID_SITE_PORC.jpg'
        # mplot.mrfile.analyse_mu_flag('flag_site',file_out)

        # file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/PLOTS/MUALL_SITE.jpg'
        # mplot.mrfile.analyse_mu_temporal_flag(False,'flag_site',file_out)
        # file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/PLOTS/MUVALID_SITE.jpg'
        # mplot.mrfile.analyse_mu_temporal_flag(True, 'flag_site',file_out)
        # file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/PLOTS/MUVALID_SITE_PORC.jpg'
        # mplot.mrfile.analyse_mu_flag('flag_site',file_out)

        # Coverage
        # file00 = os.path.join(file_base,'MUALL_SITE.jpg')
        # file10 = os.path.join(file_base,'MUVALID_SITE.jpg')
        # file_out = os.path.join(file_base,'MUDistribution.jpg')
        # pm = PlotMultiple()
        # pm.start_multiple_plot(2,1)
        # pm.plot_image(file00,0,0)
        # pm.plot_image(file10, 1, 0)
        # pm.set_text('a',1850,-50)
        # pm.set_text('b', 1850, 1400)
        # pm.save_fig(file_out)
        # pm.close_plot()

        # PorcValidSiteFlag.jpg
        # file00 = os.path.join(file_base,'MUVALID_SITE_PORC.jpg')
        # file01 = os.path.join(file_base,'FlagTotal.jpg')
        # file_out = os.path.join(file_base,'PorcValidSiteFlag.jpg')
        # pm = PlotMultiple()
        # pm.start_multiple_plot(1,2)
        # pm.plot_image(file00,0,0)
        # pm.plot_image(file01, 0, 1)
        # pm.set_text('a',-50,1400)
        # pm.set_text('b', 1850, 1400)
        # pm.save_fig(file_out)
        # pm.close_plot()

        # scatterplots
        # file00 = os.path.join(file_base,'OVERAL_SCATTERPLOT_ACOLITE.jpg')
        # file10 = os.path.join(file_base,'OVERAL_SCATTERPLOT_C2RCC.jpg')
        # file_out = os.path.join(file_base,'OverallScatterplot.jpg')
        # pm = PlotMultiple()
        # pm.start_multiple_plot(2,1)
        # pm.plot_image(file00,0,0)
        # pm.plot_image(file10, 1, 0)
        # pm.set_text('ACOLITE', 800, -1350)
        # pm.set_text('C2RCC',800,100)
        # pm.save_fig(file_out)
        # pm.close_plot()

        # stats site
        # file00 = os.path.join(file_base, 'PLOT_STAT_flag_site_RMSD.jpg')
        # file10 = os.path.join(file_base, 'PLOT_STAT_flag_site_DETER(r2).jpg')
        # file20 = os.path.join(file_base, 'PLOT_STAT_flag_site_BIAS.jpg')
        # file_out = os.path.join(file_base, 'StatsByWL.jpg')
        # pm = PlotMultiple()
        # pm.start_multiple_plot(3, 1)
        # pm.plot_image(file00, 0, 0)
        # pm.plot_image(file10, 1, 0)
        # pm.plot_image(file20, 2, 0)
        # pm.set_text('a',1800,-1475)
        # pm.set_text('b', 1800, -50)
        # pm.set_text('c', 1800, 1400)
        # pm.save_fig(file_out)
        # pm.close_plot()

        # stats ac
        # file00 = os.path.join(file_base,'PLOT_STAT_flag_ac_RMSD.jpg')
        # file10 = os.path.join(file_base,'PLOT_STAT_flag_ac_DETER(r2).jpg')
        # file20 = os.path.join(file_base, 'PLOT_STAT_flag_ac_BIAS.jpg')
        # file_out = os.path.join(file_base,'StatsByWL.jpg')
        # pm = PlotMultiple()
        # pm.start_multiple_plot(3,1)
        # pm.plot_image(file00,0,0)
        # pm.plot_image(file10, 1, 0)
        # pm.plot_image(file20, 2, 0)
        # pm.set_text('a',1800,-1475)
        # pm.set_text('b', 1800, -50)
        # pm.set_text('c', 1800, 1400)
        # pm.save_fig(file_out)
        # pm.close_plot()

        # stats gloal
        file00 = os.path.join(file_base,'PLOT_STAT_GLOBAL_RMSD.jpg')
        file10 = os.path.join(file_base,'PLOT_STAT_GLOBAL_DETER(r2).jpg')
        file20 = os.path.join(file_base, 'PLOT_STAT_GLOBAL_BIAS.jpg')
        file_out = os.path.join(file_base,'StatsByWL.jpg')
        pm = PlotMultiple()
        pm.start_multiple_plot(3,1)
        pm.plot_image(file00,0,0)
        pm.plot_image(file10, 1, 0)
        pm.plot_image(file20, 2, 0)
        pm.set_text('a',1800,-1475)
        pm.set_text('b', 1800, -50)
        pm.set_text('c', 1800, 1400)
        pm.save_fig(file_out)
        pm.close_plot()


        # sites
        # file00 = os.path.join(file_base,'PLOT_SAT_INS_COMPARISON_VEIT.jpg')
        # file01 = os.path.join(file_base, 'PLOT_SAT_INS_COMPARISON_BEFR.jpg')
        # file10 = os.path.join(file_base, 'PLOT_SAT_INS_COMPARISON_LPAR.jpg')
        # file11 = os.path.join(file_base,'PLOT_SAT_INS_COMPARISON_MAFR.jpg')
        # file20 = os.path.join(file_base, 'PLOT_SAT_INS_COMPARISON_GAIT.jpg')
        # file21 = os.path.join(file_base, 'PLOT_SAT_INS_COMPARISON_M1BE.jpg')
        # file_out = os.path.join(file_base,'PLOT_SAT_INS_COMPARISON.jpg')
        # pm = PlotMultiple()
        # pm.start_multiple_plot(3,2)
        # pm.plot_image(file00, 0, 0)
        # pm.plot_image(file01, 0, 1)
        # pm.plot_image(file10, 1, 0)
        # pm.plot_image(file11, 1, 1)
        # pm.plot_image(file20, 2, 0)
        # pm.plot_image(file21, 2, 1)
        # pm.set_text('VEIT',-200,-1475)
        # pm.set_text('BEFR', 1700, -1475)
        # pm.set_text('LPAR', -200,-45)
        # pm.set_text('MAFR', 1700, -45)
        # pm.set_text('GAIT', -200, 1380)
        # pm.set_text('M1BE', 1700, 1380)
        # pm.save_fig(file_out)
        # pm.close_plot()

        # comparison ab
        # file00 = os.path.join(file_base, 'OVERAL_SCATTERPLOT_SAT.jpg')
        # file01 = os.path.join(file_base, 'PLOT_STAT_SAT_flag_satellite_BIAS.jpg')
        # file10 = os.path.join(file_base, 'PLOT_STAT_SAT_flag_satellite_RMSD.jpg')
        # file11 = os.path.join(file_base, 'PLOT_STAT_SAT_flag_satellite_DETER(r2).jpg')
        # file_out = os.path.join(file_base, 'ComparisonAB.jpg')
        # pm = PlotMultiple()
        # pm.start_multiple_plot(2, 2)
        # pm.plot_image(file00, 0, 0)
        # pm.plot_image(file01, 0, 1)
        # pm.plot_image(file10, 1, 0)
        # pm.plot_image(file11, 1, 1)
        # pm.set_text('a',-100,-50)
        # pm.set_text('b', 1800, -50)
        # pm.set_text('c', -100, 1375)
        # pm.set_text('d', 1800, 1375)
        # pm.save_fig(file_out)
        # pm.close_plot()

        # comparison site
        # file00 = os.path.join(file_base, 'PLOT_SAT_INS_COMPARISON_VEIT_ACOLITE.jpg')
        # file01 = os.path.join(file_base, 'PLOT_SAT_INS_COMPARISON_VEIT_C2RCC.jpg')
        # file10 = os.path.join(file_base, 'PLOT_SAT_INS_COMPARISON_BEFR_ACOLITE.jpg')
        # file11 = os.path.join(file_base, 'PLOT_SAT_INS_COMPARISON_BEFR_C2RCC.jpg')
        # file_out = os.path.join(file_base, 'Sites_S2.jpg')
        # pm = PlotMultiple()
        # pm.start_multiple_plot(2, 2)
        # pm.plot_image(file00, 0, 0)
        # pm.plot_image(file01, 0, 1)
        # pm.plot_image(file10, 1, 0)
        # pm.plot_image(file11, 1, 1)
        # pm.set_text('a. VEIT ACOLITE',-400,-40)
        # pm.set_text('b. VEIT C2RCC', 1500, -40)
        # pm.set_text('c. BEFR ACOLITE', -400, 1385)
        # pm.set_text('d. BEFR C2RCC', 1500, 1385)
        # pm.save_fig(file_out)
        # pm.close_plot()


        # spectra, stats = mplot.mrfile.get_all_insitu_valid_spectra()
        # from PlotSpectra import PlotSpectra
        # pspectra = PlotSpectra()
        # wavelength = mplot.mrfile.get_insitu_wl()
        # #spectra = spectra[:, 0:990]
        # pspectra.plot_multiple_spectra(wavelength,spectra,stats,None,800)
        # # pspectra.xdata = mplot.mrfile.get_insitu_wl()[0:990]
        # # spectra = spectra[:, 0:990]
        # # pspectra.plot_data(spectra.transpose(), pspectra.line_style_default)
        # file_out_here = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/tal.jpg'
        # pspectra.save_plot(file_out_here)
        # pspectra.close_plot()


if __name__ == '__main__':
    main()
