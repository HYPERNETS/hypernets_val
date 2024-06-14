import datetime

import matplotlib.pyplot as plt
import numpy as np
import pytz
import os.path
import sys
import argparse
import warnings

import pandas as pd

warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', RuntimeWarning)

from MDBFile import MDBFile
from MDB_builder.INSITU_base import INSITUBASE

parser = argparse.ArgumentParser(
    description="Match-ups extraction from MDB files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument("-m", "--mode", help="Mode",
                    choices=["GENERATEMU", "GENERATEMU_S", "CONCATENATE", "REMOVEREP", "PLOT", "PLOT_CSV", "COMMONMU",
                             "COMMONMU_NOSAT",
                             "COMMONMU_INS", "CHECK_WL", "UPDATE_SAT_WL", "UPDATE_INSITU_WL", "CHECK_SAT_TIME",
                             "CHECK_PROTOCOLS", "TEST", "ADDFLAGBAND"],
                    required=True)
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-i', "--input_path", help="Input MDB path")
parser.add_argument('-o', "--output", help="Path to output")
parser.add_argument('-arep', "--allow_repeated", help="Set to allow match-ups on the same date. Only for transects.",
                    action="store_true")
parser.add_argument('-rmdb', "--reduce_mdbr", help="MDBr should be reduced to only one insitu_id",
                    action="store_true")

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

    def create_mdb_results_file(self, fout, reduce_mdbr):
        if not self.mfile.VALID:
            return

        if self.mfile.df_validation is None:
            nmu_valid, df_valid = self.mfile.prepare_df_validation()

        foutcsv = fout.replace('.nc', '_summary.csv')
        foutcsv = foutcsv.replace('MDBr', 'CSVr')
        df_valid.to_csv(foutcsv, sep=';')

        ibase = INSITUBASE(None)

        if reduce_mdbr:
            new_MDB = ibase.copy_nc_reduced(self.mfile.file_path, fout)
        else:
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
                'type': 'f8'
            },
            'mu_ins_time': {
                'namedf': 'Ins_Time',
                'fillvalue': -999.0,
                'type': 'f8'
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
            'mu_valid_complete': {
                'namedf': 'spectrum_complete',
                'fillvalue': -1,
                'type': 'i1'
            },
            'mu_insitu_id': {
                'namedf': 'mu_insitu_id',
                'fillvalue': -1,
                'type': 'i2'
            }
        }

        if 'insitu_latitude' in self.mfile.variables and 'insitu_longitude' in self.mfile.variables:
            new_variables_sat_mu['mu_insitu_latitude'] = {
                'namedf': 'Ins_Lat',
                'fillvalue': -999.0,
                'type': 'f4'
            }
            new_variables_sat_mu['mu_insitu_longitude'] = {
                'namedf': 'Ins_Lon',
                'fillvalue': -999.0,
                'type': 'f4'
            }

        if self.mfile.df_mu is None:
            self.mfile.prepare_df_mu()

        foutcsv = fout.replace('.nc', '.csv')
        foutcsv = foutcsv.replace('MDBr', 'CSVr')
        status = self.mfile.df_validation[self.mfile.df_validation['Index_Band'] == 0]['status']
        status_array = np.array(status)
        self.mfile.df_mu['status'][:] = status_array

        self.mfile.df_mu.to_csv(foutcsv, sep=';')

        for new_var_name in new_variables_sat_mu:
            new_var = new_MDB.createVariable(new_var_name, new_variables_sat_mu[new_var_name]['type'],
                                             ('satellite_id',),
                                             zlib=True,
                                             complevel=6)
            array = np.array(self.mfile.df_mu.loc[:, new_variables_sat_mu[new_var_name]['namedf']])

            if new_var_name == 'mu_valid':
                array_status = np.array(self.mfile.df_mu.loc[:, 'status'])

            if new_var_name == 'mu_sat_time' or new_var_name == 'mu_ins_time':
                array_t = []
                for idx in range(len(array)):
                    val = array[idx]
                    try:
                        valdt = datetime.datetime.strptime(val, '%Y-%m-%d %H:%M').replace(tzinfo=pytz.utc)
                        array_t.append(float(valdt.timestamp()))
                    except:
                        array_t.append(-999.0)
                array = np.array(array_t, dtype=np.float64)

            fillValue = new_variables_sat_mu[new_var_name]['fillvalue']
            if fillValue is not None:
                array = ma.masked_array(array, mask=array == fillValue)

            for idx in range(self.mfile.n_mu_total):

                if new_var_name == 'mu_valid':
                    status_here = array_status[idx] * (-1)
                    val = array[idx]
                    if 3 >= status_here >= 1:
                        val = -1
                    new_var[idx] = [val]
                else:
                    # print(idx, '-->', array[idx])
                    new_var[idx] = [array[idx]]

        ##validitiy of spectrums
        if args.verbose:
            print('[INFO] Check validity spectra...')
        new_var = new_MDB.createVariable('insitu_valid', 'i1', ('satellite_id', 'insitu_id'), zlib=True, complevel=6,
                                         fill_value=-1)
        if reduce_mdbr:
            array_mu_insitu_id = np.array(self.mfile.df_mu.loc[:, new_variables_sat_mu['mu_insitu_id']['namedf']])

            for var_name in new_MDB.variables:
                if var_name.startswith(
                        'insitu') and var_name != 'insitu_valid' and var_name != 'insitu_original_bands' and var_name != 'insitu_Rrs':
                    for index_mu in range(self.mfile.n_mu_total):
                        insitu_id_here = array_mu_insitu_id[index_mu]
                        new_MDB.variables[var_name][index_mu, 0] = self.mfile.nc.variables[var_name][
                            index_mu, insitu_id_here]
            for index_mu in range(self.mfile.n_mu_total):
                insitu_id_here = array_mu_insitu_id[index_mu]
                new_MDB.variables['insitu_Rrs'][index_mu, :, 0] = self.mfile.nc.variables['insitu_Rrs'][index_mu, :,
                                                                  insitu_id_here]

            # for index_mu in range(self.mfile.n_mu_total):
            #     insitu_id_here = array_mu_insitu_id[index_mu]
            #     if (index_mu % 100) == 0 and args.verbose:
            #         print(f'[INFO] Checking spectra validity for mu {index_mu} of {self.mfile.n_mu_total} -> {insitu_id_here}')
            #     validity_spectra = self.mfile.qc_insitu.check_validity_spectra_mu(index_mu)
            #     new_var[index_mu] = validity_spectra[insitu_id_here]
            #     new_MDB.variables['mu_insitu_id'][index_mu] = [0]

        else:
            for index_mu in range(self.mfile.n_mu_total):
                if (index_mu % 100) == 0 and args.verbose:
                    print(f'[INFO] Checking spectra validity for mu : {index_mu} of {self.mfile.n_mu_total}')
                validity_spectra = self.mfile.qc_insitu.check_validity_spectra_mu(index_mu)
                new_var[index_mu] = validity_spectra[:]

        new_MDB.close()
        if args.verbose:
            print(f'[INFO] Completed')

        return nmu_valid

    def create_mdb_results_single(self, fout):
        if not os.path.exists(fout):
            print(f'[INFO] Creating copy of {self.mfile.file_path}...')
            ibase = INSITUBASE(None)
            new_MDB = ibase.copy_nc(self.mfile.file_path, fout)
            new_MDB.createDimension('mu_id', None)
            print(f'[INFO] Copy completed')
        else:
            from netCDF4 import Dataset
            new_MDB = Dataset(fout, 'a')

        create = self.mfile.create_mu_single_variables(new_MDB,False,False,None)

        new_MDB.close()

        if not create:
            print(f'[ERROR] Error creating file: {fout}')
            #os.remove(fout)
        else:
            print(f'[INFO] File created: {fout}')


        # if self.mfile.df_validation is None:
        #     nmu_valid, df_valid = self.mfile.prepare_df_validation_single()

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


def get_band_list(input_path):
    band_list = {}
    for name in os.listdir(input_path):
        if not name.endswith('.nc'):
            continue
        file_in = os.path.join(input_path, name)
        reader = MDB_READER(file_in, True)
        for namevar in reader.mfile.variables:
            var = reader.mfile.variables[namevar]
            if namevar not in band_list.keys():
                # band_list.append(var)
                band_list[namevar] = {}
                atts = var.ncattrs()
                for at in atts:
                    band_list[namevar][at] = var.getncattr(at)
            else:
                atts = var.ncattrs()
                for at in atts:
                    if at not in band_list[namevar].keys():
                        band_list[namevar][at] = var.getncattr(at)

        reader.mfile.close()
    return band_list


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


def get_wl_list_from_file(fwl):
    wl_list = []
    f1 = open(fwl, 'r')
    for line in f1:
        try:
            wl = float(line.strip())
            wl_list.append(wl)
        except:
            print(f'[ERROR] {line.strip()} is not a valid wavelength')
            return None
    if len(wl_list) == 0:
        return None
    return wl_list


def update_sat_extract_file_version(input_file):
    from netCDF4 import Dataset
    input_dataset = Dataset(input_file)
    output_file = os.path.join(os.path.dirname(input_file), 'Temp.nc')
    ncout = Dataset(output_file, 'w', format='NETCDF4')

    # copy attributs
    for at in input_dataset.ncattrs():
        val = input_dataset.getncattr(at)
        if at == 'insitu_site_name':
            at = 'site'
        if at == 'insitu_lat':
            at = 'in_situ_lat'
        if at == 'insitu_lon':
            at = 'in_situ_lon'
        ncout.setncattr(at, val)

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    for name, variable in input_dataset.variables.items():
        if name == 'satellite_PDU':
            continue
        newname = name
        if name == 'OZA' or name == 'OAA' or name == 'SAA' or name == 'SZA':
            newname = f'satellite_{name}'

        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(newname, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)

        # copy variable attributes all at once via dictionary
        ncout[newname].setncatts(input_dataset[name].__dict__)

        ncout[newname][:] = input_dataset[name][:]

    ncout.close()
    input_dataset.close()

    os.rename(output_file, input_file)


def update_mdb_file_version(input_file):
    from netCDF4 import Dataset
    input_dataset = Dataset(input_file)
    output_file = os.path.join(os.path.dirname(input_file), 'Temp.nc')
    ncout = Dataset(output_file, 'w', format='NETCDF4')

    # copy attributs
    for at in input_dataset.ncattrs():
        val = input_dataset.getncattr(at)
        if at == 'insitu_site_name':
            at = 'site'
        if at == 'insitu_lat':
            at = 'in_situ_lat'
        if at == 'insitu_lon':
            at = 'in_situ_lon'
        ncout.setncattr(at, val)

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    add_mu_valid_complete = True
    for name, variable in input_dataset.variables.items():
        if name == 'mu_valid_complete':
            add_mu_valid_complete = False
        if name == 'satellite_PDU':
            continue
        if name == 'insitu_filename':
            continue
        newname = name
        if name == 'OZA' or name == 'OAA' or name == 'SAA' or name == 'SZA':
            newname = f'satellite_{name}'

        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(newname, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)

        # copy variable attributes all at once via dictionary
        ncout[newname].setncatts(input_dataset[name].__dict__)

        ncout[newname][:] = input_dataset[name][:]

    if add_mu_valid_complete:
        var = ncout.createVariable('mu_valid_complete', 'i1', ('satellite_id',), zlib=True, shuffle=True, complevel=6)
        var[:] = 1

    os.rename(output_file, input_file)


def creating_copy_with_flag_band(input_file, name_var, flag_values, flag_meanings, value_array):
    from netCDF4 import Dataset
    input_dataset = Dataset(input_file)
    output_file = os.path.join(os.path.dirname(input_file), 'Temp.nc')
    ncout = Dataset(output_file, 'w', format='NETCDF4')

    # copy attributs
    for at in input_dataset.ncattrs():
        val = input_dataset.getncattr(at)
        if at == 'insitu_site_name':
            at = 'site'
        if at == 'insitu_lat':
            at = 'in_situ_lat'
        if at == 'insitu_lon':
            at = 'in_situ_lon'
        ncout.setncattr(at, val)

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    for name, variable in input_dataset.variables.items():
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(input_dataset[name].__dict__)

        ncout[name][:] = input_dataset[name][:]

    var_flag = ncout.createVariable(name_var, 'i4', ('satellite_id',), zlib=True, shuffle=True, complevel=6)
    var_flag.flag_meanings = ' '.join(flag_meanings)
    var_flag.flag_values = flag_values
    var_flag[:] = value_array
    ncout.close()
    input_dataset.close()

    os.rename(output_file, input_file)


def creating_copy_correcting_attributes(input_file):
    from netCDF4 import Dataset
    input_dataset = Dataset(input_file)
    output_file = os.path.join(os.path.dirname(input_file), 'Temp.nc')
    ncout = Dataset(output_file, 'w', format='NETCDF4')

    # copy attributs
    for at in input_dataset.ncattrs():
        val = input_dataset.getncattr(at)
        if at == 'insitu_site_name':
            at = 'site'
        if at == 'insitu_lat':
            at = 'in_situ_lat'
        if at == 'insitu_lon':
            at = 'in_situ_lon'
        ncout.setncattr(at, val)

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    for name, variable in input_dataset.variables.items():
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(input_dataset[name].__dict__)

        ncout[name][:] = input_dataset[name][:]

    ncout.close()
    input_dataset.close()

    os.rename(output_file, input_file)


def creating_copy_changing_valid_mu(input_file, output_file, valid_mu, valid_mu_common):
    from netCDF4 import Dataset
    input_dataset = Dataset(input_file)
    ncout = Dataset(output_file, 'w', format='NETCDF4')

    # copy global attributes all at once via dictionary
    ncout.setncatts(input_dataset.__dict__)

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    for name, variable in input_dataset.variables.items():
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(input_dataset[name].__dict__)

        # copy data
        if name == 'mu_valid':
            ncout[name][:] = valid_mu[:]
        elif name == 'mu_valid_common':
            ncout[name][:] = valid_mu_common[:]
        else:
            ncout[name][:] = input_dataset[name][:]

    ncout.close()
    input_dataset.close()


def creating_copy_adding_new_flag(input_file, output_file, new_flag):
    from netCDF4 import Dataset
    input_dataset = Dataset(input_file)
    ncout = Dataset(output_file, 'w', format='NETCDF4')

    # copy global attributes all at once via dictionary
    ncout.setncatts(input_dataset.__dict__)

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    for name, variable in input_dataset.variables.items():
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(input_dataset[name].__dict__)

        # copy data
        if name == 'satellite_WQSF':
            ncout[name][:] = new_flag[:]
        else:
            ncout[name][:] = input_dataset[name][:]

    ncout.close()
    input_dataset.close()


def creating_copy_correcting_sat_time(input_file, output_file):
    from netCDF4 import Dataset
    import numpy as np
    from datetime import datetime as dt
    # from datetime import timedelta
    from datetime import timezone
    input_dataset = Dataset(input_file)

    insitu_time = np.array(input_dataset.variables['insitu_time'])
    satellite_source = np.array(input_dataset.variables['satellite_PDU'])

    satellite_time_in = np.array(input_dataset.variables['satellite_time'])
    satellite_time_out = satellite_time_in.copy()
    time_diff_in = np.array(input_dataset.variables['time_difference'])
    time_diff_out = time_diff_in.copy()

    for idx in range(time_diff_in.shape[0]):
        sat_source_date = None
        sat_source = satellite_source[idx].split('/')[-1]
        sat_source_l = sat_source.split('_')
        for sl in sat_source_l:
            try:
                sat_source_date = dt.strptime(sl, '%Y%m%dT%H%M%S')
            except:
                pass
            if sat_source_date is not None:
                break

        if sat_source_date is None:
            print('f[ERROR] Error getting the satellite date from file name')
            input_dataset.close()
            return

        sat_time_new = sat_source_date.replace(tzinfo=timezone.utc).timestamp()

        satellite_time_out[idx] = np.float32(sat_time_new)

        print(f'[INFO] Sat. source: {sat_source} Sat. time new: {dt.utcfromtimestamp(sat_time_new)}')

        for inx in range(time_diff_in.shape[1]):
            if insitu_time[idx][inx] > 10e+35:
                continue
            if insitu_time[idx][inx] < 0:
                continue
            insitu_time_val = np.float32(insitu_time[idx][inx])
            insitu_time_here = dt.utcfromtimestamp(np.int32(insitu_time[idx][inx]))
            time_diff_out[idx][inx] = np.float32(np.abs(sat_time_new - insitu_time_val))
            time_diff_here = abs((insitu_time_here - sat_source_date).total_seconds())
            diff = abs(time_diff_here - time_diff_out[idx][inx])
            if diff > 0:
                print('[ERROR] It should not be here')

    ncout = Dataset(output_file, 'w', format='NETCDF4')

    # copy global attributes all at once via dictionary
    ncout.setncatts(input_dataset.__dict__)

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    for name, variable in input_dataset.variables.items():
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(input_dataset[name].__dict__)

        # copy data
        if name == 'satellite_time':
            ncout[name][:] = satellite_time_out[:]
        elif name == 'time_difference':
            ncout[name][:] = time_diff_out[:]
        else:
            ncout[name][:] = input_dataset[name][:]

    ncout.close()

    input_dataset.close()


def creating_copy_adding_mu(input_file, ref_file, output_file, idmin, idmax):
    from netCDF4 import Dataset
    import numpy as np
    input_dataset = Dataset(input_file)
    ref_dataset = Dataset(ref_file)
    ncout = Dataset(output_file, 'w', format='NETCDF4')
    nnew = (idmax - idmin) + 1
    norig = len(input_dataset.dimensions['satellite_id'])
    ntotal = norig + nnew
    idmin_new = norig
    idmax_new = (idmin_new + nnew) - 1

    # copy global attributes all at once via dictionary
    ncout.setncatts(ref_dataset.__dict__)

    # copy dimensions
    for name, dimension in ref_dataset.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

        # copy variables
    for name, variable in ref_dataset.variables.items():
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue
        # print(name)
        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(ref_dataset[name].__dict__)
        # copy variable data
        if name == 'insitu_exact_wavelenghts' or name == 'insitu_Rrs':
            array_ref = np.array(ref_dataset[name][:])
            array_add = np.array(input_dataset[name][idmin:idmax, :, 0:30])
            print('----------------------')
            print(name)
            print(array_ref.shape)
            print(array_add.shape)
            array_end = np.concatenate((array_ref, array_add))
            print(array_end.shape)
            print('----------------------')
            ncout[name][:] = array_end[:]
        elif name == 'insitu_time' or name == 'time_difference':
            array_ref = np.array(ref_dataset[name][:])
            array_add = np.array(input_dataset[name][idmin:idmax, 0:30])
            print('----------------------')
            print(name)
            print(array_ref.shape)
            print(array_add.shape)
            array_end = np.concatenate((array_ref, array_add))
            print(array_end.shape)
            print('---------------------')
            ncout[name][:] = array_end[:]
        elif name == 'insitu_original_bands' or name == 'satellite_bands':
            ncout[name][:] = ref_dataset[name][:]
        elif name == 'satellite_latitude' or name == 'satellite_longitude':
            print('====================')
            print(name)
            array_ref = np.array(ref_dataset[name][:])
            array_add = np.array(input_dataset[name][idmin:idmax, :, :])
            print(array_ref.shape)
            print(array_add.shape)
            array_end = np.concatenate((array_ref, array_add))
            print(array_end.shape)
            ncout[name][:] = array_end[:]
            print('====================')
        elif name == 'satellite_time' or name == 'satellite_PDU':
            print('====================')
            print(name)
            array_ref = np.array(ref_dataset[name][:])
            array_add = np.array(input_dataset[name][idmin:idmax])
            print(array_ref.shape)
            print(array_add.shape)
            array_end = np.concatenate((array_ref, array_add))
            print(array_end.shape)
            ncout[name][:] = array_end[:]
            print('====================')
        elif name == 'satellite_Rrs':
            print('====================')
            print(name)
            array_ref = np.array(ref_dataset[name][:])
            array_add = np.array(input_dataset[name][idmin:idmax, :, :, :])
            print(array_ref.shape)
            print(array_add.shape)
            array_end = np.concatenate((array_ref, array_add))
            print(array_end.shape)
            ncout[name][:] = array_end[:]
            print('====================')
        else:
            print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^', name)

    ncout.close()


def creating_copy_with_mu_valid_common(reader, file_out, mu_valid_common):
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

    # creating mu_valid_common variable
    var = ncout.createVariable('mu_valid_common', 'i1', ('satellite_id',), fill_value=-1, zlib=True, shuffle=True,
                               complevel=6)
    var[:] = mu_valid_common[:]

    ncout.close()
    reader.mfile.close()
    return True


def creating_copy_with_flag_bands(reader, file_out, flag_lists, satellite_id_ref, extra_bands):
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
        if name_flag in reader.mfile.nc.variables:
            print(f'[INFO] Flag name {name_flag} is already in the original file')
            continue
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

    angle_bands = ['OZA', 'OAA', 'SZA', 'SAA']
    if len(extra_bands) > 0:
        for name_band in extra_bands:
            if name_band.startswith('satellite_') or name_band in angle_bands:
                fill_value = -999.0
                if '_FillValue' in extra_bands[name_band].keys():
                    fill_value = extra_bands[name_band]['_FillValue']
                var = ncout.createVariable(name_band, 'f4', ('satellite_id', 'rows', 'columns'), zlib=True,
                                           shuffle=True, complevel=6, fill_value=fill_value)
                for at in extra_bands[name_band]:
                    if at != '_FillValue':
                        var.setncattr(at, extra_bands[name_band][at])
                var[:] = fill_value

    ncout.close()
    reader.mfile.close()
    return True


def creating_copy_correcting_band_bis(file_in, file_out, band_to_correct, new_array):
    # reader = MDB_READER('', True)
    from netCDF4 import Dataset
    import numpy as np
    ncin = Dataset(file_in)
    ncout = Dataset(file_out, 'w', format='NETCDF4')

    # copy global attributes all at once via dictionary
    ncout.setncatts(ncin.__dict__)

    # copy dimensions
    for name, dimension in ncin.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    # copy variables
    for name, variable in ncin.variables.items():
        fill_value = None
        if '_FillValue' in list(ncin.ncattrs()):
            fill_value = variable._FillValue
        # print(name)
        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(ncin[name].__dict__)

        # copy variable data
        ncout[name][:] = ncin[name][:]
        if name == 'satellite_Rrs' and band_to_correct is None and new_array is None:
            ncout[name][:, 1, :, :] = ncin[name][:, 1, :, :] * np.pi
        if band_to_correct is not None and name == band_to_correct:
            ncout[name][:] = new_array[:]

    ncout.close()
    ncin.close()
    return True


def creating_copy_correcting_band(reader, file_out, band_to_correct, new_array):
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
        if name == 'satellite_Rrs' and band_to_correct is None and new_array is None:
            ncout[name][:, 1, :, :] = reader.mfile.nc[name][:, 1, :, :] * np.pi
        if band_to_correct is not None and name == band_to_correct:
            ncout[name][:] = new_array[:]

    ncout.close()
    reader.mfile.close()
    return True


def creating_copy_with_new_bands(reader, file_out, wl_new):
    from netCDF4 import Dataset
    import numpy as np
    wl_prev = np.array(reader.mfile.satellite_bands)
    indices = []
    for wl in wl_new:
        diff = np.abs(wl - wl_prev)
        imin = np.argmin(diff)
        if diff[imin] < 1:
            indices.append(imin)
            print(f'[INFO] Wavelength: {wl} Index: {imin} Previous wavelength: {wl_prev[imin]}')
        else:
            indices.append(-1)
            print(f'[INFO] Wavelength: {wl} Index: -1 (not previously available in the file)')

    nwl = len(wl_new)
    ncout = Dataset(file_out, 'w', format='NETCDF4')
    # copy global attributes all at once via dictionary
    if args.verbose:
        print(f'[INFO] Copying attributes...')
    ncout.setncatts(reader.mfile.nc.__dict__)
    # copy dimensions
    if args.verbose:
        print(f'[INFO] Copying dimensions...')
    for name, dimension in reader.mfile.nc.dimensions.items():
        if name == 'satellite_bands':
            ncout.createDimension(name, nwl)
        else:
            ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

    # copy variables
    for name, variable in reader.mfile.nc.variables.items():
        if args.verbose:
            print(f'[INFO] Copying variable {name}')
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(reader.mfile.nc[name].__dict__)

        # copy variable data

        if name == 'satellite_bands':
            ncout[name][:] = [wl_new[:]]
        elif name == 'satellite_Rrs':
            for idx in range(nwl):
                index = indices[idx]
                if index >= 0:
                    ncout[name][:, idx, :, :] = reader.mfile.nc[name][:, index, :, :]
                else:
                    ncout[name][:, idx, :, :] = fill_value
        else:
            ncout[name][:] = reader.mfile.nc[name][:]

    ncout.close()
    reader.mfile.close()


def creating_copy_with_new_insitu(reader, file_out, wl_new):
    from netCDF4 import Dataset
    import numpy as np
    wl_prev = np.array(reader.mfile.insitu_bands)
    indices_new = {}
    for wl in wl_new:
        indices_here = []
        diff = np.abs(wl - wl_prev)
        for idx in range(len(wl_prev)):
            if diff[idx] <= 5:
                indices_here.append(idx)
        indices_new[wl] = indices_here

    print(indices_new)

    nwl = len(wl_new)
    ncout = Dataset(file_out, 'w', format='NETCDF4')
    # copy global attributes all at once via dictionary
    if args.verbose:
        print(f'[INFO] Copying attributes...')
    ncout.setncatts(reader.mfile.nc.__dict__)

    # copy dimensions
    if args.verbose:
        print(f'[INFO] Copying dimensions...')
    for name, dimension in reader.mfile.nc.dimensions.items():
        if name == 'insitu_original_bands':
            ncout.createDimension(name, nwl)
        else:
            ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

    # copy variables
    for name, variable in reader.mfile.nc.variables.items():
        if args.verbose:
            print(f'[INFO] Copying variable {name}')
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(reader.mfile.nc[name].__dict__)

        # copy variable data

        if name == 'insitu_original_bands':
            ncout[name][:] = [wl_new[:]]
        elif name == 'insitu_Rrs' or name == 'insitu_exact_wavelenghts':
            for idx in range(nwl):
                wl_here = wl_new[idx]
                indices_here = indices_new[wl_here]
                # print(name, ':', wl_here, '---->', indices_here, len(indices_here))
                if len(indices_here) == 0:
                    ncout[name][:, idx, :] = fill_value
                elif len(indices_here) == 1:
                    index = indices_here[0]
                    ncout[name][:, idx, :] = reader.mfile.nc[name][:, index, :]
                else:
                    data_here = reader.mfile.nc[name][:, indices_here, 0]
                    ndata = data_here.shape[0]
                    for idata in range(ndata):
                        igood = np.where(~data_here[idata, :].mask)
                        if len(igood[0]) == 0:
                            ncout[name][idata, idx, :] = fill_value
                        else:
                            index = indices_here[igood[0][0]]
                            ncout[name][idata, idx, :] = reader.mfile.nc[name][idata, index, :]

        else:
            ncout[name][:] = reader.mfile.nc[name][:]

    ncout.close()
    reader.mfile.close()


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


def creating_copy_with_valid_indices(reader, file_out, indices_sat_valid):
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

        if name == 'mu_valid':
            mu_valid_array = np.array(reader.mfile.nc[name][:])
            mu_valid_array[indices_sat_valid == False] = 0
            ncout[name][:] = mu_valid_array[:]
        else:
            ncout[name][:] = reader.mfile.nc[name][:]

    ncout.close()
    reader.mfile.close()
    return True


# satellite: S3 or S2
def creating_copy_mdb_publication(file_in, file_out, satellite):
    from netCDF4 import Dataset
    input_dataset = Dataset(file_in)
    ncout = Dataset(file_out, 'w', format='NETCDF4')
    for at in input_dataset.ncattrs():
        val = input_dataset.getncattr(at)
        if val == 'olci':
            val = 'OLCI'
        if at == 'insitu_site_name':
            at = 'site'
        if at == 'insitu_lat':
            at = 'in_situ_lat'
        if at == 'insitu_lon':
            at = 'in_situ_lon'
        if satellite == 'S2' and at == 'satellite_proc_version':
            val = ''
        if at.startswith('satellite_ws'):
            continue
        if at.startswith('satellite_S'):
            continue
        if at.startswith('satellite_V'):
            continue

        ncout.setncattr(at, val)

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

    for name, variable in input_dataset.variables.items():

        if name == 'insitu_filename':
            continue
        if name == 'satellite_PDU':
            continue
        name_new = name
        if name == 'insitu_solar_azimuth_angle':
            name_new = 'insitu_SAA'
        if name == 'insitu_solar_zenith_angle':
            name_new = 'insitu_SZA'
        if name == 'insitu_viewing_azimuth_angle':
            name_new = 'insitu_OAA'
        if name == 'insitu_viewing_zenith_angle':
            name_new = 'insitu_OZA'
        if name == 'OAA':
            name_new = 'satellite_OAA'
        if name == 'OZA':
            name_new = 'satellite_OZA'
        if name == 'SAA':
            name_new = 'satellite_SAA'
        if name == 'SZA':
            name_new = 'satellite_SZA'

        if name == 'satellite_WQSF' and satellite == 'S2':
            name_new = 'satellite_IdePix'

        print('Var: ', name_new)
        datatype = variable.datatype
        if name == 'insitu_time':
            datatype = 'f8'
        if name == 'satellite_WQSF':
            datatype = 'u8'

        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        if name == 'satellite_WQSF':
            ncout.createVariable(name_new, datatype, variable.dimensions, zlib=True, shuffle=True, complevel=6)
        else:
            ncout.createVariable(name_new, datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                                 shuffle=True, complevel=6)

        # copy variable attributes all at once via dictionary
        if name_new == 'insitu_Rrs':
            ncout[name_new].short_name = 'remote_sensing_reflectance'
            ncout[name_new].long_name = 'Remote-Sensing reflectance of the water column at the surface'
            ncout[name_new].units = 'sr-1'
        elif name_new == 'insitu_Rrs_nosc':
            ncout[name_new].short_name = 'remote_sensing_reflectance_nosc'
            ncout[
                name_new].long_name = 'Remote-Sensing reflectance of the water column at the surface without correction for the NIR similarity spectrum (see Ruddick et al., 2006)'
            ncout[name_new].units = 'sr-1'
        else:
            for at in input_dataset[name].ncattrs():
                if at == '_FillValue':
                    continue

                at_new = at
                if at == 'description':
                    at_new = 'long_name'
                if at == 'standard_name':
                    at_new = 'short_name'
                if at == 'flag_mask':
                    at_new = 'flag_masks'
                if at == 'flag_values':
                    at_new = 'flag_masks'
                val = input_dataset[name].getncattr(at)
                ncout[name_new].setncattr(at_new, val)

            if name == 'insitu_original_bands':
                ncout[name_new].units = 'nm'
            if name == 'insitu_time':
                ncout[name_new].long_name = 'In situ time in ISO 8601 format.'
                ncout[name_new].units = 'Seconds since 1970-01-01 00:00:00 UTC'
            if name == 'satellite_time':
                ncout[name_new].long_name = 'Satellite time in ISO 8601 format.'
                ncout[name_new].units = 'Seconds since 1970-01-01 00:00:00 UTC'
            if name == 'satellite_bands':
                ncout[name_new].long_name = 'Satellite bands in nm.'

            if name_new == 'satellite_IdePix':
                ncout[name_new].long_name = 'IdePix Pixel Classification Flag Dataset'

            # ncout[name_new].setncatts(input_dataset[name].__dict__)

        ncout[name_new][:] = input_dataset[name][:]

    ncout.close()
    input_dataset.close()

    return True


def creating_copy_with_convolution(file_in, file_out, array_new):
    from netCDF4 import Dataset
    input_dataset = Dataset(file_in)
    ncout = Dataset(file_out, 'w', format='NETCDF4')

    # copy global attributes all at once via dictionary
    ncout.setncatts(input_dataset.__dict__)

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

    # copy variables
    for name, variable in input_dataset.variables.items():
        fill_value = None
        if '_FillValue' in list(input_dataset.ncattrs()):
            fill_value = variable._FillValue
        # print(name)
        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)

        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(input_dataset[name].__dict__)

        if name == 'mu_ins_rrs':
            ncout[name][:] = array_new[:]
        else:
            ncout[name][:] = input_dataset[name][:]
    input_dataset.close()
    ncout.close()


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

        if cm[time_here_str][sat] and used[time_here_str][sat] < 2:
            indices.append(idx)
            indices_valid[idx] = True
            used[time_here_str][sat] = used[time_here_str][sat] + 1

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

    return np.array(indices_valid, dtype=bool), np.array(indices_mu_common_valid, dtype=bool)


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


def correct_mdb_time(site):
    print('CORRECT MDB')
    dir_mdb = f'/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/{site}/NO_SAT_TIME_CORRECTED'
    dir_out = f'/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/{site}'
    from netCDF4 import Dataset
    import numpy as np
    for name in os.listdir(dir_mdb):
        if name.startswith('MDB_'):
            file_in = os.path.join(dir_mdb, name)
            file_out = os.path.join(dir_out, name)
            creating_copy_correcting_sat_time(file_in, file_out)

    return True


def do_test():
    print('TEST')
    # from netCDF4 import Dataset
    # import numpy as np
    # input_path = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/ALLSITES/MDBrc_S2AB_ALLSITES_MSI_20M_COMMONMU.nc'
    # dataset = Dataset(input_path)
    # common_valid = np.array(dataset.variables['mu_valid_common'])
    # flag_site = np.array(dataset.variables['flag_site'])
    # flag_ac = np.array(dataset.variables['flag_ac'])
    # site_values = np.unique(flag_site).tolist()
    # nsites = len(site_values)
    # ncommon_bysite = np.zeros(nsites)
    # for idx in range(len(common_valid)):
    #     site_value = flag_site[idx]
    #     index = site_values.index(site_value)
    #     if common_valid[idx]==1 and flag_ac[idx]==1:
    #         ncommon_bysite[index] = ncommon_bysite[index]+1
    # dataset.close()
    # print(ncommon_bysite)
    # from netCDF4 import Dataset
    # import numpy as np
    # from datetime import datetime as dt
    # from datetime import timedelta
    # file_data = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/MAFR/MDBrc_S2AB_MAFR_ACOLITE_MSI_20M.nc'
    # # file_data = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/BEFR/MDBrc_S2AB_BEFR_ACOLITE_MSI_20M.nc'
    # dataset = Dataset(file_data)
    # mu_valid = np.array(dataset.variables['mu_valid'])
    # mu_insitu_id = np.array(dataset.variables['mu_insitu_id'])

    # file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/BEFR/BEFR_ACOLITE_EXTRACTS.csv'
    # f1 = open(file_out, 'w')
    # satellite_source = dataset.variables['satellite_PDU']
    # for idx in range(mu_valid.shape[0]):
    #     if mu_valid[idx] == 1:
    #         f1.write(satellite_source[idx].split('/')[-1])
    #         f1.write('\n')
    # f1.close()

    # insitu_source = dataset.variables['insitu_filename']
    # satellite_time = np.array(dataset.variables['satellite_time'])
    # insitu_time = np.array(dataset.variables['insitu_time'])
    # insitu_file_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/INSITU_HYPSTAR/MAFR'
    # file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/MAFR/MAFR_Valid_MU_CORRECTED.csv'
    # fline = 'Date;SatelliteID;InsituID;SatTime;InsituTime;SatelliteSource;InsituSource;TimeDiff'
    # f1 = open(file_out, 'w')
    # f1.write(fline)
    # for idx in range(mu_valid.shape[0]):
    #     if mu_valid[idx] == 1:
    #         print(idx, '--------------------------------')
    #         stime = dt.utcfromtimestamp(np.int32(satellite_time[idx]))
    #         date = stime.strftime('%Y-%m-%d')
    #         insitu_id = mu_insitu_id[idx]
    #         itime = dt.utcfromtimestamp(np.int32(insitu_time[idx][insitu_id]))
    #         ssource = satellite_source[idx].split('/')[-1]
    #         isource = insitu_source[idx][insitu_id]
    #         # fisource = os.path.join(insitu_file_base,itime.strftime('%Y'),itime.strftime('%m'),itime.strftime('%d'),isource)
    #         # print(ssource)
    #         # print(stime)
    #         #
    #         # print(isource)
    #         # print(itime)
    #         #
    #         # print(fisource,os.path.exists(fisource))
    #         #
    #         # dinsitu = Dataset(fisource)
    #         # var = dinsitu.variables['acquisition_time']
    #         # value = int(var[0])
    #         # itimeinsitu = dt(1970,1,1)+timedelta(seconds=value)
    #         # itimeinsitu2 = dt.fromtimestamp(value)
    #         # itimeinsitu3 = dt.utcfromtimestamp(value)
    #         # print(insitu_time[idx][insitu_id],value,itimeinsitu,itimeinsitu2,itimeinsitu3)
    #         # dinsitu.close()
    #
    #         shour = stime.strftime('%H:%M:%S')
    #         ihour = itime.strftime('%H:%M:%S')
    #         tdiff = abs((stime - itime).total_seconds()) / 60
    #
    #         line = f'{date};{idx};{insitu_id};{shour};{ihour};{ssource};{isource};{tdiff}'
    #         f1.write('\n')
    #         f1.write(line)
    # f1.close()
    # dataset.close()

    ##CHECKING OLD AND NEW FLAGS
    from netCDF4 import Dataset
    import numpy as np
    dirbase = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/SAT_EXTRACTS/MAFR/ACOLITEv1'
    dirc = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/SAT_EXTRACTS/MAFR/ACOLITE'
    for name in os.listdir(dirbase):
        file_ardl = os.path.join(dirbase, name)
        file_orig = os.path.join(dirc, name)
        if os.path.exists(file_ardl) and os.path.exists(file_orig):
            dataset_ardl = Dataset(file_ardl)
            dataset_orig = Dataset(file_orig)
            flag_ardl = np.array(dataset_ardl.variables['satellite_WQSF'])
            flag_orig = np.array(dataset_orig.variables['satellite_WQSF'])
            ratio = flag_orig / flag_ardl
            print(np.min(ratio[:]), np.max(ratio[:]), np.average(ratio[:]))
            dataset_ardl.close()
            dataset_orig.close()

    ##CHECKING IF ARDL and IDEPIX are the same
    # if name.endswith('ARDL.nc'):
    #     name_orig = name[:-8]
    #     file_ardl = os.path.join(dirbase,name)
    #     file_orig = os.path.join(dirbase,f'{name_orig}.nc')
    #     if os.path.exists(file_ardl) and os.path.exists(file_orig):
    #         dataset_ardl = Dataset(file_ardl)
    #         dataset_orig = Dataset(file_orig)
    #         flag_ardl = np.array(dataset_ardl.variables['satellite_WQSF'])
    #         flag_orig = np.array(dataset_orig.variables['satellite_WQSF'])
    #         ratio = flag_orig/flag_ardl
    #         print(np.min(ratio[:]),np.max(ratio[:]),np.average(ratio[:]))
    #         dataset_ardl.close()
    #         dataset_orig.close()

    ##CHECKING IF ARDL FLAG ARE THE SAME FOR BOTH AC
    # for name in os.listdir(dirbase):
    # if name.endswith('_acolite.nc'):
    #     name_c2rcc = name.replace('_acolite.nc','_c2rcc.nc')
    #     file_acolite = os.path.join(dirbase,name)
    #     file_c2rcc = os.path.join(dirbase,name_c2rcc)
    #     dataset_acolite = Dataset(file_acolite)
    #     dataset_c2rcc = Dataset(file_c2rcc)
    #     flag_ardl = np.array(dataset_acolite.variables['satellite_WQSF'])
    #     flag_orig = np.array(dataset_c2rcc.variables['satellite_WQSF'])
    #     ratio = flag_orig/flag_ardl
    #     print(np.min(ratio[:]),np.max(ratio[:]),np.average(ratio[:]))
    #     dataset_acolite.close()
    #     dataset_c2rcc.close()

    ##CHECKING IF PIXELS C2RCC WITH VALUE==0 ARE FLAGGED.
    # if name.endswith('_c2rcc.nc'):
    #     file_c2rcc = os.path.join(dirbase,name)
    #     dataset = Dataset(file_c2rcc)
    #     flag = np.array(dataset.variables['satellite_WQSF'])
    #     rrs = np.array(dataset.variables['satellite_Rrs'])
    #     rrs_704 = np.zeros((1,25,25))
    #     rrs_704[0,:,:] = rrs[0,4,:,:]
    #
    #     indices = np.where(rrs_704==0)
    #     flag_values = flag[indices]
    #     print(np.min(flag_values),np.max(flag_values))
    #
    #     dataset.close()

    dirbase = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/SAT_EXTRACTS/extract_MAFR_Rrs_ARDL'
    from netCDF4 import Dataset
    import numpy as np
    values_end = None
    nfiles = 0
    flag_meanings = None
    flag_masks = None
    for name in os.listdir(dirbase):
        if name.endswith('_acolite_ARDL.nc'):
            nfiles = nfiles + 1
            file_acolite = os.path.join(dirbase, name)
            dataset = Dataset(file_acolite)
            flag = np.array(dataset.variables['satellite_WQSF'])
            values = flag.flatten()
            if values_end is None:
                values_end = values
            else:
                values_end = np.concatenate((values_end, values))
            if nfiles == 1:
                satellite_flag = dataset.variables['satellite_WQSF']
                flag_meanings = satellite_flag.flag_meanings
                if isinstance(flag_meanings, list):
                    flag_meanings = ' '.join(flag_meanings)
                flag_masks = satellite_flag.flag_masks.astype('uint64')
            dataset.close()

    import COMMON.Class_Flags_OLCI as flag
    print(flag_masks, flag_meanings)
    flagging = flag.Class_Flags_OLCI(flag_masks, flag_meanings)

    print(values_end.shape, nfiles)
    values_u = np.unique(values_end)
    print(values_u.shape)
    file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/SAT_EXTRACTS/flags_ardl.csv'
    f1 = open(file_out, 'w')

    meanings = flag_meanings.replace(' ', ';')
    first_line = f'Value;N;{meanings}'
    f1.write(first_line)
    second_line = ';'
    for x in flag_masks:
        print(x)
        second_line = f'{second_line};{x}'
    f1.write('\n')
    f1.write(second_line)

    meanings_list = flag_meanings.split(' ')
    for v in values_u:
        nv = np.where(values_end == v)
        nvalues = len(nv[0])
        fvalue = np.array(v).astype(np.uint64)
        print(v, nvalues, fvalue)
        line = f'{fvalue};{nvalues}'
        for idx in range(len(meanings_list)):
            here = [meanings_list[idx]]
            # herev = flag_masks[idx]
            res = flagging.Mask(fvalue, here)
            if res == 0:
                line = f'{line};0'
            else:
                line = f'{line};{nvalues}'
        f1.write('\n')
        f1.write(line)

    f1.close()

    # print(flag_masks,flag_meanings)

    return True


def do_add_new_flag(site):
    dir_base = f'/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/SAT_EXTRACTS/{site}'
    from netCDF4 import Dataset
    import numpy as np
    refs = ['ACOLITE', 'C2RCC']
    for ref in refs:
        dir_orig = os.path.join(dir_base, f'{ref}v1')
        dir_flag = os.path.join(dir_base, f'{ref}_newflag')
        dir_out = os.path.join(dir_base, f'{ref}')
        for name in os.listdir(dir_orig):
            forig = os.path.join(dir_orig, name)
            fout = os.path.join(dir_out, name)
            fnew = os.path.join(dir_flag, f'{name[:-3]}_ARDL.nc')
            dnew = Dataset(fnew)
            var_new = np.array(dnew.variables['satellite_WQSF'])
            print(forig)
            dnew.close()
            creating_copy_adding_new_flag(forig, fout, var_new)
    return True


def do_check_times(site):
    from netCDF4 import Dataset
    import numpy as np
    from datetime import datetime as dt
    dir_base = f'/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/{site}'
    for name in os.listdir(dir_base):
        if name.startswith('MDBr_'):
            fmdb = os.path.join(dir_base, name)
            print(name)
            dataset = Dataset(fmdb)
            satellite_time = np.array(dataset.variables['satellite_time'])
            satellite_PDU = dataset.variables['satellite_PDU']
            insitu_time = np.array(dataset.variables['insitu_time'])
            insitu_file = dataset.variables['insitu_filename']
            time_diff = dataset.variables['time_difference']

            nsat = satellite_time.shape[0]
            ninsitu = insitu_time.shape[1]
            ##satellite
            for idx in range(nsat):
                time_sat = dt.utcfromtimestamp(float(satellite_time[idx]))
                name_sat = satellite_PDU[idx].split('/')[-1]
                time_name_sat = dt.strptime(name_sat.split('_')[4], '%Y%m%dT%H%M%S')
                time_dif = abs((time_sat - time_name_sat).total_seconds())
                if time_dif > 120:
                    print('ERROR:', time_sat, name_sat, time_name_sat, time_dif)
                for ids in range(ninsitu):
                    if insitu_time[idx][ids] > 9e+36:
                        continue
                    time_insitu = dt.utcfromtimestamp(float(insitu_time[idx][ids]))
                    name_insitu = insitu_file[idx][ids]
                    time_name_insitu = dt.strptime(name_insitu.split('_')[5], '%Y%m%dT%H%M')
                    time_dif = abs((time_insitu - time_name_insitu).total_seconds())
                    if time_dif > 600:
                        print('ERROR:', time_insitu, name_insitu, time_name_insitu, time_dif)
                    time_diff_institu_sat = abs(float(satellite_time[idx]) - float(insitu_time[idx][ids]))

                    res = abs(time_diff[idx][ids] - time_diff_institu_sat)
                    if res > 150:
                        print('ERROR Time diff in file:', time_diff[idx][ids], ' checked it: ', time_diff_institu_sat,
                              'res ', res)

            dataset.close()
    return True


def do_check_times_S3(site):
    from netCDF4 import Dataset
    import numpy as np
    from datetime import datetime as dt
    dir_base = f'/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/{site}'
    for name in os.listdir(dir_base):
        if name.startswith('MDBr_'):
            fmdb = os.path.join(dir_base, name)
            print(name)
            dataset = Dataset(fmdb)
            satellite_time = np.array(dataset.variables['satellite_time'])
            satellite_PDU = dataset.variables['satellite_PDU']
            insitu_time = np.array(dataset.variables['insitu_time'])
            insitu_file = dataset.variables['insitu_filename']
            time_diff = dataset.variables['time_difference']

            nsat = satellite_time.shape[0]
            ninsitu = insitu_time.shape[1]
            ##satellite
            for idx in range(nsat):
                time_sat = dt.utcfromtimestamp(float(satellite_time[idx]))
                name_sat = satellite_PDU[idx].split('/')[-1]
                time_name_sat = dt.strptime(name_sat.split('_')[7], '%Y%m%dT%H%M%S')
                time_dif = abs((time_sat - time_name_sat).total_seconds())
                if time_dif > 120:
                    print('ERROR:', time_sat, name_sat, time_name_sat, time_dif)
                for ids in range(ninsitu):
                    if insitu_time[idx][ids] > 9e+36:
                        continue
                    if insitu_time[idx][ids] < 0:
                        continue
                    time_insitu = dt.utcfromtimestamp(float(insitu_time[idx][ids]))
                    name_insitu = insitu_file[idx][ids]
                    time_name_insitu = dt.strptime(name_insitu.split('_')[5], '%Y%m%dT%H%M')
                    time_dif = abs((time_insitu - time_name_insitu).total_seconds())
                    if time_dif > 600:
                        print('ERROR:', time_insitu, name_insitu, time_name_insitu, time_dif)
                    time_diff_institu_sat = abs(float(satellite_time[idx]) - float(insitu_time[idx][ids]))

                    res = abs(time_diff[idx][ids] - time_diff_institu_sat)
                    if res > 150:
                        print('ERROR Time diff in file:', time_diff[idx][ids], ' checked it: ', time_diff_institu_sat,
                              'res ', res)

            dataset.close()
    return True


def set_id_as_invalid(file_in, list_id):
    file_out = os.path.join(os.path.dirname(file_in), 'TempHere.nc')
    from netCDF4 import Dataset
    import numpy as np
    dataset = Dataset(file_in)
    mu_valid = np.array(dataset.variables['mu_valid'])
    mu_valid_common = None
    if 'mu_valid_common' in dataset.variables:
        mu_valid_common = np.array(dataset.variables['mu_valid_common'])
    dataset.close()
    for id in list_id:
        mu_valid[id] = 0
        if mu_valid_common is not None:
            mu_valid_common[id] = 0
    creating_copy_changing_valid_mu(file_in, file_out, mu_valid, mu_valid_common)
    os.rename(file_out, file_in)

    return True


def do_check_protocols(input_csv):
    print(f'[INFO] Starting check protocols from {input_csv}')
    path_output = os.path.dirname(input_csv)
    import pandas as pd
    df = pd.read_csv(input_csv, sep=';')

    # ref_values = df.iloc[:, 0].drop_duplicates().tolist()
    ref_values = ['BEFR', 'VEIT', 'MAFR', 'LPAR', 'GAIT', 'M1BE']
    param = 'R2'
    path_concatenate = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/CONCATENATE'
    # print(ref_values)

    # time_diff_columns = list(range(15, 181, 15))
    # df_temporal = pd.DataFrame(data=None, index=ref_values, columns=time_diff_columns)
    # work for each ref values
    # for ref in ref_values:
    #     df_ref = df[df.iloc[:, 0] == ref]
    #
    #     for time_diff in time_diff_columns:
    #         new_params = {'time_diff': time_diff * 60}
    #         nvalid = 0
    #         for name in os.listdir(path_concatenate):
    #             os.remove(os.path.join(path_concatenate, name))
    #
    #         for index, row in df_ref.iterrows():
    #             input_mdb = row.iat[1]
    #             input_qc = row.iat[2]
    #             n = get_nvalid(input_mdb, input_qc, new_params,path_concatenate)
    #             if n >= 0:
    #                 nvalid = nvalid + n
    #             else:
    #                 print('[INFO] No valid match-ups were found')
    #         if param=='N':
    #             df_temporal.loc[ref].at[time_diff] = nvalid
    #         else:
    #             files_n = []
    #             for name in os.listdir(path_concatenate):
    #                 files_n.append(os.path.join(path_concatenate,name))
    #             if len(files_n)==0:
    #                 print(f'[ERROR] MDBr files were not created....')
    #                 return
    #             if len(files_n)==1:
    #                 file_mdbrc = files_n[0]
    #             else:
    #                 file_mdbrc =  os.path.join(path_concatenate,'MDBrcOutput.nc')
    #                 concatenate_nc_impl(files_n,path_concatenate,file_mdbrc)
    #             from MDBPlotV2 import MDBPlot
    #             mplot = MDBPlot(file_mdbrc)
    #             mplot.set_data_scatterplot(None,None,None,None,None)
    #             mplot.compute_statistics(False,False)
    #             print(mplot.valid_stats)
    #             df_temporal.loc[ref].at[time_diff] = mplot.valid_stats['DETER(r2)']
    #             output_file_csv = os.path.join(path_output, 'TemporalVariationR2.csv')
    #             df_temporal.to_csv(output_file_csv, sep=';')
    #
    # output_file_csv = os.path.join(path_output, 'TemporalVariation.csv')
    # df_temporal.to_csv(output_file_csv, sep=';')

    # min_valid_columns = list(range(1, 10))
    # df_spatial = pd.DataFrame(data=None, index=ref_values, columns=min_valid_columns)
    # # work for each ref values
    # for ref in ref_values:
    #     df_ref = df[df.iloc[:, 0] == ref]
    #     print('------------------------------------------------------------------------------>',ref)
    #     for min_valid in min_valid_columns:
    #         new_params = {'n_min_valid': min_valid}
    #         nvalid = 0
    #
    #         for name in os.listdir(path_concatenate):
    #             os.remove(os.path.join(path_concatenate, name))
    #
    #         for index, row in df_ref.iterrows():
    #             input_mdb = row.iat[1]
    #             input_qc = row.iat[2]
    #             n = get_nvalid(input_mdb, input_qc, new_params,path_concatenate)
    #             if n >= 0:
    #                 nvalid = nvalid + n
    #             else:
    #                 print('[ERROR]')
    #         if param == 'N':
    #             df_spatial.loc[ref].at[min_valid] = nvalid
    #         else:
    #             files_n = []
    #             for name in os.listdir(path_concatenate):
    #                 files_n.append(os.path.join(path_concatenate, name))
    #             if len(files_n) == 0:
    #                 print(f'[ERROR] MDBr files were not created....')
    #                 return
    #             if len(files_n) == 1:
    #                 file_mdbrc = files_n[0]
    #             else:
    #                 file_mdbrc = os.path.join(path_concatenate, 'MDBrcOutput.nc')
    #                 concatenate_nc_impl(files_n, path_concatenate, file_mdbrc)
    #             from MDBPlotV2 import MDBPlot
    #             mplot = MDBPlot(file_mdbrc)
    #             mplot.set_data_scatterplot(None, None, None, None, None)
    #             mplot.compute_statistics(False, False)
    #             print(mplot.valid_stats)
    #             df_spatial.loc[ref].at[min_valid] = mplot.valid_stats['DETER(r2)']
    #             output_file_csv = os.path.join(path_output, 'SpatialVariationR2.csv')
    #             df_spatial.to_csv(output_file_csv, sep=';')

    # output_file_csv = os.path.join(path_output, 'SpatialVariation.csv')
    # df_spatial.to_csv(output_file_csv, sep=';')

    ##COMBINING FIGURES

    output_file_csv = os.path.join(path_output, 'TemporalVariation.csv')
    df_temporal = pd.read_csv(output_file_csv, sep=';')
    output_file_csv = os.path.join(path_output, 'SpatialVariation.csv')
    df_spatial = pd.read_csv(output_file_csv, sep=';')
    output_file_csv = os.path.join(path_output, 'TemporalVariationR2.csv')
    df_temporalr = pd.read_csv(output_file_csv, sep=';')
    output_file_csv = os.path.join(path_output, 'SpatialVariationR2.csv')
    df_spatialr = pd.read_csv(output_file_csv, sep=';')

    file_out_t = os.path.join(path_output, 'TemporalVariation.tif')
    xdefaultst = [120] * 6
    xdefaultss = [9, 9, 1, 9, 1, 9]

    plot_df_lines(df_temporal, file_out_t, 'Maximum time difference (minutes)', '# Valid match-ups', xdefaultst,
                  [0, 250])
    file_out_s = os.path.join(path_output, 'SpatialVariation.tif')
    handles, str_legend = plot_df_lines(df_spatial, file_out_s, 'Mininum number of valid pixels', '# Valid match-ups',
                                        xdefaultss, [0, 300])

    file_out_tr = os.path.join(path_output, 'TemporalVariationR2.tif')
    plot_df_lines(df_temporalr, file_out_tr, 'Maximum time difference (minutes)', f'R$^2$', xdefaultst, [0.3, 1])
    file_out_ts = os.path.join(path_output, 'SpatialVariationR2.tif')
    handles, str_legend = plot_df_lines(df_spatialr, file_out_ts, 'Mininum number of valid pixels', f'R$^2$',
                                        xdefaultss, [0.3, 1])

    file_out = os.path.join('/mnt/c/DATA_LUIS/HYPERNETS_WORK/PUBLICATION/Figures/REVIEW_ROUND_1/Figure13.tif')
    from PlotMultiple import PlotMultiple
    pm = PlotMultiple()
    pm.start_multiple_plot_advanced(2, 2, 12, 9, 0, 0, True)
    pm.plot_image(file_out_t, 0, 0)
    pm.plot_image(file_out_tr, 0, 1)
    pm.plot_image(file_out_s, 1, 0)
    pm.plot_image(file_out_ts, 1, 1)
    pm.set_text(-150, 0, '(a)')
    pm.set_text(1600, 0, '(b)')
    pm.set_text(-150, 1300, '(c)')
    pm.set_text(1600, 1300, '(d)')

    pm.set_global_legend(handles, str_legend)
    pm.save_fig(file_out)

    print('COMPLETED')


def plot_df_lines(df, file_out, xlabel, ylabel, xdefaults, yrange):
    # df = pd.DataFrame()
    sites = 'BEFR,VEIT,MAFR,LPAR,GAIT,M1BE'
    sites = sites.split(',')
    colors = 'blue,red,green,cyan,magenta,orange'
    colors = colors.split(',')

    from PlotSpectra import PlotSpectra
    import numpy as np
    ps = PlotSpectra()
    ps.close_plot()
    ps.start_plot()
    xdata = np.array(df.columns.values[1:], dtype=np.float32)
    ps.xdata = xdata
    # idx = 0
    handles = [''] * len(sites)
    str_values = [''] * len(sites)
    for index, row in df.iterrows():
        y_data = np.array(row.values[1:], dtype=np.float32)
        site = row.values[0]
        idx = sites.index(site)
        color = colors[idx]
        print(site, idx, color, y_data)
        h = ps.plot_single_line(y_data, color, 'solid', 1, 'o', 8)
        handles[idx] = h[0]
        str_values[idx] = site
        if len(xdefaults) == len(sites):
            ihighlight = np.argmin(np.abs(xdefaults[idx] - xdata))
            # print(xdefaults[idx],ihighlight,xdata[ihighlight],y_data[ihighlight],'???????????????????????')
            ps.plot_single_marker(xdata[ihighlight], y_data[ihighlight], 'o', 8, 'w', color, 1)

        # handles.append(h[0])
        # str_values.append(site)

    ps.set_y_range(yrange[0], yrange[1])
    xdata_int = np.array(xdata, dtype=np.int32).tolist()
    ps.set_xticks(xdata, xdata_int, 0, 11)
    if yrange[1] > 1:
        yticks_val = np.array(ps.get_yticks(), dtype=np.int32).tolist()
        ps.set_yticks(ps.get_yticks(), yticks_val, 0, 11)
    ps.set_grid()
    ps.set_xaxis_title_f(xlabel, 12)
    ps.set_yaxis_title_f(ylabel, 12)
    # if xdefault >= 0:
    #     plt.axvline(xdefault, color='black', linewidth=1, ls='-')

    ps.save_fig(file_out)
    ps.close_plot()

    return handles, str_values


def get_nvalid(input_path, input_qc, new_params, path_concatenate):
    print(f'[INFO] Working with file: {input_path} ---------------------------------------------------------------')
    # print(input_path,os.path.exists(input_path))
    reader = MDB_READER(input_path, True)
    if not reader.mfile.VALID:
        return -1
    if not reader.mfile.check_repeated():
        return -1
    # if args.verbose:
    #     print(f'[INFO] Using file: {args.config_file} to set quality control options...')
    import configparser
    from QC_OPTIONS import QC_OPTIONS
    options = configparser.ConfigParser()
    options.read(input_qc)
    qco = QC_OPTIONS(options)
    wllist = qco.get_wllist()
    check_sat_bands = reader.mfile.check_bands(wllist, 5)
    if check_sat_bands == -1:
        return -1

    reader.mfile.set_wl_ref(wllist)
    reader.mfile.qc_sat.ncdataset = reader.mfile.nc
    reader.mfile.qc_sat = qco.get_qcsat(reader.mfile.qc_sat, reader.mfile.nc)

    reader.mfile.qc_insitu.ncdataset = reader.mfile.nc
    reader.mfile.qc_insitu = qco.get_qc_insitu(reader.mfile.qc_insitu)
    if reader.mfile.qc_insitu is None:
        return -1
    if not reader.mfile.qc_insitu.apply_nir_correction:
        reader.mfile.qc_insitu.insitu_rrs = reader.mfile.variables['insitu_Rrs_nosc']
        # if args.verbose:
        #     print('[INFO] NIR Correction is not applied. Using insitu_Rrs_nosc as input variable')

    ##parametres to be modified with respect to quality control file
    if 'time_diff' in new_params:
        reader.mfile.qc_insitu.time_max = new_params['time_diff']

    if 'n_min_valid' in new_params:
        reader.mfile.qc_sat.min_valid_pixels = new_params['n_min_valid']

    # print('---------------------------------------->', reader.mfile.qc_insitu.time_max,reader.mfile.qc_sat.min_valid_pixels)

    if path_concatenate is None:
        nmu_valid, df_valid = reader.mfile.prepare_df_validation()
    else:
        name_output = os.path.basename(input_path).replace('MDB', 'MDBr_TEMPORAL')
        nmu_valid = reader.create_mdb_results_file(os.path.join(os.path.dirname(input_path), name_output))
        os.rename(os.path.join(os.path.dirname(input_path), name_output), os.path.join(path_concatenate, name_output))

    reader.mfile.close()

    return nmu_valid


def plot_distribution_times(date_here, title):
    dir_base = '/mnt/c/DATA_LUIS/TARA_TEST/TIME_DISTRIBUTION'
    date_here = date_here.replace(tzinfo=pytz.utc)
    date_here_str = date_here.strftime('%Y%m%d%H%M')
    file_out = os.path.join(dir_base, f'TimeDistribution_{date_here_str}.tif')
    from MDB_builder.INSITU_tara import INSITU_TARA
    itara = INSITU_TARA(None, True)
    from datetime import timedelta
    date_ini = date_here - timedelta(hours=2)
    date_fin = date_here + timedelta(hours=2)
    date_ini = date_ini.replace(tzinfo=pytz.utc)
    date_fin = date_fin.replace(tzinfo=pytz.utc)
    distributions, summary = itara.get_time_distribution_from_metadata(date_here, date_ini, date_fin)
    # x_data = np.linspace(450,85950,96)
    x_data = np.linspace(0, 95, 96)
    x_ticks = list(range(0, 24))
    x_data_ticks = [0] * 24
    idx = 0
    for hour in range(24):
        x_data_ticks[hour] = x_data[idx]
        idx = idx + 4

    y_data_total = distributions[:, 0]
    y_data_q3 = distributions[:, 1]
    y_data_q4 = distributions[:, 2]
    y_data_q5 = distributions[:, 3]

    bar1 = plt.bar(x_data, y_data_total, color='darkblue', edgecolor='black',
                   linewidth=0, width=1)  # ,marker='o',markersize=10)
    bar2 = plt.bar(x_data, y_data_q3, color='palegreen', edgecolor='black',
                   linewidth=0, width=1)  # , marker='o', markersize=10)
    bar3 = plt.bar(x_data, y_data_q4, color='lime', edgecolor='black',
                   linewidth=0, width=1)  # , marker='o', markersize=10)
    bar4 = plt.bar(x_data, y_data_q5, color='green', edgecolor='black',
                   linewidth=0, width=1)  # , marker='o', markersize=10)
    plt.xticks(x_data_ticks, x_ticks, rotation=90)

    ylim = plt.gca().get_ylim()

    date_ref = date_here.replace(hour=0, minute=0, second=0, tzinfo=pytz.utc)

    x_sat_time = (date_here - date_ref).total_seconds() / 900.0
    x_sat_time_ini = (date_ini - date_ref).total_seconds() / 900.0
    x_sat_time_fin = (date_fin - date_ref).total_seconds() / 900.0
    plt.vlines(x_sat_time, ylim[0], ylim[1], linestyles='-', colors='red', linewidths=1.5)
    plt.vlines(x_sat_time_ini, ylim[0], ylim[1], linestyles='--', colors='red', linewidths=1.5)
    plt.vlines(x_sat_time_fin, ylim[0], ylim[1], linestyles='--', colors='red', linewidths=1.5)
    plt.xlabel('Hour')
    plt.ylabel('Number of spectra')
    plt.legend([bar1, bar2, bar3, bar4],
               [f'Total({int(summary[0])})', f'q_3({int(summary[1])})', f'q_4({int(summary[2])})',
                f'q_5({int(summary[3])})'])
    plt.title(title)
    plt.xlim([20, 68])
    # plt.grid(axis='y')

    # plt.xlim([20000,60000])

    plt.savefig(file_out, dpi=300)
    plt.close()

    return summary


def convert_tara_files():
    import pandas as pd
    from datetime import datetime as dt
    dir_base = '/mnt/c/DATA_LUIS/TARA_TEST/insitu_data/HyperPro-Orig'
    dir_meta = '/mnt/c/DATA_LUIS/TARA_TEST/insitu_data/HyperPro_meta'
    dir_rad = '/mnt/c/DATA_LUIS/TARA_TEST/insitu_data/HyperPro-Rad'

    for name in os.listdir(dir_base):
        if name.endswith('b.csv'):
            continue
        fcsv = os.path.join(dir_base, name)
        df = pd.read_csv(fcsv, sep=',')
        # lat_array = np.array(df['Latitude'])
        # lon_array = np.array(df['Longitude'])
        # day = df['Day']
        # time = df['Time']
        # fmeta = os.path.join(dir_meta,name)
        # fwmeta = open(fmeta,'w')
        # fwmeta.write('lat,lon,timestamp')
        # for idx in range(len(day)):
        #     datetime_str = f'{day[idx]} {time[idx]}'
        #     datetime = dt.strptime(datetime_str,'%Y-%m-%d %H:%M:%S')
        #     datatime_n = datetime.strftime('%Y-%m-%d %H:%M:%S.%f')
        #     line = f'{lat_array[idx]},{lon_array[idx]},{datatime_n}'
        #     fwmeta.write('\n')
        #     fwmeta.write(line)
        # fwmeta.close()
        frad = os.path.join(dir_rad, name)
        fwrad = open(frad, 'w')
        col_names = df.columns.tolist()[4:]
        wl_values = [float(x[4:]) for x in col_names]
        first_line = ",".join([str(x) for x in wl_values])
        first_line = f'#{first_line}'
        fwrad.write(first_line)

        for index, row in df.iterrows():
            rad_values = np.array(row[4:]).tolist()
            rad_values_str = ",".join([str(x) for x in rad_values])
            fwrad.write('\n')
            fwrad.write(rad_values_str)

        fwrad.close()


def apply_spectra_convolution():
    print('started')
    file_nc = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/ALLSITES/MDBrc_S3AB_ALLSITES_OLCI_WFR.nc'

    # file_nc = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/M1BE/MDBrc_S3AB_M1BE_WFR.nc'
    from netCDF4 import Dataset
    from QC_INSITU import QC_INSITU
    dataset = Dataset(file_nc)
    qc_insitu = QC_INSITU(dataset.variables['insitu_Rrs'], dataset.variables['insitu_original_bands'])
    qc_insitu_nosc = QC_INSITU(dataset.variables['insitu_Rrs_nosc'], dataset.variables['insitu_original_bands'])
    sat_wl = np.array(dataset.variables['satellite_bands'])
    qc_insitu.set_wllist_using_wlref(sat_wl)
    qc_insitu_nosc.set_wllist_using_wlref(sat_wl)
    sat_id = -1
    # insitu_id = -1

    sites = ['BEFR', 'GAIT', 'LPAR', 'M1BE', 'MAFR', 'VEIT']
    sats = ['S3A', 'S3B']
    sites_dict = {}
    for site in sites:
        file_site = f'/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/{site}/MDBrc_S3AB_{site}_WFR.nc'
        dataset_site = Dataset(file_site)
        # sat_wl = np.array(dataset_site.variables['satellite_bands'])
        insitu_wl = np.array(dataset_site.variables['insitu_original_bands'])
        sat_insitu_indices = []
        for wl in sat_wl:
            index = np.argmin(np.abs(wl - insitu_wl))
            sat_insitu_indices.append(index)
            # print(wl,'->',index)
        if site == 'LPAR':
            sat_insitu_indices[10] = 809
            sat_insitu_indices[12] = 950
        if site == 'MAFR':
            sat_insitu_indices[10] = 808
            # sat_insitu_indices[12] = 950
        if site == 'M1BE':
            sat_insitu_indices[9] = 753
            sat_insitu_indices[12] = 949
        if site == 'GAIT':
            sat_insitu_indices[9] = 752
        sites_dict[site] = {
            'sat_insitu_indices': sat_insitu_indices,
            'insitu_wl': insitu_wl
        }
        dataset_site.close()

    for site in sites:
        print(site, sites_dict[site]['sat_insitu_indices'])

    array_out = dataset.variables['mu_ins_rrs'][:]
    array_new = array_out.copy()

    for idx, sat_id_here in enumerate(dataset.variables['mu_satellite_id']):
        if idx == 0 or (idx % 100) == 0:
            print('-->', idx)
        insitu_id_here = dataset.variables['mu_insitu_id'][sat_id_here]
        index_site = int(np.log2(dataset.variables['flag_site'][sat_id_here]))
        index_sat = int(np.log2(dataset.variables['flag_satellite'][sat_id_here]))
        site = sites[index_site]
        sat = sats[index_sat]
        qc_insitu.insitu_bands = sites_dict[site]['insitu_wl']
        qc_insitu_nosc.insitu_bands = sites_dict[site]['insitu_wl']
        if sat == 'S3A':
            qc_insitu.srf = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/S3A_OL_SRF_20160713_mean_rsr.nc4'
            qc_insitu_nosc.srf = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/S3A_OL_SRF_20160713_mean_rsr.nc4'
        elif sat == 'S3B':
            qc_insitu.srf = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/S3B_OL_SRF_0_20180109_mean_rsr.nc4'
            qc_insitu_nosc.srf = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/S3B_OL_SRF_0_20180109_mean_rsr.nc4'

        if sat_id_here != sat_id:
            sat_id = sat_id_here
            insitu_id = insitu_id_here
            rrs_values, indices, valid_bands = qc_insitu.get_spectrum_for_mu_and_index_insitu(sat_id, insitu_id)
            rrs_values_nosc, indices_nosc, valid_bands_nosc = qc_insitu_nosc.get_spectrum_for_mu_and_index_insitu(
                sat_id, insitu_id)

        if dataset.variables['mu_valid'][sat_id] == 1:
            wl = dataset.variables['mu_wavelength'][idx]
            iwl = np.argmin(np.abs(wl - sat_wl))
            index = iwl
            # sat_insitu_indices = sites_dict[site]['sat_insitu_indices']
            # index = sat_insitu_indices[iwl]
            if site == 'BEFR' or site == 'VEIT' or site == 'GAIT':
                rrs_here = rrs_values[index]
            else:
                rrs_here = rrs_values_nosc[index]
            # rrs_prev = dataset.variables['mu_ins_rrs'][idx]
            array_new[idx] = rrs_here
            # diff = rrs_here - rrs_prev
            # # if diff > 0.0000001:
            # #     print(site, sat, idx, index, '-->', sat_id, insitu_id, wl, rrs_prev, rrs_here, diff)
            # print(site, sat, idx, index, '-->', sat_id, insitu_id, wl, rrs_prev, rrs_here, diff)

    # file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/ALLSITES/ComparisonWithoutWithConvolution.csv'
    # f1 = open(file_out,'w')
    # f1.write('Withoout;WithConvolution')
    # for idx, rrs_prev in enumerate(dataset.variables['mu_ins_rrs']):
    #     rrs_new = array_new[idx]
    #     line=f'{rrs_prev};{rrs_new}'
    #     f1.write('\n')
    #     f1.write(line)
    # f1.close()
    # diff = abs(rrs_prev-rrs_new)

    dataset.close()

    print('Creating copy...')
    file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/ALLSITES/MDBrc_S3AB_ALLSITES_OLCI_WFR_WITHCONVOLUTION.nc'
    creating_copy_with_convolution(file_nc, file_out, array_new)
    print('Completed')


def do_image_with_centro():
    dir_img = '/mnt/c/DATA_LUIS/INFO_CNR2/JSIT'
    dir_out = '/mnt/c/DATA_LUIS/INFO_CNR2/JSIT/SUN'
    for name in os.listdir(dir_img):
        if name.startswith('Sun'):
            file_img = os.path.join(dir_img, name)
            file_out = os.path.join(dir_out, name)
            # file_img = '/mnt/c/DATA_LUIS/INFO_CNR2/JSIT/01_090_0000_0_0000.jpg'
            # file_out = '/mnt/c/DATA_LUIS/INFO_CNR2/JSIT/01_090_0000_0_0000_grid.png'
            from PIL import Image
            from matplotlib import pyplot as plt
            image = Image.open(file_img)
            # rimage = image.rotate(270, expand=True)
            rimage = image.transpose(Image.ROTATE_270)

            plt.imshow(rimage)

            w, h = rimage.size

            ##central point
            plt.axvline(w / 2, 0.48, 0.52, color='red', linewidth=0.25)
            plt.axhline(h / 2, 0.48, 0.52, color='red', linewidth=0.25)
            ##grid
            incremx = int(w / 4)
            incremy = int(h / 4)
            for x in range(0, w, incremx):
                plt.axvline(x, color='red', linewidth=0.5)
            for y in range(0, h, incremy):
                plt.axhline(y, color='red', linewidth=0.5)

            plt.xticks([])
            plt.yticks([])

            plt.savefig(file_out, bbox_inches='tight', dpi=300)


def get_certo_dates_olci_step1():
    from netCDF4 import Dataset
    from datetime import datetime as dt
    dir_base = '/mnt/c/DATA_LUIS/DOORS_WORK/MDBs'
    certo_dates = {}
    for name in os.listdir(dir_base):
        if name.find('CERTO_OLCI') > 0:
            dataset = Dataset(os.path.join(dir_base, name))
            sat_time = np.array(dataset.variables['satellite_time'])
            for stime in sat_time:
                sat_time_utc = dt.utcfromtimestamp(float(stime))
                key = sat_time_utc.strftime('%Y%m%d')
                certo_dates[key] = {
                    'time': sat_time_utc,
                    'stamp': float(stime)
                }
            dataset.close()

    file_out = os.path.join(dir_base, 'CERTO_OLCI_TIMES.csv')
    fout = open(file_out, 'w')
    fout.write('date;time;stamp')

    for name in os.listdir(dir_base):
        if name.find('CMEMS_OLCI') > 0:
            file_in = os.path.join(dir_base, name)
            dataset = Dataset(file_in)
            sat_time = np.array(dataset.variables['satellite_time'])
            for stime in sat_time:

                sat_time_utc_prev = dt.utcfromtimestamp(float(stime))
                key = sat_time_utc_prev.strftime('%Y%m%d')
                if key in certo_dates:
                    sat_time_utc_new = certo_dates[key]['time']
                    time_stamp = certo_dates[key]['stamp']
                    line = f'{key};{sat_time_utc_new.strftime("%Y-%m-%dT%H:%M:%S.%f")};{time_stamp}'
                else:
                    line = f'{key};NaN;NaN'
                fout.write('\n')
                fout.write(line)

    fout.close()
    dataset.close()


def get_certo_dates_olci_step2():
    from datetime import datetime as dt
    import subprocess
    dir_base = '/store3/DOORS/MDBs'
    dir_sources = '/store3/DOORS/CERTO_SOURCES'
    # dir_base = '/mnt/c/DATA_LUIS/DOORS_WORK/MDBs'
    # dir_sources = '/mnt/c/DATA_LUIS/DOORS_WORK/SOURCES'
    file_in = os.path.join(dir_base, 'CERTO_OLCI_TIMES.csv')
    file_out = os.path.join(dir_base, 'CERTO_OLCI_TIMES_OUT.csv')
    import pandas as pd
    df = pd.read_csv(file_in, sep=';')

    for index, row in df.iterrows():
        there = str(row['time'])
        if there == 'nan':
            dhere = str(row['date'])
            date_here = dt.strptime(dhere, '%Y%m%d')
            if date_here.year < 2017:
                continue

            yyyy = date_here.strftime('%Y')
            jjj = date_here.strftime('%j')
            if date_here.year < 2023:
                namefile = f'CERTO_blk_{dhere}_OLCI_RES300__final_l3_product.nc'
            else:
                namefile = f'CERTO_olci_blk_p2_{dhere}_OLCI_RES300__final_l3_product.nc'
            dir_year = os.path.join(dir_sources, yyyy)
            if not os.path.isdir(dir_year):
                os.mkdir(dir_year)
            dir_jjj = os.path.join(dir_year, jjj)
            if not os.path.isdir(dir_jjj):
                os.mkdir(dir_jjj)
            ofile = os.path.join(dir_jjj, namefile)

            cmd = f'wget --user=rsg_dump --password=yohlooHohw2Pa9ohv1Chi ftp://ftp.rsg.pml.ac.uk/DOORS_matchups/OLCI/{namefile} -O {ofile}'
            print(cmd)
            proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            try:
                outs, errs = proc.communicate(timeout=1800)
            except subprocess.TimeoutExpired:
                proc.kill()
                outs, errs = proc.communicate()
            if os.path.exists(ofile) and os.stat(ofile).st_size > 0:
                sat_time_new = get_satellite_time_from_global_attributes(ofile)
                print(ofile, '->', sat_time_new)
                if sat_time_new is not None:
                    df.loc[index, 'time'] = sat_time_new.strftime("%Y-%m-%dT%H:%M:%S.%f")
                    df.loc[index, 'stamp'] = sat_time_new.timestamp()

    df.to_csv(file_out, sep=';')


def get_certo_dates_msi():
    from datetime import datetime as dt
    import subprocess
    dir_base = '/store3/DOORS/MDBs'
    dir_sources = '/store3/DOORS/CERTO_SOURCES'
    # dir_base = '/mnt/c/DATA_LUIS/DOORS_WORK/MDBs'
    # dir_sources = '/mnt/c/DATA_LUIS/DOORS_WORK/SOURCES'
    file_in = os.path.join(dir_base, 'DOORS_BlackSea_insitu_cnr_iop_extract_CERTO_MSI.csv')
    file_out = os.path.join(dir_base, 'DOORS_BlackSea_insitu_cnr_iop_extract_CERTO_MSI_TIME.csv')
    fout = open(file_out, 'w')
    fout.write('date;stamp')

    import pandas as pd
    df = pd.read_csv(file_in, sep=';')

    for index, row in df.iterrows():
        date_here_str = str(row['Date'])
        date_here = dt.strptime(date_here_str, '%Y-%m-%d')
        yyyy = date_here.strftime('%Y')
        jjj = date_here.strftime('%j')
        orig_file = str(row['OrigFile'])
        if orig_file.lower() == 'nan':
            line = f'{date_here_str};NaN'
            fout.write('\n')
            fout.write(line)
            continue
        dir_year = os.path.join(dir_sources, yyyy)
        dir_jjj = os.path.join(dir_year, jjj)
        file_orig = os.path.join(dir_jjj, orig_file)
        if not os.path.exists(file_orig):
            if not os.path.exists(dir_year):
                os.mkdir(dir_year)
            if not os.path.exists(dir_jjj):
                os.mkdir(dir_jjj)

            cmd = f'wget --user=rsg_dump --password=yohlooHohw2Pa9ohv1Chi ftp://ftp.rsg.pml.ac.uk/DOORS_matchups/MSI/{orig_file} -O {file_orig}'
            print(cmd)
            proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            try:
                outs, errs = proc.communicate(timeout=1800)
            except subprocess.TimeoutExpired:
                proc.kill()
                outs, errs = proc.communicate()
        if os.path.exists(file_orig) and os.stat(file_orig).st_size > 0:
            sat_date = get_satellite_time_from_global_attributes(file_orig)
            if sat_date is not None:
                line = f'{date_here_str};{sat_date.strftime("%Y-%m-%d %H:%M:%S")}'
            else:
                line = f'{date_here_str};SatTimeError'
        else:
            line = f'{date_here_str};NaN'
        fout.write('\n')
        fout.write(line)

    fout.close()


def get_certo_dates_olci():
    from datetime import datetime as dt
    import subprocess
    dir_base = '/store3/DOORS/MDBs'
    dir_sources = '/store3/DOORS/CERTO_SOURCES'
    # dir_base = '/mnt/c/DATA_LUIS/DOORS_WORK/MDBs'
    # dir_sources = '/mnt/c/DATA_LUIS/DOORS_WORK/SOURCES'
    file_in = os.path.join(dir_base, 'DOORS_insitu_from_metadata_11102023_extract_CERTO_OLCI.csv')
    file_out = os.path.join(dir_base, 'DOORS_insitu_from_metadata_11102023_extract_CERTO_OLCI_TIME.csv')
    fout = open(file_out, 'w')
    fout.write('date;stamp')

    import pandas as pd
    df = pd.read_csv(file_in, sep=';')

    for index, row in df.iterrows():
        date_here_str = str(row['Date'])
        date_here = dt.strptime(date_here_str, '%Y-%m-%d')
        yyyy = date_here.strftime('%Y')
        jjj = date_here.strftime('%j')
        index = str(row['Index'])
        if index == -1:
            line = f'{date_here_str};NaN'
            fout.write('\n')
            fout.write(line)
            continue
        dir_year = os.path.join(dir_sources, yyyy)
        dir_jjj = os.path.join(dir_year, jjj)
        if date_here.year == 2023:
            name_file = f'CERTO_blacksea_{date_here.strftime("%Y%m%d")}_OLCI_RES300__final_l3_product.nc'
        else:
            name_file = f'CERTO_blk_{date_here.strftime("%Y%m%d")}_OLCI_RES300__final_l3_product.nc'
        file_orig = os.path.join(dir_jjj, name_file)
        if not os.path.exists(file_orig):
            if not os.path.exists(dir_year):
                os.mkdir(dir_year)
            if not os.path.exists(dir_jjj):
                os.mkdir(dir_jjj)

            cmd = f'wget --user=rsg_dump --password=yohlooHohw2Pa9ohv1Chi ftp://ftp.rsg.pml.ac.uk/DOORS_matchups/OLCI/{name_file} -O {file_orig}'
            print(cmd)
            proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            try:
                outs, errs = proc.communicate(timeout=1800)
            except subprocess.TimeoutExpired:
                proc.kill()
                outs, errs = proc.communicate()
        if os.path.exists(file_orig) and os.stat(file_orig).st_size > 0:
            sat_date = get_satellite_time_from_global_attributes(file_orig)
            if sat_date is not None:
                line = f'{date_here_str};{sat_date.strftime("%Y-%m-%d %H:%M:%S")}'
            else:
                line = f'{date_here_str};SatTimeError'
        else:
            line = f'{date_here_str};NaN'
        fout.write('\n')
        fout.write(line)

    fout.close()


def set_certo_dates_olci_mdb():
    from datetime import datetime as dt
    dir_base = '/mnt/c/DATA_LUIS/DOORS_WORK/MDBs'
    dir_out = '/mnt/c/DATA_LUIS/DOORS_WORK/MDBs/OUT'
    file_new_dates = '/mnt/c/DATA_LUIS/DOORS_WORK/MDBs/CERTO_OLCI_TIMES_OUT.csv'
    df = pd.read_csv(file_new_dates, sep=';')
    info_dates = {}
    for index, row in df.iterrows():
        key = str(row['date'])
        dtime = str(row['time'])
        if dtime == 'nan':
            date_here_str = f'{key}T08:00:00'
            date_here = dt.strptime(date_here_str, '%Y%m%dT%H:%M:%S').replace(tzinfo=pytz.utc).replace(microsecond=0)
            info_dates[key] = {
                'time': date_here,
                'stamp': date_here.timestamp()
            }
        else:
            info_dates[key] = {
                'time': dt.strptime(dtime, "%Y-%m-%dT%H:%M:%S.%f"),
                'stamp': float(row['stamp'])
            }
    from netCDF4 import Dataset
    for name in os.listdir(dir_base):
        if name.find('CMEMS_OLCI') > 0:
            print('-->', name)
            file_in = os.path.join(dir_base, name)
            dataset = Dataset(file_in)
            sat_time = np.array(dataset.variables['satellite_time'])
            sat_time_new = sat_time.copy()
            for idx in range(len(sat_time)):
                date_here = dt.utcfromtimestamp(float(sat_time[idx]))
                key = date_here.strftime('%Y%m%d')
                sat_time_new[idx] = float(info_dates[key]['stamp'])
                print(':->', sat_time[idx], '->', sat_time_new[idx])
            dataset.close()
            output_file = os.path.join(dir_out, name)
            areader = MDB_READER(file_in, True)
            creating_copy_correcting_band(areader, output_file, 'satellite_time', sat_time_new)


def set_certo_dates_msi():
    from datetime import datetime as dt
    dir_base = '/mnt/c/DATA_LUIS/DOORS_WORK/MDBs'
    dir_out = '/mnt/c/DATA_LUIS/DOORS_WORK/MDBs/OUT'
    site = 'Section-7'
    file_new_dates = '/mnt/c/DATA_LUIS/DOORS_WORK/in_situ_extracts/certo_msi/DOORS_insitu_BlackSea_AeronetOC_Section-7_Platform_extract_CERTO_MSI.csv'
    df = pd.read_csv(file_new_dates, sep=';')
    info_dates = {}
    for index, row in df.iterrows():
        sat_time_str = str(row['SatTime'])
        sat_time = dt.strptime(sat_time_str, '%Y-%m-%d %H:%M:%S').replace(microsecond=0).replace(tzinfo=pytz.utc)
        key = sat_time.strftime('%Y%m%d')
        info_dates[key] = {
            'time': sat_time,
            'stamp': sat_time.timestamp()
        }
    from netCDF4 import Dataset
    for name in os.listdir(dir_base):
        if name.find('CERTO_MSI') > 0 and name.find(site) > 0:
            print('-->', name)
            file_in = os.path.join(dir_base, name)
            dataset = Dataset(file_in)
            sat_time = np.array(dataset.variables['satellite_time'])
            sat_time_new = sat_time.copy()
            for idx in range(len(sat_time)):
                date_here = dt.utcfromtimestamp(float(sat_time[idx]))
                key = date_here.strftime('%Y%m%d')
                sat_time_new[idx] = float(info_dates[key]['stamp'])
                print(':->', sat_time[idx], '->', sat_time_new[idx])
            dataset.close()
            output_file = os.path.join(dir_out, name)
            areader = MDB_READER(file_in, True)
            creating_copy_correcting_band(areader, output_file, 'satellite_time', sat_time_new)


def set_certo_dates_extracts():
    from datetime import datetime as dt
    ##olci
    # dir_base = '/mnt/c/DATA_LUIS/DOORS_WORK/extracts_certo_olci'
    # dir_out = '/mnt/c/DATA_LUIS/DOORS_WORK/extracts_certo_olci_out'
    # file_extracts = '/mnt/c/DATA_LUIS/DOORS_WORK/in_situ_extracts/certo_olci/DOORS_BlackSea_insitu_cnr_iop_extract_CERTO_OLCI.csv'
    # msi
    # dir_base = '/mnt/c/DATA_LUIS/DOORS_WORK/extracts_certo_msi'
    # dir_out = '/mnt/c/DATA_LUIS/DOORS_WORK/extracts_certo_msi_out'
    # #file_extracts = '/mnt/c/DATA_LUIS/DOORS_WORK/in_situ_extracts/certo_msi/DOORS_BlackSea_insitu_cnr_iop_extract_CERTO_MSI.csv'
    # file_extracts = '/mnt/c/DATA_LUIS/DOORS_WORK/in_situ_extracts/certo_msi/DOORS_insitu_from_metadata_11102023_extract_CERTO_MSI.csv'
    # cmems olci
    dir_base = '/mnt/c/DATA_LUIS/DOORS_WORK/extracts/cmems_olci'
    dir_out = '/mnt/c/DATA_LUIS/DOORS_WORK/extracts/cmems_olci_out'
    # file_extracts = '/mnt/c/DATA_LUIS/DOORS_WORK/INSITU/in_situ_extracts/cmems/DOORS_insitu_from_metadata_11102023_extract_CMEMS.csv'
    file_extracts = '/mnt/c/DATA_LUIS/DOORS_WORK/INSITU/in_situ_extracts/cmems/DOORS_BlackSea_insitu_cnr_iop_extract_CMEMS.csv'
    df = pd.read_csv(file_extracts, sep=';')
    for index, row in df.iterrows():
        if int(row['IndexExtractOlci']) == -1:
            continue
        name_file = f'extract_{str(row["ExtractOlci"])}'
        input_file = os.path.join(dir_base, name_file)
        output_file = os.path.join(dir_out, name_file)
        sat_time = dt.strptime(str(row['SatTimeOlci']), '%Y-%m-%d %H:%M:%S').replace(tzinfo=pytz.utc)
        ts = sat_time.timestamp()
        if os.path.exists(output_file):
            from netCDF4 import Dataset
            dataset_out = Dataset(output_file)
            ts_out = float(dataset_out.variables['satellite_time'][0])
            dataset_out.close()
            diff = abs(ts - ts_out)
            sat_time_check = dt.utcfromtimestamp(ts_out)
            # print(sat_time.strftime('%Y-%m-%d %H:%M:%S'), sat_time_check.strftime('%Y-%m-%d %H:%M:%S'),diff)
            if diff != 0:
                print(output_file, sat_time.strftime('%Y-%m-%d %H:%M:%S'))
                # fnew = os.path.join('/mnt/c/DATA_LUIS/DOORS_WORK',name_file)
                # os.rename(output_file,fnew)
        else:
            sat_time_check = dt.utcfromtimestamp(ts)
            # print(input_file)
            print(sat_time.strftime('%Y-%m-%d %H:%M:%S'), sat_time_check.strftime('%Y-%m-%d %H:%M:%S'),
                  os.path.exists(input_file))
            array_new = np.array([ts], dtype=np.float64)
            creating_copy_correcting_band_bis(input_file, output_file, 'satellite_time', array_new)


def check_dates():
    from datetime import datetime as dt
    from netCDF4 import Dataset
    dir_base = '/mnt/c/DATA_LUIS/DOORS_WORK/MDBs'
    dir_out = '/mnt/c/DATA_LUIS/DOORS_WORK/MDBs/OUT'
    for name in os.listdir(dir_base):
        if name.find('CMEMS_MULTI') > 0:
            file_in = os.path.join(dir_base, name)
            dataset = Dataset(file_in)
            sat_time = np.array(dataset.variables['satellite_time'])
            ins_time = np.array(dataset.variables['insitu_time'])
            time_diff = np.array(dataset.variables['time_difference'])
            time_diff_new = time_diff.copy()
            for idx in range(len(sat_time)):
                sat_time_here = dt.utcfromtimestamp(float(sat_time[idx]))
                ins_time_here = dt.utcfromtimestamp(float(ins_time[idx][0]))
                time_diff_here = time_diff[idx][0]
                time_diff_c = abs((sat_time_here - ins_time_here).total_seconds())
                time_diff_c2 = abs(sat_time[idx] - ins_time[idx][0])
                time_diff_new[idx, :] = np.abs(sat_time[idx] - ins_time[idx][:])
                print(sat_time_here.strftime('%Y-%m-%d %H:%M:%S'), ins_time_here.strftime('%Y-%m-%d %H:%M:%S'),
                      time_diff_here, time_diff_c, time_diff_c2, time_diff_new[idx, 0])
            dataset.close
            # file_out = os.path.join(dir_out,name)
            # areader = MDB_READER(file_in,True)
            # creating_copy_correcting_band(areader,file_out,'time_difference',time_diff_new)


def get_satellite_time_from_global_attributes(fproduct):
    from netCDF4 import Dataset
    from datetime import datetime as dt
    dataset = Dataset(fproduct)
    if 'time_coverage_start' in dataset.ncattrs() and 'time_coverage_end' in dataset.ncattrs():
        try:
            start_date = dt.strptime(dataset.time_coverage_start, '%d-%b-%Y %H:%M:%S.%f').replace(tzinfo=pytz.utc)
            end_date = dt.strptime(dataset.time_coverage_end, '%d-%b-%Y %H:%M:%S.%f').replace(tzinfo=pytz.utc)
            sat_time = (start_date.timestamp() + end_date.timestamp()) / 2
            sat_time = dt.utcfromtimestamp(sat_time).replace(tzinfo=pytz.utc)
            return sat_time

        except:
            return None
    dataset.close()
    return None


def move_extracts():
    dir_orig = '/mnt/c/DATA_LUIS/DOORS_WORK/temp'
    dir_dest = '/mnt/c/DATA_LUIS/DOORS_WORK/extracts_certo_msi_section-7'
    path_csv = '/mnt/c/DATA_LUIS/DOORS_WORK/in_situ_extracts/certo_msi/DOORS_insitu_BlackSea_AeronetOC_Section-7_Platform_extract_CERTO_MSI.csv'
    df = pd.read_csv(path_csv, sep=';')
    for index, row in df.iterrows():
        name_file = row['Extract']
        name_file = f'extract_{name_file}'
        input_file = os.path.join(dir_orig, name_file)
        output_file = os.path.join(dir_dest, name_file)
        os.rename(input_file, output_file)
    print('DONE')


def do_temporal_owt_comparison(input_path):
    var_data = 'CMEMSVal'
    owt_list = [3, 6, 9, 13]
    file_out = f'/mnt/c/DATA_LUIS/DOORS_WORK/COMPARISON_CMEMS_CERTO/PLOTS/Spectra_{var_data}_4_OWT.tif'
    from PlotSpectra import PlotSpectra
    import pandas as pd
    from CSVPlot import CSVPLOT
    pspectra = PlotSpectra()
    wl_list_str = ['400', '412_5', '442_5', '490', '510', '560', '620', '665', '673_75', '681_25', '708_75', '753_75',
                   '778_75', '865', '885', '1020']
    wl_list = [float(x.replace('_', '.')) for x in wl_list_str]
    wl_ticks = [str(x) for x in wl_list]
    colors = ['Blue', 'Red', 'Green', 'm', 'Cyan', 'Orange', 'Yellow']
    legend = [f'OWT {x}' for x in owt_list]
    df = pd.DataFrame(index=owt_list, columns=wl_list_str)
    options_out = {
        'yvar': var_data,
        'selectBy': 'blended_dominat_owt',
        'selectValue': -1
    }
    for wls in wl_list_str:
        file_csv = os.path.join(input_path, f'rrs_{wls}_points_valid_common.csv')
        cplot = CSVPLOT(file_csv)
        for owt in owt_list:
            options_out['selectValue'] = owt
            median, p25, p75 = cplot.compute_distribution(options_out)
            df.loc[owt].at[wls] = median
    print(df)
    pspectra.xdata = wl_list
    for index, owt in enumerate(owt_list):
        data = np.array(df.loc[owt])
        pspectra.plot_single_line(data, colors[index], '-', 1, 'o', 3)
    pspectra.set_xticks(wl_list, wl_ticks, 45, 8)
    pspectra.set_xaxis_title('Wavelength (nm)')
    pspectra.set_yaxis_title('Rrs')
    pspectra.set_grid()
    pspectra.legend_options['loc'] = 'lower center'
    pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.4)
    pspectra.legend_options['ncols'] = len(owt_list)
    pspectra.set_legend(legend)
    pspectra.set_tigth_layout()

    pspectra.save_plot(file_out)


def do_temporal_spectra_stats_from_csv(input_path):
    do_temporal_spectra_stats_from_csv_impl(input_path, -1)
    for owt in range(1, 18):
        do_temporal_spectra_stats_from_csv_impl(input_path, owt)


def do_temporal_spectra_stats_from_csv_impl(input_path, owt):
    # owt = 13
    file_out = f'/mnt/c/DATA_LUIS/DOORS_WORK/COMPARISON_CMEMS_CERTO/PLOTS/SpectraComparison_OWT_{owt}.tif'
    from PlotSpectra import PlotSpectra
    import pandas as pd
    from CSVPlot import CSVPLOT
    pspectra = PlotSpectra()
    wl_list_str = ['400', '412_5', '442_5', '490', '510', '560', '620', '665', '673_75', '681_25', '708_75', '753_75',
                   '778_75', '865', '885', '1020']
    wl_list = [float(x.replace('_', '.')) for x in wl_list_str]
    wl_ticks = [str(x) for x in wl_list]
    # owt_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
    # owt_list = [3, 6, 9, 13]

    indices = ['median_CMEMS', 'p25_CMEMS', 'p75_CMEMS', 'median_CERTO', 'p25_CERTO', 'p75_CERTO']
    df = pd.DataFrame(index=indices, columns=wl_list_str)
    options_out = {
        'yvar': 'CMEMSVal',
        'selectBy': 'blended_dominat_owt',
        'selectValue': owt
    }
    for wls in wl_list_str:
        file_csv = os.path.join(input_path, f'rrs_{wls}_points_valid_common.csv')
        cplot = CSVPLOT(file_csv)
        options_out['yvar'] = 'CMEMSVal'
        median, p25, p75 = cplot.compute_distribution(options_out)
        df.loc['median_CMEMS'].at[wls] = median
        df.loc['p25_CMEMS'].at[wls] = p25
        df.loc['p75_CMEMS'].at[wls] = p75
        options_out['yvar'] = 'CERTOVal'
        median, p25, p75 = cplot.compute_distribution(options_out)
        df.loc['median_CERTO'].at[wls] = median
        df.loc['p25_CERTO'].at[wls] = p25
        df.loc['p75_CERTO'].at[wls] = p75

    print(df)
    pspectra.xdata = np.array(wl_list)
    pspectra.set_iqr_as_stats_plot()
    pspectra.set_color('red')
    stats = {
        'median': np.array(df.loc['median_CMEMS']).astype(np.float32),
        'p25': np.array(df.loc['p25_CMEMS']).astype(np.float32),
        'p75': np.array(df.loc['p75_CMEMS']).astype(np.float32)
    }
    h1 = pspectra.plot_stats(stats, 0, 16)
    pspectra.set_color('blue')
    stats = {
        'median': np.array(df.loc['median_CERTO']).astype(np.float32),
        'p25': np.array(df.loc['p25_CERTO']).astype(np.float32),
        'p75': np.array(df.loc['p75_CERTO']).astype(np.float32)
    }
    h2 = pspectra.plot_stats(stats, 0, 16)
    pspectra.set_xticks(wl_list, wl_ticks, 45, 8)
    pspectra.set_xaxis_title('Wavelenght(nm)')
    pspectra.set_yaxis_title('Rrs')
    pspectra.legend_options['loc'] = 'lower center'
    pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.4)
    pspectra.legend_options['ncols'] = 2
    pspectra.set_legend_h([h1[0], h2[0]], ['CMEMS-OLCI', 'CERTO-OLCI'])
    pspectra.set_title(f'Optical water type {owt}')
    pspectra.set_grid()
    pspectra.set_tigth_layout()
    pspectra.save_plot(file_out)


def do_temporal_wl_stats_from_csv(type, input_path):
    stats = ['N', 'slope', 'intercept', 'BIAS', 'RMSD', 'DETER(r2)', 'RPD', 'APD']
    ytitles = ['N', 'SLOPE TYPE-II REGRESSION', 'OFFSET TYPE-II REGRESSION', 'bias', 'RMSD', f'R$^2$', 'RPD', 'APD']
    for stat, ytitle in zip(stats, ytitles):
        do_temporal_wl_stats_from_csv_impl(type, input_path, stat, ytitle)


def do_temporal_wl_stats_from_csv_impl(type, input_path, stat, ytitle):
    # N, slope, intercept, BIAS, RMSD, DETER(r2)
    # stat = 'DETER(r2)'
    # f'R$^2$'
    # ytitle = f'R$^2$'
    file_out = f'/mnt/c/DATA_LUIS/DOORS_WORK/COMPARISON_CMEMS_CERTO/PLOTS/Stats_{stat}_{type}.tif'
    from PlotSpectra import PlotSpectra
    import pandas as pd
    from CSVPlot import CSVPLOT
    pspectra = PlotSpectra()
    wl_list_str = ['400', '412_5', '442_5', '490', '510', '560', '620', '665', '673_75', '681_25', '708_75', '753_75',
                   '778_75', '865', '885', '1020']
    wl_list = [float(x.replace('_', '.')) for x in wl_list_str]
    wl_ticks = [str(x) for x in wl_list]
    # owt_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
    if type == '4owt':
        owt_list = [3, 6, 9, 13]
    if type == '6owt':
        owt_list = [1, 2, 3, 6, 9, 13]

    colors = ['Blue', 'Red', 'Green', 'm', 'Cyan', 'Orange', 'Yellow']
    legend = [f'OWT {x}' for x in owt_list]
    df = pd.DataFrame(index=owt_list, columns=wl_list_str)
    pspectra.xdata = wl_list
    options_out = {
        'xvar': 'CMEMSVal',
        'yvar': 'CERTOVal',
        'selectBy': 'blended_dominat_owt',
        'selectValue': 0,
        'log_scale': False
    }
    for wls in wl_list_str:
        file_csv = os.path.join(input_path, f'rrs_{wls}_points_valid_common.csv')
        cplot = CSVPLOT(file_csv)
        for owt in owt_list:
            options_out['selectValue'] = owt
            print(wls, owt)
            valid_stats = cplot.compute_statistics(options_out)
            df.loc[owt].at[wls] = valid_stats[stat]
    print(df)

    for index, owt in enumerate(owt_list):
        data = np.array(df.loc[owt])
        pspectra.plot_single_line(data, colors[index], '-', 1, 'o', 6)
    pspectra.set_xticks(wl_list, wl_ticks, 45, 8)
    if stat == 'APD':
        pspectra.set_y_range(0, 100)
    if stat == 'RPD':
        pspectra.set_y_range(-100, 100)
    pspectra.set_xaxis_title('Wavelength (nm)')
    pspectra.set_yaxis_title(ytitle)
    pspectra.set_grid()
    pspectra.set_legend(legend)
    pspectra.set_tigth_layout()

    pspectra.save_plot(file_out)


def do_temporal_all_owt_stats(type, input_path):
    ## type: CHL or TSM
    stats = ['N', 'slope', 'intercept', 'BIAS', 'RMSD', 'DETER(r2)', 'RPD', 'APD']
    ytitles = ['N', 'SLOPE TYPE-II REGRESSION', 'OFFSET TYPE-II REGRESSION', 'bias', 'RMSD', f'R$^2$', 'RPD', 'APD']
    for stat, ytitle in zip(stats, ytitles):
        do_temporal_all_owt_stats_impl(input_path, type, stat, ytitle)


def do_temporal_all_owt_stats_impl(input_path, type, stat, ytitle):
    # N, slope, intercept, BIAS, RMSD, DETER(r2)
    import pandas as pd
    from PlotSpectra import PlotSpectra
    from CSVPlot import CSVPLOT
    param = 'Chla'
    if type == 'TSM':
        param = 'tsm'
    file_out = f'/mnt/c/DATA_LUIS/DOORS_WORK/COMPARISON_CMEMS_CERTO/PLOTS/Stats_{param}_{stat}.tif'

    cplot = CSVPLOT(input_path)
    options_out = {
        'xvar': 'CMEMS_chl',
        'yvar': 'CERTO_blended_chla',
        'selectBy': 'blended_dominant_owt',
        'selectValue': 0,
        'log_scale': True
    }

    owt = list(range(18))
    if type == 'CHL':
        indices = ['CERTO_blended_chla', 'CERTO_blended_chla_from_predominant_owt', 'CERTO_blended_chla_top_2_weighted',
                   'CERTO_blended_chla_top_3_weighted']
    if type == 'TSM':
        options_out['xvar'] = 'CMEMS_tsmnn'
        indices = ['CERTO_blended_tsm', 'CERTO_blended_tsm_top_2_weighted', 'CERTO_blended_tsm_top_3_weighted']

    colors = ['Blue', 'Red', 'Green', 'm', 'Cyan', 'Orange', 'Yellow']
    df = pd.DataFrame(index=indices, columns=owt)
    for yvar in indices:
        options_out['yvar'] = yvar
        for ow in owt:
            selectValue = ow
            if ow == 0:
                selectValue = -1
            options_out['selectValue'] = selectValue
            valid_stats = cplot.compute_statistics(options_out)
            df.loc[yvar].at[ow] = valid_stats[stat]

    pspectra = PlotSpectra()
    pspectra.xdata = owt
    for index, yvar in enumerate(indices):
        data = np.array(df.loc[yvar])
        pspectra.plot_single_line(data, colors[index], '-', 1, 'o', 5)
    owt_ticks = []
    for ow in owt:
        if ow == 0:
            tick = 'NO OWT'
        else:
            tick = f'OWT {ow}'
        owt_ticks.append(tick)
    pspectra.set_xticks(owt, owt_ticks, 90, 8)
    pspectra.set_xaxis_title('Optical Water Types')
    pspectra.set_yaxis_title(ytitle)
    pspectra.set_grid()
    legend = ['Blended', 'OWT', 'Top-2', 'Top-3']
    if type == 'TSM':
        legend = ['Blended', 'Top-2', 'Top-3']
    pspectra.legend_options['loc'] = 'lower center'
    pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.4)
    pspectra.legend_options['ncols'] = len(legend)
    pspectra.set_legend(legend)
    pspectra.set_tigth_layout()

    pspectra.save_plot(file_out)

    print(df)

def check_n_values_cmems_certo():
    file_nc = '/mnt/c/DATA_LUIS/DOORS_WORK/COMPARISON_CMEMS_CERTO/Coverage_CMEMS_CERTO.nc'
    file_out = '/mnt/c/DATA_LUIS/DOORS_WORK/SummaryCoverage.csv'
    certo_bands = ['400','412','443','490','510','560','620','665','674','681','709','754','779','865','885','1020','blended_chla_from_predominant_owt','blended_chla','blended_chla_top_2_weighted','blended_chla_top_3_weighted']
    cmems_bands = ['400','412_5','442_5','490','510','560','620','665','673_75','681_25','708_75','753_75','778_75','865','885','1020','chl','chl','chl','chl']
    fout = open(file_out,'w')
    fout.write('CERTO_Ref;CMEMS_Ref;N_CERTO;N_CMEMS')


    from netCDF4 import Dataset
    dataset = Dataset(file_nc)
    for idx in range(len(certo_bands)):
        certo_var = f'CERTO_{certo_bands[idx]}_N'
        cmems_var = f'CMEMS_{cmems_bands[idx]}_N'
        ncerto = np.ma.sum(dataset.variables[certo_var][:])
        ncmems = np.ma.sum(dataset.variables[cmems_var][:])
        line = f'{certo_bands[idx]};{cmems_bands[idx]};{ncerto};{ncmems}'
        fout.write('\n')
        fout.write(line)
    dataset.close()
    fout.close()


def main():
    mode = args.mode
    print(f'Started MDBReader with mode: {mode}')

    if args.mode == 'TEST':
        check_n_values_cmems_certo()
        # do_image_with_centro()

        # get_certo_dates_olci_step1()
        # get_certo_dates_olci_step2()
        # set_certo_dates_olci_mdb()
        # move_extracts()
        # set_certo_dates_msi()
        # get_certo_dates_msi()
        # get_certo_dates_olci()
        # check_dates()
        # set_certo_dates_extracts()

        # from BSC_QAA import bsc_qaa_EUMETSAT as qaa
        # import MDBFile
        # code_home = os.path.dirname(os.path.dirname(MDBFile.__file__))
        # # # code_home = os.path.abspath('../')
        # sys.path.append(code_home)
        # code_aeronet = os.path.join(os.path.dirname(code_home), 'aeronet')
        # sys.path.append(code_aeronet)
        # print(code_aeronet)
        # from base.anet_nc_reader import AERONETReader
        # aeronet_file = '/mnt/c/DATA_LUIS/AERONET_OC/AERONET_NC/20020101_20231111_Gloria.LWN_lev20_15.nc'
        # areader = AERONETReader(aeronet_file)
        # all_rrs = areader.extract_rrs(False)
        # all_wl = areader.extract_spectral_data('Exact_Wavelengths',False)
        #
        # print(all_rrs.shape,all_wl.shape)
        #
        # rrs_input = all_rrs[2325,:]
        # wl_input = all_wl[2325,:]
        # print(rrs_input)
        # print(wl_input)
        # rrs_in = [0.004897506441920996, 0.006124403327703476, 0.009032594971358776, 0.011597496457397938,
        #           0.01306516770273447, 0.0037191512528806925, 0.00033687229733914137, 0.00010173415648750961]
        # bands_in = [411.29998779296875, 442.0, 491.3999938964844, 530.7000122070312, 551.9000244140625, 668.5,
        #          869.4000244140625, 1017.2999877929688]
        # print(len(rrs_in),len(bands_in))
        # nominal_wl = areader.get_all_nominal_wl()
        # print(nominal_wl)
        #
        # ##for a single band out
        # bands_out = [412]
        # rrs_out = qaa.bsc_qaa(rrs_in, bands_in, bands_out)
        # print(f'Original wavelength: {bands_in[0]} Output wavelength: {bands_out[0]} Wl diff: {abs(bands_in[0]-bands_out[0])}')
        # print(f'Original rrs: {rrs_in[0]} Output rrs: {rrs_out[0]}')
        #
        # #for CMEMS multi bands
        # bands_out = [412,443,490,510,555,670]
        # rrs_out = qaa.bsc_qaa(rrs_in, bands_in, bands_out)
        # for iwl,wl in enumerate(bands_out):
        #     band_input_ref = bands_in[np.argmin(np.abs(np.array(bands_in)-wl))]
        #     diff_wl = abs(wl-band_input_ref)
        #     if diff_wl<=5:
        #         print(f'Band shifting from {band_input_ref} nm to {wl} nm. Rrs: {rrs_in[iwl]} -> {rrs_out[iwl]}')
        #     else:
        #         print(f'Band shifting from {band_input_ref} nm to {wl} nm. Rrs: {rrs_in[iwl]} -> {rrs_out[iwl]} WARNING: {diff_wl} is higher than 5 nm. It should not be used for validation')

        # print(rrs_out)
        # convert_tara_files()

        ##convert mdb format
        # satellite = 'S2'
        # path_in = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/PUBLICATION/proof/ZenodoMDB'
        # path_output = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/PUBLICATION/proof/ZenodoMDB_final'
        # for name in os.listdir(path_in):
        #     file_in = os.path.join(path_in, name)
        #     file_out = os.path.join(path_output, name)
        #     if name.startswith(f'MDB_{satellite}'):
        #         print('---------------> ', name)
        #         creating_copy_mdb_publication(file_in, file_out, 'S2')

        # from datetime import datetime as dt
        # date_tal = dt(2023,4,5,10,18)
        # file_csv = '/mnt/c/DATA_LUIS/TARA_TEST/MDBs/SatelliteDates.csv'
        # file_out = '/mnt/c/DATA_LUIS/TARA_TEST/MDBs/SatelliteDates_out.csv'
        # f1 = open(file_csv, 'r')
        # fw = open(file_out, 'w')
        # for line in f1:
        #     if line.strip().startswith('Satellite'):
        #         line_out = f'{line.strip()};NTotal;N_q_1;N_q_2;N_q_3'
        #         fw.write(line_out)
        #         continue
        #     values = [x.strip() for x in line.split(';')]
        #     date_here = dt.strptime(values[1], '%Y-%m-%d %H:%M')
        #     title = f'{values[0]} {values[1]}'
        #     summary = plot_distribution_times(date_here, title)
        #     line_out = f'{line.strip()};{int(summary[0])};{int(summary[1])};{int(summary[2])};{int(summary[3])}'
        #     fw.write('\n')
        #     fw.write(line_out)
        # f1.close()
        # fw.close()

        # input_file = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/publication/revision2/Figure7.jpg'
        # output_file = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/publication/production/Figure3.tif'
        # from matplotlib import image as img
        # from matplotlib import pyplot as plt
        # image = img.imread(input_file)
        # plt.imshow(image)
        # plt.axis(False)
        # plt.savefig(output_file, dpi=350, bbox_inches='tight', pil_kwargs={"compression": "tiff_lzw"})
        #
        # from BSC_QAA import bsc_qaa_EUMETSAT as qaa
        # import numpy as np
        # bands: B1	B2	B3	B4
        # bands_in = [442.7,492.4,559.8,664.6]
        # rrs_in  = [0.002348428, 0.003592771, 0.004582189, 0.001049218]
        # band_out = [570]
        # rrsout = qaa.bsc_qaa(rrs_in,bands_in,band_out)
        # rrsoutpolynominal = 0.000178 + (0.734 * rrs_in[2]) + (63.7 * (rrs_in[2]**2)) - (3.08e3 * (rrs_in[2]**3))
        # print(f'rrs at 560 nm: {rrs_in[2]}')
        # print(f'rrs at 570 nm using band-shifting: {rrsout[0]}')
        # print(f'rrs at 570 nm using polynomial: {rrsoutpolynominal}')

        # file_in = '/mnt/c/DATA_LUIS/S3_VALIDATION_TEAM_MEETING_2023/VEIT/MDB_S3B_OLCI_WFR_STANDARD_20230315T000000_20231115T235959_HYPSTAR_VEIT.nc'
        # creating_copy_correcting_attributes(file_in)

        # file_in = '/mnt/c/DATA_LUIS/S3_VALIDATION_TEAM_MEETING_2023/VEIT/MDBrc__S3AB_OLCI_WFR_STANDARD_20230315T000000_20231115T235959_HYPSTAR_VEIT.nc'
        # file_in = '/mnt/c/DATA_LUIS/S3_VALIDATION_TEAM_MEETING_2023/VEIT/MDBrc__S3AB_OLCI_WFR_STANDARD_HYPSTARv1_VEIT.nc'
        # name_var = 'flag_hypstar'
        # flag_values = [1,2]
        # flag_meanings = ['HYPSTARv1','HYPSTARv3']
        # creating_copy_with_flag_band(file_in,name_var,flag_values,flag_meanings,1)
        # update_mdb_file_version(file_in)
        # if do_check_times_S3('LPAR'):
        #     return

        # input_path = '/mnt/c/DATA_LUIS/S3_VALIDATION_TEAM_MEETING_2023/GAIT/extracts'
        # for name in os.listdir(input_path):
        #     input_file = os.path.join(input_path,name)
        #     print(input_file)
        #     update_sat_extract_file_version(input_file)

        # if do_add_new_flag('MAFR'):
        #     return
        # if correct_mdb_time('MAFR'):
        #     return
        # if do_test():
        #     return
        # file = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/ALLSITES/MDBrc_S2AB_ALLSITES_MSI_20M_COMMONMU.nc'
        # list_id = [771]
        # if set_id_as_invalid(file, list_id):
        #     return

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

        # list_dates, list_sat, cm = getting_common_matchups()
        # na = 0
        # nb = 0
        # for datestr in cm:
        #     if cm[datestr]['A']:
        #         na = na + 1
        #     if cm[datestr]['B']:
        #         nb = nb + 1
        # ncommon = na + nb
        # print(na, nb, ncommon)
        #
        # dir_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI'
        # file_in = os.path.join(dir_base, 'CONTATENATED_WITHOUT_COMMON_MATCHUPS', 'MDBrc_ALL.nc')
        # file_out = os.path.join(dir_base, 'MRDrc_ALL.nc')
        # indices, indices_mu = getting_common_satellite_id(file_in, cm)
        #
        # import numpy as np
        # print(np.sum(indices))
        # # print(np.sum(indices_mu))
        # reader = MDB_READER(file_in, True)
        # creating_copy_with_valid_indices(reader, file_out, indices)

        # comparison site
        # from PlotMultiple import PlotMultiple
        # # file_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/COMPARISON_OLCI_MULTI/PLOTS'
        # # format = 'jpg'
        # file_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/publication/FIGURES/'
        # format = 'png'
        # file00 = os.path.join(file_base, f'scatter_plot_443.{format}')
        # file01 = os.path.join(file_base, f'scatter_plot_490.{format}')
        # file02 = os.path.join(file_base, f'scatter_plot_510.{format}')
        # file10 = os.path.join(file_base, f'scatter_plot_560.{format}')
        # file11 = os.path.join(file_base, f'scatter_plot_670.{format}')
        # file12 = os.path.join(file_base, f'scatter_plot_chla.{format}')
        # file_out = os.path.join(file_base, 'Figure14.png')
        # pm = PlotMultiple()
        # pm.start_multiple_plot_advanced(2, 3, 12, 7.5, 0, 0, True)
        # pm.plot_image(file00, 0, 0)
        # pm.plot_image(file01, 0, 1)
        # pm.plot_image(file02, 0, 2)
        # pm.plot_image(file10, 1, 0)
        # pm.plot_image(file11, 1, 1)
        # pm.plot_image(file12, 1, 2)
        # # pm.set_text('a. VEIT ACOLITE',-400,-40)
        # # pm.set_text('b. VEIT C2RCC', 1500, -40)
        # # pm.set_text('c. BEFR ACOLITE', -400, 1385)
        # # pm.set_text('d. BEFR C2RCC', 1500, 1385)
        # pm.save_fig(file_out)
        # pm.close_plot()

        # chla
        # from PlotMultiple import PlotMultiple
        # # file_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/PUBLICATION'
        # # file00 = os.path.join(file_base, 'scatter_plot_chla_cci_all.jpg')
        # # file01 = os.path.join(file_base, 'scatter_plot_chla_cci_common.jpg')
        # # file02 = os.path.join(file_base, 'scatter_plot_chla_olci_common.jpg')
        # file_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/publication/FIGURES'
        # file00 = os.path.join(file_base, 'Figure13a.png')
        # file01 = os.path.join(file_base, 'Figure13b.png')
        # file02 = os.path.join(file_base, 'Figure13c.png')
        # file_out = os.path.join(file_base, 'Figure13.png')
        # pm = PlotMultiple()
        # pm.start_multiple_plot_advanced(1,3,12,7.5,0,0,False)
        # pm.plot_image(file00, 0, 0)
        # pm.plot_image(file01, 0, 1)
        # pm.plot_image(file02, 0, 2)
        # pm.set_text(-1670, 1350, '(a)')
        # pm.set_text(-200,1350,'(b)')
        # pm.set_text(1300, 1350, '(c)')
        # pm.save_fig(file_out)
        # pm.close_plot()

        # added MU DATA
        # input_file = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/CCIv6/MDB_2016-2022/MDB_MULTI_MULTI_1KM_CCI_20180724T000000_20220923T000000_AERONET_Irbe_Lighthouse.nc'
        # ref_file = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/CCIv6/MDB_allperiod/NO2022/MDB_MULTI_MULTI_1KM_CCIALL_20050101T000000_20220907T000000_AERONET_Irbe_Lighthouse.nc'
        # output_file = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/CCIv6/MDB_allperiod/MDB_MULTI_MULTI_1KM_CCIALL_20050101T000000_20220907T000000_AERONET_Irbe_Lighthouse.nc'
        # # GUSTAV:
        # # idini = 296
        # # idfin = 377
        # # IRBE:
        # idini = 166
        # idfin = 244
        #
        # creating_copy_adding_mu(input_file,ref_file,output_file,idini,idfin)

        # check temporal distribution of points
        # input_file = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/MDBrc_OLCI_AERONET_COMMONMU.nc'
        # #reader  = MDB_READER(input_file,True)
        # file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/MUVALID_SITE.jpg'
        # from MDBPlotV2 import MDBPlot
        # mplot = MDBPlot(input_file)
        # mplot.mrfile.analyse_mu_temporal_flag(True, 'flag_ac',file_out)

        # compasions check
        # from datetime import datetime as dt
        # from datetime import timedelta
        # import pandas as pd
        # import numpy as np
        # from scipy import stats
        # ref = 'RRS443'
        # dir_base = f'/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/COMPARISON_OLCI_MULTI/{ref}'
        # file_out = f'/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/COMPARISON_OLCI_MULTI/stats_{ref.lower()}.csv'
        # first_line = 'date;r2'
        # f1 = open(file_out,'w')
        # f1.write(first_line)
        # start_date = dt(2016,5,1)
        # end_date = dt(2022,12,31)
        # date_here = start_date
        # while date_here<=end_date:
        #     date_here_str = date_here.strftime('%Y%j')
        #     date_here_s = date_here.strftime('%Y-%m-%d')
        #     name = f'Comparison_{ref}_{date_here_str}.csv'
        #     fname = os.path.join(dir_base,name)
        #     if os.path.exists(fname):
        #         print(date_here)
        #         df = pd.read_csv(fname,';')
        #         multi_vals = np.array(df['MultiVal'])
        #         olci_vals = np.array(df['OlciVal'])
        #         if len(multi_vals)>0:
        #             slope, intercept, r_value, p_value, std_err = stats.linregress(olci_vals, multi_vals)
        #             r2 = r_value * r_value
        #             line = f'{date_here_s};{r2}'
        #             f1.write('\n')
        #             f1.write(line)
        #     date_here = date_here + timedelta(hours=240)
        # f1.close()

        # list_dates = []
        # for name in os.listdir(dir_base):
        #     datestr = name.split('_')[-1][:-4]
        #     datehere = dt.strptime(datestr,'%Y%j')
        #     #print(datehere)
        #     fcsv = os.path.join(dir_base,name)
        #     dataset = pd.read_csv(fcsv,sep=';')
        #     n_data = 0
        #     ntotal  = 0
        #     for index,row in dataset.iterrows():
        #         val_multi = row['MultiVal']*1000
        #         val_olci = row['OlciVal']*1000
        #         ntotal = ntotal +1
        #         if val_multi<=5 and val_olci>=10:
        #            n_data = n_data +1
        #     if n_data>100:
        #         list_dates.append(datehere)
        #         p = (n_data/ntotal)*100
        #         print(datehere,n_data, p)
        # print(len(list_dates))

        # import pandas as pd
        # dir_base = f'/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/COMPARISON_OLCI_MULTI/points_file_end'
        # print(os.path.exists(dir_base))
        # #file_chla = os.path.join(dir_base,'chla_points.csv')
        # file_rrs = os.path.join(dir_base, 'rrs670_points.csv')
        # #file_output =os.path.join(dir_base,'chla_error.csv')
        # file_output = os.path.join(dir_base, 'rrs670_error.csv')
        # f1 = open(file_output,'w')
        # f1.write('ERROR')
        # #df = pd.read_csv(file_chla,sep=';')
        # df = pd.read_csv(file_rrs, sep=';')
        # for index,row in df.iterrows():
        #     val = '0'
        #     # if row['OlciVal']<0.5 and (row['MultiVal']>=1.40 and row['MultiVal']<=30):
        #     #     val = '1'
        #     if row['OlciVal']<=0.001 and row['MultiVal']>=0.0025:
        #         val = '1'
        #     f1.write('\n')
        #     f1.write(val)
        #
        # f1.close()

        # file_mdb = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/MDBrc_MULTI_AERONET_L3_WITHCCIALL_COMMONMU.nc'
        # flag_band = 'satellite_WQSF'
        # flag_ac_value = 8

        return

    if args.mode == 'CHECK_PROTOCOLS' and args.input_path:
        if os.path.isfile(args.input_path):
            do_check_protocols(args.input_path)

    if args.input_path and args.config_file and (mode == 'UPDATE_SAT_WL' or mode == 'UPDATE_INSITU_WL'):
        fwl = args.config_file
        input_file = args.input_path
        if not os.path.exists(fwl):
            print(f'[ERROR] Configuration file {fwl} with wl list does not exist')
            return
        wl_list = get_wl_list_from_file(fwl)
        if wl_list is None:
            print(f'[ERROR] Error in the wavelenght list defined in file: {fwl}')
            return
        if not os.path.exists(input_file):
            print(f'[ERROR] Input MDB file {input_file} does not exist')
            return
        reader = MDB_READER(input_file, True)
        if not reader.mfile.VALID:
            print(f'[ERROR] Input file {input_file} is not a valid MDB file')
        if args.output:
            output_file = args.output
        else:
            dir_name = os.path.dirname(input_file)
            fname = os.path.basename(input_file)[:-3]
            output_file = os.path.join(dir_name, f'{fname}_OUTPUT.nc')
        if os.path.exists(output_file):
            print(f'[ERROR] Output file {output_file} already exist')
            return
        if args.verbose:
            print(f'[INFO] Input file: {input_file}')
            print(f'[INFO] Output file: {output_file}')
            print(f'[INFO] New wavelenght list: {wl_list}')
        if mode == 'UPDATE_SAT_WL':
            creating_copy_with_new_bands(reader, output_file, wl_list)
        if mode == 'UPDATE_INSITU_WL':
            creating_copy_with_new_insitu(reader, output_file, wl_list)

    if args.input_path and mode == 'CHECK_WL':
        print('Running option: Checking wavelength')
        input_path = args.input_path
        if not os.path.exists(input_path):
            print(f'[ERROR] Input path {input_path} does not exist')
            return
        reader = MDB_READER(input_path, True)
        if not reader.mfile.check_wl():
            reader.mfile.correct_wl()
        else:
            print(f'[INFO] Wavelengths in variables mu_wavelength and satellite_bands are coherent')
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

    ##ADDING FLAG BAND TO SINGLE MDB FILE
    if args.input_path and args.config_file and mode == 'ADDFLAGBAND':
        print('[INFO] Running option: add flag band')
        if not os.path.exists(args.input_path):
            print(f'[ERROR] Input file: {args.input_path} does not exist.')
            return
        if not os.path.exists(args.config_file):
            print(f'[ERROR] Flags configuration file: {args.config_file} does not exist.')
            return
        from FlagBuilder import FlagBuilder
        fbuilder = FlagBuilder(args.input_path, args.config_file, None)
        if not fbuilder.VALID:
            return
        for flag_here in fbuilder.flag_list:
            b = fbuilder.omanager.get_value_param(flag_here,'apply',True,'boolean')
            if b is not None and not b:
                continue
            print(f'[INFO] Working with flag: {flag_here}')
            fbuilder.create_flag_array(flag_here, True)

        return

    if args.input_path and mode == 'CHECK_SAT_TIME':
        print(f'[INFO] Running option: CHECK_SAT_TIME')
        input_path = args.input_path
        if not os.path.exists(input_path):
            print(f'[ERROR] Input path {input_path} does not exist')
            return
        output_path = os.path.join(os.path.dirname(input_path), 'TempCS.nc')
        creating_copy_correcting_sat_time(input_path, output_path)
        if os.path.exists(output_path):
            os.rename(output_path, input_path)

    ##CREATING MDB READER FILES WORKING WITH DEFAULT OPTIONS FROM A INPUT MDB FILE
    if args.input_path and mode == 'GENERATEMU':
        allow_repeated = False
        if args.allow_repeated:
            allow_repeated = True

        reduce_mdbr = False
        if args.reduce_mdbr:
            reduce_mdbr = True

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
        if not reader.mfile.check_repeated() and not allow_repeated:
            print('[ERROR] To remove repeated satellite ids you could run:')
            print(f'python MDB_readerV2.py -m REMOVEREP -i {input_path} -v')
            return
        wllist = None
        if args.config_file and os.path.exists(args.config_file):
            if args.verbose:
                print(f'[INFO] Using file: {args.config_file} to set quality control options...')
            import configparser
            from QC_OPTIONS import QC_OPTIONS
            options = configparser.ConfigParser()
            options.read(args.config_file)
            qco = QC_OPTIONS(options)
            wllist = qco.get_wllist()
            check_sat_bands = reader.mfile.check_bands(wllist, 5)
            if check_sat_bands == -1:
                return
            copy_with_wllist = False
            if check_sat_bands == 0:
                copy_with_wllist = qco.get_create_copy_with_band_list()

            reader.mfile.set_wl_ref(wllist)
            reader.mfile.qc_sat.ncdataset = reader.mfile.nc
            reader.mfile.qc_sat = qco.get_qcsat(reader.mfile.qc_sat, reader.mfile.nc)
            reader.mfile.qc_insitu.ncdataset = reader.mfile.nc
            reader.mfile.qc_insitu = qco.get_qc_insitu(reader.mfile.qc_insitu)
            if reader.mfile.qc_insitu is None:
                return
            reader.mfile.PI_DIVIDED = reader.mfile.qc_insitu.pi_divided
            if not reader.mfile.qc_insitu.apply_nir_correction:
                reader.mfile.qc_insitu.insitu_rrs = reader.mfile.variables['insitu_Rrs_nosc']
                if args.verbose:
                    print('[INFO] NIR Correction is not applied. Using insitu_Rrs_nosc as input variable')
        else:
            reader.set_defaults_olci_wfr_hypstar()

        reader.mfile.allow_repeated = allow_repeated
        reader.create_mdb_results_file(output_path, reduce_mdbr)
        if copy_with_wllist and wllist is not None:
            if args.verbose:
                print('[INFO] Creating copy with limited number of bands....')
            reader.mfile.close()
            reader = MDB_READER(output_path, True)
            file_temp = os.path.join(output_folder, 'TempMDB.nc')
            creating_copy_with_new_bands(reader, file_temp, wllist)
            os.rename(file_temp, output_path)

        return

    if args.input_path and mode == 'GENERATEMU_S':
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
            print(f'[INFO] No spectral variables are required for match-ups for single variables. No problem!')
        if not args.config_file:
            print(f'[ERROR] Config. file with Quality Control options (-c option) is required')
            return
        if not os.path.exists(args.config_file):
            print(f'[ERROR] Config. file {args.config_file} with Quality Control options does not exist')
            return
        if args.config_file and os.path.exists(args.config_file):
            if args.verbose:
                print(f'[INFO] Using file: {args.config_file} to set quality control options...')
            import configparser
            from QC_OPTIONS import QC_OPTIONS

            options = configparser.ConfigParser()
            options.read(args.config_file)
            qco = QC_OPTIONS(options)
            reader.mfile.qc_single = qco.get_qc_single(reader.mfile.file_path)
            reader.create_mdb_results_single(output_path)
            reader.mfile.qc_single.close_dataset()
            reader.mfile.close()

    ##CONCATENATING MDB READER FILES INCLUDING FLAGS BANDS FOR SATELLITE+PLATFORM, SENSOR, AC, SITE
    if args.input_path and mode == 'CONCATENATE':
        import pandas as pd
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

        # ats_in = [['satellite', 'platform'], 'sensor', 'satellite_aco_processor', 'insitu_site_name']
        ats_in = [['satellite', 'platform'], 'sensor', 'satellite_aco_processor', 'site']
        flag_bands = ['flag_satellite', 'flag_sensor', 'flag_ac', 'flag_site']
        flag_lists = get_flag_lists(input_path, ats_in, flag_bands)
        all_bands = get_band_list(input_path)

        idfile = 1
        satellite_id_ref = 0
        list_files = []
        csv_list_df = []
        for name in os.listdir(input_path):
            if name.endswith('.csv') and name.startswith('CSVr'):
                df = pd.read_csv(os.path.join(input_path, name), sep=';')
                csv_list_df.append(df)
                continue
            if not name.endswith('.nc'):
                continue
            if args.verbose:
                print(f'[INFO] Preparing file: {name} for concatenation...')
            file_in = os.path.join(input_path, name)
            reader = MDB_READER(file_in, True)
            bands_excluded = {}
            for band in all_bands:
                if band not in reader.mfile.variables:
                    # bands_excluded.append(band)
                    bands_excluded[band] = all_bands[band]
                    if args.verbose:
                        print(f'[INFO]-> Added not available band: {band}')
            file_out = os.path.join(input_path, f'Temp_{idfile}.nc')
            creating_copy_with_flag_bands(reader, file_out, flag_lists, satellite_id_ref, bands_excluded)
            satellite_id_ref = satellite_id_ref + reader.mfile.n_mu_total
            list_files.append(file_out)
            idfile = idfile + 1
        concatenate_nc_impl(list_files, input_path, ncout_file)
        if len(csv_list_df) > 0:
            dfend = pd.concat(csv_list_df, ignore_index=True)
            namef = os.path.basename(ncout_file)
            namef = namef.replace('.nc', '.csv')
            namef = namef.replace('MDB', 'CSV')
            fend = os.path.join(os.path.dirname(ncout_file), namef)
            dfend.to_csv(fend, sep=';')

    ##SETTIN VARIABLES INDICATING COMMON MU
    if args.input_path and mode.startswith('COMMONMU'):
        input_path = args.input_path
        if not os.path.exists(input_path):
            print(f'[ERROR] {input_path} does not exist')
            return
        if args.output:
            output_path = args.output
        else:
            output_path = f'{input_path[:-3]}_COMMONMU.nc'
        if args.verbose:
            print(f'[INFO] Input path: {input_path}')
            print(f'[INFO] Output path: {output_path}')
        reader = MDB_READER(input_path, True)
        if not reader.mfile.VALID:
            print(f'[INFO] Input path: {input_path} is not a valid MDB file')
            return
        if mode == 'COMMONMU':
            mu_valid_common = reader.mfile.obtain_info_common_mu()
        elif mode == 'COMMONMU_NOSAT':
            mu_valid_common = reader.mfile.obtain_info_common_mu_nosat()
        elif mode == 'COMMONMU_INS':
            mu_valid_common = reader.mfile.obtain_info_common_mu_ins()

        if mu_valid_common is None:
            print(f'[INFO] No valid common match-ups were defined in file: {input_path}')
            return
        if creating_copy_with_mu_valid_common(reader, output_path, mu_valid_common):
            print(f'[INFO] Output file: {output_path} completed')
        else:
            print(f'[ERROR] Error creating output file: {output_path}')

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

        ##WITH MDBPLOTV3
        from MDBPlotV3 import MDBPlot
        mplot = MDBPlot(input_path)
        mplot.plot_from_options_file(config_file)

        # mplot.output_path = output_path
        #

        ##WITH MDBPlotV2
        # from MDBPlotV2 import MDBPlot
        # import configparser
        # mplot = MDBPlot(input_path)
        # options = configparser.ConfigParser()
        # options.read(config_file)
        # # print(mplot.VALID)
        # mplot.set_global_options(options)
        # mplot.plot_from_options(options)

        # path_img = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/PLOTS'
        # from PlotMultiple import PlotMultiple
        # file_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S3OLCI/PLOTS/'

        # path_img = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/PLOTS'
        # from PlotMultiple import PlotMultiple
        # file_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/S2MSI/PLOTS/'

        # path_img = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/CMEMSMULTI/VEIT/PLOTS'
        # from PlotMultiple import PlotMultiple
        # file_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/WP7_FINAL_ANALYSIS/MDBs/CMEMSMULTI/VEIT/PLOTS'

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
        # file00 = os.path.join(file_base,'PLOT_STAT_GLOBAL_RMSD.jpg')
        # file10 = os.path.join(file_base,'PLOT_STAT_GLOBAL_DETER(r2).jpg')
        # file20 = os.path.join(file_base, 'PLOT_STAT_GLOBAL_BIAS.jpg')
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

    ##PLOTTING CSV
    if args.input_path and args.config_file and mode == 'PLOT_CSV':
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
            path_csv = os.path.dirname(input_path)
            output_path = os.path.join(path_csv, 'PLOTS')
        if not os.path.exists(output_path):
            try:
                os.mkdir(output_path)
            except:
                pass
        if not os.path.isdir(output_path):
            print(f'[ERROR] Ouput path: {output_path} does not exist or is not a directory')
            return

        ##temporal,stats plots from multiple csv files
        if os.path.isdir(input_path):
            # do_temporal_wl_stats_from_csv('4owt',input_path)
            # do_temporal_wl_stats_from_csv('6owt', input_path)
            do_temporal_spectra_stats_from_csv(input_path)
            # do_temporal_owt_comparison(input_path)
            return
        if os.path.basename(config_file) == 'skip.ini':
            # do_temporal_all_owt_stats('CHL',input_path)
            do_temporal_all_owt_stats('TSM', input_path)
            return

        from CSVPlot import CSVPLOT
        cplot = CSVPLOT(input_path)
        if not cplot.VALID:
            print(f'[ERROR] Input path {input_path} is not valid')
            return
        import configparser
        options = configparser.ConfigParser()
        options.read(config_file)
        cplot.plot_from_options(options)


if __name__ == '__main__':
    main()
