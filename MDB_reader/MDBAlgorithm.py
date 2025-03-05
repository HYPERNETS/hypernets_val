import argparse
import warnings
import os
import numpy as np
from datetime import datetime as dt

import pandas as pd
import pytz
from netCDF4 import Dataset

from MDBFile import MDBFile
from MDB_builder.INSITU_base import INSITUBASE
import __init__

warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', RuntimeWarning)

parser = argparse.ArgumentParser(description="Algorithms implementations from MDB files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument("-m", "--mode", help="Mode", choices=["CONFIGFILE", "CYANOFLAG","BALMLP202411"])
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-i', "--input_path", help="Input MDB path", required=True)
parser.add_argument('-o', "--output", help="Path to output")
args = parser.parse_args()


def do_test():
    ####TEMPORAL
    # file_nc = '/mnt/c/DATA_LUIS/OCTAC_WORK/MATCH-UPS_ANALYSIS_2024/BAL/MDBs/MDBr__MULTI_CCI_1KM_OC-CCI-BALCHL2021_19970316T000000_20240415T000000.nc'
    file_nc = '/mnt/c/DATA_LUIS/OCTAC_WORK/MATCH-UPS_ANALYSIS_2024/BAL/MDBs/MDBr__S3_OLCI_300M_CMEMS-OLCI_20160404T000000_20240415T000000.nc'
    dataset = Dataset(file_nc)
    insitu_id = dataset.variables['mu_insitu_CHLA_id'][:]
    insitu_time = dataset.variables['insitu_time'][:]
    flag_cyano_array = dataset.variables['satellite_flag_cyano'][:]
    file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/MATCH-UPS_ANALYSIS_2024/BAL/MDBs/FlagCyanoOlci.csv'
    fout = open(file_out, 'w')
    fout.write('DATE;HOUR;FLAG_CYANO;N_VALID;N_CYANO')
    for idx in range(flag_cyano_array.shape[0]):
        if idx % 100 == 0: print(idx)
        insitu_id_here = insitu_id[idx]
        if np.ma.is_masked(insitu_id_here):
            continue

        insitu_time_here = float(insitu_time[idx][insitu_id_here])
        insitu_time_here_dt = dt.utcfromtimestamp(insitu_time_here)
        date = insitu_time_here_dt.strftime('%Y-%m-%d')
        time = insitu_time_here_dt.strftime('%H:%M:%S')
        # rrs_555 = satellite_Rrs_555[idx, 12, 12]
        # rrs_670 = satellite_Rrs_670[idx, 12, 12]
        flag_cyano_here = flag_cyano_array[idx, 12, 12]
        flag_cyano_values = flag_cyano_array[idx, 11:14, 11:14]
        flag_cyano_values_good = flag_cyano_values[~flag_cyano_values.mask]
        n_all = len(flag_cyano_values_good)
        flag_cyano_centre = flag_cyano_values_good[flag_cyano_values_good == flag_cyano_here]
        n_cyano = len(flag_cyano_centre)
        # line = f'{date};{time};{rrs_555};{rrs_670};{flag_cyano_here}'
        line = f'{date};{time};{flag_cyano_here};{n_all};{n_cyano}'
        fout.write('\n')
        fout.write(line)
    fout.close()
    dataset.close()

    ####

    return True


def create_mdb_from_csv():
    file_csv = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/MATCH-UPS_ANALYSIS_2024/CSV_MATCH-UPS/MULTI/Baltic_CHLA_Valid_AllSources_1997-2023_FINAL_TIMEFILTERED_complete.csv'
    dir_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/MATCH-UPS_ANALYSIS_2024/MDBs'
    file_out = os.path.join(dir_out, 'MDBr__MULTI_CCI_1KM_OC-CCI-BALCHL202411_19970909T000000_20231219T000000.nc')
    nc_out = Dataset(file_out, 'w')
    nc_out.createDimension('mu_id', size=None)
    nc_out.createDimension('satellite_id', size=None)
    df = pd.read_csv(file_csv, sep=';')
    col_names =df.columns.tolist()

    mu_variables = ['LATITUDE', 'LONGITUDE', 'INSITU_CHL','FLAG_CYANO_PORC_MODE','DISTANCE']

    flag_variables_num = ['satellite_CHL_NVALID','FLAG_CDF','BLOOM','SUB_SURFACE','SURFACE']
    flag_variables_str = ['SOURCE_ORIG','SOURCE','FlagPrec','FlagBrando','FlagOldNew','FLAG_DISTANCE']

    for col in col_names:
        if col.startswith('satellite_RRS'):
            mu_variables.append(col)
        if col.startswith('satellite_CHL') and not col.endswith('NVALID'):
            mu_variables.append(col)
        if col.startswith('satellite_CDF'):
            mu_variables.append(col)
        if col.startswith('satellite_WEIGHT'):
            mu_variables.append(col)
        if col.startswith('FLAG_CYANO') and col!='FLAG_CYANO_PORC_MODE':
            flag_variables_num.append(col)
        if col.startswith('use'):
            flag_variables_num.append(col)
        if col.startswith('FLAG_ERROR'):
            flag_variables_num.append(col)
        if col.startswith('FLAG_NOERROR'):
            flag_variables_num.append(col)


        # if col.startswith('chl_') or col.startswith('cdf_'):
        #     mu_variables.append(col)
        # if col.endswith('_cdf') or col.startswith('use_cdf'):
        #     flag_variables_num.append(col)
        # if col.endswith('_NVALID'):
        #     flag_variables_num.append(col)




    datetime_array = np.array([dt.strptime(x,'%Y-%m-%dT%H:%M:%S').replace(tzinfo=pytz.UTC).timestamp() for x in df['DATETIME']])
    var = nc_out.createVariable('time','f4',('satellite_id',),zlib=True,complevel=6,fill_value=-999.0)
    var[:] = datetime_array

    for var_name in mu_variables:
        print('Creating mu variable:',var_name)
        var = nc_out.createVariable(var_name.upper(),'f4',('mu_id',),zlib=True,complevel=6,fill_value=-999.0)
        var[:] = df[var_name]

    for var_name in flag_variables_num:
        print('Creating flag num',var_name)
        var_array = np.array(df[var_name])
        if var_name=='FLAG_CYANO' or var_name=='FLAG_CYANO_CENTRAL' or var_name=='FLAG_CYANO_MODE':
            flag_values = [0,1,2,3]
            flag_meanings = 'NO_BLOOM SUB_SURFACE_BLOOM SURFACE_BLOOM BOTH_BLOOMS'
        elif var_name == 'FLAG_CDF':
            flag_values = [1,2,4,8,6,10,12,14]
            flag_meanings = 'NO_CDF CDF_MLP_3B CDF_MLP_4B CDF_MLP_5B CDF_MLP_3B+4B CDF_MLP_3B+5B CDF_MLP_4B+5B CDF_MLP_3B+4B+5B'
        elif var_name == 'BLOOM' or var_name == 'SUB_SURFACE' or var_name == 'SURFACE':
            flag_values = [1,2]
            flag_meanings = 'NO_BLOOM BLOOM'
        elif var_name == 'FLAG_CYANO_HOMOGENEITY':
            flag_values = [0,1]
            flag_meanings = ['CYANO_MIX','CYANO_HOMEGENEOUS']
        else:
            values_unique = np.unique(var_array).tolist()
            flag_values = values_unique
            flag_meanings = ' '.join([str(x) for x in flag_values])
        var = nc_out.createVariable(var_name.upper(), 'i4', ('satellite_id',), zlib=True, complevel=6, fill_value=-999)
        var[:] = var_array[:]
        var.flag_meanings = flag_meanings
        var.flag_values = flag_values




    for var_name in flag_variables_str:
        print('Flag str:', var_name)
        var_array_obj = np.array(df[var_name])
        obj_list = np.unique(var_array_obj).tolist()

        var_array = np.zeros((var_array_obj.shape))
        var_array[:] = -999
        flag_values = []
        for idx,obj in enumerate(obj_list):
            flag_value = 2**idx
            flag_values.append(flag_value)
            var_array[var_array_obj==obj] = flag_value
        flag_meanings = ' '.join(obj_list)
        var = nc_out.createVariable(var_name.upper(), 'i4', ('satellite_id',), zlib=True, complevel=6, fill_value=-999)
        var[:] = var_array[:]
        var.flag_meanings = flag_meanings
        var.flag_values = flag_values



    nc_out.close()


    return True


def main():
    # if create_mdb_from_csv():
    #     return
    # if do_test():
    #     return
    print('Started MDBAlgorithm')
    input_path = args.input_path
    if not os.path.exists(input_path):
        print(f'[ERROR] Input path {input_path} does not exist')
        return
    try:
        dset = Dataset(input_path)
        dset.close()
    except:
        print(f'[ERROR] Input path {input_path} is not a valid MDB (NetCDF File)')
        return
    output_path = None
    if args.output:
        output_path = args.output
        if not os.path.isdir(os.path.dirname(output_path)):
            try:
                os.mkdir(os.path.dirname(output_path))
            except:
                print(
                    f'[ERROR] Ouput path {os.path.basename(output_path)} is not valid as {os.path.dirname(output_path)} is not a valid directory')
                return
        if output_path.endswith('.nc'):
            print(f'[ERROR] Output path {output_path} should be a NC file (.nc)')
            return

    if args.mode == 'CYANOFLAG':
        create_cyano_flag(input_path, output_path)

    if args.mode == 'BALMLP202411':
        try:
            from  baltic_202411 import BALTIC_202411_PROCESSOR
        except:
            print(f'[ERROR] baltic_202411 code is not available')
        bprocessor = BALTIC_202411_PROCESSOR(None, False)
        bprocessor.run_from_mdb_file(args.input_path)
        return

def create_cyano_flag(input_path, output_path):
    if output_path is None:
        dataset_w = Dataset(input_path, 'a')
    else:
        ibase = INSITUBASE(None)
        dataset_w = ibase.copy_nc(input_path, output_path)

    ##satellite variable
    if not 'satellite_Rrs' in dataset_w.variables or not 'satellite_bands' in dataset_w.variables:
        print(f'[ERROR] satellite_Rrs and satellite_bands are required to compute CYANOFLAG')
        dataset_w.close()
        if output_path is not None:
            os.remove(output_path)
        return
    th_555_sub = 4.25e-3
    th_670_sur = 1.22e-3
    satellite_bands = dataset_w.variables['satellite_bands'][:]
    index_555 = np.ma.argmin(np.ma.abs(satellite_bands - 555.0))
    wl_555 = satellite_bands[index_555]
    diff_555 = abs(wl_555 - 555.0)
    if diff_555 > 5:
        print(f'[ERROR] Band at 555 nm is not avaialable (nearest band is {wl_555})')
        return
    else:
        print(f'[INFO] Wavelength for 555 nm (sub-surface blooms): {wl_555}')

    index_670 = np.argmin(np.abs(satellite_bands - 670.0))
    wl_670 = satellite_bands[index_670]
    diff_670 = abs(wl_670 - 670.0)
    if diff_670 > 5:
        print(f'[ERROR] Band at 670 nm is not avaialable (nearest band is {wl_670})')
        return
    else:
        if wl_670 > 670:  ##665 better than 673.75
            satellite_bands_l = satellite_bands[satellite_bands < 670]
            index_670_l = np.argmin(np.abs(satellite_bands_l - 670.0))
            wl_670_l = satellite_bands[index_670_l]
            diff_670_l = abs(wl_670_l - 670.0)
            if diff_670_l <= 5:
                index_670 = index_670_l
                wl_670 = wl_670_l
                diff_670 = diff_670_l
        print(f'[INFO] Wavelength for 670 nm (sub-surface blooms): {wl_670}')

    satellite_Rrs = dataset_w.variables['satellite_Rrs']

    if diff_555 == 0.0:
        satellite_Rrs_555 = np.ma.squeeze(satellite_Rrs[:, index_555, :, :])
    else:
        ##apply band shifting
        print(f'[INFO] Applying band shifting from {wl_555} nm to 555.0 nm')
        from BSC_QAA import bsc_qaa_EUMETSAT as bsc
        satellite_Rrs_555 = np.ma.zeros((satellite_Rrs.shape[0], satellite_Rrs.shape[2], satellite_Rrs.shape[3]))
        for index_mu in range(satellite_Rrs.shape[0]):
            rrs_in = np.ma.squeeze(satellite_Rrs[index_mu, :, :, :])
            ndata_valid = np.ma.count(rrs_in)
            if ndata_valid == 0:
                satellite_Rrs_555[index_mu, :, :] = np.ma.masked
            else:
                satellite_Rrs_555[index_mu, :, :] = bsc.bsc_qaa(rrs_in, satellite_bands, np.ma.array([555.0]))

    if diff_670 == 0.0:
        satellite_Rrs_670 = np.ma.squeeze(satellite_Rrs[:, index_670, :, :])
    else:
        ##apply band shifting
        print(f'[INFO] Applying band shifting from {wl_670} nm to 670.0 nm')
        from BSC_QAA import bsc_qaa_EUMETSAT as bsc
        satellite_Rrs_670 = np.ma.zeros((satellite_Rrs.shape[0], satellite_Rrs.shape[2], satellite_Rrs.shape[3]))
        for index_mu in range(satellite_Rrs.shape[0]):
            rrs_in = np.ma.squeeze(satellite_Rrs[index_mu, :, :, :])
            ndata_valid = np.ma.count(rrs_in)
            if ndata_valid == 0:
                satellite_Rrs_670[index_mu, :, :] = np.ma.masked
            else:
                satellite_Rrs_670[index_mu, :, :] = bsc.bsc_qaa(rrs_in, satellite_bands, np.ma.array([670.0]))

    if 'satellite_flag_cyano' not in dataset_w.variables:
        print(f'[INFO] Creating variable satellite_flag_cyano')
        satellite_cyano_array = np.zeros(satellite_Rrs_555.shape)
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 >= th_555_sub, satellite_Rrs_670 >= th_670_sur)] = 3
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 >= th_555_sub, satellite_Rrs_670 < th_670_sur)] = 1
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 < th_555_sub, satellite_Rrs_670 >= th_670_sur)] = 2
        satellite_cyano_array[np.logical_and(satellite_Rrs_555.mask, satellite_Rrs_670.mask)] = -999
        satellite_cyano_var = dataset_w.createVariable('satellite_flag_cyano', 'i2',
                                                       ('satellite_id', 'rows', 'columns'),
                                                       zlib=True, complevel=6, fill_value=-999.0)
        satellite_cyano_var.descripton = 'Satellite Cyano Flag'
        satellite_cyano_var.flag_masks = [0, 1, 2, 3]
        satellite_cyano_var.flag_meanings = "NO_BLOOM SUB-SURFACE_BLOOM SURFACE_BLOOM BOTH_BLOOMS"
        satellite_cyano_var[:] = satellite_cyano_array[:]
    else:
        print(f'[WARNING] Variable satellite_flag_cyano already exits. Skipping...')
        satellite_cyano_array = np.zeros(satellite_Rrs_555.shape)
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 >= th_555_sub, satellite_Rrs_670 >= th_670_sur)] = 3
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 >= th_555_sub, satellite_Rrs_670 < th_670_sur)] = 1
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 < th_555_sub, satellite_Rrs_670 >= th_670_sur)] = 2
        satellite_cyano_array[np.logical_and(satellite_Rrs_555.mask, satellite_Rrs_670.mask)] = -999
        # satellite_cyano_var = dataset_w.variables['satellite_flag_cyano']
        # satellite_cyano_var.flag_masks = [0, 1, 2, 3]
        # satellite_cyano_var.flag_meanings = "NO_BLOOM SUB-SURFACE_BLOOM SURFACE_BLOOM BOTH_BLOOMS"
        # satellite_cyano_var[:] = satellite_cyano_array[:]

    if 'flag_cyano' not in dataset_w.variables:
        print(f'[INFO] Creating variable flag_cyano')
        flag_cyano_array = np.ma.squeeze(satellite_cyano_array[:, 12, 12])
        flag_cyano_var = dataset_w.createVariable('flag_cyano', 'i2', ('satellite_id',), zlib=True, complevel=6,
                                                  fill_value=-999)
        flag_cyano_var.flag_values = [0, 1, 2, 3]
        flag_cyano_var.flag_meanings = "NO_BLOOM SUB-SURFACE_BLOOM SURFACE_BLOOM BOTH_BLOOMS"
        flag_cyano_var[:] = flag_cyano_array
    else:
        print(f'[WARNING] Variable flag_cyano already exits. Skipping...')
        # flag_cyano_array = np.ma.squeeze(satellite_cyano_array[:, 12, 12])

        # flag_cyano_var = dataset_w.variables['flag_cyano']
        # flag_cyano_var.flag_values = [0, 1, 2, 3]
        # flag_cyano_var.flag_meanings = "NO_BLOOM SUB-SURFACE_BLOOM SURFACE_BLOOM BOTH_BLOOMS"
        # flag_cyano_var[:] = flag_cyano_array

    # insitu_time = dataset_w.variables['insitu_time'][:]
    # insitu_id = dataset_w.variables['mu_insitu_CHLA_id'][:]

    dataset_w.close()

    if output_path is None:
        output_path = input_path
    print(f'[INFO] Completed. Output file {output_path}')


if __name__ == '__main__':
    main()
