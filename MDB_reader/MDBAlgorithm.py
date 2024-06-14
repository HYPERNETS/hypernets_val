import argparse
import warnings
import os
import numpy as np
from datetime import datetime as dt
from netCDF4 import Dataset

from MDBFile import MDBFile
from MDB_builder.INSITU_base import INSITUBASE

warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', RuntimeWarning)

parser = argparse.ArgumentParser(description="Algorithms implementations from MDB files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument("-m", "--mode", help="Mode", choices=["CONFIGFILE", "CYANOFLAG"])
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-i', "--input_path", help="Input MDB path", required=True)
parser.add_argument('-o', "--output", help="Path to output")
args = parser.parse_args()


def main():
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
        print(f'[INFO] Wavelength for 670 nm (sub-surface blooms): {wl_670}')

    satellite_Rrs = dataset_w.variables['satellite_Rrs']

    if diff_555==0.0:
        satellite_Rrs_555 = np.ma.squeeze(satellite_Rrs[:, index_555, :, :])
    else:
        ##apply band shifting
        print(f'[INFO] Applying band shifting from {wl_555} nm to 555.0 nm')
        from BSC_QAA import bsc_qaa_EUMETSAT as bsc
        satellite_Rrs_555 = np.ma.zeros((satellite_Rrs.shape[0],satellite_Rrs.shape[2],satellite_Rrs.shape[3]))
        for index_mu in range(satellite_Rrs.shape[0]):
            rrs_in = np.ma.squeeze(satellite_Rrs[index_mu,:,:,:])
            ndata_valid = np.ma.count(rrs_in)
            if ndata_valid==0:
                satellite_Rrs_555[index_mu,:,:] = np.ma.masked
            else:
                satellite_Rrs_555[index_mu,:,:]  = bsc.bsc_qaa(rrs_in, satellite_bands,np.ma.array([555.0]))
            if index_mu == 689:
                print('555',satellite_Rrs_555[index_mu, 12, 12])

    if diff_670==0.0:
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
            if index_mu==689:
                print('670',satellite_Rrs_670[index_mu,12,12])

    if 'satellite_flag_cyano' not in dataset_w.variables:
        print(f'[INFO] Creating variable satellite_flag_cyano')
        satellite_cyano_array = np.zeros(satellite_Rrs_555.shape)
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 >= th_555_sub, satellite_Rrs_670 >= th_670_sur)] = 3
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 >= th_555_sub, satellite_Rrs_670 < th_670_sur)] = 1
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 < th_555_sub, satellite_Rrs_670 >= th_670_sur)] = 2
        satellite_cyano_array[np.logical_and(satellite_Rrs_555.mask,satellite_Rrs_670.mask)] = -999
        satellite_cyano_var = dataset_w.createVariable('satellite_flag_cyano', 'i2', ('satellite_id', 'rows', 'columns'),
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
        flag_cyano_array = np.ma.squeeze(satellite_cyano_array[:,12,12])
        flag_cyano_var = dataset_w.createVariable('flag_cyano','i2',('satellite_id',),zlib=True,complevel=6,fill_value=-999)
        flag_cyano_var.flag_values = [0, 1, 2, 3]
        flag_cyano_var.flag_meanings = "NO_BLOOM SUB-SURFACE_BLOOM SURFACE_BLOOM BOTH_BLOOMS"
        flag_cyano_var[:] = flag_cyano_array
    else:
        print(f'[WARNING] Variable flag_cyano already exits. Skipping...')
        flag_cyano_array = np.ma.squeeze(satellite_cyano_array[:, 12, 12])

        # flag_cyano_var = dataset_w.variables['flag_cyano']
        # flag_cyano_var.flag_values = [0, 1, 2, 3]
        # flag_cyano_var.flag_meanings = "NO_BLOOM SUB-SURFACE_BLOOM SURFACE_BLOOM BOTH_BLOOMS"
        # flag_cyano_var[:] = flag_cyano_array


    insitu_time = dataset_w.variables['insitu_time'][:]
    insitu_id = dataset_w.variables['mu_insitu_CHLA_id'][:]

    dataset_w.close()

    ####TEMPORAL
    file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/MATCH-UPS_ANALYSIS_2024/BAL/MDBs/FlagCyano.csv'
    fout = open(file_out,'w')
    fout.write('DATE;HOUR;RRS555;RRR670;FLAG_CYANO')
    for idx in range(10294):
        insitu_id_here = insitu_id[idx]
        insitu_time_here = float(insitu_time[idx][insitu_id_here])
        insitu_time_here_dt = dt.utcfromtimestamp(insitu_time_here)
        date = insitu_time_here_dt.strftime('%Y-%m-%d')
        time = insitu_time_here_dt.strftime('%H:%M:%S')
        rrs_555 = satellite_Rrs_555[idx,12,12]
        rrs_670 = satellite_Rrs_670[idx,12,12]
        flag_cyano_here = flag_cyano_array[idx]
        line = f'{date};{time};{rrs_555};{rrs_670};{flag_cyano_here}'
        fout.write('\n')
        fout.write(line)
    fout.close()

    ####


    if output_path is None:
        output_path = input_path
    print(f'[INFO] Completed. Output file {output_path}')

if __name__ == '__main__':
    main()
