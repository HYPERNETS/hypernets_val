import argparse
import os
from datetime import datetime as dt
from datetime import timedelta
import sys
import INSITU_comparison
from INSITU_comparison import INSITUCOMPARISON

code_home = os.path.dirname(os.path.dirname(INSITU_comparison.__file__))
sys.path.append(code_home)
code_aeronet = os.path.join(os.path.dirname(code_home), 'aeronet')
sys.path.append(code_aeronet)

parser = argparse.ArgumentParser(description="Creation of insitu nc files")

parser.add_argument('-c', "--config_file", help="Config File.")
# parser.add_argument('-o', "--output",help="Output ")
# parser.add_argument('-edir', "--sat_extract_dir",
#                     help="Input sat. extract dir. Optional for --listdates, required for single concatenation")
# parser.add_argument('-site', "--sitename", help="Site name. Only required with --listdates")
# parser.add_argument('-ld', "--listdates",
#                     help="Option to obtain a date list for a specific HYPERNETS site (-site option).",
#                     action="store_true")
# parser.add_argument('-sd', "--start_date", help="Start date. Optional with --listdates (YYYY-mm-dd)")
# parser.add_argument('-ed', "--end_date", help="End date. Optional with --listdates (YYYY-mm-dd)")
# parser.add_argument('-nd', "--nodelfiles", help="Do not delete temp files.", action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")

args = parser.parse_args()


def make_concatenation():
    path_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC'
    file_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/Comparison_20230323_20231031.nc'
    files_append = []
    start_date = dt(2023, 3, 23)
    end_date = dt(2023, 10, 31)
    date_here = start_date
    while date_here <= end_date:
        file_comparison = get_file_comparison_date(path_out, date_here, False)
        if file_comparison is not None:
            files_append.append(file_comparison)
        date_here = date_here + timedelta(hours=24)
    concatenate_nc_impl(files_append, path_out, file_out)


def create_multiple_comparison_files():
    print(f'[INFO] Started in situ comparison...')
    path_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC'
    start_date = dt(2023, 3, 23)
    end_date = dt(2023, 10, 31)
    date_here = start_date
    while date_here <= end_date:
        date_here_str = date_here.strftime('%Y-%m-%d')
        print(f'[INFO] Working with date: {date_here_str}')
        create_comparison_file_date(date_here)
        # add_mu_to_file(date_here)
        date_here = date_here + timedelta(hours=24)


def create_comparison_file_date(date_here):
    print(f'[INFO] Started in situ comparison for a specific date...')
    path_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC'
    # date_here = dt(2023, 9, 26)
    aeronet_file = '/mnt/c/DATA_LUIS/AERONET_OC/AERONET_NC/20020101_20231111_AAOT.LWN_lev20_15.nc'
    path_hypstar = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT'

    file_comparison = get_file_comparison_date(path_out, date_here, True)

    path_hypstar_date = get_folder_date(path_hypstar, date_here)

    from base.anet_nc_reader import AERONETReader
    areader = AERONETReader(aeronet_file)
    ic = INSITUCOMPARISON(file_comparison)
    ic.start_file_w()
    baeronet = ic.create_aeronet_variables(areader, date_here)
    if not baeronet:
        print(f'[ERROR] Comparison file for date: {date_here} could not be done. ')
        os.remove(file_comparison)
        return
    ic.create_hypstar_variables(path_hypstar_date, date_here)
    ic.create_hypstar_to_aeronet_variables()
    ic.close_file_w()


# if create is True, create the output folder
def get_folder_date(path_out, date_here):
    yyyy = date_here.strftime('%Y')
    mm = date_here.strftime('%m')
    dd = date_here.strftime('%d')
    folder_date = os.path.join(path_out, yyyy, mm, dd)
    if os.path.isdir(folder_date):
        return folder_date
    else:
        return None

def add_mu_to_file_date(date_here):
    path_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC'
    # date_here = dt(2023, 9, 26)
    file_comparison = get_file_comparison_date(path_out, date_here, False)
    if file_comparison is not None:
        add_mu_to_file(file_comparison)

def add_mu_to_file(file_comparison):

    if file_comparison is None:
        return
    print(f'[INFO] Working with file: {file_comparison}')
    ic = INSITUCOMPARISON(file_comparison)
    ic.add_match_ups_to_file(11)


def make_plots():
    # path_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC'
    # date_here = dt(2023, 9, 26)
    # file_comparison = get_file_comparison_date(path_out, date_here, False)
    file_comparison = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/Comparison_20230323_20231031.nc'
    print(file_comparison)

    ic = INSITUCOMPARISON(file_comparison)
    #ic.plot_scatterplot(None)
    # ic.plot_all_scatterplots_wl()
    ic.plot_all_spectra()
    # ic.plot_rho_scatterplot()


def get_file_comparison_date(path_out, date_here, to_create):
    yyyy = date_here.strftime('%Y')
    mm = date_here.strftime('%m')
    dd = date_here.strftime('%d')
    date_here_str = date_here.strftime('%Y%m%d')
    path_year = os.path.join(path_out, yyyy)
    path_month = os.path.join(path_year, mm)
    path_day = os.path.join(path_month, dd)
    if to_create:
        create_dir(path_year)
        create_dir(path_month)
        create_dir(path_day)
    file_comparison = os.path.join(path_day, f'COMPARISON_{date_here_str}.nc')
    if to_create:
        return file_comparison
    if os.path.exists(file_comparison):
        return file_comparison
    else:
        return None


def create_dir(path):
    if not os.path.isdir(path):
        os.mkdir(path)

def concatenate_nc_impl(list_files, path_out, ncout_file):
    if len(list_files) == 0:
        print(f'[WARNING] No MDB sat extract files were found. Please review')
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
            if icent == 74:
                print(list_files_here)
            ncout_file_tmp = os.path.join(path_out, f'Temp_{indextmp}.nc')
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

        # if not args.nodelfiles:
        #     [os.remove(f) for f in list_files]

        for f in list_files_tmp:
            fname = f.split('/')[-1]
            if fname.startswith('Temp_'):
                os.remove(f)


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

        if not args.nodelfiles:
            [os.remove(f) for f in list_files[:-1]]
    if args.verbose:
        print(f'[INFO] Concatenated MDB file created: {ncout_file}')

def check_angles():
    print('checking')
    from netCDF4 import Dataset
    path_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT'
    start_date = dt(2023, 3, 23)
    end_date = dt(2023, 10, 31)
    date_here = start_date
    file_out = os.path.join(path_out,'Azimuth_angles.csv')
    fout = open(file_out,'w')
    fout.write('Date;Time;ViewingAzimuthAngle;PointAzimuthAngle')
    while date_here <= end_date:
        date_here_str = date_here.strftime('%Y-%m-%d')

        print(f'[INFO] Working with date: {date_here_str}')
        folder = get_folder_date(path_out,date_here)
        for name in os.listdir(folder):
            if name.find('L2A')>0:
                file = os.path.join(folder,name)
                dataset = Dataset(file)
                vaa = float(dataset.variables['viewing_azimuth_angle'][0])
                paa = float(dataset.variables['pointing_azimuth_angle'][0])
                time_stamp = float(dataset.variables['acquisition_time'][0])
                time_here = dt.utcfromtimestamp(time_stamp)
                time_here_str = time_here.strftime('%H:%M')
                #time_here_str = date_here.strftime(dt.utcfromtimestamp(float(dataset.variables['acquisition_time'][0])))
                line  = f'{date_here_str};{time_here_str};{vaa};{paa}'
                print(line)
                fout.write('\n')
                fout.write(line)
                #print(date_here_str,dataset.variables['viewing_azimuth_angle'][0],dataset.variables['pointing_azimuth_angle'][0])
                dataset.close()
        date_here = date_here + timedelta(hours=24)
    fout.close()
def check_hypstar_qf():
    print('checking')
    #file = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT/2023/09/26/HYPERNETS_W_VEIT_L1C_ALL_20230926T0700_20240116T1733_270_v2.0.nc'

    file = '/mnt/c/DATA_LUIS/AERONET_OC/AERONET_NC/20020101_20231111_AAOT.LWN_lev20_15.nc'

    #file = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/2023/04/19/COMPARISON_20230419.nc'
    from netCDF4 import Dataset
    import numpy as np
    dset = Dataset(file)
    for var in dset.variables:
        print(var)
    # var= np.array(dset.variables['HYPSTAR_Lw'][0,4])
    # print(var)
    # for idx in range(60):
    #     v = var[0,idx]
    #     print(idx,dt.utcfromtimestamp(v))
    dset.close()
    # path_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT'
    # start_date = dt(2023, 3, 23)
    # end_date = dt(2023, 10, 31)
    # date_here = start_date
    # while date_here <= end_date:
    #     date_here_str = date_here.strftime('%Y-%m-%d')
    #     path_date = get_folder_date(path_out,date_here)


# %%
if __name__ == '__main__':
    #multiple_dates
    ##1. Create files
    #create_multiple_comparison_files()
    ##2. Concatenate
    #make_concatenation()
    ##3. Add mu
    #file_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/Comparison_20230323_20231031.nc'
    #add_mu_to_file(file_out)

    ## example for a specific date
    # date_here = dt(2023, 9, 26)
    # create_comparison_file_date(date_here)
    # add_mu_to_file_date(date_here)

    ##PLOTING
    #make_plots()

    ##CHECKING HYPSTAR QF (TEST)
    check_hypstar_qf()
    #check_angles()
