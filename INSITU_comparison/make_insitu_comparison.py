import argparse
import os
import shutil
from datetime import datetime as dt
from datetime import timedelta
import sys
import INSITU_comparison
from INSITU_comparison import INSITUCOMPARISON
from INSITU_plots import INSITU_plots

code_home = os.path.dirname(os.path.dirname(INSITU_comparison.__file__))
sys.path.append(code_home)
code_aeronet = os.path.join(os.path.dirname(code_home), 'aeronet')
sys.path.append(code_aeronet)

parser = argparse.ArgumentParser(description="Creation of insitu nc files")
parser.add_argument('-m', "--mode", help="Mode.", choices=["COMPARISON", "CONCAT", "GENERATEMU", "TEST"])
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


def make_concatenation_test():
    path_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC'
    file_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/Comparison_20230425_20230625.nc'
    files_append = []
    start_date = dt(2023, 4, 25)
    end_date = dt(2023, 6, 25)
    date_here = start_date
    while date_here <= end_date:
        file_comparison = get_file_comparison_date(path_out, date_here, False)
        if file_comparison is not None:
            files_append.append(file_comparison)
        date_here = date_here + timedelta(hours=24)


def make_concatenation(path_out, file_out, start_dates, end_dates):
    files_append = []

    for idate in range(len(start_dates)):
        files_append = get_files_concatenation(files_append, path_out, start_dates[idate], end_dates[idate])

    if len(files_append)>1:
        if args.verbose:
            print(f'[INFO] Concatenating: {len(files_append)}')
        concatenate_nc_impl(files_append, path_out, file_out)
    else:
        if args.verbose:
            print(f'[INFO] Less than 1 files were retrieved for concatenation. Operation cancelled')


def get_files_concatenation(files_append, path_out, start_date, end_date):
    if args.verbose:
        print(f'[INFO] Getting files from {start_date} to {end_date}')
    date_here = start_date
    nfiles = 0
    while date_here <= end_date:
        file_comparison = get_file_comparison_date(path_out, date_here, False)
        if file_comparison is not None:
            files_append.append(file_comparison)
            nfiles = nfiles + 1
        date_here = date_here + timedelta(hours=24)
    if args.verbose:
        print(f'[INFO] Files retrieved: {nfiles} Total: {len(files_append)}')
    return files_append


def create_multiple_comparison_files(start_date, end_date, aeronet_file, path_hypstar, path_out, overwrite):
    # path_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC'
    # start_date = dt(2023, 3, 24)
    # end_date = dt(2023, 10, 31)
    if args.verbose:
        print(f'[INFO] Started in situ comparison...')
        print(f'[INFO] > AERONET-NC file: {aeronet_file}')
        print(f'[INFO] > HYPSTAR input path: {path_hypstar}')
        print(f'[INFO] > Output path: {path_out}')
        print(f'[INFO] > Start date: {start_date.strftime("%Y-%m-%d")}')
        print(f'[INFO] > End date: {start_date.strftime("%Y-%m-%d")}')
        print(f'[INFO] > Overwrite: {overwrite}')
    date_here = start_date
    while date_here <= end_date:
        date_here_str = date_here.strftime('%Y-%m-%d')
        if args.verbose:
            print(
                f'[INFO] -----------------------------------------------------------------------------------------------')
            print(f'[INFO] Working with date: {date_here_str}')
        create_comparison_file_date(date_here, aeronet_file, path_hypstar, path_out, overwrite)
        ##add_mu_to_file_date(date_here)##-> only for test a specific date
        date_here = date_here + timedelta(hours=24)


def create_comparison_file_date(date_here, aeronet_file, path_hypstar, path_out, overwrite):
    # path_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC'
    # # date_here = dt(2023, 9, 26)
    # aeronet_file = '/mnt/c/DATA_LUIS/AERONET_OC/AERONET_NC/20020101_20231111_AAOT.LWN_lev20_15.nc'
    # path_hypstar = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT'

    if not os.path.exists(aeronet_file):
        ##it shouldn't arrive here, file is checked before
        print(f'[ERROR] Aeronet file {aeronet_file} does not exist')
        return

    file_comparison = get_file_comparison_date(path_out, date_here, True)

    if os.path.exists(file_comparison) and not overwrite:
        print(f'[WARNING] Comparison file {file_comparison} already exist. Skipping...')
        return

    path_hypstar_date = get_folder_date(path_hypstar, date_here)
    if not check_hypstar_files(path_hypstar_date, date_here):
        date_here_str = date_here.strftime('%Y-%m-%d')
        print(f'[WARNING] Hyptar files are not available for date: {date_here_str}')
        return

    if args.verbose:
        print(f'[INFO] --> Started comparison file: {file_comparison}')

    from base.anet_nc_reader import AERONETReader
    areader = AERONETReader(aeronet_file)
    ic = INSITUCOMPARISON(file_comparison)

    ic.start_file_w()
    if args.verbose:
        print(f'[INFO] --> Creating AERONET_OC variables...')
    baeronet = ic.create_aeronet_variables(areader, date_here)
    if not baeronet:
        print(f'[ERROR] Comparison file for date: {date_here} could not be done. ')
        os.remove(file_comparison)
        return
    if args.verbose:
        print(f'[INFO] --> Creating HYPSTAR variables...')
    ic.create_hypstar_variables(path_hypstar_date, date_here)
    if args.verbose:
        print(f'[INFO] --> Creating HYPSTAR to AERONET variables..')
    ic.create_hypstar_to_aeronet_variables()
    ic.close_file_w()
    if args.verbose:
        print(f'[INFO] --> Completed')


def check_hypstar_files(path_day, date_here):
    import pytz
    from netCDF4 import Dataset
    time_ini = date_here.replace(tzinfo=pytz.utc, hour=0, minute=0, second=0).timestamp()
    time_end = date_here.replace(tzinfo=pytz.utc, hour=23, minute=59, second=59).timestamp()
    for name in os.listdir(path_day):
        if name.find('L2A_REF') < 0:
            continue
        file_here = os.path.join(path_day, name)
        dataset = Dataset(file_here, 'r')
        time_stamp = float(dataset.variables['acquisition_time'][0])
        if time_ini <= time_stamp <= time_end:
            dataset.close()
            return True
    print(f'[WARNING] HYPSTAR files were not found in folder {path_day} for date {date_here.strftime("%Y-%m-%d")}')
    return False


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
    print(f'[INFO] Adding match-ups to file: {file_comparison}')
    ic = INSITUCOMPARISON(file_comparison)
    ic.add_match_ups_to_file(11)


def make_plots():
    # path_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC'
    # date_here = dt(2023, 9, 26)
    # file_comparison = get_file_comparison_date(path_out, date_here, False)

    file_comparison = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/Comparison_Valid_20230425_20240331.nc'

    file_config = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/config_plot.ini'
    import configparser
    options = configparser.ConfigParser()
    options.read(file_config)
    ic = INSITUCOMPARISON(file_comparison)
    # iplots = INSITU_plots(ic)
    # iplots.plot_from_options(options)

    # if xvariable.find('Li') > 0:
    #     ylabel = r'AERONET-OC Li [μW/(cm$^2$·sr·nm)]'
    #     xlabel = r'HYPSTAR Li [μW/(cm$^2$·sr·nm)]'
    # elif xvariable.find('Lt') > 0:
    #     ylabel = r'AERONET-OC Lt [μW/(cm$^2$·sr·nm)]'
    #     xlabel = r'HYPSTAR Lt [μW/(cm$^2$·sr·nm)]'
    # elif xvariable.find('Lw') > 0:
    #     ylabel = r'AERONET-OC [μW/(cm$^2$·sr·nm)]'
    #     xlabel = r'HYPSTAR Lw [μW/(cm$^2$·sr·nm)]'



    #ic.set_spectra_stats('mu_HYPSTAR_TO_AERONET_Lt_mean','mu_AERONET_Lt_mean','mu_wavelength')

    ic.set_spectra_stats('mu_HYPSTAR_TO_AERONET_Li_mean', 'mu_AERONET_Li_mean', 'mu_wavelength')
    file_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/PLOTS/SpectraComparison_Li.tif'
    legend = ['HYPSTAR - Li', 'AERONET-OC - Li']
    ylabel = r'Li [μW/(cm$^2$·sr·nm)]'
    title = 'Downwelling radiance'
    ic.plot_spectra_stats(file_out, legend, ylabel, title)

    ic.set_spectra_stats('mu_HYPSTAR_TO_AERONET_Lt_mean', 'mu_AERONET_Lt_mean', 'mu_wavelength')
    file_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/PLOTS/SpectraComparison_Lt.tif'
    legend = ['HYPSTAR - Lt', 'AERONET-OC - Lt']
    ylabel = r'Lt [μW/(cm$^2$·sr·nm)]'
    title = 'Upwelling radiance'
    ic.plot_spectra_stats(file_out, legend, ylabel, title)

    ic.set_spectra_stats('mu_HYPSTAR_TO_AERONET_Lw', 'mu_AERONET_Lw', 'mu_wavelength')
    file_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/PLOTS/SpectraComparison_Lw.tif'
    legend = ['HYPSTAR - Lw', 'AERONET-OC - Lw']
    ylabel = r'Lw [μW/(cm$^2$·sr·nm)]'
    title = 'Water-leaving radiance'
    ic.plot_spectra_stats(file_out,legend,ylabel,title)

    #ic.plot_scatterplot(None)
    #ic.plot_all_scatterplots_wl()
    #ic.plot_all_spectra()
    #ic.plot_rho_scatterplot()


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

        # if not args.nodelfiles:
        #     [os.remove(f) for f in list_files[:-1]]
    if args.verbose:
        print(f'[INFO] Concatenated MDB file created: {ncout_file}')


def check_angles():
    print('checking')
    from netCDF4 import Dataset
    path_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT'
    start_date = dt(2023, 3, 23)
    end_date = dt(2023, 10, 31)
    date_here = start_date
    file_out = os.path.join(path_out, 'Azimuth_angles.csv')
    fout = open(file_out, 'w')
    fout.write('Date;Time;ViewingAzimuthAngle;PointAzimuthAngle')
    while date_here <= end_date:
        date_here_str = date_here.strftime('%Y-%m-%d')

        print(f'[INFO] Working with date: {date_here_str}')
        folder = get_folder_date(path_out, date_here)
        for name in os.listdir(folder):
            if name.find('L2A') > 0:
                file = os.path.join(folder, name)
                dataset = Dataset(file)
                vaa = float(dataset.variables['viewing_azimuth_angle'][0])
                paa = float(dataset.variables['pointing_azimuth_angle'][0])
                time_stamp = float(dataset.variables['acquisition_time'][0])
                time_here = dt.utcfromtimestamp(time_stamp)
                time_here_str = time_here.strftime('%H:%M')
                # time_here_str = date_here.strftime(dt.utcfromtimestamp(float(dataset.variables['acquisition_time'][0])))
                line = f'{date_here_str};{time_here_str};{vaa};{paa}'
                print(line)
                fout.write('\n')
                fout.write(line)
                # print(date_here_str,dataset.variables['viewing_azimuth_angle'][0],dataset.variables['pointing_azimuth_angle'][0])
                dataset.close()
        date_here = date_here + timedelta(hours=24)
    fout.close()


def check_hypstar_qf():
    print('checking')
    # file = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT/2023/09/26/HYPERNETS_W_VEIT_L1C_ALL_20230926T0700_20240116T1733_270_v2.0.nc'

    file = '/mnt/c/DATA_LUIS/AERONET_OC/AERONET_NC/20020101_20231111_AAOT.LWN_lev20_15.nc'

    # file = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/2023/04/19/COMPARISON_20230419.nc'
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


def do_test():
    csv_file = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/QUALITY_CONTROL/VEIT/PERIOD_VALID_1/Comparison_VEIT_20230425_20230629.csv'
    path = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/QUALITY_CONTROL/VEIT/PERIOD_VALID_1'
    import pandas as pd
    df = pd.read_csv(csv_file, sep=';')
    for index, row in df.iterrows():
        folder = os.path.join(path, row['FOLDER'])
        name_seq = row['sequence_ref']
        name_img = f'VEIT_{name_seq}_Report.png'
        file_img = os.path.join(path, name_img)
        file_out = os.path.join(folder, name_img)
        if not os.path.exists(file_img):
            print('error')
        shutil.copy(file_img, file_out)


def get_options_comparison(config_file):
    if not os.path.exists(config_file):
        print(f'[ERROR] Config file: {config_file} does not exist')
        return None
    import configparser
    try:
        options = configparser.ConfigParser()
        options.read(config_file)
    except:
        print(f'[ERROR] Parse error reading config file: {config_file}')
        return None
    section = 'comparison_files'
    if not options.has_section(section):
        print(f'[ERROR] Section {section} is required in the config file: {config_file}')
        return None

    required_options = ['aeronet_file', 'path_hypstar', 'path_out', 'start_date', 'end_date']
    options_dict = {}
    for roption in required_options:
        if not options.has_option(section, roption):
            print(f'[ERROR] Option {roption} is required in {section} section of the config file: {config_file}')
            return None
        options_dict[roption] = options[section][roption].strip()

    if not os.path.isfile(options_dict['aeronet_file']):
        print(f'[ERROR] {options_dict["aeronet_file"]} is not a valid file or does not exist')
        return None
    if not os.path.isdir(options_dict['path_hypstar']):
        print(f'[ERROR] {options_dict["path_hypstar"]} is not a valid directory or does not exist')
        return None
    if not os.path.isdir(options_dict['path_out']):
        try:
            os.mkdir(options_dict['path_out'])
        except:
            print(f'[ERROR] {options_dict["path_out"]} is not a valid directory and could not be created')
            return None
    try:
        start_date = dt.strptime(options_dict['start_date'], '%Y-%m-%d')
    except:
        print(f'[ERROR] Start date {options_dict["start_date"]} is not in the correct format %Y-%m-%d')
        return None
    try:
        end_date = dt.strptime(options_dict['end_date'], '%Y-%m-%d')
    except:
        print(f'[ERROR] End date {options_dict["end_date"]} is not in the correct format %Y-%m-%d')
        return None
    if end_date < start_date:
        print(
            f'[ERROR] End date {options_dict["end_date"]} should be greater or equal to start date {options_dict["start_date"]}')
        return None
    options_dict['start_date'] = start_date
    options_dict['end_date'] = end_date

    options_dict['overwrite'] = False
    if options.has_option(section, 'overwrite'):
        sval = options[section]['overwrite'].strip().lower()
        strue = ['true', 't', '1']
        if sval in strue:
            options_dict['overwrite'] = True

    return options_dict


def get_options_concat(config_file):
    if not os.path.exists(config_file):
        print(f'[ERROR] Config file: {config_file} does not exist')
        return None
    import configparser
    try:
        options = configparser.ConfigParser()
        options.read(config_file)
    except:
        print(f'[ERROR] Parse error reading config file: {config_file}')
        return None
    section = 'concatenation'
    required_options = ['path_comparison', 'file_out', 'start_date', 'end_date']
    options_dict = {}
    for roption in required_options:
        if not options.has_option(section, roption):
            print(f'[ERROR] Option {roption} is required in {section} section of the config file: {config_file}')
            return None
        options_dict[roption] = options[section][roption].strip()

    if not os.path.isdir(options_dict['path_comparison']):
        print(f'[ERROR] {options_dict["path_comparison"]} is not a valid directory or does not exist')
        return None

    file_out = options_dict['file_out']
    if os.path.exists(file_out):
        print(f'[ERRR] File out {file_out} already exists. Please use other output file name')
        return None
    path_out = os.path.dirname(file_out)
    name_out = os.path.basename(file_out)
    if not name_out.endswith('.nc'):
        print(f'[ERROR] file_out name {name_out} should be a .nc file')
        return None
    if not os.path.isdir(path_out):
        try:
            os.mkdir(path_out)
        except:
            print(f'[ERROR] {path_out} is not a valid directory and could not be created')
            return None

    start_dates_str = [x.strip() for x in options_dict['start_date'].split(',')]
    end_dates_str = [x.strip() for x in options_dict['end_date'].split(',')]
    if len(start_dates_str) != len(end_dates_str):
        print(f'[ERROR] Number of start_date options should be equal to number of end_date options')
        return None
    start_dates = [None] * len(start_dates_str)
    end_dates = [None] * len(end_dates_str)

    for idx in range(len(start_dates_str)):
        try:
            start_dates[idx] = dt.strptime(start_dates_str[idx], '%Y-%m-%d')
        except:
            print(f'[ERROR] Start date {start_dates_str[idx]} is not in the correct format %Y-%m-%d')
            return None
        try:
            end_dates[idx] = dt.strptime(end_dates_str[idx], '%Y-%m-%d')
        except:
            print(f'[ERROR] End date {end_dates_str[idx]} is not in the correct format %Y-%m-%d')
            return None
        if end_dates[idx] < start_dates[idx]:
            print(f'[ERROR] End date {end_dates[idx]} should be greater or equal to start date {start_dates[idx]}')
            return None

    options_dict['start_date'] = start_dates
    options_dict['end_date'] = end_dates

    return options_dict


def main():
    if args.mode == 'COMPARISON':
        if not args.config_file:
            print(f'[ERROR] Config file is compulsory for option: {args.mode}')
            return
        options = get_options_comparison(args.config_file)
        if options is None:
            return
        create_multiple_comparison_files(options['start_date'], options['end_date'], options['aeronet_file'],
                                         options['path_hypstar'], options['path_out'], options['overwrite'])

    if args.mode == 'CONCAT':
        if not args.config_file:
            print(f'[ERROR] Config file is compulsory for option: {args.mode}')
            return
        options = get_options_concat(args.config_file)
        if options is None:
            return
        make_concatenation(options['path_comparison'],options['file_out'],options['start_date'],options['end_date'])

    elif args.mode == 'TEST':
        # file = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT/2023/05/05/HYPERNETS_W_VEIT_L2A_REF_20230505T1540_20240118T1418_270_v2.0.nc'
        # file = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/2023/05/05/COMPARISON_20230505.nc'
        # file = '/mnt/c/DATA_LUIS/AERONET_OC/AERONET_NC/20020101_20240406_AAOT.LWN_lev20_15.nc'
        # from netCDF4 import Dataset
        # dataset = Dataset(file)
        # for name in dataset.variables:
        #     print(name)
        # # print(dataset.variables['HYPSTAR_epsilon'][:])
        # # print(dataset.variables['rhof'])
        # dataset.close()
        # # do_test()
        #file_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/Comparison_Valid_20230425_20240331.nc'
        #add_mu_to_file(file_out)
        make_plots()


# %%
if __name__ == '__main__':
    # multiple_dates
    ##1. Create files
    # create_multiple_comparison_files()
    ##2. Concatenate
    # make_concatenation()
    ##3. Add mu
    # file_out = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/Comparison_20230425_20230625.nc'
    # add_mu_to_file(file_out)

    ## example for a specific date
    # date_here = dt(2023, 9, 26)
    # create_comparison_file_date(date_here)
    # add_mu_to_file_date(date_here)

    ##PLOTING
    make_plots()

    ##CHECKING HYPSTAR QF (TEST)
    # check_hypstar_qf()
    # check_angles()

    #main()
