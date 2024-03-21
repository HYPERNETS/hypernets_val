import os.path
import shutil
from datetime import datetime as dt
from datetime import timedelta
import argparse

import numpy as np
import pandas as pd

from hypernets_day import HYPERNETS_DAY

parser = argparse.ArgumentParser(description="Creation of insitu nc files")
parser.add_argument('-m', "--mode",
                    choices=['GETFILES', 'CREATEDAYFILES', 'REPORTDAYFILES', 'SUMMARYFILES', 'NCFROMCSV'],
                    required=True)
parser.add_argument('-sd', "--start_date", help="Start date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-ed', "--end_date", help="End date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-st', "--start_time", help="Start time. (HH:MM)")
parser.add_argument('-et', "--end_time", help="End time. (HH:MM)")
parser.add_argument('-i', "--input_path", help="Input path")
parser.add_argument('-o', "--output_path", help="Output path")
parser.add_argument('-c', "--config_path", help="Configuration file path")
parser.add_argument('-site', "--site_name", help="Site name")
parser.add_argument('-sopt', "--summary_options", help="Summary options,separated by '_': csv,nc,copy")
parser.add_argument('-ndays', "--ndays_interval", help="Interval days between start date and end date")
parser.add_argument('-ndel', "--nodelfiles", help="Do not delete temp files.", action="store_true")
parser.add_argument("-ndw", "--nodownload", help="No download (for launching without connection with RBINS).",
                    action="store_true")
parser.add_argument("-ow", "--overwrite", help="Overwrite output file(s).", action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")

args = parser.parse_args()


def test():
    import os
    from netCDF4 import Dataset
    import numpy as np
    start_date = dt(2023, 10, 1)
    end_date = dt(2023, 10, 10)
    work_date = start_date
    path = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/QUALITY_CONTROL/VEIT'
    output_file = os.path.join(path, 'VEIT_20241001_20241010')
    f1 = open(output_file, 'w')
    f1.write('TimeSeq;sza;saa;paa;epsilon')
    while work_date <= end_date:
        yearstr = work_date.strftime('%Y')
        monthstr = work_date.strftime('%m')
        daystr = work_date.strftime('%d')
        datestr = work_date.strftime('%Y%m%d')
        name_file = f'HYPERNETES_W_DAY_{datestr}.nc'
        path_date = os.path.join(path, yearstr, monthstr, daystr, name_file)
        if os.path.exists(path_date):
            dataset = Dataset(path_date)
            seq_ref = dataset.variables['sequence_ref'][:]
            for idx in range(len(seq_ref)):
                seq = seq_ref[idx]
                if np.ma.is_masked(seq):
                    continue
                time_here = dt.utcfromtimestamp(float(seq)).strftime('%Y-%m-%dT%H:%M')
                sza = dataset.variables['l2_solar_zenith_angle'][idx]
                saa = dataset.variables['l2_solar_azimuth_angle'][idx]
                paa = dataset.variables['l2_pointing_azimuth_angle'][idx]
                epsilon = dataset.variables['l2_epsilon'][idx]
                line = f'{time_here};{sza};{saa};{paa};{epsilon}'
                print(line)
                f1.write('\n')
                f1.write(line)

            dataset.close()
        work_date = work_date + timedelta(hours=24)
    f1.close()
    return True


def test_2():
    print('test2')
    # work_date = dt(2023,10,5)
    # hdayfile = hday.get_hypernets_day_file(site, work_date)
    # hdayfile.set_path_images_date(site, work_date)
    # wimages1 = hdayfile.get_water_images(site, work_date, None, None, None)
    # work_date = dt(2023, 10, 7,6,0,0)
    # hdayfile = hday.get_hypernets_day_file(site, work_date)
    # hdayfile.set_path_images_date(site, work_date)
    # wimages2 = hdayfile.get_water_images(site, work_date, None, None, None)
    # end_date = dt(2023, 10, 7,17,0,0)
    # wimages = {}
    # while work_date <= end_date:
    #
    #     wd_str_2 = work_date.strftime('%Y%m%dT%H%M')
    #     wd_str_1 = wd_str_2.replace('20231007','20231005')
    #     if wd_str_1 not in wimages1.keys():
    #         continue
    #     if wd_str_2 not in wimages2.keys():
    #         continue
    #     print(work_date,wd_str_1,wd_str_2)
    #     if wimages1[wd_str_1]['file_img'] is not None and wimages2[wd_str_2]['file_img'] is not None:
    #         wimages[wd_str_1[9:]] = {
    #             'file_1': wimages1[wd_str_1]['file_img'],
    #             'title_1':  wimages1[wd_str_1]['title'],
    #             'file_2': wimages2[wd_str_2]['file_img'],
    #             'title_2': wimages2[wd_str_2]['title'],
    #         }
    #     work_date = work_date + timedelta(minutes=60)
    # print('============================')
    # for wd in wimages:
    #     print(wd,wimages[wd])
    # hdayfile.plot_water_images(wimages)


def make_report_files(input_path, output_path, site, start_date, end_date):
    if args.verbose:
        print(f'[INFO] Started making report')
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)
    hday = HYPERNETS_DAY(input_path, output_path)
    interval = 24
    if args.ndays_interval:
        interval = 24 * int(args.ndays_interval)

    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')
        hdayfile = hday.get_hypernets_day_file(site, work_date)
        if hdayfile is None:
            print(f'[WARNING] HYPERNETS day file for date {work_date} is not available. Skipping...')
            work_date = work_date + timedelta(hours=interval)
            continue
        hdayfile.set_path_images_date(site, work_date)

        for isequence in range(len(hdayfile.sequences)):
            hdayfile.isequence = isequence
            # if isequence==0:
            #     hdayfile.save_report_image(site,args.nodelfiles, args.overwrite)
            hdayfile.save_report_image(site, args.nodelfiles, args.overwrite)
            # hdayfile.save_angle_files(True)

        # hdayfile.save_img_files(True)
        # hdayfile.save_spectra_files(True)
        # hdayfile.save_angle_files(True)
        # hdayfile.get_flags_sequence()
        # hdayfile.get_title()
        # hdayfile.save_report_image(False)
        # create_daily_pdf_report(input_path,output_path,site,work_date)
        work_date = work_date + timedelta(hours=interval)


def create_daily_pdf_report(input_path, output_path, site, date_here):
    import os
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib import pyplot as plt
    hday = HYPERNETS_DAY(input_path, output_path)
    folder_day = hday.get_folder_date(site, date_here)
    date_here_str = date_here.strftime('%Y%m%d')
    file_pdf = os.path.join(folder_day, f'Report_{site}_{date_here_str}.pdf')
    pdf = PdfPages(file_pdf)
    index = 0
    exist_image = True
    while exist_image:
        fimage = os.path.join(folder_day, f'ReportImage_{index}.tif')
        exist_image = os.path.exists(fimage)
        if exist_image:
            print(fimage)
            plt.close()
            plt.figure(figsize=(10, 18))
            plt.imshow(plt.imread(fimage))
            plt.axis('off')
            pdf.savefig(dpi=300)

        index = index + 1

    pdf.close()


def make_get_files(input_path, site, start_date, end_date):
    if args.verbose:
        print(f'[INFO] Started getting files')
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)

    hday = HYPERNETS_DAY(input_path, None)
    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')
        hday.get_files_date(site, work_date)
        work_date = work_date + timedelta(hours=24)


def make_create_dayfiles(input_path, output_path, site, start_date, end_date):
    if args.verbose:
        print(f'[INFO] Started creating files')
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)
    interval = 24
    if args.ndays_interval:
        interval = 24 * int(args.ndays_interval)
    hday = HYPERNETS_DAY(input_path, output_path)
    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')
        if args.nodownload:
            hday.get_files_date_nodownload(site, work_date)
        else:
            hday.get_files_date(site, work_date)

        if len(hday.files_dates) == 0:
            print(f'[WARNING] No data files found for the sequences on date: {work_date}. Skipping...')
            work_date = work_date + timedelta(hours=interval)
            continue

        if args.verbose:
            print(f'[INFO] Number of sequences: {len(hday.files_dates)}')
        nseq = hday.start_file_date_complete(site, work_date, args.overwrite)
        if nseq < 0:
            work_date = work_date + timedelta(hours=interval)
            continue
        else:
            if args.verbose:
                print(f'[INFO] Valid sequences: {nseq}')
        hday.set_data(site, work_date)
        hday.close_datafile_complete()
        work_date = work_date + timedelta(hours=interval)

    # hday = HYPERNETS_DAY(None)
    # hday.get_files_date(site, date_here)
    # hday.start_file_date_complete(site, date_here, True)
    # hday.set_data(site, date_here)
    # hday.close_datafile_complete()
    # file_rgb = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT/2023/03/25/HYPERNETS_W_VEIT_IMG_20230325T1040_20240118T0338_009_40_90_v2.0.jpg'
    # from PIL import Image
    # import numpy as np
    # img = Image.open(file_rgb)
    # red_array_orig = np.array(img.getdata(0)).reshape(img.size)
    # red_array = np.array(img.getdata(0)).reshape(img.size).astype(np.uint8)


def make_flagged_nc_from_csv(config_file):
    if args.verbose:
        print(f'[INFO] Started concatenaded nc file from csv using configuration file: {config_file}')

    import configparser
    from hypernets_day_file import HYPERNETS_DAY_FILE
    options = configparser.ConfigParser()
    options.read(config_file)
    if not options.has_section('basic'):
        print(f'[ERROR] Basic section should be included in the configuration file')
        return
    input_csv = None
    if options.has_option('basic', 'input_csv'):
        input_csv = options['basic']['input_csv'].strip()
    if input_csv is None:
        print(f'[ERROR] input_csv should be included in the basic section of the configuration file')
        return
    if not os.path.exists(input_csv):
        print(f'[ERROR] File {input_csv} does not exist')
        return

    if options.has_option('basic', 'output_path'):
        output_path = options['basic']['output_path'].strip()
    else:
        output_path = os.path.dirname(input_csv)
    if not os.path.isdir(output_path):
        try:
            os.mkdir(output_path)
        except:
            print(f'[ERROR] Output path {output_path} does not exist')
            return

    col_sequence = 'sequence_ref'
    col_file_nc = 'file_nc'
    col_isequence = 'isequence'
    flags = get_flags(options)
    name_csv = input_csv.split('/')[-1]
    name_csv = name_csv[:-4]
    file_nc = os.path.join(output_path, f'CSV_NC_{name_csv}.nc')

    df = pd.read_csv(input_csv, sep=';')

    file_nc_prev = ''
    dataset_w = None
    hdayfile = None
    sindices = []
    pindices = []
    index_w = 0
    for index, row in df.iterrows():
        file_nc_here = str(row[col_file_nc])
        if not os.path.exists(file_nc_here):
            continue
        isequence = int(row[col_isequence])

        if file_nc_here != file_nc_prev:
            if dataset_w is not None:  ##set data
                dataset_w = hdayfile.set_data_dataset_w(dataset_w, sindices, index_w)
                for flag in flags:
                    print(pindices)
                    print(flags)
                    print(flag)
                    flag_array = get_flag_array(pindices, df, flags, flag)
                    print(flag_array)
                    dataset_w = hdayfile.set_data_flag(dataset_w, index_w, flag, flag_array)
                index_w = index_w + len(sindices)
                hdayfile = HYPERNETS_DAY_FILE(file_nc_here, None)
                sindices = [isequence]
                pindices = [int(index)]
            else:
                hdayfile = HYPERNETS_DAY_FILE(file_nc_here, None)
                dataset_w = hdayfile.start_dataset_w(file_nc)
                dataset_w = hdayfile.add_flag_variables(dataset_w, flags)
                sindices = [isequence]
                pindices = [int(index)]
            file_nc_prev = file_nc_here
        else:
            sindices.append(isequence)
            pindices.append(int(index))

    if dataset_w is not None:
        dataset_w.close()


def get_flags(options):
    flags = {}
    if options.has_section('flags'):
        index = 0
        while options.has_option('flags', f'flag_{index}.name'):
            name = options['flags'][f'flag_{index}.name'].strip()
            flag_values = []
            flag_meanings = []
            if options.has_option('flags', f'flag_{index}.values'):
                flag_values_str = options['flags'][f'flag_{index}.values'].strip()
                flag_values = [int(x) for x in flag_values_str.split(',')]
            if options.has_option('flags', f'flag_{index}.values'):
                flag_meanings_str = options['flags'][f'flag_{index}.meanings'].strip()
                flag_meanings = [str(x) for x in flag_meanings_str.split(',')]
            flags[name] = {
                'values': flag_values,
                'meanings': flag_meanings
            }
            index = index + 1
    return flags


def get_flag_array(indices, df, flags, flag):
    print('Columna: ',df.columns)
    print('Indices: ',indices)
    print(type(indices))
    print('Flag',flag)
    data_flag = df.loc[indices, flag]
    try:
        data_flag_array = np.array(data_flag).astype(np.int64)
    except:
        for idx in range(len(flags[flag]['values'])):
            val = flags[flag]['values'][idx]
            meaning = flags[flag]['meanings'][idx]
            data_flag[data_flag == meaning] = val
        data_flag_array = np.array(data_flag).astype(np.int64)

    return data_flag_array


def make_summary_files(input_path, output_path, site, start_date, end_date, start_time, end_time, options):
    if args.verbose:
        print(f'[INFO] Started summary files with options: {options}')

    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)

    hday = HYPERNETS_DAY(None, input_path)
    interval = 24
    if args.ndays_interval:
        interval = 24 * int(args.ndays_interval)

    if not os.path.isdir(output_path):
        try:
            os.mkdir(output_path)
        except:
            print(f'[ERROR] Output path {output_path} is not a valid directory and could not be created')
            return

    start_date_str = start_date.strftime('%Y%m%d')
    end_date_str = end_date.strftime('%Y%m%d')
    dataset_w = None
    index_w = 0
    df_csv = None
    col_names_csv = None

    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')
        hdayfile = hday.get_hypernets_day_file(site, work_date)
        if hdayfile is None:
            print(f'[WARNING] HYPERNETS day file for date {work_date} is not available. Skipping...')
            work_date = work_date + timedelta(hours=interval)
            continue
        work_time_str = work_date.strftime('%Y-%m-%d')
        start_work_time = dt.strptime(f'{work_time_str}T{start_time}', '%Y-%m-%dT%H:%M')
        end_work_time = dt.strptime(f'{work_time_str}T{end_time}', '%Y-%m-%dT%H:%M')
        sequences, sequence_indices = hdayfile.get_sequences_interval(start_work_time, end_work_time)
        if 'copy' in options:
            report_files = hdayfile.get_report_files_interval(sequences, site, start_time, end_time)
            if args.verbose:
                print(f'[INFO] -> Copy {len(report_files)} report files')
            if len(report_files) > 0:
                for rinfile in report_files:
                    routfile = os.path.join(output_path, os.path.basename(rinfile))
                    shutil.copy(rinfile, routfile)
        if 'nc' in options and sequence_indices is not None:
            if dataset_w is None:
                file_nc = os.path.join(output_path, f'Comparison_{site}_{start_date_str}_{end_date_str}.nc')
                dataset_w = hdayfile.start_dataset_w(file_nc)
            print('Setting data->', len(sequence_indices))
            dataset_w = hdayfile.set_data_dataset_w(dataset_w, sequence_indices, index_w)
            index_w = index_w + len(sequence_indices)

        if 'csv' in options and sequence_indices is not None:
            if col_names_csv is None:
                col_names_csv = hdayfile.get_csv_col_names()
                df_csv = hdayfile.get_dataframe_lines(sequence_indices, col_names_csv)
            else:
                df_here = hdayfile.get_dataframe_lines(sequence_indices, col_names_csv)
                df_csv = pd.concat([df_csv, df_here])

        work_date = work_date + timedelta(hours=interval)

    if df_csv is not None:
        file_csv = os.path.join(output_path, f'Comparison_{site}_{start_date_str}_{end_date_str}.csv')
        df_csv.to_csv(file_csv, sep=';')
    if dataset_w is not None:
        dataset_w.close()


def get_start_and_end_dates():
    start_date = dt.now()
    if args.start_date:
        try:
            start_date = dt.strptime(args.start_date, '%Y-%m-%d')
        except:
            print(f'[ERROR] Start date is not in valid format YYYY-mm-DD')
            return None, None
    if args.end_date:
        try:
            end_date = dt.strptime(args.end_date, '%Y-%m-%d')
        except:
            print(f'[ERROR] End date is not in valid format YYYY-mm-DD')
            return None, None
    else:
        end_date = start_date

    return start_date, end_date


def get_start_and_end_times():
    start_time = '00:00'
    end_time = '23:59'
    if args.start_time:
        try:
            start_time_p = dt.strptime(args.start_time, '%H:%M')
            start_time = start_time_p.strftime('%H:%M')
        except:
            pass
    if args.end_time:
        try:
            end_time_p = dt.strptime(args.end_time, '%H:%M')
            end_time = end_time_p.strftime('%H:%M')
        except:
            pass

    return start_time, end_time


def main():
    print('STARTED')
    # b = test()
    # if b:
    #     return
    start_date, end_date = get_start_and_end_dates()
    if start_date is None:
        return
    start_time, end_time = get_start_and_end_times()

    site = 'VEIT'
    if args.site_name:
        site = args.site_name
    input_path = None
    output_path = None

    if args.input_path:
        input_path = args.input_path
        if args.verbose:
            print(f'[INFO] Input path set to: {input_path}')
    if args.output_path:
        output_path = args.output_path
        if args.verbose:
            print(f'[INFO] Output path set to: {output_path}')

    summary_options = ['copy']
    if args.summary_options:
        summary_options = args.summary_options.split('_')

    if args.mode == 'GETFILES':
        make_get_files(input_path, site, start_date, end_date)

    if args.mode == 'CREATEDAYFILES':
        make_create_dayfiles(input_path, output_path, site, start_date, end_date)

    if args.mode == 'REPORTDAYFILES':
        make_report_files(input_path, output_path, site, start_date, end_date)

    if args.mode == 'SUMMARYFILES':
        make_summary_files(input_path, output_path, site, start_date, end_date, start_time, end_time, summary_options)

    if args.mode == 'NCFROMCSV':
        make_flagged_nc_from_csv(args.config_path)


# %%
if __name__ == '__main__':
    main()
