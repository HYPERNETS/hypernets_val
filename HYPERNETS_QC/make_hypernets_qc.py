from datetime import datetime as dt
from datetime import timedelta
import argparse
from hypernets_day import HYPERNETS_DAY

parser = argparse.ArgumentParser(description="Creation of insitu nc files")
parser.add_argument('-m', "--mode", choices=['GETFILES', 'CREATEDAYFILES', 'REPORTDAYFILES'], required=True)
parser.add_argument('-sd', "--start_date", help="Start date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-ed', "--end_date", help="End date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-i', "--input_path", help="Input path")
parser.add_argument('-o', "--output_path", help="Output path")
parser.add_argument('-site', "--site_name", help="Site name")
# parser.add_argument('-nd', "--nodelfiles", help="Do not delete temp files.", action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")

args = parser.parse_args()

def make_report_files(input_path,output_path,site,start_date,end_date):
    if args.verbose:
        print(f'[INFO] Started making report')
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)
    hday = HYPERNETS_DAY(input_path,output_path)
    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')
        hdayfile = hday.get_hypernets_day_file(site,work_date)
        if hdayfile is None:
            print(f'[WARNING] HYPERNETS day file for date {work_date} is not available. Skipping...')
        for isequence in range(len(hdayfile.sequences)):
            hdayfile.isequence = isequence
            hdayfile.save_report_image(False,False)
            #hdayfile.save_angle_files(True)
            #hdayfile.save_spectra_files(True)
        #hdayfile.save_img_files(True)
        #hdayfile.save_spectra_files(True)
        #hdayfile.save_angle_files(True)
        #hdayfile.get_flags_sequence()
        #hdayfile.get_title()
        #hdayfile.save_report_image(False)
        #create_daily_pdf_report(input_path,output_path,site,work_date)
        work_date = work_date + timedelta(hours=24)

def create_daily_pdf_report(input_path,output_path,site,date_here):
    import os
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib import pyplot as plt
    hday = HYPERNETS_DAY(input_path,output_path)
    folder_day = hday.get_folder_date(site,date_here)
    date_here_str = date_here.strftime('%Y%m%d')
    file_pdf = os.path.join(folder_day,f'Report_{site}_{date_here_str}.pdf')
    pdf = PdfPages(file_pdf)
    index = 0
    exist_image = True
    while exist_image:
        fimage = os.path.join(folder_day,f'ReportImage_{index}.tif')
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



def make_get_files(input_path,output_path,site, start_date, end_date):
    if args.verbose:
        print(f'[INFO] Started getting files')
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)

    hday = HYPERNETS_DAY(input_path,output_path)
    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')
        hday.get_files_date(site, work_date)
        work_date = work_date + timedelta(hours=24)


def make_create_dayfiles(input_path,output_path,site, start_date,end_date):
    if args.verbose:
        print(f'[INFO] Started creating files')
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)


    hday = HYPERNETS_DAY(input_path,output_path)
    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')
        hday.get_files_date(site, work_date)
        hday.start_file_date_complete(site, work_date, True)
        hday.set_data(site, work_date)
        hday.close_datafile_complete()
        work_date = work_date + timedelta(hours=24)

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


def main():
    print('STARTED')
    start_date, end_date = get_start_and_end_dates()
    if start_date is None:
        return
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


    if args.mode == 'GETFILES':
        make_get_files(input_path,output_path,site, start_date, end_date)

    if args.mode == 'REPORTDAYFILES':
        make_report_files(input_path,output_path,site,start_date,end_date)

    if args.mode == 'CREATEDAYFILES':
        make_create_dayfiles(input_path,output_path,site,start_date,end_date)

# %%
if __name__ == '__main__':
    main()
