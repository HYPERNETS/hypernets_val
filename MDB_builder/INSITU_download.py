import argparse
import os.path
from datetime import datetime as dt
from datetime import timedelta
from INSITU_hypernets import INSITU_HYPERNETS_DAY

parser = argparse.ArgumentParser(
    description="Download in situ files.")

parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-o', "--output",help="Output directory")
parser.add_argument('-site', "--sitename", help="Site name. Only required with --listdates")
parser.add_argument('-sd', "--start_date", help="The Start Date - format YYYY-MM-DD ")
parser.add_argument('-ed', "--end_date", help="The End Date - format YYYY-MM-DD ")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")

args = parser.parse_args()


def main():
    print('[INFO] Started donwload...')
    if args.config_file:
        print('[ERROR] Option config file is not implemented yet')
        return
    start_date, end_date = get_dates_from_arg()
    if start_date is None or end_date is None:
        return
    output_folder = args.output
    if not os.path.exists(output_folder):
        try:
            os.mkdir(output_folder)
        except:
            print(f'[ERROR] Output folder does not exists and could not be created')
            return
    site = args.sitename
    make_download(start_date,end_date,site,output_folder)


def make_download(start_date,end_date,site,output_folder):
    ih = INSITU_HYPERNETS_DAY(None,None,args.verbose)
    date_download = start_date
    while date_download<=end_date:

        files_download = ih.get_files_download(date_download,site)
        if files_download is None:
            if args.verbose:
                print(f'[WARNING] No files are available for site: {site} and date: {date_download}')
            date_download = date_download + timedelta(hours=24)
            continue

        output_folder_site = get_folder_new(output_folder,site)
        output_folder_date = get_folder_date(output_folder_site,date_download)
        if output_folder_date is None:
            print(f'[ERROR] Output folder date does not exist and could not be created')
            date_download = date_download + timedelta(hours=24)
            continue
        if args.verbose:
            print(f'[INFO] Output folder date: {output_folder_date}')
            print(f'[INFO] Files available for download: {len(files_download)}')

        ih.transfer_files_to_output_folder_via_ssh(files_download,output_folder_date)

        date_download = date_download + timedelta(hours=24)


def get_folder_new(output_folder_base,new_folder):
    if output_folder_base is None:
        return None
    folder_new = os.path.join(output_folder_base,new_folder)
    if not os.path.exists(folder_new):
        try:
            os.mkdir(folder_new)
        except:
            return None
    return folder_new

def get_folder_date(output_folder_base,date_here):
    folder_year = get_folder_new(output_folder_base,date_here.strftime('%Y'))
    folder_month = get_folder_new(folder_year,date_here.strftime('%m'))
    folder_day = get_folder_new(folder_month,date_here.strftime('%d'))
    return folder_day


def get_dates_from_arg():

    start_date = None
    end_date = None
    if args.start_date:
        try:
            start_date = dt.strptime(args.start_date, '%Y-%m-%d')
        except:
            try:
                tdelta = int(args.start_date)
                start_date = dt.now() + timedelta(days=tdelta)
                start_date = start_date.replace(hour=12, minute=0, second=0, microsecond=0)
            except:
                print(f'[ERROR] Start date {args.start_date} is not in the correct format: YYYY-mm-dd or integer')
    if args.end_date:
        try:
            end_date = dt.strptime(args.end_date, '%Y-%m-%d')
        except:
            try:
                tdelta = int(args.end_date)
                end_date = dt.now() + timedelta(days=tdelta)
                end_date = end_date.replace(hour=12, minute=0, second=0, microsecond=0)
            except:
                print(f'[ERROR] End date {args.end_date} is not in the correct format: YYYY-mm-dd or integer')
    if args.start_date and not args.end_date:
        end_date = start_date

    return start_date, end_date

# %%
if __name__ == '__main__':
    main()