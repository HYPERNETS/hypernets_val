import argparse
import os.path
from datetime import datetime as dt
from datetime import timedelta
from INSITU_hypernets import INSITU_HYPERNETS_DAY

parser = argparse.ArgumentParser(
    description="Download in situ files.")

parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-o', "--output", help="Output directory")
parser.add_argument('-site', "--sitename", help="Site name. Only required with --listdates")
parser.add_argument('-sd', "--start_date", help="The Start Date - format YYYY-MM-DD ")
parser.add_argument('-ed', "--end_date", help="The End Date - format YYYY-MM-DD ")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-check', "--check_mode", help="Check mode.", action="store_true")
parser.add_argument('-meta', "--download_metadata", help="Option do download metadata", action="store_true")

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

    if args.check_mode:
        print('[INFO] Entering in check mode....')
        make_check(start_date, end_date, site, output_folder)
        return
    if args.download_metadata:
        print('[INFO] Entering metadata mode...')
        make_download_metadata(start_date, end_date, site, output_folder)
        return

    make_download(start_date, end_date, site, output_folder)


def make_check(start_date, end_date, site, output_folder):
    ih = INSITU_HYPERNETS_DAY(None, None, args.verbose)
    # ih.find_ref = 'HYPERNETS_W_SITE_L1C_ALL*'
    # ih.find_ref = 'HYPERNETS_W_SITE_L2A_REF*'
    output_file = os.path.join(output_folder, f'NFiles_{site}.csv')
    fout = open(output_file, 'w')
    fout.write('Date;NFiles_L1A;NFiles_L1C_ALL;NFiles_L2A_REF')
    date_download = start_date
    while date_download <= end_date:
        date_download_str = date_download.strftime('%Y-%m-%d')
        # ih.find_ref = '*.jpg'
        ih.find_ref = 'HYPERNETS_W_SITE_L1A_*'
        files_l1a = ih.get_files_download(date_download, site)
        ih.find_ref = 'HYPERNETS_W_SITE_L1C_ALL*'
        files_l1c = ih.get_files_download(date_download, site)
        # ih.transfer_files_to_output_folder_via_ssh(files_l1a,output_folder)
        ih.find_ref = 'HYPERNETS_W_SITE_L2A_REF*'
        files_l2a = ih.get_files_download(date_download, site)
        line = f'{date_download_str};{len(files_l1a)};{len(files_l1c)};{len(files_l2a)}'
        fout.write('\n')
        fout.write(line)
        # if args.verbose:
        #     print(f'[INFO] Date: {date_download_str} Files available for download: {len(files_l1a)}')
        if args.verbose:
            print(
                f'[INFO] Date: {date_download_str} L1AFiles: {len(files_l1a)} L1CFiles: {len(files_l1c)} L2Files: {len(files_l2a)}')

        date_download = date_download + timedelta(hours=24)

    fout.close()


def make_download_metadata(start_date, end_date, site, output_folder):
    ih = INSITU_HYPERNETS_DAY(None, None, args.verbose)

    date_download = start_date
    while date_download <= end_date:
        sequence_folders = ih.get_sequences_day_ssh(site, date_download)
        if len(sequence_folders) == 0:
            if args.verbose:
                print(f'[WARNING] No sequences are available for site: {site} and date: {date_download}')
            date_download = date_download + timedelta(hours=24)
            continue
        output_folder_site = get_folder_new(output_folder, site)
        output_folder_date = get_folder_date(output_folder_site, date_download)
        if output_folder_date is None:
            print(f'[ERROR] Output folder date does not exist and could not be created')
            date_download = date_download + timedelta(hours=24)
            continue
        for seq in sequence_folders:
            print(f'[INFO] Donwload metadata file for {site}/{seq} to {output_folder_date}')
            ih.download_sequence_metadata(site, seq, output_folder_date)
            file_metadata = os.path.join(output_folder_date, 'metadata.txt')
            file_metadata_new = os.path.join(output_folder_date, f'{seq}_metadata.txt')
            os.rename(file_metadata, file_metadata_new)
        date_download = date_download + timedelta(hours=24)


def make_download(start_date, end_date, site, output_folder):
    ih = INSITU_HYPERNETS_DAY(None, None, args.verbose)

    date_download = start_date
    while date_download <= end_date:
        if site == 'JSIT':
            sequence_folders = ih.get_sequences_date_l2(site, date_download)
        else:
            sequence_folders = ih.get_sequences_day_ssh(site, date_download)
        if len(sequence_folders) == 0:
            if args.verbose:
                print(f'[WARNING] No sequences are available for site: {site} and date: {date_download}')
            date_download = date_download + timedelta(hours=24)
            continue
        # print(sequence_folders)



        ##water sites
        if site != 'JSIT':
            files_download_all = None
            ih.find_ref = 'HYPERNETS_W_SITE_L1C_ALL*'
            files_download_l1 = ih.get_files_download(date_download, site)
            if files_download_l1 is not None:
                files_download_all = files_download_l1
            ih.find_ref = 'HYPERNETS_W_SITE_L2A_REF*'
            files_download_l2 = ih.get_files_download(date_download, site)
            if files_download_l2 is not None:
                if files_download_all is None:
                    files_download_all = files_download_l2
                else:
                    files_download_all = files_download_all + files_download_l2

            ih.find_ref = 'HYPERNETS_W_SITE_IMG*jpg'
            files_download_img = ih.get_files_download(date_download, site)
            if files_download_img is not None:
                if files_download_all is None:
                    files_download_all = files_download_img
                else:
                    files_download_all = files_download_all + files_download_img

            if files_download_all is None:
                print(f'[WARNING] Files are not available for downloading. Skipping date: {date_download}')
                date_download = date_download + timedelta(hours=24)
                continue
            output_folder_site = get_folder_new(output_folder, site)
            output_folder_date = get_folder_date(output_folder_site, date_download)
            if output_folder_date is None:
                print(f'[ERROR] Output folder date does not exist and could not be created')
                date_download = date_download + timedelta(hours=24)
                continue
            if args.verbose:
                print(f'[INFO] Output folder date: {output_folder_date}')
                print(f'[INFO] Files available for download: {len(files_download_all)}')
            ih.transfer_files_to_output_folder_via_ssh(files_download_all, output_folder_date)

        if site == 'JSIT':
            output_folder_site = get_folder_new(output_folder, site)
            output_folder_date = get_folder_date(output_folder_site, date_download)
            if output_folder_date is None:
                print(f'[ERROR] Output folder date does not exist and could not be created')
                date_download = date_download + timedelta(hours=24)
                continue
            for seq in sequence_folders:
                print('[INFO] Downloading images for sequence:', seq)
                list_img = ih.get_files_img_download_land(site, seq)
                if list_img is not None and len(list_img) > 0:
                    ih.transfer_files_to_output_folder_via_ssh_land(list_img, output_folder_date)
                    for file_img_orig in list_img:
                        file_img = os.path.join(output_folder_date,os.path.basename(file_img_orig))
                        if os.path.exists(file_img):
                            name_new = get_name_new_file_img(file_img, site, seq)
                            file_img_new = os.path.join(output_folder_date, name_new)
                            os.rename(file_img, file_img_new)
            files_download_all = None
            ih.find_ref = 'HYPERNETS_L_SITE_L1C_ALL*'
            files_download_l1 = ih.get_files_download_land(date_download, site)
            #print('-->',files_download_l1)
            if files_download_l1 is not None:
                files_download_all = files_download_l1
            ih.find_ref = 'HYPERNETS_L_SITE_L2A_REF*'
            files_download_l2 = ih.get_files_download_land(date_download, site)
            if files_download_l2 is not None:
                if files_download_all is None:
                    files_download_all = files_download_l2
                else:
                    files_download_all = files_download_all + files_download_l2
            if files_download_all is None:
                print(f'[WARNING] Files are not available for downloading. Skipping date: {date_download}')
                date_download = date_download + timedelta(hours=24)
                continue
            if args.verbose:
                print(f'[INFO] Output folder date: {output_folder_date}')
                print(f'[INFO] Files available for download: {len(files_download_all)}')
            ih.transfer_files_to_output_folder_via_ssh_land(files_download_all, output_folder_date)


        if len(sequence_folders) > 0:
            save_sequence_list(sequence_folders, output_folder_date)

        date_download = date_download + timedelta(hours=24)


def get_name_new_file_img(file_img,site,seq):
    type = 'W'
    if site=='JSIT':
        type = 'L'
    name_split = os.path.basename(file_img)[:-4].split('_')
    picture = name_split[1]
    zenith = str(int(name_split[4]))
    azimuth = str(int(name_split[2]))
    date_img =  dt.fromtimestamp(os.path.getmtime(file_img))
    date_img_str = date_img.strftime('%Y%m%dT%H%M')
    name_new = f'HYPTERNETS_{type}_{site}_IMG_{seq[3:-2]}_{date_img_str}_{picture}_{zenith}_{azimuth}_v2.0.jpg'
    return name_new

def save_sequence_list(sequence_folders, output_folder_date):
    file_out = os.path.join(output_folder_date, 'sequence_list.txt')
    f1 = open(file_out, 'w')
    started = True
    for seq in sequence_folders:
        if not started:
            f1.write('\n')
        f1.write(seq)
        started = False
    f1.close()


def get_folder_new(output_folder_base, new_folder):
    if output_folder_base is None:
        return None
    folder_new = os.path.join(output_folder_base, new_folder)
    if not os.path.exists(folder_new):
        try:
            os.mkdir(folder_new)
        except:
            return None
    return folder_new


def get_folder_date(output_folder_base, date_here):
    folder_year = get_folder_new(output_folder_base, date_here.strftime('%Y'))
    folder_month = get_folder_new(folder_year, date_here.strftime('%m'))
    folder_day = get_folder_new(folder_month, date_here.strftime('%d'))
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
