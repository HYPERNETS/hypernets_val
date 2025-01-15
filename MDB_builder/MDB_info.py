import argparse
import os.path
from datetime import timedelta
from netCDF4 import Dataset

parser = argparse.ArgumentParser(
    description="Obtaining information for running MDB_builder.")

parser.add_argument('-m', "--mode", help='Mode option', choices=["instrument_id"], required=True)
parser.add_argument('-i', "--input_path", help="Input path.")
parser.add_argument('-o', "--output", help="Output file.")
parser.add_argument('-sd', "--start_date", help="Start date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-ed', "--end_date", help="End date. Optional with --listdates (YYYY-mm-dd)")

args = parser.parse_args()


def main():
    print(f'[INFO] Started MDB_info!')
    if args.mode == 'instrument_id':
        if not check_required_params(['input_path', 'start_date', 'end_date']):
            return
        start_date, end_date = get_dates_from_arg()
        if start_date is None or end_date is None:
            return
        if not os.path.isdir(args.input_path):
            print(f'[ERROR] {args.input_path} does not exist or is not a valid directory')
            return
        get_info_instrument_id(args.input_path, start_date, end_date)


def get_info_instrument_id(input_path, start_date, end_date):
    date_ref = start_date
    ids_list = {}
    while date_ref <= end_date:
        print(f'[INFO] Checking date: {date_ref.strftime("%Y-%m-%d")}')
        yyyy = date_ref.strftime('%Y')
        mm = date_ref.strftime('%m')
        dd = date_ref.strftime('%d')
        folder_date = os.path.join(input_path, yyyy, mm, dd)
        if os.path.isdir(folder_date):
            for name in os.listdir(folder_date):
                if name.find('L2A_REF') > 0:
                    file_nc = os.path.join(folder_date, name)
                    dataset = Dataset(file_nc)
                    if 'instrument_id' in dataset.ncattrs():
                        id_here = str(dataset.instrument_id)
                        if id_here not in ids_list.keys():
                            nwl = dataset.variables['wavelength'].shape[0]
                            ids_list[id_here] = nwl
                    dataset.close()
        date_ref = date_ref + timedelta(hours=24)
    for id_here in ids_list:
        print(f'[INFO] Instrument id: {id_here} Number of wavelengths: {ids_list[id_here]}')


def check_required_params(param_list):
    b = True
    for param in param_list:
        if not args.__dict__[param]:
            print(f'[ERROR] {param} is required for mode {args.mode}')
            b = False
    return b


def get_dates_from_arg():
    from datetime import datetime as dt
    from datetime import timedelta
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
