import argparse
import os.path
from netCDF4 import Dataset
from datetime import datetime as dt

parser = argparse.ArgumentParser(description="Adding flag band to extract files")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-i', "--input_dir", help="Input directory with sat extracts files")
parser.add_argument('-o', "--output_dir", help="Ouput directory with sat extracts files")
parser.add_argument('-f', "--flag_band", help="Bands to be added (separated by ,)")
args = parser.parse_args()


def main():
    print('[INFO] Started add flag mask')
    if not os.path.isdir(args.output_dir) or not os.path.exists(args.output_dir):
        return

    flag_bands = args.flag_band.split(',')
    if args.verbose:
        print(f'[INFO] Bands to be added:')
        for flag in flag_bands:
            print(f'[INFO]->{flag}')

    file_dates = get_file_dates(flag_bands)

    if file_dates is None or len(file_dates) == 0:
        print(f'[ERROR] Bands to be added were not found in any file in input directory')
        return
    if args.verbose:
        print('[INFO] Retrieved input file dates...')



    output_folder = os.path.join(os.path.dirname(args.output_dir), args.output_dir.split('/')[-1] + '_noflag')
    if not os.path.exists(output_folder):
        os.mkdir((output_folder))





    for fname in os.listdir(args.output_dir):
        if fname.endswith('nc'):
            if args.verbose:
                print('------------------')
                print(f'Working with file: {fname}')

            file_nc = os.path.join(args.output_dir, fname)
            dataset = Dataset(file_nc, mode='a')
            date_here = dt.fromtimestamp(int(dataset.variables['satellite_time'][0]))
            date_here_str = date_here.strftime('%Y%m%d')
            ref_here = f'{dataset.satellite.upper()}{dataset.platform.upper()}{dataset.sensor.upper()}'
            key_here = f'{date_here_str}_{ref_here}'
            if key_here in file_dates.keys():
                dataset_orig = Dataset(file_dates[key_here]['file_path'])
                for flag_band in flag_bands:
                    if flag_band in dataset.variables:
                        if args.verbose:
                            print(f'[INFO] Variable: {flag_band} already available. Skipping...')
                        continue
                    if args.verbose:
                        print(f'[INFO] Adding variable: {flag_band}')
                    satellite_flag = dataset.createVariable(flag_band, 'f4', ('satellite_id', 'rows', 'columns'),
                                                            fill_value=-999, zlib=True, complevel=6)
                    var_orig = dataset_orig.variables[flag_band]
                    satellite_flag[0, :, :] = var_orig[0, :, :]
                    if 'description' in var_orig.ncattrs():
                        satellite_flag.description = var_orig.description
                    if 'flag_masks' in var_orig.ncattrs():
                        satellite_flag.flag_masks = var_orig.flag_masks
                    if 'flag_meanings' in var_orig.ncattrs():
                        satellite_flag.flag_meanings = var_orig.flag_meanings
                dataset_orig.close()
                dataset.close()
            else:
                dataset.close()
                output_path = os.path.join(output_folder, fname)
                os.rename(file_nc, output_path)


def get_file_dates(flag_bands):
    if not os.path.isdir(args.input_dir) or not os.path.exists(args.input_dir):
        return None

    file_dates = {}

    for fname in os.listdir(args.input_dir):
        if fname.endswith('nc'):
            file_nc = os.path.join(args.input_dir, fname)

            dataset = Dataset(file_nc)
            have_bands = True
            for flag_band in flag_bands:
                if not flag_band in dataset.variables:
                    have_bands = False

            if have_bands:
                date_here = dt.fromtimestamp(int(dataset.variables['satellite_time'][0]))
                date_here_str = date_here.strftime('%Y%m%d')
                ref_here = f'{dataset.satellite.upper()}{dataset.platform.upper()}{dataset.sensor.upper()}'
                key_here = f'{date_here_str}_{ref_here}'
                file_dates[key_here] = {
                    'file_path': file_nc,
                }

    return file_dates


if __name__ == '__main__':
    main()
