import argparse
import os.path
from netCDF4 import Dataset
from datetime import datetime as dt

parser = argparse.ArgumentParser(description="Adding flag band to extract files. ")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-i', "--input_dir", help="Input directory with sat extracts files with flag to be added")
parser.add_argument('-o', "--output_dir", help="Ouput directory with sat extracts files")
parser.add_argument('-f', "--flag_band", help="Bands to be added (separated by ,)")
args = parser.parse_args()


def main():
    print('[INFO] Started add flag mask')
    if not os.path.isdir(args.output_dir) or not os.path.exists(args.output_dir):
        return
    if not args.flag_band:
        return

    if args.flag_band=='CHECK_EXTRACTS': ##MODE TO CHECK THAT THE EXTRACTS IN INPUT AND OUTPUT PATHS ARE THE SAME
        input_dir = args.input_dir
        output_dir = args.output_dir
        if not os.path.isdir(input_dir) or not os.path.exists(input_dir):
            print(f'[ERROR] input_dir: {input_dir} is not a valid directory')
            return
        if not os.path.isdir(output_dir) or not os.path.exists(output_dir):
            print(f'[ERROR] output_dir: {output_dir} is not a valid directory')
            return
        print(f'[INFO] Input directory: {input_dir}')
        print(f'[INFO] Output directory: {output_dir}')
        filesi = get_all_file_dates(input_dir)
        print(f'[INFO] Number of files in input_dir: {len(filesi)}')
        fileso = get_all_file_dates(output_dir)
        print(f'[INFO] Number of files in output_dir: {len(fileso)}')

        input_dir_output = os.path.join(os.path.dirname(input_dir),f'{os.path.basename(input_dir)}_nocommon')
        output_dir_output = os.path.join(os.path.dirname(output_dir), f'{os.path.basename(input_dir)}_nocommon')
        if not os.path.exists(input_dir_output):
            os.mkdir(input_dir_output)
        if not os.path.exists(output_dir_output):
            os.mkdir(output_dir_output)

        for key in filesi:
            if not key in fileso:
                file_here = filesi[key]['file_path']
                file_output = os.path.join(input_dir_output,filesi[key]['file_name'])
                print(f'[INFO] Moving ifile {file_here} to {input_dir_output}')
                os.rename(file_here,file_output)

        for key in fileso:
            if not key in filesi:
                file_here = fileso[key]['file_path']
                file_output = os.path.join(output_dir_output,fileso[key]['file_name'])
                print(f'[INFO] Moving ofile {file_here} to {input_dir_output}')
                os.rename(file_here,file_output)

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
            if 'insitu_time' in dataset.variables:
                insitu_time_here = dt.utcfromtimestamp(dataset.variables['insitu_time'][0][0])
                insitu_time_here_str = insitu_time_here.strftime('%H%M%S')
                ref_here = f'{ref_here}{insitu_time_here_str}'
            key_here = f'{date_here_str}_{ref_here}'
            if key_here in file_dates.keys():
                print(f'KEY HERE: {key_here}')
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


def get_all_file_dates(input_dir):
    file_dates = {}
    for fname in os.listdir(input_dir):
        if fname.endswith('nc'):
            file_nc = os.path.join(input_dir, fname)
            dataset = Dataset(file_nc)
            date_here = dt.fromtimestamp(int(dataset.variables['satellite_time'][0]))
            date_here_str = date_here.strftime('%Y%m%d')
            ref_here = f'{dataset.satellite.upper()}{dataset.platform.upper()}{dataset.sensor.upper()}'
            if 'insitu_time' in dataset.variables:
                insitu_time_here = dt.utcfromtimestamp(dataset.variables['insitu_time'][0][0])
                insitu_time_here_str = insitu_time_here.strftime('%H%M%S')
                ref_here = f'{ref_here}{insitu_time_here_str}'
            key_here = f'{date_here_str}_{ref_here}'
            if key_here in file_dates.keys():
                print(f'[INFO] Repeated key: {key_here} for file: {fname}')
            file_dates[key_here] = {
                'file_path': file_nc,
                'file_name' : fname
            }
            dataset.close()
    return file_dates

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
                if 'insitu_time' in dataset.variables:
                    insitu_time_here = dt.utcfromtimestamp(dataset.variables['insitu_time'][0][0])
                    insitu_time_here_str = insitu_time_here.strftime('%H%M%S')
                    ref_here = f'{ref_here}{insitu_time_here_str}'
                key_here = f'{date_here_str}_{ref_here}'
                file_dates[key_here] = {
                    'file_path': file_nc,
                }

    return file_dates


if __name__ == '__main__':
    main()
