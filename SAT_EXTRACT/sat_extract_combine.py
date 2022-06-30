import argparse
import os.path
import numpy as np
import numpy.ma as ma
from datetime import datetime as dt
from netCDF4 import Dataset
import sys

code_home = os.path.abspath('../')
sys.path.append(code_home)
import COMMON.Class_Flags_OLCI as flag

parser = argparse.ArgumentParser(description="Create Match-up DataBase files (MDB) files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-i', "--input_extracts", help="Input directory with extracts combining A + B.")
parser.add_argument('-ac',"--ac_processor",help="AC Processors",choices=['STANDARD','POLYMER'])
parser.add_argument('-t', "--time", help="Time to be assigned to S3A extracts HH:MM")

args = parser.parse_args()


def main():
    print('Started')
    input_path = args.input_extracts
    if not os.path.exists(input_path) or not os.path.isdir(input_path):
        return
    output_path = os.path.join(os.path.dirname(input_path), 'extractsAB')
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    date_list = {}
    for fname in os.listdir(input_path):
        if not fname.endswith('.nc'):
            continue
        dthere = get_sat_time_from_fname(fname)
        dtkey = dthere.strftime('%Y%m%d')
        if not dtkey in date_list:
            date_list[dtkey] = {
                'files': [fname],
                'nfiles': 1
            }
        else:
            list_files = date_list[dtkey]['files']
            list_files.append(fname)
            date_list[dtkey]['files'] = list_files
            date_list[dtkey]['nfiles'] = date_list[dtkey]['nfiles'] + 1

    for d in date_list:
        name_orig = date_list[d]['files'][0]
        name_output, time_output = get_name_output(name_orig, args.time)
        forig = os.path.join(input_path, name_orig)
        foutput = os.path.join(output_path, name_output)
        dst = copy_nc(forig, foutput)
        dst.plafform = 'AB'
        dst.variables['satellite_time'][0] = float(time_output.timestamp())
        if date_list[d]['nfiles'] == 2:
            dataset2 = Dataset(os.path.join(input_path, date_list[d]['files'][1]))
            dst = combine_rrs_values(dst, dataset2, args.ac_processor)
        dst.close()


def combine_rrs_values(dataset, dother, ac):
    rrs1 = ma.array(dataset.variables['satellite_Rrs'][:])
    rrs2 = ma.array(dother.variables['satellite_Rrs'][:])
    nbands = rrs1.shape[1]
    if ac=='POLYMER':
        bitmask1 = np.array(dataset.variables['satellite_bitmask'][:])
        bitmask2 = np.array(dother.variables['satellite_bitmask'][:])
        ny = bitmask1.shape[1]
        nx = bitmask1.shape[2]
        array_mean_ref = np.zeros(bitmask1.shape)
        satellite_flag = dataset.variables['satellite_bitmask']
        flagging = flag.Class_Flags_Polymer(satellite_flag.flag_masks, satellite_flag.flag_meanings)
        flag1 = flagging.MaskGeneral(bitmask1)
        flag2 = flagging.MaskGeneral(bitmask2)

    if ac=='STANDARD':
        bitmask1 = np.array(dataset.variables['satellite_WQSF'][:])
        bitmask2 = np.array(dother.variables['satellite_WQSF'][:])
        ny = bitmask1.shape[1]
        nx = bitmask1.shape[2]
        array_mean_ref = np.zeros(bitmask1.shape)
        satellite_flag = dataset.variables['satellite_WQSF']
        flagging = flag.Class_Flags_OLCI(satellite_flag.flag_masks, satellite_flag.flag_meanings)
        flag_list = 'CLOUD,CLOUD_AMBIGUOUS,CLOUD_MARGIN,INVALID,COSMETIC,SATURATED,SUSPECT,HISOLZEN,HIGHGLINT,SNOW_ICE,AC_FAIL,WHITECAPS,RWNEG_O2,RWNEG_O3,RWNEG_O4,RWNEG_O5,RWNEG_O6,RWNEG_O7,RWNEG_O8'
        flag_list = flag_list.replace(" ", "")
        flag_list = str.split(flag_list, ',')
        flag1 = flagging.Mask(bitmask1, flag_list)
        flag2 = flagging.Mask(bitmask2, flag_list)


    nvalid1 = 0
    nvalid2 = 0
    for y in range(ny):
        for x in range(nx):
            escentral = False
            if y >= 11 and y <= 13 and x >= 11 and x <= 13:
                escentral = True
            if flag1[0][y][x] == 0 and flag2[0][y][x] == 0:  ##ambos validos
                array_mean_ref[0][y][x] = 3
                if escentral:
                    nvalid1 = nvalid1 + 1
                    nvalid2 = nvalid2 + 1
            elif flag1[0][y][x] == 0 and flag2[0][y][x] > 0:  # solo 1(S3A) valido
                array_mean_ref[0][y][x] = 1
                if escentral:
                    nvalid1 = nvalid1 + 1
            elif flag1[0][y][x] > 0 and flag2[0][y][x] == 0:  # solo 2 (S3B) valido
                array_mean_ref[0][y][x] = 2
                if escentral:
                    nvalid2 = nvalid2 + 1
            else:
                array_mean_ref[0][y][x] = 0

    if nvalid1 == 9 and nvalid2 == 9:  ##Hacemos la media
        rrsfin = np.copy(rrs1)
        bitmaskfin = np.copy(bitmask1)
        for b in range(nbands):
            for y in range(ny):
                for x in range(nx):
                    bitmaskfin[0][y][x] = 0
                    if array_mean_ref[0][y][x] == 0:
                        bitmaskfin[0][y][x] = bitmask1[0][y][x]
                    if array_mean_ref[0][y][x]==3:
                        rrsfin[0][b][y][x] = (rrs1[0][b][y][x]+rrs2[0][b][y][x])/2
                    elif array_mean_ref[0][y][x]==2:
                        rrsfin[0][b][y][x] = rrs2[0][b][y][x]
    elif nvalid1 >= nvalid2:
        rrsfin = np.copy(rrs1)
        bitmaskfin = np.copy(bitmask1)
    else:
        rrsfin = np.copy(rrs2)
        bitmaskfin = np.copy(bitmask2)

    dataset.variables['satellite_Rrs'] = rrsfin[:]
    dataset.variables['satellite_bitmask'] = bitmaskfin[:]

    return dataset


# new_time as HH:MM
def get_name_output(name_orig, new_time):
    name_new = name_orig
    if name_orig.startswith('S3A'):
        name_new = name_orig.replace('S3A', 'S3AB')
    if name_orig.startswith('S3B'):
        name_new = name_orig.replace('S3B', 'S3AB')

    sat_time = get_sat_time_from_fname(name_orig)
    str_orig = sat_time.strftime('%Y%m%dT%H%M%S')

    str_new = sat_time.strftime('%Y%m%dT')
    str_new = f'{str_new}{new_time}'
    time_new = dt.strptime(str_new, '%Y%m%dT%H:%M')
    str_new = time_new.strftime('%Y%m%dT%H%M%S')
    name_new = name_new.replace(str_orig, str_new)
    return name_new, time_new


def get_sat_time_from_fname(fname):
    val_list = fname.split('_')
    sat_time = None
    for v in val_list:
        try:
            sat_time = dt.strptime(v, '%Y%m%dT%H%M%S')
            break
        except ValueError:
            continue
    return sat_time


def copy_nc(ifile, ofile):
    with Dataset(ifile) as src:
        dst = Dataset(ofile, 'w', format='NETCDF4')
        if args.verbose:
            print(f'debug: copying {ifile} -----> {ofile}')

            # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)

        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))

        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            dst.createVariable(name, variable.datatype, variable.dimensions)

            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)

            dst[name][:] = src[name][:]

    return dst


if __name__ == '__main__':
    main()
