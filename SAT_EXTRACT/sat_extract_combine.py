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
parser.add_argument('-ac', "--ac_processor", help="AC Processors", choices=['STANDARD', 'POLYMER'])
parser.add_argument('-t', "--time", help="Time to be assigned to S3A extracts HH:MM")
parser.add_argument('-c', "--config_file", help="Config File.")
args = parser.parse_args()


def main():
    print('Started')
    input_path = args.input_extracts
    if not os.path.exists(input_path) or not os.path.isdir(input_path):
        return
    output_path = os.path.join(os.path.dirname(input_path), 'extractsAB')
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    config_file = args.config_file
    if not os.path.exists(config_file):
        return
    from sat_extract_quality import QC_EXTRACT


    new_extracts = get_new_extract_list(input_path,output_path)
    from sat_extract_quality import QC_EXTRACT
    qc = QC_EXTRACT(config_file)
    create_avg_extracts(input_path,new_extracts,qc)

    ##DEPRECATED
    # date_list = {}
    # nfiles = 0
    # for fname in os.listdir(input_path):
    #     if not fname.endswith('.nc'):
    #         continue
    #     nfiles = nfiles + 1
    #     dtsathere = get_sat_time_from_fname(fname)
    #     dtinsituhere = get_insitu_time(os.path.join(input_path, fname))
    #     dtsathere_s = dtsathere.strftime('%Y%m%d')
    #     dtinsituhere_s = dtinsituhere.strftime('%H:%M:%S')
    #     dtkey = f'{dtsathere_s}T{dtinsituhere_s}'
    #     if not dtkey in date_list:
    #         date_list[dtkey] = {
    #             'files': [fname],
    #             'nfiles': 1
    #         }
    #     else:
    #         list_files = date_list[dtkey]['files']
    #         list_files.append(fname)
    #         date_list[dtkey]['files'] = list_files
    #         date_list[dtkey]['nfiles'] = date_list[dtkey]['nfiles'] + 1
    #
    # ndates = 0
    # nab = 0
    # for d in date_list:
    #     ndates = ndates + 1
    #     name_orig = date_list[d]['files'][0]
    #     name_output, time_output = get_name_output(name_orig, args.time)
    #     forig = os.path.join(input_path, name_orig)
    #     foutput = os.path.join(output_path, name_output)
    #     dst = copy_nc(forig, foutput)
    #     dst.plafform = 'AB'
    #     dst.variables['satellite_time'][0] = float(time_output.timestamp())
    #     if date_list[d]['nfiles'] == 2:
    #         print(f'A+B for: {d}')
    #         nab = nab +1
    #         dataset2 = Dataset(os.path.join(input_path, date_list[d]['files'][1]))
    #         dst = combine_rrs_values(dst, dataset2, args.ac_processor)
    #     dst.close()
    #
    # print('NDates: ', ndates, 'AB:', nab)

def get_new_extract_list(input_path,output_path):
    from datetime import datetime as dt
    from datetime import timezone as tz
    new_extracts = {}
    for name in os.listdir(input_path):
        if name.startswith('S3A_'):
            date_here = get_sat_time_from_fname(name)
            key = date_here.strftime('%Y%m%d')
            if not key in new_extracts:
                new_extracts[key] = {
                    'S3A': name,
                    'S3ADate':date_here,
                    'S3B': None,
                    'S3BDate': None,
                    'file_out': None,
                    'date_out': None
                }
            else:
                new_extracts[key]['S3A'] = name
                new_extracts[key]['S3ADate'] = date_here
        if name.startswith('S3B_'):
            date_here = get_sat_time_from_fname(name)
            key = date_here.strftime('%Y%m%d')
            if not key in new_extracts:
                new_extracts[key] = {
                    'S3A': None,
                    'S3ADate': None,
                    'S3B': name,
                    'S3BDate': date_here,
                    'file_out': None,
                    'date_out': None
                }
            else:
                new_extracts[key]['S3B'] = name
                new_extracts[key]['S3BDate'] = date_here

    for key in new_extracts:
        if new_extracts[key]['S3A'] is not None and new_extracts[key]['S3B'] is not None:

            tsavg = (new_extracts[key]['S3ADate'].timestamp()+new_extracts[key]['S3BDate'].timestamp())/2
            date_out = dt.fromtimestamp(tsavg).replace(tzinfo=tz.utc)
            date_out_str = date_out.strftime('%Y%m%dT%H%M%S')
            date_creation_str = dt.utcnow().strftime('%Y%m%dT%H%M%S')
            name_list = new_extracts[key]['S3A'].split('_')
            file_out_name = f'S3AB_{name_list[1]}_{name_list[2]}_{name_list[3]}____{date_out_str}_{date_out_str}_{date_creation_str}'
            for idx in range(10,len(name_list)):
                file_out_name = f'{file_out_name}_{name_list[idx]}'
            file_out = os.path.join(output_path,file_out_name)
            new_extracts[key]['file_out'] = file_out
            new_extracts[key]['date_out'] = date_out
        else:
            if new_extracts[key]['S3A'] is not None:
                file_out_name = new_extracts[key]['S3A'].replace('S3A','S3AB')
                file_out = os.path.join(output_path, file_out_name)
                new_extracts[key]['file_out'] = file_out
                new_extracts[key]['date_out'] = new_extracts[key]['S3ADate'].replace(tzinfo=tz.utc)
            if new_extracts[key]['S3B'] is not None:
                file_out_name = new_extracts[key]['S3B'].replace('S3B','S3AB')
                file_out = os.path.join(output_path, file_out_name)
                new_extracts[key]['file_out'] = file_out
                new_extracts[key]['date_out'] = new_extracts[key]['S3BDate'].replace(tzinfo=tz.utc)

    return new_extracts

def create_avg_extracts(input_path,new_extracts,qc):
    for key in new_extracts:
        file_out = new_extracts[key]['file_out']
        date_out = new_extracts[key]['date_out']
        if new_extracts[key]['S3A'] is not None and new_extracts[key]['S3B'] is not None:
            filea = os.path.join(input_path,new_extracts[key]['S3A'])
            valida = qc.compute_array_valid(filea)
            na = np.sum(valida)
            fileb = os.path.join(input_path, new_extracts[key]['S3B'])
            validb = qc.compute_array_valid(fileb)
            nb = np.sum(validb)
            if na>0 and nb>0:
                create_avg_extract_impl(file_out,filea,valida,fileb,validb,date_out)
            elif na>=0 and nb==0:
                create_single_extract_impl(file_out,filea,valida,date_out)
            elif na==0 and nb>0:
                create_single_extract_impl(file_out, fileb, validb, date_out)
        elif new_extracts[key]['S3A'] is not None and new_extracts[key]['S3B'] is  None:
            filea = os.path.join(input_path, new_extracts[key]['S3A'])
            valida = qc.compute_array_valid(filea)
            create_single_extract_impl(file_out, filea, valida, date_out)
        elif new_extracts[key]['S3A'] is None and new_extracts[key]['S3B'] is not None:
            fileb = os.path.join(input_path, new_extracts[key]['S3B'])
            validb = qc.compute_array_valid(fileb)
            create_single_extract_impl(file_out, fileb, validb, date_out)

def create_single_extract_impl(file_out,file,valid,date_out):
    print(f'[INFO] Creating single file out: {file_out}')
    dst = copy_nc(file, file_out)
    dst.platform = 'AB'
    dst.variables['satellite_time'][0] = float(date_out.timestamp())
    var = dst.createVariable('satellite_combine_flag', 'f4', ('satellite_id', 'rows', 'columns'), fill_value=-999,
                             zlib=True, complevel=6)
    var_array = np.zeros(valid.shape)
    var_array[valid == 0] = 1
    var[0, :, :] = var_array[:, :]
    var.description = 'Flag created after combine satellite extracts from A and B platforms'
    var.flag_masks = [1]
    var.flag_meanings = 'INVALID'

    dst.close()


def create_avg_extract_impl(file_out,filea,valida,fileb,validb,date_out):
    print(f'[INFO] Creating file out: {file_out}')
    dst = copy_nc(filea, file_out)
    dst.platform = 'AB'
    dst.variables['satellite_time'][0] = float(date_out.timestamp())

    # date_out_check = dt.utcfromtimestamp(float(dst.variables['satellite_time'][0]))
    # print(date_out,date_out_check)
    dataseta = Dataset(filea)
    rrsa = np.array(dataseta.variables['satellite_Rrs'][:])
    filla = dataseta.variables['satellite_Rrs']._FillValue
    dataseta.close()

    datasetb = Dataset(fileb)
    rrsb = np.array(datasetb.variables['satellite_Rrs'][:])
    fillb = dataseta.variables['satellite_Rrs']._FillValue
    datasetb.close()

    rrs = np.zeros(rrsa.shape)
    n = valida + validb
    nbands = rrs.shape[1]

    for iband in range(nbands):
        ahere = rrsa[0,iband,:,:]
        ahere[valida==0] = 0
        ahere[ahere==filla] = 0

        bhere = rrsb[0, iband, :, :]
        bhere[bhere==fillb] = 0
        bhere[validb == 0] = 0

        nhere = n.copy()
        nhere[ahere == filla] = nhere[ahere == filla] - 1
        nhere[ahere == fillb] = nhere[ahere == fillb] - 1
        nhere[nhere<0] = 0

        rrs[0,iband,:,:] = np.divide(np.add(ahere,bhere),nhere)

    rrs[np.isinf(rrs)] = -999.0
    rrs[np.isnan(rrs)] = -999.0

    dst.variables['satellite_Rrs'] = rrs[:]

    var = dst.createVariable('satellite_combine_flag', 'f4', ('satellite_id','rows','columns'),fill_value=-999, zlib=True, complevel=6)
    var_array = np.zeros(n.shape)
    var_array[n==0] = 1

    var[0, :, :] = var_array[:,:]
    var.description = 'Flag created after combine satellite extracts from A and B platforms'
    var.flag_masks = [1]
    var.flag_meanings = 'INVALID'

    dst.close()

def combine_rrs_values(dataset, dother, ac):
    rrs1 = ma.array(dataset.variables['satellite_Rrs'][:])
    rrs2 = ma.array(dother.variables['satellite_Rrs'][:])
    nbands = rrs1.shape[1]
    if ac == 'POLYMER':
        bitmask1 = np.array(dataset.variables['satellite_bitmask'][:])
        bitmask2 = np.array(dother.variables['satellite_bitmask'][:])
        ny = bitmask1.shape[1]
        nx = bitmask1.shape[2]
        array_mean_ref = np.zeros(bitmask1.shape)
        satellite_flag = dataset.variables['satellite_bitmask']
        flagging = flag.Class_Flags_Polymer(satellite_flag.flag_masks, satellite_flag.flag_meanings)
        flag1 = flagging.MaskGeneral(bitmask1)
        flag2 = flagging.MaskGeneral(bitmask2)

    if ac == 'STANDARD':
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
                    if array_mean_ref[0][y][x] == 3:
                        rrsfin[0][b][y][x] = (rrs1[0][b][y][x] + rrs2[0][b][y][x]) / 2
                    elif array_mean_ref[0][y][x] == 2:
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


def get_insitu_time(fname):
    nc = Dataset(fname)
    sec = float(nc.variables['insitu_time'][0])
    insitudt = dt.fromtimestamp(sec)
    nc.close()
    return insitudt


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
        # if args.verbose:
        #     print(f'debug: copying {ifile} -----> {ofile}')

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
