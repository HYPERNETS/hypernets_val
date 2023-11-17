import argparse, configparser
from MDB_builder_options import MDBBuilderOptions
import os
from SATEXTRACTS_list import SAT_EXTRACTS_LIST
import numpy as np
from netCDF4 import Dataset
from datetime import datetime as dt
from datetime import timezone
import pytz

parser = argparse.ArgumentParser(
    description="Create Match-up DataBase files (MDB) files from satellite extracts and AERONET files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-site', "--sitename", help="Site name. Only with --listdates")
parser.add_argument('-c', "--config_file", help="Config File.", required=True)
parser.add_argument('-o', "--output", help="Output file. Only with --listdates")
parser.add_argument('-edir', "--sat_extract_dir", help="Input sat. extract dir. Only with --listdates")
parser.add_argument('-ld', "--listdates",
                    help="Option to obtain a date list for a specific AERONET site (-site option).",
                    action="store_true")
parser.add_argument('-nd', "--nodelfiles", help="Do not delete temp files.", action="store_true")

args = parser.parse_args()

import MDB_builder_options
import sys

code_home = os.path.dirname(os.path.dirname(MDB_builder_options.__file__))

# # code_home = os.path.abspath('../')
sys.path.append(code_home)
code_aeronet = os.path.join(os.path.dirname(code_home), 'aeronet')
sys.path.append(code_aeronet)

# import BRDF.brdf_olci as brdf
from COMMON import common_functions as cfs


def get_time_from_extract_file(extract_path):
    dataset = Dataset(extract_path)
    satellite_time_var = None
    if 'satellite_time' in dataset.variables:
        satellite_time_var = dt.utcfromtimestamp(float(dataset.variables['satellite_time'][0]))
    dataset.close()
    satellite_time_file = None
    try:
        time_str = extract_path.split('/')[-1].split('_')[7]
        satellite_time_file = dt.strptime(time_str, '%Y%m%dT%H%M%S')
    except:
        try:
            time_str = extract_path.split('/')[-1][1:8]
            satellite_time_file = dt.strptime(time_str, '%Y%j')
        except:
            pass

    if satellite_time_file is not None and satellite_time_var is not None:
        if satellite_time_file == satellite_time_var:
            return satellite_time_var
        else:
            print(f'[WARNING] Date time in name file is different from time in satellite_time variable')
            print(f'[WARNING] Satellite time set to: {satellite_time_file}')
            return satellite_time_file
    if satellite_time_file is not None and satellite_time_var is None:
        return satellite_time_file

    return satellite_time_var


def add_insitu_aeronet(extract_info, ofile, areader, mo_options, time_list, ins_time_ini,ins_time_end):
    extract_path = extract_info['path']
    satellite_datetime =  extract_info['time']
    # satellite_datetime = get_time_from_extract_file(extract_path)
    #
    # if satellite_datetime.hour == 0 and satellite_datetime.minute == 0:
    #     satellite_datetime = dt.strptime(extract_info['time'],'%Y%m%dT%H%M%S')
    #     if args.verbose:
    #         print(f'[WARNING] Satellite time set to default hour: {satellite_datetime}')

    if satellite_datetime is None:
        print(f'[ERROR] Extract satellite datetime could not be retrieved.')
        return False
    else:
        print(f'[INFO] Satellite time -> {satellite_datetime}')

    if satellite_datetime<ins_time_ini or satellite_datetime>ins_time_end:
        print(f'[WARNING] No in situ spectra were found for date: {satellite_datetime}')
        return False
    satellite_date = satellite_datetime.strftime('%Y-%m-%d')
    check_date = areader.prepare_data_fordate(satellite_date)
    if not check_date:
        print(f'[WARNING] No in situ spectra were found for date: {satellite_date}')
        return False

    nspectra = len(time_list)
    nspectra_within_timewindow = 0
    spectra_within_timewindow = np.zeros(nspectra, dtype=bool)
    time_dif_seconds = np.zeros(nspectra, dtype=float)
    time_window_seconds = mo_options.insitu_options['time_window']  # seconds

    for i in range(nspectra):
        time_dif_seconds[i] = float(abs((time_list[i] - satellite_datetime).total_seconds()))
        #print(satellite_datetime, time_list[i], time_dif_seconds[i], '============================================')
        if time_dif_seconds[i] < time_window_seconds:
            spectra_within_timewindow[i] = True
            nspectra_within_timewindow = nspectra_within_timewindow + 1

    if args.verbose:
        print(f'[INFO] # In situ spectra within the time window: {nspectra_within_timewindow}')

    if nspectra_within_timewindow == 0:
        print(f'[WARNING] No in situ spectra were found for date: {satellite_date} within the specified time window')
        return False

    nominal_wavelengths = np.array(mo_options.get_wllist())
    subset_wl = True
    if mo_options.get_wllist() is None:
        nominal_wavelengths = areader.dataset['Nominal_Wavelenghts'][:]
        subset_wl = False


    rrs = areader.extract_rrs(False)
    exactwl = areader.extract_spectral_data('Exact_Wavelengths', False)

    ##WHEN A SUBSET IS DEFINED, BANDS ARE BAND SHIFTED
    if subset_wl:
        from BSC_QAA import bsc_qaa_EUMETSAT as qaa
        maxwldiff = mo_options.get_maxwl_diff()
        exactwl_h = exactwl[0][:]

        nominal_wavelengths = update_nominal_wl_with_sat_wl(nominal_wavelengths, extract_path, maxwldiff)
        indiceswl_valid, nominal_wavelengths_exact = get_wl_valid_indices(nominal_wavelengths, exactwl_h, maxwldiff)


        rrsh = rrs[0][:]
        indices_good = np.ma.where(rrsh.mask == False)

        ns = rrs.shape[0]
        nwl = len(nominal_wavelengths)
        rrsnew = np.ma.zeros((ns, nwl))
        exactwlnew = np.ma.zeros((ns, nwl))

        for idx in range(ns):
            if len(indices_good[0]) == 1 and indices_good[0][0] == 0:
                rrsnew[idx, :] = qaa.bsc_qaa(rrs[idx][:], exactwl[idx][:], nominal_wavelengths)
            else:
                rrsnew[idx, :] = qaa.bsc_qaa(rrs[idx][indices_good], exactwl[idx][indices_good], nominal_wavelengths)
            rrsnew[idx, indiceswl_valid < 0] = np.ma.masked
            exactwlnew[idx, :] = nominal_wavelengths_exact[:]

        rrs = rrsnew
        exactwl = exactwlnew

    # to append to nc file
    if subset_wl:
        new_MDB = copy_nc(extract_path, ofile, satellite_datetime, nominal_wavelengths, maxwldiff)
    else:
        new_MDB = copy_nc(extract_path, ofile, satellite_datetime, None, None)

    # add time window diff
    new_MDB.time_diff = f'{time_window_seconds}'  # in seconds

    if subset_wl:
        n_insitu_bands = len(nominal_wavelengths)
    else:
        n_insitu_bands = areader.nwl

    # create in situ dimensions
    new_MDB.createDimension('insitu_id', 50)
    new_MDB.createDimension('insitu_original_bands', n_insitu_bands)

    # create variables
    insitu_time = new_MDB.createVariable('insitu_time', 'f8', ('satellite_id', 'insitu_id',), zlib=True, complevel=6)
    insitu_time.units = "Seconds since 1970-1-1"
    insitu_time.description = 'In situ time in ISO 8601 format (UTC).'

    insitu_original_bands = new_MDB.createVariable('insitu_original_bands', 'f4', ('insitu_original_bands'),
                                                   fill_value=-999, zlib=True, complevel=6)
    insitu_original_bands.description = 'In situ bands in nm'

    insitu_Rrs = new_MDB.createVariable('insitu_Rrs', 'f4', ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                        fill_value=-999, zlib=True, complevel=6)
    insitu_Rrs.description = 'In situ Rrs'

    insitu_exact_wavelenghts = new_MDB.createVariable('insitu_exact_wavelenghts', 'f4',
                                                      ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                                      fill_value=-999, zlib=True, complevel=6)
    insitu_exact_wavelenghts.description = 'In situ bands in nm'

    time_difference = new_MDB.createVariable('time_difference', 'f4', ('satellite_id', 'insitu_id'), fill_value=-999,
                                             zlib=True, complevel=6)
    time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
    time_difference.units = "seconds"

    ##SETTING DATA
    insitu_original_bands[:] = nominal_wavelengths[:]

    insitu_idx = 0
    for i in range(areader.row_ini, areader.row_fin + 1):
        if not spectra_within_timewindow[i]:
            continue

        insitu_time[0, insitu_idx] = float(time_list[i].replace(tzinfo=timezone.utc).timestamp())

        time_difference[0, insitu_idx] = time_dif_seconds[i]


        idx = i - areader.row_ini
        rrs_here = np.array(rrs[idx][:])

        rrs_here[rrs[idx][:].mask] = -999
        insitu_Rrs[0, :, insitu_idx] = rrs_here

        wl_here = np.array(exactwl[idx][:])
        wl_here[exactwl[idx][:].mask] = -999
        insitu_exact_wavelenghts[0, :, insitu_idx] = wl_here

        insitu_idx = insitu_idx + 1

    new_MDB.close()

    if insitu_idx == 0:
        if os.path.exists(ofile):
            os.remove(ofile)
            if args.verbose:
                print('Not in situ measurements within the time window. MDB file deleted!')
        return False
    else:
        return True


def get_wl_valid_indices(nominal_wl, exact_wl, maxdiff):
    indices = np.zeros(nominal_wl.shape).astype(int)
    nominal_wl_exact = np.ma.zeros(nominal_wl.shape)
    indices[:] = -1
    for idx in range(len(nominal_wl)):
        wlhere = nominal_wl[idx]
        iexact = np.argmin(np.abs(wlhere - exact_wl))
        if np.abs(exact_wl[iexact] - wlhere) < maxdiff:
            indices[idx] = iexact
            nominal_wl_exact[idx] = exact_wl[iexact]
        else:
            nominal_wl_exact[idx] = np.ma.masked

    return indices, nominal_wl_exact


def update_nominal_wl_with_sat_wl(nominal_wl, extract_path, maxwldiff):
    dataset = Dataset(extract_path)
    satellite_wl = np.array(dataset.variables['satellite_bands'][:])
    dataset.close()

    for idx in range(len(nominal_wl)):
        wl = nominal_wl[idx]
        diff = np.abs(wl - satellite_wl)
        iwl = np.argmin(diff)
        # print(wl,iwl,satellite_wl[iwl],maxwldiff)
        if diff[iwl] < maxwldiff:
            nominal_wl[idx] = satellite_wl[iwl]

    return nominal_wl


def copy_nc(ifile, ofile, satellite_time, nominalwl, maxwldif):
    with Dataset(ifile) as src:
        dst = Dataset(ofile, 'w', format='NETCDF4')

        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)

        # copy dimensions
        for name, dimension in src.dimensions.items():
            if name == 'satellite_bands' and nominalwl is not None:
                dst.createDimension(name, len(nominalwl))
            else:
                dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            dst.createVariable(name, variable.datatype, variable.dimensions)

            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)
            if nominalwl is None:
                dst[name][:] = src[name][:]
            else:
                if name == 'satellite_bands':
                    dst[name][:] = nominalwl[:]
                elif name == 'satellite_Rrs':
                    dim_orig = src.variables[name].shape
                    array = np.ma.zeros((dim_orig[0], len(nominalwl), dim_orig[2], dim_orig[3]))
                    array[:] = np.ma.masked
                    wl_sat_array = np.array(src['satellite_bands'][:])
                    indices_valid, wl_sat_valid = get_wl_valid_indices(nominalwl, wl_sat_array, maxwldif)
                    for idx in range(len(nominalwl)):
                        index = indices_valid[idx]
                        if index >= 0:
                            array[0, idx, :, :] = src[name][0, index, :, :]
                    dst[name][:] = array[:]
                elif name == 'satellite_time':
                    satellite_time_real = float(satellite_time.replace(tzinfo=timezone.utc).timestamp())
                    satellite_time_prev = float(src[name][0])
                    if satellite_time_real != satellite_time_prev:
                        print(f'[WARNING] Setting time to: {satellite_time}')
                        dst[name][0] = satellite_time_real
                    else:
                        dst[name][:] = src[name][:]

                else:
                    dst[name][:] = src[name][:]
    return dst


def main():
    print('[INFO] Creating MDB files!')
    # Option to a list dates
    if args.sitename and args.output and args.listdates:
        # no implemented
        return

    if os.path.isfile(args.config_file):
        options = configparser.ConfigParser()
        options.read(args.config_file)
        mo = MDBBuilderOptions(options, args.verbose)
    else:
        print('[ERROR] Configuration file does not exist')
        return

    valid_options = mo.set_compulsory_options(cfs)
    if not valid_options:
        return

    if not mo.insitu_type == 'AERONET':
        return
    ##sat extracts options
    mo.get_param_sat_extracts()

    ##in situ options
    mo.get_insitu_options()

    aeronet_file = mo.get_insitu_file(None)
    if aeronet_file is None:
        print(f'[ERROR] Aeronet NC file is not avaialable for the site in path: {mo.insitu_path_source}')
    from base.anet_nc_reader import AERONETReader
    areader = AERONETReader(aeronet_file)
    if args.verbose:
        print(f'[INFO] Path to AERONET NC file: {aeronet_file}')

    ##dates
    if args.verbose:
        print(f'[INFO] Checking available dates--------------------------------------------------------------START')
    mo.get_dates()
    sdate, edate = areader.get_time_ini_fin()
    if sdate>mo.start_date:
        mo.start_date = sdate.replace(hour=0,minute=0,second=0,microsecond=0)

    if edate<mo.end_date:
        from datetime import timedelta
        mo.end_date = (edate+timedelta(days=1)).replace(hour=0,minute=0,second=0,microsecond=0)

    if args.verbose:
        print(f'[INFO] Start date for MDB_builder:{mo.start_date}')
        print(f'[INFO] End date for MDB_builder: {mo.end_date}')
        print(f'[INFO] Checking available dates--------------------------------------------------------------STOP')

    ins_sensor = 'AERONET'
    time_list = areader.extract_time_list()
    time_list_range = []
    for itime in time_list:
        if mo.start_date <= itime <= mo.end_date:
            time_list_range.append(itime)
    if len(time_list_range)==0:
        print(f'[WARNING] No insitu files were found in the given range')
        return
    ins_time_ini = time_list_range[0].replace(hour=0, minute=0, second=0, microsecond=0)
    from datetime import timedelta
    ins_time_end = (time_list_range[-1] + timedelta(hours=24)).replace(hour=0, minute=0, second=0, microsecond=0)

    ##retrieving sat extract list
    if args.verbose:
        print(f'[INFO] Obtaining extract list----------------------------------------------------------------START')
    slist = SAT_EXTRACTS_LIST(mo, args.verbose)
    extract_list = slist.get_list_as_dict()
    if args.verbose:
        print(f'[INFO] Obtaining extract list----------------------------------------------------------------STOP')

    mdb_extract_files = []

    # for extract in extract_list:
    date_ref = mo.start_date
    while date_ref<=mo.end_date:
        date_ref_str = date_ref.strftime('%Y%m%d')
        if date_ref_str in extract_list.keys():
            extract = extract_list[date_ref_str]
            name = extract['path'].split('/')[-1]
            if args.verbose:
                print(f'[INFO] Working with extract file: {name} *******************')
            ofile = mo.get_mdb_extract_path(extract['path'], ins_sensor)
            if args.verbose:
                print(f'[INFO] MDB extract file: {ofile}')
            if os.path.exists(ofile):
                mdb_extract_files.append(ofile)
                print(f'[WARNING] MDB extract file already exits. Skipping...')
                continue
            b = add_insitu_aeronet(extract_list[date_ref_str], ofile, areader, mo, time_list,ins_time_ini,ins_time_end)
            if b:
                mdb_extract_files.append(ofile)
    nextract_files = len(mdb_extract_files)
    if args.verbose:
        print(f'[INFO] {nextract_files} were created/added')
        print(f'[INFO] Generating MDB extract files----------------------------------------------------------STOP')
    if nextract_files == 0:
        print(f'[INFO] Completed. No MDB file was created')
        return

    if args.verbose:
        print(f'[INFO] Generating final MDB file-------------------------------------------------------------START')
    path_mdb = mo.get_mdb_path()
    concatenate_nc_impl(mdb_extract_files, mo.path_out, path_mdb)
    if args.verbose:
        print(f'[INFO] Generating final MDB file-------------------------------------------------------------STOP')


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

        if not args.nodelfiles:
            [os.remove(f) for f in list_files]

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

        if not args.nodelfiles:
            [os.remove(f) for f in list_files[:-1]]
    if args.verbose:
        print(f'[INFO] Concatenated MDB file created: {ncout_file}')


# %%
if __name__ == '__main__':
    main()
