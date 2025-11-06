import argparse
import os.path
from datetime import timedelta
from datetime import datetime as dt
from datetime import timezone

import numpy as np
import pandas as pd
from netCDF4 import Dataset

parser = argparse.ArgumentParser(
    description="Obtaining information for running MDB_builder.")

parser.add_argument('-m', "--mode", help='Mode option', choices=["instrument_id", "match-up-info","extract_info","TEST"], required=True)
parser.add_argument('-i', "--input_path", help="Input path.")
parser.add_argument('-o', "--output", help="Output file.")
parser.add_argument('-sd', "--start_date", help="Start date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-ed', "--end_date", help="End date. Optional with --listdates (YYYY-mm-dd)")

args = parser.parse_args()


def main():
    print(f'[INFO] Started MDB_info!')
    if args.mode == 'TEST':
        make_test()
        return
    if args.mode == 'extract_info':
        if not check_required_params(['input_path', 'output']):
            return
        if not os.path.isdir(args.input_path):
            print(f'[ERROR] Input path {args.input_path} is not a valid directory')
            return
        if not args.output.endswith('.csv'):
            print(f'[ERROR] Output file {args.output} should be a CSV file')
            return
        try:
            os.makedirs(os.path.dirname(args.output),exist_ok=True)
        except Exception as ex:
            print(f'[ERROR] {os.path_dirname(args.output)} is not a valid directory and could not be created')

        get_extract_info(args.input_path,args.output)

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

    if args.mode == 'match-up-info':
        if not check_required_params(['input_path']):
            return
        if not os.path.isfile(args.input_path):
            print(f'[ERROR] {args.input_path} does not exist or is not a valid directory')
            return
        get_match_up_info(args.input_path)

def get_extract_info(input_path,output_file):
    print(f'[INFO] Starting getting info from extracts:')
    df = None
    for name in os.listdir(input_path):
        if not name.endswith('.nc'):
            continue
        file_nc = os.path.join(input_path,name)
        try:
            dataset = Dataset(file_nc,'r')
        except:
            continue
        print(f'[INFO] -->{name}')
        vals = dataset.__dict__
        dims = dataset.dimensions
        for dim in dims:
            vals[dim]=dims[dim].size
        if 'satellite_time' in dataset.variables:
            vals['sat_time'] = dt.fromtimestamp(float(dataset.variables['satellite_time'][0]),timezone.utc).strftime('%Y%m%dT%H%M%S')
        df_extract = pd.DataFrame({key:[vals[key]] for key in vals},index=[name])
        if df is None:
            df = df_extract.copy()
        else:
            df = pd.concat([df,df_extract])
        dataset.close()

    print(f'[INFO] Saving to: {output_file}')
    df.to_csv(output_file,sep=';')
    print(f'[INFO] Completed')


def get_match_up_info(input_file):
    dataset = Dataset(input_file)
    date_list = []
    date_min = None
    date_max = None
    for t in dataset.variables['satellite_time'][:]:
        tobj = dt.utcfromtimestamp(t)
        if date_min is None:
            date_min = tobj
        else:
            if tobj<date_min:
                date_min = tobj
        if date_max is None:
            date_max = tobj
        else:
            if tobj > date_max:
                date_max = tobj
        tstr = tobj.strftime('%Y%m%d')
        if tstr not in date_list:
            date_list.append(tstr)

    valid_mu = dataset.variables['mu_valid'][:]
    flag_sat = dataset.variables['flag_satellite'][:]
    nvalid = [0]*3
    ntotal = [0]*3
    for idx in range(valid_mu.shape[0]):
        flag_sat_index = flag_sat[idx]-1
        ntotal[2] = ntotal[2] + 1
        ntotal[flag_sat_index] = ntotal[flag_sat_index]+1
        if valid_mu[idx]==1:
            nvalid[2] = nvalid[2] + 1
            nvalid[flag_sat_index] = nvalid[flag_sat_index] + 1
    dataset.close()

    print(f'Start date: {date_min.strftime("%Y-%m-%d")}')
    print(f'End date: {date_max.strftime("%Y-%m-%d")}')
    print(f'Number of days: {len(date_list)}')
    print(nvalid)
    print(ntotal)
    porc = (np.array(nvalid)/np.array(ntotal))*100
    print(f'Sentinel3A: {ntotal[0]} {nvalid[0]} {porc[0]:.2f}')
    print(f'Sentinel3B: {ntotal[1]} {nvalid[1]} {porc[1]:.2f}')
    print(f'TOTAL: {ntotal[2]} {nvalid[2]} {porc[2]:.2f}')

def make_test():

    #file_csv = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/TEMPORAL/info_extracts_v3_1.csv'
    file_csv = ('/store3/IvanFarace/GaiaBlue/info_extracts_v3_1.csv')
    df = pd.read_csv(file_csv,sep=';')
    times = df['sat_time']
    dir_base = '/store3/OC/OCI'

    for time in times:
        time_obj = dt.strptime(time,'%Y%m%dT%H%M%S')
        name_dt = f'PACE_OCI.{time}.L2.OC_AOP.V3_0.nc'
        name_nrt = f'PACE_OCI.{time}.L2.OC_AOP.V3_0.NRT.nc'
        dir_out = os.path.join(dir_base,f'{time_obj.strftime("%Y")}',f'{time_obj.strftime("%j")}')
        os.makedirs(dir_out,exist_ok=True)
        file_out = os.path.join(dir_out,name_dt)
        if os.path.exists(file_out):
            print(f'[INFO] File {file_out} already exist. Skipping...')
            continue
        file_dt = os.path.join(dir_base,name_dt)
        file_nrt = os.path.join(dir_base, name_nrt)
        if os.path.exists(file_dt):
            print(f'[INFO] Moving DT file from {file_dt} to {file_out}')
            os.rename(file_dt,file_out)
        elif os.path.exists(file_nrt):
            print(f'[INFO ]Moving NRT file from {file_nrt} to {file_out}')
            os.rename(file_nrt, file_out)
        else:
            print(f'[ERROR] File for {time} is not available')



    ##CHECK FLAGS
    # file_level_1_2 = '/mnt/c/DATA_LUIS/ITALIAN_SITES_VALIDATION_PUBLICATION/OLCI/VEIT/MDBs/HYPERNETS_W_VEIT_L2A_REF_20220809T1640_20230330T1309_v1.2.nc'
    # file_level_2_0 = '/mnt/c/DATA_LUIS/ITALIAN_SITES_VALIDATION_PUBLICATION/OLCI/VEIT/MDBs/HYPERNETS_W_VEIT_L2A_REF_20250103T1500_20250103T1624_090_v2.0.nc'
    # file_level_2_1 = '/mnt/c/DATA_LUIS/ITALIAN_SITES_VALIDATION_PUBLICATION/OLCI/VEIT/MDBs/HYPERNETS_W_VEIT_L2A_REF_20240710T1645_20240831T1324_090_v2.0.nc'
    # dataset = Dataset(file_level_1_2)
    # flag_meanings = dataset.variables['quality_flag'].flag_meanings.split(' ')
    # flag_values = dataset.variables['quality_flag'].flag_masks.split(',')
    #
    # for m,v in zip(flag_meanings,flag_values):
    #     print(m,'->',v)
    # dataset.close()

    ##CHECK DATES
    # file_mdb = '/mnt/c/DATA_LUIS/ITALIAN_SITES_VALIDATION_PUBLICATION/OLCI/GAIT/MDBs/MDB_S3B_OLCI_WFR_STANDARD_20220609T000000_20241031T235959_HYPSTAR_GAIT.nc'
    # file_out = '/mnt/c/DATA_LUIS/ITALIAN_SITES_VALIDATION_PUBLICATION/OLCI/VEIT/MDBs/Temp.nc'
    # import pytz
    # date_1 = dt(2021, 4, 8).replace(tzinfo=pytz.utc).timestamp()
    # date_2 = dt(2023, 4, 25).replace(tzinfo=pytz.utc).timestamp()
    # date_3 = dt(2024, 1, 1).replace(tzinfo=pytz.utc).timestamp()
    # date_4 = dt(2024, 6, 4).replace(tzinfo=pytz.utc).timestamp()
    # date_5 = dt(2024, 10, 1).replace(tzinfo=pytz.utc).timestamp()
    # date_6 = dt(2025, 1, 31).replace(tzinfo=pytz.utc).timestamp()
    # dataset = Dataset(file_mdb)
    # qf = dataset.variables['insitu_quality_flag'][:]
    # site_flag = dataset.variables['insitu_site_flag'][:]
    #
    # vals_valid = [0, 1, 2, 3, 268435456, 268435457, 268435458, 268435459]
    # site_flag[qf.mask] = np.ma.masked
    # n_non_masked = np.ma.count(site_flag)
    # for val in vals_valid:
    #     site_flag[qf == val] = 1
    #     # print(np.sum(site_flag))
    # n_valid_qf = np.sum(site_flag)
    # n_qf = n_non_masked - n_valid_qf
    #
    # if 'insitu_epsilon' in dataset.variables:
    #     epsilon = dataset.variables['insitu_epsilon'][:]
    #     site_flag[epsilon < (-0.05)] = 0
    #     site_flag[epsilon > 0.05] = 0
    #     n_epsilon = (n_valid_qf - np.sum(site_flag))
    # else:
    #     n_epsilon = 0
    # print(f'[INFO] N. Spectra: {n_non_masked}')
    # print(f'[INFO] N. QF: {n_qf}')
    # print(f'[INFO] N. EPSILON: {n_epsilon}')
    # print(f'[INFO] N. VALID: {np.sum(site_flag)}')

    # sat_time = dataset.variables['satellite_time'][:]
    # n_sat_time = sat_time.shape[0]
    # flag_period = np.zeros((n_sat_time,))  ##1,2,4,8,16
    # periods = ['A','B','C','D','E']
    # periods_v = [1,2,4,8,16]
    # flag_hypstar_system = np.zeros((n_sat_time,))  ##1,2
    # systems = ['V1', 'V3']
    # systems_v = [1,2]
    # flag_hypstar_sn = np.zeros((n_sat_time,))  ##1,2,4
    # sn = ['120242','122304','122305']
    # sn_v = [1,2,4]
    # flag_hypernets_processing = np.zeros((n_sat_time,))  # 1,2
    # processings = ['v.1.2','v.2.0','v.2.1']
    # processings_v = [1,2,4]
    # for idx in range(n_sat_time):
    #     itime = sat_time[idx]
    #     if date_1 <= itime <= date_2:
    #         flag_period[idx] = 1
    #         flag_hypstar_system[idx] = 1
    #         flag_hypstar_sn[idx] = 1
    #         flag_hypernets_processing[idx] = 1
    #     elif date_2 <= itime <= date_3:
    #         flag_period[idx] = 2
    #         flag_hypstar_system[idx] = 2
    #         flag_hypstar_sn[idx] = 2
    #         flag_hypernets_processing[idx] = 2
    #     elif date_3 <= itime <= date_4:
    #         flag_period[idx] = 4
    #         flag_hypstar_system[idx] = 2
    #         flag_hypstar_sn[idx] = 2
    #         flag_hypernets_processing[idx] = 4
    #     elif date_4 <= itime <= date_5:
    #         flag_period[idx] = 8
    #         flag_hypstar_system[idx] = 2
    #         flag_hypstar_sn[idx] = 4
    #         flag_hypernets_processing[idx] = 4
    #     elif date_5 <= itime <= date_6:
    #         flag_period[idx] = 16
    #         flag_hypstar_system[idx] = 2
    #         flag_hypstar_sn[idx] = 4
    #         flag_hypernets_processing[idx] = 2
    #     iperiod = int(np.log2(flag_period[idx]))
    #     isystem = int(np.log2(flag_hypstar_system[idx]))
    #     isn = int(np.log2(flag_hypstar_sn[idx]))
    #     iproc = int(np.log2(flag_hypernets_processing[idx]))
    #     itime_obj = dt.utcfromtimestamp(itime)
    #     print(f'{itime_obj.strftime("%Y-%m-%d")} {periods[iperiod]} {systems[isystem]} {sn[isn]} {processings[iproc]}')

    ##WRITTING OUPPUT DATASET
    # ncout = Dataset(file_out, 'w', format='NETCDF4')
    # # copy global attributes all at once via dictionary
    # ncout.setncatts(dataset.__dict__)
    # # copy dimensions
    # for name, dimension in dataset.dimensions.items():
    #     ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
    #
    # # copy variables
    # for name, variable in dataset.variables.items():
    #     fill_value = None
    #     if '_FillValue' in list(variable.ncattrs()):
    #         fill_value = variable._FillValue
    #     ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True, complevel=6)
    #     ncout[name].setncatts(dataset[name].__dict__)
    #     if name=='insitu_site_flag':
    #         ncout[name][:] = site_flag[:]
    #     else:
    #         ncout[name][:] = dataset[name][:]

    ##new flag variables
    # var = ncout.createVariable('flag_period','i2',('satellite_id',),fill_value=None,zlib=True,complevel=6)
    # var[:] = flag_period[:]
    # var.flag_meanings = " ".join(periods)
    # var.flag_values = periods_v
    #
    # var = ncout.createVariable('flag_hypstar_system', 'i2', ('satellite_id',), fill_value=None, zlib=True, complevel=6)
    # var[:] = flag_hypstar_system[:]
    # var.flag_meanings = " ".join(systems)
    # var.flag_values = systems_v
    #
    # var = ncout.createVariable('flag_hypstar_sn', 'i2', ('satellite_id',), fill_value=None, zlib=True, complevel=6)
    # var[:] = flag_hypstar_sn
    # var.flag_meanings = " ".join(sn)
    # var.flag_values = sn_v
    #
    # var = ncout.createVariable('flag_hypernets_processing', 'i2', ('satellite_id',), fill_value=None, zlib=True, complevel=6)
    # var[:] = flag_hypernets_processing[:]
    # var.flag_meanings = " ".join(processings)
    # var.flag_values = processings_v

    # ncout.close()
    #
    #
    # dataset.close()
    #
    # os.rename(file_out,file_mdb)


def get_info_instrument_id(input_path, start_date, end_date):
    date_ref = start_date
    ids_list = {}
    while date_ref <= end_date:
        # print(f'[INFO] Checking date: {date_ref.strftime("%Y-%m-%d")}')
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
                        time_here = float(dataset.variables['acquisition_time'][0])
                        if time_here == 0:
                            continue
                        if np.isnan(time_here):
                            continue
                        time_here_obj = dt.utcfromtimestamp(time_here)
                        id_here = str(dataset.instrument_id)
                        if id_here not in ids_list.keys():
                            nwl = dataset.variables['wavelength'].shape[0]
                            ids_list[id_here] = {
                                'nwl': nwl,
                                'start_time': time_here_obj,
                                'end_time': time_here_obj
                            }
                        else:
                            if time_here_obj < ids_list[id_here]['start_time']:
                                ids_list[id_here]['start_time'] = time_here_obj
                            if time_here_obj > ids_list[id_here]['end_time']:
                                ids_list[id_here]['end_time'] = time_here_obj

                    dataset.close()
        date_ref = date_ref + timedelta(hours=24)
    for id_here in ids_list.keys():
        stime = ids_list[id_here]['start_time'].strftime('%Y-%m-%d %H:%M:%S')
        etime = ids_list[id_here]['end_time'].strftime('%Y-%m-%d %H:%M:%S')
        print(
            f'[INFO] Instrument id: {id_here} Number of wavelengths: {ids_list[id_here]["nwl"]} Start time: {stime} End time: {etime} ')


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
