import datetime

import netCDF4

from MDBFile import MDBFile
from MDBInSituFile import MDBInSituFile
from MDBPlot import MDBPlot
from MDBFileList import MDBFileList
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
from PlotSpectra import PlotSpectra
import numpy as np
import math
from datetime import datetime as dt
from datetime import timedelta


class MDB_READER():
    def __init__(self, path_mdb, start_mdb):
        if path_mdb is not None:
            self.path_mdb = path_mdb
            if start_mdb:
                self.mfile = MDBFile(path_mdb)

    def compare_different_qcinsitu(self, qc_insitu_list, file_out):

        col_names = ['Param', 'Wavelength']
        for qc in qc_insitu_list:
            col_names.append(qc['name'])
        print(col_names)
        df_summary = pd.DataFrame(columns=col_names)
        started = False
        for qc in qc_insitu_list:
            mfile_here = MDBFile(self.path_mdb)
            mfile_here.qc_sat.set_eumetsat_defaults(3)
            mfile_here.qc_insitu.time_max = qc['time_max']
            mfile_here.qc_insitu.name = qc['name']
            # mfile_here.qc_insitu.apply_band_shift = False
            # mfile_here.qc_insitu.check_indices_by_mu = False
            # print(mfile_here.qc_insitu.name, mfile_here.qc_insitu.time_max)
            n = mfile_here.prepare_df_validation()
            if n > 0:
                mplot = MDBPlot(None, mfile_here.df_validation_valid)
                mplot.compute_all_statistics(mfile_here.wlref)
                if not started:
                    started = True
                    col_names = list(mplot.df_valid_stats.columns)[2:]
                    param_names = mplot.df_valid_stats['Param']

                    param_series = pd.Series(param_names)
                    nparam = param_series.size
                    wl_series = pd.Series('All', range(nparam))
                    val_series = mplot.df_valid_stats['All']
                    for col in col_names:
                        param_series = pd.concat([param_series, pd.Series(param_names)], ignore_index=True)
                        wl_series = pd.concat([wl_series, pd.Series(col, range(nparam))], ignore_index=True)
                        val_series = pd.concat([val_series, mplot.df_valid_stats[col]], ignore_index=True)
                    df_summary['Param'] = param_series
                    df_summary['Wavelength'] = wl_series
                    df_summary[qc['name']] = val_series
                else:
                    col_names = list(mplot.df_valid_stats.columns)[2:]
                    val_series = mplot.df_valid_stats['All']
                    for col in col_names:
                        val_series = pd.concat([val_series, mplot.df_valid_stats[col]], ignore_index=True)
                    df_summary[qc['name']] = val_series
            mfile_here.close()

        print(df_summary)
        if file_out is not None:
            df_summary.to_csv(file_out, sep=';')


def main():
    # fmdb = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/MDBs/MDB_S3A_B_OLCI_POLYMER_INSITU_20160401_20220531.nc'
    # do_check_mdb_times_impl(fmdb)
    # do_check_mdb_times()
    # do_check_extract_times()

    do_final_results()
    # do_final_results_l3()
    # do_final_results_CCI()
    # do_final_results_CNR()

    # path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/MDBs'
    # # name_mdb = 'MDB_CCIv6_INSITU_19970101_20221231.nc'
    # for name_mdb in os.listdir(path_base):
    #     if not name_mdb.endswith('nc'):
    #         continue
    #     do_chla(path_base, name_mdb)

    # path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/MDBs'
    # for name_csv in os.listdir(path_base):
    #     if not name_csv.endswith('csv'):
    #         continue
    #     # name_csv = 'MDB_CCIv6_INSITU_19970101_20221231.csv'
    #     do_chla_validation(path_base, name_csv)

    # retrieve_insitu_chla_from_paper_files()
    # compare_chlainsitu_with_paper()

    # check_dates_MDB()

    # make_validation_single_MDB()
    # make_variations_qc_single_MDB()

    # make_validation_list_MDB()
    # acnames = ['POLYMER', 'STANDARD', 'FUB', 'C2RCC']
    # platforms = ['A','B']
    # for ac in acnames:
    #     for p in platforms:
    #         make_validation_from_dfvalid(ac,p)

    # make_togheter_ab('SYKE')

    # bands = [400, 412, 443, 490, 510, 560, 620, 667, 779]
    # acnames = ['POLYMER', 'STANDARD', 'FUB', 'C2RCC']
    # name_out = 'SYKECOMBINATIONS'
    # bands = [412,443,490,510,560,667]
    # acnames = ['POLYMER', 'CCI']
    # name_out = 'POLYMER_CCI'

    # for band in bands:
    #     make_together_atm(band,acnames,name_out)
    # make_together_atm_params(bands,acnames,name_out)

    # acnames = ['POLYMER', 'STANDARD', 'FUB', 'C2RCC']
    # for ac in acnames:
    #     print(ac)
    #     make_average_spectra(ac, 'A')
    # for ac in acnames:
    #     launch_complete_sam(ac)

    # graficos()

    # get_table_stats(-1)

    # plot_sam()

    # check_c2rcc_flags()

    # COMPARISON AERONET-PANTHYR
    # do_comparison_aeronet_panthyr()
    # do_example()
    # do_chla()
    # do_chla_sat()
    # do_check_reflectances()
    # do_check_extract_times()

    # do_lps()

    # dir_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/OLCI_VEIT_UPDATED/WFR'
    # name_extract = 'S3A_OL_2_WFR____20210102T095605_20210102T095610_20210906T082804__trim_EXT_067_022_MAR_R_NT_003_SEN3_extract_VEIT.nc'
    # name_result = 'S3A_OL_2_WFR____20210102T095605_20210102T095610_20210906T082804__trim_EXT_067_022_MAR_R_NT_003.SEN3/Oa01_reflectance.nc'
    #
    # file_extract = os.path.join(dir_base, 'extracts', name_extract)
    # file_result = os.path.join(dir_base, 'results', name_result)
    #
    # from netCDF4 import Dataset
    # from datetime import datetime
    # from datetime import timedelta
    #
    # nc_extract = Dataset(file_extract)
    # dtfin = datetime(1970, 1, 1) + timedelta(seconds=float(nc_extract.variables['satellite_time'][0]))
    # print('Extract: ',dtfin)
    # print('Extract value: ',float(nc_extract.variables['satellite_time'][0]))
    # dtfin_alt = datetime.fromtimestamp(float(nc_extract.variables['satellite_time'][0]))
    # print('Extract alt ',dtfin_alt)
    #
    #
    # nc_sat = Dataset(file_result)
    # print('Sat: ',nc_sat.start_time)
    # print('Sat time py ', datetime.strptime(nc_sat.start_time, "%Y-%m-%dT%H:%M:%S.%fZ"))
    # val = float(datetime.strptime(nc_sat.start_time, "%Y-%m-%dT%H:%M:%S.%fZ").timestamp())
    # print('Sat time value: ', val)
    # print('Conversion: ',datetime.fromtimestamp(val))
    #
    # dir_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/OLCI_VEIT_UPDATED/MDBs_20052022'
    # file_s3a = os.path.join(dir_base, 'MDB_S3A_OLCI_WFR_STANDARD_L2_HYPERNETS_VEIT.nc')
    # mdb = MDBFile(file_s3a)
    # mdb.load_mu_datav2(4)

def do_final_results_CNR():
    print('CCI Results')
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/MED_MATCHUPS/MDBs'
    name_mdb = f'MDB___1KM_MULTI_L2_AERONET_Casablanca_Platform.nc'
    fmdb = os.path.join(path_base,name_mdb)
    if os.path.exists(fmdb):
        make_validation_single_MDB(path_base, name_mdb)

def do_final_results_CCI():
    print('CCI Results')
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    sites = ['Gustav_Dalen_Tower', 'Helsinki_Lighthouse', 'Irbe_Lighthouse']

    # SINGLE VALIDATION FOR SITE
    for site in sites:
        path_mdb = os.path.join(path_base,'CCIv6')
        name_mdb = f'MDB___1KM_CCI_L2_AERONET_{site}.nc'
        fmdb = os.path.join(path_mdb,name_mdb)
        if os.path.exists(fmdb):
            make_validation_single_MDB(path_mdb, name_mdb)

    # PREPARE DF CSV COMBINING ALL THE STATIONS AND VALIDATING
    # make_validation_list_MDB('CCIv6', '')
    make_validation_from_dfvalid('CCIv6', '')



    # # SCATTER PLOT FOR BANDS COMBINING AC PROCESSORS
    # bands = [412, 443, 490, 510, 560, 665]
    # for band in bands:
    #     make_together_atm(band,acnames,'AB','ACCombinations_StandardAB_PolymerAB_CCIv6')

    # # PARAMETERS BY BAND
    acnames = ['STANDARD', 'POLYMER', 'CCIv6']
    bands = [412, 443, 490, 510, 560, 665]
    make_together_atm_params(bands,acnames,'AB','Params_StandardAB_PolymerAB_CCIv6')
    acnames = ['POLYMER', 'CCIv6']
    make_together_atm_params(bands, acnames, 'AB', 'Params_PolymerAB_CCIv6')



    # AVERAGE SPECTRA
    # for ac in acnames:
    #     print(ac)
    #     make_average_spectra(ac, 'AB', 'SPECTRA')

    # STATS TABLE
    bands = [412, 443, 490, 510, 560, 665]
    acs = ['STANDARD', 'POLYMER', 'CCIv6']
    params = {
        'N': 0,
        'SLOPE': 11,
        'OFFSET': 12,
        'XAVG': 13,
        'YAVG': 14,
        'DETER(r2)': 10,
        'RMSE': 6,
        'CPRMSE': 15,
        'BIAS': 9,
        'PCC(r)': 3,
        'RPD': 7,
        'APD': 8,
        'MAE': 16
    }
    get_table_stats(-1, 'AB', acs, params)
    for wl in bands:
        get_table_stats(wl, 'AB', acs, params)
    get_table_stats_complete(bands, acs, 'AB', params)

    # MU INFO
    acnames = ['CCIv6']
    # make_mu_info(path_base, acnames, '')
    # make_mu_info_bytower(path_base, 'CCIv6', 2005, 2021)


def do_final_results_l3():
    print('L3 Results')
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    sites = ['Gustav_Dalen_Tower', 'Helsinki_Lighthouse', 'Irbe_Lighthouse']
    # acnames = ['STANDARD', 'POLYMER', 'C2RCC', 'FUB']
    acnames = ['STANDARD', 'POLYMER']
    platform = 'AB'
    # SINGLE VALIDATIONS FOR  AC/SITE
    for ac in acnames:
        for site in sites:
            path_mdb, name_mdb, fmdb = get_file_mdb(path_base, platform, site, ac,3)
            print(f'MDB: {fmdb}')
            if os.path.exists(fmdb):
                make_validation_single_MDB(path_mdb, name_mdb)

    # PREPARE DF CSV COMBINING ALL THE STATIONS AND VALIDATING
    for ac in acnames:
        # make_validation_list_MDB(ac, 'AB')
        make_validation_from_dfvalid(ac, 'AB')

    # SCATTER PLOT FOR BANDS COMBINING AC PROCESSORS
    # bands = [400, 412, 443, 490, 510, 560, 620, 667, 779]
    # for band in bands:
    #     make_together_atm(band,acnames,'AB','ACCombinations_S3AB')

    # PARAMETERS BY BAND
    bands = [400, 412, 443, 490, 510, 560, 620, 667, 779]
    make_together_atm_params(bands,acnames,'AB','ParamsS3AB')

    # AVERAGE SPECTRA
    # for ac in acnames:
    #     print(ac)
    #     make_average_spectra(ac, 'AB', 'SPECTRA')

    # STATS TABLE
    bands = [400, 412, 443, 490, 510, 560, 620, 667, 779]
    acs = ['STANDARD', 'POLYMER']
    params = {
        'N': 0,
        'SLOPE': 11,
        'OFFSET': 12,
        'XAVG': 13,
        'YAVG': 14,
        'DETER(r2)': 10,
        'RMSE': 6,
        'CPRMSE': 15,
        'BIAS': 9,
        'PCC(r)': 3,
        'RPD': 7,
        'APD': 8,
        'MAE': 16
    }
    get_table_stats(-1, 'AB', acs, params)
    for wl in bands:
        get_table_stats(wl, 'AB', acs, params)
    get_table_stats_complete(bands, acs, 'AB', params)

    # MUINfo
    # make_mu_info(path_base, acnames, 'AB')
    # make_mu_info_bytower(path_base, 'POLYMER', 2016, 2021)


def do_final_results():
    print('FINAL RESULTS')
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    platforms = ['A', 'B']
    sites = ['Gustav_Dalen_Tower', 'Helsinki_Lighthouse', 'Irbe_Lighthouse']
    acnames = ['STANDARD', 'POLYMER', 'C2RCC', 'FUB']
    acnames = ['STANDARD', 'POLYMER']


    # SINGLE VALIDATIONS FOR PLATFORM/AC/SITE
    # for ac in acnames:
    #     for platform in platforms:
    #         for site in sites:
    #             path_mdb, name_mdb, fmdb = get_file_mdb(path_base, platform, site, ac, 2)
    #             print(f'MDB: {fmdb}')
    #             if os.path.exists(fmdb):
    #                 make_validation_single_MDB(path_mdb, name_mdb)

    # PREPARE DF CSV COMBINING ALL THE STATIONS AND VALIDATING
    platforms = ['A']
    acnames = ['POLYMER']
    for ac in acnames:
        for platform in platforms:
            #make_validation_list_MDB(ac, platform)
            make_validation_from_dfvalid(ac, platform)

    # SCATTER PLOT FOR BANDS COMBINING AC PROCESSORS
    # bands = [412, 443, 490, 510, 560, 620, 667]
    # for band in bands:
    #     make_together_atm(band,acnames,'A','ACCombinations_S3A')
    #     make_together_atm(band,acnames,'B', 'ACCombinations_S3B')
    # acnamesnofub = ['STANDARD', 'POLYMER', 'C2RCC']
    # bandsnofub = [400,779]
    # for band in bandsnofub:
    #     make_together_atm(band,acnamesnofub,'A','ACCombinations_S3A')
    #     make_together_atm(band,acnamesnofub,'B', 'ACCombinations_S3B')

    # PARAMETERS BY BAND
    # bands = [412, 443, 490, 510, 560, 620, 667]
    # make_together_atm_params(bands,acnames,'A','ParamsS3A')
    # make_together_atm_params(bands, acnames,'B','ParamsS3B')

    # AVERAGE SPECTRA
    # acnames = ['POLYMER', 'STANDARD', 'FUB', 'C2RCC']
    # for ac in acnames:
    #     for platform in platforms:
    #         print(ac)
    #         make_average_spectra(ac, platform, 'SPECTRA')

    # AB COMPARISON
    # for ac in acnames:
    #     make_togheter_ab(ac)

    # STATS TABLE
    # bands = [400, 412, 443, 490, 510, 560, 620, 667, 779]
    # acs = ['STANDARD', 'POLYMER']
    # params = {
    #     'N': 0,
    #     'SLOPE': 11,
    #     'OFFSET': 12,
    #     'XAVG': 13,
    #     'YAVG': 14,
    #     'DETER(r2)': 10,
    #     'RMSE': 6,
    #     'CPRMSE': 15,
    #     'BIAS': 9,
    #     'PCC(r)': 3,
    #     'RPD': 7,
    #     'APD': 8,
    #     'MAE': 16
    # }
    # get_table_stats(-1, 'A', acs, params)
    # get_table_stats(-1, 'B', acs, params)
    # for wl in bands:
    #     get_table_stats(wl, 'A', acs, params)
    #     get_table_stats(wl, 'B', acs, params)
    # get_table_stats_complete(bands, acs, 'A', params)
    # get_table_stats_complete(bands, acs, 'B', params)

    # MUINFO
    # make_mu_info(path_base, acnames, 'AB')


def get_file_mdb(path_base, platform, site, acname, level):
    path_mdb = os.path.join(path_base, acname)
    res = 'EFR'
    if acname == 'STANDARD':
        res = 'WFR'
    name_mdb = f'MDB_S3{platform}_OLCI_{res}_{acname}_L{level}_AERONET_{site}.nc'
    fmdb = os.path.join(path_mdb, name_mdb)
    return path_mdb, name_mdb, fmdb


def check_dates_MDB():
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    platforms = ['A', 'B']
    sites = ['Gustav_Dalen_Tower', 'Helsinki_Lighthouse', 'Irbe_Lighthouse']
    acnames = ['STANDARD', 'POLYMER', 'C2RCC', 'FUB']
    date_ini = dt(2016, 4, 29)
    date_fin = dt(2022, 5, 31)
    delta = date_fin - date_ini
    ndays = delta.days + 1

    # for platform in platforms:
    #     for site in sites:
    #
    #         fout = os.path.join(path_base, 'INFO_DATES', f'Dates_S3{platform}_{site}.csv')
    #
    #         # acnames = ['STANDARD', 'POLYMER', 'C2RCC']
    #         df = pd.DataFrame(index=list(range(ndays)),
    #                           columns=['Date', 'STANDARD', 'POLYMER', 'C2RCC', 'FUB', 'Total'])
    #         for idx in range(ndays):
    #             df.loc[idx, 'Date'] = (date_ini + timedelta(days=idx)).strftime('%Y-%m-%d')
    #             for ac in acnames:
    #                 df.loc[idx, ac] = 0
    #             df.loc[idx, 'Total'] = 0
    #
    #         for ac in acnames:
    #             res = 'EFR'
    #             if ac == 'STANDARD':
    #                 res = 'WFR'
    #             path_mdb = os.path.join(path_base, ac)
    #             name_mdb = f'MDB_S3{platform}_OLCI_{res}_{ac}_L2_AERONET_{site}.nc'
    #             fmdb = os.path.join(path_mdb, name_mdb)
    #             print(fmdb, os.path.exists(fmdb))
    #             mdbfile = MDBFile(fmdb)
    #             for time in mdbfile.sat_times:
    #                 date_here = time
    #                 idx = (date_here - date_ini).days
    #                 df.loc[idx, ac] = df.loc[idx, ac] + 1
    #
    #         for idx in range(ndays):
    #             for ac in acnames:
    #                 df.loc[idx, 'Total'] = df.loc[idx, 'Total'] + df.loc[idx, ac]
    #
    #         df.to_csv(fout, sep=';')

    for site in sites:
        fcsv_a = os.path.join(path_base, 'INFO_DATES', f'Dates_S3A_{site}.csv')
        fcsv_b = os.path.join(path_base, 'INFO_DATES', f'Dates_S3B_{site}.csv')
        df_a = pd.read_csv(fcsv_a, sep=';')
        df_b = pd.read_csv(fcsv_b, sep=';')
        df = pd.DataFrame(index=list(range(ndays)), columns=['Date', 'S3A', 'S3B', 'Total'])
        for idx in range(ndays):
            df.loc[idx, 'Date'] = (date_ini + timedelta(days=idx)).strftime('%Y-%m-%d')
            df.loc[idx, 'S3A'] = 0
            df.loc[idx, 'S3B'] = 0
            if df_a.loc[idx, 'Total'] == 4:
                df.loc[idx, 'S3A'] = 1
            if df_b.loc[idx, 'Total'] == 4:
                df.loc[idx, 'S3B'] = 1
            df.loc[idx, 'Total'] = df.loc[idx, 'S3A'] + df.loc[idx, 'S3B']
        fout = os.path.join(path_base, 'INFO_DATES', f'Dates_S3_{site}.csv')
        df.to_csv(fout, sep=';')

    for platform in platforms:
        fcsv_gdt = os.path.join(path_base, 'INFO_DATES', f'Dates_S3{platform}_Gustav_Dalen_Tower.csv')
        fcsv_hlh = os.path.join(path_base, 'INFO_DATES', f'Dates_S3{platform}_Helsinki_Lighthouse.csv')
        fcsv_ilh = os.path.join(path_base, 'INFO_DATES', f'Dates_S3{platform}_Irbe_Lighthouse.csv')
        df_gdt = pd.read_csv(fcsv_gdt, sep=';')
        df_hlh = pd.read_csv(fcsv_hlh, sep=';')
        df_ilh = pd.read_csv(fcsv_ilh, sep=';')
        df = pd.DataFrame(index=list(range(ndays)),
                          columns=['Date', 'Gustav_Dalen_Tower', 'Helsinki_Lighthouse', 'Irbe_Lighthouse'])
        for idx in range(ndays):
            df.loc[idx, 'Date'] = (date_ini + timedelta(days=idx)).strftime('%Y-%m-%d')
            df.loc[idx, 'Gustav_Dalen_Tower'] = 0
            df.loc[idx, 'Helsinki_Lighthouse'] = 0
            df.loc[idx, 'Irbe_Lighthouse'] = 0
            if df_gdt.loc[idx, 'Total'] == 4:
                df.loc[idx, 'Gustav_Dalen_Tower'] = 1
            if df_hlh.loc[idx, 'Total'] == 4:
                df.loc[idx, 'Helsinki_Lighthouse'] = 1
            if df_ilh.loc[idx, 'Total'] == 4:
                df.loc[idx, 'Irbe_Lighthouse'] = 1
            df.loc[idx, 'Total'] = df.loc[idx, 'Gustav_Dalen_Tower'] + df.loc[idx, 'Helsinki_Lighthouse'] + df.loc[
                idx, 'Irbe_Lighthouse']
        fout = os.path.join(path_base, 'INFO_DATES', f'Dates_S3{platform}.csv')
        df.to_csv(fout, sep=';')


def do_check_mdb_times():
    # path_mdb = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDB'
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    # platforms = ['A', 'B']

    # acnames = ['STANDARD', 'POLYMER', 'C2RCC', 'FUB']
    acnames = ['POLYMER']
    platforms = ['A', 'B']
    sites = ['Irbe_Lighthouse']
    hours_array = []
    for platform in platforms:
        for site in sites:
            for ac in acnames:
                res = 'EFR'
                if ac == 'STANDARD':
                    res = 'WFR'
                path_mdb = os.path.join(path_base, ac)
                name_mdb = f'MDB_S3{platform}_OLCI_{res}_{ac}_L2_AERONET_{site}.nc'
                fmdb = os.path.join(path_mdb, name_mdb)
                print(fmdb, os.path.exists(fmdb))
                mdbfile = MDBFile(fmdb)
                for time in mdbfile.sat_times:
                    date_ref = time.replace(hour=0, minute=0, second=0, microsecond=0)
                    hours = (time - date_ref).total_seconds() / 3600
                    print(hours)
                    hours_array.append(hours)
    mean_hour = np.mean(np.asarray(hours_array, dtype=float))
    print('--------')
    date_mean_site = dt.now().replace(hour=0, minute=0, second=0, microsecond=0) + timedelta(hours=mean_hour)
    print(date_mean_site)


def do_check_mdb_times_impl(fmdb):
    hours_array = []
    mdbfile = MDBInSituFile(fmdb)
    for time in mdbfile.sat_times:
        date_ref = time.replace(hour=0, minute=0, second=0, microsecond=0)
        hours = (time - date_ref).total_seconds() / 3600
        print(hours)
        hours_array.append(hours)
    mean_hour = np.mean(np.asarray(hours_array, dtype=float))
    print('--------')
    date_mean_site = dt.now().replace(hour=0, minute=0, second=0, microsecond=0) + timedelta(hours=mean_hour)
    print(date_mean_site)


def do_check_extract_times():
    path_extracts = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/Gustav_Dalen_Tower/polymer/extractsAB'
    # path_extracts = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/CCI/extractsv4'
    from netCDF4 import Dataset
    from datetime import timedelta
    for name in os.listdir(path_extracts):
        if not name.endswith('nc'):
            continue
        fname = os.path.join(path_extracts, name)
        nc = Dataset(fname)
        dthere = dt.utcfromtimestamp(float(nc.variables['satellite_time'][0]))
        dthereotro = dt(1970, 1, 1) + timedelta(seconds=float(nc.variables['satellite_time'][0]))
        dtherebis = get_sat_time_from_fname(name)
        # dtherebis = get_sat_time_fromcci_fname(name)
        dif = abs((dthere - dtherebis).total_seconds())
        print(dthere, dtherebis, dif)
        if dif > 1:
            print('ERROR', dthere, dthereotro, dtherebis)
        # if dthere.year == 2016:
        #     print(dthere, dtherebis, dif)
        #     if abs(dif) < 1:
        #         print('---------------------------------------->GOOD CUA')


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


def get_sat_time_fromcci_fname(fname):
    sat_time = None
    try:
        sat_time = dt.strptime(fname[1:8], '%Y%j')
    except ValueError:
        pass
    return sat_time


def do_lps():
    dir_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/OLCI_VEIT_UPDATED/MDBs_20052022'
    dir_ab = os.path.join(dir_base, 'MDB_S3AB')

    # Single validation (scaterplot for each platform)
    file_s3a = os.path.join(dir_base, 'MDB_S3A_OLCI_WFR_STANDARD_L2_HYPERNETS_VEIT.nc')
    file_s3b = os.path.join(dir_base, 'MDB_S3B_OLCI_WFR_STANDARD_L2_HYPERNETS_VEIT.nc')
    # make_validation_single_MDB(dir_base,file_s3a)
    # make_validation_single_MDB(dir_base, file_s3b)

    ##Concatenate dataframes A + B
    df_all_a = pd.read_csv(os.path.join(dir_ab, 'DataS3A.csv'), sep=';')
    df_all_b = pd.read_csv(os.path.join(dir_ab, 'DataS3B.csv'), sep=';')
    df_valid_a = pd.read_csv(os.path.join(dir_ab, 'DataValidS3A.csv'), sep=';')
    df_valid_b = pd.read_csv(os.path.join(dir_ab, 'DataValidS3B.csv'), sep=';')
    df_all = pd.concat([df_all_a, df_all_b], ignore_index=True)
    df_all.to_csv(os.path.join(dir_ab, 'DataAll.csv'), sep=';')
    df_valid = pd.concat([df_valid_a, df_valid_b], ignore_index=True)
    df_valid.to_csv(os.path.join(dir_ab, 'DataValid.csv'), sep=';')

    ##Scaterplot A + B
    do_example(os.path.join(dir_ab, 'DataValid.csv'), os.path.join(dir_ab, 'ScatterPlotAB.jpg'))

    ## Param graphics
    path_a = os.path.join(dir_base, 'MDB_S3A_OLCI_WFR_STANDARD_L2_HYPERNETS_VEIT')
    path_b = os.path.join(dir_base, 'MDB_S3B_OLCI_WFR_STANDARD_L2_HYPERNETS_VEIT')
    make_graphics_param_ab(path_a, path_b, dir_ab)

    filein = os.path.join(dir_ab, 'MUByMonth.csv')
    fileout = os.path.join(dir_ab, 'MUByMonth.jpg')
    make_bargraphics_matchups(filein, fileout)


def do_chla(path_base, name_mdb):
    name_out = name_mdb[:-3]
    file_out = os.path.join(path_base, f'{name_out}.csv')
    if os.path.exists(file_out):
        return
    file_mdb = os.path.join(path_base, name_mdb)
    baltic_mlp_code = '/home/lois/PycharmProjects/aeronet'
    mfile = MDBInSituFile(file_mdb)
    mfile.delta_t = (2*3600)

    if not mfile.VALID:
        return
    mfile.start_baltic_chla(baltic_mlp_code, 'insitu_CHLA')
    mfile.qc_sat.max_diff_wl = 6
    mfile.qc_sat.window_size = 3
    mfile.qc_sat.min_valid_pixels = 9
    mfile.qc_sat.apply_outliers = False
    mfile.qc_sat.stat_value = 'median'

    mfile.set_wlsatlist_aswlref([412, 443, 490, 510, 555, 670])

    # for imu in range(mfile.n_mu_total):
    #     mfile.compute_baltic_chla_ensemble_mu(imu)
    mfile.prepare_df_validation()

    mfile.df_validation_valid.to_csv(file_out, sep=';')


def do_chla_validation(path_base, name_csv):
    fcsv = os.path.join(path_base, name_csv)
    path_out = os.path.join(path_base, name_csv[:-4])
    if not os.path.exists(path_out):
        os.mkdir((path_out))

    lcsv = name_csv[:-4].split('_')

    df = pd.read_csv(fcsv, sep=';')
    col_names = df.columns
    var_names = []
    col_names_stats = ['Param']
    for col in col_names:
        if col.startswith('ens') or col.startswith('chl'):
            var_names.append(col)
            col_names_stats.append(col)
    params = ['N', 'slope', 'intercept', 'r_value', 'p_value', 'std_err', 'rmse_val',
              'mean_rel_diff', 'mean_abs_rel_diff', 'bias', 'r2', 'slope_typeII', 'offset_typeII', 'XAVG', 'YAVG',
              'CPRMSE', 'MAE']
    df_valid_stats = pd.DataFrame(index=params, columns=col_names_stats)
    mdbplot = MDBPlot(None, None)
    mdbplot.set_units('chla')
    insitu_chla = df.loc[:, 'insitu_CHLA']
    for var in var_names:
        print(var)
        sat_chla = df.loc[:, var]
        sat_chla_here = sat_chla[~np.isnan(sat_chla)]
        insitu_chla_here = insitu_chla[~np.isnan(sat_chla)]

        mdbplot.xdata = insitu_chla_here
        mdbplot.ydata = sat_chla_here
        mdbplot.wldata = pd.Index(np.zeros(insitu_chla.shape))
        mdbplot.xlabel = r'Chl-a $_R$$_E$$_F$'
        mdbplot.ylabel = r'Chl-a $_S$$_A$$_T$'
        mdbplot.ylabel = f'{mdbplot.ylabel} ({var})'
        mdbplot.title = f'In situ chl-a vs. {lcsv[1]} ({var}) {lcsv[3]}-{lcsv[4]}'
        mdbplot.log_scale = True
        mdbplot.plot_scatter_plot(True, False, True, os.path.join(path_out, f'ScatterPlot_InSituChla_{var}.jpg'))
        mdbplot.compute_statistics()
        for param in mdbplot.valid_stats:
            df_valid_stats.loc[param, var] = mdbplot.valid_stats[param]

    file_results = os.path.join(path_out, f'Params.csv')
    df_valid_stats.to_csv(file_results, sep=';')


def do_chla_sat():
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/MDBs'
    # file_mdb = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/MDBs/MDB_S3A_B_OLCI_POLYMER_INSITU_20160401_20220531.nc'
    name_mdb = 'MDB___CCIv5_INSITU_19970101_20191231.nc'
    file_mdb = os.path.join(path_base, name_mdb)
    name_out = 'MDB___CCIv5_INSITU_CHLASAT_19970101_20191231.nc'
    file_out = os.path.join(path_base, name_out)
    mfile = MDBInSituFile(file_mdb)
    if not mfile.VALID:
        return

    mfile.add_baltic_chla(file_out, None, True)


def do_check_reflectances():
    path_extracts = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/CCI/extractsv4'
    path_extracts_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/CCI/extractsv4yyyyjjj'
    from netCDF4 import Dataset
    from datetime import datetime as dt
    import shutil

    # for name in os.listdir(path_extracts):
    #      if not name.endswith('nc'):
    #          continue
    #
    #      yearstr = name[1:5]
    #      jdaystr = name[5:8]
    #      yearpath = os.path.join(path_extracts_out,yearstr)
    #      if not os.path.exists(yearpath):
    #          os.mkdir(yearpath)
    #      jpath = os.path.join(yearpath,jdaystr)
    #      if not os.path.exists(jpath):
    #          os.mkdir(jpath)
    #      file_orig = os.path.join(path_extracts,name)
    #      file_dest = os.path.join(jpath,name)
    #      print(file_orig,file_dest)
    #      shutil.copy(file_orig,file_dest)

    # file_orig = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/CHLA_DATA/matchup_bal_cci_chl_surf6_orig__rrsmatch_w12CHL__FM.csv'
    # df = pd.read_csv(file_orig,sep=';')
    # var_names = ['date','time','lat','lon','CHL_in','rrs_412','rrs_443','rrs_490','rrs_510','rrs_555','rrs_670']
    # var_names_rrs = ['rrs_412','rrs_443','rrs_490','rrs_510','rrs_555','rrs_670']
    # dfrec = df.loc[:,var_names]
    # dfrec.loc[:,var_names_rrs] = -999
    # pix_pos = [0] * len(dfrec.index)
    # chlaref = -999
    # iref = -1
    # for index, row in dfrec.iterrows():
    #     chlastr = row['CHL_in']
    #     if chlastr!=chlaref:
    #         iref = 0
    #         pix_pos[index] = iref
    #     else:
    #         iref = iref +1
    #         pix_pos[index] = iref
    #     chlaref = chlastr
    # file_new = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/CHLA_DATA/matchup_bal_cci_chl_surf6_orig__rrsmatch_w12CHL__FM_pix.csv'
    # dfrec = pd.concat([dfrec,pd.Series(pix_pos)],axis=1)
    # dfrec.to_csv(file_new,sep=';')

    # file_orig = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/CHLA_DATA/matchup_bal_cci_chl_surf6_orig__rrsmatch_w12CHL__FM_pix.csv'
    # file_new = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/CHLA_DATA/matchup_bal_cci_chl_surf6_orig__rrsmatch_w12CHL__FM_comp.csv'
    # dfrec = pd.read_csv(file_orig, sep=';')
    # dset = None
    # rrs412 = None
    # rrs443 = None
    # rrs490 = None
    # rrs510 = None
    # rrs555 = None
    # rrs670 = None
    # sat_lat = None
    # sat_lon = None
    #
    # for index, row in dfrec.iterrows():
    #     ppos = int(row['PixPos'])
    #     chla_here = float(row['CHL_in'])
    #     chla_here_str = "{:.2f}".format(chla_here)
    #     datestr = str(int(row['date']))
    #
    #     if ppos == 0:
    #         print('-------------------------------')
    #         if dset is not None:
    #             print('Closing previous dset...')
    #             dset.close()
    #         dset = None
    #         dateref = dt.strptime(datestr, '%Y%m%d').strftime('%Y%j')
    #         yaearstr = dateref[0:4]
    #         jdaystr = dateref[4:7]
    #         path_here = os.path.join(path_extracts_out, yaearstr, jdaystr)
    #
    #         for name in os.listdir(path_here):
    #             fpath = os.path.join(path_here, name)
    #             dset = netCDF4.Dataset(fpath)
    #             chla = float(dset.variables['insitu_CHLA'][0])
    #             chlastr = "{:.2f}".format(chla)
    #             if chlastr == chla_here_str:
    #                 rrs412 = np.array(dset.variables['satellite_Rrs'][0, 0, 11:14, 11:14])
    #                 rrs412 = list(rrs412.flatten())
    #                 rrs443 = np.array(dset.variables['satellite_Rrs'][0, 1, 11:14, 11:14])
    #                 rrs443 = list(rrs443.flatten())
    #                 rrs490 = np.array(dset.variables['satellite_Rrs'][0, 2, 11:14, 11:14])
    #                 rrs490 = list(rrs490.flatten())
    #                 rrs510 = np.array(dset.variables['satellite_Rrs'][0, 3, 11:14, 11:14])
    #                 rrs510 = list(rrs510.flatten())
    #                 rrs555 = np.array(dset.variables['satellite_Rrs'][0, 4, 11:14, 11:14])
    #                 rrs555 = list(rrs555.flatten())
    #                 rrs670 = np.array(dset.variables['satellite_Rrs'][0, 5, 11:14, 11:14])
    #                 rrs670 = list(rrs670.flatten())
    #                 sat_lat = dset.variables['satellite_latitude'][0, 12, 12]
    #                 sat_lon = dset.variables['satellite_longitude'][0, 12, 12]
    #                 break
    #             else:
    #                 dset.close()
    #                 dset = None
    #     dsetvalid = dset is not None
    #     if dsetvalid:
    #         dfrec.loc[index, 'rrs_412'] = rrs412[ppos]
    #         dfrec.loc[index, 'rrs_443'] = rrs443[ppos]
    #         dfrec.loc[index, 'rrs_490'] = rrs490[ppos]
    #         dfrec.loc[index, 'rrs_510'] = rrs510[ppos]
    #         dfrec.loc[index, 'rrs_555'] = rrs555[ppos]
    #         dfrec.loc[index, 'rrs_670'] = rrs670[ppos]
    #         dfrec.loc[index, 'sat_lat'] = sat_lat
    #         dfrec.loc[index, 'sat_lon'] = sat_lon

    # for idx in range(6):
    #     rrs = np.array(nc.variables['satellite_Rrs'][0,idx,11:14,11:14])

    # central_r, central_c, r_s, r_e, c_s, c_e = get_dimensions(nc.variables['satellite_Rrs'],3)
    # print(central_r,central_c,r_s,r_e,c_s,c_e)

    # dfrec.to_csv(file_new, sep=';')


def retrieve_insitu_chla_from_paper_files():
    path = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/CHLA_DATA/BRANDO2021_CHLA'
    # file_orig = os.path.join(path, 'matchup_bal_cci_chl_surf6_orig__rrsmatch_w12CHL__AlgLine.csv')
    # file_out = os.path.join(path, 'AlgLine_Paper_chlainsitu.csv')
    file_orig = os.path.join(path, 'matchup_bal_cci_chl_surf6_orig__rrsmatch_w12CHL__COMBINE.csv')
    file_out = os.path.join(path, 'Combine_Paper_chlainsitu.csv')
    df = pd.read_csv(file_orig, sep=';')
    dfnew = df.copy()
    chlaref = -999
    latref = -999
    lonref = -999
    index_df = 0
    for index, row in df.iterrows():
        chlastr = row['CHL_in']
        latstr = row['lat']
        lonstr = row['lon']
        if chlastr != chlaref and latstr != latref and lonstr != lonref:
            dfnew.iloc[index_df, :] = row[:]
            chlaref = chlastr
            latref = latstr
            lonref = lonstr
            index_df = index_df + 1
    dfnew = dfnew[:index_df]
    col_names = ['date', 'time', 'lat', 'lon', 'CHL_in']
    dfnew = dfnew.loc[:, col_names]
    for index, row in dfnew.iterrows():
        # strdate = str(int(row['date']))
        # datehere = dt.strptime(strdate, '%Y%m%d')
        datehere = dt.strptime(row['date'], '%d/%m/%Y')
        datehere = datehere.replace(hour=0, minute=0, second=0, microsecond=0) + timedelta(hours=row['time'])
        strdatenew = datehere.strftime('%Y-%m-%d')
        strtimenew = datehere.strftime('%H:%M')
        dfnew.loc[index, 'date'] = strdatenew
        dfnew.loc[index, 'time'] = strtimenew

    dfnew.to_csv(file_out, sep=';')


def compare_chlainsitu_with_paper():
    path = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/CHLA_DATA'
    filebal = os.path.join(path, 'Baltic_CHLA_All.csv')
    filepaper = os.path.join(path, 'BRANDO2021_CHLA', 'All_chlainsitu.csv')
    fileout = os.path.join(path, 'BRANDO2021_CHLA', 'All_chlainsitu_allidx.csv')

    datelist = []
    chlalist = []
    idxlist = []
    dfbal = pd.read_csv(filebal, sep=';')
    for index, row in dfbal.iterrows():
        datestr = row['DATE']
        timestr = row['HOUR']
        strdate = f'{datestr}T{timestr}'
        datehere = dt.strptime(strdate, '%d/%m/%YT%H:%M:%S')
        datelist.append(datehere)
        chlalist.append(float(row['CHLA']))
        idxlist.append(int(row['IDX']))

    dfpaper = pd.read_csv(filepaper, sep=';')
    dfnew = dfpaper.copy()
    idxbal = [-1] * len(dfpaper.index)
    for index, row in dfpaper.iterrows():
        datestr = row['date']
        timestr = row['time']
        strdate = f'{datestr}T{timestr}'
        datepaper = dt.strptime(strdate, '%d/%m/%YT%H:%M')
        chlapaper = float(row['CHL_in'])

        for idx in range(len(datelist)):
            datebal = datelist[idx]
            chlabal = chlalist[idx]
            tdif = abs((datepaper - datebal).total_seconds())
            cdif = abs(chlapaper - chlabal)
            if tdif < 120 and cdif < 0.000001:
                idxbal[index] = idxlist[idx]
    dfnew['IDX_BAL'] = pd.Series(idxbal)
    dfnew.to_csv(fileout,sep=';')


def get_dimensions(satellite_rrs, window_size):
    # Dimensions
    nrows = satellite_rrs.shape[2]
    ncols = satellite_rrs.shape[3]
    central_r = int(np.floor(nrows / 2))
    central_c = int(np.floor(ncols / 2))
    r_s = central_r - int(np.floor(window_size / 2))  # starting row
    r_e = central_r + int(np.floor(window_size / 2)) + 1  # ending row
    c_s = central_c - int(np.floor(window_size / 2))  # starting col
    c_e = central_c + int(np.floor(window_size / 2)) + 1  # ending col
    return central_r, central_c, r_s, r_e, c_s, c_e


def do_example(file_csv, file_out):
    # file_csv = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/OLCI_VEIT_UPDATED/MDBs/DataValid.csv'
    df_valid = pd.read_csv(file_csv, sep=';')
    # file_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/OLCI_VEIT_UPDATED/MDBs/ComparisonAB.jpg'
    h = plt.figure()
    sns.set_theme()
    sns.lmplot(data=df_valid, x="Ins_Rrs", y="Sat_Rrs", hue="Platform", truncate=False, scatter_kws={"s": 10})
    xmin, xmax = plt.gca().get_xlim()
    ymin, ymax = plt.gca().get_ylim()
    xmin = np.min([xmin, ymin])
    xmax = np.min([xmax, ymax])
    if xmin < 0:
        xmin = 0
    plt.gca().set_xlim([xmin, xmax])
    plt.gca().set_ylim([xmin, xmax])
    plt.plot([xmin, xmax], [xmin, xmax], '--k')

    plt.title(f'Venise')
    plt.savefig(file_out, dpi=300)
    plt.close(h)


def do_comparison_aeronet_panthyr():
    file_csv = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/OLCI_VEIT_UPDATED/PANTHYR/VALID_SPECTRA/Validation_Hypstar_Panthir_270.csv'
    path_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/OLCI_VEIT_UPDATED/PANTHYR/VALID_SPECTRA/'
    df_val = pd.read_csv(file_csv, sep=';')
    wldata = df_val[:]['Wavelength']
    wllist = list(wldata.unique())
    mplot = MDBPlot(None, df_val)
    mplot.make_validation_dfval_insitu(path_out, 'PANTHYR-HYPSTAR Apr-Dec 2021', 'PANTHYR-HYPSTAR', wllist)

    # file_csv = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/OLCI_VEIT_UPDATED/PANTHYR/VALID_SPECTRA/Validation_Hypstar_Panthir_270_AvgSpectra.csv'
    # df_new = pd.read_csv(file_csv,sep=';')
    # h = plt.figure()
    # sns.set_theme()
    # sns.relplot(data=df_new, kind='line', x="Wavelength", y="RRS", hue="TYPE", style="TYPE", ci='sd')
    # plt.title(f'PANTHYR-HYPSTAR Apr-Dec 2021')
    # plt.savefig(os.path.join(path_out, f'PANTHYR-HYPSTAR_AvgSpectra.jpg'), dpi=300)
    # plt.close(h)


def check_c2rcc_flags():
    file_mdb = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/C2RCC/MDB_S3A_OLCI_EFR_C2RCC_L2_AERONET_Gustav_Dalen_Tower.nc'
    mfile = MDBFile(file_mdb)
    mfile.prepare_df_validation()


def plot_sam():
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/SAM/'
    file_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/SAM/ALLSAM.csv'
    dfsam = pd.read_csv(file_base, sep=';')

    h = plt.figure()
    print(dfsam.shape)

    acnames = ['STANDARD', 'POLYMER', 'C2RCC', 'FUB']
    platforms = ['A', 'B']
    col_names = []
    allarray = []
    for ac in acnames:
        for p in platforms:
            cname = f'{ac} S3{p}'
            array = np.array(dfsam.loc[:, cname])
            array = array[~np.isnan(array)]
            allarray.append(array)
            col_names.append(f'{ac} {p}')

    plt.boxplot(allarray)
    plt.xticks(ticks=range(1, 9), labels=col_names, fontsize='xx-small')
    plt.xlabel('Atmospheric correction')
    plt.ylabel('SAM')
    plt.savefig(os.path.join(path_base, 'SamBoxPlot.jpg'), dpi=300)
    plt.close(h)


def get_table_stats(wlref, platform, acs, params):
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    path_out = os.path.join(path_base, 'TABLE_PARAM')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    # acs = ['POLYMER', 'CCI']
    # params = ['N', 'r2', 'RMSD', 'APD', 'RPD', 'bias']
    table = pd.DataFrame(columns=acs, index=params.keys())
    for ac in acs:
        postname = f'_S3{platform}'
        if ac == 'CCIv6':
            postname = ''
        path_ac = os.path.join(path_base, ac, f'MDB_{ac}{postname}', 'Params.csv')
        table_ac = pd.read_csv(path_ac, sep=';')
        index = -1
        if wlref == -1:
            index = 2
            name_out = f'ParamAll_S3{platform}.csv'
        else:

            wllist = list(table_ac.columns[3:].astype(dtype=np.float64))
            index_wl, wl_aca = get_index_wl_list(wlref, wllist)
            if index_wl >= 0:
                index = index_wl + 3
                name_out = f'Param{wlref}_{platform}.csv'
        if index >= 0:
            for param in params.keys():
                irow_ac = params[param]
                table.loc[param, ac] = table_ac.loc[irow_ac].iat[index]

            file_out = os.path.join(path_out, name_out)
            table.to_csv(file_out, sep=';')


def get_table_stats_complete(bands, acs, platform, params):
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    path_out = os.path.join(path_base, 'TABLE_PARAM')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    nrows = len(params) * len(acs)
    col_names = ['Param', 'AC', 'All']
    for wl in bands:
        col_names.append(str(wl))
    table = pd.DataFrame(columns=col_names, index=range(nrows))

    idx = 0
    idxs = {}
    for param in params.keys():
        for ac in acs:
            table.loc[idx, 'Param'] = param
            table.loc[idx, 'AC'] = ac
            if param not in idxs.keys():
                idxs[param] = {
                    ac: idx
                }
            else:
                idxs[param][ac] = idx
            idx = idx + 1

    for ac in acs:
        postname = f'_S3{platform}'
        if ac == 'CCIv6':
            postname = ''
        path_ac = os.path.join(path_base, ac, f'MDB_{ac}{postname}', 'Params.csv')
        table_ac = pd.read_csv(path_ac, sep=';')
        for param in params.keys():
            idx = idxs[param][ac]
            irow_ac = params[param]
            icol_ac = 2
            table.loc[idx, 'All'] = table_ac.loc[irow_ac].iat[icol_ac]
            wllist = list(table_ac.columns[3:].astype(dtype=np.float64))
            for wl in bands:
                wls = str(wl)
                icol_ac = -1
                index_wl, wl_aca = get_index_wl_list(wl, wllist)
                if index_wl >= 0:
                    icol_ac = index_wl + 3
                if icol_ac >= 0:
                    table.loc[idx, wls] = table_ac.loc[irow_ac].iat[icol_ac]

    file_out = os.path.join(path_out, f'TableAllParams_S3{platform}.csv')
    table.to_csv(file_out, sep=';')


def make_validation_list_MDB(acname, platform):
    path_ini = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    mdblist = MDBFileList()
    path_base = os.path.join(path_ini, acname)
    if not os.path.exists(path_base):
        return
    wcard = f'S3{platform}'
    if acname == 'CCIv6':
        wcard = 'CCI'

    nfiles = 0
    for f in os.listdir(path_base):
        if (wcard is not None and f.find(wcard) > 0) and f.endswith('.nc'):
            path_mdb = os.path.join(path_base, f)
            mdblist.add_mdb_file(path_mdb)
            nfiles = nfiles + 1
    if nfiles == 0:
        print(f'No files found for: {acname} {platform}')
        return

    if acname == 'CCIv6':
        path_out = os.path.join(path_base, f'MDB_{acname}')
    else:
        path_out = os.path.join(path_base, f'MDB_{acname}_S3{platform}')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    wllist = [400, 412, 443, 490, 510, 560, 620, 665, 779]  # -> POLYMER, STANDARD, C2RCC
    if acname == 'FUB':
        wllist = [412, 443, 490, 510, 560, 620, 665]
    if acname == 'CCIv6':
        wllist = [412, 443, 490, 510, 560, 665]

    mdblist.set_wl_ref(wllist)
    ##qc_base not implemented. See prepared_df_for_validation in MDBFileList to check QC
    # mdblist.qc_base.set_sat_eumetsat_defaults(3)

    mdblist.prepare_df_for_validation()
    mdblist.save_df_validation_to_file(path_out)

    dfval, groups, start_date, end_date = mdblist.get_df_validation(None, None, False)

    mplot = MDBPlot(None, dfval)
    mplot.save_validation_dfval_data(path_out)

    # mplot.obtain_mu_info(mdblist.df_validation,dfval,path_out)
    # mplot.plot_all_scatter_plot(path_out)

    # dates = start_date.strftime('%Y-%m-%d') + '_' + end_date.strftime('%Y-%m-%d')
    # file_name_base = f'CCI_BalTower_{dates}'
    # mplot.make_validation_dfval(path_out, None, file_name_base, mdblist.wlref)


def make_validation_single_MDB(path_base, name_mdb):
    # path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/CCI/MDBs'
    # name_mdb = 'MDB___1KM_CCI_L2_AERONET_Gustav_Dalen_Tower.nc'
    # name_mdb = 'MDB___1KM_CCI_L2_AERONET_Helsinki_Lighthouse.nc'

    path_mdb = os.path.join(path_base, name_mdb)
    reader = MDB_READER(path_mdb, True)
    ##AERONET WAVELENGHTS
    wllist = [400, 412, 443, 490, 510, 560, 620, 665, 779]
    if name_mdb.find('FUB') > 0:
        wllist = [412, 443, 490, 510, 560, 620, 665]
    if name_mdb.find('CCI') > 0 or name_mdb.find('MULTI')>0:
        wllist = [412, 443, 490, 510, 560, 665]
        reader.mfile.set_hour_sat_time(11, 0)

    reader.mfile.set_wl_ref(wllist)
    reader.mfile.qc_insitu.set_wllist_using_wlref(reader.mfile.wlref)

    # IN SITU QUALITY CONTROL
    reader.mfile.qc_insitu.check_indices_by_mu = True
    reader.mfile.qc_insitu.set_thershold(0, None, 0, 600)
    reader.mfile.qc_insitu.set_thershold(None, 0.005, 615, 625)
    reader.mfile.qc_insitu.set_thershold(None, 0.01, 410, 415)
    if name_mdb.find('CCI') > 0:
        reader.mfile.qc_insitu.set_thershold(None, 0.003, 650, 670)  ##CCI
        reader.mfile.qc_insitu.set_thershold(None, 0.004, 440, 450)  ##CCI
    reader.mfile.qc_insitu.apply_band_shift = True

    # SATELLITE QUALITY CONTROL
    reader.mfile.qc_sat.set_eumetsat_defaults(3)
    if 'satellite_pixel_classif_flags' in reader.mfile.nc.variables:
        idepix_flag = reader.mfile.nc.variables['satellite_pixel_classif_flags']
        reader.mfile.qc_sat.set_idepix_as_flag(idepix_flag)
    reader.mfile.qc_sat.add_band_statistics(-1, 400, 'avg', True, 0.003, 'greater')

    reader.mfile.prepare_df_validation()
    mplot = MDBPlot(reader.mfile, None)
    path_out = os.path.join(path_base, f'{name_mdb[:-3]}')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    mplot.make_validation_mdbfile(path_out)


def make_validation_from_dfvalid(ac, platform):
    # path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/CCI/MDBs/CCI_BalTower_Results_5'
    path_orig = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    path_base_b = os.path.join(path_orig, ac)
    if ac == 'CCIv6':
        name = f'MDB_{ac}'
    else:
        name = f'MDB_{ac}_S3{platform}'
    path_base = os.path.join(path_base_b, name)
    file_path = os.path.join(path_base, 'DataValid.csv')
    if not os.path.exists(file_path):
        print(f'File: {file_path} does not exist')
        return
    dfval = pd.read_csv(file_path, sep=';')
    wldata = dfval[dfval['Valid']]['Wavelenght']
    wllist = list(wldata.unique())
    mplot = MDBPlot(None, dfval)
    if ac == 'CCIv6':
        title = ac
        filenamebase = ac
    else:
        title = f'S3{platform} {ac}'
        filenamebase = f'S3{platform}_{ac}'
    print(path_base)
    mplot.make_validation_dfval(path_base, title, filenamebase, wllist)


def make_variations_qc_single_MDB():
    # path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/CCI/MDBs'
    path_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/OLCI_VEIT_UPDATED/MDBs'

    # name_mdb = 'MDB___1KM_CCI_L2_AERONET_Gustav_Dalen_Tower.nc'
    # name_mdb = 'MDB___1KM_CCI_L2_AERONET_Helsinki_Lighthouse.nc'
    name_mdb = 'MDB_S3A_OLCI_WFR_STANDARD_L2_HYPERNETS_VEIT.nc'
    path_mdb = os.path.join(path_base, name_mdb)
    reader = MDB_READER(path_mdb, False)

    # time window variation
    file_out = os.path.join(path_base, 'TimeEffect.csv')
    list_qc_ref = []
    qc_ref = {
        'name': '',
        'time_max': 0
    }
    timew = [30, 60, 90, 120, 150, 180]
    for t in timew:
        tseconds = t * 60
        qc_here = qc_ref.copy()
        qc_here['time_max'] = tseconds
        qc_here['name'] = f':{t} min.'
        list_qc_ref.append(qc_here)

    reader.compare_different_qcinsitu(list_qc_ref, file_out)


def make_bargraphics(df, xlabel, ylabel, title, fileout):
    h = plt.figure()
    plt.bar(df.index, df.loc[:, 'Total'], color=[0.8, 0.8, 0.8, 1])
    plt.bar(df.index, df.loc[:, 'Valid'], color=[0.65, 0.8, 0.55, 1])
    plt.xticks(df.index)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    str_legend = ['# Total match-ups', '# Valid match-ups']
    plt.legend(str_legend, loc='upper right')  # bbox_to_anchor=)
    plt.title(title)
    plt.savefig(fileout, dpi=300)
    plt.close(h)


def make_bargraphics_matchups(filein, fileout):
    df = pd.read_csv(filein, sep=';')
    print(df)
    h = plt.figure()
    plt.bar(df.loc[:, 'Month'], df.loc[:, 'Total'], color=[0.8, 0.8, 0.8, 1])
    plt.bar(df.loc[:, 'Month'], df.loc[:, 'Valid'], color=[0.65, 0.8, 0.55, 1])
    plt.xticks(df.loc[:, 'Month'])
    plt.xlabel('Month')
    plt.ylabel('# Match-ups')
    str_legend = ['# Total match-ups', '# Valid match-ups']
    plt.legend(str_legend, loc='upper right')  # bbox_to_anchor=)
    plt.title('Hypstar-OLCI WFR Match-ups (08-04-2021 to 08-05-2022)')
    plt.savefig(fileout, dpi=300)
    plt.close(h)


def make_graphics_param_ab(path_a, path_b, path_out):
    df_param_a = pd.read_csv(os.path.join(path_a, 'Params.csv'), sep=';')
    df_param_b = pd.read_csv(os.path.join(path_b, 'Params.csv'), sep=';')

    wl = list(df_param_a.columns)[3:]
    index_names = ['S3A', 'S3B']
    params = ['N', 'slope', 'intercept', 'rmse_val', 'mean_rel_diff', 'mean_abs_rel_diff', 'r2', 'bias']
    params_file = ['N', 'slope', 'intercept', 'RMSD', 'RPD', 'APD', 'r2', 'bias']
    for idx in range(len(params)):
        param = params[idx]
        param_file = params_file[idx]
        param_values_a = df_param_a[df_param_a['Param'] == param].iloc[0, 3:]
        param_values_b = df_param_b[df_param_b['Param'] == param].iloc[0, 3:]
        df_param = pd.DataFrame(index=index_names, columns=wl)
        df_param.loc['S3A'] = param_values_a
        df_param.loc['S3B'] = param_values_b
        df_param.to_csv(os.path.join(path_out, params_file[idx] + '.csv'), sep=';')
        df_param_new = get_df_from_pivot_df(df_param, 'platform', 'wavelength', param)

        h = plt.figure()
        sns.set_theme()
        dashes = [''] * len(index_names)
        markers = ['o'] * len(index_names)
        sns.relplot(kind='line', data=df_param_new, x="wavelength", y=param, hue="platform", markers=markers,
                    dashes=dashes,
                    style="platform")
        plt.savefig(os.path.join(path_out, param_file + '.jpg'), dpi=300)
        plt.close(h)


def make_togheter_ab(ac):
    file_path = f'/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/{ac}'
    prename = f'MDB_{ac}_S3'

    path_out = os.path.join(file_path, prename + 'AB')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    path_a = os.path.join(file_path, prename + 'A')
    path_b = os.path.join(file_path, prename + 'B')

    df_all_a = pd.read_csv(os.path.join(path_a, 'Data.csv'), sep=';')
    df_all_b = pd.read_csv(os.path.join(path_b, 'Data.csv'), sep=';')
    df_all = pd.concat([df_all_a, df_all_b], ignore_index=True)
    df_all.to_csv(os.path.join(path_out, 'Data.csv'), sep=';')

    df_valid_a = pd.read_csv(os.path.join(path_a, 'DataValid.csv'), sep=';')
    df_valid_b = pd.read_csv(os.path.join(path_b, 'DataValid.csv'), sep=';')
    df_valid = pd.concat([df_valid_a, df_valid_b], ignore_index=True)
    df_valid.to_csv(os.path.join(path_out, 'DataValid.csv'), sep=';')

    mplot = MDBPlot(None, df_valid)
    mplot.obtain_mu_info(df_all, df_valid, path_out)

    make_graphics_param_ab(path_a, path_b, path_out)

    h = plt.figure()
    sns.set_theme()
    sns.lmplot(data=df_valid, x="Ins_Rrs", y="Sat_Rrs", hue="platform", truncate=False, scatter_kws={"s": 10})
    xmin, xmax = plt.gca().get_xlim()
    ymin, ymax = plt.gca().get_ylim()
    xmin = np.min([xmin, ymin])
    xmax = np.min([xmax, ymax])
    if xmin < 0:
        xmin = 0
    plt.gca().set_xlim([xmin, xmax])
    plt.gca().set_ylim([xmin, xmax])
    plt.plot([xmin, xmax], [xmin, xmax], '--k')
    plt.savefig(os.path.join(path_out, prename + 'AB' + '.jpg'), dpi=300)
    plt.close(h)


def make_mu_info_bytower(path_base, ac, year_min, year_max):
    if not ac.startswith('CCI'):
        name_out = f'MUINFO_{ac}_AB'
        path_out = os.path.join(path_base, name_out)
    else:
        name_out = f'MUINFO_{ac}'
        path_out = os.path.join(path_base, name_out)
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    year = list(range(year_min, year_max + 1))
    year.reverse()
    towers_df = ['Gustav Dalen Tower', 'Irbe Lighthouse', 'Helsinki Lighthouse']
    towers = ['Gustav_Dalen_Tower', 'Irbe_Lighthouse', 'Helsinki_Lighthouse']
    col_names = ['Total', 'NoValid', 'Valid', '%Valid']
    dftower = pd.DataFrame(index=towers_df, columns=col_names, dtype=float)
    dfall = pd.DataFrame(index=year, columns=towers_df, dtype=float)
    dfvalid = pd.DataFrame(index=year, columns=towers_df, dtype=float)
    for year in range(year_min, year_max + 1):
        dfall.loc[year, :] = 0
        dfvalid.loc[year, :] = 0
    for tower in towers_df:
        dftower.loc[tower, :] = 0

    path_ac = os.path.join(path_base, ac)
    for tower in towers:
        tower_here = tower.replace('_', ' ')
        lista_paths = []
        for name in os.listdir(path_ac):
            path_mdb = os.path.join(path_ac, name)
            if os.path.isdir(path_mdb) and name.find(tower) >= 0 and (name.find('S3AB_') >= 0 or name.find('CCI') >= 0):
                lista_paths.append(path_mdb)
        for path_mdb in lista_paths:
            fall = os.path.join(path_mdb, 'Data.csv')
            fvalid = os.path.join(path_mdb, 'DataValid.csv')
            dfall_data = pd.read_csv(fall, sep=';')
            dfvalid_data = pd.read_csv(fvalid, sep=';')
            for index, row in dfall_data.iterrows():
                if row['Wavelenght'] == 490:
                    date_here = dt.strptime(row['Sat_Time'], '%Y-%m-%d %H:%M')
                    year_here = date_here.year
                    dfall.loc[year_here, tower_here] = dfall.loc[year_here, tower_here] + 1
                    dftower.loc[tower_here, 'Total'] = dftower.loc[tower_here, 'Total'] + 1
            for index, row in dfvalid_data.iterrows():
                if row['Wavelenght'] == 490:
                    date_here = dt.strptime(row['Sat_Time'], '%Y-%m-%d %H:%M')
                    year_here = date_here.year
                    dfvalid.loc[year_here, tower_here] = dfvalid.loc[year_here, tower_here] + 1
                    dftower.loc[tower_here, 'Valid'] = dftower.loc[tower_here, 'Valid'] + 1

    dfall.to_csv(os.path.join(path_out, f'{name_out}_ByYearTower_All.csv'), sep=';')
    dfvalid.to_csv(os.path.join(path_out, f'{name_out}_ByYearTower_Valid.csv'), sep=';')
    for tower in towers_df:
        dftower.loc[tower, 'NoValid'] = dftower.loc[tower, 'Total'] - dftower.loc[tower, 'Valid']
        dftower.loc[tower, '%Valid'] = (dftower.loc[tower, 'Valid'] / dftower.loc[tower, 'Total']) * 100
    dftower.to_csv(os.path.join(path_out, f'{name_out}_Tower.csv'), sep=';')

    title = f'# Potential match-ups ({ac})'
    plot_heatmap(dfall, 'Tower', 'Year', title, os.path.join(path_out, f'{name_out}_ByYearTower_All.jpg'))
    title = f'# Valid match-ups ({ac})'
    plot_heatmap(dfvalid, 'Tower', 'Year', title, os.path.join(path_out, f'{name_out}_ByYearTower_Valid.jpg'))

    title = f'# Match-ups ({ac})'
    make_bargraphics(dftower, 'Tower', '# Match-ups', title, os.path.join(path_out, f'{name_out}_Tower.jpg'))


def plot_heatmap(df, xlabel, ylabel, title, fjpg):
    h = plt.Figure()
    cmap = cm.get_cmap('RdYlBu_r')
    dfall_withnan = df
    dfall_withnan[df == 0] = np.nan
    sns.heatmap(dfall_withnan, annot=False, cmap=cmap, linewidths=1, linecolor='black')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(fjpg, dpi=300)
    plt.close(h)
    plt.close('all')


def make_mu_info(path_base, acnames, platform):
    # path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'

    for ac in acnames:
        postname = f'_S3{platform}'
        if ac == 'CCIv6':
            postname = ''
        file_valid = os.path.join(path_base, ac, 'MDB_' + ac + postname, 'DataValid.csv')
        file_all = os.path.join(path_base, ac, 'MDB_' + ac + postname, 'Data.csv')
        dfall = pd.read_csv(file_all, sep=';')
        dfvalid = pd.read_csv(file_valid, sep=';')
        mdbplot = MDBPlot(None, None)
        path_out = os.path.join(path_base, f'MUINFO_{ac}_{platform}')
        if not os.path.exists(path_out):
            os.mkdir((path_out))
        mdbplot.obtain_mu_info(dfall, dfvalid, path_out)


def make_together_atm(wlref, acnames, platform, name_out):
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    path_out = os.path.join(path_base, name_out)
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    df_wl = None
    name_acs = None

    for ac in acnames:
        postname = f'_S3{platform}'
        if ac == 'CCIv6':
            postname = ''

        file_ac = os.path.join(path_base, ac, 'MDB_' + ac + postname, 'DataValid.csv')
        print(file_ac)
        df_ac = pd.read_csv(file_ac, sep=';')
        wldata = df_ac[df_ac['Valid']]['Wavelenght']
        wllist = wldata.unique()
        index, wlref_f = get_index_wl_list(wlref, wllist)
        if index >= 0:
            print(wlref_f)
            dfwl_here = df_ac[(df_ac['Valid']) & (df_ac['Wavelenght'] == wlref_f)][:]
            if df_wl is None:
                df_wl = dfwl_here
            else:
                df_wl = pd.concat([df_wl, dfwl_here], ignore_index=True)

    df_wl.to_csv(os.path.join(path_out, 'WL_' + str(wlref) + '.csv'), sep=';')

    h = plt.figure()
    sns.set_theme()
    sns.lmplot(data=df_wl, x="Ins_Rrs", y="Sat_Rrs", hue="ac", truncate=False, scatter_kws={"s": 10})
    # plt.gca().set_aspect('equal', adjustable='box')
    xmin, xmax = plt.gca().get_xlim()
    ymin, ymax = plt.gca().get_ylim()
    xmin = np.min([xmin, ymin])
    xmax = np.min([xmax, ymax])
    if xmin < 0:
        xmin = 0
    plt.gca().set_xlim([xmin, xmax])
    plt.gca().set_ylim([xmin, xmax])
    plt.plot([xmin, xmax], [xmin, xmax], '--k')
    plt.savefig(os.path.join(path_out, 'WL_' + str(wlref) + '.jpg'), dpi=300)
    plt.close(h)


def make_together_atm_params(wlreflist, acnames, platform, name_out):
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    path_out = os.path.join(path_base, name_out)
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    params = ['N', 'slope_typeII', 'offset_typeII', 'rmse_val', 'mean_rel_diff', 'mean_abs_rel_diff', 'r2', 'bias',
              'r_value', 'XAVG', 'YAVG', 'CPRMSE', 'MAE']
    params_file = ['N', 'SLOPE', 'OFFSET', 'RMSE', 'RPD', 'APD', 'r2', 'BIAS', 'r', 'XAVG', 'YAVG', 'CPRMSE', 'MAE']

    includeCCI = False

    for idx in range(len(params)):
        param = params[idx]
        param_file = params_file[idx]
        df_param = pd.DataFrame(index=acnames, columns=wlreflist)
        name_acs = None
        for ac in acnames:
            postname = f'_S3{platform}'
            if ac == 'CCIv6':
                postname = ''
                includeCCI = True
            if name_acs is None:
                name_acs = ac
            else:
                name_acs = f'{name_acs}_{ac}'
            file_ac = os.path.join(path_base, ac, 'MDB_' + ac + postname, 'Params.csv')
            df_param_here = pd.read_csv(file_ac, sep=';')
            df_param_here = df_param_here[df_param_here['Param'] == param]
            wlparam = list(df_param_here.columns[3:].astype(dtype=np.float64))
            indices = get_indices_wl_list(wlreflist, wlparam)
            for index_wl in range(len(indices)):
                if indices[index_wl] >= 0:
                    index_df = indices[index_wl] + 3
                    # print(df_param_here.iat[0,index_df])
                    df_param.loc[ac].iat[index_wl] = df_param_here.iat[0, index_df]
        namefile = f'S3{platform}_{param_file}'
        if includeCCI:
            namefile = f'{name_acs}_{param_file}'
        df_param.to_csv(os.path.join(path_out, f'{namefile}.csv'), sep=';')
        df_param_new = get_df_from_pivot_df(df_param, 'ac', 'wavelength', param)

        df_param_new.rename(columns={param: param_file}, inplace=True)
        df_param_new.rename(columns={'wavelength': 'Wavelength (nm)'}, inplace=True)

        h = plt.figure()
        sns.set_theme()
        dashes = [''] * len(acnames)
        markers = ['o'] * len(acnames)
        sns.relplot(kind='line', data=df_param_new, x="Wavelength (nm)", y=param_file, hue="ac", markers=markers,
                    dashes=dashes,
                    style="ac")
        # plt.xlabel = 'Wavelength (nm)'
        plt.savefig(os.path.join(path_out, f'{namefile}.jpg'), dpi=300)
        plt.close(h)


def make_average_spectra(ac, platform, name_out):
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    path_out = os.path.join(path_base, name_out)
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    postname = f'_S3{platform}'
    if ac == 'CCIv6':
        postname = ''
    path_valid = os.path.join(path_base, ac, 'MDB_' + ac + postname, 'DataValid.csv')
    df_valid = pd.read_csv(path_valid, sep=';')

    # indices_mu = get_indicesmu_withcompletespectra_fromdfvalid(df_valid)

    ntotal = len(df_valid.index) * 2
    col_names = ['Wavelength', 'Rrs', 'Type', 'Valid']
    df_new = pd.DataFrame(columns=col_names, index=range(ntotal))

    nvalid = 0
    for idx in range(len(df_valid.index)):
        # index_mu_here = df_valid.iloc[idx].at['Index_MU']
        index_mu_here = df_valid.iloc[idx].index
        # valid_idx = False
        # if index_mu_here in indices_mu:
        #     nvalid = nvalid + 2
        #     valid_idx = True
        idx_insitu = idx
        idx_sat = idx + len(df_valid.index)
        df_new.loc[idx_insitu, 'Wavelength'] = df_valid.iloc[idx].at['Wavelenght']
        df_new.loc[idx_sat, 'Wavelength'] = df_valid.iloc[idx].at['Wavelenght']
        df_new.loc[idx_insitu, 'Type'] = 'In situ Rrs'
        df_new.loc[idx_sat, 'Type'] = 'Satellite Rrs'
        df_new.loc[idx_insitu, 'Rrs'] = df_valid.iloc[idx].at['Ins_Rrs']
        df_new.loc[idx_sat, 'Rrs'] = df_valid.iloc[idx].at['Sat_Rrs']

    # nmuvalid = nvalid/18
    # print('Antes',len(df_new.index))
    # df_new = df_new[df_new['Valid']==True][:]
    # print('Despues', len(df_new.index))
    # prevista = len(indices_mu)*18
    # print('Prevista', prevista, 'NValid: ', nvalid, 'Nmuvalid', nmuvalid)
    df_new.to_csv(os.path.join(path_out, f'Spectra_{ac}_S3{platform}.csv'), sep=';')

    h = plt.figure()
    sns.set_theme()
    sns.relplot(data=df_new, kind='line', x="Wavelength", y="Rrs", hue="Type", style="Type", markers=True, ci='sd')
    plt.title(f'S3{platform} {ac}')
    plt.savefig(os.path.join(path_out, f'Spectra_{ac}_S3{platform}.jpg'), dpi=300)
    plt.close(h)


def get_indicesmu_withcompletespectra_fromdfvalid(dfvalid):
    wllist = dfvalid.loc[:, 'Wavelenght']
    wllist_complete = list(np.unique(np.array(wllist)))

    indices_mu = []
    wl_here = []
    index_mu_ref = -1
    for index, row in dfvalid.iterrows():
        index_mu_here = row['Index_MU']
        wl_val_here = row['Wavelenght']
        if index_mu_here != index_mu_ref:
            if index_mu_ref != -1:
                if wl_here == wllist_complete:
                    indices_mu.append(index_mu_ref)
                    print(len(wl_here), index_mu_ref)
            wl_here = [wl_val_here]
        else:
            wl_here.append(wl_val_here)
        index_mu_ref = index_mu_here

    return indices_mu


def get_df_from_pivot_df(dforig, name_index, name_col, name_values):
    col_names = [name_index, name_col, name_values]
    nvalues = len(dforig.columns) * len(dforig.index)
    dfnew = pd.DataFrame(columns=col_names, index=range(nvalues))
    idx = 0
    for index in dforig.index:
        for col in dforig.columns:
            dfnew.iat[idx, 0] = index
            dfnew.iat[idx, 1] = np.float64(col)
            dfnew.iat[idx, 2] = np.float64(dforig.at[index, col])
            idx = idx + 1
    return dfnew


def get_indices_wl_list(wlreflist, wllist):
    indices = []
    for wlref in wlreflist:
        index, wlref_f = get_index_wl_list(wlref, wllist)
        indices.append(index)

    return indices


def get_index_wl_list(wlref, wllist):
    index = -1
    for idx in range(len(wllist)):
        wll = wllist[idx]
        dif = abs(wll - wlref)
        if dif <= 3:
            index = idx
    if index >= 0:
        wlref_f = wllist[index]
    else:
        wlref_f = None

    return index, wlref_f


def launch_complete_sam(ac):
    platforms = ['A', 'B']
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/'
    dfsam_fin = None
    for platform in platforms:
        fvalid = os.path.join(path_base, ac, f'MDB_{ac}_S3{platform}', 'DataValid.csv')
        dfsam, sam_values = compute_sam(fvalid, platform, ac)
        df_sam_check = pd.DataFrame(np.array(sam_values))
        if dfsam_fin is None:
            dfsam_fin = df_sam_check
        else:
            # dfsam_fin = pd.concat([dfsam_fin, dfsam], ignore_index=True)
            dfsam_fin = pd.concat([dfsam_fin, df_sam_check], axis=1)

    columns = []
    for idx in range(len(platforms)):
        p = platforms[idx]
        columns.append(f'{ac} S3{p}')

    dfsam_fin.columns = columns
    file_out = os.path.join(path_base, 'SAM', f'SAM{ac}.csv')

    # dfsam_fin.to_csv(file_out, sep=';')
    dfsam_fin.to_csv(file_out, sep=';')


def compute_sam(fvalid, platform, ac):
    # fvalid = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/POLYMER/MDB_POLYMER_S3B/DataValid.csv'
    dfvalid = pd.read_csv(fvalid, sep=';')
    wlpotential_lists = get_spectra_from_dfvalid(dfvalid)
    print(len(wlpotential_lists))
    wlcomplete = wlpotential_lists[0]
    col_names = ['AC', 'Platform', 'Complete', 'NBands']
    for idx in range(len(wlcomplete)):
        name = f'Band_{idx + 1}'
        col_names.append(name)
    col_names.append('NMu')
    col_names.append('SAM_sum')
    col_names.append('SAM_value')
    df_sam = pd.DataFrame(columns=col_names, index=range(len(wlpotential_lists)))
    for idx in range(len(df_sam.index)):
        complete_value = 0
        if idx == 0:
            complete_value = 1
        df_sam.loc[idx].at['AC'] = ac
        df_sam.loc[idx].at['Platform'] = platform
        df_sam.loc[idx].at['Complete'] = complete_value
        df_sam.loc[idx].at['NMu'] = 0
        df_sam.loc[idx].at['SAM_sum'] = 0
        df_sam.loc[idx].at['SAM_value'] = 0
        wl_list_here = wlpotential_lists[idx]
        df_sam.loc[idx].at['NBands'] = len(wl_list_here)
        for iband in range(len(wl_list_here)):
            name = f'Band_{iband + 1}'
            df_sam.loc[idx].at[name] = wl_list_here[iband]

    df_sam, sam_values = compute_sam_impl(dfvalid, df_sam, wlpotential_lists)

    return df_sam, sam_values


def get_spectra_from_dfvalid(dfvalid):
    wllist = dfvalid.loc[:, 'Wavelenght']
    wllist_complete = list(np.unique(np.array(wllist)))

    wlpotential_list = [wllist_complete]

    wl_here = []
    index_mu_ref = -1
    for index, row in dfvalid.iterrows():
        index_mu_here = row['Index_MU']
        wl_val_here = row['Wavelenght']
        if index_mu_here != index_mu_ref:
            if index_mu_ref != -1:
                exists = False
                for wlexisting in wlpotential_list:
                    if wlexisting == wl_here:
                        exists = True
                if not exists:
                    wlpotential_list.append(wl_here)
            wl_here = [wl_val_here]
        else:
            wl_here.append(wl_val_here)
        index_mu_ref = index_mu_here

    for wlexisting in wlpotential_list:
        print(wlexisting)
    return wlpotential_list


def compute_sam_impl(dfvalid, dfsam, wl_potentiallist):
    wl_here = []
    insitu_here = []
    sat_here = []
    index_mu_ref = -1
    sam_values = []
    for index, row in dfvalid.iterrows():
        index_mu_here = row['Index_MU']
        insitu_val_here = row['Ins_Rrs']
        sat_val_here = row['Sat_Rrs']
        wl_val_here = row['Wavelenght']
        if index_mu_here != index_mu_ref:
            if index_mu_ref != -1:
                Topline = np.dot(np.array(insitu_here),
                                 np.array(sat_here))  # TOTAL(double(class_rrs) * double(input_rrs))
                Botline1 = np.sqrt(np.sum(np.array(insitu_here) * np.array(insitu_here)))
                Botline2 = np.sqrt(np.sum(np.array(sat_here) * np.array(sat_here)))
                val_mu = math.acos(Topline / (Botline1 * Botline2))
                index_dfsam = -1
                for idx in range(len(wl_potentiallist)):
                    wl_potential = wl_potentiallist[idx]
                    if wl_potential == wl_here:
                        index_dfsam = idx
                if index_dfsam >= 0:
                    dfsam.loc[index_dfsam].at['NMu'] = dfsam.loc[index_dfsam].at['NMu'] + 1
                    dfsam.loc[index_dfsam].at['SAM_sum'] = dfsam.loc[index_dfsam].at['SAM_sum'] + val_mu
                    sam_values.append(val_mu)
            wl_here = [wl_val_here]
            insitu_here = [insitu_val_here]
            sat_here = [sat_val_here]
        else:
            wl_here.append(wl_val_here)
            insitu_here.append(insitu_val_here)
            sat_here.append(sat_val_here)
        index_mu_ref = index_mu_here

    for idx in range(len(dfsam.index)):
        dfsam.loc[idx].at['SAM_value'] = dfsam.loc[idx].at['SAM_sum'] / dfsam.loc[idx].at['NMu']

    # print(dfsam)
    return dfsam, sam_values


def compute_sam_impl_deprecated(dfvalid, dfsam, wl_potentiallist):
    wl_here = []
    insitu_here = []
    sat_here = []
    index_mu_ref = -1
    for index, row in dfvalid.iterrows():
        index_mu_here = row['Index_MU']
        insitu_val_here = row['Ins_Rrs']
        sat_val_here = row['Sat_Rrs']
        wl_val_here = row['Wavelenght']
        if index_mu_here != index_mu_ref:
            if index_mu_ref != -1:
                dot_p = np.dot(np.array(insitu_here), np.array(sat_here))
                nor_p = np.linalg.norm(insitu_here) * np.linalg.norm(sat_here)
                val_mu = math.acos(dot_p / nor_p)
                index_dfsam = -1
                for idx in range(len(wl_potentiallist)):
                    wl_potential = wl_potentiallist[idx]
                    if wl_potential == wl_here:
                        index_dfsam = idx
                if index_dfsam >= 0:
                    dfsam.loc[index_dfsam].at['NMu'] = dfsam.loc[index_dfsam].at['NMu'] + 1
                    dfsam.loc[index_dfsam].at['SAM_sum'] = dfsam.loc[index_dfsam].at['SAM_sum'] + val_mu
            wl_here = [wl_val_here]
            insitu_here = [insitu_val_here]
            sat_here = [sat_val_here]
        else:
            wl_here.append(wl_val_here)
            insitu_here.append(insitu_val_here)
            sat_here.append(sat_val_here)
        index_mu_ref = index_mu_here

    for idx in range(len(dfsam.index)):
        dfsam.loc[idx].at['SAM_value'] = dfsam.loc[idx].at['SAM_sum'] / dfsam.loc[idx].at['NMu']

    print(dfsam)
    return dfsam


def graficos():
    fbase = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/PRESENTATION_02_05'

    # aeronet bands
    file_out = os.path.join(fbase, 'AeronetBands.jpg')
    fbands = os.path.join(fbase, 'AeronetBands.csv')
    dfbands = pd.read_csv(fbands, sep=';', header=0, index_col=0)
    h = plt.figure()
    colormap = [
        [1, 1, 1, 1],
        [0.8, 0.8, 0.8],
        [0.8, 0.8, 0.8],
        [0.8, 0.8, 0.8],
        [0.8, 0.8, 0.8],
        [0.8, 0.8, 0.8],
        [0.8, 0.8, 0.8],
        [0.8, 0.8, 0.8]
    ]
    sns.heatmap(dfbands, annot=False, cmap=colormap, linewidths=1, linecolor='black', cbar=False)
    plt.xlabel('Wavelength (nm)')
    plt.savefig(file_out, dpi=300)
    plt.close(h)
    file_out = os.path.join(fbase, 'AeronetStations.jpg')
    fbands = os.path.join(fbase, 'AeronetStations.csv')
    dfbands = pd.read_csv(fbands, sep=';', header=0, index_col=0)
    plt.figure()
    sns.heatmap(dfbands, annot=False, cmap=colormap, linewidths=1, linecolor='black', cbar=False)
    # plt.xlabel('Wavelength (nm)')
    plt.savefig(file_out, dpi=300)
    plt.close()


def test():
    print('test')
    # ngood = 0
    # for imu in range(reader.mfile.qc_insitu.nmu):
    #     print('imu',imu)
    #     ins_time_index, mu_insitu_time, time_condition, spectrum_complete, rrs_values =reader.mfile.retrieve_ins_info_mu_spectra(imu)
    #     if time_condition and spectrum_complete:
    #         ngood = ngood + 1
    #
    # print(ngood)

    # print(reader.mfile.qc_insitu.mu_spectra_complete[imu],reader.mfile.qc_insitu.wl_indices[imu],reader.mfile.qc_insitu.wl_list[imu])

    # reader.mfile.qc_insitu.apply_band_shift = True
    # reader.mfile.qc_sat.set_eumetsat_defaults(3)
    # reader.mfile.prepare_df_validation()
    # mplot = MDBPlot(reader.mfile, None)
    # path_out = os.path.join(path_base, f'{name_mdb[:-3]}')
    # if not os.path.exists(path_out):
    #     os.mkdir(path_out)
    # mplot.make_validation_mdbfile(path_out)

    # print(reader.mfile.qc_insitu.wl_list)
    # print(reader.mfile.qc_insitu.wl_indices)
    # spectra = reader.mfile.qc_insitu.get_all_good_spectra()
    # print(spectra.shape)
    # print(reader.mfile.qc_insitu.wl_list)
    # print(reader.mfile.qc_insitu.nbands, reader.mfile.qc_insitu.nmu)
    # for i in range(spectra.shape[1]):
    #     print(i,spectra[:,i].count())

    # reader.mfile.prepare_df_validation()
    # mplot = MDBPlot(reader.mfile, None)
    # path_out = os.path.join(path_base, name_mdb[:-3])
    # if not os.path.exists(path_out):
    #     os.mkdir(path_out)
    # mplot.make_validation_mdbfile(path_out)

    # reader.mfile.qc_sat.compute_invalid_masks(0)
    # reader.mfile.qc_sat.compute_flag_stats(5)
    # reader.mfile.qc_sat.compute_flag_masks(5)
    # reader.mfile.qc_sat.add_theshold_mask(0, -1, 0.001, 'greater')
    # reader.mfile.qc_sat.add_band_statistics(-1, 560,'CV', 20, 'greater')
    # b = reader.mfile.qc_sat.compute_statistics(0)
    # print(b)
    # b = reader.mfile.qc_sat.do_check_statistics()
    # print(b)


if __name__ == '__main__':
    main()
