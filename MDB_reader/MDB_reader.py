from MDBFile import MDBFile
from MDBInSituFile import MDBInSituFile
from MDBPlot import MDBPlot
from MDBFileList import MDBFileList
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from PlotSpectra import PlotSpectra
import numpy as np
import math


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
            #print(mfile_here.qc_insitu.name, mfile_here.qc_insitu.time_max)
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

    do_lps()

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

def do_lps():
    dir_base = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/OLCI_VEIT_UPDATED/MDBs_20052022'
    dir_ab = os.path.join(dir_base, 'MDB_S3AB')

    # Single validation (scaterplot for each platform)
    file_s3a = os.path.join(dir_base, 'MDB_S3A_OLCI_WFR_STANDARD_L2_HYPERNETS_VEIT.nc')
    file_s3b = os.path.join(dir_base, 'MDB_S3B_OLCI_WFR_STANDARD_L2_HYPERNETS_VEIT.nc')
    #make_validation_single_MDB(dir_base,file_s3a)
    #make_validation_single_MDB(dir_base, file_s3b)

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
    path_a = os.path.join(dir_base,'MDB_S3A_OLCI_WFR_STANDARD_L2_HYPERNETS_VEIT')
    path_b = os.path.join(dir_base, 'MDB_S3B_OLCI_WFR_STANDARD_L2_HYPERNETS_VEIT')
    make_graphics_param_ab(path_a, path_b, dir_ab)

    filein = os.path.join(dir_ab, 'MUByMonth.csv')
    fileout = os.path.join(dir_ab, 'MUByMonth.jpg')
    make_bargraphics_matchups(filein, fileout)




def do_chla():
    file_mdb = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/MDBs/MDB_S3A_B_OLCI_POLYMER_INSITU_20160401_20220531.nc'
    baltic_mlp_code = '/home/lois/PycharmProjects/aeronet'
    mfile = MDBInSituFile(file_mdb)
    mfile.start_baltic_chla(baltic_mlp_code, 'insitu_CHLA')
    mfile.qc_sat.max_diff_wl = 6
    mfile.set_wlsatlist_aswlref([412, 443, 490, 510, 555, 670])

    # for imu in range(mfile.n_mu_total):
    #     mfile.compute_baltic_chla_ensemble_mu(imu)
    mfile.prepare_df_validation()
    print(mfile.df_validation_valid)
    file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/CHLA/MDBs/ChlaBalMatchUps.csv'
    mfile.df_validation_valid.to_csv(file_out, sep=';')


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


def get_table_stats(wlref):
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/SYKE'
    path_out = os.path.join(path_base, 'TABLE_PARAM')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    acs = ['STANDARD', 'POLYMER', 'C2RCC', 'FUB']
    # acs = ['POLYMER', 'CCI']
    params = ['N', 'r2', 'RMSD', 'APD', 'RPD', 'bias']
    table = pd.DataFrame(columns=acs, index=params)
    for ac in acs:
        postname = '_S3A'
        if ac == 'CCI':
            postname = ''
        path_ac = os.path.join(path_base, ac, f'MDB_{ac}{postname}', 'Params.csv')
        table_ac = pd.read_csv(path_ac, sep=';')
        index = -1
        if wlref == -1:
            index = 2
            name_out = 'ParamAll_SYKE.csv'
        else:

            wllist = list(table_ac.columns[3:].astype(dtype=np.float64))
            index_wl, wl_aca = get_index_wl_list(wlref, wllist)
            if index_wl >= 0:
                index = index_wl + 3
                name_out = f'Param{wlref}.csv'
        if index >= 0:
            table.loc['N', ac] = table_ac.loc[0].iat[index]
            table.loc['r2', ac] = table_ac.loc[11].iat[index]
            table.loc['RMSD', ac] = table_ac.loc[6].iat[index]
            table.loc['APD', ac] = table_ac.loc[8].iat[index]
            table.loc['RPD', ac] = table_ac.loc[7].iat[index]
            table.loc['bias', ac] = table_ac.loc[9].iat[index]

            file_out = os.path.join(path_out, name_out)
            table.to_csv(file_out, sep=';')


def make_validation_list_MDB():
    mdblist = MDBFileList()
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/C2RCC'
    wcard = 'S3A'
    for f in os.listdir(path_base):
        if (wcard is not None and f.find(wcard) > 0) and f.endswith('.nc'):
            path_mdb = os.path.join(path_base, f)
            mdblist.add_mdb_file(path_mdb)

    path_out = os.path.join(path_base, 'MDB_C2RCC_S3A')
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    wllist = [400, 412, 443, 490, 510, 560, 620, 665, 779]  # -> POLYMER, STANDARD, FUB (NO INCLUIRIA 412)
    mdblist.set_wlsatlist_from_wlreflist_asref(wllist)
    print(mdblist.wlref)

    mdblist.qc_base.set_sat_eumetsat_defaults(3)
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
    reader.mfile.set_wlsatrange_aswlref(400, 800)
    reader.mfile.qc_insitu.set_wllist_using_wlref(reader.mfile.wlref)
    reader.mfile.qc_insitu.check_indices_by_mu = False
    reader.mfile.qc_insitu.set_thershold(0, None, 0, 600)
    reader.mfile.qc_insitu.set_thershold(None, 0.006, 615, 625)
    reader.mfile.qc_insitu.set_thershold(None, 0.01, 410, 415)
    # reader.mfile.qc_insitu.apply_band_shift = True
    reader.mfile.qc_sat.set_eumetsat_defaults(3)
    reader.mfile.prepare_df_validation()
    mplot = MDBPlot(reader.mfile, None)
    path_out = os.path.join(path_base, f'{name_mdb[:-3]}')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    mplot.make_validation_mdbfile(path_out)


def make_validation_from_dfvalid(ac, platform):
    # path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/CCI/MDBs/CCI_BalTower_Results_5'
    path_base_b = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/SYKE/'
    name = f'MDB_{ac}_S3{platform}'
    path_base = os.path.join(path_base_b, name)
    file_path = os.path.join(path_base, 'DataValid.csv')
    dfval = pd.read_csv(file_path, sep=';')
    wldata = dfval[dfval['Valid']]['Wavelenght']
    wllist = list(wldata.unique())
    mplot = MDBPlot(None, dfval)
    mplot.make_validation_dfval(path_base, ac, ac, wllist)


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


def make_together_atm(wlref, acnames, name_out):
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/SYKE'
    path_out = os.path.join(path_base, name_out)
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    df_wl = None
    for ac in acnames:
        postname = '_S3A'
        if ac == 'CCI':
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


def make_together_atm_params(wlreflist, acnames, name_out):
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs/SYKE'
    path_out = os.path.join(path_base, name_out)
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    params = ['N', 'slope', 'intercept', 'rmse_val', 'mean_rel_diff', 'mean_abs_rel_diff', 'r2', 'bias']
    params_file = ['N', 'slope', 'intercept', 'RMSD', 'RPD', 'APD', 'r2', 'bias']
    for idx in range(len(params)):
        param = params[idx]
        param_file = params_file[idx]
        df_param = pd.DataFrame(index=acnames, columns=wlreflist)
        for ac in acnames:
            postname = '_S3A'
            if ac == 'CCI':
                postname = ''
            file_ac = os.path.join(path_base, ac, 'MDB_' + ac + postname, 'Params.csv')
            df_param_here = pd.read_csv(file_ac, sep=';')
            df_param_here = df_param_here[df_param_here['Param'] == param]
            wlparam = list(df_param_here.columns[3:].astype(dtype=np.float64))
            indices = get_indices_wl_list(wlreflist, wlparam)
            print(ac, indices)
            for index_wl in range(len(indices)):
                if indices[index_wl] >= 0:
                    index_df = indices[index_wl] + 3
                    # print(df_param_here.iat[0,index_df])
                    df_param.loc[ac].iat[index_wl] = df_param_here.iat[0, index_df]
        df_param.to_csv(os.path.join(path_out, param_file + '.csv'), sep=';')
        df_param_new = get_df_from_pivot_df(df_param, 'ac', 'wavelength', param)
        h = plt.figure()
        sns.set_theme()
        dashes = [''] * len(acnames)
        markers = ['o'] * len(acnames)
        sns.relplot(kind='line', data=df_param_new, x="wavelength", y=param, hue="ac", markers=markers, dashes=dashes,
                    style="ac")
        plt.savefig(os.path.join(path_out, param_file + '.jpg'), dpi=300)
        plt.close(h)


def make_average_spectra(ac, platform):
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    path_out = os.path.join(path_base, 'SPECTRA')
    if not os.path.exists(path_out):
        os.mkdir(path_out)
    postname = f'_S3{platform}'
    if ac == 'CCI':
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
