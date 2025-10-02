import argparse
import os.path
import shutil, os
from datetime import timedelta
from datetime import datetime as dt
import numpy as np
from netCDF4 import Dataset
import warnings
import COMMON.args_functions as arf
import zipfile as zp

warnings.filterwarnings("ignore", category=UserWarning)

parser = argparse.ArgumentParser(
    description="MDB Utilities")

parser.add_argument('-m', "--mode", help='Mode option',
                    choices=["unzip_s3","insitu_brdf","add_instrument_id", "hypstar_check", "correct_neg_values", "TEST"],
                    required=True)
parser.add_argument('-i', "--input_path", help="Input path.")
parser.add_argument('-o', "--output", help="Output file or path")
parser.add_argument('-s', "--source_path", help="Source path.", default="/dst04-data1/OC/OLCI/daily_v202311_bc")
parser.add_argument('-p', "--param", help="Param for TEST")
parser.add_argument('-sd', "--start_date", help="The Start Date - format YYYY-MM-DD ")
parser.add_argument('-ed', "--end_date", help="The End Date - format YYYY-MM-DD ")
args = parser.parse_args()


def main():
    print(f'[INFO] Started MDB_utils!')
    if args.mode == 'unzip_s3':
        if not args.input_path:
            print(f'[ERROR] Input path is required')
            return
        if not os.path.isdir(args.input_path):
            print(f'[ERROR] Input path is not available or is not a valid directory')
            return
        start_date,end_date = arf.get_start_end_date_from_args(args)
        if start_date is None:
            print(f'[ERROR] -sd (--start_date) and -ed (--end_date) are required')
            return
        output_path = args.input_path if not args.output else args.output
        run_unzip_s3(input_path,output_path,start_date,end_date)

    if args.mode == 'insitu_brdf':
        if not check_required_params(['input_path']):
            return
        if not os.path.isfile(args.input_path):
            print(f'[ERROR] {args.input_path} does not exist or is not a valid file')
            return
        output_path = args.output if args.output else None
        run_insitu_brdf()

    if args.mode == 'add_instrument_id':
        if not check_required_params(['input_path']):
            return
        if not os.path.isfile(args.input_path):
            print(f'[ERROR] {args.input_path} does not exist or is not a valid file')
            return
        output_path = args.output if args.output else None
        add_instrument_id(args.input_path, output_path)

    if args.mode == 'correct_neg_values':
        if not check_required_params(['input_path', 'output']):
            return
        if os.path.isdir(args.input_path) and os.path.isdir(args.output):
            correct_negative_values_from_extracts(args.input_path, args.output)
        if os.path.isfile(args.input_path):
            correct_negative_values_from_mdb(args.input_path, args.output)
        else:
            if not os.path.exists(args.input_path):
                print(f'[ERROR]{args.input_path} does not exist. It should be a valid file or directory')
                return
            if not os.path.exists(args.output) and os.path.isdir(args.input_path):
                print(f'[ERROR]{args.output} does not exist. It should be a valid file or directory')
                return

    if args.mode == 'hypstar_check':
        run_hypstar_check()

    if args.mode == 'TEST':
        # year = -1
        # if args.param:
        #     year = int(args.param)
        # run_test(year)
        #make_test()
        site = 'TRIT'
        path_base = f'/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/ITALIAN_SITES_VALIDATION_PUBLICATION/OCI/{site}'
        file_in = os.path.join(path_base,f'MDB_rc_PACE_OCI_1KM_HYPSTAR_{site}_COMMONMU.nc')
        file_out  = os.path.join(path_base,f'MDB_rc_PACE_OCI_1KM_HYPSTAR_{site}_COMMONMU_V3WL.nc')
        wl_list = [591,593,596,598,601,603,605,608,610]
        remove_wl_from_mu_variables(file_in,file_out,wl_list)


def run_unzip_s3(input_path,output_path,start_date,end_date):
    work_date = start_date
    while work_date<=end_date:
        input_path_date = os.path.join(input_path,work_date.strftime('%Y'),work_date.strftime('%j'))
        output_path_date = os.path.join(output_path,work_date.strftime('%Y'),work_date.strftime('%j'))
        if not os.path.isdir(output_path_date):
            if not create_output_path_date(output_path,work_date):
                print(f'[WARNING] Output path for date {work_date.strftime("%Y-%m-%d")}: {output_path_date} is not available and could not be created. Please review permisssions.')
                work_date = work_date + timedelta(days=1)
                continue
        if os.path.isdir(input_path_date):
            for name in os.listdir(input_path_date):
                if name.startswith('S3') and name.endswith('.zip'):
                    file_zip = os.path.join(input_path_date,name)
                    if not zp.is_zipfile(file_zip):
                        print(f'[WARNING] File {file_zip} is not a valid zip file. Skipping...')
                        continue
                    name_base = name[0:-4]
                    path_prod_u = os.path.join(output_path_date, name_base)
                    with zp.ZipFile(file_zip, 'r') as zprod:
                        print(f'[INFO] Unziping {name_base} to {output_path_date}')
                        zprod.extractall(path=output_path_date)

                    if os.path.isdir(path_prod_u):
                        os.remove(file_zip)

        work_date = work_date + timedelta(days=1)

def create_output_path_date(output_path,work_date):
    output_path_year = os.path.join(output_path,work_date.strftime('%Y'))
    if not os.path.isdir(output_path_year):
        try:
            os.mkdir(output_path_year)
        except:
            return False
    output_path_date = os.path.join(output_path_year,work_date.strftime('%j'))
    if not os.path.isdir(output_path_date):
        try:
            os.mkdir(output_path_date)
        except:
            return False
    return True


def run_insitu_brdf():
    import __init__
    import sys
    from netCDF4 import Dataset
    dset = Dataset(args.input_path)
    code_home = os.path.dirname(__init__.__file__)
    code_025 = os.path.join(os.path.dirname(os.path.dirname(code_home)),'O25')

    sys.path.append(code_025)
    print(code_025,os.path.exists(code_025))
    from o25 import O25_hyp
    # Define inputs
    l = dset.variables['insitu_original_bands'][1,:]  # Wavelengths in nm
    Rrs = dset.variables['insitu_Rrs'][5,:,0] # Remote sensing reflectance
    print(l.shape,Rrs.shape)
    sza = dset.variables['insitu_solar_zenith_angle'][5,0]
    vza = dset.variables['insitu_viewing_zenith_angle'][5,0]
    raa = dset.variables['insitu_viewing_azimuth_angle'][5,0]-dset.variables['insitu_solar_azimuth_angle'][5,0]
    print(sza,vza,raa)
    geom_old = [sza,vza, raa]
    geom_new = [0,0,0]
    #geom_old = [...]  # Original geometry: [solar zenith angle, viewing zenith angle, relative azimuth angle]
    #geom_new = [...]  # New geometry
    # Run O25
    a, bb, Rrs_N = O25_hyp.O25_hyp(l, Rrs, geom_old, geom_new)
    fout = '/mnt/c/DATA/FICE2025/BRDF_TEST/test_brdf.csv'
    fw = open(fout,'w')
    fw.write('WL;RRS_In;RRS_out')
    for idx in range(len(Rrs)):
        fw.write('\n')
        fw.write(f'{l[idx]};{Rrs[idx]};{Rrs_N[idx]}')
    fw.close()


def remove_wl_from_mu_variables(input_file,output_file,wl_list):
    from netCDF4 import Dataset
    input_dataset = Dataset(input_file)
    ncout = Dataset(output_file, 'w', format='NETCDF4')
    wl_array = input_dataset.variables['mu_wavelength'][:]
    data_include = np.ones(wl_array.shape)
    for wl in wl_list:
        data_include[wl_array==wl]=0
    nmu_new = np.sum(data_include)
    print('New number of mu: ',nmu_new)
    variables_to_modify = ['mu_ins_rrs','mu_ins_rrs_unc','mu_sat_rrs','mu_sat_rrs_unc','mu_satellite_id','mu_wavelength']

    # copy global attributes all at once via dictionary
    ncout.setncatts(input_dataset.__dict__)

    #dimensions
    for name, dimension in input_dataset.dimensions.items():
        if name=='mu_id':
            ncout.createDimension(name,nmu_new)
        else:
            ncout.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

    for name, variable in input_dataset.variables.items():
        print('======================>', name)
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             shuffle=True, complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(input_dataset[name].__dict__)

        if name in variables_to_modify:
            array_old = input_dataset[name][:]
            array_new = array_old[data_include==1]
            ncout[name][:] = array_new[:]
        else:
            ncout[name][:] = input_dataset[name][:]

    input_dataset.close()
    ncout.close()
print('COMPLETED')

def make_test():
    print('STARTED')
    import sys
    sys.path.append('/home/lois/PycharmProjects/pyOWT')
    from pyowt.OWT import OWT
    from pyowt.OpticalVariables import OpticalVariables
    from matplotlib import pyplot as plt
    from matplotlib.colors import Normalize
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    file_nc = '/mnt/c/DATA_LUIS/TARA_WORK/MDBs/EUMETSAT_L2/wide/MDBrc_S3AB_OLCI_WFR_STANDARD_20230101T000000_20231231T235959_HYPERBOOST_wide.nc'
    dir_output = '/mnt/c/DATA_LUIS/TARA_WORK'
    from netCDF4 import Dataset
    dataset = Dataset(file_nc, 'r')
    insitu_wl = dataset.variables['insitu_original_bands'][:]
    insitu_rrs = np.ma.squeeze(dataset.variables['insitu_Rrs'][:])
    insitu_latitude = np.ma.squeeze(dataset.variables['insitu_latitude'][:])
    insitu_longitude = np.ma.squeeze(dataset.variables['insitu_longitude'][:])


    ##SATELLITE
    # satellite_wl = dataset.variables['satellite_bands'][0:14]
    # satellite_Rrs = dataset.variables['satellite_Rrs'][:,0:14,:,:]
    # satellite_Rrs = np.ma.transpose(satellite_Rrs,(0,2,3,1))
    # satellite_lat = dataset.variables['satellite_latitude'][:]
    # satellite_lon = dataset.variables['satellite_longitude'][:]
    # n_all_data = satellite_Rrs.shape[0]*satellite_Rrs.shape[1]*satellite_Rrs.shape[2]
    # n_days = satellite_Rrs.shape[0]
    # satellite_Rrs = satellite_Rrs.reshape((n_all_data,14))

    dataset.close()

    ov_insitu = OpticalVariables(Rrs=insitu_rrs, band=insitu_wl, sensor=None)
    ##feed data into classification
    owt = OWT(AVW=ov_insitu.AVW, Area=ov_insitu.Area, NDI=ov_insitu.NDI)
    owt_classes = {
        -1: 'NaN',
        0: '1',
        1: '2',
        2: '3a',
        3: '3b',
        4: '4a',
        5: '4b',
        6: '5a',
        7: '5b',
        8: '6',
        9: '7'
    }
    owt_res = np.squeeze(np.array(list(owt.type_idx)))

    print(owt_res.shape[0])


    color_OWT = owt.classInfo.typeColHex


    ##SPECTRA

    from PlotSpectra import PlotSpectra

    pspectra = PlotSpectra()
    style = pspectra.line_style_default
    pspectra.xdata = insitu_wl

    owt_make = 'Mix'

    for owt in owt_classes:
        if owt==-1:
            indices = np.where(owt_res == owt)
            rrs_here = np.ma.squeeze(insitu_rrs[indices, :]).transpose()
            print(owt, rrs_here.shape)
            continue
        # if owt!=owt_make:
        #     continue

        indices = np.where(owt_res == owt)
        color = color_OWT[owt] if owt >= 0 else 'gray'
        rrs_here = np.ma.squeeze(insitu_rrs[indices,:]).transpose()
        print(owt,rrs_here.shape)

        style_here =style.copy()
        style_here['color'] = color
        pspectra.plot_data(rrs_here,style_here)

    pspectra.set_xaxis_title('Wavelength(nm)')
    xlabel = r'R$_r$$_s$ (sr$^-$$^1$)'
    pspectra.set_yaxis_title(xlabel)
    pspectra.set_grid_horizontal()
    #pspectra.set_y_range(-0.025,0.175)
    pspectra.set_tigth_layout()

    pspectra.save_fig(os.path.join(dir_output,f'All_Spectra_{owt_make}.tif'))

    ###MAP
    # geo_limits = [35, 65, -15, 30]
    # extent = (geo_limits[2], geo_limits[3], geo_limits[0], geo_limits[1])
    # ax = plt.axes(projection=ccrs.PlateCarree(), extent=extent)
    # ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black', linewidth=0.5)
    # gl = ax.gridlines(linewidth=0.5, linestyle='dotted', draw_labels=True)
    # gl.xlabels_top = False
    # gl.ylabels_right = False
    # handles = []
    # for owt in owt_classes:
    #     indices = np.where(owt_res == owt)
    #     color = color_OWT[owt] if owt >= 0 else 'gray'
    #     lat_points = insitu_latitude[indices]
    #     lon_points = insitu_longitude[indices]
    #     h = plt.plot(lon_points.tolist(), lat_points.tolist(), color=color, marker='o',markersize=3, linestyle='-',linewidth=0, transform=ccrs.PlateCarree())
    #     handles.append(h[0])
    # str_legend = [f'OWT_{owt_classes[owt]}' for owt in owt_classes]
    # plt.legend(handles, str_legend, framealpha=1, markerscale=1)
    # plt.tight_layout()
    # file_out = os.path.join(dir_output, 'OWT_Map_legend.tif')
    # plt.savefig(file_out, dpi=300, bbox_inches='tight', pil_kwargs={"compression": "tiff_lzw"})

    # for idx in range(owt_res):
    #     lat_p = insitu_latitude[idx]
    #     lon_p = insitu_longitude[idx]

    # all_owt = []
    # for idx,o_res in enumerate(owt_res):
    #     all_owt.append(o_res[0])
    #     print(ov_insitu.AVW[idx])
    # all_owt = np.unique(np.array(all_owt))
    # print(all_owt)

    ##SATELLITE
    # ov_sat = OpticalVariables(Rrs=satellite_Rrs,band=satellite_wl,sensor='olci-s3a')
    # avw = np.squeeze(np.array(list(ov_sat.AVW)))
    # avw_array = np.ma.array(avw.reshape((n_days,25,25)))
    # avw_array[avw_array<400] = np.ma.masked
    # avw_array[avw_array>700] = np.ma.masked
    # print(avw_array.shape,type(avw_array))
    # avw_array = np.ma.masked_where(avw_array,np.logical_or(avw_array<400,avw_array>700))

    # land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face',
    #                                         facecolor=cfeature.COLORS['land'])
    # for iday in range(n_days):
    #     data_array = np.ma.squeeze(avw_array[iday,:,:])
    #     if np.ma.count(data_array)>0:
    #         lat_array = satellite_lat[iday, :, :]
    #         lon_array = satellite_lon[iday, :, :]
    #         extent = (np.min(lon_array[:]),np.max(lon_array[:]), np.min(lat_array[:]),np.max(lat_array[:]))
    #         projection = ccrs.Mercator()
    #         ax = plt.axes(projection=projection, extent=extent)
    #         print(iday,np.ma.min(data_array),np.ma.max(data_array))
    #         plt.pcolormesh(lon_array, lat_array, data_array, transform=ccrs.PlateCarree(), cmap='jet',
    #                        norm=Normalize(vmin=400, vmax=700))
    #         plt.colorbar(shrink=0.6)
    #         ax.add_feature(land_10m, zorder=0, edgecolor='black', linewidth=0.5)
    #         plt.tight_layout()
    #         file_out = os.path.join(dir_output,f'Image_{iday}.png')
    #         plt.savefig(file_out)
    #         plt.close()

    # for idx in range(len(ov_sat.AVW)):
    #     print(satellite_Rrs[idx,1],ov_sat.AVW[idx][0])

    # show classification results
    ##print(owt.type_str)

    # file_original = '/mnt/c/DATA_LUIS/OCTAC_WORK/INC_NEG_OLCI_VALUES/original/2022/139/O2022139-rrs442_5-bs-fr.nc'
    # file_corrected = '/mnt/c/DATA_LUIS/OCTAC_WORK/INC_NEG_OLCI_VALUES/slurm/2022/139/O2022139-rrs442_5-bs-fr.nc'
    # from netCDF4 import Dataset
    # dataset_o = Dataset(file_original,'r')
    # dataset_c = Dataset(file_corrected,'r')
    # array_o = dataset_o.variables['RRS442_5'][:]
    # print(type(array_o))
    # array_c = dataset_c.variables['RRS442_5'][:]
    # print(f'ORIGINAL: {np.ma.min(array_o)} {np.ma.max(array_o)}')
    # print(f'CORRECTED: {np.ma.min(array_c)} {np.ma.max(array_c)}' )
    # indices_bad = np.where(array_o<(-400))
    # print(f'NBAD: {len(indices_bad[0])}/{array_o.shape[1]*array_o.shape[2]}')
    # dataset_o.close()
    # dataset_c.close()


def run_hypstar_check():
    dir_data = '/store3/HYPERNETS/INSITU_HYPSTARv2.0bis/JSIT'
    dir_out = '/store3/HYPERNETS/DATA_CHECK/JSIT'
    dir_zip_level1c = '/store3/HYPERNETS/DATA_CHECK/JSIT/level1c'
    year = 2024
    start_date_impl = None
    end_date_impl = None
    n_operational = 0
    work_date = dt(year, 1, 1)
    end_date = dt(year, 12, 31)
    while work_date <= end_date:
        print(f'[INFO] {work_date}')
        yyyy = work_date.strftime('%Y')
        mm = work_date.strftime('%m')
        dd = work_date.strftime('%d')
        dir_data_date = os.path.join(dir_data, yyyy, mm, dd)
        if os.path.isdir(dir_data_date):
            if start_date_impl is None: start_date_impl = work_date
            end_date_impl = work_date
            n_operational = n_operational + 1
            for name in os.listdir(dir_data_date):
                file_in = os.path.join(dir_data_date, name)
                if name.endswith('.nc'):
                    if name.find('L2A_REF') > 0:
                        file_out = os.path.join(dir_out, name)
                    if name.find('L1C_ALL') > 0:
                        file_out = os.path.join(dir_zip_level1c, name)
                    shutil.copy(file_in, file_out)
        work_date = work_date + timedelta(hours=24)

    # dir_out = '/mnt/c/DATA_LUIS/HYPERNETS_WORK/temp'
    # dir_zip_level1c = os.path.join(dir_out,'level1c')
    n_valid = 0
    n_total = 0
    file_zip_level1c = os.path.join(dir_out, 'level1c.zip')
    file_zip_level2a = os.path.join(dir_out, 'level2a.zip')
    file_csv_level2a = os.path.join(dir_out, 'level2a.csv')
    fout = None
    start_date_qc = None
    end_date_qc = None
    for name in os.listdir(dir_out):
        if name.find('L2A_REF') > 0:
            file_h = os.path.join(dir_out, name)
            dataset = Dataset(file_h)
            if fout is None:
                fout = open(file_csv_level2a, 'w')
                wl = dataset.variables['wavelength'][:]
                wl = [str(x) for x in wl.tolist()]
                wl_list = ';'.join(wl)
                fout.write(f'Date;{wl_list}')
            rrs = np.ma.squeeze(dataset.variables['reflectance'][:])
            ts = float(dataset.variables['acquisition_time'][0])
            qf = float(dataset.variables['quality_flag'][0])
            # epsilon = float(dataset.variables['epsilon'][0])
            n_total = n_total + 1
            if qf == 0:
                n_valid = n_valid + 1
                therevalid = dt.utcfromtimestamp(ts)
                if start_date_qc is None:
                    start_date_qc = therevalid
                if end_date_qc is None:
                    end_date_qc = therevalid
                else:
                    if therevalid > end_date_qc:
                        end_date_qc = therevalid

            # if qf == 0 or qf == 268435456:
            #     if (-0.005) <= epsilon <= 0.005:
            #         n_valid = n_valid + 1

            line = dt.utcfromtimestamp(ts).strftime('%Y%m%dT%H%M%S')
            rrs_line = ';'.join([str(x) for x in rrs.tolist()])
            line = f'{line};{rrs_line}'
            fout.write('\n')
            fout.write(line)

            dataset.close()

    fout.close()

    cmd_level1c = f'zip -r {file_zip_level1c} {dir_zip_level1c}'
    cmd_level2a = f'zip {file_zip_level2a} {file_csv_level2a}'

    import subprocess
    if os.path.isfile(file_zip_level1c):
        os.remove(file_zip_level1c)
    if os.path.isfile(file_zip_level2a):
        os.remove(file_zip_level2a)
    prog = subprocess.Popen(cmd_level1c, shell=True, stderr=subprocess.PIPE)
    prog.communicate()
    stats_1c = os.stat(file_zip_level1c)

    prog = subprocess.Popen(cmd_level2a, shell=True, stderr=subprocess.PIPE)
    prog.communicate()
    stats_2a = os.stat(file_zip_level2a)

    print(f'START DATE: {start_date_impl.strftime("%Y-%m-%d")}')
    print(f'END DATE: {end_date_impl.strftime("%Y-%m-%d")}')
    print(f'N OPERATIONAL DAYS: {n_operational}')
    print(f'LEVEL 1C: {stats_1c.st_size} bytes ({stats_1c.st_size / (1024 * 1024)} Mb)')
    print(f'LEVEL 2A: {stats_2a.st_size} bytes ({stats_2a.st_size / (1024 * 1024)} Mb)')
    print(f'#VALID MEASURMENTS: {n_valid} / {n_total}')
    print(f'START DATE QC: {start_date_qc.strftime("%Y-%m-%d")}')
    print(f'END DATE QC: {end_date_qc.strftime("%Y-%m-%d")}')


def run_test(year):
    if year == -2:  ##test mdb files
        old_dir = '/mnt/c/DATA_LUIS/DOORS_WORK/Extracts_2024/AERONET_OC'
        # new_dir = '/mnt/c/DATA_LUIS/DOORS_WORK/Extracts_2024/'
        names = [
            'MDB_CMEMS_OLCI_300M_CMEMS_OBS-OC_BLK_BGC_20190827T000000_20240818T000000_AERONET_Section-7_Platform.nc',
            'MDB_CMEMS_OLCI_300M_CMEMS_OBS-OC_BLK_BGC_20160401T000000_20231220T000000_AERONET_Galata_Platform.nc',
            'MDB_CMEMS_OLCI_300M_CMEMS_OBS-OC_BLK_BGC_20160401T000000_20190808T000000_AERONET_Gloria.nc']

        for name in names:
            print('--------------------------------------------------------------------------------------------------')
            print(name)
            file_old = os.path.join(old_dir, name)
            # file_new = os.path.join(new_dir,name)
            dataset_old = Dataset(file_old, 'r')
            # dataset_new = Dataset(file_new,'r')
            rrs_old = dataset_old.variables['satellite_Rrs'][:]
            # rrs_new = dataset_new.variables['satellite_Rrs'][:]
            print('DATA: ', rrs_old.size, np.ma.count(rrs_old), np.ma.min(rrs_old), np.ma.max(rrs_old))
            # print('NEW: ', rrs_new.size, np.ma.count(rrs_new), np.ma.min(rrs_new), np.ma.max(rrs_new))

            # indices_neg = np.logical_and(rrs_old.mask == False, rrs_old < (-10))
            # rrs_old_bad = rrs_old[indices_neg]
            # rrs_new_bad = rrs_new[indices_neg]
            # print('OLD BAD: ',np.ma.min(rrs_old_bad),'<->',np.ma.max(rrs_old_bad))
            # print('NEW BAD: ', np.ma.min(rrs_new_bad), '<->', np.ma.max(rrs_new_bad))
            #
            # indices_pos = np.logical_and(rrs_old.mask==False,rrs_old>(-10))
            # rrs_old_good = rrs_old[indices_pos]
            # rrs_new_good = rrs_new[indices_pos]
            # check_rrs = rrs_new_good/rrs_old_good
            # print('CHECK RRS: ',np.ma.min(check_rrs),np.ma.max(check_rrs))
            # for name_var in dataset_old.variables:
            #     if name_var=='satellite_Rrs':
            #         continue
            #     var_old = dataset_old.variables[name_var][:]
            #     var_new = dataset_new.variables[name_var][:]
            #     check_var = var_old/var_new
            #     print(f'CHECK  {name_var}', np.ma.min(check_var), np.ma.max(check_var))

            dataset_old.close()
            # dataset_new.close()

        return

    if year == -3:  ##basic test for extract file
        dir_extracts = '/mnt/c/DATA_LUIS/DOORS_WORK/Extracts_2024/extracts_cmems_olci'

        nvalid = 0
        for name in os.listdir(dir_extracts):
            file_extract = os.path.join(dir_extracts, name)
            dataset = Dataset(file_extract, 'r')
            rrs = dataset.variables['satellite_Rrs'][:]
            if np.ma.count(rrs) > 0:
                nvalid = nvalid + 1
                if np.ma.min(rrs) < (-10):
                    print(name, '->', np.ma.min(rrs), '<------------------>', np.ma.max(rrs))
            dataset.close()
        print('nvalid', nvalid)
        return

    # dir_extracts = '/mnt/c/DATA_LUIS/DOORS_WORK/Extracts_2024/extracts_cmems_olci'
    # source_folder = '/mnt/c/DATA_LUIS/DOORS_WORK/SOURCES'
    # file_out = f'/mnt/c/DATA_LUIS/DOORS_WORK/NegData_{year}.csv'

    dir_extracts = '/store3/DOORS/extracts/cmems_olci'
    source_folder = '/dst04-data1/OC/OLCI/daily_v202311_bc'
    file_out = f'/store3/DOORS/extracts/NegData_{year}.csv'

    fcsv = open(file_out, 'w')
    fcsv.write('Date;WL;RRS_O;RRS_Oa;RRS_Ob;CHL_O;CHL_Oa;CHL_Ob')
    bands = None
    bands_str = None
    for name in os.listdir(dir_extracts):

        file_extract = os.path.join(dir_extracts, name)
        time_obj = dt.strptime(name.split('_')[-3], '%Y%m%d')
        if 0 < year != time_obj.year:
            continue
        print(f'[INFO] Extract file: {name}')
        if bands is None:
            dataset = Dataset(file_extract, 'r')
            bands = dataset.variables['satellite_bands'][:]
            bands_str = [f'{b:.2f}'.replace('.', '_').replace('_00', '').replace('_50', '_5') for b in bands]
            dataset.close()
        fcsv = check_file_extract(file_extract, source_folder, time_obj, bands, bands_str, fcsv)

    # file_extract = '/mnt/c/DATA_LUIS/DOORS_WORK/Extracts_2024/extracts_cmems_olci/extract_CMEMS_OLCI_300m_20240619_1404_635.nc'

    fcsv.close()


def check_file_extract(file_extract, source_folder, time_obj, bands, bands_str, fcsv):
    yyyy = time_obj.strftime('%Y')
    jjj = time_obj.strftime('%j')
    source_folder_date = os.path.join(source_folder, yyyy, jjj)
    if not os.path.isdir(source_folder_date):
        return fcsv

    y_point = int(file_extract[:-3].split('_')[-2])
    x_point = int(file_extract[:-3].split('_')[-1])
    rmin = y_point - 12
    rmax = y_point + 13
    cmin = x_point - 12
    cmax = x_point + 13

    ##CHL
    file_chl_b = os.path.join(source_folder_date, f'Ob{yyyy}{jjj}-chl-bs-fr.nc')
    if not os.path.isfile(file_chl_b):
        return fcsv
    dataset_chl_b = Dataset(file_chl_b)
    chl_here_b = dataset_chl_b['CHL'][0, rmin:rmax, cmin:cmax]
    dataset_chl_b.close()
    file_chl = os.path.join(source_folder_date, f'O{yyyy}{jjj}-chl-bs-fr.nc')
    dataset_chl = Dataset(file_chl)
    chl_here = dataset_chl['CHL'][0, rmin:rmax, cmin:cmax]
    dataset_chl.close()
    file_chl_a = os.path.join(source_folder_date, f'Oa{yyyy}{jjj}-chl-bs-fr.nc')
    dataset_chl_a = Dataset(file_chl_a)
    chl_here_a = dataset_chl_a['CHL'][0, rmin:rmax, cmin:cmax]
    dataset_chl_a.close()

    dataset = Dataset(file_extract, 'r')
    rrs = dataset.variables['satellite_Rrs'][:]
    dataset.close()

    rrs_o = np.ma.masked_all(rrs.shape)

    for idx, b in enumerate(bands_str):
        rrs_here = np.ma.squeeze(rrs[0, idx, :, :])
        file_a = os.path.join(source_folder_date, f'Oa{yyyy}{jjj}-rrs{b}-bs-fr.nc')
        file_b = os.path.join(source_folder_date, f'Ob{yyyy}{jjj}-rrs{b}-bs-fr.nc')
        if os.path.exists(file_a) and os.path.exists(file_b):
            dataset_a = Dataset(file_a)
            rrs_a_here = dataset_a.variables[f'RRS{b}'][0, rmin:rmax, cmin:cmax]
            dataset_a.close()
            dataset_b = Dataset(file_b)
            rrs_b_here = dataset_b.variables[f'RRS{b}'][0, rmin:rmax, cmin:cmax]
            dataset_b.close()

            indices_neg = np.logical_and(rrs_here.mask == False, rrs_here < (-10))
            n_neg = np.count_nonzero(indices_neg)
            if n_neg > 0:
                rrs_here_neg = rrs_here[indices_neg]
                rrs_a_here_neg = rrs_a_here[indices_neg]
                rrs_b_here_neg = rrs_b_here[indices_neg]
                chl_here_neg = chl_here[indices_neg]
                chl_here_a_neg = chl_here_a[indices_neg]
                chl_here_b_neg = chl_here_b[indices_neg]

                for ihere in range(n_neg):
                    line = f'{time_obj.strftime("%Y-%m-%d")};{bands[idx]};{rrs_here_neg[ihere]};{rrs_a_here_neg[ihere]};{rrs_b_here_neg[ihere]};{chl_here_neg[ihere]};{chl_here_a_neg[ihere]};{chl_here_b_neg[ihere]}'
                    fcsv.write('\n')
                    fcsv.write(line)

        file_o = os.path.join(source_folder_date, f'O{yyyy}{jjj}-rrs{b}-bs-fr.nc')
        if os.path.exists(file_o):
            dataset_o = Dataset(file_o)
            var_dat = dataset_o.variables[f'RRS{b}'][0, rmin:rmax, cmin:cmax]
            dataset_o.close()
            rrs_o[0, idx, :, :] = var_dat[:, :]

    check = rrs / rrs_o
    print(
        f'[INFO] Check for: {time_obj.strftime("%Y-%m-%d")} Min. ratio: {np.ma.min(check)} Max. ratio: {np.ma.min(check)}')

    return fcsv


def correct_negative_values_from_mdb(input_path, output_path):
    source_folder = args.source_path
    if not os.path.isdir(source_folder):
        print(f'[ERROR] Source folder: {source_folder} is not a valid directory')
        return
    dataset = Dataset(input_path, 'r')
    bands = dataset.variables['satellite_bands'][:]
    bands_str = [f'{b:.2f}'.replace('.', '_').replace('_00', '').replace('_50', '_5') for b in bands]
    rrs = dataset.variables['satellite_Rrs'][:]
    rrs_final = rrs.copy()
    time = dataset.variables['satellite_time'][:]
    dataset.close()

    nmu = rrs.shape[0]
    dims = None
    for imu in range(nmu):
        time_obj = dt.utcfromtimestamp(float(time[imu]))
        print(f'[INFO] Index MU: {imu} Date: {time_obj.strftime("%Y-%m-%d")}')
        if dims is None:
            lat_array_source, lon_array_source = get_latlon_arrays_from_source(source_folder, time_obj)
            if lat_array_source is not None and lon_array_source is not None:
                dims = get_dims(input_path, lat_array_source, lon_array_source)
                print(f'[INFO] Dims have been defined as: {dims}')
        if dims is not None:
            rrs_new = get_rrs_new(source_folder, rrs, imu, time_obj, bands_str, dims)
            rrs_final[imu, :, :, :] = rrs_new[:, :, :]

    create_new_file_with_corrected_rrs(input_path, output_path, rrs_final)


def correct_negative_values_from_extracts(input_path, output_path):
    source_folder = args.source_path
    if not os.path.isdir(source_folder):
        print(f'[ERROR] Source folder: {source_folder} is not a valid directory')
    bands = None
    bands_str = None
    lat_array_source = None
    lon_array_source = None
    for name in os.listdir(input_path):
        if not name.endswith('.nc'): continue
        print(f'[INFO] Working with file: {name}')
        input_file = os.path.join(input_path, name)
        output_file = os.path.join(output_path, name)
        dataset = Dataset(input_file, 'r')
        rrs = dataset.variables['satellite_Rrs'][:]
        time = float(dataset.variables['satellite_time'][0])
        time_obj = dt.utcfromtimestamp(time)
        if bands is None:
            bands = dataset.variables['satellite_bands'][:]
            bands_str = [f'{b:.2f}'.replace('.', '_').replace('_00', '').replace('_50', '_5') for b in bands]
        if lat_array_source is None and lon_array_source is None:
            lat_array_source, lon_array_source = get_latlon_arrays_from_source(source_folder, time_obj)
        dataset.close()
        dims = get_dims(input_file, lat_array_source, lon_array_source)
        rrs_new = get_rrs_new(source_folder, rrs, 0, time_obj, bands_str, dims)
        rrs[0, :, :, :] = rrs_new[:, :, :]
        create_new_file_with_corrected_rrs(input_file, output_file, rrs)


def create_new_file_with_corrected_rrs(input_file, output_file, rrs):
    from netCDF4 import Dataset
    input_dataset = Dataset(input_file)
    ncout = Dataset(output_file, 'w', format='NETCDF4')

    # copy global attributes all at once via dictionary
    ncout.setncatts(input_dataset.__dict__)

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    for name, variable in input_dataset.variables.items():
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                             complevel=6)
        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(input_dataset[name].__dict__)
        if name == 'satellite_Rrs':
            ncout[name][:] = rrs[:]
        else:
            ncout[name][:] = input_dataset[name][:]

    ncout.close()
    input_dataset.close()


def get_dims(file_extract, lat_array_source, lon_array_source):
    if lat_array_source is not None and lon_array_source is not None:
        dataset = Dataset(file_extract)
        lat_c = float(dataset.variables['satellite_latitude'][0, 12, 12])
        lon_c = float(dataset.variables['satellite_longitude'][0, 12, 12])
        dataset.close()
        y_point = np.argmin(np.abs(lat_c - lat_array_source))
        x_point = np.argmin(np.abs(lon_c - lon_array_source))
    else:
        y_point = int(file_extract[:-3].split('_')[-2])
        x_point = int(file_extract[:-3].split('_')[-1])
    rmin = y_point - 12
    rmax = y_point + 13
    cmin = x_point - 12
    cmax = x_point + 13
    return [y_point, x_point, rmin, rmax, cmin, cmax]


def get_latlon_arrays_from_source(source_folder, time_obj):
    yyyy = time_obj.strftime('%Y')
    jjj = time_obj.strftime('%j')
    source_folder_date = os.path.join(source_folder, yyyy, jjj)
    if not os.path.isdir(source_folder_date):
        return None, None
    file_ref = os.path.join(source_folder_date, f'O{yyyy}{jjj}-rrs400-bs-fr.nc')
    if not os.path.exists(file_ref):
        return None, None
    dataset_ref = Dataset(file_ref)
    lat_array = dataset_ref.variables['lat'][:]
    lon_array = dataset_ref.variables['lon'][:]
    dataset_ref.close()
    return lat_array, lon_array


def get_rrs_new(source_folder, rrs, index_rrs, time_obj, bands_str, dims):
    rmin = dims[2]
    rmax = dims[3]
    cmin = dims[4]
    cmax = dims[5]
    rrs_new = np.ma.squeeze(rrs[index_rrs, :, :, :])

    yyyy = time_obj.strftime('%Y')
    jjj = time_obj.strftime('%j')
    source_folder_date = os.path.join(source_folder, yyyy, jjj)
    if not os.path.isdir(source_folder_date):
        return rrs_new

    print(
        f'[INFO][BEFORE] --> RRS for index {index_rrs} Shape (nbandsx25x25) {rrs_new.shape}. Min value: {np.ma.min(rrs_new)}')

    for idx, b in enumerate(bands_str):
        rrs_here = np.ma.squeeze(rrs[index_rrs, idx, :, :])
        file_a = os.path.join(source_folder_date, f'Oa{yyyy}{jjj}-rrs{b}-bs-fr.nc')
        file_b = os.path.join(source_folder_date, f'Ob{yyyy}{jjj}-rrs{b}-bs-fr.nc')
        if os.path.exists(file_a) and os.path.exists(file_b):
            dataset_a = Dataset(file_a)
            rrs_a_here = dataset_a.variables[f'RRS{b}'][0, rmin:rmax, cmin:cmax]
            dataset_a.close()
            dataset_b = Dataset(file_b)
            rrs_b_here = dataset_b.variables[f'RRS{b}'][0, rmin:rmax, cmin:cmax]
            dataset_b.close()

            indices_neg = np.logical_and(rrs_here.mask == False, rrs_here < (-10))
            n_neg = np.count_nonzero(indices_neg)
            if n_neg > 0:
                indices_neg_a = np.logical_and(rrs_a_here.mask == False, rrs_here < (-10))
                indices_neg_b = np.logical_and(rrs_b_here.mask == False, rrs_here < (-10))
                # print('n_neg',n_neg,' a: ',np.count_nonzero(indices_neg_a), 'b:', np.count_nonzero(indices_neg_b))
                if np.count_nonzero(indices_neg_a) > 0:
                    rrs_here[indices_neg_a] = rrs_a_here[indices_neg_a]
                if np.count_nonzero(indices_neg_b) > 0:
                    rrs_here[indices_neg_b] = rrs_b_here[indices_neg_b]
                rrs_new[idx, :, :] = rrs_here[:, :]

    print(
        f'[INFO][AFTER] --> RRS for index {index_rrs} Shape (nbandsx25x25) {rrs_new.shape}. Min value: {np.ma.min(rrs_new)}')
    return rrs_new


def add_instrument_id(input_path, output_path):
    rename_file = False
    if output_path is None:
        rename_file = True
        output_path = os.path.join(os.path.dirname(input_path), 'Temp.nc')
    ##WRITTING OUPPUT DATASET
    dataset = Dataset(input_path)
    ncout = Dataset(output_path, 'w', format='NETCDF4')
    # copy global attributes all at once via dictionary
    ncout.setncatts(dataset.__dict__)
    # copy dimensions
    for name, dimension in dataset.dimensions.items():
        name_new = 'insitu_bands' if name == 'insitu_original_bands' else name
        ncout.createDimension(name_new, (len(dimension) if not dimension.isunlimited() else None))
    ncout.createDimension('instrument_id', 1)

    # copy variables
    for name, variable in dataset.variables.items():
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue
        dims_l = list(variable.dimensions)
        if 'insitu_original_bands' in dims_l:
            dims_l[dims_l.index('insitu_original_bands')] = 'insitu_bands'
        if name == 'insitu_original_bands':
            dims_l = ['instrument_id', 'insitu_bands']
        print(f'[INFO] Variable: {name} Dimensions: {dims_l}')
        ncout.createVariable(name, variable.datatype, tuple(dims_l), fill_value=fill_value, zlib=True, complevel=6)
        ncout[name].setncatts(dataset[name].__dict__)
        ncout[name][:] = dataset[name][:]

    # new variable
    var_instrument_id = ncout.createVariable('institu_instrument_id', 'i4', ('satellite_id', 'insitu_id'),
                                             fill_value=-999, zlib=True, complevel=6)
    var_instrument_id.description = 'Instrument id'
    var_instrument_id.flag_values = [0]
    var_instrument_id.flag_meanings = ['N/A']
    instrument_id = dataset.variables['insitu_time'][:]
    instrument_id[instrument_id.mask == False] = 0
    var_instrument_id[:] = instrument_id[:]

    ncout.close()
    dataset.close()
    if rename_file:
        os.rename(output_path, input_path)


def check_required_params(param_list):
    b = True
    for param in param_list:
        if not args.__dict__[param]:
            print(f'[ERROR] {param} is required for mode {args.mode}')
            b = False
    return b


# %%
if __name__ == '__main__':
    main()
