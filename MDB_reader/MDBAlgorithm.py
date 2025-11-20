import argparse
import warnings
import os
import numpy as np
from datetime import datetime as dt

import pandas as pd
import pytz
from netCDF4 import Dataset

from MDBFile import MDBFile
from CommonMu import COMMON_MU
from MDB_builder.INSITU_base import INSITUBASE
import __init__,sys
from COMMON import args_functions
import MDBWritter
from OPTIONS.OptionsManager import OptionsManager
code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)

warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', RuntimeWarning)

parser = argparse.ArgumentParser(description="Algorithms implementations from MDB files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument("-m", "--mode", help="Mode", choices=["CONFIGFILE", "CYANOFLAG","BALMLP202411"])
parser.add_argument('-c', "--config_file", help="Config File.")
parser.add_argument('-i', "--input_path", help="Input MDB path", required=True)
parser.add_argument('-o', "--output", help="Path to output")
parser.add_argument('-s', "--section",help="Section to be processed for CONFIGFILE mode")
args = parser.parse_args()


class AlgorithOptions:

    def __init__(self,config_file,section,verbose):
        self.verbose = verbose
        self.omanager = None
        self.gmanager = None

        general_options_file = os.path.join(code_home, 'OPTIONS', 'algorithm_options.ini')
        self.gmanager = OptionsManager(general_options_file, None)
        self.omanager = OptionsManager(config_file, None)

        self.is_valid = self.gmanager.is_valid() and self.omanager.is_valid()
        self.section = None
        if self.is_valid:
            if section is None:
                print(f'[ERROR] --section (-s) argument is required')
                self.is_valid = False
            else:
                if not self.omanager.options.has_section(section):
                    print(f'[ERROR] Section {section} is not available in the configuration file')
                    self.is_valid = False
                else:
                    self.section = section




    def get_general_options(self,section):
        if not self.is_valid:
            return None
        soptions, required = self.gmanager.get_retrieve_options(section)

        if soptions is not None:
            options = self.omanager.get_options_as_dict(section, soptions,required)
            return options
        else:
            return None

    def get_options(self):
        if not self.is_valid:
            return None
        soptions, required = self.gmanager.get_retrieve_options('GLOBAL')
        if soptions is None:
            return None
        options_global = self.omanager.get_options_as_dict(self.section, soptions, required)
        if not 'type_algo' in options_global:
            return None
        if options_global['type_algo'] is None:
            return None
        soptions, required = self.gmanager.get_retrieve_options(options_global['type_algo'])
        options_specific = self.omanager.get_options_as_dict(self.section, soptions, required)
        options_specific['required_args'] = self.gmanager.get_required_args(options_global['type_algo'])
        if options_global is None or options_specific is None:
            return None
        return options_global | options_specific

    def get_algo_types(self):
        soptions, required = self.gmanager.get_retrieve_options('GLOBAL')
        return soptions['type_algo']['list_values']

class MDBProcessing:
    def __init__(self,path_mdb):
        self.path_mdb = path_mdb
        self.mfile = MDBFile(path_mdb)

    def close_mfile(self):
        self.mfile.close()

    def run(self,reference,options,output):
        type_algo = options['type_algo']
        if type_algo=='expression':
            self.run_expression(reference,options,output)
        elif type_algo=='common_mu':
            self.run_common_mu(reference,options,output)
        elif type_algo=='subset':
            self.run_subset(options,output)
        else:
            print(f'[ERROR] {type_algo} processing is not implemented yet. Please add the corresponding function in the run() method in the MDBProcessing.py class (file MDBAlgorithm.py)')

    def run_subset(self,options,output):
        type_subset = options['type_subset']


        array_subset = None

        if type_subset=='basic_filter':
            if options['basic_filter'] is None:
                print(f'[ERROR] basic_filter option with a valid expression is required for type_subset: basic_filter')
                return array_subset
            expression, dims, shape = self.get_expression_from_mfile(options['basic_filter'])
            array_subset, dims = self.eval_expression(expression,{},dims,shape)

        if array_subset is None:
            return
        if os.path.exists(output):
            os.remove(output)
        writer = MDBWritter.MDBWritter(self.mfile, output)
        writer.create_subset(array_subset)


    def run_common_mu(self,reference,options,output):
        output_variable = reference if options['output_variable'] is None else options['output_variable']
        output_variable_mu = f'{output_variable}_mu'
        exist_output_variable = output_variable in self.mfile.variables
        exist_output_variable_mu  = output_variable_mu in self.mfile.variables
        self.close_mfile()##using path to mdb
        options['type_cmu']='single_mdb'
        options['input_file'] = self.path_mdb
        if options['output_mode']=='subset_file' and not args_functions.check_arg_type_impl(output,'output_file_nc'):
            print(f'[ERROR] Argument output should be a valid output NC file. Common mu processing stopped.')
            return

        #options['output'] = output


        cmu = COMMON_MU(options)
        if not exist_output_variable:
            cmu_array = cmu.run_single_mdb()
            self.write_output_variable(output,output_variable,cmu_array,None,'i2',-999)
        else:
            print(f'[WARNING] Output variable {output_variable} already exists. To repeat the process, please choose another output variable name')

        if options['create_mu_variable']:
            if not exist_output_variable_mu:
                cmu_array_mu = cmu.get_mu_variable(output_variable)
                self.write_output_variable(output, output_variable_mu, cmu_array_mu, None, 'i2', -999)
            else:
                print(
                    f'[WARNING] Output mu variable {output_variable} already exists. To repeat the process, please choose another output variable name')


    def run_expression(self,reference,options,output):
        n_exp = 0
        while f'expression_{n_exp}' in options:
            if options[f'expression_{n_exp}'] is not None:

                n_exp = n_exp + 1

        if n_exp==0:
            print(f'[ERROR] At least one valid expression is required for the "expression" processing.')
            return
        output_variable = reference if options['output_variable'] is None else options['output_variable']
        arrays = {}
        dims = None
        for i_exp in range(n_exp):
            expression, dims, shape = self.get_expression_from_mfile(options[f'expression_{i_exp}'])
            local_dict = {}
            if i_exp>=1:
                expression,local_dict = self.get_local_dict(expression,arrays)
            arrays[i_exp],dims = self.eval_expression(expression, local_dict, dims, shape)

        self.close_mfile()
        array = arrays[n_exp-1]##last array
        if array is None:
            return
        self.write_output_variable(output,output_variable,array,dims,None,-999.0)

    def get_local_dict(self,expression,arrays):
        loc_dict = locals()
        for idx in range(0,len(arrays)):
            r = expression.find(f'$result_{idx}$')
            if r>=0:
                expression = expression.replace(f'$result_{idx}$',f'result_{idx}')
                loc_dict[f'result_{idx}'] = arrays[idx]

        return expression,loc_dict

    def get_expression_from_mfile(self,expression):
        vars = expression.split('$')
        var_names = []
        for var in vars:
            if var in self.mfile.variables:
                expression = expression.replace(f'${var}$', f"self.mfile.variables['{var}']")
                var_names.append(var)
        dims = None
        shape = None
        if len(var_names) > 0:
            for var_name in var_names:
                dims_h = self.mfile.variables[var_name].dimensions
                shape_h = self.mfile.variables[var_name].shape
                if dims is None:
                    dims = dims_h
                    shape = shape_h
                else:
                    if dims_h != dims or shape_h != shape:
                        print(f'[WARNING] Inconsistency in the dimension and shape of the input variables:')
                        print(f'[WARNING] ->Variable {var_name}: {dims_h}/{shape_h} but expected {dims}/{shape}')


        return expression, dims, shape

    def eval_expression(self,expression,local_dict,dims,shape):
        print(f'[INFO] Evaluating expression:')
        print(f'[INFO] {expression}')
        try:
            if len(local_dict)>0:
                array = eval(expression,globals(),local_dict)
            else:
                array = eval(expression)
            if isinstance(array,np.ndarray) and array.shape!=shape:
                print(f'[WARNING] Inconsistency in the output array. Obtained shape: {array.shape} but expected {shape}')
                return array,None
            if isinstance(array, np.ndarray):
                print(f'[INFO] Output array shape: {array.shape} Dimensions: {dims}')

            return array,dims
        except Exception as ex:
            print(f'[ERROR] Exception evaluating expression: {ex}')
            return None,None


    def write_output_variable(self,output,output_variable,array,dims,dtype,fill_value):
        if output is None:
            writer = MDBWritter.MDBWritter(None, self.path_mdb)
        else:
            MDBWritter.copy_nc(self.mfile.file_path, output)
            writer = MDBWritter.MDBWritter(None, output)

        print(f'[INFO] Output variable: {output_variable}')
        writer.add_variable(output_variable, array, dims, None, None)
        writer.close()



def do_test():
    ####TEMPORAL
    # file_nc = '/mnt/c/DATA_LUIS/OCTAC_WORK/MATCH-UPS_ANALYSIS_2024/BAL/MDBs/MDBr__MULTI_CCI_1KM_OC-CCI-BALCHL2021_19970316T000000_20240415T000000.nc'
    file_nc = '/mnt/c/DATA_LUIS/OCTAC_WORK/MATCH-UPS_ANALYSIS_2024/BAL/MDBs/MDBr__S3_OLCI_300M_CMEMS-OLCI_20160404T000000_20240415T000000.nc'
    dataset = Dataset(file_nc)
    insitu_id = dataset.variables['mu_insitu_CHLA_id'][:]
    insitu_time = dataset.variables['insitu_time'][:]
    flag_cyano_array = dataset.variables['satellite_flag_cyano'][:]
    file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/MATCH-UPS_ANALYSIS_2024/BAL/MDBs/FlagCyanoOlci.csv'
    fout = open(file_out, 'w')
    fout.write('DATE;HOUR;FLAG_CYANO;N_VALID;N_CYANO')
    for idx in range(flag_cyano_array.shape[0]):
        if idx % 100 == 0: print(idx)
        insitu_id_here = insitu_id[idx]
        if np.ma.is_masked(insitu_id_here):
            continue

        insitu_time_here = float(insitu_time[idx][insitu_id_here])
        insitu_time_here_dt = dt.utcfromtimestamp(insitu_time_here)
        date = insitu_time_here_dt.strftime('%Y-%m-%d')
        time = insitu_time_here_dt.strftime('%H:%M:%S')
        # rrs_555 = satellite_Rrs_555[idx, 12, 12]
        # rrs_670 = satellite_Rrs_670[idx, 12, 12]
        flag_cyano_here = flag_cyano_array[idx, 12, 12]
        flag_cyano_values = flag_cyano_array[idx, 11:14, 11:14]
        flag_cyano_values_good = flag_cyano_values[~flag_cyano_values.mask]
        n_all = len(flag_cyano_values_good)
        flag_cyano_centre = flag_cyano_values_good[flag_cyano_values_good == flag_cyano_here]
        n_cyano = len(flag_cyano_centre)
        # line = f'{date};{time};{rrs_555};{rrs_670};{flag_cyano_here}'
        line = f'{date};{time};{flag_cyano_here};{n_all};{n_cyano}'
        fout.write('\n')
        fout.write(line)
    fout.close()
    dataset.close()

    ####

    return True


def create_mdb_from_csv(sensor):
    sensor = sensor.upper()
    dir_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/MATCH-UPS_ANALYSIS_2024/MDBs'
    file_out = None
    if sensor=='MULTI':
        file_csv = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/MATCH-UPS_ANALYSIS_2024/CSV_MATCH-UPS/MULTI/Baltic_CHLA_Valid_AllSources_1997-2023_FINAL_TIMEFILTERED_complete.csv'
        file_out = os.path.join(dir_out, 'MDBr__MULTI_CCI_1KM_OC-CCI-BALCHL202411_19970909T000000_20231219T000000.nc')

    if sensor=='OLCI':
        file_csv = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION_202411/MATCH-UPS_ANALYSIS_2024/CSV_MATCH-UPS/OLCI/Baltic_CHLA_Valid_AllSources_2016-2023_FINAL_OLCI_rrs_chl_3x3_filtered_match-ups_complete.csv'
        file_out = os.path.join(dir_out, 'MDBr__S3_OLCI_300M_CMEMS-OLCI_20160501T000000_20240415T000000_BAL202411.nc')

    if file_out is None:
        return

    nc_out = Dataset(file_out, 'w')
    nc_out.createDimension('mu_id', size=None)
    nc_out.createDimension('satellite_id', size=None)
    df = pd.read_csv(file_csv, sep=';')
    col_names =df.columns.tolist()

    if sensor=='MULTI':
        mu_variables = ['LATITUDE', 'LONGITUDE', 'INSITU_CHL','FLAG_CYANO_PORC_MODE','DISTANCE']
        flag_variables_str = ['SOURCE_ORIG', 'SOURCE', 'FlagPrec', 'FlagBrando', 'FlagOldNew', 'FLAG_DISTANCE']

    if sensor=='OLCI':
        mu_variables = ['LATITUDE', 'LONGITUDE', 'INSITU_CHL','FLAG_CYANO_PORC_MODE']
        flag_variables_str = ['SOURCE_ORIG', 'SOURCE', 'FlagPrec', 'FlagBrando', 'FlagOldNew']

    flag_variables_num = ['satellite_CHL_NVALID','FLAG_CDF','BLOOM','SUB_SURFACE','SURFACE']


    for col in col_names:
        if col.startswith('satellite_RRS'):
            mu_variables.append(col)
        if col.startswith('satellite_CHL') and not col.endswith('NVALID'):
            mu_variables.append(col)
        if col.startswith('satellite_CDF'):
            mu_variables.append(col)
        if col.startswith('satellite_WEIGHT'):
            mu_variables.append(col)
        if col.startswith('FLAG_CYANO') and col!='FLAG_CYANO_PORC_MODE':
            flag_variables_num.append(col)
        if col.startswith('use'):
            flag_variables_num.append(col)
        if col.startswith('FLAG_ERROR'):
            flag_variables_num.append(col)
        if col.startswith('FLAG_NOERROR'):
            flag_variables_num.append(col)


        # if col.startswith('chl_') or col.startswith('cdf_'):
        #     mu_variables.append(col)
        # if col.endswith('_cdf') or col.startswith('use_cdf'):
        #     flag_variables_num.append(col)
        # if col.endswith('_NVALID'):
        #     flag_variables_num.append(col)




    datetime_array = np.array([dt.strptime(x,'%Y-%m-%dT%H:%M:%S').replace(tzinfo=pytz.UTC).timestamp() for x in df['DATETIME']])
    var = nc_out.createVariable('time','f4',('satellite_id',),zlib=True,complevel=6,fill_value=-999.0)
    var[:] = datetime_array

    for var_name in mu_variables:
        print('Creating mu variable:',var_name)
        var = nc_out.createVariable(var_name.upper(),'f4',('mu_id',),zlib=True,complevel=6,fill_value=-999.0)
        var[:] = df[var_name]

    for var_name in flag_variables_num:
        print('Creating flag num',var_name)
        var_array = np.array(df[var_name])
        if var_name=='FLAG_CYANO' or var_name=='FLAG_CYANO_CENTRAL' or var_name=='FLAG_CYANO_MODE':
            flag_values = [0,1,2,3]
            flag_meanings = 'NO_BLOOM SUB_SURFACE_BLOOM SURFACE_BLOOM BOTH_BLOOMS'
        elif var_name == 'FLAG_CDF':
            flag_values = [1,2,4,8,6,10,12,14]
            flag_meanings = 'NO_CDF CDF_MLP_3B CDF_MLP_4B CDF_MLP_5B CDF_MLP_3B+4B CDF_MLP_3B+5B CDF_MLP_4B+5B CDF_MLP_3B+4B+5B'
        elif var_name == 'BLOOM' or var_name == 'SUB_SURFACE' or var_name == 'SURFACE':
            flag_values = [1,2]
            flag_meanings = 'NO_BLOOM BLOOM'
        elif var_name == 'FLAG_CYANO_HOMOGENEITY':
            flag_values = [0,1]
            flag_meanings = ['CYANO_MIX','CYANO_HOMEGENEOUS']
        else:
            values_unique = np.unique(var_array).tolist()
            flag_values = values_unique
            flag_meanings = ' '.join([str(x) for x in flag_values])
        var = nc_out.createVariable(var_name.upper(), 'i4', ('satellite_id',), zlib=True, complevel=6, fill_value=-999)
        var[:] = var_array[:]
        var.flag_meanings = flag_meanings
        var.flag_values = flag_values




    for var_name in flag_variables_str:
        print('Flag str:', var_name)
        var_array_obj = np.array(df[var_name])
        obj_list = np.unique(var_array_obj).tolist()

        var_array = np.zeros((var_array_obj.shape))
        var_array[:] = -999
        flag_values = []
        for idx,obj in enumerate(obj_list):
            flag_value = 2**idx
            flag_values.append(flag_value)
            var_array[var_array_obj==obj] = flag_value
        flag_meanings = ' '.join(obj_list)
        var = nc_out.createVariable(var_name.upper(), 'i4', ('satellite_id',), zlib=True, complevel=6, fill_value=-999)
        var[:] = var_array[:]
        var.flag_meanings = flag_meanings
        var.flag_values = flag_values



    nc_out.close()


    return True


def main():
    # if create_mdb_from_csv('OLCI'):
    #     return
    # if do_test():
    #     return
    print('Started MDBAlgorithm')
    if args.mode=='CONFIGFILE':
        run_from_config_file()
        return
    input_path = args.input_path
    if not os.path.exists(input_path):
        print(f'[ERROR] Input path {input_path} does not exist')
        return
    try:
        dset = Dataset(input_path)
        dset.close()
    except:
        print(f'[ERROR] Input path {input_path} is not a valid MDB (NetCDF File)')
        return
    output_path = None
    if args.output:
        output_path = args.output
        if not os.path.isdir(os.path.dirname(output_path)):
            try:
                os.mkdir(os.path.dirname(output_path))
            except:
                print(
                    f'[ERROR] Ouput path {os.path.basename(output_path)} is not valid as {os.path.dirname(output_path)} is not a valid directory')
                return
        if output_path.endswith('.nc'):
            print(f'[ERROR] Output path {output_path} should be a NC file (.nc)')
            return

    if args.mode == 'CYANOFLAG':
        create_cyano_flag(input_path, output_path)

    if args.mode == 'BALMLP202411':
        try:
            from  baltic_202411 import BALTIC_202411_PROCESSOR
        except:
            print(f'[ERROR] baltic_202411 code is not available')
        if output_path is None:
            input_dir = os.path.dirname(args.input_path)
            name = os.path.basename(args.input_path)[:-3]+'_BAL202411.nc'
            output_path = os.path.join(input_dir,name)
            print(f'[INFO] Ouput path set to: {output_path}')
        bprocessor = BALTIC_202411_PROCESSOR(None, False)
        bprocessor.run_from_mdb_file(args.input_path,output_path)
        return

def run_from_config_file():
    algo_options = AlgorithOptions(args.config_file,args.section,args.verbose)
    if not algo_options.is_valid:
        print(f'[ERROR] Problem retrieving algorithm options for {args.section}. Please review algorithm_options.ini')
        return
    options = algo_options.get_options()
    print(options)
    args_d = args_functions.get_args_as_dict(args,options['required_args'],False)
    if args_d is None:
        return
    mprocessing = MDBProcessing(args_d['input_path'])
    mprocessing.run(args.section,options,args.output)

def create_cyano_flag(input_path, output_path):
    if output_path is None:
        dataset_w = Dataset(input_path, 'a')
    else:
        ibase = INSITUBASE(None)
        dataset_w = ibase.copy_nc(input_path, output_path)

    ##satellite variable
    if not 'satellite_Rrs' in dataset_w.variables or not 'satellite_bands' in dataset_w.variables:
        print(f'[ERROR] satellite_Rrs and satellite_bands are required to compute CYANOFLAG')
        dataset_w.close()
        if output_path is not None:
            os.remove(output_path)
        return
    th_555_sub = 4.25e-3
    th_670_sur = 1.22e-3
    satellite_bands = dataset_w.variables['satellite_bands'][:]
    index_555 = np.ma.argmin(np.ma.abs(satellite_bands - 555.0))
    wl_555 = satellite_bands[index_555]
    diff_555 = abs(wl_555 - 555.0)
    if diff_555 > 5:
        print(f'[ERROR] Band at 555 nm is not avaialable (nearest band is {wl_555})')
        return
    else:
        print(f'[INFO] Wavelength for 555 nm (sub-surface blooms): {wl_555}')

    index_670 = np.argmin(np.abs(satellite_bands - 670.0))
    wl_670 = satellite_bands[index_670]
    diff_670 = abs(wl_670 - 670.0)
    if diff_670 > 5:
        print(f'[ERROR] Band at 670 nm is not avaialable (nearest band is {wl_670})')
        return
    else:
        if wl_670 > 670:  ##665 better than 673.75
            satellite_bands_l = satellite_bands[satellite_bands < 670]
            index_670_l = np.argmin(np.abs(satellite_bands_l - 670.0))
            wl_670_l = satellite_bands[index_670_l]
            diff_670_l = abs(wl_670_l - 670.0)
            if diff_670_l <= 5:
                index_670 = index_670_l
                wl_670 = wl_670_l
                diff_670 = diff_670_l
        print(f'[INFO] Wavelength for 670 nm (sub-surface blooms): {wl_670}')

    satellite_Rrs = dataset_w.variables['satellite_Rrs']

    if diff_555 == 0.0:
        satellite_Rrs_555 = np.ma.squeeze(satellite_Rrs[:, index_555, :, :])
    else:
        ##apply band shifting
        print(f'[INFO] Applying band shifting from {wl_555} nm to 555.0 nm')
        from BSC_QAA import bsc_qaa_EUMETSAT as bsc
        satellite_Rrs_555 = np.ma.zeros((satellite_Rrs.shape[0], satellite_Rrs.shape[2], satellite_Rrs.shape[3]))
        for index_mu in range(satellite_Rrs.shape[0]):
            rrs_in = np.ma.squeeze(satellite_Rrs[index_mu, :, :, :])
            ndata_valid = np.ma.count(rrs_in)
            if ndata_valid == 0:
                satellite_Rrs_555[index_mu, :, :] = np.ma.masked
            else:
                satellite_Rrs_555[index_mu, :, :] = bsc.bsc_qaa(rrs_in, satellite_bands, np.ma.array([555.0]))

    if diff_670 == 0.0:
        satellite_Rrs_670 = np.ma.squeeze(satellite_Rrs[:, index_670, :, :])
    else:
        ##apply band shifting
        print(f'[INFO] Applying band shifting from {wl_670} nm to 670.0 nm')
        from BSC_QAA import bsc_qaa_EUMETSAT as bsc
        satellite_Rrs_670 = np.ma.zeros((satellite_Rrs.shape[0], satellite_Rrs.shape[2], satellite_Rrs.shape[3]))
        for index_mu in range(satellite_Rrs.shape[0]):
            rrs_in = np.ma.squeeze(satellite_Rrs[index_mu, :, :, :])
            ndata_valid = np.ma.count(rrs_in)
            if ndata_valid == 0:
                satellite_Rrs_670[index_mu, :, :] = np.ma.masked
            else:
                satellite_Rrs_670[index_mu, :, :] = bsc.bsc_qaa(rrs_in, satellite_bands, np.ma.array([670.0]))

    if 'satellite_flag_cyano' not in dataset_w.variables:
        print(f'[INFO] Creating variable satellite_flag_cyano')
        satellite_cyano_array = np.zeros(satellite_Rrs_555.shape)
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 >= th_555_sub, satellite_Rrs_670 >= th_670_sur)] = 3
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 >= th_555_sub, satellite_Rrs_670 < th_670_sur)] = 1
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 < th_555_sub, satellite_Rrs_670 >= th_670_sur)] = 2
        satellite_cyano_array[np.logical_and(satellite_Rrs_555.mask, satellite_Rrs_670.mask)] = -999
        satellite_cyano_var = dataset_w.createVariable('satellite_flag_cyano', 'i2',
                                                       ('satellite_id', 'rows', 'columns'),
                                                       zlib=True, complevel=6, fill_value=-999.0)
        satellite_cyano_var.descripton = 'Satellite Cyano Flag'
        satellite_cyano_var.flag_masks = [0, 1, 2, 3]
        satellite_cyano_var.flag_meanings = "NO_BLOOM SUB-SURFACE_BLOOM SURFACE_BLOOM BOTH_BLOOMS"
        satellite_cyano_var[:] = satellite_cyano_array[:]
    else:
        print(f'[WARNING] Variable satellite_flag_cyano already exits. Skipping...')
        satellite_cyano_array = np.zeros(satellite_Rrs_555.shape)
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 >= th_555_sub, satellite_Rrs_670 >= th_670_sur)] = 3
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 >= th_555_sub, satellite_Rrs_670 < th_670_sur)] = 1
        satellite_cyano_array[np.logical_and(satellite_Rrs_555 < th_555_sub, satellite_Rrs_670 >= th_670_sur)] = 2
        satellite_cyano_array[np.logical_and(satellite_Rrs_555.mask, satellite_Rrs_670.mask)] = -999
        # satellite_cyano_var = dataset_w.variables['satellite_flag_cyano']
        # satellite_cyano_var.flag_masks = [0, 1, 2, 3]
        # satellite_cyano_var.flag_meanings = "NO_BLOOM SUB-SURFACE_BLOOM SURFACE_BLOOM BOTH_BLOOMS"
        # satellite_cyano_var[:] = satellite_cyano_array[:]

    if 'flag_cyano' not in dataset_w.variables:
        print(f'[INFO] Creating variable flag_cyano')
        flag_cyano_array = np.ma.squeeze(satellite_cyano_array[:, 12, 12])
        flag_cyano_var = dataset_w.createVariable('flag_cyano', 'i2', ('satellite_id',), zlib=True, complevel=6,
                                                  fill_value=-999)
        flag_cyano_var.flag_values = [0, 1, 2, 3]
        flag_cyano_var.flag_meanings = "NO_BLOOM SUB-SURFACE_BLOOM SURFACE_BLOOM BOTH_BLOOMS"
        flag_cyano_var[:] = flag_cyano_array
    else:
        print(f'[WARNING] Variable flag_cyano already exits. Skipping...')
        # flag_cyano_array = np.ma.squeeze(satellite_cyano_array[:, 12, 12])

        # flag_cyano_var = dataset_w.variables['flag_cyano']
        # flag_cyano_var.flag_values = [0, 1, 2, 3]
        # flag_cyano_var.flag_meanings = "NO_BLOOM SUB-SURFACE_BLOOM SURFACE_BLOOM BOTH_BLOOMS"
        # flag_cyano_var[:] = flag_cyano_array

    # insitu_time = dataset_w.variables['insitu_time'][:]
    # insitu_id = dataset_w.variables['mu_insitu_CHLA_id'][:]

    dataset_w.close()

    if output_path is None:
        output_path = input_path
    print(f'[INFO] Completed. Output file {output_path}')


if __name__ == '__main__':
    main()
