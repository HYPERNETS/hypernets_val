import sys

import __init__,argparse,os,warnings
import numpy as np
import pandas as pd
from COMMON import args_functions
from OPTIONS.OptionsManager import OptionsManager
from MDBFile import MDBFile
import MDBWritter as mdbw
from datetime import datetime as dt
from datetime import timezone

warnings.simplefilter('ignore', UserWarning)
code_home = os.path.dirname(os.path.dirname(__init__.__file__))

parser = argparse.ArgumentParser(description="Algorithms implementations from MDB files.")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
parser.add_argument('-c', "--config_file", help="Config File.",required=True)
parser.add_argument('-i', "--input_path", help="Input MDB path",required=True)
parser.add_argument('-o', "--output_path", help="Path to output")
parser.add_argument('-s', "--section",help="Section to be processed",required=True)
args = parser.parse_args()

class FlagOptions:

    def __init__(self,config_file,section,verbose):
        self.verbose = verbose
        self.omanager = None
        self.gmanager = None

        general_options_file = os.path.join(code_home, 'OPTIONS', 'flag_options.ini')
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
        if not 'type_flag' in options_global:
            return None
        if options_global['type_flag'] is None:
            return None
        soptions, required = self.gmanager.get_retrieve_options(options_global['type_flag'])
        if soptions is None:
            return None

        options_specific = self.omanager.get_options_as_dict(self.section, soptions, required)
        if options_global is None or options_specific is None:
            return None
        return options_global | options_specific

    def get_flag_types(self):
        soptions, required = self.gmanager.get_retrieve_options('GLOBAL')
        return soptions['type_flag']['list_values']


def get_output_path(ext):
    if not args.output_path:
        ref = (-1)*len(ext)
        output_path = os.path.join(os.path.dirname(args.input_path),f'{os.path.basename(args.input_path)[:ref]}_{args.section}{ext}')
    else:
        output_path = args.output_path
    try:
        os.makedirs(os.path.dirname(output_path),exist_ok=True)
    except Exception as ex:
        print(f'[ERROR] Output directory {os.path.dirname(output_path)} does not exist and could not be created. Exception {ex}')
        return None
    return output_path

def run_polygon_shp_flag(options):
    try:
        import geopandas
    except:
        print(f'[ERROR] geopandas package is required to run run_polygon_shp_flag. ')
        return
    try:
        from shapely import Point
    except:
        print(f'[ERROR] shapely package is required to run run_polygon_shp_flag')
        return

    file_shp = options['input_shape']
    if file_shp is None:
        print(f'[ERROR] input_shape option with the input file shape in required in the configuration file, section {args.section}')
        return
    if not os.path.exists(file_shp):
        print(f'[ERROR] {file_shp} is not available')
        return
    try:
        df_geo = geopandas.read_file(file_shp)
    except Exception as ex:
        print(f'[ERROR] Error reading {file_shp} using geopandas')
        return
    crs_lat_lon = "EPSG:4326"
    crs_original = str(df_geo.crs)
    if crs_original!=crs_lat_lon:
        try:
            df_geo = df_geo.to_crs(crs_lat_lon)
        except:
            print(f'[ERROR] Original crs {crs_original} of the shape file is not lat/lon WGS84 and could not be converted')
            return

    ngeoms = len(df_geo.index)
    col_flag = options['col_flag']
    if col_flag is None:
        all_flags = [f'Flag_{x+1}' for x in range(ngeoms)]
    else:
        if col_flag not in df_geo.columns:
            print(f'[ERROR] {col_flag} is not a valid column in the dbf file associated with the shape file')
            return
        all_flags = df_geo[col_flag][:].tolist()

    flag_list_options = options['flag_list']
    flag_list_implemented = {}
    if flag_list_options is None:
        iflag_value = 1
        for idx,flag_meaning in enumerate(all_flags):
            flag_meaning = str(flag_meaning)##make sure flag meaning is a string
            if flag_meaning not in flag_list_implemented:
                flag_list_implemented[flag_meaning]={
                    'index_geom': [idx],
                    'flag_value': iflag_value
                }
                iflag_value  = iflag_value + 1
            else:
                flag_list_implemented[flag_meaning]['index_geom'].append(idx)

    # else: ##to be implemented
    #     for flag_meaning in flag_list_options:
    #         idx = all_flags.index(idx)
    for flag_meaning in flag_list_implemented:
        print(f'[INFO] Flag meaning: {flag_meaning} Flag value: {flag_list_implemented[flag_meaning]["flag_value"]} Geometries: {flag_list_implemented[flag_meaning]["index_geom"]}')

    n_other_col = 0
    if options['other_col'] is not None:
        n_other_col = len(options['other_col'])
        for flag_meaning in flag_list_implemented:
            indexes_geom = flag_list_implemented[flag_meaning]['index_geom']
            other_col_obj = ['NAv']*n_other_col
            for icol in range(n_other_col):
                col_name = options['other_col'][icol]
                list_v = df_geo[col_name][indexes_geom].tolist()
                other_col_obj[icol] = ','.join(list_v)
            flag_list_implemented[flag_meaning]['other_col_obj'] = other_col_obj

    is_csv = True if args.input_path.endswith('.csv') else False

    if is_csv:
        try:
            df_csv = pd.read_csv(args.input_path,sep=';')
        except:
            print(f'[ERROR] {args.input_path} is not a valid CSV file')
    else:
        mfile = MDBFile(args.input_path)
        if not mfile.VALID_NC:
            print(f'[ERROR] {args.input_path} is not a valid NetCDF file')
            return
    ext = '.csv' if is_csv else '.nc'
    output_path = get_output_path(ext)
    if output_path is None:
        return

    if args.verbose:
        print(f'[INFO] Output file: {output_path}')

    if not os.path.exists(output_path) and not is_csv:
        mdbw.copy_nc(args.input_path,output_path)

    lat_variable = options['input_var_lat']
    lon_variable = options['input_var_lon']

    if is_csv:
        if not lat_variable in df_csv.columns:
            print(f'[ERROR] {lat_variable} is not available in the file {args.input_path}')
            return
        if not lon_variable in df_csv.columns:
            print(f'[ERROR] {lon_variable} is not available in the file {args.input_path}')
            return
        lat_array = np.array(df_csv[lat_variable][:])
        lon_array = np.array(df_csv[lon_variable][:])
        lon_lat_array = np.ma.stack([lon_array, lat_array], axis=1)
    else:##MDB
        if not lat_variable in mfile.variables:
            print(f'[ERROR] {lat_variable} is not available in the file {args.input_path}')
            return
        if not lon_variable in mfile.variables:
            print(f'[ERROR] {lon_variable} is not available in the file {args.input_path}')
            return
        lat_array = mfile.variables[lat_variable][:]
        lon_array = mfile.variables[lon_variable][:]
        if lat_array.shape!=lon_array.shape:
            print(f'[ERROR] {lat_variable} and {lon_variable} does not have the same shapes')
            return
        shape_orig = lat_array.shape
        dims_orig = mfile.variables[lat_variable].dimensions
        lon_lat_array = np.ma.stack([lon_array.flatten(),lat_array.flatten()],axis=1)

    ndata = lon_lat_array.shape[0]
    indices_valid = np.where(np.ma.count(lon_lat_array,axis=1) == 2)
    lon_lat_array = lon_lat_array[indices_valid[0],:]

    point_list = [Point(lon_lat_array[idx,0],lon_lat_array[idx,1]) for idx in range(lon_lat_array.shape[0])]

    flag_array_valid = np.zeros(lon_lat_array.shape[0])
    flag_meanings =  list(flag_list_implemented.keys())
    nobj = len(flag_meanings)
    flag_values = [flag_list_implemented[meaning]['flag_value'] for meaning in flag_meanings]
    ndata_by_flag = np.zeros(nobj+1)
    other_col = None
    if n_other_col>0:
        other_col = np.zeros((nobj+1,n_other_col)).astype(np.object_)

    for iobj in range(nobj):
        meaning = flag_meanings[iobj]
        flag_value = flag_values[iobj]
        indexes_geom = flag_list_implemented[meaning]['index_geom']
        if len(indexes_geom)==0:
            geom = df_geo.loc[indexes_geom[0], 'geometry']
            res = geom.contains(point_list)
            flag_array_valid[res==True] = flag_value
        else:
            for igeom in indexes_geom:
                geom = df_geo.loc[igeom, 'geometry']
                res = geom.contains(point_list)
                flag_array_valid[res == True] = flag_value


        ndata_by_flag[iobj] = np.sum(flag_array_valid==flag_value)
        if other_col is not None:
            other_col[iobj,:] = np.array(flag_list_implemented[meaning]['other_col_obj'])

    ndata_by_flag[-1] = lon_lat_array.shape[0]-np.sum(ndata_by_flag[:-1])

    if options['other_col'] is not None:
        col_info = options['other_col'] + ['NPoints']
    else:
        col_info = ['NPoints']
    df_info = pd.DataFrame(index=flag_meanings+['No Classified'],columns=col_info)


    for icol,col_name in enumerate(options['other_col']):
        data_col = other_col[:,icol]
        data_col[-1] = 'NaN'
        df_info[col_name] = data_col
    df_info['NPoints'] = ndata_by_flag
    outfile_csv = output_path.replace(ext, '_summary.csv')

    name_var = args.section if options['name_var'] is None else options['name_var']
    name_var = f'flag_{name_var}' if not name_var.startswith('flag') else name_var

    df_info.to_csv(outfile_csv, sep=';', index_label=name_var)

    if args.verbose:
        print(f'[INFO] Creating variable: {name_var}')

    if is_csv:
        flag_array = np.ma.masked_all((ndata,))
        flag_array[indices_valid[0]] = flag_array_valid[:]
        flag_array = np.ma.filled(flag_array,-999.0)
        flag_array_meanings = ['NaN']*ndata
        flag_array_other_col = None
        if n_other_col>0:
            flag_array_other_col = np.zeros((ndata,n_other_col)).astype(np.object_)
            flag_array_other_col[:] = 'NaN'
        for idata in range(ndata):
            if flag_array[idata]!=-999.0:
                v = int(flag_array[idata])-1
                if v>=0:
                    meaning_here = flag_meanings[v]
                    flag_array_meanings[idata] = meaning_here
                    if flag_array_other_col is not None:
                        flag_array_other_col[idata,:] = np.array(flag_list_implemented[meaning_here]['other_col_obj'])
                else:
                    flag_array_meanings[idata] = 'NotClassified'
                    if flag_array_other_col is not None:
                        flag_array_other_col[idata,:] = np.array(['NotClassified']*n_other_col)

        df_csv[name_var] = flag_array_meanings[:]
        if flag_array_other_col is not None:
            for icol in range(n_other_col):
                col_name = options['other_col'][icol]
                df_csv[col_name] = flag_array_other_col[:,icol]
        df_csv.to_csv(output_path,sep=';',index=None)
        print(f'[INFO] Completed csv file with flag: {output_path}')
    else:##MDB
        mw = mdbw.MDBWritter(args.input_path,output_path)


        flag_array = np.ma.masked_all((ndata,))
        flag_array[indices_valid[0]] = flag_array_valid[:]
        flag_array = np.reshape(flag_array,shape_orig)

        mw.add_variable(name_var,flag_array,dims_orig,'i4',-999)
        attrs = {
            'flag_values': flag_values,
            'flag_meanings': " ".join(flag_meanings)
        }
        for col_name in  options['other_col']:
            col_name_list = df_geo[col_name][:]
            str_val = None
            for flag_meaning in flag_list_implemented:
                indexes_geom = flag_list_implemented[flag_meaning]['index_geom']
                str_here = "+".join(col_name_list[indexes_geom].tolist())

                if str_val is None:
                    str_val = str_here
                else:
                    str_val = f'{str_val},{str_here}'
            attrs[col_name] = str_val
        mw.add_attrs_to_variable(name_var,attrs)
        mw.close()

def run_csv_label_flag(options):
    mfile = MDBFile(args.input_path)
    if not mfile.VALID_NC:
        print(f'[ERROR] {args.input_path} is not a valid NetCDF file')
        return
    nc_time_variable =options['time_variable']
    if not nc_time_variable in mfile.variables:
        print(f'[ERROR] time_variable {nc_time_variable} is not available in the NetCDF (MDB) file')
        return
    try:
        df_csv = pd.read_csv(options['input_csv'],sep=';')
    except:
        print(f'[ERROR] input_csv {options["input_csv"]} is not a valid CSV file separated by ;')
        return
    time_csv = get_time_from_dataframe(df_csv,options)
    if time_csv is None:
        return
    if not options['col_flag'] in df_csv:
        print(f'[ERROR] col_flag {options["col_flag"]} is not a valid column flag')
        return
    flag_array = df_csv[options['col_flag']][:]

    if options['flag_meanings'] is not None:
        flag_array_check = np.unique(flag_array).tolist()
        flag_meanings = options['flag_meanings']
        any_flag_present = False
        for fm in flag_meanings:
            any_flag_present = fm in flag_array_check
            if any_flag_present:
                break
        if not any_flag_present:
            print(f'[ERROR] None of the flags given in flag_meanings {flag_meanings} is present in the column {options["col_flag"]} of the CSV file.')
            return
    else:
        flag_meanings = np.unique(flag_array).tolist()

    if options['flag_values'] is not None:
        flag_values = options['flag_values']
        if len(flag_values)!=len(flag_meanings):
            print(f'[ERROR] Number of flag_meanings ({len(flag_meanings)}) {flag_meanings} is not equal to the number of flag_values ({len(flag_values)}) {flag_values}')
            return
    else:
        if options['use_pow2']:
            flag_values = np.power(2,np.arange(len(flag_meanings))).tolist()
        else:
            flag_values = np.arange(1,len(flag_meanings)+1,1).tolist()

    all_flags_csv = {}
    for t, fvalue in zip(time_csv, flag_array):
        if fvalue in flag_meanings:
            all_flags_csv[t] = flag_values[flag_meanings.index(fvalue)]
        else:
            all_flags_csv[t] = 0

    time_nc = mfile.variables[nc_time_variable][:]
    dims_orig  = mfile.variables[nc_time_variable].dimensions
    shape_orig = time_nc.shape
    time_nc = time_nc.flatten()
    ndata = time_nc.shape[0]
    indices_valid = np.where(time_nc.mask==False)
    time_nc = time_nc[indices_valid]
    flag_array_valid = np.zeros(time_nc.shape)
    for it,t in enumerate(time_nc):
        time_here = dt.fromtimestamp(t).astimezone(timezone.utc).strftime('%Y%m%dT%H%M%S')
        if time_here in all_flags_csv:
            flag_array_valid[it] = all_flags_csv[time_here]
    flag_array = np.zeros(ndata)
    flag_array[indices_valid] = flag_array_valid[:]
    flag_array = np.reshape(flag_array,shape_orig)

    name_var = args.section if options['name_var'] is None else options['name_var']
    name_var = f'flag_{name_var}' if not name_var.startswith('flag') else name_var
    output_path = get_output_path('.nc')
    if output_path is None:
        return
    print(f'[INFO] Adding variable {name_var} to output path {output_path}')

    mw = mdbw.MDBWritter(args.input_path, output_path)
    mw.add_variable(name_var, flag_array, dims_orig, 'i4', -999)
    attrs = {
        'flag_values': flag_values,
        'flag_meanings': " ".join(flag_meanings)
    }
    mw.add_attrs_to_variable(name_var, attrs)
    mw.close()
    print(f'[INFO] Completed')

def run_csv_int_flag(options):
    print(f'[INFO] Starting CSV int flag')
    mfile = MDBFile(args.input_path)
    if not mfile.VALID_NC:
        print(f'[ERROR] {args.input_path} is not a valid NetCDF file')
        return
    nc_time_variable =options['time_variable']
    if not nc_time_variable in mfile.variables:
        print(f'[ERROR] time_variable {nc_time_variable} is not available in the NetCDF (MDB) file')
        return
    try:
        df_csv = pd.read_csv(options['input_csv'],sep=';')
    except:
        print(f'[ERROR] input_csv {options["input_csv"]} is not a valid CSV file separated by ;')
        return
    time_csv = get_time_from_dataframe(df_csv,options)
    if time_csv is None:
        return
    if not options['col_flag'] in df_csv:
        print(f'[ERROR] col_flag {options["col_flag"]} is not a valid column flag')
        return
    flag_array = df_csv[options['col_flag']][:]

    if options['flag_values'] is not None:
        flag_values_check = np.unique(flag_array).tolist()
        flag_values = options['flag_values']
        any_flag_present = False
        for fv in flag_values:
            any_flag_present = fv in flag_values_check
            if any_flag_present:
                break
        if not any_flag_present:
            print(f'[ERROR] None of the flag values given in flag_values {flag_values} is present in the column {options["col_flag"]} of the CSV file, with values {flag_values_check}')
            return


    else:
        flag_values = np.unique(flag_array).tolist()

    if options['flag_meanings'] is not None:
        flag_meanings = options['flag_meanings']
        if len(flag_values)!=len(flag_meanings):
            print(f'[ERROR] Number of flag_meanings ({len(flag_meanings)}) {flag_meanings} is not equal to the number of flag_values ({len(flag_values)}) {flag_values}')
            return
    else:
        flag_meanings = [f'flag_{x}' for x in flag_values]


    all_flags_csv = {}
    for t, fvalue in zip(time_csv, flag_array):
        if fvalue in flag_values:
            all_flags_csv[t] = fvalue
        else:
            all_flags_csv[t] = 0

    time_nc = mfile.variables[nc_time_variable][:]
    dims_orig  = mfile.variables[nc_time_variable].dimensions
    shape_orig = time_nc.shape
    time_nc = time_nc.flatten()
    ndata = time_nc.shape[0]
    indices_valid = np.where(time_nc.mask==False)
    time_nc = time_nc[indices_valid]
    flag_array_valid = np.zeros(time_nc.shape)
    for it,t in enumerate(time_nc):
        time_here = dt.fromtimestamp(t).astimezone(timezone.utc).strftime('%Y%m%dT%H%M%S')
        if time_here in all_flags_csv:
            flag_array_valid[it] = all_flags_csv[time_here]
    flag_array = np.zeros(ndata)
    flag_array[indices_valid] = flag_array_valid[:]
    flag_array = np.reshape(flag_array,shape_orig)

    name_var = args.section if options['name_var'] is None else options['name_var']
    name_var = f'flag_{name_var}' if not name_var.startswith('flag') else name_var
    output_path = get_output_path('.nc')
    if output_path is None:
        return
    print(f'[INFO] Adding variable {name_var} to output path {output_path}')

    mw = mdbw.MDBWritter(args.input_path, output_path)
    mw.add_variable(name_var, flag_array, dims_orig, 'i4', -999)
    attrs = {
        'flag_values': flag_values,
        'flag_meanings': " ".join(flag_meanings)
    }
    mw.add_attrs_to_variable(name_var, attrs)
    mw.close()
    print(f'[INFO] Completed')

def get_time_from_dataframe(df,options):
    if options['col_datetime'] is not None:
        col_dt = options['col_datetime']
        format_dt = options['format_col_datetime']
        if col_dt not in df.columns:
            print(f'[ERROR] col_datetime {col_dt} is not a column in the CSV file')
            return None
        if format_dt is None:
            print(f'[ERROR] format_col_datetime is required to parse datetime objects in column {col_dt}. First value: {df.loc[0,col_dt]}')
            return None

        time_array = df[col_dt]
        try:
            time_list = [dt.strptime(x,format_dt).strftime('%Y%m%dT%H%M%S') for x in time_array]
        except Exception as ex:
            print(f'[ERROR] Error parsing datetime objects from CSV file: Exception {ex}')
            return None
        return time_list
    elif options['col_date'] is not None:
        col_d = options['col_date']
        format_d = options['format_col_date']
        col_t = options['col_time']
        format_t = options['format_col_time']
        if col_d not in df.columns:
            print(f'[ERROR] col_date {col_d} is not a column in the CSV file')
            return None
        if format_d is None:
            print(f'[ERROR] format_col_date is required to parse date objects in column {col_d}. First value: {df.loc[0, col_d]}')
            return None
        if col_t is None:##use only format date
            time_array = df[col_d]
            try:
                time_list = [dt.strptime(x, format_d).strftime('%Y%m%dT%H%M%S') for x in time_array]
            except Exception as ex:
                print(f'[ERROR] Error parsing date objects from CSV file: Exception {ex}')
                return None
            return time_list
        else:
            if col_t not in df.columns:
                print(f'[ERROR] col_time {col_t} is not a column in the CSV file')
                return None
            if format_t is None:
                print(f'[ERROR] format_col_time is required to parse time objects in column {col_t}. First value: {df.loc[0, col_t]}')
                return None
            format_dt = f'{format_d}T{format_t}'
            d_array = df[col_d]
            t_array = df[col_t]
            try:
                time_list = [dt.strptime(f'{x}T{y}',format_dt).strftime('%Y%m%dT%H%M%S') for x,y in zip(d_array,t_array)]
            except Exception as ex:
                print(f'[ERROR] Error parsing time objects from CSV file: Exception {ex}')
                return None
            return time_list

def run_flag_multiple_insitu(options):
    mfile = MDBFile(args.input_path)
    if not mfile.VALID_NC:
        print(f'[ERROR] {args.input_path} is not a valid NetCDF file')
        return

    var_ref = options['var_ref']
    if not var_ref in mfile.variables:
        print(f'[ERROR] Reference variable {var_ref} is not available in the NetCDF (MDB) file {args.input_path}')

    array = mfile.variables[var_ref][:]
    dims_orig = mfile.variables[var_ref].dimensions
    nvalid = np.ma.count(array,axis=1)

    flag_values = [int(x) for x in np.unique(nvalid)]
    flag_meanings = [str(x) for x in flag_values]
    name_var = args.section if options['name_var'] is None else options['name_var']
    name_var = f'flag_{name_var}' if not name_var.startswith('flag') else name_var
    output_path = get_output_path('.nc')
    if output_path is None:
        return
    print(f'[INFO] Adding variable {name_var} to output path {output_path}')

    mw = mdbw.MDBWritter(args.input_path, output_path)
    mw.add_variable(name_var, nvalid, (dims_orig[0],), 'i4', -999)
    attrs = {
        'flag_values': flag_values,
        'flag_meanings': " ".join(flag_meanings)
    }
    mw.add_attrs_to_variable(name_var, attrs)
    mw.close()
    print(f'[INFO] Completed')

def get_var_ref_for_dimensions(options,mfile,param_index):
    if options['var_dim'] is not None:
        var_dim = options['var_dim']
        if var_dim not in mfile.variables:
            print(f'[ERROR] Reference variable for dimensions {var_dim} is not available in the NetCDF (MDB) file')
            return None
    else:
        var_dim = None
        ndims_prev = sys.maxsize
        index = 0
        while f'range_{index}' in options:
            var_here = options[f'{param_index}{index}'][0]
            if var_here in mfile.variables:
                ndims = len(mfile.variables[var_here].dimensions)
                if ndims<ndims_prev:
                    var_dim = var_here
                    ndims_prev = ndims
            index = index + 1
        if var_dim is None:
            print(f'[ERROR] Reference variable for dimensions could not be retrieved, please use the option var_dim to set it')
            return None
    return var_dim

def get_range_array(mfile,options,general_options):
    var_name = options[0]
    if not var_name in mfile.variables:
        print(f'[ERROR] {var_name} is not available in the netCDF (MDB) file. Range array could not be creted')
        return None
    array_here = mfile.variables[var_name][:]

    minV, maxV = np.nan,np.nan
    if options[1]!='MIN':
        try:
            minV = float(options[1])
        except:
            print(f'[ERROR] Maximum value for range {options[1]} is not a valid number. Range array could not be created')
            return None
    if options[2]!='MAX':
        try:
            maxV = float(options[2])
        except:
            print(f'[ERROR] Maximum value for range {options[2]} is not a valid number. Range array could not be created')
            return None

    dim_c_info = None
    if len(options) == 4:  ##change of dimension
        dim_c = options[3]
        if dim_c not in general_options:
            print(f'[ERROR] Dimension change  {dim_c} is not defined in the configuration file')
            return None
        else:
            dim_c_info = general_options[dim_c]

    if options[1]=='MIN' and options[2]!='MAX':
        array_range = array_here < maxV
    elif options[1]!='MIN' and options[2]=='MAX':
        array_range = array_here >= minV
    elif options[1]!='MIN' and options[2]!='MAX':
        if minV==maxV:
            array_range = array_here==minV
        else:
            array_range = (array_here>=minV) & (array_here<maxV)
    else:
        print(f'[ERROR] Error defining the ranges, is should be: var,MIN,maxV or var,minV,MAX or var,minV,maxV, with minV and maxV being numerical values')
        print(f'[ERROR] But options are: {options}')
        return None

    if dim_c_info is not None and dim_c_info[0]=='all': #the condition must be in all the valid points
        nvalid = np.ma.count(array_range,axis=1)
        svalid = np.ma.sum(array_range,axis=1)
        array_range = nvalid == svalid
    if dim_c_info is not None and dim_c_info[0]=='any': #the condition must be in any of the valid points
        svalid = np.ma.sum(array_range,axis=1)
        array_range = svalid > 0

    return array_range

def get_combine_array(options,array_dict):
    condition = options[0]
    array_list = options[1:]
    if len(array_list)==1:
        return array_dict[array_list[0]]
    if condition=='AND':
        final_array = np.logical_and(array_dict[array_list[0]],array_dict[array_list[1]])
        if len(array_list)>=3:
            for iarray in range(2,len(array_list)):
                final_array = np.logical_and(final_array,array_dict[array_list[iarray]])

        return final_array




def run_flag_ranges(options):
    mfile = MDBFile(args.input_path)
    if not mfile.VALID_NC:
        print(f'[ERROR] {args.input_path} is not a valid NetCDF file')
        return
    var_dim  = get_var_ref_for_dimensions(options,mfile,'range_')
    if var_dim is None:
        return


    print(f'[INFO] Getting range arrays...')
    ranges_arrays = {}
    index = 0
    while f'range_{index}' in options:
        ranges_arrays[f'range_{index}'] = get_range_array(mfile,options[f'range_{index}'],options)
        array_temp = ranges_arrays[f'range_{index}']
        index = index + 1

    print(f'[INFO] Getting flag array...')
    dims_orig = mfile.variables[var_dim].dimensions
    shape_orig = mfile.variables[var_dim].shape
    array_flag = np.zeros(shape_orig)
    flag_values = []
    flag_meanings = []
    index = 0
    while f'flag_{index}' in options:
        fvalues = options[f'flag_{index}']
        try:
            flag_meanings.append(fvalues[0].strip())
            flag_value = int(fvalues[1].strip())
            flag_values.append(flag_value)
            array_combine = get_combine_array(fvalues[2:],ranges_arrays)
            array_flag[array_combine] = flag_value
        except Exception as ex:
            print(f'[ERROR] Error creating flag_{index}. Exception {ex}')

        index = index + 1

    if len(flag_values)>1:
        indices_sort = np.argsort(np.array(flag_values))
        flag_values = flag_values[indices_sort[0]]
        flag_meanings = flag_meanings[indices_sort[0]]
    print(f'[INFO] Flag values: {flag_values}')
    print(f'[INFO] Flag meanings: {flag_meanings}')

    name_var = args.section if options['name_var'] is None else options['name_var']
    name_var = f'flag_{name_var}' if not name_var.startswith('flag') else name_var
    output_path = get_output_path('.nc')
    if output_path is None:
        return

    print(f'[INFO] Adding variable {name_var} to output path {output_path}')
    mw = mdbw.MDBWritter(args.input_path, output_path)
    mw.add_variable(name_var, array_flag, dims_orig, 'i4', -999)
    attrs = {
        'flag_values': flag_values,
        'flag_meanings': " ".join(flag_meanings)
    }
    mw.add_attrs_to_variable(name_var, attrs)
    mw.close()
    print(f'[INFO] Completed')

def create_new_flag(mfile,options):
    var_flag = options['var_flag']
    output_path = get_output_path('.nc')
    if output_path is None:
        return
    flag_values_new = options['flag_values']
    flag_meanings_new = options['flag_meanings']

    flag_dims = options['flag_dims']
    flag_dims_var = tuple(flag_dims)
    if flag_dims==['satellite_id']:
        n_mu_total = len(mfile.dimensions['satellite_id'])
        array_flag_new = np.ma.masked_all(n_mu_total)
        array_flag_new[:] = options['default_value']
    elif flag_dims==['satellite_id','insitu_id']:
        n_mu_total = len(mfile.dimensions['satellite_id'])
        n_insitu_day = len(mfile.dimensions['insitu_id'])
        array_flag_new = np.ma.masked_all((n_mu_total,n_insitu_day))
        array_flag_new[:,0] = options['default_value']
    else:
        print(f'[ERROR] {flag_dims} is not a valid dimension combination')
        return

    print(f'[INFO] Creating variable {var_flag} in output path {output_path}')
    mw = mdbw.MDBWritter(args.input_path, output_path)
    mw.add_variable(var_flag,array_flag_new,flag_dims_var,'i4',-999)
    attrs = {
        'flag_values': flag_values_new,
        'flag_meanings': " ".join(flag_meanings_new)
    }
    index = 0
    while f'attr_{index}' in options:
        edit_at = options[f'attr_{index}']
        at = edit_at[0]
        at_val = ",".join(edit_at[1:])
        attrs[at] = at_val
        index = index + 1
    mw.add_attrs_to_variable(var_flag, attrs)
    mw.close()
    print(f'[INFO] Completed')

def run_edit_flag(options):
    mfile = MDBFile(args.input_path)
    if not mfile.VALID_NC:
        print(f'[ERROR] {args.input_path} is not a valid NetCDF file')
        return
    var_flag = options['var_flag']
    if options['flag_values'] is not None and options['flag_meanings'] is not None and var_flag not in mfile.variables:
        ##spetial case it creates a new flag
        create_new_flag(mfile,options)
        return

    if var_flag not in mfile.variables:
        print(f'[ERROR] Flag variable {var_flag} is not in the NetCDF (MDB) file {args.input_path}')
        return

    if options['flag_values'] is not None and options['flag_meanings'] is not None:
        ##rewritte of flag_values and flag_meanings
        array_flag = mfile.variables[var_flag][:]
        array_flag_new = array_flag.copy()
        flag_values_prev = mfile.variables[var_flag].flag_values
        flag_meanings_prev = mfile.variables[var_flag].flag_meanings.split(' ')
        flags_prevs = {m:v for m,v in zip(flag_meanings_prev,flag_values_prev)}
        flag_values_new = options['flag_values']
        flag_meanings_new = options['flag_meanings']

        for meaning,value in zip(flag_meanings_new,flag_values_new):
            if meaning in flag_meanings_prev:
                value_prev = flags_prevs[meaning]
                print(f'[INFO] Reassigning flag meaning {meaning} from {value_prev} to {value}')
                array_flag_new[array_flag==value_prev] = value
    else:
        array_flag_new = mfile.variables[var_flag][:]
        flag_values_new = mfile.variables[var_flag].flag_values
        flag_meanings_new = mfile.variables[var_flag].flag_meanings.split(' ')


    output_path = get_output_path('.nc')
    if output_path is None:
        return
    print(f'[INFO] Updating variable {var_flag} in output path {output_path}')
    mw = mdbw.MDBWritter(args.input_path, output_path)
    mw.update_data_variable(var_flag,array_flag_new)
    attrs = {
        'flag_values': flag_values_new,
        'flag_meanings': " ".join(flag_meanings_new)
    }
    index = 0
    while f'attr_{index}' in options:
        edit_at = options[f'attr_{index}']
        at = edit_at[0]
        at_val = ",".join(edit_at[1:])
        attrs[at] = at_val
        index = index + 1
    mw.add_attrs_to_variable(var_flag, attrs)
    mw.close()
    print(f'[INFO] Completed')


def test():
    dir_base = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/OCTAC_WORK/BAL_EVOLUTION_202411/MATCH-UPS_ANALYSIS_2024/MDBsV3'
    #file_bal = os.path.join(dir_base,'MDB_CMEMS_MULTI_BAL_STATIONS_19970909_20241006.nc')
    file_kat = os.path.join(dir_base,'MDB_CMEMS_MULTI_KATTEGAT_STATIONS_19970904_20231206.nc')
    from netCDF4 import Dataset
    # for file in [file_bal,file_kat]:
    #     dset = Dataset(file)
    #     longitude = dset.variables['insitu_longitude'][:]
    #     print(np.ma.min(longitude),np.ma.max(longitude))
    #     dset.close()
    dset = Dataset(file_kat)
    insitu_time = dset.variables['insitu_time'][:]
    for idx in range(100):
        itime = insitu_time[idx,0]
        print(dt.fromtimestamp(itime).astimezone(timezone.utc).strftime('%Y-%m-%d %H:%M:%S'))
    dset.close()
    return True
def main():
    # if test():
    #     return

    required_args = {
        'input_path':{'type':'input_file'},
        'config_file':{'type':'input_file'},
        'section':{'type':'str'}
    }
    args_d = args_functions.get_args_as_dict(args,required_args,False)
    if args_d is None:
        return
    print(f'[INFO] Working from configuration file: {args.config_file}')
    print(f'[INFO] Section: {args.section}')
    flag_options = FlagOptions(args.config_file, args.section, args.verbose)
    if not flag_options.is_valid:
        print(f'[ERROR] Problem retrieving flags options for {args.section}. Please review flag_options.ini')
        return
    options = flag_options.get_options()
    if options is None:
        return


    if options['type_flag']=='polygons_shp':
        run_polygon_shp_flag(options)

    if options['type_flag']=='csv_label':
        run_csv_label_flag(options)

    if options['type_flag']=='csv_int':
        run_csv_int_flag(options)

    if options['type_flag']=='multiple_insitu':
        run_flag_multiple_insitu(options)

    if options['type_flag']=='flag_ranges':
        run_flag_ranges(options)

    if options['type_flag']=='edit_flag':
        run_edit_flag(options)



if __name__ == '__main__':
    main()