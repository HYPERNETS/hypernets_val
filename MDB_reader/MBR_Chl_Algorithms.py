import os,math,__init__
import pandas as pd
import numpy as np
from netCDF4 import Dataset

class MBR_Chl_Algorithms:

    def __init__(self,algorithm):
        self.algo = algorithm
        self.df_algo = self.get_info_algorithm()
        self.valid_algo = True if self.df_algo is not None else False
        self.max_diff = 2
        self.CI_METHOD = 'NASA'
        self.MERGE_METHOD = 'NASA'
        if self.valid_algo:
            print(f'[INFO] Loaded MBR algorithm: {self.algo}')


    def get_wl(self):
        if not self.valid_algo:
            return None
        blue_wl = np.ma.masked_equal(self.df_algo[['Blue_1','Blue_2','Blue_3','Blue_4']].to_numpy()[0],-999).compressed()
        green_wl = np.ma.masked_equal(self.df_algo[['Green_1','Green_2']].to_numpy()[0],-999).compressed()
        return blue_wl,green_wl

    def get_indices_wl(self,wl_required,wl_array):
        nrequired = len(wl_required)
        narray = len(wl_array)
        wl_required_2d = np.repeat(wl_required, narray).reshape((nrequired, narray))
        wl_array_2d = np.tile(wl_array,nrequired).reshape((nrequired,narray))
        diff = np.abs(wl_required_2d-wl_array_2d)
        min_diff = np.min(diff, axis=1)
        min_diff_indices = np.argmin(diff, axis=1)
        if np.max(min_diff)>self.max_diff:
            for iwl in range(nrequired):
                if min_diff[iwl]>self.max_diff:
                    print(f'[ERROR] Data for wl {wl_required[iwl]} is not available as difference with the nearest band ({wl_array[min_diff_indices[iwl]]}) is greater than the allowed maximum of {self.max_diff} nm')
            return None
        else:
            return min_diff_indices

    def get_coeffs(self):
        if not self.valid_algo:
            return None
        return self.df_algo[['a0', 'a1', 'a2', 'a3', 'a4']].to_numpy()[0]

    def get_rss555_for_ci_nasa(self,rrs,wl_ref):
        min_range = [553,543,548,558,563]
        max_range = [557,547,552,562,567]
        irange = np.where((wl_ref>=min_range) & (wl_ref<=max_range))[0]
        if len(irange)==0:
            print(f'[ERROR] Rrs for green band at 555 nm could not be extracted')
            return None
        iref = str(int(irange[0]))
        rrs = np.ma.masked_less(rrs,0)
        if iref=='0':
            return rrs
        transform = {
            '1':{
                'th':0.001723,
                'a1':0.986,
                'b1':0.081495,
                'a2':1.031,
                'b2':0.000216
            },
            '2': {
                'th': 0.001597,
                'a1': 0.988,
                'b1': 0.062195,
                'a2': 1.014,
                'b2': 0.000128
            },
            '3': {
                'th': 0.001148,
                'a1': 1.023,
                'b1': -0.103624,
                'a2': 0.979,
                'b2': -0.000121
            },
            '4': {
                'th': 0.000891,
                'a1': 1.039,
                'b1': -0.183044,
                'a2': 0.971,
                'b2': -0.000170
            }
        }


        rrs[rrs <  transform[iref]['th']] = np.power(10,transform[iref]['a1']*np.ma.log10(rrs[rrs<transform[iref]['th']])-transform[iref]['b1'])
        rrs[rrs >= transform[iref]['th']] = transform[iref]['a2'] * rrs[rrs >= transform[iref]['th']] - transform[iref]['b2']

        return rrs

    def get_input_array_impl(self,input_array,nbands,iwl):

        if len(input_array.shape)==2 and input_array.shape[1]==nbands:
            iwl_orig = 1
            iwl_dest = 1
        else:
            shape_orig = input_array.shape
            iwl_dest = len(shape_orig)-1
            if iwl==-1:
                iwl = np.where(np.array(shape_orig)==nbands)[0]
                if len(iwl)==0:
                    print(f'[ERROR] Number of bands ({nbands}) was not found in the input array shape: {input_array.shape}')
                    return None
                elif len(iwl)>1:
                    print(f'[ERROR] Number of bands ({nbands}) is available in more than one dimension in the input array shape: {shape_orig}')
                    return None
                else:
                    iwl_orig = iwl[0]
            elif 0 <= iwl <= iwl_dest:
                if shape_orig[iwl]==nbands:
                    iwl_orig = iwl
                else:
                    print(f'[ERROR] iwl {iwl} does not correspond with the number of bands ({nbands}) in the input array shape: {shape_orig}')
                    return None
            else:
                print(f'[ERROR] iwl {iwl} defining the band dimension is not valid, it should be between 0 and {iwl_dest}')
                return None

        if iwl_dest>iwl_orig:
            input_array_impl = np.moveaxis(input_array,iwl_orig,iwl_dest)
        else:
            input_array_impl = input_array

        shape_impl = input_array_impl.shape

        print(f'[INFO] - Index band orig in array: {iwl_orig}. Index band dest in array: {iwl_dest}. Input array shape after moving axis: {shape_impl}')

        if len(shape_impl)>2:

            ndata = np.prod(np.array(shape_impl)[:-1])
            output_shape = (ndata,nbands)

            input_array_impl = np.reshape(input_array_impl,output_shape)

        print(f'[INFO] - Implementation array shape: {input_array_impl.shape}')

        return input_array_impl,shape_impl

    def run_merge_mbr_ci(self,chl_mbr,chl_ci):
        print(f'[INFO] Merging MBR and CI chl-a with method: {self.MERGE_METHOD}')
        if self.MERGE_METHOD=='NASA':
            t1 = 0.25
            t2 = 0.35
            output_array = chl_mbr.copy()
            output_array[chl_ci<=t1] = chl_ci[chl_ci<=t1]
            indices_merge = np.where((chl_ci>t1) & (chl_ci<t2))
            output_array[indices_merge] = (chl_ci[indices_merge]*((t2-chl_ci[indices_merge])/(t2-t1))) + (chl_mbr[indices_merge]*((chl_ci[indices_merge]-t1)/(t2-t1)))
            return output_array
        else:
            print(f'[ERROR] Merge method {self.MERGE_METHOD} is not available. Choose among: NASA')
            return None

    #rrs_blue and rrs_red are  arrays with the same shape. wl_blue, wl_green and wl_red could be 1D n arrays or single values
    def run_color_index_impl(self,rrs_blue,rrs_green,rrs_red,wl_blue,wl_green,wl_red):
        if not isinstance(rrs_blue,np.ndarray) or not isinstance(rrs_red,np.ndarray) or not isinstance(rrs_green,np.ndarray):
            print(f'[ERROR] rrs_blue and rrs_red should be numpy arrays')
            return [None]*2
        if rrs_blue.shape!=rrs_red.shape or rrs_blue.shape!=rrs_green.shape:
            print(f'[ERROR] rrs_blue, rrs_green, and rrs_red should have the same shape but rrs_blue: {rrs_blue.shape} and rrs_red: {rrs_red.shape}')
            return [None] * 2

        shape_orig = rrs_blue.shape
        wl_term = (wl_green-wl_blue)/(wl_red-wl_blue)

        if isinstance(wl_term,(np.float32,np.float64,float)):
            wl_term = np.reshape(np.repeat(wl_term,rrs_blue.size),shape_orig)
        elif isinstance(wl_term,np.ndarray):
            if wl_term.shape != shape_orig:
                print(f'[ERROR] Shape of blue, red and green wavelength arrays ({wl_term.shape}) should be the same as Rrs arrays: {shape_orig}')
                return [None]*2
        else:
            print(f'[ERROR] wl_blue, wl_green and wl_green should be numerical values or arrays with the same shape as rrs_blue and rrs_red')
            return [None] * 2


        valid_indices = np.where((rrs_blue.mask==False) & (rrs_green.mask==False) & (rrs_red.mask==False))
        rrs_blue = rrs_blue[valid_indices]
        rrs_green = rrs_green[valid_indices]
        rrs_red = rrs_red[valid_indices]
        wl_term = wl_term[valid_indices]

        CI = rrs_green - (rrs_blue + wl_term * (rrs_red-rrs_blue))

        if self.CI_METHOD == 'NASA':
            chl_CI = np.power(10,(-0.4287 + 230.47   * CI))
        elif self.CI_METHOD == 'ESA':
            chl_CI = np.power(10,(-0.5379 + 180.9642 * CI))
        else:
            print(f'[ERROR] Color index method {self.CI_METHOD} is not available. Please choose between NASA and ESA')
            return [None]*2

        CI_output = np.ma.masked_all(shape_orig)
        chl_CI_output = np.ma.masked_all(shape_orig)

        CI_output[valid_indices] = CI
        chl_CI_output[valid_indices] = chl_CI

        return CI_output,chl_CI_output



    def run_color_index(self,input_array,wl_array,iwl=-1,output_as_implemented = False):
        print(f'[INFO] Running Color Index algorithm {self.algo}')

        index_blue = int(np.argmin(np.abs(wl_array-433)))
        index_green = int(np.argmin(np.abs(wl_array-555)))
        index_red = int(np.argmin(np.abs(wl_array-670)))
        wl_blue = wl_array[index_blue]
        wl_green = wl_array[index_green]
        wl_red = wl_array[index_red]
        print(f'[INFO] - Blue band: {wl_blue}. Index: {index_blue}.')
        print(f'[INFO] - Green band: {wl_green}. Index: {index_green}.')
        print(f'[INFO] - Red band: {wl_red}. Index: {index_red}.')

        nbands = len(wl_array)
        input_array_impl,shape_impl = self.get_input_array_impl(input_array,nbands,iwl)
        rrs_blue = input_array_impl[:,index_blue]
        rrs_green = input_array_impl[:,index_green]
        rrs_red = input_array_impl[:,index_red]

        rrs_green = self.get_rss555_for_ci_nasa(rrs_green,wl_green)

        CI, chl_CI = self.run_color_index_impl(rrs_blue,rrs_green,rrs_red,wl_blue,wl_green,wl_red)
        if chl_CI is None:
            return
        print(f'[INFO] Created CI and chl_CI arrays with shape {CI.shape}')

        if len(shape_impl)>2 and not output_as_implemented:
            final_output_shape = tuple(np.array(shape_impl)[:-1])
            CI = np.reshape(CI,final_output_shape)
            chl_CI = np.reshape(chl_CI,final_output_shape)

        if not output_as_implemented:
            print(f'[INFO] - Final output  shape: {CI.shape}')

        return CI, chl_CI

    #input_array: ndata x nbands
    #wl_array: nbands
    def run_mbr_algorithm_impl(self,input_array,indices_blue=None,indices_green=None,wl_array=None):
        if not self.valid_algo:
            return None

        if indices_blue is None or indices_green is None:
            if wl_array is None:
                print(f'[ERROR] Indices for blue and green bands could not be obtained without the wavelengths for the bands')
                return None
            blue_wl, green_wl = self.get_wl()
            if indices_blue is None:
                indices_blue = self.get_indices_wl(blue_wl, wl_array)
                if indices_blue is None:
                    return None
            if indices_green is None:
                indices_green = self.get_indices_wl(green_wl, wl_array)
                if indices_green is None:
                    return None

        coefs = self.get_coeffs()
        if coefs is None:
            return None

        ##getting blue and green Rrs
        blue_data = np.ma.array(input_array[:,indices_blue])
        green_data = np.ma.array(input_array[:,indices_green])
        ##negative values are considered non-valid to get the maximum blue and mean green
        blue_data = np.ma.masked_less(blue_data,0)
        green_data = np.ma.masked_less(green_data, 0)
        ##max blue (if needed)
        if blue_data.shape[1]>1:
            blue_data = np.ma.max(blue_data,axis=1)
        else:
            blue_data = np.squeeze(blue_data)
        ##mean green (if needed)
        if green_data.shape[1]>1:
            green_data = np.ma.mean(green_data,axis=1)
        else:
            green_data = np.squeeze(green_data)

        ndata = blue_data.shape[0]


        ##processing only using valid data
        indices_valid = np.where((blue_data.mask == False) & (green_data.mask == False))
        blue_data = blue_data[indices_valid]
        green_data = green_data[indices_valid]
        X = np.ma.log10(blue_data/green_data)
        nvalid = X.shape[0]
        X = np.repeat(X,5).reshape((nvalid,5))
        indices = np.tile(np.arange(5),nvalid).reshape((nvalid,5))
        coefs = np.tile(coefs,nvalid).reshape((nvalid,5))
        res = np.power(X,indices)*coefs
        log_chl = np.sum(res,axis=1)
        chl = np.power(10,log_chl)

        #creating output array
        output_chl = np.ma.masked_all((ndata,))
        output_chl[indices_valid] = chl[:]

        print(f'[INFO] - Processed chl-a values: {chl.shape[0]}/{ndata}')


        return output_chl

    def run_mbr_algorithm(self,input_array,wl_array=None,iwl=-1):
        if not self.valid_algo:
            return None


        blue_wl, green_wl = self.get_wl()

        if wl_array is None:#in this case, input_array should contain the blue and green bands in the same order
            wl_array = np.concat([blue_wl,green_wl])
            indices_blue = np.arange(0,len(blue_wl)).astype(np.int8)
            indices_green = np.arange(len(blue_wl),len(wl_array)).astype(np.int8)
        else:
            indices_blue = self.get_indices_wl(blue_wl, wl_array)
            indices_green = self.get_indices_wl(green_wl, wl_array)


        nbands = len(wl_array)

        print(f'[INFO] Running MRBr Algorithm: {self.algo}')
        print(f'[INFO] - Blue bands: {blue_wl}')
        print(f'[INFO] - Green bands: {green_wl}')
        print(f'[INFO] - Input array shape: {input_array.shape}')
        print(f'[INFO] - Number of bands: {nbands}. Blue indices: {indices_blue}. Green indices: {indices_green}')



        input_array_impl, shape_impl = self.get_input_array_impl(input_array,nbands,iwl)
        output_array = self.run_mbr_algorithm_impl(input_array_impl,indices_blue=indices_blue,indices_green=indices_green)
        print(f'[INFO] - Implementation output array shape: {output_array.shape}')

        if self.algo.endswith('_CI'):
            print(f'[INFO] Getting chl-a from Color Index (CI) algorithm to be used with low chl-a estimations')
            _,chl_ci = self.run_color_index(input_array,wl_array,iwl,output_as_implemented=True)
            output_array = self.run_merge_mbr_ci(output_array,chl_ci)
            if output_array is None:
                return

        if len(shape_impl)>2:
            final_output_shape = tuple(np.array(shape_impl)[:-1])
            output_array = np.reshape(output_array,final_output_shape)
        print(f'[INFO] - Final output array shape: {output_array.shape}')

        # if iwl_dest > iwl_orig:
        #     output_array = np.moveaxis(output_array,iwl_dest,iwl_orig)

        return output_array



    def run_from_mdb(self,path_mdb_file,options):
        ncout = Dataset(path_mdb_file,'r+')
        input_array = ncout.variables['satellite_Rrs'][:]
        wl_array = ncout.variables['satellite_bands'][:]
        if self.algo=='CI_NASA':
            CI,output_array = self.run_color_index(input_array,wl_array,iwl=1)
        elif self.algo=='CI_ESA':
            self.CI_METHOD = 'ESA'
            CI,output_array = self.run_color_index(input_array,wl_array,iwl=1)
        else:
            CI = None
            output_array = self.run_mbr_algorithm(input_array,wl_array=wl_array,iwl=1)
        if output_array is None:
            return

        sat_var_name = f'satellite_CHL_{self.algo}' if options['satellite_var_name'] is None else options['satellite_var_name']
        ncout = self.create_chl_variable(ncout,sat_var_name,output_array,('satellite_id','rows','columns'),options['overwrite_vars'])

        if CI is not None:
            ci_var_name = f'satellite_color_index'
            ncout = self.create_chl_variable(ncout, ci_var_name,CI, ('satellite_id', 'rows', 'columns'),options['overwrite_vars'],comment='Color index (CI)')


        if options['create_mu_variable']:
            mu_var_name = f'mu_satellite_CHL_{self.algo}' if options['mu_var_name'] is None else options['mu_var_name']
            stat, min_y, max_y, min_x, max_x = self.get_info_create_mu(options['mu_var_method'],output_array)
            if stat is None:
                return
            output_mu = None
            if stat=='median':
                output_mu = np.ma.median(output_array[:,min_y:max_y,min_x:max_x],axis=[1,2])
            elif stat=='avg':
                output_mu = np.ma.mean(output_array[:, min_y:max_y, min_x:max_x], axis=[1, 2])
            print(f'[INFO] Mu variable: {mu_var_name}')
            print(f'[INFO] Number of valid chl-a values in the mu variable: {np.ma.count(output_mu)}')
            if options['mu_var_ref'] is not None:
                array_ref = ncout.variables[options['mu_var_ref']][:]
                output_mu[array_ref.mask] = np.ma.masked
                print(f'[INFO] Number of valid chl-a values in the mu variable after masking based on a reference variable: {np.ma.count(output_mu)}')

            ncout = self.create_chl_variable(ncout,mu_var_name,output_mu,('mu_id',),options['overwrite_vars'])



        ncout.close()

        return

    def get_info_create_mu(self,mu_method,output_array):
        if len(mu_method) < 2:
            print(
                f'[ERROR] mu_var_method should at least contains two elements-> mu_method: metric,window_size. Default: median,3')
            return [None]*5
        stat = mu_method[0]
        stat_list = ['median', 'avg']
        if not stat in stat_list:
            print(
                f'[ERROR] Error in mu_var_method (comma-separated parameter). First element should be one of the list: {stat_list}')
            return [None]*5
        try:
            wsize = int(mu_method[1])
        except:
            print(f'[ERROR] Window size {mu_method[1]} is not a valid integer value')
            return [None]*5
        mid_y = int(np.round(output_array.shape[1] / 2))
        mid_x = int(np.round(output_array.shape[2] / 2))
        if wsize > min(mid_y, mid_x):
            print(f'[ERROR] Window size {wsize} should be lower than {min(mid_y, mid_x)}')
            return [None]*5
        min_y = mid_y - int(np.floor(wsize / 2))
        max_y = mid_y + int(np.ceil(wsize / 2))
        min_x = mid_x - int(np.floor(wsize / 2))
        max_x = mid_x + int(np.ceil(wsize / 2))

        return stat,min_y,max_y,min_x,max_x

    def create_chl_variable(self,ncout,var_name,output_array,dims_var,overwrite,comment=None):
        if not var_name in ncout.variables:
            var = ncout.createVariable(var_name, 'f4', dims_var, zlib=True, complevel=6,fill_value=-999.0)
            if comment is None:
                var.comment = f'Chl-a computed using MBR algorithm {self.algo}'
            else:
                var.comment = comment
            var[:] = output_array[:]
        else:
            old_output_array = ncout[var_name][:]
            diff = old_output_array / output_array
            are_equal = math.isclose(np.min(diff), 1, abs_tol=1e-6) and math.isclose(np.max(diff), 1, abs_tol=1e-6)

            if not are_equal:
                if overwrite:
                    print(f'[INFO] Overwriting variable: {var_name}')
                    ncout[var_name][:] = output_array[:]
                else:
                    print(
                        f'[WARNING] Variable {var_name} already exists in the file but with different values. To overwrite it, set overwrite_vars: True, and repeat.')
            else:
                print(f'[INFO] The variable {var_name} already exists in the file with the same expected values. Skipping...')

        return ncout

    def get_info_algorithm(self):
        file_algo = os.path.join(os.path.dirname(__init__.__file__),'MBR_Algorithms.csv')
        if not os.path.exists(file_algo):
            print(f'[ERROR] MBR algorithms data frame {file_algo} is not available')
            return None
        df_all = pd.read_csv(file_algo,sep=';')
        algo_list = df_all['Algo'].to_list()
        if self.algo not in algo_list:
            print(f'[ERROR] MBR algorithm {self.algo} is not in the list. Choose between: ')
            for algo_name in algo_list:
                print(f'   - {algo_name}')
            return None

        df_algo = df_all[df_all['Algo']==self.algo]
        return df_algo

