import numpy as np
import pandas
#import BSC_QAA.bsc_qaa_EUMETSAT as bsc_qaa
import warnings
import COMMON.common_functions as cfs
import COMMON.flag_functions as ffs

warnings.simplefilter(action='ignore', category=FutureWarning)


class QC_INSITU:

    def __init__(self, dataset):
        self.dataset = dataset
        self.wl_list = None
        self.indices_wl = None
        self.n_instruments = None
        self.spectral_convolution = None

        ##basic information
        self.bands_variable = 'insitu_original_bands'
        self.spectral_variable = 'insitu_Rrs'
        self.spectral_variable_unc = 'insitu_Rrs_unc'
        self.instrument_variable = 'insitu_instrument_id'
        self.time_max = 7200
        self.max_diff_wl = 5.0
        self.spectral_convolution_methods = ["nearest","interp","srf","band_shifting"]
        ##filters
        self.filter_invalid = None
        self.filter_flag = []
        self.filter_spectral_range = []
        self.filter_band_range = []

        self.spectral_stats = ['all', 'any', 'avg', 'median', 'std', 'iqr', 'min', 'max', 'CV']
        self.actions = ['keep','remove']
        self.filter_invalid_options = ['all', 'any']

        ##DIMENSIONS
        ##dimensions
        self.n_bands = None  ##defined for wl_list
        self.n_original_bands = None
        self.n_mu = None
        self.n_insitu_id = None
        # self.name = ''
        #
        # self.insitu_rrs = insitu_rrs
        # self.insitu_rrs_unc = None
        # self.insitu_bands = insitu_bands
        # self.n_instrument = -1
        # if len(self.insitu_bands.shape) == 2:
        #     self.n_instrument = self.insitu_bands.shape[0]
        # self.instrument_id_array = None
        #
        # self.nmu = self.insitu_rrs.shape[0]
        # self.nbands = self.insitu_rrs.shape[1]
        # self.nid = self.insitu_rrs.shape[2]
        #
        # self.maxdifwl = 10
        #
        # self.time_max = 7200
        #
        # self.band_stats = None
        #
        # self.wl_original = None
        # self.wl_list = None
        # self.wl_indices = None
        #
        #
        # self.thersholds = None
        #
        # self.apply_band_shift = False
        # self.apply_nir_correction = True
        # self.srf = None
        #
        # self.check_indices_by_mu = False
        # self.only_complete_spectra = True
        # self.pi_divided = False
        #
        # ##checking flags and other bands
        # self.ncdataset = None
        # self.check_flags = {}
        # self.check_th_other_bands = {}
        #
        # # for idendifying bands
        # self.msibands = self.get_msi_bands_dict()
        # self.olcibands = self.get_olci_bands_dict()
        #
        # self.insitu_valid_rrs = None
        #
        # # temporal for checking
        # self.insitu_time = None
        # self.sat_time = None



    def set_basic_info(self,options_config):
        self.bands_variable = options_config['bands_variable']
        self.spectral_variable = options_config['spectral_variable']
        self.spectral_variable_unc = options_config['spectral_variable_unc']
        self.instrument_variable = options_config['instrument_variable']

        self.max_diff_wl  = options_config['max_diff_wl']

        max_time_units = options_config['time_diff_max_units']
        max_time_value = options_config['time_diff_max']
        self.time_max = None
        if not isinstance(max_time_units, str):
            print(f'[ERROR][QC_INS] time_diff_max_units should be a string, not {max_time_units}')
        elif not isinstance(max_time_value, float):
            print(f'[ERROR][QC_INS] time_diff_max should be a float, not {max_time_value}')
        else:
            units = ['hours','minutes','seconds']
            u_factor = [3600, 60, 1]
            if not max_time_units in units:
                print(f'[ERROR][QC_INS] time_diff_max_units should be one of {units}')
            else:
                self.time_max = max_time_value * u_factor[units.index(max_time_units)]

        self.spectral_convolution = options_config['spectral_convolution']
        self.filter_invalid = options_config['filter_invalid']

    def set_wl_list(self,options_config):
        if not self.bands_variable in self.dataset.variables:
            print(f'[ERROR][QC_INS] {self.bands_variable} variable is not available in the NetCDF dataset')
            return

        original_bands = self.dataset.variables[self.bands_variable][:]
        if len(original_bands.shape)==1:
            original_bands = np.expand_dims(original_bands,axis=0)
            self.n_instruments = 1
        elif len(original_bands.shape)==2:
            self.n_instruments = original_bands.shape[0]
        else:
            print(f'[ERROR][QC_INS] Wavelength {self.bands_variable} variable should have 1 or 2 dimensions')
            return


        wl_list = options_config['wl_list']
        wl_min = options_config['wl_min']
        wl_max = options_config['wl_max']
        if wl_list is not None and (wl_min is not None or wl_max is not None):
            print(f'[WARNING] As wl_list is given, wl_min and wl_max values are not used for setting the wavelength list')

        if wl_list is None:
            if wl_min is not None or wl_max is not None:
                if wl_min is None:
                    wl_min = np.ma.min(original_bands)
                if wl_max  is None:
                    wl_max = np.ma.max(original_bands)
                if wl_max<wl_min:
                    print(f'[ERROR][QC_INS] wl_max ({wl_max}) must be greater or equal than wl_min ({wl_min})')

                valid_bands = np.where(np.sum((original_bands>=wl_min) & (original_bands<=wl_max),axis=0)>=1)[0]
                index_valid_min, index_valid_max = min(valid_bands), max(valid_bands)

                wl_list = original_bands[:,index_valid_min:index_valid_max+1]


                print(f'[INFO][QC_INS] wl_list set to {wl_list.shape[1]} bands between {np.ma.min(wl_list)} and {np.ma.max(wl_list)}')
            elif wl_min is None and wl_max is None:
                wl_list = original_bands
                print(f'[INFO][QC_INS] wl_list set to the original in situ bands with {wl_list.shape[1]} bands between {np.ma.min(wl_list)} and {np.ma.max(wl_list)}')

        else:##set the same wl_list for all the instruments
            n_wl_here = len(wl_list)
            wl_list = np.tile(np.ma.array(wl_list),self.n_instruments).reshape(self.n_instruments,n_wl_here)

        n_original_bands = original_bands.shape[1]
        n_list = wl_list.shape[1]

        self.wl_list = np.ma.masked_all((self.n_instruments,n_list))
        self.indices_wl = np.ma.masked_all((self.n_instruments,n_list),np.int32)
        self.n_bands = n_list
        self.n_original_bands = n_original_bands

        ##check maximum wavelength difference by instrument
        for i_instrument in range(self.n_instruments):
            original_bands_i = original_bands[i_instrument,:]
            if original_bands_i.mask.all():
                print(f'[WARNING][QC_INS] No bands are available for instrument {i_instrument}. It could cause problems if validation of spectra for this instrument is required. Skipping to the next instrument...')
                continue
            wl_list_i = wl_list[i_instrument,:]
            wl_list_m = np.repeat(wl_list_i,n_original_bands).reshape((n_list,n_original_bands))
            wl_original_m = np.tile(original_bands_i,n_list).reshape((n_list,n_original_bands))

            diff_wl = np.abs(wl_list_m-wl_original_m)
            min_diff_wl = np.min(diff_wl,axis=1)

            if np.max(min_diff_wl)>=self.max_diff_wl:
                print(f'[WARNING][QC_INS][START]=======================================================================')
                print(f'Some bands given in the parameter wl_list are not available as wavelength difference with the in situ original bands is greater than the allowed maximum of {self.max_diff_wl} nm')
                print(f'Options:')
                min_diff_wl_f = np.ma.filled(min_diff_wl,self.max_diff_wl)  ##no_valid_bands are defined only for non-masked
                no_valid_bands = wl_list_i[min_diff_wl_f>self.max_diff_wl]
                print(f' 1. Remove the following bands and try it again: ')
                for no_valid_band in no_valid_bands:
                    print(f'    * {no_valid_band} nm')
                print(f' 2. Modify the allowed  maximum difference using the max_diff_wl parameter in QC_INS and try it again')
                print(f' 3. Use these results, considering that in situ spectra could be incomplete for some of the required bands')
                print(f'[WARNING][QC_INS][STOP]========================================================================')

            indices_min_diff_wl = np.ma.array(np.argmin(diff_wl,axis=1),mask=min_diff_wl.mask)
            indices_min_diff_wl[min_diff_wl>self.max_diff_wl] = np.ma.masked
            self.wl_list[i_instrument,:] = wl_list_i[:]
            self.wl_list[i_instrument,indices_min_diff_wl.mask] = np.ma.masked
            self.indices_wl[i_instrument,:] = indices_min_diff_wl[:]


    def set_flag_filter(self,options_config, key_values = None):
        self.filter_flag = cfs.get_filter_list(options_config,'filter_flag_',key_values=key_values)


    def set_spectral_range_filter(self,options_config, key_values = None):
        self.filter_spectral_range = cfs.get_filter_list(options_config, 'filter_spectral_range_', key_values=key_values)
        #print(self.spectral_range_filter)

    def set_band_range_filter(self,options_config, key_values = None):
        self.filter_band_range = cfs.get_filter_list(options_config, 'filter_band_range_', key_values=key_values)
        #print(self.band_range_filter)

    def check_parameters(self):

        check_qc = self.check_spectral_variable(self.spectral_variable,self.bands_variable)

        if self.wl_list is None or self.indices_wl is None:
            check_qc = False

        if self.time_max is None:
            check_qc = False

        if self.n_mu is None:
            print(f'[ERROR][QC_INS] Problem retrieving the number of potential match-ups')
            check_qc = False
        if self.n_insitu_id is None:
            print(f'[ERROR][QC_INS] Problem retrieving the number of potential spectra')
            check_qc = False
        if self.n_bands is None:
            print(f'[ERROR][QC_INS] Problem retrieving the number of bands for wl_list')
            check_qc = False
        if self.n_original_bands is None:
            print(f'[ERROR][QC_INS] Problem retrieving the number of original bands')
            check_qc = False

        if not self.spectral_convolution  in self.spectral_convolution_methods:
            print(f'[ERROR][QC_INS] spectral_convolution {self.spectral_convolution} method is not available, it should be one of {self.spectral_convolution_methods}')
            check_qc = False

        if self.filter_invalid not in self.filter_invalid_options:
            print(f'[ERROR][QC_INS] filter_invalid {self.filter_invalid} should be one of {self.filter_invalid_options}')
            check_qc = False

        if self.spectral_variable_unc is not None and not self.check_spectral_variable(self.spectral_variable_unc,self.bands_variable):
            print(f'[ERROR][QC_INS] spectral_variable_unc {self.spectral_variable_unc} is not available, it should not be used in the quality control scripts')
            check_qc = False

        if self.n_instruments>1 and not self.check_non_spectral_variable(self.instrument_variable):
            print(f'[ERROR][QC_INS] instruments_variable {self.instrument_variable} is not available, but it is required as the number of instruments {self.n_instruments} > 1:wq')
            check_qc = False

        if len(self.filter_flag)>0:
            for idx in range(len(self.filter_flag)):
                if self.filter_flag[idx]['name_var'] is None and self.filter_flag[idx]['insitu_instrument'] is not None:
                    self.filter_flag[idx]['name_var'] = None ##not implemented, to get name_var from default insitu_instrument

                flag_list = self.filter_flag[idx]['flag_list']
                flag_list_valid = self.filter_flag[idx]['flag_list_valid']
                if flag_list is None and flag_list_valid is None:
                    print('NOT IMPLEMENTED, TO GET FLAG LISTS FROM insitu_instrument')
                if flag_list is None and flag_list_valid is None:
                    print(f'[ERROR][QC_INS] flag_list or flag_list_valid are required for filter_flag_{idx}')
                    check_qc = False
                else:
                    if flag_list is not None:
                        if '$' in flag_list:
                            print(f'[ERROR][QC_INS] Key flag \"$NONE$\" for empty flags is not allowed for flag_list, but only for flag_list_valid (filter_flag_{idx})')
                            check_qc = False
                        if not self.check_flag_list(self.filter_flag[idx]['name_var'],flag_list):
                            check_qc = False
                    if flag_list_valid is not None:
                        if not self.check_flag_list(self.filter_flag[idx]['name_var'],flag_list_valid):
                            check_qc = False
                    if flag_list_valid is not None and flag_list is not None:
                        irs = set(flag_list_valid).intersection(set(flag_list))
                        if len(irs)>0:
                            print(f'[ERROR][QC_INS] Flags \"{",".join(list(irs))}\" in flag_list and flag_list_valid are duplicated for filter_flag_{idx}')

        if len(self.filter_spectral_range)>0:
            for idx in range(len(self.filter_spectral_range)):
                if self.filter_spectral_range[idx]['name_var']=="$spectral_variable$":
                    self.filter_spectral_range[idx]['name_var'] = self.spectral_variable
                if self.filter_spectral_range[idx]['name_var_wl']=="$bands_variable$":
                    self.filter_spectral_range[idx]['name_var_wl'] = self.bands_variable
                check_b = self.check_spectral_variable(self.filter_spectral_range[idx]['name_var'],self.filter_spectral_range[idx]['name_var_wl'])
                if check_b:
                    wl_min_abs = np.min(self.dataset.variables[self.filter_spectral_range[idx]['name_var_wl']][:]) - self.max_diff_wl
                    wl_max_abs = np.max(self.dataset.variables[self.filter_spectral_range[idx]['name_var_wl']][:]) + self.max_diff_wl
                    wl_min_here = self.filter_spectral_range[idx]['wl_min']
                    wl_max_here = self.filter_spectral_range[idx]['wl_min']
                    if not isinstance(wl_min_here, float):
                        print(f'[ERROR][QC_INS] wl_min {wl_min_here} for filter_spectral_range_{idx} should be a float value')
                        check_qc = False
                    if not isinstance(wl_max_here, float):
                        print(f'[ERROR][QC_INS] wl_max {wl_max_here} for filter_spectral_range_{idx} should be a float value')
                        check_qc = False
                    if isinstance(wl_min_here, float) and isinstance(wl_max_here, float):
                        if wl_min_here > wl_max_here:
                            print(f'[ERROR][QC_INS] wl_max {wl_max_here} should be greater or equal to wl_min {wl_min_here} for filter_spectral_range_{idx}')
                            check_qc = False
                        else:
                            if wl_min_here < wl_min_abs or wl_min_here > wl_max_abs:
                                print(f'[ERROR][QC_INS] wl_min {wl_min_here} should be in the spectral range {wl_min_abs} - {wl_max_abs} for filter_spectral_range_{idx}')
                                check_qc = False
                            if wl_max_here < wl_min_abs or wl_max_here > wl_max_abs:
                                print(f'[ERROR][QC_INS] wl_max {wl_max_here} should be in the spectral range {wl_min_abs} - {wl_max_abs} for filter_spectral_range_{idx}')
                                check_qc = False
                else:
                    check_qc = False

                th_min_valid = isinstance(self.filter_spectral_range[idx]['th_min'], float)
                th_max_valid = isinstance(self.filter_spectral_range[idx]['th_max'], float)
                if th_min_valid or th_max_valid:
                    pass
                else:
                    if not th_min_valid:
                        print(f'[ERROR][QC_INS] Minimum threshold th_min {self.filter_spectral_range[idx]['th_min']} for filter_spectral_range_{idx} should be a float')
                        check_qc = False
                    if not th_max_valid:
                        print(f'[ERROR][QC_INS] Minimum threshold th_max {self.filter_spectral_range[idx]['th_max']} for filter_spectral_range_{idx} should be a float')
                        check_qc = False
                if not self.filter_spectral_range[idx]['spectral_stat'] in self.spectral_stats:
                    print(f'[ERROR][QC_INS] spectral_stat {self.filter_spectral_range[idx]['spectral_stat']} for filter_spectral_range_{idx} should one of {self.spectral_stats}')
                    check_qc = False
                if not self.filter_spectral_range[idx]['action'] in self.actions:
                    print(f'[ERROR][QC_INS] action {self.filter_spectral_range[idx]['action']} for filter_spectral_range_{idx} should one of {self.actions}')
                    check_qc = False

        if len(self.filter_band_range) > 0:
            for idx in range(len(self.filter_band_range)):
                check_b = self.check_non_spectral_variable(self.filter_band_range[idx]['name_var'])
                if not check_b:
                    print(f'[ERROR][QC_INS] Error with name_var in filter_band_range {idx}')
                    check_qc = False
                th_min_valid = isinstance(self.filter_band_range[idx]['th_min'], float)
                th_max_valid = isinstance(self.filter_band_range[idx]['th_max'], float)
                if th_min_valid or th_max_valid:
                    pass
                else:
                    if not th_min_valid:
                        print(f'[ERROR][QC_INS] Minimum threshold th_min {self.filter_band_range[idx]['th_min']} for filter_band_range_{idx} should be a float')
                        check_qc = False
                    if not th_max_valid:
                        print(f'[ERROR][QC_INS] Minimum threshold th_max {self.filter_band_range[idx]['th_max']} for filter_band_range_{idx} should be a float')
                        check_qc = False

                if not self.filter_band_range[idx]['action'] in self.actions:
                    print(
                        f'[ERROR][QC_INS]action {self.filter_band_range[idx]['action']} for filter_band_range_{idx} should one of {self.actions}')
                    check_qc = False

        return check_qc


    def get_wl_list_subset_from_main_wl_list(self,wl_list_here):
        wl_list_here = np.array(wl_list_here)
        n_list_here = wl_list_here.shape[0]
        wl_list_here = np.tile(wl_list_here, self.n_instruments).reshape(self.n_instruments, n_list_here)
        n_original_bands = self.wl_list.shape[1]
        wl_list_final = np.ma.masked_all((self.n_instruments, n_list_here), dtype=self.wl_list.dtype)
        indices_wl_final = np.ma.masked_all((self.n_instruments, n_list_here), dtype=self.indices_wl.dtype)
        for instrument_id in range(self.n_instruments):
            wl_list_here_m = np.repeat(wl_list_here[instrument_id, :], n_original_bands).reshape(
                (n_list_here, n_original_bands))
            wl_original_m = np.tile(self.wl_list[instrument_id], n_list_here).reshape((n_list_here, n_original_bands))
            diff_wl = np.abs(wl_list_here_m - wl_original_m)
            min_diff_wl = np.min(diff_wl, axis=1)
            indices_min_diff_wl = np.ma.array(np.argmin(diff_wl, axis=1), mask=min_diff_wl.mask)
            indices_min_diff_wl[min_diff_wl > self.max_diff_wl] = np.ma.masked
            wl_list_final[instrument_id, indices_min_diff_wl.mask == False] = self.wl_list[instrument_id, indices_min_diff_wl.compressed()]
            indices_wl_final[instrument_id,indices_min_diff_wl.mask == False] = self.indices_wl[instrument_id, indices_min_diff_wl.compressed()]

        return wl_list_final,indices_wl_final

    def get_all_spectra(self,name_var,wl_list_here = None,to_original_shape = True,limit_to=None):
        if wl_list_here is None:
            wl_list_here = self.wl_list
            indices_wl_here = self.indices_wl
        else:
            wl_list_here,indices_wl_here = self.get_wl_list_subset_from_main_wl_list(wl_list_here)





        # if np.sum(np.diff(self.indices_wl,axis=0))==0:##indices wl are equal for all the instruments
        #     indices_check = np.arange(self.dataset.variables[name_var].shape[1]).astype(np.int32).all()==self.indices_wl[0,:].all()
        #     if indices_check:##data don't need to subsetted using the indices
        #         data_dest = self.dataset.variables[name_var][:]
        #     else:
        #         data_dest = self.dataset.variables[name_var][:, self.indices_wl[0, :], :]
        #     return data_dest

        ##data are subsetted considering different instruments
        n_sat = self.dataset.variables[name_var].shape[0]
        n_ins = self.dataset.variables[name_var].shape[2]
        n_all = n_sat*n_ins
        n_bands = self.dataset.variables[name_var].shape[1]

        data_orig = np.reshape(np.moveaxis(self.dataset.variables[name_var][:],1,2),(n_all,n_bands))
        instrument_array = np.reshape(self.dataset.variables[self.instrument_variable][:],(n_all,))
        n_bands_dest = wl_list_here.shape[1]
        data_dest = np.ma.masked_all((n_all,n_bands_dest),dtype=data_orig.dtype)


        for instrument_id in range(self.n_instruments):
            indices_wl_instrument = indices_wl_here[instrument_id,:].compressed()
            wl_instrument_mask = indices_wl_here[instrument_id,:].mask
            is_sat_instrument = instrument_array==instrument_id
            if np.sum(is_sat_instrument)>0:
                data_dest_instrument = data_dest[is_sat_instrument,:]
                data_orig_instrument = data_orig[:,indices_wl_instrument]
                data_dest_instrument[:,wl_instrument_mask==False] = data_orig_instrument[is_sat_instrument,:]
                data_dest[is_sat_instrument,:] = data_dest_instrument[:,:]


        if limit_to is not None:
            if limit_to=='VALID_ALL':
                is_valid = np.ma.count(data_dest,axis=1)==data_dest.shape[1]
                data_dest = data_dest[is_valid,:]
            elif limit_to=='VALID_ANY':
                is_valid = np.ma.sum(data_dest, axis=1) > 0
                data_dest = data_dest[is_valid, :]

            return data_dest

        if to_original_shape:
            data_dest = np.moveaxis(np.reshape(data_dest,(n_sat,n_ins,n_bands_dest)),2,1)


        return data_dest

    def compute_validity(self):
        original_shape = self.dataset.variables[self.spectral_variable].shape
        spectral_data = self.get_all_spectra(self.spectral_variable,to_original_shape = False)
        if self.n_mu!=original_shape[0]:
            print(f'[ERROR] Number of match-ups {self.n_mu} does not agree with the shape of the spectral variable {self.spectral_variable}: {original_shape[0]} ')
            return
        if self.n_insitu_id!=original_shape[2]:
            print(f'[ERROR] Number of spectra {self.n_insitu_id} does not agree with the shape of the spectral variable {self.spectral_variable}: {original_shape.shape[1]} ')
            return
        n_all_spectra = self.n_mu * self.n_insitu_id
        if n_all_spectra != spectral_data.shape[0]:
            print(
                f'[ERROR] Total number of spectra {n_all_spectra} does not agree with the shape of the transformed spectral variable {self.spectral_variable}: {spectral_data.shape[0]} ')
            return

        if self.verbose:
            print(f'[INFO][QC_INS] Started spectral validity for {n_all_spectra} spectra ({self.n_mu} potential match-ups x {self.n_insitu_id} potential spectra by match-up)')

        invalid_filter_array = self.compute_invalid_filter(spectral_data)
        if self.verbose:
            print(f'[INFO][QC_INS] Number of non available or incomplete spectra: {np.sum(invalid_filter_array)}')

        mask_flag = self.compute_flag_mask_array()
        if self.verbose:
            print(f'[INFO][QC_INS] Number of filtered spectra using flag variables: {int(np.sum(mask_flag))}')


    def check_flag_list(self,name_variable,flag_list):
        if name_variable is None:
            print(f'[ERROR][QC_INS] name_variable is required for flagging filters, it could not be None')
            return False
        if not name_variable in self.dataset.variables:
            print(f'[ERROR][QC_INS] {name_variable} is not available in the dataset')
            return False
        if flag_list is None:
            print(f'[ERROR][QC_INS] {flag_list} is required for flagging filters, it could not be None')
            return False

        flag_list_var, flag_values_var = ffs.get_info_from_flag_variable(self.dataset.variables[name_variable],key_error='QC_INS')
        if flag_values_var is None or flag_list_var is None:
            return False
        flag_list_check = flag_list.copy()
        if '$NONE$' in flag_list_check:
            flag_list_check.remove('$NONE$')

        check = set(flag_list_check).issubset(set(flag_list_var))
        if not check:
            print(f'[ERROR][QC_INS] {flag_list_check} flags are not available in the variable {name_variable} flag list: {flag_list_var}')
        return check

    def check_spectral_variable(self,spectral_variable,bands_variable):
        check_qc = True
        if bands_variable is None or not bands_variable in self.dataset.variables:
            print(f'[ERROR][QC_INS] Band wavelengths variable  {bands_variable} is not available in the dataset. Choose among: ')
            print(f'[ERROR][QC_INS]{list(self.dataset.variables)}')
            return False
        if spectral_variable is None or not spectral_variable in self.dataset.variables:
            print(f'[ERROR][QC_INS] Spectral variable {spectral_variable} is not available in the dataset. Choose among:')
            print(f'[ERROR][QC_INS]{list(self.dataset.variables)}')
            return False

        if self.n_original_bands is None:
            original_bands = self.dataset.variables[self.bands_variable][:]
            if len(original_bands.shape) == 1:
                n_bands = original_bands.shape[0]
                self.n_original_bands = n_bands
            elif len(original_bands.shape) == 2:
                n_bands = original_bands.shape[1]
                self.n_original_bands = n_bands
            else:
                print(f'[ERROR][QC_INS] Wavelength {self.bands_variable} variable should have 1 or 2 dimensions')
                check_qc = False
                n_bands = -1
        else:
            n_bands = self.n_original_bands

        if len(self.dataset.variables[spectral_variable].shape) != 3:
            print(f'[ERROR][QC_INS] Spectral variable {self.spectral_variable} should have 3 dimensions: satellite_id,insitu_bands,insitu_id')
            check_qc = False
        else:
            n_bands_spectral = self.dataset.variables[spectral_variable].shape[1]
            n_mu = self.dataset.variables[spectral_variable].shape[0]
            n_id = self.dataset.variables[spectral_variable].shape[2]
            if n_bands_spectral != n_bands and n_bands>0:
                print(f'[ERROR][QC_INS] The number of bands in the spectral variable {spectral_variable} ({n_bands_spectral}) should be equal to the number of bands in the band wavelengths variable ({n_bands}) ')
                check_qc = False
            if self.n_mu is not None and self.n_mu!=n_mu:
                print(f'[ERROR][QC_INS] The number of match-ups in the spectral variable {spectral_variable} ({n_mu}) should be equal to the general number of match-ups ({self.n_mu}) ')
                check_qc = False
            if self.n_mu is None:
                self.n_mu = n_mu
            if self.n_insitu_id is not None and self.n_insitu_id!=n_id:
                print(f'[ERROR][QC_INS] The number of spectra in the spectral variable {spectral_variable} ({n_id}) should be equal to the general number of spectra ({self.n_insitu_id}) ')
                check_qc = False
            if self.n_insitu_id is None:
                self.n_insitu_id = n_id


        return check_qc

    def check_non_spectral_variable(self,name_var):
        check_qc = True
        if name_var is None or not name_var in self.dataset.variables:
            print(f'[ERROR][QC_INS] Non-spectral variable {name_var} is not available in the dataset. Choose among: ')
            print(f'[ERROR][QC_INS] {list(self.dataset.variables)}')
            return False

        if len(self.dataset.variables[name_var].shape) != 2:
            print(f'[ERROR][QC_INS] Spectral variable {name_var} should have 2 dimensions: satellite_id, insitu_id')
            check_qc = False
        n_mu = self.dataset.variables[name_var].shape[0]
        n_id = self.dataset.variables[name_var].shape[1]
        if self.n_mu is not None and self.n_mu != n_mu:
            print(f'[ERROR][QC_INS] The number of match-ups in the non-spectral variable {name_var} ({n_mu}) should be equal to the general number of match-ups ({self.n_mu}) ')
            check_qc = False
        if self.n_mu is None:
            self.n_mu = n_mu
        if self.n_insitu_id is not None and self.n_insitu_id != n_id:
            print(f'[ERROR][QC_INS] The number of spectra in the non-spectral variable {name_var} ({n_id}) should be equal to the general number of spectra ({self.n_insitu_id}) ')
            check_qc = False
        if self.n_insitu_id is None:
            self.n_insitu_id = n_id



        return check_qc


    def compute_invalid_filter(self,spectral_data):

        invalid_filter = None
        if self.filter_invalid.lower() == 'any':
            invalid_filter = np.where(np.ma.count_masked(spectral_data,axis=1)>0,1,0)

        elif self.filter_invalid.lower() == 'all':
            invalid_filter = np.where(np.ma.count_masked(spectral_data, axis=1) == spectral_data.shape[1], 1, 0)

        return invalid_filter

    def compute_flag_mask_array(self):
        n_all = self.n_mu * self.n_insitu_id
        flag_mask = np.zeros((n_all,), dtype=np.uint64)
        for idx in range(len(self.filter_flag)):
            flag_mask_here = self.compute_flag_mask_array_impl(self.filter_flag[idx])
            if self.verbose:
                print(f'[INFO][QC_INS] Number of removed spectra using filter_flag_{idx} ({self.filter_flag[idx]["name_var"]}): {np.sum(flag_mask_here)}')
            if flag_mask_here is not None:
                flag_mask = flag_mask + flag_mask_here
        flag_mask[flag_mask > 0] = 1

        return flag_mask

    def compute_flag_mask_array_impl(self, f_flag):
        flag_array_all = self.dataset.variables[f_flag['name_var']][:].flatten()
        if np.issubdtype(flag_array_all.dtype, np.floating): ##floating points are not allowed
            flag_array_all = flag_array_all.astype('uint64')

        flag_array = flag_array_all.compressed()
        
        fw = ffs.start_flag_work_from_variable(self.dataset.variables[f_flag['name_var']])
        if fw is None:
            print(f'[ERROR][QC_INS] Flag work object for variable {f_flag["name_var"]} could not be started')
            return np.zeros(flag_array.shape, dtype=np.int8)
        if self.verbose:
            print(f'[INFO][QC_INS] Starting FlagWork object with variable {f_flag["name_var"]}. Data type: {fw.dType}')

        ##invalid
        mask_array = None
        if f_flag['flag_list'] is not None:
            mask_array = fw.mask(flag_array,f_flag['flag_list'])
            mask_array[mask_array >= 1] = 1
        if f_flag['flag_list_valid'] is not None:
            f_list = f_flag['flag_list_valid'].copy()

            if '$NONE$' in f_list:
                mask_array_v = np.where(flag_array==0,0,1)
                if len(f_list)>1:
                    f_list.remove('$NONE$')
                    m_array = fw.mask(flag_array,f_list)
                    mask_array_v[m_array>=1] = 0
            else:
                m_array = fw.mask(flag_array, f_list)
                mask_array_v = np.where(m_array>=1,0,1)
            if mask_array is None:
                mask_array = mask_array_v
            else:
                mask_array = np.where(np.logical_and(mask_array_v==0,ask_array==0,mask_array),0,1)

        mask_array_all = np.ma.masked_all(flag_array_all.shape,dtype=mask_array.dtype)
        mask_array_all[flag_array_all.mask==False] = mask_array[:]

        return mask_array_all


    # def check_validity(self):
    #
    #
    #     nmu = self.insitu_rrs.shape[0]
    #     nid = self.insitu_rrs.shape[2]
    #     nspectra = nmu * nid
    #     nbands = self.insitu_rrs.shape[1]
    #
    #     self.insitu_valid_rrs = np.ma.zeros((nmu,nid))
    #
    #     ##check first masked values or nan,inf,-inf, etc
    #     insitu_rrs_check = np.ma.masked_invalid(self.insitu_rrs)
    #     insitu_rrs_check = np.moveaxis(insitu_rrs_check,1,2)
    #     insitu_rrs_check = np.reshape(insitu_rrs_check,(nspectra,nbands))
    #
    #     # print('check validity line 83',insitu_rrs_check.shape)
    #     # print(self.wl_list,self.wl_indices,self.instrument_id_array is None)
    #     #
    #     # if self.instrument_id_array is not None:
    #     #     instrument_id_array = np.reshape(self.instrument_id_array,(nspectra,))
    #     #     print(instrument_id_array.shape)
    #
    #
    #     #insitu_rrs_check_sum = np.ma.count_masked(insitu_rrs_check, axis=1) > 0
    #     #insitu_rrs_check_sum = np.ma.count_masked(insitu_rrs_check, axis=1) > 3 ##for LAIT
    #     #insitu_rrs_check_sum = np.ma.count_masked(insitu_rrs_check, axis=1) > 1 ##for TRIT
    #     insitu_rrs_check_sum = np.ma.count_masked(insitu_rrs_check,axis=1)>2 ##for VEIT with instrument change
    #
    #
    #
    #     insitu_rrs_check_sum = np.reshape(insitu_rrs_check_sum,(nmu,nid))
    #     self.insitu_valid_rrs[insitu_rrs_check_sum==False]=self.insitu_valid_rrs[insitu_rrs_check_sum==False]+1
    #
    #
    #
    #     nconditions = 1
    #
    #     # checking flag: TO BE IMPLEMENTED
    #     if len(self.check_flags) > 0:
    #         for flag_name in self.check_flags:
    #             flag_var = self.check_flags[flag_name]['variable']
    #             oflag = self.check_flags[flag_name]['oflags']
    #             if oflag is not None:
    #                 nconditions = nconditions + 1
    #                 flag_value = flag_var[:] #np.array(flag_var[index_mu, insitu_id])
    #                 flag_list = self.check_flags[flag_name]['flag_list']
    #                 m = oflag.Mask(flag_value, flag_list)
    #                 if self.check_flags[flag_name]['remove_spectra'] and m > 0:
    #                     check = False
    #                 if not self.check_flags[flag_name]['remove_spectra'] and m == 0:
    #                     check = False
    #
    #
    #     # checking threshold other bands
    #     if len(self.check_th_other_bands) > 0:
    #         for band_name in self.check_th_other_bands:
    #             nconditions = nconditions+1
    #             val_here = self.check_th_other_bands[band_name]['variable'][:]
    #             th_type = self.check_th_other_bands[band_name]['th_type']
    #             val_min = self.check_th_other_bands[band_name]['value_min']
    #             val_max = self.check_th_other_bands[band_name]['value_max']
    #             is_angle = self.check_th_other_bands[band_name]['isangle']
    #
    #
    #             if val_max >= val_min:
    #                 check_condition = np.logical_and(val_min <= val_here,val_here<=val_max)
    #             elif val_max < val_min and is_angle:
    #                 check_condition = val_here >= val_min or val_here <= val_max
    #             else:
    #                 check_condition = None
    #             if check_condition is not None:
    #                 if th_type == 'keep':
    #                     self.insitu_valid_rrs[check_condition==True]=self.insitu_valid_rrs[check_condition==True]+1
    #                 if th_type == 'remove':
    #                     self.insitu_valid_rrs[check_condition==False]=self.insitu_valid_rrs[check_condition==False]+1
    #
    #     # checking rrs thresholds
    #     if self.thersholds is not None:
    #         for idx in range(len(self.wl_list)):
    #             wl = self.wl_list[idx]
    #             #val = np.ma.squeeze(self.insitu_rrs[:,idx,:])
    #             rrs_min = None
    #             rrs_max = None
    #             wls = str(wl)
    #             if self.thersholds[wls]['min_th']['apply'] or self.thersholds[wls]['max_th']['apply']:
    #                 min_wl_real = self.thersholds[wls]['min_wl']
    #                 max_wl_real = self.thersholds[wls]['max_wl']
    #                 rrs_min,rrs_max = self.get_min_values_for_spectral_range(min_wl_real,max_wl_real)
    #
    #
    #             if self.thersholds[wls]['min_th']['apply'] and rrs_min is not None:
    #                 nconditions = nconditions + 1
    #                 check_condition = rrs_min>=self.thersholds[wls]['min_th']['value']
    #                 print(f'[INFO][QC_INS] Condition: {nconditions} Wavelength: {wls} Lower thershold: {self.thersholds[wls]['min_th']['value']} Number of bad spectra: {np.ma.count(check_condition) - np.sum(check_condition)}')
    #                 self.insitu_valid_rrs[check_condition == True] = self.insitu_valid_rrs[check_condition == True] + 1
    #
    #             if self.thersholds[wls]['max_th']['apply'] and rrs_max is not None:
    #                 nconditions = nconditions + 1
    #                 check_condition = rrs_max<=self.thersholds[wls]['max_th']['value']
    #                 print(f'[INFO][QC_INS] Condition: {nconditions} Wavelength: {wls} Upper thershold: {self.thersholds[wls]['max_th']['value']} Number of bad spectra: {np.ma.count(check_condition)-np.sum(check_condition)}')
    #                 self.insitu_valid_rrs[check_condition == True] = self.insitu_valid_rrs[check_condition == True] + 1
    #
    #
    #
    #
    #     print(f'[INFO][QC_INS]->Number of conditions analysed: {nconditions}')
    #     self.insitu_valid_rrs = np.where(self.insitu_valid_rrs==nconditions,True,False)
    #     print(f'[INFO][QC_INS]->Number of valid spectra: {np.ma.sum(self.insitu_valid_rrs)}')
    #
    # def get_spectra_rrs(self, insitu_id, valid_sat_id=None):
    #     return self.get_spectra_impl(self.insitu_rrs,insitu_id,valid_sat_id)
    #
    # def get_spectra_unc(self, insitu_id, valid_sat_id=None):
    #     return self.get_spectra_impl(self.insitu_rrs_unc,insitu_id,valid_sat_id)
    #
    # def get_spectra_impl(self,insitu_rrs_impl,insitu_id,valid_sat_id=None):
    #     insitu_rrs_impl = insitu_rrs_impl[:]
    #     insitu_rrs_all = np.moveaxis(insitu_rrs_impl, 1, 2)
    #     if valid_sat_id is not None:
    #         insitu_rrs_nearest = insitu_rrs_all[valid_sat_id==1, insitu_id, :]
    #     else:
    #         insitu_rrs_nearest = insitu_rrs_all[np.arange(insitu_rrs_all.shape[0]).astype(np.int32), insitu_id, :]
    #     if len(self.wl_indices.shape) == 1:
    #         insitu_rrs_nearest = insitu_rrs_nearest[:, self.wl_indices]
    #         return insitu_rrs_nearest
    #
    #     if self.wl_indices.shape[0]==1:
    #         insitu_rrs_nearest = insitu_rrs_nearest[:, self.wl_indices[0]]
    #         return insitu_rrs_nearest
    #
    #
    #     indices_equal = (self.wl_indices[0]==self.wl_indices).all()
    #     if indices_equal:
    #         insitu_rrs_nearest = insitu_rrs_nearest[:,self.wl_indices[0]]
    #         return insitu_rrs_nearest
    #     else:
    #         if self.instrument_id_array is None:
    #             print(f'[ERROR][QC_INS] {self.wl_indices.shape[0]} different wl arrays  for different instruments are given, but instrument id array is not available')
    #             return None
    #         if valid_sat_id is not None:
    #             instrument_id_array = self.instrument_id_array[valid_sat_id == 1, insitu_id]
    #         else:
    #             instrument_id_array = self.instrument_id_array[np.arange(self.instrument_id_array.shape[0]).astype(np.int32), insitu_id]
    #         instruments_here = np.unique(instrument_id_array).compressed()
    #         insitu_rrs_end = np.ma.masked_all((insitu_rrs_nearest.shape[0],self.wl_indices.shape[1]))
    #         for instrument_id in instruments_here:
    #             insitu_rrs_nearest_here = insitu_rrs_nearest[instrument_id_array==instrument_id, :]
    #             insitu_rrs_end[instrument_id_array==instrument_id,:] = insitu_rrs_nearest_here[:,self.wl_indices[instrument_id]]
    #
    #         return insitu_rrs_end
    #
    #
    #
    #
    # def get_olci_bands_dict(self):
    #     olcibands = {
    #         'Oa01': {'wl': 400, 'apply': False},
    #         'Oa02': {'wl': 412.5, 'apply': False},
    #         'Oa03': {'wl': 442.5, 'apply': False},
    #         'Oa04': {'wl': 490, 'apply': False},
    #         'Oa05': {'wl': 510, 'apply': False},
    #         'Oa06': {'wl': 560, 'apply': False},
    #         'Oa07': {'wl': 620, 'apply': False},
    #         'Oa08': {'wl': 665, 'apply': False},
    #         'Oa09': {'wl': 673.75, 'apply': False},
    #         'Oa10': {'wl': 681.25, 'apply': False},
    #         'Oa11': {'wl': 708.75, 'apply': False},
    #         'Oa12': {'wl': 753.75, 'apply': False},
    #         'Oa13': {'wl': 761.25, 'apply': False},
    #         'Oa14': {'wl': 764.375, 'apply': False},
    #         'Oa15': {'wl': 767.5, 'apply': False},
    #         'Oa16': {'wl': 778.75, 'apply': False},
    #         'Oa17': {'wl': 865, 'apply': False},
    #         'Oa18': {'wl': 885, 'apply': False},
    #         'Oa19': {'wl': 900, 'apply': False},
    #         'Oa20': {'wl': 940, 'apply': False},
    #         'Oa21': {'wl': 1020, 'apply': False},
    #     }
    #     return olcibands
    #
    # def get_msi_bands_dict(self):
    #     msibands = {
    #         'B1': {'wl': 443, 'apply': False},
    #         'B2': {'wl': 490, 'apply': False},
    #         'B3': {'wl': 560, 'apply': False},
    #         'B4': {'wl': 665, 'apply': False},
    #         'B5': {'wl': 705, 'apply': False},
    #         'B6': {'wl': 740, 'apply': False},
    #         'B7': {'wl': 783, 'apply': False},
    #         'B8': {'wl': 842, 'apply': False},
    #         'B8A': {'wl': 865, 'apply': False},
    #         'B9': {'wl': 945, 'apply': False},
    #         'B10': {'wl': 1375, 'apply': False},
    #         'B11': {'wl': 1610, 'apply': False},
    #         'B12': {'wl': 2190, 'apply': False}
    #     }
    #     return msibands
    #
    # ##Method to make the subset of the spectra
    # def set_wllist_min_max(self, wlmin, wlmax):
    #     self.wl_list = []
    #     iwl_min = 0
    #     iwl_max = self.nbands
    #     started = False
    #     for iwl in range(len(self.insitu_bands)):
    #         wl = self.insitu_bands[iwl]
    #         if wlmin is not None and wlmax is not None and wlmax > wlmin:
    #             if wl >= wlmin and not started:
    #                 iwl_min = iwl
    #                 started = True
    #             if started and wl <= wlmax:
    #                 iwl_max = iwl
    #                 self.wl_list.append(wl)
    #         else:
    #             self.wl_list.append(wl)
    #
    #     if iwl_max < self.nbands:
    #         iwl_max = iwl_max + 1
    #
    #     self.wl_indices = list(range(iwl_min, iwl_max))
    #
    # def set_wllist_using_wlref(self, wlreflist):
    #
    #     self.wl_list = wlreflist
    #     nwl_max = len(wlreflist)
    #     n_check = 1 if self.n_instrument == -1 else self.n_instrument
    #     self.wl_indices = np.ma.masked_all((n_check, nwl_max)).astype(np.int32)
    #     self.msibands = self.get_msi_bands_dict()
    #     self.olcibands = self.get_olci_bands_dict()
    #
    #     for idx in range(n_check):
    #         for iwl, wl in enumerate(self.wl_list):
    #             index, wl_index = self.get_insitu_index_instrument(wl, idx)
    #             if index >= 0:
    #                 self.wl_indices[idx, iwl] = index
    #                 for b in self.msibands:
    #                     if abs(wl - self.msibands[b]['wl']) < 10:
    #                         self.msibands[b]['apply'] = True
    #                 for b in self.olcibands:
    #                     if abs(wl - self.olcibands[b]['wl']) < 2:
    #                         self.olcibands[b]['apply'] = True
    #
    #
    # def start_quality_control(self):
    #     self.thersholds = {}
    #     for wl in self.wl_list:
    #         wls = str(wl)
    #         self.thersholds[wls] = {
    #             'min_th': {
    #                 'apply': False,
    #                 'value': 0
    #             },
    #             'max_th': {
    #                 'apply': False,
    #                 'value': 0
    #             },
    #             'min_wl': None,
    #             'max_wl': None
    #         }
    #
    # def set_thershold(self, valuemin, valuemax, wlmin, wlmax):
    #     if self.thersholds is None:
    #         self.start_quality_control()
    #
    #     for wl in self.wl_list:
    #         if wl < wlmin or wl > wlmax:
    #             continue
    #         wls = str(wl)
    #         if valuemin is not None:
    #             self.thersholds[wls]['min_th']['apply'] = True
    #             self.thersholds[wls]['min_th']['value'] = valuemin
    #         if valuemax is not None:
    #             self.thersholds[wls]['max_th']['apply'] = True
    #             self.thersholds[wls]['max_th']['value'] = valuemax
    #         self.thersholds[wls]['min_wl'] = wlmin
    #         self.thersholds[wls]['max_wl'] = wlmax
    #
    # def add_flag_expression(self, flag_band, flag_list, remove_spectra):
    #     if self.ncdataset is None:
    #         return
    #
    #     if flag_band not in self.ncdataset.variables:
    #         return
    #
    #     flag_variable = self.ncdataset.variables[flag_band]
    #
    #     if flag_list == 'ALL':  # APPLY ALL THE FLAGS
    #         flag_list = flag_variable.flag_meanings.split()
    #     else:
    #         flag_list = [x.strip() for x in flag_list.split(',')]
    #
    #     oflags = None
    #     try:
    #         from COMMON.Class_Flags_OLCI import Class_Flags_OLCI
    #         flag_values = [np.uint32(x.strip()) for x in flag_variable.flag_mask.split(',')]
    #         oflags = Class_Flags_OLCI(flag_values, flag_variable.flag_meanings)
    #     except:
    #         print(f'[WARNING] Flag class could not be defined for variable: {flag_band}')
    #         pass
    #
    #     self.check_flags[flag_band] = {
    #         'variable': flag_variable,
    #         'flag_list': flag_list,
    #         'remove_spectra': remove_spectra,
    #         'oflags': oflags
    #     }
    #
    # def add_other_band_thersholds(self, band_name, type_th, value_min, value_max, isangle):
    #     if self.ncdataset is None:
    #         return
    #
    #     if band_name not in self.ncdataset.variables:
    #         return
    #
    #     band_variable = self.ncdataset.variables[band_name]
    #
    #     self.check_th_other_bands[band_name] = {
    #         'variable': band_variable,
    #         'th_type': type_th,
    #         'value_min': value_min,
    #         'value_max': value_max,
    #         'isangle': isangle
    #     }
    #
    #     # print('---->',self.check_th_other_bands)
    #
    # def check_validity_spectrum(self, rrs_values, index_mu, insitu_id):
    #
    #     if rrs_values is None:
    #         return False
    #
    #     if len(rrs_values)==1:
    #         rrs_values = rrs_values[0]
    #
    #     # if np.sum(np.isnan(rrs_values)) > 0:
    #     #     return False
    #
    #     if self.only_complete_spectra and (rrs_values.count() != len(rrs_values)):
    #         return False
    #     check = True
    #
    #     # checking flag
    #     if len(self.check_flags) > 0:
    #         for flag_name in self.check_flags:
    #             flag_var = self.check_flags[flag_name]['variable']
    #             oflag = self.check_flags[flag_name]['oflags']
    #             if oflag is not None:
    #                 flag_value = np.array(flag_var[index_mu, insitu_id])
    #                 flag_list = self.check_flags[flag_name]['flag_list']
    #                 m = oflag.Mask(flag_value, flag_list)
    #                 if self.check_flags[flag_name]['remove_spectra'] and m > 0:
    #                     check = False
    #                 if not self.check_flags[flag_name]['remove_spectra'] and m == 0:
    #                     check = False
    #
    #     if not check:
    #         return False
    #
    #     # checking threshold other bands
    #     if len(self.check_th_other_bands) > 0:
    #
    #         for band_name in self.check_th_other_bands:
    #
    #             var_here = self.check_th_other_bands[band_name]['variable']
    #             val_here = np.array(var_here[index_mu, insitu_id])
    #
    #             th_type = self.check_th_other_bands[band_name]['th_type']
    #             val_min = self.check_th_other_bands[band_name]['value_min']
    #             val_max = self.check_th_other_bands[band_name]['value_max']
    #             is_angle = self.check_th_other_bands[band_name]['isangle']
    #
    #             check_condition = False
    #             if val_max >= val_min:
    #                 check_condition = val_min <= val_here <= val_max
    #             elif val_max < val_min and is_angle:
    #                 check_condition = val_here >= val_min or val_here <= val_max
    #
    #             if th_type == 'keep' and not check_condition:
    #                 check = False
    #
    #             if th_type == 'remove' and check_condition:
    #                 check = False
    #
    #     if not check:
    #         return False
    #
    #     # checking thresholds
    #     if self.thersholds is not None:
    #
    #         for idx in range(len(self.wl_list)):
    #             wl = self.wl_list[idx]
    #             val = rrs_values[idx]
    #             wls = str(wl)
    #
    #
    #             if self.thersholds[wls]['min_th']['apply'] and val < self.thersholds[wls]['min_th']['value']:
    #                 check = False
    #                 break
    #             if self.thersholds[wls]['max_th']['apply'] and val > self.thersholds[wls]['max_th']['value']:
    #                 check = False
    #                 break
    #
    #     return check
    #
    # def check_validity_spectra_mu(self, index_mu):
    #     val_array = np.zeros((self.nid))
    #     for idx in range(self.nid):
    #         rrs_values, indices, valid_bands = self.get_spectrum_for_mu_and_index_insitu(index_mu, idx)
    #         valid = self.check_validity_spectrum(rrs_values, index_mu, idx)
    #         if valid:
    #             val_array[idx] = 1
    #             break
    #     return val_array
    #
    #
    # def get_insitu_index_instrument(self, wl, idx):
    #     if self.n_instrument == -1:
    #         all_wl = self.insitu_bands[:]
    #     else:
    #         all_wl = self.insitu_bands[idx, :]
    #     dif_wl = np.abs(all_wl - wl)
    #     index = np.argmin(dif_wl)
    #     wl_index = all_wl[index]
    #     if np.ma.is_masked(wl_index):
    #         index = -1
    #     else:
    #         if dif_wl[index] > self.maxdifwl:
    #             index = -1
    #     return index, wl_index
    #
    # def get_min_values_for_spectral_range(self,wl_min,wl_max):
    #     insitu_bands = self.insitu_bands
    #     if len(self.insitu_bands.shape)==1:
    #         insitu_bands = np.expand_dims(self.insitu_bands,0)
    #
    #     insitu_rrs_min = None
    #     insitu_rrs_max = None
    #     if self.instrument_id_array is None:
    #         all_wl = insitu_bands[0, :]
    #         indices = np.where(np.logical_and(all_wl>=wl_min,all_wl<=wl_max))
    #         insitu_rrs_indices = self.insitu_rrs[(np.arange(self.insitu_rrs.shape[0]),indices[0],np.arange(self.insitu_rrs.shape[2]))]
    #         insitu_rrs_min = np.ma.min(insitu_rrs_indices,axis=1)
    #         insitu_rrs_max = np.ma.max(insitu_rrs_indices,axis=1)
    #
    #     return insitu_rrs_min,insitu_rrs_max
    #
    # def get_insitu_index(self, wl):
    #     all_wl = self.insitu_bands[:]
    #     dif_wl = np.abs(all_wl - wl)
    #     index = np.argmin(dif_wl)
    #     wl_index = all_wl[index]
    #     if dif_wl[index] > self.maxdifwl:
    #         index = -1
    #     return index, wl_index
    #
    # def get_insitu_index_mu(self, imu, wl):
    #     dif_ref = self.maxdifwl
    #     index = -1
    #     for iband in range(self.nbands):
    #         dif = np.abs(wl - self.insitu_bands[iband])
    #         valarray = ma.array(self.insitu_rrs[imu, iband, 0])
    #         if dif <= dif_ref and not valarray.mask:
    #             index = iband
    #             dif_ref = dif
    #     return index
    #
    # def get_insitu_indices_mu(self, imu):
    #     indices = []
    #     valid_bands = [False] * len(self.wl_list)
    #     for idx in range(len(self.wl_list)):
    #         wl = self.wl_list[idx]
    #         index = self.get_insitu_index_mu(imu, wl)
    #         if index >= 0:
    #             valid_bands[idx] = True
    #             indices.append(index)
    #         else:
    #             if self.only_complete_spectra:
    #                 return None, None
    #             else:
    #                 first_index_valid = -1
    #                 for iband in range(self.nbands):
    #                     valarray = ma.array(self.insitu_rrs[imu, iband, 0])
    #                     if not valarray.mask:
    #                         first_index_valid = iband
    #                         break
    #                 indices.append(
    #                     first_index_valid)  # note that this index is not used, rrs value in invalid bands will be masked in a later step
    #     return indices, valid_bands
    #
    # # retrieve stats
    # # see compute_band_statistics for correct values of stat
    # def get_stat_spectra(self, wlmin, wlmax, stat):
    #     if self.band_stats is None:
    #         self.compute_good_spectra_statistics()
    #     stat_spectra = []
    #     for wls in self.band_stats:
    #         wl = float(wls)
    #         if wlmin is not None and wl < wlmin:
    #             continue
    #         if wlmax is not None and wl > wlmax:
    #             continue
    #         stat_spectra.append(self.band_stats[wls][stat])
    #
    #     return np.array(stat_spectra)
    #
    # def get_finalspectrum_mu(self, index_mu, dif_time_array, exact_wl_array, wl_ref):
    #
    #     time_condition = False
    #     spectrum_complete = False
    #     valid_values = False
    #     id_min_time = -1
    #     rrs_values = None
    #     #spectra_with_time_condition = False
    #
    #     ngood = np.ma.count(dif_time_array)
    #     if ngood<dif_time_array.shape[0]:
    #         dif_time_good = dif_time_array[~dif_time_array.mask]
    #     else:
    #         dif_time_good = dif_time_array
    #
    #     if ngood == 0:
    #         return id_min_time, time_condition, valid_values, spectrum_complete, rrs_values
    #
    #     indices_good = np.argsort(dif_time_good)
    #
    #     time_condition_array = dif_time_good[indices_good]<self.time_max
    #
    #     spectra_array = np.squeeze(self.insitu_rrs[index_mu, :, :])
    #     spectra_validity = np.squeeze(self.insitu_valid_rrs[index_mu, :])
    #
    #     if ngood < dif_time_array.shape[0]:
    #         spectra_array = spectra_array[:,~dif_time_array.mask]
    #         spectra_validity = spectra_validity[~dif_time_array.mask]
    #
    #     spectra_array = spectra_array[:,indices_good]
    #     spectra_validity = spectra_validity[indices_good]
    #
    #     if self.instrument_id_array is not None:
    #         instrument_valid = self.instrument_id_array[~dif_time_array.mask]
    #         instrument_valid = instrument_valid[indices_good]
    #         ins_here = np.unique(instrument_valid)
    #         if len(ins_here)==1:
    #             if self.wl_indices.shape[0]==1:
    #                 indices = self.wl_indices[0,:]
    #             else:
    #                 indices = self.wl_indices[ins_here[0],:]
    #             rrs_values = spectra_array[indices,:]
    #         else:
    #             n_bands = len(self.wl_indices[ins_here[0],:])
    #             n_here = spectra_array.shape[1]
    #             rrs_values = np.ma.masked_all((n_bands,n_here))
    #             for ihere in range(n_here):
    #                 indices = self.wl_indices[instrument_valid[ihere],:]
    #                 rrs_values[:,ihere]= spectra_array[indices,ihere]
    #     else:
    #         indices = self.wl_indices[0, :]
    #         rrs_values = spectra_array[indices,:]
    #
    #
    #
    #     n_valid_bands = np.sum(np.invert(ma.getmaskarray(rrs_values)),axis=0)
    #     complete_spectra = n_valid_bands == len(self.wl_list)
    #
    #     all_valid = np.logical_and(time_condition_array,spectra_validity)
    #     if self.only_complete_spectra:
    #         all_valid = np.logical_and(all_valid,complete_spectra)
    #
    #     if np.count_nonzero(all_valid)>0:
    #         ins_valid = np.where(all_valid)[0][0]
    #         id_min_time = indices_good[ins_valid]
    #         rrs_values = rrs_values[:, ins_valid]
    #         insitu_valid = True
    #     else:
    #         ins_valid = 0
    #         id_min_time = indices_good[ins_valid]
    #         rrs_values = None
    #         insitu_valid = False
    #
    #     # if rrs_values is not None and rrs_values[-1]>0.005:
    #     #     print(self.wl_list)
    #     #     print('me llega aqui',index_mu,id_min_time,rrs_values[-1],rrs_values[-1]*1000)
    #
    #     return id_min_time,time_condition_array[ins_valid],insitu_valid,complete_spectra[ins_valid],rrs_values
    #
    #
    #     # if index_mu==0:
    #     #     print(n_valid_bands)
    #     #     print(complete_spectra)
    #     # print('==================')
    #
    #
    #
    #
    #     # for idx in indices_good:
    #     #     id_min_time = idx
    #     #     time_dif = dif_time_array[idx]
    #     #     time_condition = time_dif < self.time_max
    #     #
    #     #
    #     #     if time_condition:
    #     #         spectra_with_time_condition = True
    #     #         rrs_values, indices, valid_bands = self.get_spectrum_for_mu_and_index_insitu(index_mu, idx)
    #     #         valid_bands_array = np.array(valid_bands, dtype=bool)
    #     #         rrs_values = np.ma.masked_where(valid_bands_array == False, rrs_values)
    #     #         valid_values = self.check_validity_spectrum(rrs_values, index_mu, idx)
    #     #
    #     #         spectrum_complete = np.sum(valid_bands_array) == len(self.wl_list)
    #     #         if valid_values and self.apply_band_shift and exact_wl_array is not None and wl_ref is not None:
    #     #             if len(exact_wl_array.shape) == 1:
    #     #                 exact_wl = exact_wl_array[indices]
    #     #             else:
    #     #                 exact_wl = exact_wl_array[indices, id_min_time]
    #     #             rrs_values = bsc_qaa.bsc_qaa(rrs_values, exact_wl, wl_ref)
    #     #         if valid_values:
    #     #             break
    #     #
    #     # if not valid_values and spectra_with_time_condition:
    #     #     time_condition = True
    #
    #     #return id_min_time, time_condition, valid_values, spectrum_complete, rrs_values
    #
    # ##OPERATIONS WITH GOOD SPECTRA
    # def compute_good_spectra_statistics(self):
    #     print(f'[INFO][QC_INS] Getting valid spectra...')
    #     spectra_res = self.get_all_good_spectra()
    #     band_stats = {}
    #     for index in range(len(self.wl_list)):
    #         wl = self.wl_list[index]
    #         wls = str(wl)
    #         print(f'[INFO][QC_INS] Computing statistics for band: {wls} ({index})')
    #         allvalid = spectra_res[:, index]
    #         band_stats[wls] = {
    #             'nvalid': 0,
    #             'min_val': 0,
    #             'max_val': 0,
    #             'mean_val': 0,
    #             'std_val': 0,
    #             'median': 0,
    #             'p25': 0,
    #             'p75': 0
    #         }
    #         # print(f'[INFO][QC_INS]Computing statistics...')
    #         if allvalid is not None:
    #             band_stats[wls]['nvalid'] = len(allvalid)
    #             band_stats[wls]['min_val'] = np.min(allvalid)
    #             band_stats[wls]['max_val'] = np.max(allvalid)
    #             band_stats[wls]['mean_val'] = np.mean(allvalid)
    #             band_stats[wls]['std_val'] = np.std(allvalid)
    #             band_stats[wls]['median'] = np.median(allvalid)
    #             band_stats[wls]['p25'] = np.percentile(allvalid, 25)
    #             band_stats[wls]['p75'] = np.percentile(allvalid, 75)
    #
    #     self.band_stats = band_stats
    #     return self.band_stats
    #
    # def get_spectrum_for_mu_and_index_insitu_unc(self, index_mu, index_insitu):
    #     if self.insitu_rrs_unc is None:
    #         return None, None, None
    #
    #     if index_mu < 0 or index_mu >= self.nmu:
    #         return None, None, None
    #     spectra_unc = ma.array(self.insitu_rrs_unc[index_mu, :, index_insitu])
    #     return self.get_spectrum_for_mu_and_index_insitu_impl(spectra_unc, index_mu, index_insitu)
    #
    # def get_spectrum_for_mu_and_index_insitu(self, index_mu, index_insitu):
    #     if index_mu < 0 or index_mu >= self.nmu:
    #         return None, None, None
    #     spectra = ma.array(self.insitu_rrs[index_mu, :, index_insitu])
    #     return self.get_spectrum_for_mu_and_index_insitu_impl(spectra, index_mu, index_insitu)
    #
    # def get_spectrum_for_mu_and_index_insitu_impl(self, spectra, index_mu,index_insitu):
    #     if self.check_indices_by_mu:
    #         indices, valid_bands = self.get_insitu_indices_mu(index_mu)
    #         rrs_values = spectra[indices]
    #         rrs_values[valid_bands == False] = np.ma.masked
    #     else:
    #         if self.n_instrument==-1:
    #             indices = self.wl_indices
    #         elif self.n_instrument==1:
    #             indices = self.wl_indices[0, :]
    #         else:
    #             iins = 0
    #             if self.instrument_id_array is not None:
    #                 iins = self.instrument_id_array[index_insitu]
    #                 #iins = self.instrument_id_array[index_mu,index_insitu]
    #             if np.ma.is_masked(iins):
    #                 return None,None,None
    #             indices = self.wl_indices[iins,:]
    #
    #
    #
    #
    #         rrs_values = spectra[indices]
    #         valid_bands = np.invert(ma.getmaskarray(rrs_values)).tolist()
    #
    #
    #     # implementation of spectral response function here
    #
    #     ##csv sentinel2 srf
    #     if self.srf is not None and self.srf.endswith('.csv'):
    #         import pandas as pd
    #         S2srf = pd.read_csv(self.srf)
    #         bands = self.insitu_bands[:]
    #         df = pd.DataFrame(data=spectra, index=bands)
    #         nnT = self.reindex_and_interpolate(df, S2srf["SR_WL"].values)
    #         rrs_values = []
    #         for a in np.arange(1, S2srf.shape[1]):
    #             band = str(S2srf.columns[a]).split('_')[-1]
    #             if self.msibands[band]['apply']:
    #                 rrs_values.append(
    #                     nnT.mul(np.array(S2srf.iloc[:, a]), axis=0).sum(axis=0) / (sum(np.array(S2srf.iloc[:, a]))))
    #         rrs_values = np.ma.array(rrs_values)
    #         rrs_values = rrs_values.flatten()
    #
    #     ##olci3 nc4 srf
    #     if self.srf is not None and self.srf.endswith('.nc4'):
    #         from netCDF4 import Dataset
    #         import pandas as pd
    #         bands = self.insitu_bands[:]
    #         olci_bands = list(self.olcibands.keys())
    #         df = pd.DataFrame(data=spectra, index=bands)
    #
    #         dataset = Dataset(self.srf)
    #
    #         wl_values = np.array(dataset.variables['mean_spectral_response_function_wavelength'])
    #         msrf = np.array(dataset.variables['mean_spectral_response_function'])
    #         nwl = wl_values.shape[0]
    #
    #         rrs_values = []
    #         for iwl in range(nwl):
    #             band = olci_bands[iwl]
    #             if self.olcibands[band]['apply']:
    #                 nnT = self.reindex_and_interpolate(df, wl_values[iwl, :])
    #                 nnT_new = nnT.mul(np.array(msrf[iwl, :]), axis=0).sum(axis=0) / np.sum(np.array(msrf[iwl, :]))
    #                 rrs_values.append(nnT_new.loc[0])
    #
    #         dataset.close()
    #
    #     return rrs_values, indices, valid_bands
    #
    # def reindex_and_interpolate(self, df, new_index):
    #     return df.reindex(df.index | new_index).interpolate(method='index', limit_direction='both').loc[new_index]
    #
    # # ngood is only for checking (assing -1 for not using it)
    # def get_good_spectrum_for_mu(self, index_mu, id_min_time, ngood):
    #
    #     spectra, indices, valid_bands = self.get_all_good_spectra_for_mu(index_mu)
    #
    #     rrs_values = None
    #     if ngood == -1:
    #         ngood = spectra.shape[0]
    #     if spectra is not None and spectra.shape[0] == ngood:
    #         rrs_values = spectra[id_min_time, :]
    #
    #     return rrs_values, indices, valid_bands
    #
    # def get_all_good_spectra_for_mu(self, index_mu):
    #     if index_mu < 0 or index_mu >= self.nmu:
    #         return None
    #
    #     spectra = np.transpose(ma.array(self.insitu_rrs[index_mu]))
    #
    #     vspectra = np.all(spectra == -999, axis=1)
    #
    #     # ids_good = np.where(vspectra == False)[0]
    #
    #     spectra_here = spectra[vspectra == False]
    #
    #     if self.check_indices_by_mu:
    #         indices, valid_bands = self.get_insitu_indices_mu(index_mu)
    #     else:
    #         indices = self.wl_indices
    #         valid_bands = [True] * len(self.wl_indices)
    #
    #     if indices is None:
    #         spectra_valid = None
    #     else:
    #         spectra_valid = spectra_here[:, indices]
    #
    #     return spectra_valid, indices, valid_bands
    #
    # def get_all_good_spectra(self, check_indices):
    #     ntotal = self.nmu * self.nid
    #
    #     nbands_here = len(self.wl_indices)
    #
    #     valid_spectra = ma.zeros((ntotal, nbands_here))
    #     index_valid = 0
    #     for index_mu in range(self.nmu):
    #         spectra_here = self.get_good_spectra_for_mu(index_mu, check_indices)
    #         nspectra = spectra_here.shape[0]
    #         iini = index_valid
    #         ifin = index_valid + nspectra
    #         # valid_spectra[iini:ifin, :] = spectra_here[0:nspectra, self.iwl_min:self.iwl_max]
    #         valid_spectra[iini:ifin, :] = spectra_here[:, :]
    #         index_valid = ifin
    #
    #     spectra_rec = valid_spectra[0:index_valid, :]
    #
    #     return spectra_rec
