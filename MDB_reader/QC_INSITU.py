import os.path

import numpy as np
import numpy.ma as ma
import pandas

import BSC_QAA.bsc_qaa_EUMETSAT as bsc_qaa
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


class QC_INSITU:

    def __init__(self, insitu_rrs, insitu_bands):
        self.name = ''

        self.insitu_rrs = insitu_rrs
        self.insitu_rrs_unc = None
        self.insitu_bands = insitu_bands
        self.n_instrument = -1
        if len(self.insitu_bands.shape) == 2:
            self.n_instrument = self.insitu_bands.shape[0]
        self.instrument_id_array = None

        self.nmu = self.insitu_rrs.shape[0]
        self.nbands = self.insitu_rrs.shape[1]
        self.nid = self.insitu_rrs.shape[2]

        self.maxdifwl = 10

        self.time_max = 7200

        self.band_stats = None

        self.wl_list = None
        self.wl_indices = None


        self.thersholds = None

        self.apply_band_shift = False
        self.apply_nir_correction = True
        self.srf = None

        self.check_indices_by_mu = False
        self.only_complete_spectra = True
        self.pi_divided = False

        ##checking flags and other bands
        self.ncdataset = None
        self.check_flags = {}
        self.check_th_other_bands = {}

        # for idendifying bands
        self.msibands = self.get_msi_bands_dict()
        self.olcibands = self.get_olci_bands_dict()

        self.insitu_valid_rrs = None

    def check_validity(self):
        self.insitu_valid_rrs = np.ma.zeros((self.insitu_rrs.shape[0],self.insitu_rrs.shape[2]))
        nconditions = 0

        # checking flag: TO BE IMPLEMENTED
        if len(self.check_flags) > 0:
            for flag_name in self.check_flags:
                flag_var = self.check_flags[flag_name]['variable']
                oflag = self.check_flags[flag_name]['oflags']
                if oflag is not None:
                    nconditions = nconditions + 1
                    flag_value = flag_var[:] #np.array(flag_var[index_mu, insitu_id])
                    flag_list = self.check_flags[flag_name]['flag_list']
                    m = oflag.Mask(flag_value, flag_list)
                    if self.check_flags[flag_name]['remove_spectra'] and m > 0:
                        check = False
                    if not self.check_flags[flag_name]['remove_spectra'] and m == 0:
                        check = False

        # checking threshold other bands
        if len(self.check_th_other_bands) > 0:
            for band_name in self.check_th_other_bands:
                nconditions = nconditions+1
                val_here = self.check_th_other_bands[band_name]['variable'][:]
                th_type = self.check_th_other_bands[band_name]['th_type']
                val_min = self.check_th_other_bands[band_name]['value_min']
                val_max = self.check_th_other_bands[band_name]['value_max']
                is_angle = self.check_th_other_bands[band_name]['isangle']


                if val_max >= val_min:
                    check_condition = np.logical_and(val_min <= val_here,val_here<=val_max)
                elif val_max < val_min and is_angle:
                    check_condition = val_here >= val_min or val_here <= val_max

                if th_type == 'keep':
                    self.insitu_valid_rrs[check_condition==True]=self.insitu_valid_rrs[check_condition==True]+1
                if th_type == 'remove':
                    self.insitu_valid_rrs[check_condition==False]=self.insitu_valid_rrs[check_condition==False]+1

        # checking rrs thresholds
        if self.thersholds is not None:

            for idx in range(len(self.wl_list)):
                wl = self.wl_list[idx]
                val = np.ma.squeeze(self.insitu_rrs[:,idx,:])
                wls = str(wl)
                if self.thersholds[wls]['min_th']['apply']:
                    nconditions = nconditions + 1
                    check_condition = val>=self.thersholds[wls]['min_th']['value']
                    self.insitu_valid_rrs[check_condition == True] = self.insitu_valid_rrs[check_condition == True] + 1
                if self.thersholds[wls]['max_th']['apply']:
                    nconditions = nconditions + 1
                    check_condition = val<=self.thersholds[wls]['max_th']['value']
                    self.insitu_valid_rrs[check_condition == True] = self.insitu_valid_rrs[check_condition == True] + 1

                # if self.thersholds[wls]['min_th']['apply'] and val < self.thersholds[wls]['min_th']['value']:
                #     check = False
                #     break
                # if self.thersholds[wls]['max_th']['apply'] and val > self.thersholds[wls]['max_th']['value']:
                #     check = False
                #     break

        print(f'[INFO]->Number of conditions analysed: {nconditions}')
        self.insitu_valid_rrs = np.where(self.insitu_valid_rrs==nconditions,True,False)
        print(f'[INFO]->Number of valid spectra: {np.ma.sum(self.insitu_valid_rrs)}')




    def get_olci_bands_dict(self):
        olcibands = {
            'Oa01': {'wl': 400, 'apply': False},
            'Oa02': {'wl': 412.5, 'apply': False},
            'Oa03': {'wl': 442.5, 'apply': False},
            'Oa04': {'wl': 490, 'apply': False},
            'Oa05': {'wl': 510, 'apply': False},
            'Oa06': {'wl': 560, 'apply': False},
            'Oa07': {'wl': 620, 'apply': False},
            'Oa08': {'wl': 665, 'apply': False},
            'Oa09': {'wl': 673.75, 'apply': False},
            'Oa10': {'wl': 681.25, 'apply': False},
            'Oa11': {'wl': 708.75, 'apply': False},
            'Oa12': {'wl': 753.75, 'apply': False},
            'Oa13': {'wl': 761.25, 'apply': False},
            'Oa14': {'wl': 764.375, 'apply': False},
            'Oa15': {'wl': 767.5, 'apply': False},
            'Oa16': {'wl': 778.75, 'apply': False},
            'Oa17': {'wl': 865, 'apply': False},
            'Oa18': {'wl': 885, 'apply': False},
            'Oa19': {'wl': 900, 'apply': False},
            'Oa20': {'wl': 940, 'apply': False},
            'Oa21': {'wl': 1020, 'apply': False},
        }
        return olcibands

    def get_msi_bands_dict(self):
        msibands = {
            'B1': {'wl': 443, 'apply': False},
            'B2': {'wl': 490, 'apply': False},
            'B3': {'wl': 560, 'apply': False},
            'B4': {'wl': 665, 'apply': False},
            'B5': {'wl': 705, 'apply': False},
            'B6': {'wl': 740, 'apply': False},
            'B7': {'wl': 783, 'apply': False},
            'B8': {'wl': 842, 'apply': False},
            'B8A': {'wl': 865, 'apply': False},
            'B9': {'wl': 945, 'apply': False},
            'B10': {'wl': 1375, 'apply': False},
            'B11': {'wl': 1610, 'apply': False},
            'B12': {'wl': 2190, 'apply': False}
        }
        return msibands

    ##Method to make the subset of the spectra
    def set_wllist_min_max(self, wlmin, wlmax):
        self.wl_list = []
        iwl_min = 0
        iwl_max = self.nbands
        started = False
        for iwl in range(len(self.insitu_bands)):
            wl = self.insitu_bands[iwl]
            if wlmin is not None and wlmax is not None and wlmax > wlmin:
                if wl >= wlmin and not started:
                    iwl_min = iwl
                    started = True
                if started and wl <= wlmax:
                    iwl_max = iwl
                    self.wl_list.append(wl)
            else:
                self.wl_list.append(wl)

        if iwl_max < self.nbands:
            iwl_max = iwl_max + 1

        self.wl_indices = list(range(iwl_min, iwl_max))

    def set_wllist_using_wlref(self, wlreflist):

        self.wl_list = wlreflist
        nwl_max = len(wlreflist)
        n_check = 1 if self.n_instrument == -1 else self.n_instrument
        self.wl_indices = np.ma.masked_all((n_check, nwl_max)).astype(np.int32)
        self.msibands = self.get_msi_bands_dict()
        self.olcibands = self.get_olci_bands_dict()

        for idx in range(n_check):
            for iwl, wl in enumerate(self.wl_list):
                index, wl_index = self.get_insitu_index_instrument(wl, idx)
                if index > 0:
                    self.wl_indices[idx, iwl] = index
                    for b in self.msibands:
                        if abs(wl - self.msibands[b]['wl']) < 10:
                            self.msibands[b]['apply'] = True
                    for b in self.olcibands:
                        if abs(wl - self.olcibands[b]['wl']) < 2:
                            self.olcibands[b]['apply'] = True

    def start_quality_control(self):
        self.thersholds = {}
        for wl in self.wl_list:
            wls = str(wl)
            self.thersholds[wls] = {
                'min_th': {
                    'apply': False,
                    'value': 0
                },
                'max_th': {
                    'apply': False,
                    'value': 0
                }
            }

    def set_thershold(self, valuemin, valuemax, wlmin, wlmax):
        if self.thersholds is None:
            self.start_quality_control()
        for wl in self.wl_list:
            if wl < wlmin or wl > wlmax:
                continue
            wls = str(wl)
            if valuemin is not None:
                self.thersholds[wls]['min_th']['apply'] = True
                self.thersholds[wls]['min_th']['value'] = valuemin
            if valuemax is not None:
                self.thersholds[wls]['max_th']['apply'] = True
                self.thersholds[wls]['max_th']['value'] = valuemax

    def add_flag_expression(self, flag_band, flag_list, remove_spectra):
        if self.ncdataset is None:
            return

        if flag_band not in self.ncdataset.variables:
            return

        flag_variable = self.ncdataset.variables[flag_band]

        if flag_list == 'ALL':  # APPLY ALL THE FLAGS
            flag_list = flag_variable.flag_meanings.split()
        else:
            flag_list = [x.strip() for x in flag_list.split(',')]

        oflags = None
        try:
            from COMMON.Class_Flags_OLCI import Class_Flags_OLCI
            flag_values = [np.uint32(x.strip()) for x in flag_variable.flag_mask.split(',')]
            oflags = Class_Flags_OLCI(flag_values, flag_variable.flag_meanings)
        except:
            print(f'[WARNING] Flag class could not be defined for variable: {flag_band}')
            pass

        self.check_flags[flag_band] = {
            'variable': flag_variable,
            'flag_list': flag_list,
            'remove_spectra': remove_spectra,
            'oflags': oflags
        }

    def add_other_band_thersholds(self, band_name, type_th, value_min, value_max, isangle):
        if self.ncdataset is None:
            return

        if band_name not in self.ncdataset.variables:
            return

        band_variable = self.ncdataset.variables[band_name]

        self.check_th_other_bands[band_name] = {
            'variable': band_variable,
            'th_type': type_th,
            'value_min': value_min,
            'value_max': value_max,
            'isangle': isangle
        }

        # print('---->',self.check_th_other_bands)

    def check_validity_spectrum(self, rrs_values, index_mu, insitu_id):

        if rrs_values is None:
            return False

        if len(rrs_values)==1:
            rrs_values = rrs_values[0]

        # if np.sum(np.isnan(rrs_values)) > 0:
        #     return False

        if self.only_complete_spectra and (rrs_values.count() != len(rrs_values)):
            return False
        check = True

        # checking flag
        if len(self.check_flags) > 0:
            for flag_name in self.check_flags:
                flag_var = self.check_flags[flag_name]['variable']
                oflag = self.check_flags[flag_name]['oflags']
                if oflag is not None:
                    flag_value = np.array(flag_var[index_mu, insitu_id])
                    flag_list = self.check_flags[flag_name]['flag_list']
                    m = oflag.Mask(flag_value, flag_list)
                    if self.check_flags[flag_name]['remove_spectra'] and m > 0:
                        check = False
                    if not self.check_flags[flag_name]['remove_spectra'] and m == 0:
                        check = False

        if not check:
            return False

        # checking threshold other bands
        if len(self.check_th_other_bands) > 0:

            for band_name in self.check_th_other_bands:

                var_here = self.check_th_other_bands[band_name]['variable']
                val_here = np.array(var_here[index_mu, insitu_id])

                th_type = self.check_th_other_bands[band_name]['th_type']
                val_min = self.check_th_other_bands[band_name]['value_min']
                val_max = self.check_th_other_bands[band_name]['value_max']
                is_angle = self.check_th_other_bands[band_name]['isangle']

                check_condition = False
                if val_max >= val_min:
                    check_condition = val_min <= val_here <= val_max
                elif val_max < val_min and is_angle:
                    check_condition = val_here >= val_min or val_here <= val_max

                if th_type == 'keep' and not check_condition:
                    check = False

                if th_type == 'remove' and check_condition:
                    check = False

        if not check:
            return False

        # checking thresholds
        if self.thersholds is not None:

            for idx in range(len(self.wl_list)):
                wl = self.wl_list[idx]
                val = rrs_values[idx]
                wls = str(wl)


                if self.thersholds[wls]['min_th']['apply'] and val < self.thersholds[wls]['min_th']['value']:
                    check = False
                    break
                if self.thersholds[wls]['max_th']['apply'] and val > self.thersholds[wls]['max_th']['value']:
                    check = False
                    break

        return check

    def check_validity_spectra_mu(self, index_mu):
        val_array = np.zeros((self.nid))
        for idx in range(self.nid):
            rrs_values, indices, valid_bands = self.get_spectrum_for_mu_and_index_insitu(index_mu, idx)
            valid = self.check_validity_spectrum(rrs_values, index_mu, idx)
            if valid:
                val_array[idx] = 1
                break
        return val_array

    def get_insitu_index_instrument(self, wl, idx):
        if self.n_instrument == -1:
            all_wl = self.insitu_bands[:]
        else:
            all_wl = self.insitu_bands[idx, :]
        dif_wl = np.abs(all_wl - wl)
        index = np.argmin(dif_wl)
        wl_index = all_wl[index]
        if np.ma.is_masked(wl_index):
            index = -1
        else:
            if dif_wl[index] > self.maxdifwl:
                index = -1
        return index, wl_index

    def get_insitu_index(self, wl):
        all_wl = self.insitu_bands[:]
        dif_wl = np.abs(all_wl - wl)
        index = np.argmin(dif_wl)
        wl_index = all_wl[index]
        if dif_wl[index] > self.maxdifwl:
            index = -1
        return index, wl_index

    def get_insitu_index_mu(self, imu, wl):
        dif_ref = self.maxdifwl
        index = -1
        for iband in range(self.nbands):
            dif = np.abs(wl - self.insitu_bands[iband])
            valarray = ma.array(self.insitu_rrs[imu, iband, 0])
            if dif <= dif_ref and not valarray.mask:
                index = iband
                dif_ref = dif
        return index

    def get_insitu_indices_mu(self, imu):
        indices = []
        valid_bands = [False] * len(self.wl_list)
        for idx in range(len(self.wl_list)):
            wl = self.wl_list[idx]
            index = self.get_insitu_index_mu(imu, wl)
            if index >= 0:
                valid_bands[idx] = True
                indices.append(index)
            else:
                if self.only_complete_spectra:
                    return None, None
                else:
                    first_index_valid = -1
                    for iband in range(self.nbands):
                        valarray = ma.array(self.insitu_rrs[imu, iband, 0])
                        if not valarray.mask:
                            first_index_valid = iband
                            break
                    indices.append(
                        first_index_valid)  # note that this index is not used, rrs value in invalid bands will be masked in a later step
        return indices, valid_bands

    # retrieve stats
    # see compute_band_statistics for correct values of stat
    def get_stat_spectra(self, wlmin, wlmax, stat):
        if self.band_stats is None:
            self.compute_good_spectra_statistics()
        stat_spectra = []
        for wls in self.band_stats:
            wl = float(wls)
            if wlmin is not None and wl < wlmin:
                continue
            if wlmax is not None and wl > wlmax:
                continue
            stat_spectra.append(self.band_stats[wls][stat])

        return np.array(stat_spectra)



    def get_finalspectrum_mu(self, index_mu, dif_time_array, exact_wl_array, wl_ref):

        time_condition = False
        spectrum_complete = False
        valid_values = False
        id_min_time = -1
        rrs_values = None
        #spectra_with_time_condition = False

        ngood = np.ma.count(dif_time_array)
        if ngood<dif_time_array.shape[0]:
            dif_time_good = dif_time_array[~dif_time_array.mask]
        else:
            dif_time_good = dif_time_array

        if ngood == 0:
            return id_min_time, time_condition, valid_values, spectrum_complete, rrs_values

        indices_good = np.argsort(dif_time_good)

        time_condition_array = dif_time_good[indices_good]<self.time_max

        spectra_array = np.squeeze(self.insitu_rrs[index_mu, :, :])
        spectra_validity = np.squeeze(self.insitu_valid_rrs[index_mu, :])
        if ngood < dif_time_array.shape[0]:
            spectra_array = spectra_array[:,~dif_time_array.mask]
            spectra_validity = spectra_validity[~dif_time_array.mask]

        spectra_array = spectra_array[:,indices_good]
        spectra_validity = spectra_validity[indices_good]

        indices = self.wl_indices[0, :]

        rrs_values = spectra_array[indices,:]

        n_valid_bands = np.sum(np.invert(ma.getmaskarray(rrs_values)),axis=0)
        complete_spectra = n_valid_bands == len(self.wl_list)

        all_valid = np.logical_and(time_condition_array,spectra_validity)
        if self.only_complete_spectra:
            all_valid = np.logical_and(all_valid,complete_spectra)

        if np.count_nonzero(all_valid)>0:
            ins_valid = np.where(all_valid)[0][0]
            tal = np.where(all_valid)[0]
            ital = np.min(tal)
            if ital!=ins_valid:
                print('ojo',ital,ins_valid)
            id_min_time = indices_good[ins_valid]
            rrs_values = rrs_values[:, ins_valid]
            insitu_valid = True
        else:
            ins_valid = 0
            id_min_time = indices_good[ins_valid]
            rrs_values = None
            insitu_valid = False



        return id_min_time,time_condition_array[ins_valid],insitu_valid,complete_spectra[ins_valid],rrs_values


        # if index_mu==0:
        #     print(n_valid_bands)
        #     print(complete_spectra)
        # print('==================')




        # for idx in indices_good:
        #     id_min_time = idx
        #     time_dif = dif_time_array[idx]
        #     time_condition = time_dif < self.time_max
        #
        #
        #     if time_condition:
        #         spectra_with_time_condition = True
        #         rrs_values, indices, valid_bands = self.get_spectrum_for_mu_and_index_insitu(index_mu, idx)
        #         valid_bands_array = np.array(valid_bands, dtype=bool)
        #         rrs_values = np.ma.masked_where(valid_bands_array == False, rrs_values)
        #         valid_values = self.check_validity_spectrum(rrs_values, index_mu, idx)
        #
        #         spectrum_complete = np.sum(valid_bands_array) == len(self.wl_list)
        #         if valid_values and self.apply_band_shift and exact_wl_array is not None and wl_ref is not None:
        #             if len(exact_wl_array.shape) == 1:
        #                 exact_wl = exact_wl_array[indices]
        #             else:
        #                 exact_wl = exact_wl_array[indices, id_min_time]
        #             rrs_values = bsc_qaa.bsc_qaa(rrs_values, exact_wl, wl_ref)
        #         if valid_values:
        #             break
        #
        # if not valid_values and spectra_with_time_condition:
        #     time_condition = True

        #return id_min_time, time_condition, valid_values, spectrum_complete, rrs_values

    ##OPERATIONS WITH GOOD SPECTRA
    def compute_good_spectra_statistics(self):
        print(f'[INFO] Getting valid spectra...')
        spectra_res = self.get_all_good_spectra()
        band_stats = {}
        for index in range(len(self.wl_list)):
            wl = self.wl_list[index]
            wls = str(wl)
            print(f'[INFO] Computing statistics for band: {wls} ({index})')
            allvalid = spectra_res[:, index]
            band_stats[wls] = {
                'nvalid': 0,
                'min_val': 0,
                'max_val': 0,
                'mean_val': 0,
                'std_val': 0,
                'median': 0,
                'p25': 0,
                'p75': 0
            }
            # print(f'[INFO]Computing statistics...')
            if allvalid is not None:
                band_stats[wls]['nvalid'] = len(allvalid)
                band_stats[wls]['min_val'] = np.min(allvalid)
                band_stats[wls]['max_val'] = np.max(allvalid)
                band_stats[wls]['mean_val'] = np.mean(allvalid)
                band_stats[wls]['std_val'] = np.std(allvalid)
                band_stats[wls]['median'] = np.median(allvalid)
                band_stats[wls]['p25'] = np.percentile(allvalid, 25)
                band_stats[wls]['p75'] = np.percentile(allvalid, 75)

        self.band_stats = band_stats
        return self.band_stats

    def get_spectrum_for_mu_and_index_insitu_unc(self, index_mu, index_insitu):
        if self.insitu_rrs_unc is None:
            return None, None, None

        if index_mu < 0 or index_mu >= self.nmu:
            return None, None, None
        spectra_unc = ma.array(self.insitu_rrs_unc[index_mu, :, index_insitu])
        return self.get_spectrum_for_mu_and_index_insitu_impl(spectra_unc, index_mu, index_insitu)

    def get_spectrum_for_mu_and_index_insitu(self, index_mu, index_insitu):
        if index_mu < 0 or index_mu >= self.nmu:
            return None, None, None
        spectra = ma.array(self.insitu_rrs[index_mu, :, index_insitu])
        return self.get_spectrum_for_mu_and_index_insitu_impl(spectra, index_mu, index_insitu)

    def get_spectrum_for_mu_and_index_insitu_impl(self, spectra, index_mu,index_insitu):
        if self.check_indices_by_mu:
            indices, valid_bands = self.get_insitu_indices_mu(index_mu)
            rrs_values = spectra[indices]
            rrs_values[valid_bands == False] = np.ma.masked

        else:
            if self.n_instrument==-1:
                indices = self.wl_indices
            elif self.n_instrument==1:
                indices = self.wl_indices[0, :]
            else:
                iins = 0
                if self.instrument_id_array is not None:
                    iins = self.instrument_id_array[index_mu,index_insitu]
                if np.ma.is_masked(iins):
                    return None,None,None
                indices = self.wl_indices[iins,:]



            rrs_values = spectra[indices]
            valid_bands = np.invert(ma.getmaskarray(rrs_values)).tolist()


        # implementation of spectral response function here

        ##csv sentinel2 srf
        if self.srf is not None and self.srf.endswith('.csv'):
            import pandas as pd
            S2srf = pd.read_csv(self.srf)
            bands = self.insitu_bands[:]
            df = pd.DataFrame(data=spectra, index=bands)
            nnT = self.reindex_and_interpolate(df, S2srf["SR_WL"].values)
            rrs_values = []
            for a in np.arange(1, S2srf.shape[1]):
                band = str(S2srf.columns[a]).split('_')[-1]
                if self.msibands[band]['apply']:
                    rrs_values.append(
                        nnT.mul(np.array(S2srf.iloc[:, a]), axis=0).sum(axis=0) / (sum(np.array(S2srf.iloc[:, a]))))
            rrs_values = np.ma.array(rrs_values)
            rrs_values = rrs_values.flatten()

        ##olci3 nc4 srf
        if self.srf is not None and self.srf.endswith('.nc4'):
            from netCDF4 import Dataset
            import pandas as pd
            bands = self.insitu_bands[:]
            olci_bands = list(self.olcibands.keys())
            df = pd.DataFrame(data=spectra, index=bands)

            dataset = Dataset(self.srf)

            wl_values = np.array(dataset.variables['mean_spectral_response_function_wavelength'])
            msrf = np.array(dataset.variables['mean_spectral_response_function'])
            nwl = wl_values.shape[0]

            rrs_values = []
            for iwl in range(nwl):
                band = olci_bands[iwl]
                if self.olcibands[band]['apply']:
                    nnT = self.reindex_and_interpolate(df, wl_values[iwl, :])
                    nnT_new = nnT.mul(np.array(msrf[iwl, :]), axis=0).sum(axis=0) / np.sum(np.array(msrf[iwl, :]))
                    rrs_values.append(nnT_new.loc[0])

            dataset.close()

        return rrs_values, indices, valid_bands

    def reindex_and_interpolate(self, df, new_index):
        return df.reindex(df.index | new_index).interpolate(method='index', limit_direction='both').loc[new_index]

    # ngood is only for checking (assing -1 for not using it)
    def get_good_spectrum_for_mu(self, index_mu, id_min_time, ngood):

        spectra, indices, valid_bands = self.get_all_good_spectra_for_mu(index_mu)

        rrs_values = None
        if ngood == -1:
            ngood = spectra.shape[0]
        if spectra is not None and spectra.shape[0] == ngood:
            rrs_values = spectra[id_min_time, :]

        return rrs_values, indices, valid_bands

    def get_all_good_spectra_for_mu(self, index_mu):
        if index_mu < 0 or index_mu >= self.nmu:
            return None

        spectra = np.transpose(ma.array(self.insitu_rrs[index_mu]))

        vspectra = np.all(spectra == -999, axis=1)

        # ids_good = np.where(vspectra == False)[0]

        spectra_here = spectra[vspectra == False]

        if self.check_indices_by_mu:
            indices, valid_bands = self.get_insitu_indices_mu(index_mu)
        else:
            indices = self.wl_indices
            valid_bands = [True] * len(self.wl_indices)

        if indices is None:
            spectra_valid = None
        else:
            spectra_valid = spectra_here[:, indices]

        return spectra_valid, indices, valid_bands

    def get_all_good_spectra(self, check_indices):
        ntotal = self.nmu * self.nid

        nbands_here = len(self.wl_indices)

        valid_spectra = ma.zeros((ntotal, nbands_here))
        index_valid = 0
        for index_mu in range(self.nmu):
            spectra_here = self.get_good_spectra_for_mu(index_mu, check_indices)
            nspectra = spectra_here.shape[0]
            iini = index_valid
            ifin = index_valid + nspectra
            # valid_spectra[iini:ifin, :] = spectra_here[0:nspectra, self.iwl_min:self.iwl_max]
            valid_spectra[iini:ifin, :] = spectra_here[:, :]
            index_valid = ifin

        spectra_rec = valid_spectra[0:index_valid, :]

        return spectra_rec
