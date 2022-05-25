import numpy as np
import numpy.ma as ma
import BSC_QAA.bsc_qaa_EUMETSAT as bsc_qaa


class QC_INSITU:

    def __init__(self, insitu_rrs, insitu_bands):
        self.name = ''
        self.insitu_rrs = insitu_rrs
        self.insitu_bands = insitu_bands

        self.nmu = self.insitu_rrs.shape[0]
        self.nbands = self.insitu_rrs.shape[1]
        self.nid = self.insitu_rrs.shape[2]

        self.maxdifwl = 10

        self.time_max = 7200

        self.band_stats = None

        # self.wl_list = np.zeros((self.nmu, self.nbands))
        # self.wl_indices = np.zeros((self.nmu, self.nbands))
        # for imu in range(self.nmu):
        #     self.wl_list[imu, :] = self.insitu_bands[:]
        #     self.wl_indices[imu, :] = np.array(range(self.nbands))
        # self.mu_spectra_complete = None
        self.wl_list = list(self.insitu_bands[:])
        self.wl_indices = list(range(self.nbands))

        self.thersholds = None

        self.apply_band_shift = True

        self.check_indices_by_mu = True
        self.only_complete_spectra = False

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
        self.wl_indices = []
        for wl in self.wl_list:
            index, wl_index = self.get_insitu_index(wl)
            self.wl_indices.append(index)

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

    def check_validity_spectrum(self, rrs_values, index_mu):
        if rrs_values is None:
            return False
        if rrs_values.count() != len(rrs_values):
            return False
        if self.thersholds is None:
            return True
        check = True
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
        time_dif = self.time_max
        id_min_time = -1
        time_condition = False
        spectrum_complete = False
        ngood = 0


        for idx in range(len(dif_time_array)):
            t = dif_time_array[idx]
            if not ma.is_masked(t):
                ngood = ngood + 1
                if t < time_dif:
                    time_dif = t
                    id_min_time = idx
                    time_condition = True

        rrs_values = None
        valid_values = False
        if time_condition:

            rrs_values, indices, valid_bands = self.get_good_spectrum_for_mu(index_mu, id_min_time, ngood)

            # if wl_ref is not None:
            #     for wlr in wl_ref:
            #         index_wlr = self.set_wllist_using_wlref()
            # print(self.wl_list)
            # print(indices)
            # print(valid_bands)
            # wl_insitu_here = self.wl_list[indices]
            # print(wl_insitu_here)
            # print(wl_ref)

            if self.check_validity_spectrum(rrs_values,index_mu):

                # print(spectrum_complete, valid_bands)
                # print('-------------> ',len(rrs_values))
                valid_values = True
                if self.apply_band_shift and exact_wl_array is not None and wl_ref is not None:
                    if len(exact_wl_array.shape)==1:
                        exact_wl = exact_wl_array[indices]
                    else:
                        exact_wl = exact_wl_array[indices, id_min_time]
                    rrs_values = bsc_qaa.bsc_qaa(rrs_values, exact_wl, wl_ref)
                # print('*************> ', len(rrs_values))
                spectrum_complete = sum(valid_bands) == len(self.wl_list)
                if not spectrum_complete:
                    rrs_values[np.array(valid_bands)==False] = ma.masked

        return id_min_time, time_condition, valid_values, spectrum_complete, rrs_values

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
