import numpy as np
import numpy.ma as ma


class QC_INSITU:

    def __init__(self, insitu_rrs, insitu_bands):
        self.insitu_rrs = insitu_rrs
        self.insitu_bands = insitu_bands

        self.nmu = self.insitu_rrs.shape[0]
        self.nbands = self.insitu_rrs.shape[1]
        self.nid = self.insitu_rrs.shape[2]

        self.maxdifwl = 5

        self.band_stats = None

        self.wl_list = None
        self.thersholds = None

    def start_quality_control(self, wllist):
        if wllist is None:
            self.wl_list = self.insitu_bands
        else:
            self.wl_list = wllist

        self.thersholds = {}
        for wl in self.wl_list:
            self.thersholds[wl] = {
                'min_th': {
                    'type': 'None',
                    'value': 0
                },
                'max_th': {
                    'type': 'None',
                    'value': 0
                }
            }

    def compute_band_statistics(self,wlmin,wlmax):
        print(f'[INFO] Getting valid spectra...')
        spectra_res, self.wl_list = self.get_all_valid_spectra(wlmin,wlmax)
        band_stats = {}
        for index in range(len(self.wl_list)):
            wl = self.wl_list[index]
            wls = str(wl)
            print(f'[INFO] Computing statistics for band: {wls} ({index})')
            allvalid = spectra_res[:,index]
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
            #print(f'[INFO]Computing statistics...')
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



        # if self.wl_list is None:
        #     self.wl_list = self.insitu_bands
        #
        # for wl in self.wl_list:
        #     wls = str(wl)
        #     ntotal = self.nmu * self.nid
        #     alldata = np.arange(ntotal, dtype=float)
        #     iband_nearest, valvalid = self.get_nearest_insitu_band(0, wl)
        #     print(f'[INFO] Getting data for band: {wls} ({iband_nearest})')
        #     ivalid = 0
        #     for imu in range(self.nmu):
        #         # iband_nearest, valvalid = self.get_nearest_insitu_band(imu, wl)
        #         valvalid = self.get_values_insitu_band(imu,iband_nearest)
        #         if len(valvalid)>0:
        #             for val in valvalid:
        #                 alldata[ivalid] = val
        #                 ivalid = ivalid + 1
        #             # if allvalid is None:
        #             #     allvalid = valvalid
        #             # else:
        #             #     allvalid = np.concatenate([allvalid, valvalid])
        #     allvalid = alldata[0:ivalid]
        #     band_stats[wls] = {
        #         'nvalid': 0,
        #         'min_val': 0,
        #         'max_val': 0,
        #         'mean_val': 0,
        #         'std_val': 0,
        #         'median': 0,
        #         'p25': 0,
        #         'p75': 0
        #     }
        #     print(f'[INFO]Computing statistics...')
        #     if allvalid is not None:
        #         band_stats[wls]['nvalid'] = len(allvalid)
        #         band_stats[wls]['min_val'] = np.min(allvalid)
        #         band_stats[wls]['max_val'] = np.max(allvalid)
        #         band_stats[wls]['mean_val'] = np.mean(allvalid)
        #         band_stats[wls]['std_val'] = np.std(allvalid)
        #         band_stats[wls]['median'] = np.median(allvalid)
        #         band_stats[wls]['p25'] = np.percentile(allvalid, 25)
        #         band_stats[wls]['p75'] = np.percentile(allvalid, 75)
        #
        # self.band_stats = band_stats
        # return self.band_stats

    def get_nearest_insitu_band(self, imu, wl):
        dif_here = self.maxdifwl
        iband_nearest = -1
        valvalid = None
        for iband in range(self.nbands):
            dif = np.abs(wl - self.insitu_bands[iband])
            if dif <= dif_here:
                valarray = ma.array(self.insitu_rrs[imu, iband])
                valvalid = np.array(valarray[valarray.mask == False])
                if len(valvalid) > 0:
                    iband_nearest = iband
                    dif_here = dif
        return iband_nearest, valvalid

    def get_values_insitu_band(self, imu, iband):
        valarray = ma.array(self.insitu_rrs[imu, iband])
        valvalid = np.array(valarray[valarray.mask == False])
        return valvalid

    def get_spectra_for_mu(self, index_mu):
        if index_mu < 0 or index_mu >= self.nmu:
            return None

        spectra = np.transpose(np.array(self.insitu_rrs[index_mu]))
        vspectra = np.all(spectra == -999, axis=1)

        spectra_valid = spectra[vspectra == False]

        return spectra_valid

    # see compute_band_statistics for correct values of stat
    def get_stat_spectra(self, wlmin, wlmax, stat):
        if self.band_stats is None:
            self.compute_band_statistics()
        stat_spectra = []
        for wls in self.band_stats:
            wl = float(wls)
            if wlmin is not None and wl < wlmin:
                continue
            if wlmax is not None and wl > wlmax:
                continue
            stat_spectra.append(self.band_stats[wls][stat])

        return np.array(stat_spectra)

    def get_wl_min_max(self, wlmin, wlmax):
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

        return iwl_min, iwl_max

    def get_all_valid_spectra(self, wlmin, wlmax):
        ntotal = self.nmu * self.nid
        iwl_min, iwl_max = self.get_wl_min_max(wlmin, wlmax)

        nbands_here = iwl_max - iwl_min

        valid_spectra = np.zeros((ntotal, nbands_here))
        index_valid = 0
        for index_mu in range(self.nmu):
            spectra_here = self.get_spectra_for_mu(index_mu)
            nspectra = spectra_here.shape[0]
            iini = index_valid
            ifin = index_valid + nspectra
            valid_spectra[iini:ifin, :] = spectra_here[0:nspectra, iwl_min:iwl_max]
            index_valid = ifin

        spectra_rec = valid_spectra[0:index_valid, :]

        return spectra_rec, np.array(self.wl_list)
