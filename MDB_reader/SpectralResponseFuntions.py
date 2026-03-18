import os,__init__

import numpy as np
import pandas as pd
from netCDF4 import Dataset

class SRF:

    def __init__(self,path_srf,sensor):
        if path_srf is None:
            path_hval = os.path.dirname(os.path.dirname(__init__.__file__))
            if os.path.basename(path_hval)=='hypernets_val':
                path_hval = os.path.dirname(path_hval)
            path_srf = os.path.join(path_hval,'SpectralResponseFunctions')


        self.sensor = sensor
        self.file_srf = None
        if os.path.isdir(path_srf):
            file_srf  = None
            if sensor.upper() =='S3A' or sensor.upper() =='OLCIA':
                file_srf = os.path.join(path_srf,'S3A_OL_SRF_20160713_mean_rsr.nc4')
            elif sensor.upper()=='S3B' or sensor.upper() =='OLCIB':
                file_srf = os.path.join(path_srf,'S3B_OL_SRF_0_20180109_mean_rsr.nc4')

            if file_srf is not None and os.path.isfile(file_srf):
                self.file_srf = file_srf

        if self.file_srf is not None:
            print(f'[INFO] Started Spectral Response Function for sensor {self.sensor} with file: {self.file_srf}')

    def reindex_and_interpolate(self, df, new_index):

        new_indices = np.sort(np.unique(np.concat([df.index,new_index])))
        return df.reindex(new_indices).interpolate(method='index',limit_direction='both').loc[new_index]
        #return df.reindex(df.index | new_index).interpolate(method='index', limit_direction='both').loc[new_index]

    def get_olci_bands_info(self,wl_array,wl_sat_ref):
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
        olci_wl_all = np.array([float(olcibands[b]['wl']) for b in olcibands])


        olci_wl_list = []
        for b in olcibands:
            if wl_sat_ref is not None:
                wl = olcibands[b]['wl']
                # tal = olci_wl_all[int(np.argmin(np.abs(wl-olci_wl_all)))]
                # print('checking wl-->',wl,tal,)
                if np.abs(wl-wl_sat_ref[int(np.argmin(np.abs(wl-wl_sat_ref)))])>1.0:
                    continue

            wl_min = olcibands[b]['wl']-2
            wl_max = olcibands[b]['wl']+2
            indices = np.where((wl_array>=wl_min) & (wl_array<=wl_max))
            if len(indices[0])>=1:
                olcibands[b]['apply'] = True
                olci_wl_list.append(olcibands[b]['wl'])




        return olcibands,olci_wl_list


    def check_srf(self):
        if self.sensor is None:
            print(f'[ERROR] Sensor for spectral response function class is not defined')
            return None
        if self.file_srf is None:
            print(f'[ERROR] Spectral response function is not defined for sensor {self.sensor}')
            return None
        if self.sensor.upper().startswith('S3') or self.sensor.upper().startswith('OLCI'):
            return 'OLCI'

    ##apply convolutions. spectra_wl: nspectra x nwl, wl_array: nwl
    def apply_convolution_impl(self,spectra,wl_array,wl_sat_ref=None):
        check_srf = self.check_srf()
        if check_srf is None:
            return [None]*2
        if check_srf=='OLCI':
            return self.apply_convolution_olci_impl(spectra,wl_array,wl_sat_ref)

    def apply_convolution_olci_impl(self,spectra,wl_array,wl_sat_ref = None):
        df = pd.DataFrame(data=spectra.transpose(), index=wl_array)
        olci_bands_dict,wl_olci_apply = self.get_olci_bands_info(wl_array,wl_sat_ref)
        olci_bands_list = list(olci_bands_dict.keys())
        dataset = Dataset(self.file_srf)

        wl_values = np.array(dataset.variables['mean_spectral_response_function_wavelength'])
        msrf = np.array(dataset.variables['mean_spectral_response_function'])
        nwl = wl_values.shape[0]

        nwl_apply = len(wl_olci_apply)
        rrs_values = np.ma.masked_all((spectra.shape[0],nwl_apply))
        print('rrs values shape',rrs_values.shape)
        iwl_apply  = 0
        for iwl in range(nwl):
            band = olci_bands_list[iwl]
            if olci_bands_dict[band]['apply']:
                nnT = self.reindex_and_interpolate(df, wl_values[iwl, :])
                nnT_new = nnT.mul(np.array(msrf[iwl, :]), axis=0).sum(axis=0) / np.sum(np.array(msrf[iwl, :]))
                rrs_values[:,iwl_apply] = np.array(nnT_new)
                iwl_apply = iwl_apply + 1

        dataset.close()

        return rrs_values, wl_olci_apply

    #wl_array: nwl or 1 x nwl
    def apply_convolution_array(self,spectra,wl_array, dim_wl, wl_sat_ref = None):
        wl_array = np.squeeze(wl_array)
        nwl = wl_array.shape[0]

        wl_indices_valid = np.where(wl_array.mask==False)
        nwl_valid = len(wl_indices_valid[0])
        print(f'[INFO] Number of wavelengths: {nwl} Valid: {nwl_valid}')

        shape_orig = spectra.shape
        if len(shape_orig)==1 and len(spectra)==nwl:
            spectra = np.expand_dims(spectra,0)
            dim_wl = 1

        if dim_wl==-1:
            try:
                dim_wl = shape_orig.index(nwl)
            except:
                print(f'[ERROR] Convolution could not be started. Number of wavelengths {nwl} is not a dimension available in the spectra array: {shape_orig}')
                return
        else:
            if shape_orig[dim_wl]!=nwl:
                print(f'[ERROR] Dimension {dim_wl} has a length of {shape_orig[dim_wl]}, different from the number of detected wavelengths')
                return

        move_axis = True if dim_wl<len(shape_orig)-1 else False
        if move_axis:
            spectra = np.moveaxis(spectra, dim_wl, len(shape_orig) - 1)
            shape_proc = list(spectra.shape)
        else:
            shape_proc = list(shape_orig)
            
        reshape = True if len(shape_orig)>=3 else False
        if reshape:
            nspectra = np.prod(spectra.shape[:-1])
            spectra  = np.reshape(spectra,(nspectra,nwl))
        else:
            nspectra = spectra.shape[0]


        spectra_indices_valid = np.where(np.ma.count(spectra,axis=1)==nwl_valid)
        spectra_to_process = spectra[spectra_indices_valid[0],:]
        spectra_to_process = spectra_to_process[:,wl_indices_valid[0]]
        spectra_values,wl_output = self.apply_convolution_impl(spectra_to_process,wl_array[wl_indices_valid[0]],wl_sat_ref=wl_sat_ref)
        nwl_sat = spectra_values.shape[1]
        spectra_out = np.ma.masked_all((nspectra,nwl_sat))
        spectra_out[spectra_indices_valid[0]]= spectra_values
        if reshape:
            shape_proc[-1] = nwl_sat
            shape_proc = tuple(shape_proc)
            spectra_out = np.reshape(spectra_out,shape_proc)
        if move_axis:
            spectra_out = np.moveaxis(spectra_out,len(spectra_out.shape)-1,dim_wl)

        return spectra_out,wl_output
        





