import os
from netCDF4 import Dataset
from datetime import datetime as dt
import numpy as np
import __init__
import sys

code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
from MDB_reader.PlotMultiple import PlotMultiple
from COMMON import Class_Flags_OLCI


class HYPERNETS_DAY_FILE():

    def __init__(self, file_nc):
        self.VALID = True
        if not os.path.exists(file_nc):
            self.VALID = False
        if file_nc is None:
            self.VALID = False
        self.file_nc = file_nc
        self.sequences = []
        self.format_img = '.png'
        if self.VALID:
            self.sequences = self.get_sequences()
            self.isequence = 0
            self.valid_sequences = None

        self.ref_wl = 800
        self.ref_wl_idx = -1
        self.rho = r'ρ$_w$'

    def get_sequences(self):
        sequences = []
        dataset = Dataset(self.file_nc)
        seq_var = dataset.variables['sequence_ref']
        for val in seq_var:
            time = dt.utcfromtimestamp(float(val))
            sequences.append(time.strftime('%Y%m%dT%H%M'))
        dataset.close()
        return sequences

    def get_ref_wl_idx(self):
        dataset = Dataset(self.file_nc)
        wavelength = np.array(dataset.variables['wavelength'][:])
        self.ref_wl_idx = np.argmin(np.abs(wavelength - self.ref_wl))
        dataset.close()

    def reflectance_ref(self, apply_nosc):
        if self.ref_wl_idx == -1:
            self.get_ref_wl_idx()
        variable = 'l2_reflectance'
        if apply_nosc:
            variable = 'l2_reflectance_nosc'
        dataset = Dataset(self.file_nc)
        reflectance = np.array(dataset.variables[variable][:, self.ref_wl_idx])
        dataset.close()
        return reflectance

    def reflectance_ref_l1(self, apply_nosc):
        if self.ref_wl_idx == -1:
            self.get_ref_wl_idx()
        variable = 'l1_reflectance'
        if apply_nosc:
            variable = 'l1_reflectance_nosc'
        dataset = Dataset(self.file_nc)
        reflectance = np.array(dataset.variables[variable][self.isequence, self.ref_wl_idx,:])
        dataset.close()
        return reflectance


    def angle_rad_level1(self, flag):
        if flag == 'sza':
            variable = 'l1_solar_zenith_angle'
        elif flag == 'saa':
            variable = 'l1_solar_azimuth_angle'
        elif flag == 'paa':
            variable = 'l1_pointing_azimuth_angle'
        dataset = Dataset(self.file_nc)
        angle = np.array(dataset.variables[variable][self.isequence, :])
        dataset.close()
        angle = (angle * np.pi) / 180
        return angle

    def angle_rad(self, flag):
        if flag == 'sza':
            variable = 'l2_solar_zenith_angle'
        elif flag == 'saa':
            variable = 'l2_solar_azimuth_angle'
        elif flag == 'paa':
            variable = 'l2_pointing_azimuth_angle'
        dataset = Dataset(self.file_nc)
        angle = np.array(dataset.variables[variable][:])
        dataset.close()
        median_op_angle = np.median(angle) + 180
        if median_op_angle > 360:
            median_op_angle = median_op_angle - 360
        angle_labels = np.arange(0, 361, 45)
        angle_label = angle_labels[int(np.argmin(np.abs(angle_labels - median_op_angle)))]
        if angle_label == 360 or angle_label == 0:
            angle_label = 315
        angle = (angle * np.pi) / 180
        return angle, angle_label

    def get_valid_flags(self):
        dataset = Dataset(self.file_nc)
        self.valid_sequences = np.array(dataset.variables['l2_quality_flag'][:])
        dataset.close()

    def get_flags_sequence(self):
        dataset = Dataset(self.file_nc)
        flag_value = np.uint64(dataset.variables['l2_quality_flag'][self.isequence])
        all_flag_values = [np.uint64(x) for x in dataset.variables['l2_quality_flag'].flag_masks.split(',')]
        all_flag_meaninigs = dataset.variables['l2_quality_flag'].flag_meanings
        cflags = Class_Flags_OLCI.Class_Flags_OLCI(all_flag_values, all_flag_meaninigs)
        list, mask = cflags.Decode(flag_value)
        dataset.close()
        return list

    def get_info_l2(self):
        dataset = Dataset(self.file_nc)
        epsilon = dataset.variables['l2_epsilon'][self.isequence]
        rho = dataset.variables['l2_rhof'][self.isequence]
        raa = dataset.variables['l1_rhof_raa'][self.isequence,0]
        sza = dataset.variables['l1_rhof_sza'][self.isequence,0]
        vza = dataset.variables['l1_rhof_vza'][self.isequence, 0]
        ws = dataset.variables['l1_rhof_wind'][self.isequence,0]
        dataset.close()
        str = f'epsilon={epsilon:.4f};rho={rho:.4f}(raa={raa:.1f};sza={sza:.1f};vza={vza:.1f};ws={ws:.2f})'
        return str

    def get_title(self):
        date_time_here = dt.strptime(self.sequences[self.isequence], '%Y%m%dT%H%M')
        date_str = date_time_here.strftime('%Y-%m-%d')
        time_str = date_time_here.strftime('%H:%M')
        title = f'SEQ{self.sequences[self.isequence]} - {date_str} {time_str} - {self.isequence + 1}/{len(self.sequences)}'
        return title

    # flag: name_variable, sky_irr_1, sky_irr_2, sky_rad_1, sky_rad_1, water_rad, sun
    def get_img_file(self, flag):
        if flag.startswith('pictures_'):
            name_var = flag
        else:
            name_var = f'pictures_{flag}'
        dataset = Dataset(self.file_nc)
        if name_var not in dataset.variables:
            return None

        var = dataset.variables[name_var]
        prefix = var.prefix
        suffix = var.suffix
        seq_here = self.sequences[self.isequence]

        val = var[self.isequence]
        file_img = None
        if not np.ma.is_masked(val):
            time = dt.utcfromtimestamp(float(val))
            time_str = time.strftime('%Y%m%dT%H%M')
            name_file_img = f'{prefix}_{seq_here}_{time_str}_{suffix}'
            file_img = os.path.join(os.path.dirname(self.file_nc), name_file_img)
            if not os.path.exists(file_img):
                file_img = None

        title = f'{flag}(az={var.oaa};zn={var.oza})'
        dataset.close()
        return file_img, title

    def save_img_files(self, multiple_plot):
        flags = ['sky_irr_1', 'sky_rad_1', 'water_rad', 'sky_rad_2', 'sky_irr_2', 'sun']
        file_out_base = os.path.join(os.path.dirname(self.file_nc), f'CameraImages_{self.isequence}')
        if multiple_plot:
            pm = PlotMultiple()
            nrow = 2
            ncol = 3
            pm.start_multiple_plot_advanced(nrow, ncol, 10, 5.2, 0, 0.15, True)
        index_row = 0
        index_col = 0
        for flag in flags:
            file_img, title = self.get_img_file(flag)
            if index_col == ncol:
                index_col = 0
                index_row = index_row + 1
            # print(f'{flag}->{index_row} {index_col}')
            if multiple_plot:
                if file_img is not None:
                    pm.plot_image_title(file_img, index_row, index_col, title)
                else:
                    pm.plot_blank_with_title(index_row, index_col, title)

            index_col = index_col + 1
        if multiple_plot:
            pm.save_fig(f'{file_out_base}_all{self.format_img}')
            pm.close_plot()

    def save_spectra_files(self, multiple_plot):
        file_out_base = os.path.join(os.path.dirname(self.file_nc), f'Spectra_{self.isequence}')
        flags = ['irradiance', 'downwelling_radiance', 'upwelling_radiance', 'water_leaving_radiance', 'reflectance_nosc','reflectance']
        if multiple_plot:
            pm = PlotMultiple()
            nrow = 2
            ncol = 3
            pm.start_multiple_plot_advanced(nrow, ncol, 10, 5.5, 0.25, 0.30, True)
        index_row = 0
        index_col = 0
        for flag in flags:
            if multiple_plot:
                if index_col == ncol:
                    index_col = 0
                    index_row = index_row + 1
                ax_here = pm.get_axes(index_row, index_col)
                self.plot_spectra(flag, ax_here)
                index_col = index_col + 1

        if multiple_plot:
            pm.save_fig(f'{file_out_base}_all{self.format_img}')
            pm.close_plot()

    def save_angle_files(self, multiple_plot):
        file_out_base = os.path.join(os.path.dirname(self.file_nc), f'Angles_{self.isequence}')
        flags = ['sza_nosc', 'saa_nosc', 'paa_nosc','sza', 'saa', 'paa']
        if multiple_plot:
            pm = PlotMultiple()
            nrow = 2
            ncol = 3
            pm.start_multiple_plot_polar(nrow, ncol, 10, 6, 0.25, 0.40, True)
        index_row = 0
        index_col = 0
        for flag in flags:
            if multiple_plot:
                if index_col == ncol:
                    index_col = 0
                    index_row = index_row + 1
                ax_here = pm.get_axes(index_row, index_col)
                self.plot_angle(flag, ax_here)
                index_col = index_col + 1
        if multiple_plot:
            pm.save_fig(f'{file_out_base}_all{self.format_img}')
            pm.close_plot()

    def save_report_image(self, delete_images, overwrite):
        print(f'[INFO] Sequence {self.isequence}: SEQ{self.sequences[self.isequence]}')
        file_out = os.path.join(os.path.dirname(self.file_nc), f'ReportImage_{self.isequence}{self.format_img}')
        if os.path.exists(file_out) and not overwrite:
            return
        file_angle = os.path.join(os.path.dirname(self.file_nc), f'Angles_{self.isequence}_all{self.format_img}')
        if not os.path.isfile(file_angle):
            self.save_angle_files(True)
        file_img = os.path.join(os.path.dirname(self.file_nc), f'CameraImages_{self.isequence}_all{self.format_img}')
        if not os.path.isfile(file_img):
            self.save_img_files(True)
        file_spectra = os.path.join(os.path.dirname(self.file_nc), f'Spectra_{self.isequence}_all{self.format_img}')
        if not os.path.exists(file_spectra):
            self.save_spectra_files(True)
        pm = PlotMultiple()
        nrow = 3
        ncol = 1
        pm.start_multiple_plot_advanced(nrow, ncol, 10, 18, 0.1, 0.1, True)
        pm.get_axes(0, 0).set_title(self.get_title(), fontsize=20)
        pm.plot_image(file_img, 0, 0)
        flag_list = self.get_flags_sequence()
        if len(flag_list) == 0:
            str_flag_list = 'NO FLAGGED'
        else:
            str_list = ','.join(flag_list)
            str_flag_list = f'FLAGS:{str_list}'
        pm.get_axes(1, 0).set_title(str_flag_list, fontsize=12)
        pm.plot_image(file_angle, 1, 0)
        info_l2 = self.get_info_l2()
        pm.get_axes(2,0).set_title(info_l2,fontsize=12)
        pm.plot_image(file_spectra, 2, 0)

        pm.save_fig(file_out)
        pm.close_plot()

    def plot_angle(self, flag, ax_here):
        angle_flag = flag.split('_')[0]
        apply_nosc = False
        if flag.endswith('_nosc'):
            apply_nosc = True

        reflectance_ref = self.reflectance_ref(apply_nosc)
        reflectance_ref_l1 = self.reflectance_ref_l1(apply_nosc)
        angle_rad, angle_label = self.angle_rad(angle_flag)
        angle_rad_l1 = self.angle_rad_level1(angle_flag)

        if self.valid_sequences is None:
            self.get_valid_flags()

        ax_here.scatter(angle_rad_l1, reflectance_ref_l1, marker='s', color='gray', s=8)

        ax_here.scatter(angle_rad[self.valid_sequences == 0], reflectance_ref[self.valid_sequences == 0], marker='o',
                        s=8, color='green')
        ax_here.scatter(angle_rad[self.valid_sequences > 0], reflectance_ref[self.valid_sequences > 0], marker='o', s=8,
                        color='red')

        ax_here.scatter(angle_rad[self.isequence], reflectance_ref[self.isequence], marker='s', color='blue')

        val = (angle_rad[self.isequence] * 180) / np.pi
        if apply_nosc:
            title = f'{self.rho}({self.ref_wl})/NOSC vs.{angle_flag} ({val:.2f}°)'
        else:
            title = f'{self.rho}({self.ref_wl}) vs.{angle_flag} ({val:.2f}°)'
        ax_here.set_title(title, fontsize=10)

        if angle_flag == 'sza':
            ax_here.set_thetamin(0)
            ax_here.set_thetamax(90)
        if angle_flag == 'saa':
            ax_here.set_rlabel_position(angle_label)
        if angle_flag == 'paa':
            ax_here.set_rlabel_position(angle_label)
        ax_here.tick_params(axis='x', labelsize=10)
        ax_here.tick_params(axis='y', labelsize=10)

    def plot_spectra(self, flag, ax_here):
        info = {
            'irradiance': {
                'ylabel': 'Ed',
                'title': 'Downwelling Irradiance'
            },
            'downwelling_radiance': {
                'ylabel': 'Li',
                'title': 'Downwelling radiance'
            },
            'upwelling_radiance': {
                'ylabel': 'Lt',
                'title': 'Upwelling radiance'
            },
            'water_leaving_radiance': {
                'ylabel': 'Lw',
                'title': 'Water-leaving radiance'
            },
            'reflectance': {
                'ylabel': r'ρ$_w$',
                'title': 'Reflectance'
            },
            'reflectance_nosc': {
                'ylabel': r'ρ$_w$',
                'title': 'Reflectance Nosc'
            },
        }
        l1_variable = f'l1_{flag}'
        l2_variable = f'l2_{flag}'
        dataset = Dataset(self.file_nc)
        # print(flag)
        spectra = np.array(dataset.variables[l1_variable][self.isequence, :, :]).transpose()
        spectra_l2 =None
        if l2_variable in dataset.variables:
            spectra_l2 = np.array(dataset.variables[l2_variable][self.isequence,:]).transpose()
        wavelength = np.array(dataset.variables['wavelength'])
        for ispectra in range(spectra.shape[0]):
            ax_here.plot(wavelength, spectra[ispectra, :], color='gray', linewidth=0.5)
        if spectra_l2 is not None:
            ax_here.plot(wavelength, spectra_l2, color='black', linewidth=0.25)



        ax_here.set_xlabel('Wavelength(nm)', fontsize=7)
        ax_here.set_ylabel(info[flag]['ylabel'], fontsize=7)
        ax_here.tick_params(axis='x', labelsize=7)
        ax_here.tick_params(axis='y', labelsize=7)
        ax_here.set_title(info[flag]['title'], fontsize=7)
        dataset.close()

    def check_img_files(self):
        flags = ['sky_irr_1', 'sky_irr_2', 'sky_rad_1', 'sky_rad_1', 'water_rad', 'sun']
        for flag in flags:
            file_img, title = self.get_img_file(flag)
            print(f'[INFO] Sequence: {self.sequences[self.isequence]}')
            if file_img is None:
                print(f'[INFO] {flag}-> N/Av')
            else:
                print(f'[INFO] {flag}->{file_img}->{os.path.exists(file_img)}')
