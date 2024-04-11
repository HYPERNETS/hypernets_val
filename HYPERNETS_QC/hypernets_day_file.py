import os
from netCDF4 import Dataset
from datetime import datetime as dt
import numpy as np
#import __init__
#import sys

#code_home = os.path.dirname(os.path.dirname(__init__.__file__))
#sys.path.append(code_home)
from MDB_reader.PlotMultiple import PlotMultiple
from COMMON import Class_Flags_OLCI


class HYPERNETS_DAY_FILE():

    def __init__(self, file_nc, path_images):
        self.VALID = True
        if not os.path.exists(file_nc):
            self.VALID = False
        if file_nc is None:
            self.VALID = False
        self.file_nc = file_nc
        self.path_images = path_images
        if self.path_images is None:
            self.path_images = os.path.dirname(self.file_nc)
        self.path_images_date = None
        self.sequences = []
        self.format_img = '.png'
        if self.VALID:
            self.sequences = self.get_sequences()
            self.isequence = 0
            self.valid_sequences = None

        self.ref_wl = 800
        self.ref_wl_idx = -1
        self.rho = r'ρ$_w$'

        self.flag_builder = None

    def get_sequences(self):
        sequences = []
        dataset = Dataset(self.file_nc)
        seq_var = dataset.variables['sequence_ref']
        for idx in range(len(seq_var)):
            val = seq_var[idx]
            if np.ma.is_masked(val):
                sequences.append(None)
                continue
            time = dt.utcfromtimestamp(float(val))
            sequences.append(time.strftime('%Y%m%dT%H%M'))
        dataset.close()
        return sequences

    ##methods for multiple dates. start_date and end_date are dt objects, start_time and end_time with format %H:%M
    ##these parameters are retrieved from options_figures using check_dates_times
    def get_sequences_range(self, start_date, end_date, start_time, end_time):
        if start_date is None and end_date is None and start_time is None and end_time is None:
            return None, None
        if self.sequences is None:
            self.get_sequences()
        sequences_here = []
        sequences_indices = None
        for iseq in range(len(self.sequences)):
            seq = self.sequences[iseq]
            if seq is None:
                continue
            time_seq = dt.strptime(seq, '%Y%m%dT%H%M')
            if start_time is not None or end_time is not None:
                date_seq_str = time_seq.strftime('%Y-%m-%d')
            if start_date is not None:
                if time_seq < start_date:
                    continue
            if end_date is not None:
                end_date = end_date.replace(hour=23, minute=59, second=59)
                if time_seq > end_date:
                    continue
            if start_time is not None:
                start_time_ref = dt.strptime(f'{date_seq_str}T{start_time}', '%Y-%m-%dT%H:%M')
                if time_seq < start_time_ref:
                    continue
            if end_time is not None:
                end_time_ref = dt.strptime(f'{date_seq_str}T{end_time}', '%Y-%m-%dT%H:%M')
                if time_seq > end_time_ref:
                    continue
            sequences_here.append(seq)
            if sequences_indices is None:
                sequences_indices = [iseq]
            else:
                sequences_indices.append(iseq)

        return sequences_here, sequences_indices

    ##method only if file is for a specific data
    def get_sequences_interval(self, start_time, end_time):
        if self.sequences is None:
            self.get_sequences()
        sequences_here = []
        sequences_indices = None
        if len(self.sequences) == 0:
            return sequences_here
        for iseq in range(len(self.sequences)):
            seq = self.sequences[iseq]
            if seq is None:
                continue
            time_seq = dt.strptime(seq, '%Y%m%dT%H%M')
            if start_time <= time_seq <= end_time:
                sequences_here.append(seq)
                if sequences_indices is None:
                    sequences_indices = [iseq]
                else:
                    sequences_indices.append(iseq)

        return sequences_here, sequences_indices

    def get_report_files_interval(self, sequences_here, site, start_time, end_time):
        if sequences_here is None:
            sequences_here, range_here = self.get_sequences_interval(start_time, end_time)
        report_files = []
        if len(sequences_here) == 0:
            return report_files
        for seq in sequences_here:
            file_out = os.path.join(os.path.dirname(self.file_nc), f'{site}_{seq}_Report{self.format_img}')
            if os.path.exists(file_out):
                report_files.append(file_out)
        return report_files

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
        reflectance = np.array(dataset.variables[variable][self.isequence, self.ref_wl_idx, :])
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
        if angle_label == 270:
            angle_label = 225
        if angle_label == 90:
            angle_label = 135
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
        raa = dataset.variables['l1_rhof_raa'][self.isequence, 0]
        sza = dataset.variables['l1_rhof_sza'][self.isequence, 0]
        vza = dataset.variables['l1_rhof_vza'][self.isequence, 0]
        ws = dataset.variables['l1_rhof_wind'][self.isequence, 0]
        dataset.close()
        str = f'epsilon={epsilon:.4f};rho={rho:.4f}(raa={raa:.1f};sza={sza:.1f};vza={vza:.1f};ws={ws:.2f})'
        return str

    def get_title(self, site):
        date_time_here = dt.strptime(self.sequences[self.isequence], '%Y%m%dT%H%M')
        date_str = date_time_here.strftime('%Y-%m-%d')
        time_str = date_time_here.strftime('%H:%M')
        title = f'{site} {self.sequences[self.isequence]} - {date_str} {time_str} - {self.isequence + 1}/{len(self.sequences)}'
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
            if self.path_images_date is not None:
                file_img = os.path.join(self.path_images_date, name_file_img)
                if not os.path.exists(file_img):
                    file_img = None

        title = f'{flag}(az={var.oaa};zn={var.oza})'
        dataset.close()
        return file_img, title

    def set_path_images_date(self, site, date_here):
        folder_date = os.path.join(self.path_images, site, date_here.strftime('%Y'), date_here.strftime('%m'),
                                   date_here.strftime('%d'))
        if os.path.exists(folder_date):
            self.path_images_date = folder_date
        else:
            self.path_images_date = folder_date

    def get_water_images(self, site, date_here, time_min, time_max, interval_minutes):
        if time_min is None:
            time_min = '0600'
        if time_max is None:
            time_max = '1700'
        if interval_minutes is None:
            interval_minutes = 20
        self.set_path_images_date(site, date_here)
        from datetime import timedelta
        date_here_str = date_here.strftime('%Y%m%d')
        date_here_min = dt.strptime(f'{date_here_str}T{time_min}', '%Y%m%dT%H%M')
        date_here_max = dt.strptime(f'{date_here_str}T{time_max}', '%Y%m%dT%H%M')
        date_here = date_here_min
        water_images = {}
        while date_here <= date_here_max:
            date_here_str = date_here.strftime('%Y%m%dT%H%M')
            water_images[date_here_str] = {
                'file_img': None,
                'title': None
            }
            date_here = date_here + timedelta(minutes=interval_minutes)

        dataset = Dataset(self.file_nc)
        img_var = dataset.variables['pictures_water_rad']
        for idx in range(len(img_var)):

            val = img_var[idx]
            if np.ma.is_masked(val):
                continue
            self.isequence = idx
            file_img, title = self.get_img_file('water_rad')
            seq = self.sequences[self.isequence]
            print(idx, seq, file_img)
            if seq in water_images.keys():
                water_images[seq]['file_img'] = file_img
                paa = float(dataset.variables['l2_pointing_azimuth_angle'][idx])
                sza = float(dataset.variables['l2_solar_zenith_angle'][idx])
                saa = float(dataset.variables['l2_solar_azimuth_angle'][idx])
                title = f'{seq} (sza={sza:.1f};paa={paa:.1f};saa={saa:.1f})'
                title = f'{seq} (paa={paa:.1f})'
                water_images[seq]['title'] = title
        dataset.close()
        for tal in water_images:
            print(tal, '->', water_images[tal]['file_img'], '->', water_images[tal]['title'])
        return water_images

    def plot_water_images(self, wimages):
        file_out_base = os.path.join(os.path.dirname(self.file_nc), f'DayComparison_202304023_20230425')
        nimages = len(wimages)
        pm = PlotMultiple()
        nrow = nimages
        ncol = 2
        pm.start_multiple_plot_advanced(nrow, ncol, 5, 12, 0, 0.35, True)
        index_row = 0
        for wimage in wimages:
            pm.plot_image_title(wimages[wimage]['file_1'], index_row, 0, wimages[wimage]['title_1'])
            pm.plot_image_title(wimages[wimage]['file_2'], index_row, 1, wimages[wimage]['title_2'])
            index_row = index_row + 1
        pm.save_fig(f'{file_out_base}.{self.format_img}')
        pm.close_plot()

    def save_img_files(self, multiple_plot):
        flags = ['sky_irr_1', 'sky_rad_1', 'water_rad', 'sky_rad_2', 'sky_irr_2', 'sun']
        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)
        file_out_base = os.path.join(dir_img, f'CameraImages_{self.isequence}')
        if multiple_plot:
            pm = PlotMultiple()
            nrow = 2
            ncol = 3
            pm.start_multiple_plot_advanced(nrow, ncol, 10, 7.0, 0.1, 0.15, True)
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
                    # pm.plot_image_title(file_img,index_row,index_col,title)
                    pm.plot_image_hypernets(file_img, index_row, index_col, title)
                else:
                    pm.plot_blank_with_title(index_row, index_col, title[:9])

            index_col = index_col + 1
        if multiple_plot:
            pm.save_fig(f'{file_out_base}_all{self.format_img}')
            pm.close_plot()

    def save_spectra_files(self, multiple_plot):
        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)
        file_out_base = os.path.join(dir_img, f'Spectra_{self.isequence}')
        flags = ['irradiance', 'downwelling_radiance', 'upwelling_radiance', 'water_leaving_radiance',
                 'reflectance_nosc', 'reflectance']
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
        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)
        file_out_base = os.path.join(dir_img, f'Angles_{self.isequence}')
        flags = ['sza_nosc', 'saa_nosc', 'paa_nosc', 'sza', 'saa', 'paa']

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

    def save_report_image(self, site, delete_images, overwrite):
        print(f'[INFO] Sequence {self.isequence}: SEQ{self.sequences[self.isequence]}')
        if self.sequences[self.isequence] is None:
            return
        file_out = os.path.join(os.path.dirname(self.file_nc),
                                f'{site}_{self.sequences[self.isequence]}_Report{self.format_img}')
        if os.path.exists(file_out) and not overwrite:
            return
        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)
        file_angle = os.path.join(dir_img, f'Angles_{self.isequence}_all{self.format_img}')
        if not os.path.isfile(file_angle) or (os.path.isfile(file_angle) and overwrite):
            self.save_angle_files(True)
        file_img = os.path.join(dir_img, f'CameraImages_{self.isequence}_all{self.format_img}')
        if not os.path.isfile(file_img) or (os.path.isfile(file_img) and overwrite):
            self.save_img_files(True)
        file_spectra = os.path.join(dir_img, f'Spectra_{self.isequence}_all{self.format_img}')
        if not os.path.exists(file_spectra) or (os.path.isfile(file_spectra) and overwrite):
            self.save_spectra_files(True)
        pm = PlotMultiple()
        nrow = 3
        ncol = 1
        pm.start_multiple_plot_advanced(nrow, ncol, 10, 18, 0.1, 0.1, True)
        pm.get_axes(0, 0).set_title(self.get_title(site), fontsize=20)
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
        pm.get_axes(2, 0).set_title(info_l2, fontsize=12)
        pm.plot_image(file_spectra, 2, 0)

        pm.save_fig(file_out)
        pm.close_plot()

        if delete_images:
            for name in os.listdir(dir_img):
                file_here = os.path.join(dir_img, name)
                os.remove(file_here)
            os.remove(dir_img)

    def plot_angle(self, flag, ax_here):
        angle_flag = flag.split('_')[0]
        apply_nosc = False
        if flag.endswith('_nosc'):
            apply_nosc = True

        reflectance_ref = self.reflectance_ref(apply_nosc)
        reflectance_ref_l1 = self.reflectance_ref_l1(apply_nosc)
        angle_rad, angle_label = self.angle_rad(angle_flag)
        angle_rad_l1 = self.angle_rad_level1(angle_flag)

        reflectance_ref[reflectance_ref < 1e-5] = 1e-5
        reflectance_ref_l1[reflectance_ref_l1 < 1e-5] = 1e-5

        if self.valid_sequences is None:
            self.get_valid_flags()

        ax_here.set_rscale('log')
        ax_here.set_rlim((1e-5, 1))
        ax_here.set_rticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
        ax_here.set_theta_zero_location("N")
        ax_here.set_theta_direction(-1)
        ax_here.scatter(angle_rad[self.valid_sequences == 0], reflectance_ref[self.valid_sequences == 0], marker='o',
                        s=8, color='green')

        ax_here.scatter(angle_rad[self.valid_sequences > 0], reflectance_ref[self.valid_sequences > 0], marker='o', s=8,
                        color='red')

        ax_here.scatter(angle_rad_l1[:], reflectance_ref_l1[:], marker='s', s=8, color='gray')
        ax_here.scatter(angle_rad[self.isequence], reflectance_ref[self.isequence], marker='s', s=12, color='blue')

        val = (angle_rad[self.isequence] * 180) / np.pi
        if apply_nosc:
            title = f'{self.rho}({self.ref_wl})/NOSC vs.{angle_flag} ({val:.2f}°)'
        else:
            title = f'{self.rho}({self.ref_wl}) vs.{angle_flag} ({val:.2f}°)'
        ax_here.set_title(title, fontsize=10)

        if angle_flag == 'sza':
            ax_here.set_thetamin(0)
            ax_here.set_thetamax(90)
            ax_here.set_xticks((np.pi * np.array([0, 15, 30, 45, 60, 75, 90]) / 180))
            ax_here.set_rticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
        if angle_flag == 'saa' or angle_flag == 'paa':
            ax_here.set_rlabel_position(angle_label)
            # if angle_label==225:
            #     labels = ax_here.yticklabels()
            #     print(labels)

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
        spectra_l2 = None
        if l2_variable in dataset.variables:
            spectra_l2 = np.array(dataset.variables[l2_variable][self.isequence, :]).transpose()
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

    def start_dataset_w(self, file_out):
        dataset_w = Dataset(file_out, 'w', format='NETCDF4')
        dataset_r = Dataset(self.file_nc)
        # copy attributes
        dataset_w.setncatts(dataset_r.__dict__)

        # copy dimensions
        for name, dimension in dataset_r.dimensions.items():
            dataset_w.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

        # copy all file data except for the excluded
        for name, variable in dataset_r.variables.items():
            dataset_w.createVariable(name, variable.datatype, variable.dimensions)
            # copy variable attributes all at once via dictionary
            dataset_w[name].setncatts(dataset_r[name].__dict__)

            if name == 'wavelength' or name == 'bandwidth':
                dataset_w[name][:] = dataset_r[name][:]

        dataset_r.close()
        return dataset_w

    def add_flag_variables(self, dataset_w, flags):
        for flag in flags:
            var_name = flag
            if not var_name.startswith('flag_'):
                var_name = f'flag_{flag}'
            var = dataset_w.createVariable(var_name, 'i8', ('series',))
            var[:] = 0
            var.flag_values = flags[flag]['values']
            var.flag_meanings = ' '.join(flags[flag]['meanings'])
        return dataset_w

    def set_data_dataset_w(self, dataset_w, sindices, index_w):
        dini = index_w
        dfin = index_w + len(sindices)
        dataset_r = Dataset(self.file_nc)
        for variable in dataset_w.variables:
            if variable == 'wavelength' or variable == 'bandwidth':
                continue
            if variable not in dataset_r.variables:  ##for flag variables
                continue
            ndim = len(dataset_w[variable].shape)
            if ndim == 1:
                dataset_w[variable][dini:dfin] = dataset_r[variable][sindices]
            elif ndim == 2:
                dataset_w[variable][dini:dfin, :] = dataset_r[variable][sindices, :]
            elif ndim == 3:
                dataset_w[variable][dini:dfin, :, :] = dataset_r[variable][sindices, :, :]

        dataset_r.close()
        return dataset_w

    def set_data_flag(self, dataset_w, index_w, flag, flag_array):
        dini = index_w
        dfin = index_w + len(flag_array)
        var_name = flag
        if not var_name.startswith('flag_'):
            var_name = f'flag_{flag}'
        if var_name in dataset_w.variables:
            dataset_w[var_name][dini:dfin] = flag_array[:]
        return dataset_w

    def get_csv_col_names(self):
        col_names = ['sequence_ref']
        dataset_r = Dataset(self.file_nc)
        for var in dataset_r.variables:
            if var.startswith('l2'):
                ndim = len(dataset_r.variables[var].shape)
                if ndim == 1:
                    col_names.append(var[3:])

        col_names = col_names + ['rhow_nosc_800', 'rhow_800', 'rhow_nosc_350_450', 'rhow_350_450', 'isequence',
                                 'file_nc']

        dataset_r.close()

        return col_names

    def get_dataframe_lines(self, sindices, col_names):
        import pandas as pd
        if self.sequences is None:
            self.sequences = self.get_sequences()
        dataset_r = Dataset(self.file_nc)
        data = {}
        for c in col_names:
            if c == 'file_nc':
                data[c] = [self.file_nc] * len(sindices)
            elif c == 'isequence':
                data[c] = sindices
            elif c == 'sequence_ref':
                data[c] = [self.sequences[idx] for idx in sindices]
            elif c.startswith('rhow'):
                lc = c.split('_')
                if lc[1] == 'nosc':
                    wl_ini = float(lc[2])
                    variable = 'l2_reflectance_nosc'
                else:
                    wl_ini = float(lc[1])
                    variable = 'l2_reflectance'
                wl_fin = float(lc[-1])
                wls = np.array(dataset_r.variables['wavelength'][:])
                i_ini = np.argmin(np.abs(wl_ini - wls))
                i_fin = np.argmin(np.abs(wl_fin - wls)) + 1
                reflectance = dataset_r.variables[variable][sindices, i_ini:i_fin]
                data[c] = np.mean(reflectance, axis=1)
            else:
                variable = f'l2_{c}'
                data[c] = np.array(dataset_r.variables[variable][sindices])

        dataset_r.close()

        df = pd.DataFrame(data)
        return df

    def plot_from_options_impl(self, options_figure):
        if not options_figure['apply']:
            return
        if options_figure['type'] == 'spectraplot' and options_figure['type_rrs'] == 'user_defined':
            self.plot_spectra_plot_from_options(options_figure)

    def plot_angle_plot_from_options(self, options_figure):
        options_figure = self.check_gs_options_impl(options_figure, 'groupBy', 'groupType', 'groupValues')
        dataset = Dataset(self.file_nc)
        angle_variable = np.array(dataset.variables[options_figure['wl_variable']])

    def reduce_l1_dimensions(self, array):
        if len(array.shape) == 3:  ##spectral variables
            nspectra_total = array.shape[2] * array.shape[0]
            new_shape = (nspectra_total, array.shape[1])
            array_new = np.zeros(new_shape)
            for iscan in range(array.shape[2]):
                iini = iscan * array.shape[0]
                ifin = iini + array.shape[0]
                array_new[iini:ifin, :] = array[:, :, iscan]
            return array_new

    def multiply_by_scan(self, nseries, nscan, sequences_here, sequence_indices):
        nhere = len(sequences_here)
        ntotal = nhere * nscan
        sequence_here_total = [None] * ntotal
        sequence_indices_total = [None] * ntotal
        for iscan in range(nscan):
            iini = iscan * nhere
            ifin = iini + nhere
            sequence_here_total[iini:ifin] = sequences_here[:]
            sequence_indices_total[iini:ifin] = sequence_indices[:] + (iscan * nseries)

        return sequence_here_total, sequence_indices_total

    def multiply_array_by_scan(self, array, nseries, nscan):
        if array.shape[0] != nseries:
            return None
        ntotal = nseries * nscan
        array_new = np.zeros((ntotal,), dtype=array.dtype)
        for iscan in range(nscan):
            iini = iscan * nseries
            ifin = iini + nseries
            array_new[iini:ifin] = array[:]
        return array_new

    ##spectra without _fILLValue
    def get_spectra_stats(self, spectra_good):
        import statistics as st

        spectra_avg = np.mean(spectra_good, axis=0)
        spectra_std = np.std(spectra_good, axis=0)
        indices_max = np.argmax(spectra_good, axis=0)
        imax = st.mode(indices_max)
        spectra_max_real = spectra_good[imax, :]
        spectra_max = np.max(spectra_good, axis=0)

        indices_min = np.argmin(spectra_good, axis=0)
        imin = st.mode(indices_min)
        spectra_mim_real = spectra_good[imin, :]
        spectra_min = np.min(spectra_good, axis=0)

        spectra_median = np.median(spectra_good, axis=0)


        spectra_p25 = np.percentile(spectra_good, 25, axis=0)
        spectra_p75 = np.percentile(spectra_good, 75, axis=0)


        spectra_stats = {
            'avg': spectra_avg,
            'std': spectra_std,
            'spectra_min_real': spectra_mim_real,
            'spectra_max_real': spectra_max_real,
            'spectra_min': spectra_min,
            'spectra_max': spectra_max,
            'median': spectra_median,
            'p25': spectra_p25,
            'p75': spectra_p75
        }
        return spectra_stats


    def plot_spectra_plot_from_options(self, options_figure):
        plot_spectra = True
        if options_figure['plot_spectra'][0].lower() == 'none':
            plot_spectra = False
        plot_stats = options_figure['plot_stats']
        if not plot_spectra and not plot_stats:
            print('[WARNING] Please active plot_spectra or plot_stats')
            return

        ##start potential virtual flags
        options_figure = self.check_gs_options_impl(options_figure, 'groupBy', 'groupType', 'groupValues')

        ##DATA SELECTION
        dataset = Dataset(self.file_nc)
        wavelength = np.array(dataset.variables[options_figure['wl_variable']])
        y_variable = options_figure['y_variable']
        var_y = dataset.variables[y_variable]
        ndim = len(var_y.shape)
        if '_FillValue' in var_y.ncattrs():
            fill_value = var_y._FillValue
        else:
            from netCDF4 import default_fillvals
            fill_value = default_fillvals[var_y.dtype.str[1:]]
        nscan = 0
        nseries = 0
        if ndim == 3:  # l1 variable
            spectra = np.array(dataset.variables[y_variable][:, :, :])
            nseries = spectra.shape[0]
            nscan = spectra.shape[2]
            spectra = self.reduce_l1_dimensions(spectra)
        elif ndim == 2:
            spectra = np.array(dataset.variables[y_variable][:, :])
            nseries = spectra.shape[0]
        ##Checking spectra with fill values
        spectra_check = np.where(spectra == fill_value, 1, 0)
        spectra_check = np.sum(spectra_check, axis=1)

        sequences_here, sequence_indices = self.get_sequences_range(options_figure['start_date'],
                                                                    options_figure['end_date'],
                                                                    options_figure['start_time'],
                                                                    options_figure['end_time'])
        if ndim == 3 and sequence_indices is not None:
            sequences_here, sequence_indices = self.multiply_by_scan(nseries, nscan, sequences_here, sequence_indices)

        if sequence_indices is not None:
            spectra_check[sequence_indices] = spectra_check[sequence_indices] + 100
            spectra_check = np.where(spectra_check == 100, 0, 1)
        spectra = spectra[spectra_check == 0, :]
        dataset.close()
        ##DATA SELECTION

        ##CHECKING GROUP AND LEGEND
        ngroup = 1
        groupValues = None
        groupArray = None
        if 'groupValues' in options_figure.keys():
            groupValues = options_figure['groupValues']
        if groupValues is not None:
            ngroup = len(groupValues)
        str_legend = []
        handles = []
        if ngroup > 1 and options_figure['legend']:
            str_legend = self.get_str_legend(options_figure)
        if ngroup > 1:
            groupArray, all_flag_values, all_flag_meanings = self.get_gs_array(options_figure,
                                                                               options_figure['groupBy'],
                                                                               options_figure['groupType'])
            if ndim == 3:
                if len(groupArray.shape) == 1 and groupArray.shape[0] == nseries:
                    groupArray = self.multiply_array_by_scan(groupArray, nseries, nscan)
            groupArray = groupArray[spectra_check == 0]




        ##PLOTTING
        from MDB_reader.PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        pspectra.xdata = wavelength

        #single spectra plotting
        if plot_spectra:
            line_color = options_figure['line_color']
            marker = options_figure['marker']
            marker_size = options_figure['marker_size']
            line_type = options_figure['line_type']
            line_size = options_figure['line_size']
            if ngroup > 1:
                if len(line_color) != ngroup:  ##assing default colors
                    line_color = [line_color[0]] * ngroup
                    for idx in range(ngroup):
                        line_color[idx] = self.get_color_default(idx, 0, ngroup)
                if len(marker) != ngroup:
                    marker = [marker[0]] * ngroup
                if len(marker_size) != ngroup:
                    marker_size = [marker_size[0]] * ngroup
                if len(line_type) != ngroup:
                    line_type = [line_type[0]] * ngroup
                if len(line_size) != ngroup:
                    line_size = [line_size[0]] * ngroup

                for idx in range(ngroup):
                    val = groupValues[idx]
                    hline = pspectra.plot_single_line(spectra[groupArray == val, :].transpose(), line_color[idx],
                                                  line_type[idx], line_size[idx], marker[idx], marker_size[idx])
                    if len(hline) > 0:
                        handles.append(hline[0])
            else:
                pspectra.plot_single_line(spectra.tranpose(),line_color[0],line_type[0],line_size[0],marker[0],marker_size[0])

        if plot_stats:
            line_color = options_figure['line_color']
            if ngroup > 1:
                if len(line_color) != ngroup:  ##assing default colors
                    line_color = [line_color[0]] * ngroup
                    for idx in range(ngroup):
                        line_color[idx] = self.get_color_default(idx, 0, ngroup)
                for idx in range(ngroup):
                    val = groupValues[idx]
                    spectra_here = spectra[groupArray == val, :]
                    stats_here = self.get_spectra_stats(spectra_here)
                    pspectra.stats_style['fill']['color'] = line_color[idx]
                    pspectra.stats_style['central']['color'] = line_color[idx]
                    hline = pspectra.plot_stats(stats_here,None,None)
                    handles.append(hline[0])
            else:
                stats = self.get_spectra_stats(spectra)
                pspectra.plot_stats(stats, None, None)


        ##legend
        if len(str_legend) > 0:
            if len(handles) == 0:
                pspectra.set_legend(str_legend)
            else:

                pspectra.set_legend_h(handles, str_legend)

        ##y-range
        pspectra.set_y_range(options_figure['y_min'], options_figure['y_max'])

        #saveing to file
        if options_figure['file_out'] is not None:
            file_out = options_figure['file_out']
            pspectra.save_fig(file_out)

        pspectra.close_plot()

    def get_str_legend(self, options):
        if options['legend_values'] is not None:
            return options['legend_values']
        str_legend = []
        groupValues = options['groupValues']
        if groupValues is not None:
            ngroup = len(groupValues)
            if ngroup > 1:
                if options['groupType'] == 'float' or options['groupType'] == 'wavelength':
                    for g in groupValues:
                        str_legend.append(f'{g:.2f}')
                if options['groupType'] == 'flag':
                    flag_name = options['groupBy']
                    str_legend = self.get_flag_list(groupValues, options[flag_name]['flag_values'],
                                                    options[flag_name]['flag_meanings'])
                    if 'FUB' in str_legend:
                        index = str_legend.index('FUB')
                        if index >= 0:
                            str_legend[index] = 'S3 FUB-CSIRO'
                    if 'STANDARD' in str_legend:
                        index = str_legend.index('STANDARD')
                        if index >= 0:
                            str_legend[index] = 'WFR'
                    if 'POLYMER' in str_legend:
                        index = str_legend.index('POLYMER')
                        if index >= 0:
                            str_legend[index] = 'CMEMS-OLCI'
                    if 'CCIALL' in str_legend:
                        index = str_legend.index('CCIALL')
                        if index >= 0:
                            str_legend[index] = 'OC-CCI v.6 (complete time series)'
                    if 'CCI' in str_legend:
                        index = str_legend.index('CCI')
                        if index >= 0:
                            str_legend[index] = 'OC-CCI v.6 (OLCI period)'

        return str_legend

    def get_flag_list(self, values, allValues, allFlags):
        flag_list = []
        for val in values:
            if val == -1:
                flag_list.append('GLOBAL')
            indext = np.where(np.array(allValues) == val)
            index = indext[0]
            if len(index) == 1:
                indexf = index[0]
                flag_list.append(allFlags[indexf])
        return flag_list

    def get_gs_array(self, options_figure, by, type):
        dataset = Dataset(self.file_nc)
        if type == 'float':
            array_flag = np.array(dataset.variables[by])
            all_flag_values = None
            all_flag_meanings = None
        else:
            if by in dataset.variables:
                array_flag = np.array(dataset.variables[by][:])
                all_flag_values = dataset.variables[by].flag_values
                all_flag_meanings = dataset.variables[by].flag_meanings.split(' ')
            else:  ##previously built as a virtual flag in method check_gs_options_impl, line 904
                array_flag = options_figure[by]['flag_array']
                all_flag_values = options_figure[by]['flag_values']
                all_flag_meanings = options_figure[by]['flag_meanings']
        dataset.close()

        return array_flag, all_flag_values, all_flag_meanings

    def check_gs_options_impl(self, options_figure, by, type, values):
        dataset = Dataset(self.file_nc)
        var_group_name = options_figure[by]
        if options_figure[type] == 'flag':
            if var_group_name in dataset.variables:
                flag_values = dataset.variables[var_group_name].flag_values
                flag_meanings_list = dataset.variables[var_group_name].flag_meanings.split(' ')
                flag_meanings = [x.strip() for x in flag_meanings_list]
                options_figure[var_group_name] = {
                    'flag_values': flag_values,
                    'flag_meanings': flag_meanings
                }
            else:  ##virtual flag
                virtual_flags_options = self.flag_builder.get_virtual_flags_options()
                array, flag_meanings, flag_values = self.flag_builder.create_flag_array_ranges_v2(
                    virtual_flags_options[var_group_name])
                options_figure[var_group_name] = {
                    'flag_values': flag_values,
                    'flag_meanings': flag_meanings,
                    'flag_array': array
                }

            if options_figure[values] is None:
                options_figure[values] = flag_values
            else:
                flag_list_config = options_figure[values]
                flag_values_config = []
                for flag_config in flag_list_config:
                    if flag_config.strip() == 'GLOBAL':
                        flag_values_config.append(-1)
                        continue
                    try:
                        iflag = flag_meanings.index(flag_config.strip())
                        flag_values_config.append(flag_values[iflag])
                    except:
                        print(f'[WARNING] Flag {flag_config.strip()} is not in the list')
                        return None
                options_figure[values] = flag_values_config

        if options_figure[type] == 'float':
            if var_group_name not in dataset.variables:
                return None
            all_group_values = np.unique(np.array(dataset.variables[var_group_name]))
            if options_figure[values] is None:
                options_figure[values] = list(all_group_values)
            else:
                group_values_given = options_figure[values]
                group_values = []
                for val in group_values_given:
                    imin = np.argmin(np.abs(val - all_group_values))
                    if abs(val - all_group_values[imin]) < 0.1:
                        group_values.append((all_group_values[imin]))
                    else:
                        print(f'[WARNING] Value {val} is not in the variable {var_group_name}')
                        return None
                options_figure[values] = group_values
        dataset.close()
        return options_figure

    def get_color_default(self, value, min, max):
        nvalues = (max - min) + 1
        colors_default = ['Blue', 'Red', 'Green', 'm', 'Cyan', 'Orange', 'Yellow']
        if nvalues < 6:
            index = value - min
            print('index', index, '-->', colors_default[index])
            return colors_default[index]
        import matplotlib as mpl
        cm = mpl.colormaps['jet']
        return cm((value - min) / (max - min))

    def get_virtual_flag(self, options, var_name):
        flag_values = None
        flag_meanings = None
        flag_array = None
        # type_f = self.get_value_param(options, var_name, 'type', None, 'str')
        # type_v = self.get_value_param(options, var_name, 'typevirtual', 'flags_c', 'str')
        # if type_f != 'virtual_flag':
        #     return flag_values, flag_meanings, flag_array
        #
        # if type_v == 'spatial' or type_v == 'temporal' or type_v == 'ranges':
        #     #from OptionsManager import OptionsManager
        #     from MDB_reader.FlagBuilder import FlagBuilder
        #     fbuilder = FlagBuilder(self.path_nc, None, options)
        #     # options = fbuilder.get_options_dict(var_name)
        #     flag_values, flag_meanings, flag_array = fbuilder.create_flag_array(var_name, False)
        #
        #     return flag_values, flag_meanings, flag_array
