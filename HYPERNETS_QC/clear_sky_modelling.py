import os.path
import numpy as np
import __init__
from MDB_reader.PlotSpectra import PlotSpectra

class ClearSkyModel():
    def __init__(self, path_model):
        self.path_model = path_model
        if path_model is None:
            self.path_model = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__init__.__file__))),
                                           'LRTModel')

        self.prename = 'uvspec_clear_sza'
        self.output_suffix = ['00', '10', '20', '30', '40', '50', '60', '70', '80', '82', '84', '86', '88']
        self.min_values = np.array([float(x) for x in self.output_suffix])
        self.max_values = self.min_values + 10
        self.max_values[8:] = self.min_values[8:] + 2


    def check_model_validity(self):
        valid = True
        for suffix in self.output_suffix:
            file_model = os.path.join(self.path_model,f'{self.prename}{suffix}.out')
            if not os.path.isfile(file_model):
                print(f'[ERROR] Error loading clear sky model. File: {file_model} is not available')
                valid = False
        return valid

    def get_lrt_model_geometry(self, sza, oaa):
        index_sza = np.where(np.logical_and(sza >= self.min_values, sza < self.max_values))
        if len(index_sza[0]) == 0:
            print(f'[ERROR] Solar zenith angle should be in the range 0 - 90°')
            return None
        oaa_array = np.arange(0, 315 + 1, 45)  # [0 45 90 135 180 225 270 315]
        index_oaa = np.where(oaa_array == oaa)
        if len(index_oaa[0]) == 0:
            print(f'[ERROR] Observation azimuth angle should have on the following values: {oaa_array}')
            return None
        index_file_ini = int(index_sza[0])
        file_lrt_ini = os.path.join(self.path_model, f'{self.prename}{self.output_suffix[index_file_ini]}.out')
        LRT_wave, LRT_eddir_ini, LRT_eddif_ini, LRT_edtot_ini, LRT_ld40_ini = self.get_lrt_model(file_lrt_ini)
        nwl = LRT_wave.shape[0]
        index_file_end = index_file_ini + 1
        last_point = False
        if index_file_end == len(self.min_values):
            index_file_end = len(self.min_values) - 1
            file_lrt_end = file_lrt_ini
            LRT_eddir_end = LRT_eddir_ini
            LRT_eddif_end = LRT_eddif_ini
            LRT_edtot_end  = LRT_edtot_ini
            LRT_ld40_end = LRT_ld40_ini
            index_file_ini = index_file_end-1
            file_lrt_ini = os.path.join(self.path_model, f'{self.prename}{self.output_suffix[index_file_ini]}.out')
            LRT_wave, LRT_eddir_ini, LRT_eddif_ini, LRT_edtot_ini, LRT_ld40_ini = self.get_lrt_model(file_lrt_ini)
            last_point = True

        if not last_point:
            file_lrt_end = os.path.join(self.path_model, f'{self.prename}{self.output_suffix[index_file_end]}.out')
            LRT_wave, LRT_eddir_end, LRT_eddif_end, LRT_edtot_end, LRT_ld40_end = self.get_lrt_model(file_lrt_end)

        x_interpolate = np.cos(np.deg2rad(sza))
        x_known = np.array(
            [np.cos(np.deg2rad(self.min_values[index_file_end])), np.cos(np.deg2rad(self.min_values[index_file_ini]))])
        x_known = x_known.repeat(nwl).reshape((2, nwl))



        ##Ld
        ld_ini = np.squeeze(LRT_ld40_ini[:, index_oaa])
        ld_end = np.squeeze(LRT_ld40_end[:, index_oaa])
        y_known = np.array([ld_end[:], ld_ini[:]])
        ld_array = self.multiple_interpolate(x_interpolate, x_known, y_known)

        ##Ed dir
        y_known = np.array([LRT_eddir_end, LRT_eddir_ini])
        eddir_array = self.multiple_interpolate(x_interpolate, x_known, y_known)

        ##Ed dif
        y_known = np.array([LRT_eddif_end, LRT_eddif_ini])
        eddif_array = self.multiple_interpolate(x_interpolate, x_known, y_known)

        ##Ed tot
        edtot_array = eddif_array + eddir_array

        ##TEST: interpolation using a loop based on np.interp to check that multiple_interpolate produce the same result
        # y_interpolate = self.multiple_interpolate(x_interpolate, x_known, y_known)
        # for idx in range(nwl):
        #     x_known_test = [np.cos(np.deg2rad(self.min_values[index_file_end])),np.cos(np.deg2rad(self.min_values[index_file_ini]))]
        #     y_known_test = [ld_end[idx],ld_ini[idx]]
        #     y_interpolate_test = np.interp(x_interpolate,x_known_test,y_known_test)
        #     print(f'{y_interpolate[idx]} / {y_interpolate_test} = {y_interpolate[idx]/y_interpolate_test}')

        return LRT_wave, eddir_array, eddif_array, edtot_array, ld_array

    def multiple_interpolate(self, x_interpolate, x_know, y_know):
        return y_know[0, :] + (x_interpolate - x_know[0, :]) * (
                    (y_know[1, :] - y_know[0, :]) / (x_know[1, :] - x_know[0, :]))

    def get_lrt_model(self, LRTfile):
        ifind = LRTfile.find("sza")
        # sza = LRTfile[ifind + 3:ifind + 5]  # finds sza (as float??)
        encoding = 'utf-8'
        LRT_nwave = 1050 - 350 + 1  ##701
        LRT_nphi = 8  # 0° to 315° step 45°
        # LRT_phi = np.arange(0,315+1,45) ##sensor relative azimuth angle [ 0 45  90 135 180 225 270 315]

        LRT_wave = np.zeros(LRT_nwave)  ##wavelength
        LRT_eddir = np.zeros(LRT_nwave)  ##Ed direct
        LRT_eddif = np.zeros(LRT_nwave)  ## Ed  diffuse
        LRT_edtot = np.zeros(LRT_nwave)  ## Ed tot
        LRT_ld40 = np.zeros((LRT_nwave, LRT_nphi))  ##Ld(Li) with oza=40 and oaa=[ 0  45  90 135 180 225 270 315]

        ##linear interpolation, better using cosine(sza). Use oaa = 90 o 135 (to be checked)

        with open(LRTfile, 'r', encoding=encoding) as f:
            for i, line in enumerate(f.readlines()):
                line = line.strip()
                if len(line) == 0: continue
                iwave, jline = divmod(i, 3)  # alternatively a%b gives remainder of a/b
                if jline == 0:
                    # wavelength and ed output
                    LRT_wave[iwave] = line.split()[0]
                    LRT_eddir[iwave] = line.split()[1]
                    LRT_eddif[iwave] = line.split()[2]
                    LRT_edtot[iwave] = LRT_eddir[iwave] + LRT_eddif[iwave]
                elif jline == 1:
                    pass
                else:
                    # radiances
                    LRT_ld40[iwave, :] = line.split()[2:]

        return LRT_wave, LRT_eddir, LRT_eddif, LRT_edtot, LRT_ld40

    ##get lrt model results interpolated for a new wavelenght array
    def get_lrt_model_geometry_wl(self,sza, ooa, wl_array):
        LRT_wave, eddir_array, eddif_array, edtot_array, ld_array = self.get_lrt_model_geometry(sza,ooa)

        ld_array_interp = np.interp(wl_array,LRT_wave,ld_array)
        eddir_array_interp = np.interp(wl_array, LRT_wave, eddir_array)
        eddif_array_interp = np.interp(wl_array, LRT_wave, eddif_array)
        edtot_array_interp = eddif_array_interp + eddir_array_interp

        return ld_array_interp,eddir_array_interp,eddif_array_interp,edtot_array_interp

def make_test():
    # Ed clear sky modelling
    # KGR 2024-10-30
    # run from base not acolite conda environment ?
    # set here directories for local machine
    import sys
    import os, glob
    ifpng = True
    ifpdf = True
    # clear sky model data from HYPERNETS processor github
    CSMpath = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/CLEAR_SKY_MODEL'
    # specify output directory for plots
    outpath = os.path.join(CSMpath, 'PLOTS')
    if not os.path.isdir(outpath):
        os.mkdir(outpath)

    # using output from enhydra libradtran simulations
    LRTpath = CSMpath
    # specify directory for plots
    figpath = LRTpath + "/figures"

    LRTfiles = glob.glob('{}/uvspec*.out'.format(LRTpath))
    LRTfiles.sort()
    print('Found {} LRT files'.format(len(LRTfiles)))
    print(LRTfiles)

    # read LRT output files
    import numpy as np

    LRT_nsza = 13  # 0° to 80° step 10° then 82 to 88° step 2°
    LRT_nmu = 1  # only 40° nadir for radiance
    LRT_nphi = 8  # 0° to 315° step 45°
    LRT_nwave = 1050 - 350 + 1
    print("nwaves=", LRT_nwave)

    LRT_ld40 = np.zeros((LRT_nwave, LRT_nsza, LRT_nphi, LRT_nmu))  ##Ld or Li
    LRT_eddir = np.zeros((LRT_nwave, LRT_nsza))  ##Ed dir
    LRT_eddif = np.zeros((LRT_nwave, LRT_nsza))  ##Ed dif
    LRT_edtot = np.zeros((LRT_nwave, LRT_nsza))  ##Ed tot

    LRT_sza = np.zeros(LRT_nsza)
    LRT_wave = np.zeros(LRT_nwave)
    LRT_phi = np.zeros(LRT_nphi)
    LRT_mu = np.zeros(LRT_nmu)

    encoding = 'utf-8'

    nLRTfiles = len(LRTfiles)
    print(nLRTfiles)

    for isza, LRTfile in enumerate(LRTfiles):

        with open(LRTfile, 'r', encoding=encoding) as f:
            print(LRTfile)
            ifind = LRTfile.find("sza")
            LRT_sza[isza] = LRTfile[ifind + 3:ifind + 5]  # finds sza (as float??)
            print("file sza=", LRT_sza[isza])
            for i, line in enumerate(f.readlines()):
                line = line.strip()
                if len(line) == 0: continue

                # print(i,line)
                iwave, jline = divmod(i, 3)  # alternatively a%b gives remainder of a/b
                # print(iwave,jline)

                if jline == 0:
                    # wavelength and ed output
                    LRT_wave[iwave] = line.split()[0]
                    LRT_eddir[iwave, isza] = line.split()[1]
                    LRT_eddif[iwave, isza] = line.split()[2]
                    LRT_edtot[iwave, isza] = LRT_eddir[iwave, isza] + LRT_eddif[iwave, isza]
                    # print(LRT_wave[iwave])
                elif jline == 1:
                    # phi colum headers
                    # print("col headers",line)
                    dummy = 0
                else:
                    # radiances
                    LRT_ld40[iwave, isza, :, 0] = line.split()[2:]
                    # print(line.split()[2:])
                    if iwave == 50:  # check outputs for 400 nm
                        print(LRT_wave[iwave])
                        print(LRT_ld40[iwave, isza, :, 0])

    # now plot the LRT sim files that have been read

    # wcol=tab["blue","green","red"]
    from matplotlib import pyplot as plt
    cm = plt.cm.get_cmap('tab20')
    wcol = [cm.colors[i] for i in range(LRT_nsza)]
    print(wcol)

    myfig = plt.figure(figsize=(10, 30), facecolor='white')  # set the size of the figures for spectral plots
    myfig, ((ax1, ax4), (ax2, ax3)) = plt.subplots(2, 2)  # nrows, ncols
    myfig.tight_layout()

    aziout = 2  # element 2 of index3 is 90° azimuth; element 3 of index3 is 135° azimuth; element 0 of index3 is 0° azimuth (towards sun)

    # loop over sims (sza)
    # for isza in range(nLRTfiles):
    for isza in range(9):  # limit to 80°
        # for isza in [6]:
        ax1.plot(LRT_wave[:], LRT_edtot[:, isza], label="sza=" + str(LRT_sza[isza]) + "°",
                 color=wcol[isza])  # ,color=wcol[iwave],label=wplots[iwave])
        ax2.plot(LRT_wave[:], LRT_ld40[:, isza, aziout, 0], label="sza=" + str(LRT_sza[isza]) + "°", linestyle="solid",
                 color=wcol[isza])  # ,color=wcol[iwave],label=wplots[iwave]) element 2 of index3 is 90° azimuth
        ax3.plot(LRT_wave[:], LRT_ld40[:, isza, aziout, 0] / LRT_edtot[:, isza],
                 label="sza=" + str(LRT_sza[isza]) + "°", linestyle="solid",
                 color=wcol[isza])  # ,color=wcol[iwave],label=wplots[iwave])
        ax4.plot(LRT_wave[:], LRT_eddir[:, isza] / LRT_edtot[:, isza], label="sza=" + str(LRT_sza[isza]) + "°",
                 linestyle="solid", color=wcol[isza])  # ,color=wcol[iwave],label=wplots[iwave])

    ax1.plot([400, 1050], [0, 0], linestyle="dotted", color="grey")  # ,color=wcol[iwave],label=wplots[iwave])
    ax2.plot([400, 1050], [0, 0], linestyle="dotted", color="grey")  # ,color=wcol[iwave],label=wplots[iwave])
    ax3.plot([400, 1050], [0.05, 0.05], linestyle="dashed", color="grey")  # ,color=wcol[iwave],label=wplots[iwave])
    ax3.plot([400, 1050], [0, 0], linestyle="dotted", color="grey")  # ,color=wcol[iwave],label=wplots[iwave])
    ax4.plot([400, 1050], [0, 0], linestyle="dotted", color="grey")  # ,color=wcol[iwave],label=wplots[iwave])

    wavemin_plt = 400
    wavemax_plt = 1050

    ax1.set_xlim([wavemin_plt, wavemax_plt])
    ax1.set_ylabel("Edtot (mW m-2 nm-1)")
    ax1.set_xlabel("Wavelength (nm)")

    ax2.set_xlim([wavemin_plt, wavemax_plt])
    ax2.set_ylabel("Ld (mW m-2 nm-1 sr-1)")
    ax2.set_xlabel("Wavelength (nm)")

    ax3.set_xlim([wavemin_plt, wavemax_plt])
    ax3.set_ylabel("Ld/Ed (sr-1)")
    ax3.set_xlabel("Wavelength (nm)")

    ax4.set_xlim([wavemin_plt, wavemax_plt])
    ax4.set_ylabel("Eddir/Edtot (-)")
    ax4.set_xlabel("Wavelength (nm)")

    ax4.legend(title="LRT Clear sky model", loc='center left', bbox_to_anchor=(1, 0.0))

    myfig.suptitle("Relative azimuth angle=" + str(aziout * 45) + "°")
    # plt.show

    if ifpng:
        plt.savefig(outpath + "/LRT_ClearSkyModel_azi" + str(45 * aziout) + ".png", bbox_inches='tight')
    if ifpdf:
        plt.savefig(outpath + "/LRT_ClearSkyModel_azi" + str(45 * aziout) + ".pdf", bbox_inches='tight')

    wavemin_plt = 600
    wavemax_plt = 900

    ax1.set_xlim([wavemin_plt, wavemax_plt])
    ax1.set_ylabel("Edtot (mW m-2 nm-1)")
    ax1.set_xlabel("Wavelength (nm)")

    ax2.set_xlim([wavemin_plt, wavemax_plt])
    ax2.set_ylim([0, 100])
    ax2.set_ylabel("Ld (mW m-2 nm-1 sr-1)")
    ax2.set_xlabel("Wavelength (nm)")

    ax3.set_xlim([wavemin_plt, wavemax_plt])
    ax3.set_ylim([0, 0.15])
    ax3.set_ylabel("Ld/Ed (sr-1)")
    ax3.set_xlabel("Wavelength (nm)")

    ax4.set_xlim([wavemin_plt, wavemax_plt])
    ax4.set_ylabel("Eddir/Edtot (-)")
    ax4.set_xlabel("Wavelength (nm)")

    ax4.legend(title="LRT Clear sky model", loc='center left', bbox_to_anchor=(1, 0.0))
    # plt.show

    myfig.tight_layout()

    if ifpng:
        plt.savefig(outpath + "/LRT_ClearSkyModel_azi" + str(45 * aziout) + "_zoom.png", bbox_inches='tight')
    if ifpdf:
        plt.savefig(outpath + "/LRT_ClearSkyModel_azi" + str(45 * aziout) + "_zoom.pdf", bbox_inches='tight')


def main():
    print('[INFO] Started clear sky modelling...')
    # make_test()
    # get_lrt_model('/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/CLEAR_SKY_MODEL/uvspec_clear_sza40.out')

    file_hypstar = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/Comparison_Valid_2024.nc'
    from netCDF4 import Dataset
    dataset = Dataset(file_hypstar)
    wl_array = dataset.variables['HYPSTAR_Nominal_Wavelengths'][:]


    dataset.close()
    sza_list = [0,10,20,30,40,50,60,70,80]
    ckm = ClearSkyModel('/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/CLEAR_SKY_MODEL')
    ps = PlotSpectra()
    ps.xdata = wl_array
    style = {'color': 'blue', 'linestyle': '-', 'linewidth': 1, 'marker': 's', 'markersize': 0}
    from matplotlib import pyplot as plt
    import matplotlib as mpl
    cm  =  mpl.colormaps['tab20']
    colors = [cm.colors[i] for i in range(len(sza_list))]
    str_legend = []
    for idx,sza in enumerate(sza_list):
        print('Working with sza: ',sza)
        ld_array_interp,eddir_array_interp,eddif_array_interp,edtot_array_interp,ld_array = ckm.get_lrt_model_geometry_wl(sza,90,wl_array)
        style_here = style.copy()
        style_here['color'] = colors[idx]
        ps.plot_data(ld_array_interp,style_here)
        str_legend.append(f'{sza:.1f}')
    ps.set_legend(str_legend)
    ps.set_tigth_layout()
    ps.save_plot('/mnt/c/DATA_LUIS/INSITU_HYPSTAR/VEIT_HYPSTAR_AERONET_OC/CLEAR_SKY_MODEL/ele.tif')
# %%
if __name__ == '__main__':
    main()
