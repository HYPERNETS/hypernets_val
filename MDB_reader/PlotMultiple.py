from matplotlib import pyplot as plt
from matplotlib import image as img


class PlotMultiple():

    def __init__(self):
        self.ax = None
        self.fig = None
        self.index_row = 0
        self.index_col = 0
        self.nrow = 1
        self.ncol = 1
        self.axhere = None

    def start_multiple_plot(self, nrow, ncol):
        self.nrow = nrow
        self.ncol = ncol
        # self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(8, 6), frameon=False,
        #                                  gridspec_kw={'wspace': 0, 'hspace': 0})

        # self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(15.9, 18), frameon=False,gridspec_kw={'wspace': 0, 'hspace': 0})

        # self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(16, 12), frameon=False,gridspec_kw={'wspace': 0, 'hspace': 0})
        self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(16, 12), frameon=False,gridspec_kw={'wspace': 0, 'hspace': 0})

        # self.ax.set_axis_off()
        # self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(8, 6), frameon=False,
        #                                  gridspec_kw={'wspace': 0, 'hspace': 0})

    def start_multiple_plot_advanced(self, nrow, ncol, xfigsize, yfigsize, wspace, hspace, frameon):
        self.nrow = nrow
        self.ncol = ncol
        self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(xfigsize, yfigsize), frameon=frameon,
                                         gridspec_kw={'wspace': wspace, 'hspace': hspace})
        # self.fig.set_axis_off()

    def set_text(self, x, y,s):
        plt.text(x, y, s, fontsize=8,backgroundcolor='w')
        #plt.text(x, y, s, fontsize=12)

    def plot_image(self, file_img, index_row, index_col):
        image = img.imread(file_img)
        if self.nrow > 1 and self.ncol > 1:
            axhere = self.ax[index_row, index_col]
        elif self.nrow == 1 and self.ncol > 1:
            axhere = self.ax[index_col]
        if self.nrow > 1 and self.ncol == 1:
            axhere = self.ax[index_row]

        # import numpy as np
        # image_end = np.zeros((1384,1487,3),dtype=image.dtype)
        # image_end[:] = 255
        # height = image.shape[0]
        # width = image.shape[1]
        # yini = 0
        # xini = 0
        # if height<1384:
        #     yini = int(np.floor((1384-height)/2))
        # if width<1487:
        #     xini = int(np.floor(1487-width)/2)
        # yfin = yini + height
        # xfin = xini + width
        # print(yini,yfin,xini,xfin)
        #
        # image_end[yini:yfin,xini:xfin,:] = image[:,:,:]
        image_end = image

        # axhere.set_axis_off()
        # axhere.set_frame_on
        # self.fig.add_axes(axhere)

        axhere.imshow(image_end)

        axhere.axis(False)

        # self.fig.patch.set_visible(False)
        # axhere.axis('off')

        # axhere.grid(None)
        # axhere.spines["top"].set_visible(False)
        # axhere.spines["bottom"].set_visible(False)
        # axhere.spines["right"].set_visible(False)
        # axhere.spines["left"].set_visible(False)

        # axhere.set_xticks([], color='w')
        # axhere.set_yticks([], color='w')
        #
        # axhere.xaxis.set_ticklabels([])
        # axhere.yaxis.set_ticklabels([])



    def set_global_legend(self, handles, str_legend):
        #self.fig.legend(handles, str_legend, fontsize = 8, loc='lower center', ncol=len(str_legend), markerscale=1.5,bbox_to_anchor=(0.55,0.08))
        self.fig.legend(handles, str_legend, fontsize=8, loc='lower center', ncol=len(str_legend), markerscale=1.0,
                        bbox_to_anchor=(0.55, 0.06))

    def save_fig(self, file_out):
        if file_out.endswith('.tif'):
            plt.savefig(file_out, dpi=300, bbox_inches='tight', transparency=False, facecolor='white',pil_kwargs={"compression": "tiff_lzw"})
        else:
            plt.savefig(file_out, dpi=300, bbox_inches='tight',transparency=False,facecolor='white')
        # plt.savefig(file_out, dpi=300)

    def close_plot(self):
        plt.close()
