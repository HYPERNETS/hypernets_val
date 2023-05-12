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

        #self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(15.9, 18), frameon=False,gridspec_kw={'wspace': 0, 'hspace': 0})

        #self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(16, 12), frameon=False,gridspec_kw={'wspace': 0, 'hspace': 0})
        self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(16, 12), frameon=False,gridspec_kw={'wspace': 0, 'hspace': 0})
        # self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(8, 6), frameon=False,
        #                                  gridspec_kw={'wspace': 0, 'hspace': 0})
    def set_text(self, str_text, x, y):
        plt.text(x, y, str_text, fontsize=10)

    def plot_image(self, file_img, index_row, index_col):
        image = img.imread(file_img)
        if self.nrow > 1 and self.ncol > 1:
            axhere = self.ax[index_row, index_col]
        elif self.nrow == 1 and self.ncol > 1:
            axhere = self.ax[index_col]
        if self.nrow > 1 and self.ncol == 1:
            axhere = self.ax[index_row]
        axhere.imshow(image)

        # axhere.axis('off')
        axhere.set_xticks([], color='w')
        axhere.set_yticks([], color='w')

        axhere.xaxis.set_ticklabels([])
        axhere.yaxis.set_ticklabels([])

    def save_fig(self, file_out):
        plt.savefig(file_out, dpi=300, bbox_inches='tight')

    def close_plot(self):
        plt.close()
