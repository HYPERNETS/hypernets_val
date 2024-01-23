from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
#import seaborn as sns
from matplotlib.ticker import FormatStrFormatter


class PlotScatter():

    def __init__(self):
        self.ax = None
        self.fig = None
        self.index_row = 0
        self.index_col = 0
        self.nrow = 1
        self.ncol = 1
        self.axhere = None
        self.start_plot()

        self.style_default = {
            'c': 'k',
            's': None,
            'marker': '.',
            'linewidths': None,
            'edgecolors': None,
            'markeredgecolor': None,
            'markeredgewidth': 0
        }

        self.xtitle_options = {
            'fontsize': 12
        }

        self.ytitle_options = {
            'fontsize': 12
        }

        self.legend_options = {
            'loc': 'upper left',
            'bbox_to_anchor': (1.0, 1.0)
        }

        self.plot_text_options = {
            'horizontalalignment': 'left',
            'fontsize': 12
        }

    def start_plot(self):
        self.fig, self.ax = plt.subplots()

    def start_multiple_plot(self,nrow,ncol):
        self.nrow = nrow
        self.ncol = ncol
        #self.fig, self.ax = plt.subplots(nrow, ncol,figsize=(7,7),gridspec_kw={'wspace':0.1,'hspace':0.1})
        self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(7, 4), gridspec_kw={'wspace': 0.1, 'hspace': 0.1})

    def start_multiple_plot_advanced(self,nrow,ncol,xfigsize,yfigsize,widthspace,heightspace):
        self.nrow = nrow
        self.ncol = ncol
        # self.fig, self.ax = plt.subplots(nrow, ncol,figsize=(7,7),gridspec_kw={'wspace':0.1,'hspace':0.1})
        self.fig, self.ax = plt.subplots(nrow, ncol, figsize=(xfigsize, yfigsize), gridspec_kw={'wspace': widthspace, 'hspace': heightspace})

    def close_plot(self):
        plt.close()

    def set_axhere(self):
        self.axhere = self.ax
        if self.nrow > 1 and self.ncol == 1:
            self.axhere = self.ax[self.index_row]
        elif self.nrow == 1 and self.ncol > 1:
            self.axhere = self.ax[self.index_col]
        elif self.nrow > 1 and self.ncol > 1:
            self.axhere = self.ax[self.index_row, self.index_col]

    def set_axhere_rc(self,index_row,index_col):
        self.axhere = self.ax[index_row,index_col]

    def set_axhere_index(self,index):
        index_row = np.floor(index/self.ncol)
        index_col = index-(index_row*self.ncol)
        self.index_row = int(index_row)
        self.index_col = int(index_col)
        self.set_axhere_rc(self.index_row,self.index_col)

    def set_log_scale(self):
        if self.axhere is None:
            self.set_axhere()
        self.axhere.set_yscale('log')
        self.axhere.set_xscale('log')
        self.axhere.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        self.axhere.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    def plot_data(self, xdata, ydata, marker, markersize, color, edgecolor, linewidth):
        style = self.style_default
        if marker is not None:
            style['marker'] = marker
        if markersize is not None:
            style['s'] = markersize
        if color is not None:
            style['c'] = color
        if edgecolor is not None:
            style['edgecolors'] = edgecolor
        if linewidth is not None:
            style['linewidths'] = linewidth
        if self.axhere is None:
            self.set_axhere()
        self.axhere.scatter(xdata, ydata,
                            marker=style['marker'],
                            s=style['s'],
                            c=style['c'],
                            edgecolors=style['edgecolors'],
                            linewidths=style['linewidths'],
                            alpha = 1.0)

    # def plot_reg_line(self, xdata, ydata, color):
    #     data_plot = pd.concat([xdata, ydata], axis=1).astype(dtype=np.float)
    #     sns.lmplot(data=data_plot, x='Ins_Rrs', y='Sat_Rrs', line_kws={'color': color})

    def colorbar(self):
        if self.axhere is None:
            self.set_axhere()
        self.axhere.colorbar()

    def set_cmap(self, cmap):
        plt.set_cmap(cmap)
        # if self.axhere is None:
        #     self.set_axhere()
        # self.axhere.set_cmap(cmap)

    def set_equal_apect(self):
        if self.axhere is None:
            self.set_axhere()
        self.axhere.set_aspect('equal', adjustable='box')


    # def tight_layout(self):
    #     if self.axhere is None:
    #         self.set_axhere()
    #     self.fig.tight_layout()

    def set_xaxis_title(self, xaxis_title):
        if self.axhere is None:
            self.set_axhere()
        xaxis_title = xaxis_title.replace('$R$',u'\u00AE')
        self.axhere.set_xlabel(xaxis_title, fontsize=self.xtitle_options['fontsize'])

    def set_yaxis_title(self, yaxis_title):
        if self.axhere is None:
            self.set_axhere()
        self.axhere.set_ylabel(yaxis_title, fontsize=self.ytitle_options['fontsize'])

    def set_limits(self, minV, maxV):
        if self.axhere is None:
            self.set_axhere()
        self.axhere.set_xlim([minV, maxV])
        self.axhere.set_ylim([minV, maxV])

    def set_legend(self, str_legend):
        if self.axhere is None:
            self.set_axhere()
        self.axhere.legend(str_legend, loc=self.legend_options['loc'],
                           bbox_to_anchor=self.legend_options['bbox_to_anchor'])

    def set_global_legend(self,str_legend):
        #self.fig.legend(str_legend, loc='lower center', ncol=len(str_legend),markerscale=2.0,bbox_to_anchor=(0.5,0.04))
        #self.fig.legend(str_legend, loc='lower center', ncol=len(str_legend), markerscale=2.0)
        #self.fig.legend(str_legend, loc='upper center', ncol=len(str_legend))
        if len(str_legend)==8:

            self.fig.legend(str_legend, fontsize=11, loc='lower center', ncol=4, markerscale=1.5,
                            bbox_to_anchor=(0.5, -0.015))
        else:

            self.fig.legend(str_legend, fontsize=11, loc='lower center', ncol=len(str_legend), markerscale=1.5,bbox_to_anchor=(0.5, -0.015))

    def set_title(self, title):
        if self.axhere is None:
            self.set_axhere()
        self.axhere.set_title(title)

    def plot_identity_line(self):
        if self.axhere is None:
            self.set_axhere()
        xmin, xmax = self.axhere.get_xlim()
        ymin, ymax = self.axhere.get_ylim()
        xmin = np.min([xmin, ymin])
        xmax = np.max([xmax, ymax])
        self.axhere.plot([xmin, xmax], [xmin, xmax], '--k', linewidth=0.75)

    def plot_regress_line(self, xdata, ydata, color):
        if self.axhere is None:
            self.set_axhere()
        self.axhere.plot(xdata, ydata, color=color, linestyle='-', linewidth=2, marker=None)

    def plot_text(self, xpos, ypos, str):
        if self.axhere is None:
            self.set_axhere()
        self.axhere.text(xpos, ypos, str,
                         horizontalalignment=self.plot_text_options['horizontalalignment'],
                         fontsize=self.plot_text_options['fontsize'],
                         transform=self.axhere.transAxes)

    def set_ticks(self,ticks,fontsize):
        if self.axhere is None:
            self.set_axhere()
        self.axhere.set_xticks(ticks)
        self.axhere.set_yticks(ticks)
        if fontsize>0:
            self.axhere.tick_params(axis='both',labelsize=fontsize)

    def set_ticks_and_labels(self,ticks,labels,fontsize):
        if self.axhere is None:
            self.set_axhere()
        self.axhere.set_xticks(ticks,labels=labels)
        self.axhere.set_yticks(ticks,labels=labels)
        if fontsize>0:
            self.axhere.tick_params(axis='both',labelsize=fontsize)

    def set_xticks_labels_off(self,ticks):
        if self.axhere is None:
            self.set_axhere()
        if ticks is not None:
            self.axhere.set_xticks(ticks)
            self.axhere.set_yticks(ticks)
        self.axhere.tick_params(axis='x',labelbottom=False,labeltop=False)

    def set_yticks_labels_off(self,ticks):
        if self.axhere is None:
            self.set_axhere()
        if ticks is not None:
            self.axhere.set_xticks(ticks)
            self.axhere.set_yticks(ticks)
        self.axhere.tick_params(axis='y',labelleft=False,labelright=False)

    def plot_blanck(self,index):
        if index>=0:
            self.set_axhere_index(index)
        if self.axhere is None:
            self.set_axhere()
        self.axhere.axis('off')

    def save_fig(self, file_out):
        if file_out.endswith('.tif'):
            plt.savefig(file_out, dpi=300, bbox_inches='tight',pil_kwargs={"compression": "tiff_lzw"})
        else:
            plt.savefig(file_out, dpi=300, bbox_inches='tight')
        #plt.savefig(file_out,dpi = 300, bbox_inches = 'tight')
