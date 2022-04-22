from matplotlib import pyplot as plt
import numpy as np

class PlotScatter():

    def __init__(self):
        self.start_plot()

        self.style_default = {
            'c': 'k',
            's': None,
            'marker': 'o',
            'linewidths': None,
            'edgecolors': None
        }

        self.xtitle_options = {
            'fontsize' : 12
        }

        self.ytitle_options = {
            'fontsize' : 12
        }

        self.legend_options = {
            'loc':'upper left',
            'bbox_to_anchor': (1.0, 1.0)
        }

        self.plot_text_options = {
            'horizontalalignment':'left',
            'fontsize':12
        }

    def start_plot(self):
        plt.figure()

    def close_plot(self):
        plt.close()

    def plot_data(self,xdata,ydata,marker,markersize,color,edgecolor,linewidth):
        style = self.style_default
        if marker is not None:
            style['marker']=marker
        if markersize is not None:
            style['s'] = markersize
        if color is not None:
            style['c'] = color
        if edgecolor is not None:
            style['edgecolors'] = edgecolor
        if linewidth is not None:
            style['linewidths'] = linewidth

        plt.scatter(xdata, ydata,
                    marker = style['marker'],
                    s = style['s'],
                    c = style['c'],
                    edgecolors=style['edgecolors'],
                    linewidths=style['linewidths'])

    def set_equal_apect(self):
        plt.gca().set_aspect('equal', adjustable='box')

    def set_xaxis_title(self,xaxis_title):
        plt.xlabel(xaxis_title,fontsize=self.xtitle_options['fontsize'])

    def set_yaxis_title(self,yaxis_title):
        plt.ylabel(yaxis_title,fontsize=self.ytitle_options['fontsize'])

    def set_legend(self,str_legend):
        plt.legend(str_legend, loc=self.legend_options['loc'], bbox_to_anchor=self.legend_options['bbox_to_anchor'])

    def set_title(self,title):
        plt.title(title)

    def plot_identity_line(self):
        xmin, xmax = plt.gca().get_xlim()
        ymin, ymax = plt.gca().get_ylim()
        xmin = np.min([xmin, ymin])
        xmax = np.max([xmax, ymax])
        plt.plot([xmin, xmax], [xmin, xmax], '--k')

    def plot_text(self,xpos, ypos, str):
        plt.text(xpos, ypos, str,
                 horizontalalignment=self.plot_text_options['horizontalalignment'],
                 fontsize=self.plot_text_options['fontsize'],
                 transform=plt.gca().transAxes)

    def save_fig(self,file_out):
        plt.savefig(file_out, dpi=300)