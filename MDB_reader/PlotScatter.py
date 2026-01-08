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
        self.ax_multiple = None
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
            'bbox_to_anchor': (1.0, 1.0),
            'framealpha': 0.8,
            'ncols': 1,
            'markerscale': 1
        }

        self.plot_text_options = {
            'horizontalalignment': 'left',
            'fontsize': 12
        }

    def prepare_poster(self):
        import matplotlib.ticker as ticker
        self.ax.tick_params(axis='both',which='major',labelsize=14)
        plt.xlabel(self.ax.get_xlabel(), fontsize=16)
        plt.ylabel(self.ax.get_ylabel(), fontsize=16)
        self.ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
        self.ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))


    def set_grid(self):
        # plt.grid(b=True, which='major', color='gray', linestyle='--')
        plt.grid(which='major', color='gray', linestyle='--', axis='both')


    def start_plot(self):
        self.fig, self.ax = plt.subplots()

    def start_plot_joliff(self):
        self.fig, self.ax = plt.subplots(figsize=(8,8))

    def start_plot_polar(self):
        self.fig, self.ax = plt.subplots(figsize=(6, 6), subplot_kw={'projection': 'polar'})

    def start_multiple_plot(self,nrow,ncol):
        self.nrow = nrow
        self.ncol = ncol
        #self.fig, self.ax = plt.subplots(nrow, ncol,figsize=(7,7),gridspec_kw={'wspace':0.1,'hspace':0.1})
        self.fig, self.ax_multiple = plt.subplots(nrow, ncol, figsize=(7, 4), gridspec_kw={'wspace': 0.1, 'hspace': 0.1})

    def start_multiple_plot_advanced(self,nrow,ncol,xfigsize,yfigsize,widthspace,heightspace):
        self.nrow = nrow
        self.ncol = ncol
        # self.fig, self.ax = plt.subplots(nrow, ncol,figsize=(7,7),gridspec_kw={'wspace':0.1,'hspace':0.1})
        self.fig, self.ax_multiple = plt.subplots(nrow, ncol, figsize=(xfigsize, yfigsize), gridspec_kw={'wspace': widthspace, 'hspace': heightspace})

    def close_plot(self):
        plt.close()

    def set_axhere(self):
        self.ax = self.ax
        if self.nrow > 1 and self.ncol == 1:
            self.ax = self.ax_multiple[self.index_row]
        elif self.nrow == 1 and self.ncol > 1:
            self.ax = self.ax_multiple[self.index_col]
        elif self.nrow > 1 and self.ncol > 1:
            self.ax = self.ax_multiple[self.index_row, self.index_col]

    def set_axhere_rc(self,index_row,index_col):
        self.ax = self.ax_multiple[index_row,index_col]

    def set_axhere_index(self,index):
        index_row = np.floor(index/self.ncol)
        index_col = index-(index_row*self.ncol)
        self.index_row = int(index_row)
        self.index_col = int(index_col)
        #print(index,self.index_row,self.index_col)
        self.set_axhere_rc(self.index_row,self.index_col)

    def set_log_scale(self):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_yscale('log')
        self.ax.set_xscale('log')
        self.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        self.ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    def plot_data(self, xdata, ydata, marker, markersize, color, edgecolor, linewidth):
        if len(color.split(';'))==3:
            color = tuple([float(x.strip()) for x in color.split(';')])
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
        if self.ax is None:
            self.set_axhere()



        hscatter = self.ax.scatter(xdata, ydata,
                            marker=style['marker'],
                            s=style['s'],
                            c=style['c'],
                            edgecolors=style['edgecolors'],
                            linewidths=style['linewidths'],
                            alpha = 1.0)
        return hscatter
    # def plot_reg_line(self, xdata, ydata, color):
    #     data_plot = pd.concat([xdata, ydata], axis=1).astype(dtype=np.float)
    #     sns.lmplot(data=data_plot, x='Ins_Rrs', y='Sat_Rrs', line_kws={'color': color})

    def colorbar(self,hscatter):
        # if self.ax is None:
        #     self.set_axhere()
        import matplotlib.pyplot as plt
        plt.colorbar(hscatter)
        #self.ax.colorbar()

    def set_cmap(self, cmap):
        plt.set_cmap(cmap)
        # if self.ax is None:
        #     self.set_axhere()
        # self.ax.set_cmap(cmap)

    def set_equal_apect(self):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_aspect('equal', adjustable='box')

    def set_label_fontsize(self,fontsize):
        self.ax.tick_params(axis='x', labelsize=fontsize)
        self.ax.tick_params(axis='y', labelsize=fontsize)

    def tight_layout(self):
        plt.tight_layout()

    def set_xaxis_title(self, xaxis_title):
        if self.ax is None:
            self.set_axhere()
        xaxis_title = xaxis_title.replace('$R$',u'\u00AE')
        self.ax.set_xlabel(xaxis_title, fontsize=self.xtitle_options['fontsize'])

    def set_yaxis_title(self, yaxis_title):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_ylabel(yaxis_title, fontsize=self.ytitle_options['fontsize'])

    def set_limits(self, minV, maxV):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_xlim([minV, maxV])
        self.ax.set_ylim([minV, maxV])

    def set_limits_X(self,minV,maxV):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_xlim([minV, maxV])

    def set_limits_Y(self,minV,maxV):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_ylim([minV, maxV])

    def set_legend(self, str_legend):
        print(self.legend_options)
        plt.legend(str_legend, loc=self.legend_options['loc'], bbox_to_anchor=self.legend_options['bbox_to_anchor'],
                   framealpha=self.legend_options['framealpha'], ncol=self.legend_options['ncols'],markerscale=self.legend_options['markerscale'])

    def set_legend_h(self, handles, str_legend):
        plt.legend(handles, str_legend, loc=self.legend_options['loc'],
                   bbox_to_anchor=self.legend_options['bbox_to_anchor'], framealpha=self.legend_options['framealpha'],
                   ncol=self.legend_options['ncols'],markerscale=self.legend_options['markerscale'])

    def set_global_legend(self,str_legend):
        #self.fig.legend(str_legend, loc='lower center', ncol=len(str_legend),markerscale=2.0,bbox_to_anchor=(0.5,0.04))
        #self.fig.legend(str_legend, loc='lower center', ncol=len(str_legend), markerscale=2.0)
        #self.fig.legend(str_legend, loc='upper center', ncol=len(str_legend))
        if len(str_legend)==8:

            self.fig.legend(str_legend, fontsize=11, loc='lower center', ncol=4, markerscale=1.5,
                            bbox_to_anchor=(0.5, -0.015))
        else:
            #self.fig.legend(str_legend, fontsize=11, loc='lower center', ncol=len(str_legend), markerscale=1.5,bbox_to_anchor=(0.5, -0.015))
            self.fig.legend(str_legend, fontsize=11, loc='lower center', ncol=len(str_legend), markerscale=1.5,bbox_to_anchor=(0.5, 0.015))
    def set_title(self, title):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_title(title)

    def plot_identity_line(self):
        if self.ax is None:
            self.set_axhere()
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        xmin = np.min([xmin, ymin])
        xmax = np.max([xmax, ymax])
        self.ax.plot([xmin, xmax], [xmin, xmax], '--k', linewidth=0.75)

    def plot_regress_line(self, xdata, ydata, color):
        if self.ax is None:
            self.set_axhere()
        self.ax.plot(xdata, ydata, color=color, linestyle='-', linewidth=2, marker=None)

    def plot_regress_line_options(self, xdata, ydata, color,linestyle,linewidth):
        if self.ax is None:
            self.set_axhere()
        self.ax.plot(xdata, ydata, color=color, linestyle=linestyle, linewidth=linewidth, marker=None)



    def plot_text(self, xpos, ypos, str):
        if self.ax is None:
            self.set_axhere()
        self.ax.text(xpos, ypos, str,
                         horizontalalignment=self.plot_text_options['horizontalalignment'],
                         fontsize=self.plot_text_options['fontsize'],
                         transform=self.ax.transAxes)

    def set_ticks(self,ticks,fontsize):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_xticks(ticks)
        self.ax.set_yticks(ticks)
        if fontsize>0:
            self.ax.tick_params(axis='both',labelsize=fontsize)

    def set_ticks_x(self, ticks, fontsize):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_xticks(ticks)
        if fontsize > 0:
            self.ax.tick_params(axis='x', labelsize=fontsize)

    def set_ticks_y(self, ticks, fontsize):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_yticks(ticks)
        if fontsize > 0:
            self.ax.tick_params(axis='y', labelsize=fontsize)

    def set_ticks_and_labels(self,ticks,labels,fontsize):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_xticks(ticks,labels=labels)
        self.ax.set_yticks(ticks,labels=labels)
        if fontsize>0:
            self.ax.tick_params(axis='both',labelsize=fontsize)

    def set_ticks_and_labels_x(self,ticks,labels,fontsize):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_xticks(ticks,labels=labels)
        if fontsize>0:
            self.ax.tick_params(axis='x',labelsize=fontsize)

    def set_ticks_and_labels_y(self,ticks,labels,fontsize):
        if self.ax is None:
            self.set_axhere()
        self.ax.set_yticks(ticks,labels=labels)
        if fontsize>0:
            self.ax.tick_params(axis='y',labelsize=fontsize)

    def set_xticks_labels_off(self,ticks):
        if self.ax is None:
            self.set_axhere()
        if ticks is not None:
            self.ax.set_xticks(ticks)
            self.ax.set_yticks(ticks)
        self.ax.tick_params(axis='x',labelbottom=False,labeltop=False)

    def set_yticks_labels_off(self,ticks):
        if self.ax is None:
            self.set_axhere()
        if ticks is not None:
            self.ax.set_xticks(ticks)
            self.ax.set_yticks(ticks)
        self.ax.tick_params(axis='y',labelleft=False,labelright=False)


    def set_as_joliff(self,x_ticks):
        self.ax.tick_params(axis='x', labelbottom=False, labeltop=False)
        self.ax.tick_params(axis='y', labelsize=10)

        ##Show circular spines
        if x_ticks is not None:
            for cspine in x_ticks:
                self.ax.add_patch(plt.Circle((0.0, 0.0), cspine, fill=False, facecolor=None, edgecolor='lightgray', linestyle='--',linewidth=0.75))

        self.ax.spines['left'].set_position('zero')
        self.ax.spines['bottom'].set_position('zero')

        # Eliminate upper and right axes
        self.ax.spines['right'].set_color('none')
        self.ax.spines['top'].set_color('none')

    def plot_blanck(self,index):
        if index>=0:
            self.set_axhere_index(index)
        if self.ax is None:
            self.set_axhere()
        self.ax.axis('off')



    ##methods for polar plots
    def plot_polar(self,angle, radius,convert_to_rad, marker, markersize, color,edgecolor, linewidth):
        if convert_to_rad:
            angle = np.deg2rad(angle)
        return self.plot_data(angle,radius,marker,markersize,color,edgecolor,linewidth)

    def set_rscale(self,rscale):
       self.ax.set_rscale(rscale)

    def set_rlim(self,rlim):
        self.ax.set_rlim(rlim)

    def set_rticks(self,rticks):
        self.ax.set_rticks(rticks)

    def set_rticks_and_labels(self,rticks,rlabels,use_minor):
        return self.ax.set_rticks(rticks,rlabels,minor=use_minor)

    def set_theta_zero_location(self,tzl):
        self.ax.set_theta_zero_location(tzl)

    def set_theta_direction(self,td):
        self.ax.set_theta_direction(td)

    def set_theta_range(self,tmin,tmax):
        self.ax.set_thetamin(tmin)
        self.ax.set_thetamax(tmax)

    def set_polar_ticks(self,pticks,convert_to_rad):
        if convert_to_rad:
            pticks = np.deg2rad(pticks)
        self.ax.set_xticks(pticks)

    def set_rlabel_position(self,rpos):
        self.ax.set_rlabel_position(rpos)

    def set_text_size(self, x, y, s,fontsize):
        plt.text(x, y, s, fontsize=fontsize, backgroundcolor='w',transform=plt.gcf().transFigure)


    def save_fig(self, file_out):
        if file_out.endswith('.tif'):
            plt.savefig(file_out, dpi=300, bbox_inches='tight',pil_kwargs={"compression": "tiff_lzw"})
        else:
            plt.savefig(file_out, dpi=300, bbox_inches='tight')
        #plt.savefig(file_out,dpi = 300, bbox_inches = 'tight')
