import numpy as np
from matplotlib import pyplot as plt


class PlotSpectra():

    def __init__(self):
        self.start_plot()
        self.xdata = []

        line_style_default = {
            'color': 'b',
            'marker': None,
            'linestyle': 'solid',
            'linewidth': 1,
            'markersize': 25
        }
        self.legend_options = {
            'loc': 'upper left',
            'bbox_to_anchor': (1.0, 1.0),
            'framealpha': 0.8
        }

        self.line_style_default = line_style_default.copy()
        self.spectra_style = line_style_default.copy()
        self.spectra_style['color'] = 'gray'

        self.stats_style = {
            'avg': line_style_default.copy(),
            'std': line_style_default.copy(),
            'minmax': line_style_default.copy(),
            'fill': {
                'color': 'gray',
                'alpha': 0.2
            }
        }
        self.stats_style['avg']['color'] = 'black'
        self.stats_style['avg']['linewidth'] = 1

        self.stats_style['std']['color'] = 'black'
        self.stats_style['std']['linewidth'] = 0
        self.stats_style['std']['linestyle'] = 'dashed'

        self.stats_style['minmax']['color'] = 'black'
        self.stats_style['minmax']['linewidth'] = 0
        self.stats_style['minmax']['linestyle'] = 'dashed'

    def start_plot(self):
        plt.close()
        plt.figure()

    def close_plot(self):
        plt.close()

    def save_plot(self, file_out):
        plt.savefig(file_out, dpi=300)

    def plot_data(self, ydata, style):
        h, = plt.plot(self.xdata, ydata,
                 color=style['color'],
                 linestyle=style['linestyle'],
                 linewidth=style['linewidth'],
                 marker=style['marker'],
                 markersize=style['markersize'])
        return h

    def plot_single_line(self, ydata, line_color, line_type, line_width, marker, marker_size):
        style = self.line_style_default
        if line_color is not None:
            style['color'] = line_color
        if line_type is not None:
            style['linestyle'] = line_type
        if line_width is not None:
            style['linewidth'] = line_width
        if marker is not None:
            style['marker'] = marker
        if marker_size is not None:
            style['markersize'] = marker_size

        h = self.plot_data(ydata, style)

        return h

    def set_legend(self, str_legend):
        plt.legend(str_legend, loc=self.legend_options['loc'], bbox_to_anchor=self.legend_options['bbox_to_anchor'], framealpha=self.legend_options['framealpha'])

    def set_legend_h(self,handles,str_legend):
        plt.legend(handles,str_legend, loc=self.legend_options['loc'], bbox_to_anchor=self.legend_options['bbox_to_anchor'],framealpha=self.legend_options['framealpha'])

    def set_title(self, title):
        plt.title(title)

    def set_xticks(self, xticks, xtickvalues, rotation, fontsize):
        if rotation is None:
            rotation = 0
        if rotation < 0 or rotation > 90:
            rotation = 0
        if xtickvalues is None:
            xtickvalues = xticks
        if fontsize is None:
            fontsize = 9
        plt.xticks(xticks, xtickvalues, rotation=rotation, fontsize=fontsize)

    def set_xaxis_title(self, xaxis_title):
        plt.xlabel(xaxis_title, fontsize=12)

    def set_yaxis_title(self, yaxis_title):
        plt.ylabel(yaxis_title, fontsize=12)

    def save_fig(self, file_out):
        plt.savefig(file_out, dpi=300)

    def set_equal_apect(self):
        plt.gca().set_aspect('equal', adjustable='box')

    def set_grid(self):
        plt.grid(b=True, which='major', color='gray', linestyle='--')

    def set_tigth_layout(self):
        plt.gcf().tight_layout()

    def set_y_range(self, ymin, ymax):
        plt.ylim(ymin, ymax)

    def plot_multiple_spectra(self, wavelength, spectra, stats, wlmin, wlmax):
        imin, imax = self.get_imin_imax_from_wavelength(wavelength,wlmin,wlmax)

        self.xdata = wavelength[imin:imax]
        if spectra is not None:
            spectra = spectra[:, imin:imax]
            self.plot_data(spectra.transpose(), self.spectra_style)

        self.plot_data(stats['avg'][imin:imax], self.stats_style['avg'])
        y1 = stats['avg'][imin:imax] - stats['std'][imin:imax]
        y2 = stats['avg'][imin:imax] + stats['std'][imin:imax]
        self.plot_data(y1, self.stats_style['std'])
        self.plot_data(y2, self.stats_style['std'])
        plt.fill_between(self.xdata, y1, y2, color=self.stats_style['fill']['color'],
                         alpha=self.stats_style['fill']['alpha'])
        self.plot_data(stats['spectra_min'][imin:imax], self.stats_style['minmax'])
        self.plot_data(stats['spectra_max'][imin:imax], self.stats_style['minmax'])

        ymin, ymax = self.get_ymin_ymax_from_stats(stats,imin,imax)
        self.set_y_range(ymin, ymax)

        return imin,imax

    def plot_stats(self,stats,imin, imax):
        h = self.plot_data(stats['avg'][imin:imax], self.stats_style['avg'])
        y1 = stats['avg'][imin:imax] - stats['std'][imin:imax]
        y2 = stats['avg'][imin:imax] + stats['std'][imin:imax]
        if self.stats_style['std']['linewidth']>0:
            self.plot_data(y1, self.stats_style['std'])
            self.plot_data(y2, self.stats_style['std'])
        plt.fill_between(self.xdata, y1, y2, facecolor=self.stats_style['fill']['color'],
                         alpha=self.stats_style['fill']['alpha'])
        if  self.stats_style['minmax']['linewidth'] > 0:
            self.plot_data(stats['spectra_min'][imin:imax], self.stats_style['minmax'])
            self.plot_data(stats['spectra_max'][imin:imax], self.stats_style['minmax'])

        return h
    def get_ymin_ymax_from_stats(self,stats,imin,imax):
        y1 = stats['avg'][imin:imax] - (2 * stats['std'][imin:imax])
        y2 = stats['avg'][imin:imax] + (2 * stats['std'][imin:imax])
        ymin = np.min(stats['spectra_min'])
        yminstd = np.min(y1)
        if yminstd > ymin:
            ymin = yminstd
        ymax = np.max(stats['spectra_max'])
        ymaxstd = np.max(y2)
        if ymaxstd < ymax:
            ymax = ymaxstd
        return ymin,ymax

    def get_imin_imax_from_wavelength(self,wavelength, wlmin, wlmax):
        imin = 0
        imax = len(wavelength) - 1
        if wlmin is not None and wavelength[imin] < wlmin < wavelength[imax]:
            imin = np.argmin(np.abs(wavelength - wlmin))
        if wlmax is not None and wavelength[imin] < wlmax < wavelength[imax]:
            imax = np.argmin(np.abs(wavelength - wlmax))
        return imin, imax