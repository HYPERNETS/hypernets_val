from matplotlib import pyplot as plt


class PlotSpectra():

    def __init__(self):
        self.start_plot()
        self.xdata = []

        self.line_style_default = {
            'color': 'b',
            'marker': None,
            'linestyle': '-',
            'linewidth': 1,
        }
        self.legend_options = {
            'loc': 'upper left',
            'bbox_to_anchor': (1.0, 1.0)
        }

    def start_plot(self):
        plt.close()
        plt.figure()
    def close_plot(self):
        plt.close()

    def save_plot(self,file_out):
        plt.savefig(file_out, dpi=300)

    def plot_data(self, ydata, style):
        plt.plot(self.xdata, ydata,
                 color=style['color'],
                 linestyle=style['linestyle'],
                 linewidth=style['linewidth'],
                 marker=style['marker'])

    def plot_single_line(self, ydata, line_color, line_type, line_width, marker):
        style = self.line_style_default
        if line_color is not None:
            style['color'] = line_color
        if line_type is not None:
            style['linestyle'] = line_type
        if line_width is not None:
            style['linewidth'] = line_width
        if marker is not None:
            style['marker'] = marker

        self.plot_data(ydata, style)

    def set_legend(self,str_legend):
        plt.legend(str_legend, loc=self.legend_options['loc'], bbox_to_anchor=self.legend_options['bbox_to_anchor'])

    def set_title(self,title):
        plt.title(title)

    def set_xaxis_title(self,xaxis_title):
        plt.xlabel(xaxis_title,fontsize=12)

    def set_yaxis_title(self,yaxis_title):
        plt.ylabel(yaxis_title,fontsize=12)

    def save_fig(self,file_out):
        plt.savefig(file_out, dpi=300)

    def set_equal_apect(self):
        plt.gca().set_aspect('equal', adjustable='box')

    def set_grid(self):
        plt.grid(b=True, which='major', color='gray', linestyle='--')

    def set_tigth_layout(self):
        plt.gcf().tight_layout()