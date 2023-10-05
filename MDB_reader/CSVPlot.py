import os.path

import pandas as pd
from MDBPlotV2 import MDBPlot

class CSVPLOT:

    def __init__(self, path_csv):
        self.dataset = None
        self.VALID = False
        if os.path.exists(path_csv):
            self.VALID = True
            try:
                self.dataset = pd.read_csv(path_csv,sep=';')
            except:
                self.VALID = False


    def plot_from_options(self,options):
        plot_list = list(options.sections())
        from PlotOptions import PlotOptions
        poptions = PlotOptions(options,None)
        poptions.set_global_options()
        for plot in plot_list:
            if plot == 'GLOBAL_OPTIONS':
                continue
            options_out = poptions.get_options(plot)

            if options_out['apply']:
                self.plot_from_options_impl(options_out)

    def plot_from_options_impl(self,options_out):
        if options_out['type'] == 'scatterplot':
            self.plot_scatter_plot(options_out)
        if options_out['type'] == 'statstable':
            self.create_table_stats(options_out)

    def create_table_stats(self,options_out):
        from MDBPlotV2 import MDBPlot
        import numpy as np
        mplot = MDBPlot(None)
        xdata = np.array(self.dataset[options_out['xvar']])
        ydata = np.array(self.dataset[options_out['yvar']])

        # file_out = options_out['file_out']

        # xdata = xdata[1:100000]
        # ydata = ydata[1:100000]
        mplot.xdata = xdata
        mplot.ydata = ydata

        use_log_scale = options_out['log_scale']
        params = options_out['params']
        flag = options_out['flag']
        mplot.compute_statistics(use_log_scale)
        table,indices = mplot.start_table([flag],params)
        table = mplot.assign_table(table, indices, params, flag)
        if not options_out['file_out'] is None:
            file_out = options_out['file_out']
            table.to_csv(file_out, sep=';')


    def plot_scatter_plot(self,options_out):
        from MDBPlotV2 import MDBPlot
        import numpy as np
        mplot = MDBPlot(None)
        xdata = np.array(self.dataset[options_out['xvar']])
        ydata = np.array(self.dataset[options_out['yvar']])



        #file_out = options_out['file_out']


        # xdata = xdata[1:100000]
        # ydata = ydata[1:100000]
        mplot.xdata = xdata
        mplot.ydata = ydata
        #options_out['file_out'] = None
        mplot.plot_scatter_plot(options_out,None,-1,-1)

        # xdata = xdata[1000:2000]
        # ydata = ydata[1000:2000]
        # mplot.xdata = xdata
        # mplot.ydata = ydata
        # options_out['file_out'] = file_out
        # options_out['include_stats'] = False
        # options_out['regression_line'] = False
        # options_out['marker'] = '+'
        # mplot.plot_scatter_plot(options_out, plot, -1, -1)