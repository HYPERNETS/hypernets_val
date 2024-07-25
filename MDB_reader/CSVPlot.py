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

            if options_out is not None:
                self.plot_from_options_impl(options_out)

    def plot_from_options_impl(self,options_out):
        if options_out['type'] == 'scatterplot':
            self.plot_scatter_plot(options_out)
        if options_out['type'] == 'statstable':
            self.create_table_stats(options_out)
        if options_out['type'] == 'histogram':
            self.plot_histogram(options_out)






    def plot_histogram(self,options_out):
        if options_out['hvar'] is None:
            print(f'[ERROR] Variable hvar was not set for histogram')
            return

        import numpy as np
        print(options_out)
        hdata = np.array(self.dataset[options_out['hvar']]).astype(np.int32)
        if options_out['type_histo']=='int':
            int_values = options_out['int_values']
            if int_values is None:
                int_values = [int(x) for x in np.unique(hdata)]
            hticks = options_out['hticks']
            if hticks is None:
                hticks = [str(x) for x in int_values]
            xlabel = options_out['xlabel']
            ylabel = options_out['ylabel']
            ntotal = 0
            nvalues = [0] * len(int_values)
            for index,ivalue in enumerate(int_values):
                nvalues[index] = np.sum(hdata==ivalue)
                ntotal = ntotal + nvalues[index]
            print(ntotal,len(hdata))

        from PlotSpectra import PlotSpectra
        pspectra = PlotSpectra()
        pspectra.xdata = np.array(int_values).astype(np.int32)
        pspectra.plot_single_bar_series(np.array(nvalues).astype(np.int32),'blue',1,0)
        pspectra.set_xticks(pspectra.xdata,hticks,45,12)
        pspectra.set_yaxis_title(ylabel)
        pspectra.set_xaxis_title(xlabel)
        pspectra.set_tigth_layout()
        pspectra.save_plot(options_out['file_out'])

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
        use_rhow = options_out['use_rhow']
        params = options_out['params']
        flag = options_out['flag']
        mplot.compute_statistics(use_log_scale,use_rhow)
        table,indices = mplot.start_table([flag],params)
        table = mplot.assign_table(table, indices, params, flag)
        if not options_out['file_out'] is None:
            file_out = options_out['file_out']
            table.to_csv(file_out, sep=';')


    def plot_scatter_plot(self,options_out):
        from MDBPlotV3 import MDBPlot
        import numpy as np
        mplot = MDBPlot(None)
        xdata = np.array(self.dataset[options_out['xvar']])
        ydata = np.array(self.dataset[options_out['yvar']])
        valid = np.logical_and(xdata!=-999.0,ydata!=-999.0)
        xdata = xdata[valid]
        ydata = ydata[valid]
        # xdata = xdata[1:100000]
        # ydata = ydata[1:100000]
        mplot.xdata = xdata
        mplot.ydata = ydata


        if options_out['groupBy'] is not None:
            groupdata = np.array(self.dataset[options_out['groupBy']])
            gData = groupdata[valid]
            mplot.groupdata = gData
            if options_out['groupBy']=='blended_dominant_owt':
                values = list(range(1,19))
                #values = [-1] + values
                flags = [f'OWT_{x}' for x in values]
                options_out['blended_dominant_owt'] = {
                    'flag_values': values,
                    'flag_meanings': flags
                }
                if options_out['groupValues'] is not None:
                    gvalues = options_out['groupValues']
                    options_out['groupValues'] = [int(x) for x in gvalues]
                else:
                    options_out['groupValues'] = values

                options_out['groupValues'] = [-1]
                #options_out['groupValues'] = [1,2,3,6,9,13]

                #print('==========>',options_out['groupValues'])

                selectData = np.zeros(xdata.shape)
                #gData = mplot.groupdata
                for gvalue in options_out['groupValues']:
                    selectData[gData==gvalue] = 1

                mplot.xdata = xdata[selectData==1]
                mplot.ydata = ydata[selectData==1]
                mplot.groupdata = gData[selectData==1]


                # title = options_out['title']
                # wls = title.split(' ')[0]
                # vals = {
                #     '400': {
                #         'min_xy': 2,
                #         'max_xy': 10,
                #         'ticks': [2, 4, 6, 8, 10]
                #     },
                #     '412_5': {
                #          'min_xy': 2,
                #         'max_xy': 10,
                #         'ticks': [2, 4, 6, 8, 10]
                #     },
                #     '442_5': {
                #         'min_xy': 2,
                #         'max_xy': 10,
                #         'ticks': [2, 4, 6, 8, 10]
                #     },
                #     '490': {
                #         'min_xy': 4,
                #         'max_xy': 12,
                #         'ticks': [4, 6, 8, 10,12]
                #     },
                #     '510':{
                #         'min_xy': 4,
                #         'max_xy': 12,
                #         'ticks': [4, 6, 8, 10,12]
                #     },
                #     '560': {
                #         'min_xy': 2,
                #         'max_xy': 8,
                #         'ticks': [2, 4, 6, 8]
                #     },
                #     '620': {
                #         'min_xy': 0,
                #         'max_xy': 10,
                #         'ticks': [0, 2, 4, 6,8,10]
                #     },
                #     '665': {
                #          'min_xy': 0,
                #         'max_xy': 10,
                #         'ticks': [0, 2, 4, 6,8,10]
                #     },
                #     '673_75': {
                #         'min_xy': 0,
                #         'max_xy': 10,
                #         'ticks': [0, 2, 4, 6,8,10]
                #     },
                #     '681_25': {
                #         'min_xy': 0,
                #         'max_xy': 10,
                #         'ticks': [0, 2, 4, 6,8,10]
                #     },
                #     '708_75': {
                #         'min_xy': 0,
                #         'max_xy': 3,
                #         'ticks': [0, 1, 2, 3]
                #     },
                #     '753_75': {
                #         'min_xy':-0.25,
                #         'max_xy': 1,
                #         'ticks': [-0.25,0,0.25,0.5,0.75,1]
                #     },
                #     '778_75': {
                #         'min_xy':-0.25,
                #         'max_xy': 1,
                #         'ticks': [-0.25,0,0.25,0.5,0.75,1]
                #     },
                #     '865':{
                #          'min_xy':-0.2,
                #         'max_xy': 0.5,
                #         'ticks': [-0.2,-0.1,0, 0.1, 0.2,0.3,0.4,0.5]
                #     },
                #     '885': {
                #          'min_xy':-0.3,
                #         'max_xy': 0.2,
                #         'ticks': [-0.3,-0.2,-0.1,0, 0.1, 0.2]
                #     },
                #     '1020': {
                #          'min_xy':-1,
                #         'max_xy': 1,
                #         'ticks': [-1,-0.5, 0, 0.5, 1]
                #     }
                # }
                # if wls in vals.keys():
                #     options_out['min_xy'] = vals[wls]['min_xy']
                #     options_out['max_xy'] = vals[wls]['max_xy']
                #     options_out['ticks'] = vals[wls]['ticks']
                # print('-------------->',wls,options_out['max_xy'])

        mplot.plot_scatter_plot(options_out,None,-1,-1,-1)


    def compute_statistics(self,options_out):
        from MDBPlotV2 import MDBPlot
        import numpy as np
        mplot = MDBPlot(None)
        xdata = np.array(self.dataset[options_out['xvar']])
        ydata = np.array(self.dataset[options_out['yvar']])
        selectData = np.array(self.dataset[options_out['selectBy']]).astype(np.int32)
        selectValue = options_out['selectValue']
        valid = np.logical_and(xdata!=-999.0, ydata!=-999.0)
        xdata = xdata[valid]
        ydata = ydata[valid]
        selectData = selectData[valid]
        mplot.xdata = xdata[selectData==selectValue]
        mplot.ydata = ydata[selectData==selectValue]
        mplot.compute_statistics(options_out['log_scale'],False)
        #print(mplot.valid_stats)
        return mplot.valid_stats

    def compute_distribution(self,options_out):
        #from MDBPlotV2 import MDBPlot
        import numpy as np
        #mplot = MDBPlot(None)
        ydata = np.array(self.dataset[options_out['yvar']])
        selectData = np.array(self.dataset[options_out['selectBy']]).astype(np.int32)
        selectValue = options_out['selectValue']
        ydata = ydata[selectData==selectValue]
        median = np.median(ydata)
        p25 = np.percentile(ydata,25)
        p75 = np.percentile(ydata,75)
        return median,p25,p75

