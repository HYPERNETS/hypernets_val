import MDBPlotDefaults as defaults
from MDBFile import MDBFile


class MDBPlot:

    def __init__(self, path_mdbr_file):
        self.mrfile = MDBFile(path_mdbr_file)

        self.VALID = self.mrfile.VALID

        # plot general options
        self.title = ''
        self.file_name_base = ''
        self.format_image = 'jpg'

        # validation scatterplot
        self.xdata = []
        self.ydata = []
        self.wldata = []
        self.xregress = []
        self.yregress = []
        self.xlabel = defaults.xlabel_default
        self.ylabel = defaults.ylabel_default

        # validation stats
        self.valid_stats = {
            'N': 0,
            'slope': 0.0,
            'intercept': 0.0,
            'r_value': 0.0,
            'p_value': 0.0,
            'std_err': 0.0,
            'rmse_val': 0.0,
            'mean_rel_diff': 0.0,
            'mean_abs_rel_diff': 0.0,
            'bias': 0.0,
            'r2': 0.0,
            'slope_typeII': 0.0,
            'offset_typeII': 0.0,
            'XAVG': 0.0,
            'YAVG': 0.0,
            'CPRMSE': 0.0,
            'MAE': 0.0
        }
        self.df_valid_stats = None

        self.units = r'sr$^-$$^1$'
        self.log_scale = False

        self.satellite = 'S3'
        self.platform = 'AB'

    def plot_from_options(self,options):
        plot_list = list(options.sections())
        for plot in plot_list:
            if options[plot]['apply']=='true':
                self.plot_from_options_impl(plot,options[plot])

    def plot_from_options_impl(self,plot_name,options):
        options_dict = dict(options)
        print(options_dict)
        if options_dict['type']=='scatterplot':
            plot_scatter_plot()