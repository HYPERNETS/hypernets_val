from MDBFile import MDBFile
from PlotOptions import PlotOptions
import os

class MDBPlot:

    def __init__(self, path_mdbr_file):

        self.mrfile = None
        self.path_mdbr_file = path_mdbr_file
        self.VALID = False
        if path_mdbr_file is not None:
            self.mrfile = MDBFile(path_mdbr_file)
            self.VALID = self.mrfile.VALID


    def plot_from_options_file(self,file_config):
        if not os.path.isfile(file_config):
            print(f'[ERROR] File config {file_config} is not a valid config file')
            return
        import configparser
        try:
            options = configparser.ConfigParser()
            options.read(file_config)
        except:
            print(f'[ERROR] Error reading file_config: {file_config}')
        self.plot_from_options(options)

    def plot_from_options(self,options):
        poptions = PlotOptions(options,None)
        poptions.set_global_options()
        print(poptions.global_options)
        # list_figures = poptions.get_list_figures()
        # for figure in list_figures:
        #     options_figure = poptions.get_options(figure)

