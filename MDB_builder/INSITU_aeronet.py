import __init__,os,sys
try:
    from INSITU_base import INSITUBASE
except:
    from MDB_builder.INSITU_base import INSITUBASE

code_home = os.path.dirname(os.path.dirname(__init__.__file__))
sys.path.append(code_home)
code_aeronet = os.path.join(os.path.dirname(code_home), 'aeronet')
sys.path.append(code_aeronet)

class INSITU_AERONET:

    def __init__(self,insitu_options,verbose):

        self.path_aeronet_nc = insitu_options['path_aeronet']
        self.verbose = verbose
        #self.insitu_options = insitu_options
        self.site = insitu_options['site']
        self.file_list = {}
        self.date_list = []
        self.start_date = None
        self.end_date = None


        self.VALID = self.check_path_aeronet_nc()


    def check_path_aeronet_nc(self):
        try:
            from base.anet_nc_reader import AERONETReader
        except Exception as ex:
            print(f'[ERROR] aeronet code is not available. Exception: {ex}')
            if not os.path.isdir(code_aeronet):
                print(f'[ERROR] Folder with the code {code_aeronet} does not exist. You can download it by using GIT:')
                print(f'git clone https://github.com/luiscnr/aeronet.git')
            return False
        try:
            areader = AERONETReader(self.path_aeronet_nc)
            areader.dataset.close()
            return True
        except Exception as ex:
            print(f'[ERROR] AERONET-OC file {self.path_aeronet_nc} could not be started. Exception: {ex}')
            return False
