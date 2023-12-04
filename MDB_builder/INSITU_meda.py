from datetime import datetime as dt
from datetime import timedelta
import os

class INSITU_MEDA:
    def __init__(self):
        self.path_meda = '/store2/data/meda/binary'
        self.file_name = 'meda_lam_opt_$DATE$_L1v2.nc'
        self.format_date = '%y%m%d'

    def check_file(self,date):

        path = os.path.join(self.path_meda,date.strftime('%Y'),date.strftime('%m'),date.strftime('%d'))
        if not os.path.isdir(path):
            return None
        date_name_file = date.strftime(self.format_date)
        name_file = self.file_name.replace('$DATE$',date_name_file)
        meda_file = os.path.join(path,name_file)
        if os.path.isfile(meda_file):
            return meda_file
        else:
            return None

    def get_datelist_file(self,output_file,start_date,end_date):
        f1 = open(output_file,'w')
        date_ref = start_date
        while date_ref<=end_date:
            meda_file = self.check_file(date_ref)
            if meda_file is not None:
                date_ref_str = date_ref.strftime('%Y-%m-%d')
                f1.write(date_ref_str)
                f1.write('\n')
            date_ref = date_ref + timedelta(hours=24)
        f1.close()
