from INSITU_base import INSITUBASE
import subprocess
from datetime import datetime as dt


class INSITU_HYPERNETS_DAY(INSITUBASE):

    def __init__(self, mdb_options):
        self.mdb_options = mdb_options
        self.new_mdb = None

        self.url_base = 'hypstar@enhydra.naturalsciences.be'
        self.ssh_base = 'ssh -X -Y -p 9022'
        self.ls_base = 'ls processed_data/'

        self.CHECK_SSH = self.check_ssh()

    def add_insitu(self,extract_path,ofile):
        self.start_add_insitu(extract_path,ofile)
        print('NEW MDB NO DEBERIA SER NONE',self.new_mdb)

    def check_ssh(self):
        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}'
        try:
            subprocess.check_output(cmd, shell = True,stderr=subprocess.STDOUT, timeout=10)
            return True
        except:
            print(f'[ERROR] Access to {self.url_base} via ssh is not allowed')
            return False


    def get_start_and_end_dates(self,sitename):
        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}'
        list_year = self.get_list_folder_dates(cmd)
        start_date = None
        end_date = None
        if len(list_year)>0:
            for y in list_year:
                cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{y}'
                list_month = self.get_list_folder_dates(cmd)
                if len(list_month)>0:
                    for m in list_month:
                        print(f'[INFO] Checking dates via SSH. Year: {y} Month: {m}')
                        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{y}/{m}'
                        list_days = self.get_list_folder_dates(cmd)
                        if len(list_days)>0:
                            for d in list_days:
                                datehere_str = f'{y}{m}{d}'
                                datehere = dt.strptime(datehere_str,'%Y%m%d')
                                if start_date is None:
                                    start_date = datehere
                                    end_date = datehere
                                else:
                                    if datehere<start_date:
                                        start_date = datehere
                                    if datehere>end_date:
                                        end_date = datehere

        return start_date,end_date

    def get_list_folder_dates(self,cmd):
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = prog.communicate()
        list = out.decode('utf-8').split('\n')
        listd = []
        for l in list:
            try:
                int(l)
                listd.append(l)
            except:
                pass
        return listd
