from INSITU_base import INSITUBASE
import subprocess, os
from datetime import datetime as dt
from datetime import timedelta


class INSITU_HYPERNETS_DAY(INSITUBASE):

    def __init__(self, mdb_options):
        self.mdb_options = mdb_options
        self.new_mdb = None

        self.url_base = 'hypstar@enhydra.naturalsciences.be'
        self.base_folder = '/home/hypstar/'
        self.ssh_base = 'ssh -X -Y -p 9022'
        self.ls_base = 'ls processed_data/'
        self.rsync_base = f'rsync -a -e \'ssh -p 9022\' {self.url_base}:{self.base_folder}'

        self.CHECK_SSH = self.check_ssh()

    def add_insitu(self, extract_path, ofile):
        self.start_add_insitu(extract_path, ofile)
        print('NEW MDB NO DEBERIA SER NONE', self.new_mdb)

    def get_files(self,site, sat_time):
        pathbase = self.mdb_options.insitu_path_source
        if pathbase.split('/')[-1]!=site:
            pathbase = os.path.join(pathbase,site)
        year_str = sat_time.strftime('%Y')
        month_str = sat_time.strftime('%m')
        day_str = sat_time.strftime('%d')
        path_day = os.path.join(pathbase, year_str, month_str, day_str)
        if not os.path.exists(path_day):
            return None

    def create_path_day(self,site,time):
        pathbase = self.mdb_options.insitu_path_source
        if pathbase.split('/')[-1] == site:
            path_site = pathbase
        else:
            path_site = os.path.join(pathbase,site)
            if not os.path.exists(path_site):
                os.mkdir(path_site)
        year_str = time.strftime('%Y')
        path_year = os.path.join(path_site,year_str)
        if not os.path.exists(path_year):
            os.mkdir(path_year)
        month_str = time.strftime('%m')
        path_month = os.path.join(path_year,month_str)
        if not os.path.exists(path_month):
            os.mkdir(path_month)
        day_str = time.strftime('%d')
        path_day = os.path.join(path_month, day_str)
        if not os.path.exists(path_day):
            os.mkdir(path_day)
        return path_day



    def get_files_day_ssh(self, sitename, sat_time):

        level = 'L2A'
        year_str = sat_time.strftime('%Y')
        month_str = sat_time.strftime('%m')
        day_str = sat_time.strftime('%d')
        sat_time_min = sat_time - timedelta(hours=3)
        sat_time_max = sat_time + timedelta(hours=3)

        list_files_d = {}

        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{year_str}/{month_str}/{day_str}'
        sequence_list = self.get_list_sequence_folders(cmd)

        if len(sequence_list) > 0:
            for sequence in sequence_list:
                print(f'[INFO] Date: {year_str}-{month_str}-{day_str} Checking sequence folder: {sequence}')
                insitu_time_str = sequence[3:]
                insitu_time = dt.strptime(insitu_time_str, '%Y%m%dT%H%M%S')
                if sat_time_min <= insitu_time <= sat_time_max:
                    cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{year_str}/{month_str}/{day_str}/{sequence}/*.nc'
                    list_files = self.get_list_files(cmd)
                    for file in list_files:
                        if file.find(level) > 0:
                            list_files_d[insitu_time_str] = file

        return list_files_d

    def transfer_files_ssh(self,list_files_d):
        for insitu_time_str in list_files_d:
            insitu_time = dt.strptime(insitu_time_str, '%Y%m%dT%H%M%S')
            path_day = self.create_path_day(insitu_time)
            input_path = list_files_d[insitu_time_str]
            name = input_path.split('/')[-1]
            output_path = os.path.join(path_day,name)
            print(f'[INFO] Transfering file: {input_path} to {output_path}')
            cmd = f'{self.rsync_base}{input_path} {output_path}'
            prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            out, err = prog.communicate()
            if err:
                print(err)
        print(f'Transfering files completed')

    def check_ssh(self):
        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}'
        try:
            subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT, timeout=10)
            return True
        except:
            print(f'[ERROR] Access to {self.url_base} via ssh is not allowed')
            return False

    def get_start_and_end_dates(self, sitename):
        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}'
        list_year = self.get_list_folder_dates(cmd)
        start_date = None
        end_date = None
        if len(list_year) > 0:
            for y in list_year:
                cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{y}'
                list_month = self.get_list_folder_dates(cmd)
                if len(list_month) > 0:
                    for m in list_month:
                        print(f'[INFO] Checking dates via SSH. Year: {y} Month: {m}')
                        cmd = f'{self.ssh_base} {self.url_base} {self.ls_base}{sitename}/{y}/{m}'
                        list_days = self.get_list_folder_dates(cmd)
                        if len(list_days) > 0:
                            for d in list_days:
                                datehere_str = f'{y}{m}{d}'
                                datehere = dt.strptime(datehere_str, '%Y%m%d')
                                if start_date is None:
                                    start_date = datehere
                                    end_date = datehere
                                else:
                                    if datehere < start_date:
                                        start_date = datehere
                                    if datehere > end_date:
                                        end_date = datehere

        return start_date, end_date

    def get_list_folder_dates(self, cmd):
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

    def get_list_sequence_folders(self, cmd):
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = prog.communicate()
        list = out.decode('utf-8').split('\n')
        listd = []
        for l in list:
            try:
                if l.startswith('SEQ'):
                    listd.append(l)
            except:
                pass
        return listd

    def get_list_files(self, cmd):
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = prog.communicate()
        list = out.decode('utf-8').split('\n')
        listd = []
        for l in list:
            try:
                listd.append(l)
            except:
                pass
        return listd
