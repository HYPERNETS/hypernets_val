

class SBATCH_SCRIPTER(object):

    def __init__(self,output_file):
        self.output_file = output_file
        self.fw = None

    def close_script(self):
        if self.fw is not None:
            self.fw.close()

    def start_script(self,sb_options,add_conda):
        if self.fw is not None:
            self.close_script()
        self.fw = open(self.output_file,'w')
        first_line = '#!/bin/bash'
        self.fw.write(first_line)
        self.add_line('#SBATCH --nodes=1')
        self.add_line('#SBATCH --ntasks=1')
        self.add_line(f'#SBATCH -p {sb_options["sbatch_partition"]}')
        if 'sbatch_email' in sb_options:
            self.add_line(f'#SBATCH --mail-user {sb_options["sbatch_email"]}')
            if 'sbatch_email_type' in sb_options:
                self.add_line(f'#SBATCH --mail-type {sb_options["sbatch_email_type"]}')

        if add_conda:
            self.add_blank_lines(3)
            self.add_conda_environment(sb_options)
            self.add_blank_lines(3)

    def add_conda_environment(self,sb_options):
        if self.fw is None:
            return
        if 'conda_source' in sb_options.keys() and 'conda_env' in sb_options.keys():
            self.add_line(f'source /home/$USER/{sb_options["conda_source"]}')
            self.add_line(f'conda activate {sb_options["conda_env"]}')

    def add_line(self,line):
        if self.fw is not None:
            self.fw.write('\n')
            self.fw.write(line)

    def add_blank_lines(self,num):
        if self.fw is not None:
            for idx in range(num):
                self.add_line('')



class SH_SCRIPTER():
    def __init__(self,output_file):
        self.output_file = output_file
        self.fw = None

    def close_script(self):
        if self.fw is not None:
            self.fw.close()

    def start_script(self):
        if self.fw is not None:
            self.close_script()
        self.fw = open(self.output_file,'w')
        first_line = '#!/bin/bash'
        self.fw.write(first_line)

    def add_line(self,line):
        if self.fw is not None:
            self.fw.write('\n')
            self.fw.write(line)

    def add_blank_lines(self,num):
        if self.fw is not None:
            for idx in range(num):
                self.add_line('')

    def add_sbatchs(self,script_list,log_list,nmax_cores):

        use_log = True if log_list is not None else False
        for ifile in range(len(script_list)):
            iprev = ifile-nmax_cores
            if iprev<0:
                line = f'job{ifile}=$(sbatch --output={log_list[ifile]} {script_list[ifile]})' if use_log else f'job{ifile}=$(sbatch {script_list[ifile]})'
                self.add_line(line)
                line = f'jobid{ifile}=$(echo "$job{ifile}" | awk \'{{print $NF}}\')'
                self.add_line(line)
            else:
                line = f'job{ifile}=$(sbatch --dependency=afterany:$jobid{iprev} --output={log_list[ifile]} {script_list[ifile]})' if use_log else f'job{ifile}=$(sbatch --dependency=afterany:$jobid{iprev} {script_list[ifile]})'
                self.add_line(line)
                line = f'jobid{ifile}=$(echo "$job{ifile}" | awk \'{{print $NF}}\')'
                self.add_line(line)

def prepare_sh_script_with_multiple_sbatch(file_sh,script_list,log_list,nmax_cores):
    shs = SH_SCRIPTER(file_sh)
    shs.start_script()
    shs.add_blank_lines(2)
    shs.add_sbatchs(script_list,log_list,nmax_cores)
    shs.add_blank_lines(2)
    shs.close_script()
