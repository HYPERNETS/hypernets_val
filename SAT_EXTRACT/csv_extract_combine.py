import argparse,os
import pandas as pd

parser = argparse.ArgumentParser(description="Combining CSV files with Extract and Index columns from different satellite sources.")
parser.add_argument('-c', "--config_file", help="Config file, list of key;csv_file lines",required=True)
parser.add_argument('-o', "--output_file", help="Output CSV file",required=True)
args = parser.parse_args()

def main():
    config_file = args.config_file
    output_file = args.output_file
    if not os.path.isfile(config_file):
        print(f'[ERROR] Config. file {config_file} does not exist or is a valid file')
        return

    csv_files = {}
    fc = open(config_file,'r')
    for line in fc:
        sline = line.strip().split(';')
        if len(sline)==2:
            key = sline[0].strip()
            file_csv = sline[1].strip()
            if os.path.isfile(file_csv):
                csv_files[key] = file_csv
    fc.close()
    if len(csv_files)==0:
        print(f'[ERROR] CSV files were not identified in the config file')
        return
    if not os.path.isdir(os.path.dirname(output_file)):
        print(f'[ERROR] {output_file} could not be created as parent directory does not exist')
        return
    if output_file.endswith('.csv'):
        print(f'[ERROR] {output_file} shoud be a CSV file')
        return
    df = None
    for key in csv_files:
        try:
            df_here = pd.read_csv(csv_files[key],sep=';')
        except:
            print(f'[ERROR] {csv_files[key]} is not a valid CSV file (separated by ;)')
            return





if __name__ == '__main__':
    main()