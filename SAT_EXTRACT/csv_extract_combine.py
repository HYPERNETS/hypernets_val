import argparse,os
import pandas as pd

parser = argparse.ArgumentParser(description="Combine CSV files with Extract and Index columns from different satellite sources. Other options using -mode")
parser.add_argument('-m', "--mode", help="Work options. See -sh for explanations",choices=['combine','split_s3'])
parser.add_argument('-c', "--config_file", help="Config file")
parser.add_argument('-o', "--output_file", help="Output CSV file")
parser.add_argument('-i', "--input_file", help="Input CSV file")
parser.add_argument('-rf', "--ref_file", help="Reference CSV file")
parser.add_argument('-rc', "--ref_col",help="Reference (index) column in CSV file")
parser.add_argument('-sh', "--show_help", help="Show help for options",action="store_true")
args = parser.parse_args()


def combine_extract_and_index_columns():
    if not args.config_file:
        print(f'[ERROR] Argument config_file is required')
        return
    if not args.output_file:
        print(f'[ERROR] Argument output_file is required')
        return
    config_file = args.config_file
    output_file = args.output_file
    if not os.path.isfile(config_file):
        print(f'[ERROR] Config. file {config_file} does not exist or is a valid file')
        return

    csv_files = {}
    fc = open(config_file, 'r')
    for line in fc:
        sline = line.strip().split(';')
        if len(sline) == 2:
            key = sline[0].strip()
            file_csv = sline[1].strip()
            if os.path.isfile(file_csv):
                csv_files[key] = file_csv
    fc.close()
    if len(csv_files) == 0:
        print(f'[ERROR] CSV files were not identified in the config file')
        return
    if not os.path.isdir(os.path.dirname(output_file)):
        print(f'[ERROR] {output_file} could not be created as parent directory does not exist')
        return
    if not output_file.endswith('.csv'):
        print(f'[ERROR] {output_file} shoud be a CSV file')
        return
    df = None
    for key in csv_files:
        try:
            df_here = pd.read_csv(csv_files[key], sep=';')
        except:
            print(f'[ERROR] {csv_files[key]} is not a valid CSV file (separated by ;)')
            return
        c_index = f'Index_{key}'
        c_extract = f'Extract_{key}'
        if df is None:
            df = df_here.copy()
            df.rename(columns={'Extract': c_extract, 'Index': c_index}, inplace=True)
        else:
            df[c_extract] = df_here['Extract']
            df[c_index] = df_here['Index']
    if df is not None:
        df.to_csv(output_file, sep=';')
        print(f'[INFO] Combination of csv data completed in file: {output_file}')
    else:
        print(f'[ERROR] New data frame could not be created')

def split_s3():
    if not args.input_file:
        print(f'[ERROR] Argument input file is required')
        return
    if not args.ref_file:
        print(f'[ERROR] Argument ref_file is required')
        return
    if not args.ref_col:
        print(f'[ERROR] Argument ref_col is required')
        return
    input_file = args.input_file
    ref_file = args.ref_file
    ref_col = args.ref_col
    if not os.path.isfile(input_file):
        print(f'[ERROR] {input_file} is not a valid input file')
        return
    if not os.path.isfile(ref_file):
        print(f'[ERROR] {input_file} is not a valid ref file')
        return

    try:
        df_input = pd.read_csv(input_file,sep=';')
    except:
        print(f'[ERROR] Input file {input_file} is not a valid CSV file separated by ;')
        return
    try:
        df_ref = pd.read_csv(ref_file,sep=';')
    except:
        print(f'[ERROR] Reference file {ref_file} is not a valid CSV file separated by ;')
        return
    if ref_col not in df_input.columns.tolist():
        print(f'[ERROR] Reference column {ref_col} is not available in the input csv file')
        return
    if ref_col not in df_ref.columns.tolist():
        print(f'[ERROR] Reference column {ref_col} is not available in the reference csv file')
        return
    if 'Extract' not in df_input.columns.tolist():
        print(f'[ERROR] Extract column  is not available in the input csv file')
        return
    if 'Index' not in df_input.columns.tolist():
        print(f'[ERROR] Index column  is not available in the input csv file')
        return

    print(f'[INFO] Started...')
    n_values = len(df_ref[ref_col])
    print(f'[INFO] Number of stations: {n_values}')
    extracts_s3a = ['NaN']*n_values
    index_s3a = [-1] * n_values
    extracts_s3b = ['NaN'] * n_values
    index_s3b = [-1] * n_values
    for index,row in df_input.iterrows():
        ref_here = row[ref_col]
        try:
            platform = str(row['Extract'])[:3]
            index_here = df_ref[df_ref[ref_col]==ref_here].index.tolist()[0]
            if platform=='S3A':
                extracts_s3a[index_here] = row['Extract']
                index_s3a[index_here] = row['Index']
            if platform=='S3B':
                extracts_s3b[index_here] = row['Extract']
                index_s3b[index_here] = row['Index']
        except:
            continue

    df_output_s3a = df_ref.copy()
    df_output_s3a['Extract'] = extracts_s3a
    df_output_s3a['Index'] =index_s3a
    file_output_s3a = os.path.join(os.path.dirname(input_file),os.path.basename(input_file)[:-4]+'_S3A.csv')
    print(f'[INFO] Saving Sentine-3A file: {file_output_s3a}')
    df_output_s3a.to_csv(file_output_s3a,sep=';',index=False)
    df_output_s3b = df_ref.copy()
    df_output_s3b['Extract'] = extracts_s3b
    df_output_s3b['Index'] = index_s3b
    file_output_s3b = os.path.join(os.path.dirname(input_file), os.path.basename(input_file)[:-4] + '_S3B.csv')
    print(f'[INFO] Saving Sentine-3B file: {file_output_s3b}')
    df_output_s3b.to_csv(file_output_s3b, sep=';',index=False)




    return
def main():
    if not args.mode:
        combine_extract_and_index_columns()
        return
    if args.mode:
        if args.mode=='combine':
            combine_extract_and_index_columns()
            return
        if args.mode=='split_s3':
            split_s3()
            return






if __name__ == '__main__':
    if args.mode  and args.show_help:
        if args.mode=='combine':
            print('Combine CSV files with Extract and Index columns from different satellite sources')
            print('Required:')
            print('    -c/--config_file: list of key;csv_file lines')
            print('    -o/--output_file: Output CSV file')
        if args.mode=='split_s3':
            print('Split csv file with EUMETSAT WFR extracts from Sentinel-3A and 3B in two csv Extract file for each mission')
            print('Required:')
            print('    -i/--input_file: Input CSV file')
            print('    -rf/--ref_file: Reference CSV file with all the potential stations')
            print('    -i/--ref_col: Reference column to the stations')
    else:
        main()