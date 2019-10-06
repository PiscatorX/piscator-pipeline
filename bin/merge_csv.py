#!/usr/bin/env python

import pandas as pd
import argparse
import glob
import sys
import os

def merge_csv(data_directory, ext_glob, outfname):

    csv_files  = glob.glob(os.path.join(data_directory, ext_glob))
    csv_dfs  = []
    for  csv in csv_files:
        df = pd.read_csv(csv, index_col=None, header=0)
        csv_dfs.append(df)
    cat_df  = pd.concat(csv_dfs)
    cat_df.to_csv(outfname, sep = "\t")

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="merge csv files")
    parser.add_argument('data_directory', action='store')
    parser.add_argument('-o','--outfname',  required=True)
    parser.add_argument('-e','--ext_glob', default = "*csv")
    args = parser.parse_args()
    clip = merge_csv(args.data_directory, args.ext_glob, args.outfname)
    
