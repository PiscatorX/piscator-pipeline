#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys
import os

def  parse_data(filename, select, outfname):

    select_list = open(args.select).read().splitlines()
    seq_data  = SeqIO.parse(open(filename), "fasta")
    select_seq_data =  []
    for rec in seq_data:
        if rec.id in  select_list:
            sys.stderr.write(rec.id+"\n")
            select_seq_data.append(rec)
            
    if  outfname:
        path, fname = os.path.split(outfname)
        new_filename = os.path.join(path, fname)
        if os.path.exists(new_filename):
            raise Exception("File exists, try another file name")
    else:
        print("~~~")
        path, fname = os.path.split(filename)
        new_filename = 'extract_' + fname
         
    with open(new_filename, "w") as new_fobj:   
        SeqIO.write(select_seq_data, new_fobj, 'fasta')
        


if __name__ == '__main__':
    parser = argparse.ArgumentParser("get sequences with id similar to provided file")
    parser.add_argument('reference', help = 'fasta reference file')
    parser.add_argument('-s','--select', help = 'list of words to look for in fasta id info', required=True )
    parser.add_argument('-o','--outfname', help = 'list of words to look for in fasta id info')
    args = parser.parse_args()
    parse_data(args.reference, args.select, args.outfname)
