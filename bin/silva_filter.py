#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys




def  parse_data(reference, select, names,  outfname):
    global select_list
    
    if select:
        select_list =  select.read().lower().split()
    elif names:
        select_list = [ entry.lower() for entry in names ]
        
    seq_data  = SeqIO.parse(reference, "fasta")
    select_seq_data =  []
    for rec in seq_data:
        taxa_data = rec.description.lower().split(';')
        if is_select(taxa_data):
            select_seq_data.append(rec)
    SeqIO.write(select_seq_data, outfname,'fasta')
    outfname.close()


def is_select(taxa_data):

    for rank in select_list:
        if rank in taxa_data:
            return True
        
    return False



if __name__ == '__main__':
    parser = argparse.ArgumentParser("get select silva sequences based on taxonomy provided")
    parser.add_argument('reference', type=argparse.FileType('r'),  help = 'fasta reference file')
    parser.add_argument('-s','--select', type=argparse.FileType('r'), help = 'file containing list of taxonomic names/words on each to filter for')
    parser.add_argument('-n','--names', nargs='+', help = 'space separated taxonomic ranks')
    parser.add_argument('-o','--outfname', type=argparse.FileType('w'), help = 'list of words/taxonomic names to look for in silva taxonomic', required = True)
    args = parser.parse_args()
    if not (args.select or args.names):
        parser.error('At least one type of input must be provied, chose either a file for  --select or  command line arguments to  --names. Try --help for more.')
    parse_data(args.reference, args.select, args.names,  args.outfname)
