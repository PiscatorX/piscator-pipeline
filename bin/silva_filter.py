#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys




def  parse_data(reference, select, outfname):
    global select_list
    select_list =  [ line.strip().lower() for line in select.read().splitlines() ]
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
    parser.add_argument('-s','--select', type=argparse.FileType('r'), help = 'list ranks to look for', required=True )
    parser.add_argument('-o','--outfname', type=argparse.FileType('w'), help = 'list of words to look for in silva taxonomic', required = True)
    args = parser.parse_args()
    parse_data(args.reference, args.select, args.outfname)
