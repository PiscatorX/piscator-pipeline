#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO
import argparse




def parse_ref(ref_fasta,  outfile, depth=6):

    sequences = SeqIO.parse(ref_fasta, "fasta")
    ref_file_fp  = open(outfile, "a")
    map_file_fp  = open(outfile+'.map', "a")
    seq_ids = {}
    for seq in sequences:
        if seq.id in seq_ids:
            print("**skipping duplicate {} ==> {}".format(seq.id, seq_ids[seq.id]))
            continue
        seq_ids[seq.id] = len(seq)
        taxonomy = seq.description.split(" ", 1)[1:]
        taxa_ranks = taxonomy[0].split(";")
        ranks = len(taxa_ranks)
        if  ranks < 6:
            continue
        map_ref = ';'.join(taxa_ranks[:depth])
        seq.description = ''
        SeqIO.write(seq, ref_file_fp, "fasta")
        print("{}\t{}".format(seq.id, map_ref), file=map_file_fp)
        
    ref_file_fp.close()
    map_file_fp.close()



if  __name__ == '__main__':
     parser = argparse.ArgumentParser(description="""trim fast reference database""")
     parser.add_argument('ref_fasta')
     parser.add_argument('-o','--outfile',  required = True)
     args, unknown = parser.parse_known_args()
     parse_ref(args.ref_fasta, args.outfile)
