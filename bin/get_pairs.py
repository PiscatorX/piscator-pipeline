#!/usr/bin/env python
import sys
import os

def get_pair(primer_pairs):

    get_data = lambda primer_tsv: open(primer_tsv).readline().split('\t')
    primer_data  = map(get_data, primer_pairs)
    seq_ids = []
    seqs = []
    [ (seq_ids.append(d[0]), seqs.append(d[1])) for d in primer_data]
   
    fname = '_'.join(seq_ids)
    with open(fname+'.tsv','w') as fp:
        fp.writelines(['\t'.join([fname,'\t'.join(seqs)])])
    
    #print os.getcwd()

get_pair(sys.argv[1:])

