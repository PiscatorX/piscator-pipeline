#!/usr/bin/python 

from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys

def get_args():

    parser = argparse.ArgumentParser(description="compile amplicon csv  data")
    parser.add_argument('-p','--primer', dest='primer_input',action='store', required=True, type=str,
                        help="primer file with id<TAB>seq")
    parser.add_argument('-f','--format', dest='seq_format',default="fasta", action='store_true', required=False)

    return parser.parse_args()


def write_fasta(primer_input, seq_format):

    primers = map(mk_fasta, open(primer_input).read().splitlines())
    for rec in primers:
        SeqIO.write(rec,'.'.join([rec.id, seq_format]), seq_format)

def mk_fasta(line):

    seq_out = lambda x:SeqRecord(Seq(x[1], IUPAC.ambiguous_dna), id=x[0], description='')    
    return seq_out(line.split('\t'))

args = get_args()
write_fasta(args.primer_input,args.seq_format)
