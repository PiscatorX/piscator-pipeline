#!/usr/bin/python
from Bio import SeqIO
import sys
import pprint

def parse_ref_seq(ref_seq):

    records = SeqIO.parse(open(ref_seq),'fasta')
    return dict((rec.id, rec) for rec in  records)


def parse_results(ref_seq_dict, search_results):

    col_id =('primer_id','amplimer','fw','fw_mis',\
           'rev','rev_mis','len','seq_id')
    
    seq_data  = None
    with open(search_results)  as fp:

        fobj = False
        curr = None
        last = None
        amplicon_records = []
        for line in  fp:
            print len(amplicon_records),curr,last
            rec = dict(zip(col_id, line.strip().split('\t')))
            try:
                seq_data = ref_seq_dict.get(rec['seq_id'],False)
            except KeyError:
                 print line
            if not seq_data:
                print 'Not found!'
                continue
            init_len = len(seq_data.seq)
            curr = rec['primer_id']
            if curr != last:
                if last:
                    print len(amplicon_records)
                    SeqIO.write(amplicon_records, open(last+'.fasta','w'),'fasta')
                amplicon_records = []
                
            fwd = int(rec['fw'])  - 1
            rev = init_len - int(rec['rev']) + 1 
            amplicon = seq_data[fwd:rev]
            amplicon_records.append(amplicon)
            last =  curr 
        SeqIO.write(amplicon_records, open(last+'.fasta','w'),'fasta')

ref_seq_dict  = parse_ref_seq(sys.argv[1])
parse_results(ref_seq_dict, sys.argv[2])

                
