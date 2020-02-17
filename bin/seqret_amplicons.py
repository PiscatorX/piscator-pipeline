#!/usr/bin/env python

from init_primerDB import PrimerDB
from Bio import SeqIO
import argparse
import pprint
import sys


class GetAmplicons(PrimerDB):

    def __init__(self):

        super(GetAmplicons, self).__init__()
        
        parser = argparse.ArgumentParser(description="""Save primer data populate mysql database""")
        parser.add_argument('-R','--reference', help="fasta reference file ", type = str, required = True)
        parser.add_argument('-f','--fwd_mis', help="maximum mismatches on forward primer", type = int, default = 3)
        parser.add_argument('-r','--rev_mis', help="maximum mismatches on reverse primer", type = int, default = 3)
        args, unknown = parser.parse_known_args()
        self.ref_seq = args.reference
        self.cnx.database = self.DB_NAME  
        self.rev_mis = args.rev_mis
        self.fwd_mis = args.fwd_mis
        
    def parse_ref_seq(self):
        
        records = SeqIO.parse(open(self.ref_seq),'fasta')
        self.ref_seq_dict  = dict((rec.id, rec) for rec in  records)

    def get_amplicons(self):
        
        col_ids =[ 'primer_id','amplimer','fwd','fwd_mis','rev','rev_mis','len','seq_id']
        self.cursor.execute("""SELECT {} FROM amplicons 
                               WHERE rev_mis <= {} AND rev_mis <= {}  
                               ORDER BY primer_id ASC""".format(",".join(col_ids),self.rev_mis, self.fwd_mis))
        amplicon_data = self.cursor.fetchall()
        
        fobj = False
        curr = None
        last = None
        seq_data  = None

        amplicon_records = []
        
        for row in amplicon_data:
            #print len(amplicon_records),curr,last
            rec = dict(zip(col_ids, row))
            #print rec
            try:
                seq_data = self.ref_seq_dict.get(rec['seq_id'],False)
            except KeyError as e:
                 print row
                 raise Exception, e 
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

            fwd = int(rec['fwd'])  - 1
            rev = init_len - int(rec['rev']) + 1 
            amplicon = seq_data[fwd:rev]
            amplicon_records.append(amplicon)
            last =  curr
            
        SeqIO.write(amplicon_records, open(last+'.fasta','w'),'fasta')

            
if __name__ ==  '__main__':
    amplicons = GetAmplicons()
    amplicons.parse_ref_seq()
    amplicons.get_amplicons()
    

                
