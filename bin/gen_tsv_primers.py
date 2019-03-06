#!/usr/bin/env python
from init_primerDB import PrimerDB
import mysql.connector
import argparse
import sys
import os


class generateTSV(PrimerDB):

    def __init__(self):

        super(generateTSV, self).__init__()

        parser = argparse.ArgumentParser(description="Generate tsv files for the primers",
        epilog='NOTE: This utility expects database to be have been populated')
        parser.add_argument('-o','--output-dir', dest='output_dir')
        
        args, unknown = parser.parse_known_args()
        
        self.output_dir = args.output_dir if args.output_dir else '.'

        self.cnx.database = self.DB_NAME  

        

    def get_pairs(self, ext = '.tsv'):
        
        sql = """SELECT Fwd_id, Fwd_Primer, Rev_id,
                 Rev_Primer from primers;"""
        
        self.cursor.execute(sql)
        
        self. primer_pairs  =   ( (fwd_id, fwd, rev_id, rev) for fwd_id, fwd, rev_id, rev\
                in self.cursor.fetchall())


    def generateTSV(self):
        
        ext = '.tsv'  
        fix_id = lambda primer_id, Dir: primer_id if primer_id.endswith(Dir) else primer_id+Dir
        save   = lambda primer_id, primer: open(os.path.join(self.output_dir, primer_id) + ext,'w').write(primer)
        out = lambda primer_id:  primer_id + ext
        
        output = []
        #pair_data = []
        
        for fwd_id, fwd, rev_id, rev in self.primer_pairs:
            fwd_id, rev_id  = map(lambda x:\
                              x.replace('_','.'),  (fix_id(fwd_id.lower(),'f'), fix_id(rev_id.lower(),'r')))
            fwd =  '\t'.join([fwd_id,fwd])
            rev =   '\t'.join([rev_id,rev])
            for primer_id, direct in [(fwd_id, fwd), (rev_id, rev)]:
                save(primer_id, direct)
            to_nextflow = map(out, [fwd_id, rev_id])
            #pair_data.append([fwd_id, rev_id])
            output.append(','.join(to_nextflow)+'\n')
        sys.stdout.writelines(output)
    

if __name__ == '__main__':
    get_tsv = generateTSV()
    get_tsv.get_pairs()
    get_tsv.generateTSV()


    
