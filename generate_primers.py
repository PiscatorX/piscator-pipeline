#!/bin/python
import itertools as it
import  tempfile
import  sqlite3
import  subprocess
import  shlex
import  shutil
import  sys
import  pprint
import  argparse
import  os


class PiscatorDB_Primers(object):
    
    def  __init__(self):
        
        """
           Connect to the PiscatorDB and generate primer pairs

        """
        
        parser = argparse.ArgumentParser(description="compile amplicon csv  data")
        parser.add_argument('-d','--db', dest='db_name',
                            action='store', required=True, type=str)
        parser.add_argument('-o','--output', dest='output_dir',
                            action='store', required=True, type=str)
        args = parser.parse_args()

        self.primerDB_conn  = sqlite3.connect(args.db_name)
        self.primerDB_c = self.primerDB_conn.cursor()
        self.output_dir = args.output_dir

        
    def get_primer_tsv(self):
        
        sql = """SELECT Fwd_id, Forward_Primer, Rev_id, 
                 Rev_Primer from primers;"""
        
        return  ((fwd_id, fwd, rev_id, rev)\
                 for fwd_id, fwd, rev_id, rev\
                 in self.get_results(sql))

            
    def get_primers_pairs(self):
        
        sql = """SELECT Fwd_id, Forward_Primer, Rev_id,
                 Rev_Primer from primers;"""
        
        return  ((fwd_id, fwd, rev_id, rev)\
                 for fwd_id, fwd, rev_id, rev\
                 in self.get_results(sql))
        
    def get_results(self, sql):
       
        """
           Return the results of an SQL query 

        """
        
        self.primerDB_c.execute(sql)

        return iter(self.primerDB_c.fetchall())
    

    def custom_sql(self, sql):
        
        """
           return results of custom sql   

        """

        return self.get_results(sql) 
        
    
    
    def Primer_pair_amplicons(self):
         
        """
           Generate a csv of primer pairs            

        """
        
        primer_pairs = self.get_primers_pairs()
        
        for fwd_id, fwd, rev_id, rev in primer_pairs:
            
            primer_id = '~'.join([fwd_id, rev_id])
            primer_data  = '\t'.join([primer_id, fwd, rev])

            with open(primer_id+'.tsv','w') as fp:
                fp.write(primer_data)

if __name__ == '__main__':
    Piscator_primers = PiscatorDB_Primers()
    Piscator_primers.Primer_pair_amplicons() 
    
