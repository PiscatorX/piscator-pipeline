#!/usr/bin/env python
from init_primerDB import PrimerDB
from itertools import product
import mysql.connector
import numpy as np
from Bio import Seq
import  collections
import  subprocess
import  argparse
import  sqlite3
import  pprint
import  shlex
import  sys




class PhysProp(PrimerDB):

    def __init__ (self):

        super(PhysProp, self).__init__()
        
        parser = argparse.ArgumentParser(description="""Save primer data""")
        parser.add_argument('primers_file')
        parser.add_argument('-w','--windowsize',dest='windowsize',action='store',required=False, type=int)
        parser.add_argument('-s','--shiftincrement', dest='shiftincrement', action='store',required=False,type=int)
        parser.add_argument('-D','--dnaconc', dest='dnaconc', action='store',required=False,default=50,type=int)
        parser.add_argument('-S','--salt', dest='saltconc', action='store', required=False, default=50,type=int)
        parser.add_argument('-t','--temparature', dest='temperature', action='store',default=55,required=False, type=int)

        
        args = parser.parse_args()
        
        self.cnx.database = self.DB_NAME  
    
        self.primer_id, self.primer = open(args.primers_file).readline().split()
        
        self.args = vars(args)
        n = len(self.primer)
        
        self.args.update({'windowsize': n, 'shiftincrement': n})
        self.emboss_dan_args = """dan -filter -windowsize {windowsize} -shiftincrement {shiftincrement} -dnaconc {dnaconc}
          -saltconc {saltconc} -product -thermo -temperature {temperature} -outfile stdout -rformat simple""".format(**self.args)        

        
    def run_dan(self, primer):

        args  = shlex.split(self.emboss_dan_args)

        proc = subprocess.Popen(args, stdin= subprocess.PIPE,
               stdout = subprocess.PIPE, stderr= subprocess.PIPE)
        stdout, stderr = proc.communicate(input=primer)
            
        #print stdout
        
        d = Seq.IUPAC.IUPACData.ambiguous_dna_values

        primer_combinations = []

        for i in product(*[d[j] for j in primer]):
            primer_combinations.append("".join(i))

        #pprint.pprint(primer_combinations)
        
        self.tm_keys = ['Tm','TmProd','GC','DeltaG','DeltaH','DeltaS']
        
        Tm_dat = dict( (col, []) for col in self.tm_keys )
        n = len(primer_combinations)
        for seq in primer_combinations:
            
            proc = subprocess.Popen(args, stdin= subprocess.PIPE,
               stdout = subprocess.PIPE, stderr= subprocess.PIPE)
            stdout, stderr = proc.communicate(input=seq)

            for param in self.tm_keys:
                field_tmp =  dict([ line.split(':') for line in stdout.split('\n') if line.startswith(param+':') ])
                Tm_dat[param].append(field_tmp[param])
            
        return stdout, Tm_dat

    
    def populate_results(self, primer_id, primer_results):
        
        primer_results.update({'primer_ID':primer_id})
        sql = """INSERT INTO primerprop
                  (primer, TmProd, DeltaS, Length, GC, DeltaG, DeltaH, Tm)
                   values("{primer_ID}", "{TmProd}", "{DeltaS}", "{Length}", "{GC}", "{DeltaG}", "{DeltaH}", "{Tm}")
               """.format(**primer_results)
        try:
            self.cursor.execute(sql)
        except mysql.connector.errors.IntegrityError as err:
            print(err)
            
        self.cnx.commit()

    def get_physprop(self):

        stdout, Tm_dat =  self.run_dan(self.primer)
        #print stdout
        if stdout:
            primer_results = dict( tuple(rec.split(': ')) for rec in 
                [ line for line in stdout.split('\n')\
                if not line.startswith('#') ] if rec )
            primer_results['Length'] = str([primer_results['Length']])
            
            for k in self.tm_keys:
                primer_results[k] = str(Tm_dat[k])
                
            self.populate_results(self.primer_id, primer_results)

if __name__ == '__main__':            
    PhysProp = PhysProp()
    PhysProp.get_physprop()
