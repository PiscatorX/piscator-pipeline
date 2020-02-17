#!/usr/bin/env python
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import argparse
import pprint
import pandas as pd
import sys
import os




class PrimerBlastn(object):

    def __init__(self):
        
        parser = argparse.ArgumentParser(description="""get primers positions on subject sequence""")
        parser.add_argument('-q','--query', dest='query', action='store', required=True, type=str,
                            help="primer file in fasta format")
        parser.add_argument('-s','--subject', dest='subject', action='store', required=True, type=str,
                            help="subject sequence file in fasta format")
        parser.add_argument('-m','--max_hsps', type=int, default = 1,
                            help="Set maximum number of HSPs per subject sequence to save for each query")
        parser.add_argument('-f','--fmt', dest='outfmt', action='store',default=5,required=False,type=int,
                            help="output file format")
        parser.add_argument('-o','--out', dest='out', action='store',required=False, type=str,
                            help="output filename (optional)")
        parser.add_argument('-i','--identity', dest='perc_identity', action='store',default=30, required=False, type=int)
        parser.add_argument('-w','--word_size', dest='word_size', action='store',default=5, required=False, type=int)
        parser.add_argument('-a','--add',dest='add',nargs='*',  
                            help="Additional unvalidated blastn args in the form keyword=value e.g \"evalue=0.001\"")

        
        self.parser_args = parser.parse_args()
        self.args = vars(self.parser_args)
        self.primer_name = os.path.basename(self.args['query']).split('.')[0]
        

    def  Blast_it(self):

         if not self.args['out']:
             xml_out  = lambda f: '.'.join([os.path.basename(f).split('.')[0],'xml'])                            
             self.args['out'] = xml_out(self.args['query'])

         self.blast_out = self.args['out']
         if self.args['add']:
            self.args.update(dict( (var.split('=')) for var in self.args['add']))             
         del self.args['add']
         
         #print(self.args)
         Blast_cline = NcbiblastnCommandline(**self.args) 
         #print Blast_cline
         Blast_cline()


        
    def parse_xml(self):

        blast_results = NCBIXML.parse(open(self.blast_out))
        for results in blast_results:
            for alignment in results.alignments:
                for rec in alignment.hsps:
                     self.primer_data = vars(rec)
                     self.primer_data.update({ "primer": self.primer_name })
                     print(self.primer_data)
                     primer_df = pd.DataFrame.from_dict(self.primer_data)
                     csv_fname = '.'.join([self.primer_name, 'csv'])
                     primer_df.to_csv(csv_fname, index = False)
                     

                 
                 
if  __name__  ==  '__main__':
    primer_blast = PrimerBlastn()
    primer_blast.Blast_it()
    primer_blast.parse_xml()
