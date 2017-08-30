#!/usr/bin/python
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import argparse
import pprint
from primer_map import PrimerDB
import os

class PrimerBlastn(object):

    def __init__(self):
        
        parser = argparse.ArgumentParser(description="""get primers positions on subject sequence""")
        parser.add_argument('-q','--query', dest='query', action='store', required=True, type=str,
                            help="primer file in fasta format")
        parser.add_argument('-s','--subject', dest='subject', action='store', required=True, type=str,
                            help="subject sequence file in fasta format")
        parser.add_argument('-f','--fmt', dest='outfmt', action='store',default=5,required=False,type=int,
                            help="output file format")
        parser.add_argument('-o','--out', dest='out', action='store',required=False, type=str,
                            help="output filename (optional)")
        parser.add_argument('-i','--identity', dest='perc_identity', action='store',default=50, required=False, type=int)
        parser.add_argument('-w','--word_size', dest='word_size', action='store',default=5, required=False, type=int)
        parser.add_argument('-a','--add',dest='add',nargs='*',  
                            help="Additional unvalidated blastn args in the form keyword=value e.g \"evalue=0.001\"")
        parser.add_argument('-d','--db-name', dest='db_name', action='store', required=True,  type=str,
                            help = "Database where primer data will saved.")        

        self.parser_args = parser.parse_args()
        self.args = vars(self.parser_args)
        self.db_name = self.args['db_name']
        self.primer_name = os.path.basename(self.args['query']).split('.')[0]
        

    def  Blast_it(self):
       
         if not self.args['out']:
             xml_out  = lambda f: '.'.join([f.split('.')[0],'xml'])                            
             self.args['out'] = xml_out(self.args['query'])

         self.blast_out = self.args['out']

         if self.args['add']:
            self.args.update(dict( (var.split('=')) for var in self.args['add']))             
         del self.args['add']
         del self.args['db_name']

         Blast_cline = NcbiblastnCommandline(**self.args) 
         #print Blast_cline
         Blast_cline()


        
    def parse_xml(self):

        blast_results = NCBIXML.parse(open(self.blast_out))
        for results in blast_results:
            for alignment in results.alignments:
                for rec in alignment.hsps:
                     self.primer_data = vars(rec)
                     self.primer_data['primer_ID'] = self.primer_name
                     return 
                     
if  __name__  ==  '__main__':
    primer_blast = PrimerBlastn()
    primer_map = PrimerDB(primer_blast.args['db_name'])
    primer_blast.Blast_it()
    primer_blast.parse_xml()
    pprint.pprint(primer_blast.primer_data)
    primer_map.DB_Insert(primer_blast.primer_data)
