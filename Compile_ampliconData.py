#!/usr/bin/python
from matplotlib import pyplot as plt
from matplotlib_venn import venn3_unweighted, venn2_unweighted
import matplotlib.cm as cm
from Bio import SeqIO
import numpy as np
import subprocess
import argparse
import sqlite3
import pprint
import os


class load_SQL_csv(object):
    
    def __init__(self):

        parser = argparse.ArgumentParser(description="compile amplicon csv  data")
        
        parser.add_argument('-t','--tsv', dest='tsv_file', 
                            action='store', required=True, type=str)
        parser.add_argument('-d','--db', dest='db_name',
                            action='store', required=True, type=str)
        parser.add_argument('-o','--output', dest='output', 
                        action='store', required=True, type=str)
        parser.add_argument('-r','--ref', dest='ref_file', 
                        action='store', required=True, type=str)
        args = parser.parse_args()
        self.db_name = args.db_name
        self.tsv_file = args.tsv_file
        self.amplicons_loaded = True
        self.output_dir = args.output
        self.ref_len = len(list(SeqIO.parse(args.ref_file,'fasta')))


        
    def import_tsv(self):

        sql_cmd = """.echo ON\n.separator "\t"\n.import {} amplicons\nselect * from amplicons;""".format(self.tsv_file)

        args = 'sqlite3 {}'.format(self.db_name)
        
        proc = subprocess.Popen(args, stdin= subprocess.PIPE,
               stdout = subprocess.PIPE, stderr= subprocess.PIPE, shell = True)

        stdout, stderr = proc.communicate(input=sql_cmd)             

        print stdout
        
        if stderr:
            raise Exception,"Sqlite3 Insert error\n{}".format(stderr)
        
        else:
            self.amplicons_loaded = True
         
    def AmpliconPlots(self):

        primerDB_conn  = sqlite3.connect(self.db_name)
        primerDB_c = primerDB_conn.cursor()
        tup2str = lambda tup : tup[0]
        primerDB_c.execute("SELECT distinct primer_id FROM amplicons")
        primer_ids = map(tup2str, primerDB_c.fetchall())
    
        for mismatch in range(5):
            amplicons = {}
            for primer_id in primer_ids:
                primerDB_c.execute("""SELECT seq_id FROM amplicons
                                      WHERE primer_id = "{0}" 
                                      AND rev_mis <= {1} 
                                      AND fw_mis <= {1}""".format(primer_id, mismatch))
                amplicons[primer_id] = map(tup2str, primerDB_c.fetchall())
            primers = [ k  for  k  in sorted(amplicons, key = lambda primer : len(amplicons[primer]), reverse = True)]
            self.barchart(amplicons, primers, mismatch)
            self.venn(amplicons, primers, mismatch)
            
    def barchart(self, amplicons, primers, mismatch):

        fig, ax = plt.subplots() 
        get_perc = lambda primer : round(len(amplicons[primer])/float(self.ref_len) * 100, 1)
        y_pos = np.arange(len(primers))
        coverage = np.array([ get_perc(k) for k in primers ])
        colors = cm.autumn(coverage / float(max(coverage)))
        
        plt.bar(y_pos, coverage, align='center', alpha=0.5, color = colors)
        plt.xticks(y_pos, primers)
        font = {'weight': 'bold', 'size': 14}
        plt.xlabel('Primer pair ID', fontdict = font)
        plt.ylabel('Amplicon coverage (%)',fontdict = font)
        plt.setp(plt.xticks()[1], rotation=30, ha='right')
        fig_fname = '.'.join(['coverage_mismatch{}'.format(mismatch),'pdf'])
        fig.tight_layout()
        plt.savefig(os.path.join(self.output_dir, fig_fname), dpi=1000)
        plt.close()
        plt.clf()
        
                   
    def venn(self, amplicons, primers, mismatch):

        fig, ax = plt.subplots()
        
        n = len(primers)
        if  n > 2:
            primers = primers[:3]
            venn_data = venn3_unweighted(map(lambda k: set(amplicons.get(k)), primers), set_labels = primers)
        elif n ==  2:
            venn_data = venn2_unweighted(map(lambda k: set(amplicons.get(k)), primers), set_labels = primers)
        else:
            return

        [ text.set_fontsize(16) for text in venn_data.set_labels ]
        fig_fname = '.'.join(['venn_mismatch{}'.format(mismatch),'pdf'])
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, fig_fname), dpi=1000)
        plt.close()
        plt.clf()

        
Db_import = load_SQL_csv()
Db_import.import_tsv()
if Db_import.amplicons_loaded:
    Db_import.AmpliconPlots()
    
