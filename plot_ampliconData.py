#!/usr/bin/python

from matplotlib_venn import venn3_unweighted, venn2_unweighted
from init_primerDB import PrimerDB
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from Bio import SeqIO
import numpy as np
import subprocess
import argparse
import sqlite3
import pprint
import os


class PlotData(PrimerDB):
    
    def __init__(self):

        super(PlotData, self).__init__()
        parser = argparse.ArgumentParser(description="Plot and analyse amplicon data")
        
        parser.add_argument('-r','--ref', dest='ref_file', action='store', required=True, type=str)
        
        args = parser.parse_args()
        self.cnx.database = self.DB_NAME  
        self.cursor =  self.cnx.cursor()
        self.ref_len = len(list(SeqIO.parse(args.ref_file,'fasta')))
        
         
    def AmpliconPlots(self):

        tup2str = lambda tup : tup[0]
        self.cursor.execute("SELECT distinct primer_id FROM amplicons")
        primer_ids = map(tup2str, self.cursor.fetchall())
    
        for mismatch in range(5):
            amplicons = {}
            for primer_id in primer_ids:
                self.cursor.execute("""SELECT seq_id FROM amplicons
                                      WHERE primer_id = "{0}" 
                                      AND rev_mis <= {1} 
                                      AND fwd_mis <= {1}""".format(primer_id, mismatch))
                amplicons[primer_id] = map(tup2str, self.cursor.fetchall())
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
        plt.savefig(fig_fname, dpi=1000)
        plt.close()
        plt.clf()
        
                   
    def venn(self, amplicons, primers, mismatch):

        fig, ax = plt.subplots()
        
        n = len(primers)

        print  "aFA?sfA?fa",n
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
        plt.savefig(fig_fname, dpi=1000)
        plt.close()
        plt.clf()


if  __name__ ==  '__main__':    
    plots = PlotData()
    plots.AmpliconPlots()
    
