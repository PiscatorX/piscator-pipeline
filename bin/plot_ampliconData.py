#!/usr/bin/python

from matplotlib_venn import venn3_unweighted, venn2_unweighted
from init_primerDB import PrimerDB
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from Bio import SeqIO
import seaborn as sns
import pandas as pd
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
        self.tup2str = lambda tup : tup[0]
        self.cursor.execute("SELECT distinct primer_id FROM amplicons")
        self.primer_ids = map(self.tup2str, self.cursor.fetchall())
    
        
        
    def PlotsLens(self):

        set_style("whitegrid")
        for mismatch in range(5):
            ampli_lens = {}
            for primer_id in self.primer_ids:
                print primer_id
                self.cursor.execute("""SELECT len FROM amplicons
                                      WHERE primer_id = "{0}" 
                                      AND rev_mis <= {1} 
                                      AND fwd_mis <= {1}""".format(primer_id, mismatch))
                ampli_lens[primer_id] = pd.Series(map(self.tup2str, self.cursor.fetchall()))
                ax = sns.countplot(ampli_lens[primer_id])
                ax.set_title(primer_id)
                ax.set(ylabel="Number of amplicons", xlabel="Length (bp)")
                plt.xticks(fontsize=10, rotation=45, ha="right")
                plt.tight_layout()
                plt.savefig('Len_dist_'+'.'.join([primer_id,'pdf']))
                plt.close()
                plt.clf()
        len_df = pd.DataFrame(ampli_lens)
        ax = sns.boxplot(data=len_df, orient='h',  palette="Set2")
        ax.set(ylabel="Primer ID", xlabel="Length (bp)")
        plt.xticks(fontsize=10, rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig("Amplicon_Lengths.pdf")
        plt.close()
        plt.clf()

        
    def AmpliconPlots(self):

        
        for mismatch in range(5):
            amplicons = {}
            for self.primer_id in primer_ids:
                self.cursor.execute("""SELECT seq_id FROM amplicons
                                      WHERE primer_id = "{0}" 
                                      AND rev_mis <= {1} 
                                      AND fwd_mis <= {1}""".format(primer_id, mismatch))
                amplicons[primer_id] = map(self.tup2str, self.cursor.fetchall())
            primers = [ k  for  k  in sorted(amplicons, key = lambda primer : len(amplicons[primer]), reverse = True)]
            
            self.barchart(amplicons, primers, mismatch)
            
            self.venn(amplicons, primers, mismatch)

            
    def barchart(self, amplicons, primers, mismatch):
        
        fig, ax = plt.subplots() 
        get_perc = lambda primer : round(len(amplicons[primer])/float(self.ref_len) * 100, 1)
        y_pos = np.arange(len(primers))
        coverage = np.array([ get_perc(k) for k in primers ])
        colors = cm.autumn(1 - coverage/float(max(coverage)))
        
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
    #plots.AmpliconPlots()
    plots.PlotsLens()