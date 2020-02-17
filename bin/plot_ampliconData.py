#!/usr/bin/env python
import os
import matplotlib as mpl
if not os.environ.get('DISPLAY',''):
    mpl.use("Agg")
from matplotlib_venn import venn3_unweighted, venn2_unweighted
import matplotlib.ticker as ticker
import itertools as it
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
import sys
import os



class PlotData(PrimerDB):
    
    def __init__(self):

        super(PlotData, self).__init__()
        parser = argparse.ArgumentParser(description="Plot and analyse amplicon data")
        
        parser.add_argument('-r','--ref', dest='ref_file', action='store', required=True, type=str)
        args, unknown = parser.parse_known_args() 
        self.cnx.database = self.DB_NAME  
        self.cursor =  self.cnx.cursor()
        self.ref_len = len(list(SeqIO.parse(args.ref_file,'fasta')))
        self.tup2str = lambda tup : tup[0]
        self.cursor.execute("SELECT distinct primer_id FROM amplicons")
        self.primer_ids = map(self.tup2str, self.cursor.fetchall())
    
        
        
    def PlotsLens(self):

        sns.set_style("white")
        step  = 10
        font = {'weight': 'bold', 'size': 14}
        for mismatch in range(4):
            ampli_lens = {}
            for primer_id in self.primer_ids:
                print primer_id
                self.cursor.execute("""SELECT len FROM amplicons
                                      WHERE primer_id = "{0}" 
                                      AND rev_mis <= {1} 
                                      AND fwd_mis <= {1}""".format(primer_id, mismatch))
                ampli_lens[primer_id] = pd.Series(map(self.tup2str, self.cursor.fetchall()))
                data = ampli_lens[primer_id]
                xticks = np.arange(min(data)-step,  max(data) + step,  step)
                ax = sns.distplot(data, kde=False, norm_hist = False)
                ax.set_title(primer_id)
                ax.set(xticks=xticks)
                ax.set_ylabel("Number of amplicons",  fontname = 'sans-serif', fontsize = 14, weight='bold')
                ax.set_xlabel("Length (bp)", fontsize = 14, fontname = 'sans-serif', weight='bold')
                plt.xticks(fontsize=10, rotation=45, ha="right")
                plt.tight_layout()
                
                plt.savefig('len_dist_'+'.'.join([primer_id,'pdf']))
                plt.close()
                plt.clf()
            len_df = pd.DataFrame(ampli_lens)
            len_df.describe().transpose().to_csv("amplicon_Lengths_mismatch{}.tsv".format(mismatch),
                                     float_format= "%.2f",
                                     sep= "\t")
            ax = sns.boxplot(data=len_df, orient='h',  linewidth = 0.5,  fliersize = 0.75,)

            plt.xlabel('Length (bp)', fontsize=14, fontname = 'sans-serif',  fontweight='bold')
            plt.ylabel('Primer pair', fontsize=14, fontname = 'sans-serif',  fontweight='bold')
            plt.xticks(fontsize=10, ha="right")
            plt.yticks(fontsize=10)
            plt.grid()
            #plt.yticks(fontsize=8, rotation=-45, ha="right")
            
            plt.tight_layout()
            plt.savefig("amplicon_Lengths_mismatch{}.pdf".format(mismatch))
            plt.close()
            plt.clf()

               
        
    def AmpliconPlots(self):

        
        for mismatch in range(5):
            amplicons = {}
            for primer_id in self.primer_ids:
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
        
        get_perc = lambda primer : 100 * round(len(amplicons[primer])/float(self.ref_len), 3)
        y_pos = np.arange(len(primers))
        coverage = np.array([ get_perc(k) for k in primers ])
        coverage_df = pd.DataFrame(index = primers, data ={'coverage': coverage})
        coverage_df.to_csv("coverage_mismatch{}.tsv".format(mismatch),
                                     float_format= "%.2f",
                                     sep= "\t")

        
        plt.bar(y_pos, coverage, align='center', alpha=0.5, color = 'grey',  edgecolor ='black',  linewidth = 1)
        plt.xticks(y_pos, primers)
        
        font = {'weight': 'bold', 'size': 14}
        plt.xlabel('Primer pair', fontdict = font)
        plt.ylabel('Amplicon coverage (%)',fontdict = font)
        plt.setp(plt.xticks()[1], rotation=30, ha='right')
        fig_fname = 'coverage_mismatch{}.pdf'.format(mismatch)
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
    plots.AmpliconPlots()
    plots.PlotsLens()
