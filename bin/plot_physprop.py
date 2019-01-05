#!/usr/bin/python

from matplotlib import pyplot as plt
from init_primerDB import PrimerDB
import seaborn as sns
import pandas as pd
import argparse
import sqlite3
import pprint
import ast
import os


class PhyspropPlot(PrimerDB):

    def  __init__(self):

        super(PhyspropPlot, self).__init__()
        self.cnx.database = self.DB_NAME  
        self.cursor =  self.cnx.cursor()

        
    def get_data(self):

        self.col_ids = ['primer', 'Length', 'Tm', 'TmProd', 'GC', 'DeltaS', 'DeltaG', 'DeltaH' ]
        self.cursor.execute('select {} from primerprop'.format(",".join(self.col_ids)))
        results = iter(self.cursor.fetchall())
        self.primer_data = {}
        for primer_prop in results:
             primer_id = primer_prop[:1]
             primer_fields = list(primer_id) +[ map(float, field)\
                            for field in map(ast.literal_eval, primer_prop[1:])]
             self.primer_data[primer_id[0]] = dict(zip(self.col_ids, primer_fields))

             
    def PlotData(self):
        
        df_list = []
        for prop in self.col_ids[1:]:
            data = {}
            for primer_id, primer_prop in self.primer_data.items():
                prop_data = primer_prop[prop]
                tuples  =  [ (prop, i) for i in range(len(prop_data)) ]
                index = pd.MultiIndex.from_tuples(tuples, names=['Propty', 'Ref'])
                data[primer_id] = pd.Series(prop_data, index = index)
            df_list.append(pd.DataFrame(data))
        df_cat = pd.concat(df_list)
        df_x = df_cat.T

        prop_list = self.col_ids[1:]
        for prop in prop_list:
            f_name = ''.join(['physchem_',prop, '.pdf'])
            self.draw(df_x , prop, f_name)
            

    def draw(self, df, prop, fname):

        fig, ax = plt.subplots(1,1)
        fig.subplots_adjust(hspace=.5)
        sns.set_style('whitegrid')
        
        ax = sns.barplot(data=df[prop].T, ax=ax, palette=sns.color_palette('pastel'), capsize=.1)
        for p in ax.patches:    
            ax.annotate("%.1f" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()),
                        ha='center', va='center', xytext=(0, 10), color='blue', size=12, textcoords='offset points')
            
        plt.xticks(fontsize=16, rotation=45, ha="right")
        plt.xlabel('Primer ID', fontsize=18, fontweight='bold')
        plt.ylabel(prop, fontsize=18, fontweight='bold')
        
        plt.tight_layout()
        figure = ax.get_figure()
        figure.savefig(fname)
        
        
if __name__ == '__main__':
    pltx = PhyspropPlot()
    pltx.get_data()
    pltx.PlotData()
