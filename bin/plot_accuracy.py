#!/usr/bin/env python
from matplotlib.ticker import StrMethodFormatter
from init_primerDB import PrimerDB
import matplotlib.pyplot as plt
import mysql.connector
import pandas as pd



class PlotAccuracy(PrimerDB):

    def __init__(self):
        
        super(PlotAccuracy, self).__init__()
        self.cnx.database = self.DB_NAME  


    def get_data(self):

        self.cursor.execute('SELECT primer_pair,accuracy FROM taxa_assignment')
        table_rows  = self.cursor.fetchall()
        df  = pd.DataFrame(table_rows, columns=["primer_pair", "accuracy"])
        fig, axs = plt.subplots(figsize=(8,4))
        flierprops = dict(markersize=1,
                          linestyle='none',
                          markeredgecolor='g')
        plot = df.boxplot(column = 'accuracy',
                                by = 'primer_pair',
                                grid = False,
                                showmeans=True,
                                notch = True,
                                fontsize = 10,
                                flierprops=flierprops)

        plt.title("")
        plt.suptitle("")
        plt.xticks(rotation=45, horizontalalignment="right")
    
        plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
        plt.ylabel('Confidence score', fontsize=14, fontname = 'sans-serif',  weight = 'bold')
        plt.xlabel('Primer pair', fontsize=14, fontname = 'sans-serif', weight = 'bold')
        plt.tight_layout()
        plt.savefig("assignment_accuracy.pdf")
        df.groupby('primer_pair').describe().to_csv("assignment_accuracy.tsv",
                                     float_format= "%.2f",
                                     sep= "\t")



if  __name__ == '__main__':
    plot_accuracy  = PlotAccuracy()
    plot_accuracy.get_data()
