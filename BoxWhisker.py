#!/usr/bin/python
import matplotlib.pyplot as plt
from scipy.stats import mstats
import itertools as it
import numpy as np
import random
import glob
import sys
import pprint

class AlignDist(object):
    
    def __init__ (self, f_list, top = 3):
        
        self.f_list = f_list
        
        self.f_data = {}
        
        for f_name in self.f_list:
            self.f_data[f_name.replace('.dat','')] = map(float, open(f_name).read().splitlines())
        
        self.data_len = dict((name, len(data)) for name, data in self.f_data.items() )
        self.labels = sorted(self.data_len, key=self.data_len.get, reverse=True)
        self.dist_data = [ self.f_data[k] for k in self.labels ]
        self.locations = dict((i,j) for j,i in enumerate(self.labels,1))
        
        data = []
        for array in self.dist_data:
             data.extend(array)

        self.y_max, self.y_min = map(lambda func: func(data), [max, min])
        
        self.y_max = self.y_max + 0.075 * top
        
        self.top = self.labels[:top]


    def Init_plt(self):
                
        fig,  axs = plt.subplots(1,1,figsize=(6,4))
        plot = axs.boxplot(self.dist_data,
                           positions=self.locations.values(),
                           patch_artist=True)
        axs.grid(axis='y',  # set y-axis grid lines
        linestyle='--',     # use dashed lines
        which='major',      # only major ticks
        color='lightgrey',  # line colour
        alpha=0.7)          # make lines semi-translucent

        box_colors = ['#008000', 'LightGreen', 'Plum',
               'SkyBlue', 'cyan', 'Plum',
               'SkyBlue', 'LightGreen', 'Plum']

        for box, colour in zip(plot['boxes'], box_colors):
            plt.setp(box, color='DarkMagenta', linewidth=0.5, facecolor=colour)
        
        for flier in plot['fliers']:
            flier.set(marker='o', color='blue', alpha=0.5)

        axs.set_xticklabels(self.labels, fontsize=10, rotation=-20, ha='left')
        axs.yaxis.set_tick_params(labelsize=10)
        plt.ylabel('Pairwise genetic distance')
        plt.xlabel('Primer pair ID')
        #fig.savefig('AlignDist.pdf', bbox_inches='tight')
        axs.set_ylim([0, self.y_max])
        self.axs = axs
        
    def annotate(self, plt_fname = 'primer_var.pdf'):
        
        for i,j in it.combinations(self.top,2):
            data_i = self.f_data[i]
            data_j = self.f_data[j]
            stars = '*'*random.randint(0,5)
            H_value, p_value = mstats.kruskalwallis(data_i, data_j)
            self.y_max = self.y_max - 0.075
            i = self.locations[i]
            j = self.locations[j]
            self.axs.annotate("",
            xy=(i, self.y_max), xycoords='data',
            xytext=(j, self.y_max), textcoords='data',
            color='green',
            arrowprops=dict(arrowstyle="-",
                            connectionstyle= "bar, fraction=0.075",
                            linewidth=2,
                            color='red'))
            self.axs.text((i + j)/2.0 ,
                self.y_max+0.065,
                self.get_stars(p_value),
                horizontalalignment='center',
                verticalalignment='center')
        plt.savefig("{}".format(plt_fname), bbox_inches='tight')
            

    def get_stars(self, p):

        if p < 0.0001:
            return "****"
        elif (p < 0.001):
            return "***"
        elif (p < 0.01):
            return "**"
        elif (p < 0.05):
            return "*"
        else:
            return "-"
                          
if __name__ == '__main__':
        var_dist  = AlignDist()
        var_dist.Init_plt()
        var_dist.annotate()
