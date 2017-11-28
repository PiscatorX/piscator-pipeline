#! /usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import itertools as it
import operator as o
import numpy as np
import argparse
import glob
import os


class ClusterAnalysis(object):

    def  __init__(self):
        
        parser = argparse.ArgumentParser(description="Analyse cdhit clusters")    
        parser.add_argument('-w','--working-dir', dest='working_dir', action='store',
                            help="directory location of CDHit files, default is current directory", required=False, default='', type=str)
        parser.add_argument('-p','--plot-type', dest='plot_type', action='store',default='pdf',
                            help="matplotlib supported image format (default=pdf)", required=False, type=str)

        args = parser.parse_args()
        self.working_dir = args.working_dir
        self.file_patt = os.path.join(self.working_dir, '*.cd_hits.clstr')
        self.plot_type = args.plot_type
    

        
    def plot_clusters(self):

        self.cluster_data = []
        for clstr_file in glob.iglob(self.file_patt):
            f_name = os.path.basename(clstr_file)
            perc = f_name.split('.cd_hits.clstr')[0].split('_')[-1]
            clstr_len = len(filter(lambda line : line.startswith('>C'),\
                                   open(clstr_file).read().splitlines()))
            self.cluster_data.append([f_name.rsplit('_',2)[0],
                                 int(float(perc)*100),#remove trailing zeros
                                 clstr_len])
        if not self.cluster_data:
            raise Exception, "No cd_hit '*.cd_hits.clstr' files found"
        
        dpoints = np.array(self.cluster_data)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        space = 0.3

        conditions = [(c, np.mean(dpoints[dpoints[:,0] == c][:,2].astype(float)))
                                                for c in np.unique(dpoints[:,0])]
        categories = [(c, np.mean(dpoints[dpoints[:,1] == c][:,2].astype(float)))
                                                for c in np.unique(dpoints[:,1])]
        conditions = [c[0] for c in sorted(conditions, key=o.itemgetter(1))]
        #pprint.pprint(conditions)
        categories = [c[0] for c in sorted(categories, key=o.itemgetter(1))]
        #pprint.pprint(categories)
        dpoints = np.array(sorted(dpoints, key=lambda x: categories.index(x[1])))
        
        n = len(conditions)

        width = (1 - space)/(len(conditions))
        patterns = it.cycle(["+", "x" ".", "*" , "/" , "\\" , "|" , "-"  , "o", "O"])
        for i,cond in enumerate(conditions):
            indices = range(1, len(categories)+1)        
            vals = dpoints[dpoints[:,0] == cond][:,2].astype(np.float)
            pos = [j - (1 - space) / 2. + i * width for j in range(1,len(categories)+1) ]    
            ax.bar(pos, vals, width=width, label=cond, color=cm.Accent(float(i) / n - i/n), hatch=patterns.next())

        ax.set_xticks(indices)
        ax.set_xticklabels(categories)
        plt.setp(plt.xticks()[1], rotation=0)
        plt.rc('legend',**{'fontsize':8})
        font = { 'fontname':'sans-serif', 'weight': 'bold', 'size': 14}
        ax.set_ylabel("Number of CD-HIT clusters", fontdict=font)
        ax.set_xlabel("Clustering percent (%)", fontdict=font)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='upper left')
        plt.savefig(os.path.join(self.working_dir, '.'.join(['CDhit_counts', self.plot_type])))              


if __name__ ==  '__main__':        
    clusters = ClusterAnalysis()
    clusters.plot_clusters()
    
    
    
