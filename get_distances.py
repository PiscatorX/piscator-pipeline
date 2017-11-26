#! /usr/bin/env python
import matplotlib.pyplot as plt
from scipy.stats import mstats
import itertools as it
import numpy as np
import subprocess
import argparse
import glob
import os 




class GetDistances():

    def  __init__(self):

        parser = argparse.ArgumentParser(description="compile a distance matrix")    
        parser.add_argument('-t','--threads', dest='threads', default=1, help = "ClustalO threads",
                            action='store', required=False, type=int)
        parser.add_argument('-w','--working-dir', dest='working_dir', action='store',
                            help="directory location of CDHit files, default is current directory", required=False, default=os.getcwd(), type=str)
        parser.add_argument('-p','--plot-type', dest='plot_type', action='store',default='pdf',
                            help="matplotlib supported image format (default=pdf)", required=False, type=str)

        args = parser.parse_args()
        self.threads = args.threads
        self.working_dir = args.working_dir
        #print self.working_dir
        self.file_patt = os.path.join(self.working_dir, '*.cd_hits.clstr')
        self.plot_type = args.plot_type
        args = parser.parse_args()
        
        self.dist_list  = []
        self.distance_arrays = {}

        
       
    def init_get_dist(self):

        file_patt  = self.file_patt.replace('.clstr','')
        map(self.get_dist, glob.iglob(file_patt))
         
        
        
    def get_dist(self, file_in):
        
        aln_fname =  file_in.replace('cd_hits','aln')
        dist_fname = file_in.replace('cd_hits','dist')
        args = "clustalo --in {0} --out {1} --force --distmat-out {2}\
        --full --threads {3}".format(file_in, aln_fname, dist_fname, self.threads)     
        
        proc = subprocess.Popen(args,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = proc.communicate()
        
        if  stderr:
            print stderr
            return
        self.dist_list.append(dist_fname.replace('dist','dat'))
        self.parse_dist(dist_fname)


        
    def parse_dist(self, dist_fname):

        
        split_row  =  lambda  row: row.split()
        seq_ids = {}
        
        try:
            with open(dist_fname) as dist_file:
                print 'Parse',dist_fname
                n = int(dist_file.next())    
                data_matrix = np.zeros((n,n))
                for i, row in enumerate(dist_file, 0):
                    data_matrix[i] = split_row(row)[1:]
                    
                self.distance_arrays[dist_fname] = data_matrix[np.triu_indices(n)]
        except IOError as err:
            print err
            return

        
    def compile_data(self):

        
        self.f_list = self.distance_arrays.keys()
        
        print 'Compile data',self.f_list
        
        cluster_ids = [cluster_fname.split('r_')[-1].replace('.dist','') for cluster_fname in self.f_list ]
        self.cluster_fnames = {}
        
        for perc_id  in set(cluster_ids):
            self.cluster_fnames[perc_id] =\
            [f_name  for f_name in  self.f_list if perc_id in f_name  ]
        
        for perc_id in self.cluster_fnames:
            plt_fname = ''.join(['PrimerVar_',perc_id,'.', self.plot_type])
            data_files = self.cluster_fnames[perc_id]
            if len(data_files) == 1:
                continue
            cluster_data = {}
            for f_name in self.f_list:
                cluster_data[os.path.basename(f_name.rsplit('_', 1)[0])] = self.distance_arrays[f_name]
            self.Init_plot(cluster_data, plt_fname) 
        
        
    def Init_plot(self, cluster_data, plt_fname, top = 3):

        
        if not cluster_data:
            return
        
        cluster_data = cluster_data
        data_len = dict((name, len(data)) for name, data in cluster_data.items() )

        
        labels = sorted(data_len, key=data_len.get, reverse=True)
        dist_data = [ cluster_data[k] for k in labels ]
        locations = dict((i,j) for j,i in enumerate(labels,1))    

        for x in cluster_data.values():
            print x.shape
        merged_data  = np.concatenate(dist_data)
        y_max, y_min = map(lambda func: func(merged_data), [np.amax, np.amin])
        y_max = y_max + 0.075 * top
        top = labels[:top]

        fig,  axs = plt.subplots(1,1,figsize=(6,4))
        plot = axs.boxplot(dist_data,
                           positions=locations.values(),
                           patch_artist=True,
                           flierprops={'alpha':0.6,
                           'markersize': 3,
                           'color':'#e7298a',
                           'markeredgecolor': 'black',
                           'marker': 'o',
                           'alpha': 0.5},
                            sym="k.")

        axs.grid(axis='y',  # set y-axis grid lines
        linestyle='-.',     # use dashed lines
        which='major',      # only major ticks
        color='lightgrey',  # line colour
        alpha=0.7)          # make lines semi-translucent

        # for box, colour in zip(plot['boxes'], box_colors):
        #     plt.setp(box, color='turquoise', linewidth=0.5, facecolor='blue')
        

        axs.set_xticklabels(labels, fontsize=10, rotation=-20, ha='left')
        axs.yaxis.set_tick_params(labelsize=10)
        
        plt.ylabel('Pairwise genetic distance')
        plt.xlabel('Primer pair ID')

        axs.set_ylim([0, y_max])

        
        for box in plot['boxes']:
            # change outline color
            box.set( color='black', linewidth=1)
            # change fill color
            box.set( facecolor = '#1b9e77' )

        ## change color and linewidth of the whiskers
        for whisker in plot['whiskers']:
            whisker.set(color='black', linewidth=1)

        ## change color and linewidth of the caps
        for cap in plot['caps']:
            cap.set(color='black', linewidth=2)

        ## change color and linewidth of the medians
        for median in plot['medians']:
            median.set(color='red', linewidth=2)

        ## change the style of fliers and their fill
        for flier in plot['fliers']:
            flier.set(marker='o', color='#e7298a', alpha=0.5)

                
        if data_len > 2:
            top = data_len if top < 3 else top
            for i,j in it.combinations(top, 2):
                print i,j
                H_value, p_value = mstats.kruskalwallis(cluster_data[i], cluster_data[j])
                print 'next'
                y_max = y_max - 0.075
                i = locations[i]
                j = locations[j]
                axs.annotate("",
                xy=(i, y_max), xycoords='data',
                xytext=(j, y_max), textcoords='data',
                color='green',
                arrowprops=dict(arrowstyle="-",
                                connectionstyle= "bar, fraction=0.075",
                                linewidth=1,
                                color='blue'))
                axs.text((i + j)/2.0 ,
                    y_max+0.065,
                    self.get_stars(p_value),
                    horizontalalignment='center',
                    verticalalignment='center')

                
        plt.savefig("{}".format(plt_fname), bbox_inches='tight')
        plt.close()

        
        
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


if __name__ ==  '__main__':        
    clusters = GetDistances()
    clusters.init_get_dist()
    clusters.compile_data()
