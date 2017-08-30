#! /usr/bin/python

from BoxWhisker import AlignDist
import argparse
import glob
import os

class GenPlots(object):

    def __init__(self, file_patt= "*.dat"):
               
        parser = argparse.ArgumentParser(description="Generate plots")
        
        parser.add_argument('-w','--working-dir', dest='working_dir', 
                            action='store', required=True, type=str)
        args = parser.parse_args()
        
        print args.working_dir
        
        os.chdir(args.working_dir)
        f_list = (glob.glob(file_patt))
        cluster_id = [cluster.split('r_')[-1].replace('.dat','') for cluster in f_list ]
        self.cluster_fnames = {}
        for perc_id  in cluster_id:
            self.cluster_fnames[perc_id] =\
            [f_name  for f_name in  f_list if perc_id in f_name  ]
    
    def gen_plots(self, AlignDist, plot_type ='pdf' ):
        
        for perc_id in self.cluster_fnames:
            plt_fname = ''.join(['PrimerVar_',perc_id,'.',plot_type])
            data_files = self.cluster_fnames[perc_id]
            if len(data_files) == 1:
                continue
            box_plots = AlignDist(data_files)
            box_plots.Init_plt()
            box_plots.annotate(plt_fname)

if __name__  == '__main__':
    primer_var  = GenPlots()
    primer_var.gen_plots(AlignDist)
