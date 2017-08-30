#! /usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import itertools as it
import operator as o
import numpy as np
import subprocess
import argparse
import pprint
import glob
import sys
import os



def glob_files():

    parser = argparse.ArgumentParser(description="compile a distance matrix")
    
    parser.add_argument('-o','--output', dest='output', 
                        action='store', required=True, type=str)

    parser.add_argument('-t','--threads', dest='threads', default=1, 
                        action='store', required=False, type=int)

    
    args = parser.parse_args()
    
    patt = os.path.join(args.output,'*.cd_hits')

    count_clusters(args.output)
    
    global threads
    
    threads = args.threads
    
    map(Align_Dist,  glob.iglob(patt))


def count_clusters(output_dir):

    patt = os.path.join(output_dir, '*.cd_hits.clstr')
    cluster_data = []
    for clstr_file in glob.iglob(patt):
        f_name = os.path.basename(clstr_file)
        perc = f_name.split('.cd_hits.clstr')[0].split('_')[-1]
        clstr_len = len(filter(lambda line : line.startswith('>C'),\
                               open(clstr_file).read().splitlines()))
        cluster_data.append([f_name.rsplit('_',2)[0],
                             float(perc),
                             clstr_len])
        
    pprint.pprint(cluster_data)    
    plot_counts(cluster_data, output_dir)

    
def plot_counts(cluster_data, output_dir):

    #pprint.pprint(cluster_data)    
    dpoints = np.array(cluster_data)
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
        indeces = range(1, len(categories)+1)        
        vals = dpoints[dpoints[:,0] == cond][:,2].astype(np.float)
        pos = [j - (1 - space) / 2. + i * width for j in range(1,len(categories)+1) ]    
        ax.bar(pos, vals, width=width, label=cond, color=cm.Accent(float(i) / n), hatch=patterns.next())
     
    ax.set_xticks(indeces)
    ax.set_xticklabels(categories)
    plt.setp(plt.xticks()[1], rotation=0)
    plt.rc('legend',**{'fontsize':8})
    font = {'weight': 'bold', 'size': 14}
    ax.set_ylabel("Number of CD-HIT clusters", fontdict=font)
    ax.set_xlabel("Clustering percent (%)", fontdict=font)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper left')
    plt.savefig(os.path.join(output_dir, 'CDhit_counts.pdf'))

    
def Align_Dist(file_in):

    aln_fname =  file_in.replace('cd_hits','aln')
    dist_fname = file_in.replace('cd_hits','dist')
    data_f  = file_in.replace('cd_hits','dat')
    args = "clustalo --in {0} --out {1} --force --distmat-out {2}\
    --full --threads {3}".format(file_in, aln_fname, dist_fname, threads)     
    #print args
    proc = subprocess.Popen(args,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    stdout, stderr = proc.communicate()

    if  stderr:
        sys.stderr.write(stderr)
       
    stats = parse_dist(dist_fname, data_f)

    

def parse_dist(dist_in, data_f):

    if not os.path.exists(dist_in):
        #print os.getcwd()
        print "failed to locate %s"%dist_in
        return

    data = []
    seq_dict = {}
    with open(dist_in) as fp:
        n = int(fp.next())   
        for i,line in enumerate(fp):
            line = line.split()
            seq_dict[i]=line[0]
            data.append(line[1:])
            
    f_obj = open(data_f,'w')
    dist_data  = [] 
    for i,j in it.combinations(range(n),2):
        #print'\t\t\t\t{0},{1},{2}'.format(seq_dict[i], seq_dict[j],data[i][j])
        print >>f_obj,'{0}'.format(data[i][j])
        #print '{0}'.format(data[i][j])
        dist_data.append(data[i][j])




glob_files()
