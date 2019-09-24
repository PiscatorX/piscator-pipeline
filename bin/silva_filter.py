#!/usr/bin/env python

from Bio import SeqIO
import argparse
import pprint
import sys



class  SilvaFilter(object):

    def  __init__(self, args):

        self.seq_data  = SeqIO.parse(args.reference, args.informat) 
        self.select    = args.select
        self.names     = args.names
        self.outfname  = args.outfname
        self.remove    = [ taxon.lower() for taxon in args.remove ]
        self.outformat = args.outformat
        self.taxon_limit = args.taxon_limit
        self.reject   = args.reject
        if args.select:
            self.select_list =  args.select.read().lower().split()
            
        elif args.names:
            self.select_list = [ entry.lower() for entry in args.names ]

        
    def  parse_data(self):

        select_seq_data =  []
        for rec in self.seq_data:
             taxon_data = rec.description.split(';')
             clean_taxon = [ taxon for taxon in taxon_data if taxon.lower() not in self.remove ]
             if self.is_select(clean_taxon):
                 if set(self.reject).intersection(clean_taxon):
                     continue
                 n = len(clean_taxon)
                 #N=6 domain, phylum, class, order, family, and genus.
                 if n < self.taxon_limit:
                     
                     continue
                 rec.description = ';'.join(clean_taxon[:6])
                 select_seq_data.append(rec)
        SeqIO.write(select_seq_data, self.outfname, self.outformat)
        self.outfname.close()

                 
    def is_select(self, taxa_data):

        for taxon in self.select_list:
            if taxon in map(lambda x:x.lower(), taxa_data):
                return True
            
        return False
            
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser("get select silva sequences based on taxonomy provided")
    parser.add_argument('reference', type=argparse.FileType('r'),  help = 'fasta reference file')
    parser.add_argument('-s','--select', type=argparse.FileType('r'), help = 'file containing list of taxonomic names/words on each to filter for')
    parser.add_argument('-f','--informat', default = "fasta")
    parser.add_argument('-F','--outformat', default = "fasta")
    parser.add_argument('-t','--taxon_limit',type = int, default = 6)
    parser.add_argument('-n','--names', nargs='+', help = 'space separated taxonomic ranks')
    parser.add_argument('-o','--outfname',    default = "filter.fasta", type=argparse.FileType('w'), help = 'list of words/taxonomic names to look for in silva taxonomic', required = True)
    parser.add_argument('-r','--remove', nargs='+', default = ['SAR', 'Alveolata', 'Diatomea', 'Excavata','Discoba','Discicristata','Euglenida'] , help = 'list of taxonomic to remove sequences', required = False)
    parser.add_argument('-x','--reject', nargs='+', default = ["uncultured"] , help = 'list of taxonomic to reject sequences', required = False)
    args = parser.parse_args()
    if not (args.select or args.names):
        parser.error('At least one type of input must be provied, chose either a file for  --select or  command line arguments to  --names. Try --help for more.')
    silva = SilvaFilter(args)
    silva.parse_data()

        #parse_data(args.reference, args.select, args.names,  args.outfname)
