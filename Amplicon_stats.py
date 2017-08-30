#!/usr/bin/python

import argparse


class Amplicon_stats(object):

    def __init__(self):

        parser = argparse.ArgumentParser(description="""get a primersearch output file
         and analyse the results""")

        parser.add_argument('-a','--amplicons', dest='amplicons', action='store', required=True, type=str)
        args = parser.parse_args()
        self.amplicons_output = open(args.amplicons).read()
        self.results_tsv = args.amplicons.replace('rst','tsv')


    def get_amplicon_stats(self):

        amplicon_ref  = {}
        amplicon_data = iter(self.amplicons_output.splitlines())
        data = [('sequence', lambda x: x.split(':'))]

        # Lay out of amplicon entries in file
        # Sequence: Z52466 Z52466 
        #    H.sapiens (D1S2660) DNA segment containing (CA) repeat; clone AFMa203yc1; single read.
        #    CACACATGCACATGCAC hits forward strand at 27 with 0 mismatches
        #    AGTGACACCAGCAGGG hits reverse strand at [103] with 0 mismatches
        #    Amplimer length: 261 bp

        for line in amplicon_data:
            if line.startswith('Primer name'):
                primer_name = line.split()[-1]
            elif line.startswith('Amplimer'):
                amplicon_vars = self.amplicon_analyse(amplicon_data)
                amplicon_vars.update({'amplimer':line,
                                      'primer_id':primer_name})
                amplicon_ref[line] = amplicon_vars
                fwd, fwd_mis = amplicon_ref[line]['fwd']
                amplicon_ref[line]['fwd'] = fwd
                amplicon_ref[line]['fwd_mis'] = fwd_mis
                rev, rev_mis = amplicon_ref[line]['rev']
                amplicon_ref[line]['rev'] = rev
                amplicon_ref[line]['rev_mis'] = rev_mis

        self.write2csv(amplicon_ref)


    def amplicon_analyse(self, amplicon_data):

        data = [('seq_id',  lambda x: x.split(':')[-1].strip()),
                ('taxa', lambda x: x ),
                ('fwd', lambda x: (int(x.split()[5]), int(x.split()[7]))), 
                ('rev', lambda x: (int(x.split()[5].translate(None,'[]')), int(x.split()[7]))),
                  ('len', lambda x: int(x.split()[2]))]
        amplicon_vars = {}
        for i, line in enumerate(amplicon_data, 0):
            key, func = data[i]
            amplicon_vars[key] = func(line.strip())
            if key == 'len':
                return amplicon_vars


    def write2csv(self, amplicon_ref):
        
        with open(self.results_tsv,'w') as amplicon_tsv:
            for ref in amplicon_ref:
                row = "{primer_id}\t{amplimer}\t{fwd}\t{fwd_mis}\t{rev}\t{rev_mis}\t{len}\t{seq_id}\n"
                row = row.format(**amplicon_ref[ref])
                amplicon_tsv.write(row)
                amplicon_tsv.flush()

if __name__ ==  '__main__':
    Amplicons = Amplicon_stats()
    Amplicons.get_amplicon_stats()
    Amplicons.get_amplicon_stats()
