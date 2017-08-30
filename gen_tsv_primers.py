#!/usr/bin/python

from generate_primers import PiscatorDB_Primers
import sys
import os

class PrimerProspector(PiscatorDB_Primers):
      
    def Primer_pair_amplicons(self, ext = '.tsv'):
         
        """
           Generate a tsv of primer for primer_prospector            

        """
        
        primer_pairs = self.get_primer_tsv()
        fix_id = lambda primer_id, Dir: primer_id if primer_id.endswith(Dir) else primer_id+Dir
        save   = lambda primer_id, primer: open(primer_id + ext,'w').write(primer)
        output_dir = self.output_dir
        out = lambda primer_id: os.path.join(output_dir, primer_id + ext)
        self.output = []
        self.pair_data = []
        for fwd_id, fwd, rev_id, rev in primer_pairs:
            fwd_id, rev_id  = map(lambda x:\
                                  x.replace('_','.'),  (fix_id(fwd_id.lower(),'f'), fix_id(rev_id.lower(),'r')))
            fwd =  '\t'.join([fwd_id,fwd])
            rev =   '\t'.join([rev_id,rev])
            save(fwd_id, fwd)
            save(rev_id, rev)
            to_nextflow = map(out, [fwd_id, rev_id])
            self.pair_data.append([fwd_id, rev_id])
            self.output.append(','.join(to_nextflow)+'\n')

    def send2stdout(self):
        
        """
            simple utility to send to print out to nextflow

        """
        sys.stdout.writelines(self.output)
    
if __name__ == '__main__':
    PP = PrimerProspector()
    PP.Primer_pair_amplicons() 
    PP.send2stdout()
