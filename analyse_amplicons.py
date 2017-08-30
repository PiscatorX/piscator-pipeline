#!/usr/bin/python
from generate_primers import PiscatorDB_Primers

class AnalyseAmplicons(PiscatorDB_Primers):

    def get_amplicons(self):
        sql = "select distinct primer_id from amplicons;"
        primer_ids = self.get_results(sql)
        self.primer_ids = ( x[0] for x in primer_ids)
    
    def get_amplicon_data(self):
        
        fix_data = lambda var : ','.join(map(str, var))
        sql = "select len from amplicons where primer_id is '{}';"
        with open('docs/amplicons.len.csv','w') as fp:
            for primer_id  in self.primer_ids:
                primer_id = primer_id.encode('utf8')            
                len_data = self.get_results(sql.format(primer_id))
                data = map(fix_data, len_data)
                data.insert(0, primer_id)
                print data
                #print len_data
                #print primer_id
                #fp.writelines(map(fix_data, len_data))
                break

Amplicons = AnalyseAmplicons()
Amplicons.get_amplicons()
Amplicons.get_amplicon_data()
