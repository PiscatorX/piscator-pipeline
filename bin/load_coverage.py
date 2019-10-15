#!/usr/bin/env python
from init_primerDB import PrimerDB
import mysql.connector
import argparse
import csv
import sys



class TaxonomyCoverage(PrimerDB):
    
    """
       Get taxonomy coverage produced by primer prospector and load into database

   """
    
    def __init__(self):
        
        super(TaxonomyCoverage, self).__init__()
        parser = argparse.ArgumentParser(description="Get taxonomy coverage produced by primer prospector and load into database")
        parser.add_argument('-r', '--report', type=argparse.FileType('r'),  required = True)
        self.args, unknown = parser.parse_known_args()
        self.cnx.database = self.DB_NAME
        #fieldnames = ['level','first_level_taxonomy','current_taxonomy','total_seqs','percent_passing']
        #taxonomy level, first level taxonomy, taxonomy classification for current level, total seqs for given classification, percent seqs passing score
        #self.report_csv = csv.DictReader(self.args.report, fieldnames = fieldnames, delimiter = "\t")
        self.report_csv = csv.reader(self.args.report)
        self.primer_pair = "_".join(self.args.report.name.split("_")[0:4:3])
        

    def get_report(self):
        
         for i,row in enumerate(self.report_csv, 1):
             if row[0].startswith('#'):
                 continue
             values = [self.primer_pair] + row
             
             sql = """INSERT IGNORE  INTO taxa_coverage (primer_pair, level, first_level, current_taxonomy, total_seqs, percent_passing) VALUES (%s, %s, %s, %s, %s, %s)"""
             try:
                 self.cursor.execute(sql, values)
             except mysql.connector.errors.IntegrityError as err:
                print(err)
         self.cnx.commit()


         
if  __name__ == '__main__':
    taxonomy_coverage = TaxonomyCoverage()
    taxonomy_coverage.get_report()
    


