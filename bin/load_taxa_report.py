#!/usr/bin/env python
from init_primerDB import PrimerDB
import mysql.connector
import argparse
import csv
import sys



class TaxonomyReport(PrimerDB):
    
    """
       Get taxonomy report produced by primer prospector and load into database

   """
    
    def __init__(self):
        
        super(TaxonomyReport, self).__init__()
        parser = argparse.ArgumentParser(description="show number of sequences")
        parser.add_argument('-r', '--report', type=argparse.FileType('r'),  required = True)
        self.args, unknown = parser.parse_known_args()
        self.cnx.database = self.DB_NAME
        fieldnames = ['amplicon_id','taxonomy','accuracy'] 
        self.report_csv = csv.DictReader(self.args.report, fieldnames = fieldnames, delimiter = "\t")
        self.primer_pair = self.args.report.name.replace("_amplicons_assignments.txt","")
        

    def get_report(self):

        
         for i,row in enumerate(self.report_csv, 1):
             amplicon_id, accuracy = (row['amplicon_id'], float(row['accuracy']),)
             sql = """INSERT IGNORE   INTO taxa_assignment (primer_pair, amplicon_id, accuracy) VALUES (%s, %s, %s)"""
             try:
                 self.cursor.execute(sql,(self.primer_pair, amplicon_id, accuracy))
             except mysql.connector.errors.IntegrityError as err:
                print(err)
         self.cnx.commit()


         
if  __name__ == '__main__':
    taxonomy_report = TaxonomyReport()
    taxonomy_report.get_report()
    


