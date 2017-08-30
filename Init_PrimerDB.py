#!/usr/bin/python
import subprocess
import argparse
import sqlite3
import csv



class PrimerDB(object):

    
    def __init__(self):

        """
            Get the filename of csv file with primers and initialise venatorDB database
        """
        
        parser = argparse.ArgumentParser(description="""Compile High-Throughput Sequencing 
        (HTS) Primer database using a csv primer file.""",
        epilog='NOTE: This utility expects sqlite3 to be installed.')

        parser.add_argument('-p','--primers-file', dest='primers_file', action='store', 
                            required=True, type=str)
        parser.add_argument('-d','--db-name', dest='db_name', action='store', required=True, type=str)
        args = parser.parse_args()
        self.csv_reader = csv.reader(open(args.primers_file))
        self.table_col = ["Fwd_id","Forward_Primer","Rev_id","Rev_Primer","Target",
                          "gene","Amplicon_length","technology",
                          "Reference","Cross_ref","notes"]
        
        self.primerDB_conn  = sqlite3.connect(args.db_name)
        self.primerDB_c = self.primerDB_conn.cursor()
        self.primerDB_c.execute("""CREATE TABLE IF NOT EXISTS primers
                                (gene VARCHAR,
                                 Fwd_id  VARCHAR,
                                 Forward_Primer VARCHAR,
                                 Rev_id VARCHAR,
                                 Rev_Primer VARCHAR,
                                 Amplicon_length VARCHAR,
                                 technology VARCHAR,
                                 Reference  VARCHAR,
                                 Cross_ref VARCHAR,
                                 notes TEXT,
                                 primary key (Fwd_id, Fwd_id));""")
        self.primerDB_conn.commit()
        
    def compileDB(self):
        
        #must skip the first line of csv file
        self.csv_reader.next()
        for row in self.csv_reader:
            row_dict = dict(zip(self.table_col,row))
            self.DB_Insert(row_dict)
        self.primerDB_conn.commit()  


    def DB_Insert(self, row_dict):

        #Inserts rows into the primers table one row at time
        #slow for large datasets, sqlite3 .import csv function is faster.
        #TO DO implement sqlite csv import
        sql = """INSERT OR IGNORE INTO primers (gene,Fwd_id,Forward_Primer,
                 Rev_id,Rev_Primer,Amplicon_length, technology,
                 Reference,Cross_ref, notes)
              values ('{gene}','{Fwd_id}','{Forward_Primer}',
                 '{Rev_id}','{Rev_Primer}','{Amplicon_length}','{technology}',
                 '{Reference}','{Cross_ref}','{notes}');""".format(**row_dict)

        self.primerDB_c.execute(sql)




PiscatorDB = PrimerDB()
PiscatorDB.compileDB()
