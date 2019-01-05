#!/usr/bin/python

from init_primerDB import PrimerDB
import mysql.connector
import argparse
import csv


class LoadDB(PrimerDB):

    def __init__(self):
        
        super(LoadDB, self).__init__()

        parser = argparse.ArgumentParser(description="""Compile High-Throughput Sequencing 
        (HTS) Primer database using a csv primer file.""",
        epilog='NOTE: This utility expects database to be created and cvs files to have headers')

        parser.add_argument('-p','--primers-file', dest='primers_file', action='store', 
                            required=True, type=str)
        parser.add_argument('-c','--custom-header', dest='headers', default=False, action='store_true', 
                            required=False, help="""use cvs headers (default=False). Used when order of column header is different from expected. Column/header are expected in the following order: "Fwd_id, Fwd_Primer, Rev_id, Rev_Primer, Target, gene, Amplicon_length, technology, Reference, Cross_ref, notes". Important! header reading is case sensitive, if provided headers must match the case those provided here""")

        self.cnx.database = self.DB_NAME  
        args = parser.parse_args()
        self.headers = args.headers
        self.csv_reader = csv.reader(open(args.primers_file))
        self.table_cols = ["Fwd_id","Fwd_Primer","Rev_id","Rev_Primer",
                          "gene","Amplicon_length","technology",
                          "Reference","Cross_ref","notes"]
        
    def LoadCSV(self):
        
        #must skip the first line of csv file
        top_line = self.csv_reader.next()
        csv_headers = top_line if self.headers else self.table_cols
        
        insert_dict = {}
        for row in self.csv_reader:
            row_dict = dict(zip(csv_headers, map(lambda w: w.strip(), row)))
            if self.headers:
                for k in row_dict:
                    if k in self.table_cols:
                        insert_dict[k] = row_dict[k]
                    else:
                        print "column {} removed from data not inserted into DB (use -h/--help for help)".format(k)
                if  insert_dict:
                    row_dict = insert_dict
                    insert_dict = {}
            self.DB_Insert(row_dict)
        self.cnx.commit()
        
    def DB_Insert(self, row_dict):
        
        cols = ','.join(row_dict.keys()) 
        values = ','.join( '"'+val.translate(None,""""'""")+'"'  for val in  row_dict.values())
        sql = """INSERT INTO primers ({}) VALUES ({})""".format(cols, values)
        
        try:
            self.cursor.execute(sql)
        except mysql.connector.errors.IntegrityError as err:
            #print err
            pass
















LoadDB = LoadDB()
LoadDB.LoadCSV()