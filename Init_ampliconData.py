#!/usr/bin/python
from Bio import SeqIO
import argparse
import sqlite3
import csv


def get_args():

    parser = argparse.ArgumentParser(description="""Initialise Amplicon database
    and database taxonomic reference data""")

    parser.add_argument('-r','--ref-fasta', dest='ref_file', 
                        action='store', required=True, type=str)
    parser.add_argument('-d','--db-name', dest='db_name', 
                        action='store', required=True, type=str)
        
    return parser.parse_args()
        

def Init_Amplicon_Tables(db_name):
    
    global primerDB_conn, primerDB_c 
    primerDB_conn  = sqlite3.connect(db_name)
    primerDB_c = primerDB_conn.cursor()
    
    primerDB_c.execute("""CREATE TABLE IF NOT EXISTS amplicons
                               (primer_id  VARCHAR,	
                                amplimer   VARCHAR,
                                fw         INT,
                                fw_mis     INT,
                                rev	       INT,
                                rev_mis    INT,
                                len	       INT,
                                seq_id     VARCHAR);""")
    table_col = ["primer_id","amplimer","fwd","fwd_mis",
                      "rev","rev_mis","len","seq_id"]

    primerDB_c.execute("""CREATE TABLE IF NOT EXISTS  Ref_taxa_data
                            (seq_id_ref VARCHAR  PRIMARY KEY,
                             taxon VARCHAR);""")
    primerDB_conn.commit()


def compile_Ref_taxa_data(ref_fasta):
    
    ref_data  = SeqIO.parse(open(ref_fasta),'fasta')
    
    data_col = ["seq_id","taxon"]    

    
    for rec in ref_data:
        data = [rec.id, rec.description ]
        rec_dict = dict(zip(data_col,data))
        sql = """INSERT OR IGNORE INTO Ref_taxa_data (seq_id_ref, taxon)
              values ("{seq_id}","{taxon}") """.format(**rec_dict)
        
        primerDB_c.execute(sql)
        
    primerDB_conn.commit()
        
if __name__ == '__main__':
    args = get_args()
    Init_Amplicon_Tables(args.db_name)
    compile_Ref_taxa_data(args.ref_file)
    
    
