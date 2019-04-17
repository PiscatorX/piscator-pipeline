#!/usr/bin/env python

from mysql.connector import errorcode
import mysql.connector
import argparse
import pprint
import sys






class PrimerDB(object):

    def __init__(self):

        """

             Get the filename of csv file with primers and initialise Piscator database

        """
        config = {'user':'root',
             'password' :'Pythonic',
              'database':'piscator'}

       # 'host': "localhost",
  
        #config file stores all the database
        parser = argparse.ArgumentParser(description="""Compile High-Throughput Sequencing 
        (HTS) Primer database using a csv primer file.""",
        epilog='NOTE: This utility expects database to be created and cvs files to have headers')

        parser.add_argument('host')

        self.args, unknown = parser.parse_known_args()
        
        config.update({'host': self.args.host } )
        
        self.DB_NAME = config.pop('database')

        self.cnx = mysql.connector.connect(user=config['user'], password=config['password'], host=config['host'])

        self.cursor = self.cnx.cursor()

        self.config = config
        

        
    def create_DB(self):

        
        try:
            self.cursor.execute("CREATE DATABASE {} DEFAULT CHARACTER SET 'utf8'".format(self.DB_NAME))
        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_DB_CREATE_EXISTS:
                print('Database "{}" already exists'.format(self.DB_NAME)) 
            else:
                raise Exception(err)

            
    def init_tables(self):

        TABLES = {}
        
        TABLES['primers'] = (
             " CREATE TABLE `primers` ("
             " `gene` CHAR(255) NOT NULL,"
             " `Fwd_id`  CHAR(30) NOT NULL,"
             " `Fwd_Primer` CHAR(255) NOT NULL,"
             " `Rev_id` CHAR(30) NOT NULL,"
             " `Rev_Primer` CHAR(255) NOT NULL,"
             " `Amplicon_length` CHAR(255),"
             " `technology` CHAR(255),"
             " `Reference`  CHAR(255),"
             " `Cross_ref`  CHAR(255),"
             " `notes` LONGTEXT,"
             " primary key (`Fwd_id`, `Rev_id`)"
             " ) ENGINE=InnoDB")
        
        
        TABLES['primerprop'] = (
             " CREATE TABLE `primerprop` ("
             " `primer` CHAR(30) PRIMARY KEY,"
             " `TmProd` LONGTEXT NOT NULL," 
             " `DeltaS` LONGTEXT NOT NULL," 
             " `Length` LONGTEXT NOT NULL," 
             " `GC`     LONGTEXT NOT NULL," 
             " `DeltaG` LONGTEXT NOT NULL," 
             " `DeltaH` LONGTEXT NOT NULL," 
             " `Tm`     LONGTEXT NOT NULL"
             " ) ENGINE=InnoDB")

        TABLES['amplicons'] = (
             " CREATE TABLE `amplicons` ("
             " `primer_ID` CHAR(255) NOT NULL,"	
             " `amplimer`  CHAR(255) NOT NULL,"	
             " `fwd`        INT NOT NULL,"	
             " `fwd_mis`    INT NOT NULL,"	
             " `rev`       INT NOT NULL,"	
             " `rev_mis`   INT NOT NULL,"	
             " `len`	   INT NOT NULL,"
             " `seq_id`    VARCHAR(1000) NOT NULL"
             " ) ENGINE=InnoDB")

        
        TABLES['taxa_data'] = (
             " CREATE TABLE `taxa_data` ("
             " `seq_id` CHAR(30)  PRIMARY KEY NOT NULL,"
             " `taxonomy` LONGTEXT NOT NULL"
             " ) ENGINE=InnoDB")

        self.TABLES = TABLES
        self.DB_tables = list(TABLES.keys())


    def create_tables(self):
        
        try:
            self.cnx.database = self.DB_NAME  
        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_BAD_DB_ERROR:
                print('Database Error')
                raise Exception(err)
            
        for table, ddl in list(self.TABLES.items()):
            try:
                self.cursor.execute(ddl)
            except mysql.connector.Error as err:
                
                if err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
                    print('Table "{}" already exists.'.format(table))
                    
                else:
                    raise Exception(err)
                
if __name__ == '__main__':
    primer_db = PrimerDB()
    primer_db.init_tables()
    primer_db.create_DB()
    primer_db.create_tables()

    
