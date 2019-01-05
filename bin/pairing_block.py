#!/usr/bin/python
from datetime import datetime
import time
import argparse
import os

class Pairing_block(object):

    def __init__(self):
               
        parser = argparse.ArgumentParser(description="Block untill all files are available")
        parser.add_argument('-p', '--file-pair', dest='file_pair', action='store',
                            required=True, type=str, nargs=2)
        
        self.args = parser.parse_args()

        
    def blocker(self):
        
        t1 = datetime.now()
        time.sleep(1)
        while len(self.args.file_pair) != 0:
            for mate_file in self.args.file_pair:
                #print mate_file
                if os.path.exists(mate_file):
                     t2 = datetime.now()
                     delta = t2 - t1
                     if delta.seconds > 1800:
                         raise Exception,"Pairing of hits took too long"
                     self.args.file_pair.remove(mate_file)
                     f_name = os.path.basename(mate_file)
                     f_name = os.path.join(os.getcwd(), f_name)
                     #print mate_file
                     #print f_name
                     os.symlink(mate_file, f_name)

                    
block = Pairing_block()
block.blocker()
