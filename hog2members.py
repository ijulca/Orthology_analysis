#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 10:39:05 2024

@author: ijulcach
"""

import argparse
import sys,os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import general_modules as gmo
import pyoma.browser.db
from pyoma.browser.models import ProteinEntry
db = pyoma.browser.db.Database('/work/FAC/FBM/DBC/cdessim2/oma/oma-browser/All.Jul2023/data/OmaServer.h5')
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

### main
parser = argparse.ArgumentParser(description="download fasta file of hogs (root hogs)")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="list of hogs or list of mnemonic")
args = parser.parse_args()

inFile = args.inFile

hogs = gmo.load_list(inFile)
outfile = open(inFile+'.members.tsv','w')
for hog in hogs:
    data = pyoma.browser.models.HOG(db, hog)
    level = data.level
    members = [x.omaid for x in data.members]
    string = hog+'\t'+level+'\t'+'; '.join(members)
    print(string,file=outfile)
outfile.close()
    
