#!/usr/bin/env python
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
#db = pyoma.browser.db.Database('/home/ijulcach/Programs/DataBases/OmaServer.h5')
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

### main
parser = argparse.ArgumentParser(description="download fasta file of hogs (root hogs)")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="hog name, starting with D")
args = parser.parse_args()

inFile = args.inFile
hog = 'HOG:'+inFile
# hogs = gmo.load_list(inFile)
outfile = open(inFile+'.members.tsv','w')
data = pyoma.browser.models.HOG(db, hog)
# for hog in hogs:
#     data = pyoma.browser.models.HOG(db, hog)
level = data.level
members = [x.omaid for x in data.members]
string = hog+'\t'+level+'\t'+'; '.join(members)
print(string,file=outfile)
outfile.close()
    
