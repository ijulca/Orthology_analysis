#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 11:29:48 2021

@author: ijulca
"""

import argparse
import subprocess

api_key='67b44932b5a9b8eb90b873da1c6129637508'

### main

parser = argparse.ArgumentParser(description="get the data type of the bioproject")
parser.add_argument("-i", "--id", dest="id", required=True, help="accession number of the bioproject")
args = parser.parse_args()

key = args.id

cmd = "efetch -db BioProject -id "+key+" -format xml | xtract -pattern DocumentSummary -element DataType -first Title"
out = subprocess.getoutput(cmd)

outfile = open(key+'.info','w')
print(key+'\t'+out, file=outfile)
outfile.close()
print('End...')