#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 19:42:41 2024

@author: ijulcach
"""

import argparse
import sys,os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import general_modules as gmo
from omadb import Client #https://dessimozlab.github.io/pyomadb/build/html/
c = Client()

def get_list_hogs(inFile):
    hogs = set()
    for line in open(inFile):
        line = line.strip()
        data = line.split('.')[0]
        hogs.add(data)
    print('Number of hogs',len(hogs))
    return hogs

def get_fasta_hogs(inFile):
    hogs = get_list_hogs(inFile)
    for g in hogs:
        level = c.hogs[g].level
        page = 'https://oma-stage.vital-it.ch/oma/hog/'
        page += g+'/'+level+'/fasta/'
        cmd = 'wget '+page
        gmo.run_command(cmd)
        cmd = 'mv index.html '+g.split(':')[1]+'.fa'
        gmo.run_command(cmd)

### main
parser = argparse.ArgumentParser(description="download fasta file of hogs (root hogs)")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="list of hogs")
args = parser.parse_args()

inFile = args.inFile

get_fasta_hogs(inFile)
print('End...')
