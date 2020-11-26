#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 09:51:58 2020

@author: ijulca
"""
import argparse 
from Bio import Entrez
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import general_modules as gmo

Entrez.email = 'irene.julcac@gmail.com'

def get_lineage(inFile):
    outfile = open(inFile+'.lineage_class', 'w')
    tax_ids = set(gmo.load_list(inFile, sep='\t', remove='taxid'))
    for tax in tax_ids:
        handle = Entrez.efetch('taxonomy', id=tax, rettype='xml')
        response = Entrez.read(handle)

        for entry in response:
            sci_name = entry.get('ScientificName')
            lineage_taxa = entry.get('Lineage').split(';')
            lineage_taxa = [x.strip() for x in lineage_taxa]
            print(tax, sci_name, '\t'.join(lineage_taxa), sep='\t',file=outfile)
    outfile.close()

### main
parser = argparse.ArgumentParser(description="get the taxonomic classification of samples")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="list of taxaId")
args = parser.parse_args()

inFile = args.inFile

get_lineage(inFile)
print('End...')