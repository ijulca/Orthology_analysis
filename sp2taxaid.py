#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 17:15:45 2021

@author: ijulca
"""
https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi

import argparse 
from Bio import Entrez
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import general_modules as gmo

Entrez.email = 'irene.julcac@gmail.com'

def get_taxaid(inFile):
    sp = 'Cycas circinalis'
    handle = Entrez.esearch(db="taxonomy", retmax=10, term=sp, idtype="acc")
    record = Entrez.read(handle)
    handle.close()
    print(record)
    
    # outfile = open(inFile+'.lineage_class', 'w')
    # for line in open(inFile):
    #     line = line.strip()
        
    #     handle = Entrez.efetch('taxonomy', id=tax, rettype='xml')
    
    # for tax in tax_ids:
    #     handle = Entrez.efetch('taxonomy', id=tax, rettype='xml')
    #     response = Entrez.read(handle)

    #     for entry in response:
    #         sci_name = entry.get('ScientificName')
    #         lineage_taxa = entry.get('Lineage').split(';')
    #         lineage_taxa = [x.strip() for x in lineage_taxa]
    #         print(tax, sci_name, '\t'.join(lineage_taxa), sep='\t',file=outfile)
    # outfile.close()

### main
# parser = argparse.ArgumentParser(description="get the taxID of a list of species")
# parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="list of species name")
# args = parser.parse_args()

# inFile = args.inFile
inFile = 1
get_taxaid(inFile)
print('End...')