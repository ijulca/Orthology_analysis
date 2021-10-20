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

#handle = Entrez.esearch(db="taxonomy", term="2759")
search = Entrez.efetch(id = "2759", db = "taxonomy", retmode = "xml")
record = Entrez.read(search)
print(record)


# search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
# esearch -db sra -query "PRJNA517527" | efetch -format runinfo | grep -v "Run" | cut -d ',' -f1
# efetch -db taxonomy -id 9606,7227,10090 -format xml

## search in the SRA web page: txid2759[Organism:exp] 

print('End...')







