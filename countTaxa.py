#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 10:46:24 2020

@author: ijulca
"""
import argparse 
from Bio import Entrez

# Entrez.email = 'irene.julcac@gmail.com'

# ### main
# # parser = argparse.ArgumentParser(description="count the number of species per taxonomica level")
# # parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="classification file")
# # parser.add_argument("-t", "--taxa", dest="taxa", required=True, help="higher taxonomic level. Ex: Viridiplantae")
# # args = parser.parse_args()

# # inFile = args.inFile
# # taxa = args.taxa

# inFile = '/home/ijulca/Archeoplastida_expression/taxa.txt.lineage_class'
# taxa = "Viridiplantae"

# handle = Entrez.esearch(db="Taxonomy", term=taxa)
# record = Entrez.read(handle)
# idrec = record["IdList"][0]
# handle = Entrez.efetch(db="Taxonomy", id=idrec, retmode="xml")
# records = Entrez.read(handle)
# rank = records[0]['Rank']
# print(records[0]["Division"])

from ete3 import NCBITaxa

ncbi = NCBITaxa()
ncbi.update_taxonomy_database()
print(ncbi)