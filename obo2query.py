#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 14:16:50 2024

@author: ijulcach
"""

import argparse

def get_list(inFile):
    dades = set()
    for line in open(inFile):
        line = line.strip()
        dades.add(line)
    return dades

def get_uberon_info(dbFile):
    toprint = False
    table = {}
    general = set()
    for line in open(dbFile):
        line = line.strip()
        if line.startswith('id:'):
            if line.startswith('id: UBERON:'):
                term = line.split(' ')[1].strip()
                toprint = True
                table[term] = {'n':'none','p':[],'i':[]}
            else:
                toprint = False
        else:
            if toprint == True:
                if line.startswith('name'):
                    name = line.split('name:')[1]
                    table[term]['n'] = name
                elif line.startswith('relationship: part_of'):
                    part = line.split('part_of')[1].strip()
                    if part.startswith('UBERON'):
                        table[term]['p'].append(part)
                elif line.startswith('is_a:'):
                    part = line.split('is_a:')[1].strip()
                    if part.startswith('UBERON'):
                        table[term]['i'].append(part)
                        general.add(part)
    return table

def get_all_terms(inFile, dbFile):
    table = get_uberon_info(dbFile)
    data = get_list(inFile)
    for u in data:
        if u in table:
            print(table[u])
   
        
    
### main
# parser = argparse.ArgumentParser(description="get the parents of the uberon terms")
# parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="list of uberons")
# args = parser.parse_args()

# inFile = args.inFile

dbFile = '/home/ijulcach/projects/ldo_project/data/animal_entity/uberon-full.obo'
inFile = '/home/ijulcach/projects/ldo_project/data/animal_entity/organ_ube.txt'

get_all_terms(inFile, dbFile)
