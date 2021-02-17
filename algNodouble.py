#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 12:08:41 2021

@author: ijulca
"""

import argparse
import sys, os
import glob
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM
import random

def get_doubles(data):
    species = {}
    for g in data:
        sp = g.split('-')[-1]
        if sp not in species:
            species[sp]=set([])
        species[sp].add(g)
    seqs = []
    i = 0
    for sp in species:
        if len(species[sp]) !=1:
            i+=1
            s = random.choice(list(species[sp]))
            seqs += [x for x in species[sp] if x != s]
    return seqs,i
            


### main
parser = argparse.ArgumentParser(description="remove duplicates from aligment")
parser.add_argument("-p", "--path", dest="path", required=True, help="path to the aligned genes in fasta format")
args = parser.parse_args()

path = args.path

files = glob.glob(path+'/*')

for f in files:
    seq = GM.load_sequences(f)
    doubles,sp = get_doubles(seq)
    if sp == 1:
        if len(doubles) != 0:
            print(f)
            outfile = open(f+'.2','w')
            for s in seq:
                if s not in doubles:
                    GM.print_sequence(s,''.join(seq[s]),outfile)
            outfile.close()

print('End...')