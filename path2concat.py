#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 17:02:45 2020

@author: ijulca
"""
import argparse
import glob
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM

def load_seq_sp(files, delimiter, pos):
    seqs = {}
    for f in files:
        family = f.split('/')[-1].split('.')[0]
        genes = {}
        for line in open(f):
            line = line.strip()
            if '>' in line:
                sp = line.split(delimiter)[pos]
                if '>' in sp:
                    sp = sp.replace('>','')
                if sp not in genes:
                    genes[sp] = []
                else:
                    print('ERROR...', 'duplicated genes for', sp)
            else:
                for b in line:
                    genes[sp].append(b)
        seqs[family] = genes
    return seqs

def get_species(seqs):
    species = set([])
    for s in seqs:
        for n in seqs[s]:
            species.add(n)
    return species

def concat_genes(files,delimiter,pos,outname):
    seqs = load_seq_sp(files, delimiter, pos)
    species = get_species(seqs)
    print('')
    
    

    


### main
parser = argparse.ArgumentParser(description="get the concatenation of aligned genes")
parser.add_argument("-p", "--path", dest="path", required=True, help="path to the aligned genes in fasta format")
args = parser.parse_args()

path = args.path

genes = glob.glob(path+'/*')