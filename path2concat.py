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
    species = set([])
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
                species.add(sp)
            else:
                for b in line:
                    genes[sp].append(b)
        seqs[family] = genes
    return seqs,species

def concat_genes(files,delimiter,pos,outname):
    seqs,species = load_seq_sp(files, delimiter, pos)
    print('concatenation of ', len(seqs),' genes and ', len(species), ' species...')
    species =list(species)
    print(species)
    concat = {}
    for g in seqs:
        for s in species:
            if s not in concat:
                concat[s] = []
            concat[s] += seqs[g][s]
    outfile = open(outname, 'w')
    for sp in concat:
        GM.print_sequence(sp,''.join(concat[sp]),outfile)
    outfile.close()
            

### main
parser = argparse.ArgumentParser(description="get the concatenation of aligned genes")
parser.add_argument("-p", "--path", dest="path", required=True, help="path to the aligned genes in fasta format")
parser.add_argument("-d", "--delimiter", dest="delimiter", default='-', help="delimiter to get the name. Default='-'")
parser.add_argument("-n", "--number", dest="number", default='-1', help="position number of the species name starting from 0. Defalut -1")
parser.add_argument("-o", "--outFile", dest="outFile", required=True, help="outfile name")
args = parser.parse_args()

path = args.path
d = args.delimiter
num = int(float(args.number))
outname = args.outFile

files = glob.glob(path+'/*')

concat_genes(files,d,num,outname)
print('End...')
