#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 14:47:31 2021

@author: ijulca
"""
import argparse
import glob
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM
import general_modules as gmo

#trimal = '/home/irene.julca/anaconda/envs/phylo/bin/trimal'
trimal = 'trimal'

def get_genes(inFile):
    genes = []
    for line in open(inFile):
        line = line.strip()
        if '>' in line:
            g = line.split('>')[1]
            genes.append(g)
    return genes


### main
parser = argparse.ArgumentParser(description="change order of aligments")
parser.add_argument("-p", "--path", dest="path", required=True, help="path where the clean_cds aligment is")
args = parser.parse_args()

path = args.path+'/'
alg = glob.glob(path+'*.alg.metalig')
algFile = alg[0]

tempfasta = path +'test.fasta'
cmd = trimal +' -in '+algFile +' -out '+tempfasta+' -fasta'
gmo.run_command(cmd)
genes = get_genes(tempfasta)

pathFile = glob.glob(path +'*.alg.paths')
pathFile = pathFile[0]

files = gmo.load_list(pathFile)
print(files)
for f in files:
    print(f)
    seqs = GM.load_sequences(f)
    outfile = open(f,'w')
    for g in genes:
        GM.print_sequence(g,''.join(seqs[g]), outfile)
    outfile.close()

cmd = 'rm '+tempfasta
gmo.run_command(cmd)
print('End')

