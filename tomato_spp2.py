#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 12:13:57 2023

@author: ijulca
"""
import glob
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM

path = "/home/ijulca/projects/tomato_project/"

blastFiles = glob.glob(path+"blast_results/*.blast")

genes = set()
for f in blastFiles:
    for line in open(f):
        line = line.strip()
        data = line.split("\t")
        e = float(data[10])
        if e <0.005:
            genes.add(data[1])
print("number of genes:", len(genes))

pepFiles = glob.glob(path+"protein_DB/*.fa")
outfile = open(path +"spp_gene.fa", "w")
for f in pepFiles:
    seqs = GM.load_sequences(f)
    for s in seqs:
        if s in genes:
            GM.print_sequence(s.replace(":","."),seqs[s],outfile)
outfile.close()