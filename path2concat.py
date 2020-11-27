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

def concat_genes(files, delimiter, pos):
    seqs = {}
    for f in files:
        
    


### main
parser = argparse.ArgumentParser(description="get the concatenation of aligned genes")
parser.add_argument("-p", "--path", dest="path", required=True, help="path to the aligned genes in fasta format")
args = parser.parse_args()

path = args.path

genes = glob.glob(path+'/*')