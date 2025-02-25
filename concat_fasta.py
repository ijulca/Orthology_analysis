#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 14:26:21 2025

@author: ijulcach
"""

import argparse
import sys, os
import glob
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import general_modules as GM


def get_species(fastaFiles, tag):
    species = []
    for file in fastaFiles:
        seq = GM.load_sequences(file)
        if tag == 'oma':
            names = [x[:5] for x in list(seq.keys())]
        elif tag == 'no':
            names = list(seq.keys())
        else:
            print('ERROR..., wrong fasta name format')
        species += names
    return list(set(species))

def get_len_alg(seq):
    names = list(seq.keys())
    num = len(seq[names[0]])
    return num

def read_alignment(fastaFiles, tag):
    species_codes = get_species(fastaFiles, tag)
    print('Total number of species in the alignments:', len(species_codes))
    algs = {x:'' for x in species_codes}
    for file in fastaFiles:
        seq = GM.load_sequences(file)
        if tag == 'oma':
            seq = {x[:5]:seq[x] for x in seq}
        if len(seq) != len(species_codes):
            len_alg = get_len_alg(seq)
        for s in species_codes:
            if s in algs[s]:
                algs[s] += ''.join(seq[s])
            else:
                algs[s] += '-'*len_alg
    return algs
    

def concatenate_alignments_from_path(base_path,outfile,tag):
    if base_path[-1:] != "/":
        base_path = base_path+"/"
    fastaFiles = glob.glob(f'{base_path}*')
    algs = read_alignment(fastaFiles, tag)
    with open(outfile, 'w') as fout:
        for sp in algs:
            GM.print_sequence(sp, algs[sp], fout)


##############
#### main ####
##############

parser = argparse.ArgumentParser(description="concat fasta files")
parser.add_argument("-i", "--inPath", dest="inPath", required=True, help="path where the fasta files are")
parser.add_argument("-o", "--outname", dest="outname", required=True, help="outname")
parser.add_argument("-t", "--tag", dest="tag", default='oma', help="tag of the format. default='oma', no=no change")
args = parser.parse_args()

if __name__=="__main__":
    print('Start...')
    concatenate_alignments_from_path(args.inPath,args.outname,args.tag)
    print('End...')