#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 16:51:35 2022

@author: irene
"""

import argparse
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM

def get_chromosomes(gffFile,chn):
    scaff = []
    for line in open(gffFile):
        line.strip()
        if line.startswith("#"):
            pass
        else:
            data = line.split('\t')
            if data[0] not in scaff:
                scaff.append(data[0])
    scaffolds = scaff[:chn]
    return scaffolds

## main
parser = argparse.ArgumentParser(description="create the gff file for mcscan")
parser.add_argument("-g", "--gffFile", dest="gffFile", required=True, help="gffFile")
parser.add_argument("-p", "--pepFile", dest="pepFile", required=True, help="pepFile")
parser.add_argument("-t", "--taxa", dest="taxa", default='no', help="mnemonic of the species if you want to add a tag")
parser.add_argument("-n", "--name", dest="name", required=True, help="three letters names of the species tag")
parser.add_argument("--chn", dest="chn", required=True, help="maximum number of chromosomes to be used")
args = parser.parse_args()

gffFile = args.gffFile
pepFile = args.pepFile
taxa = args.taxa
name = args.name
chn = int(args.chn)

pep = GM.load_sequences(pepFile)
scaffolds = get_chromosomes(gffFile,chn)
outfile1 = open(taxa+'.mcs.gff','w')
chromo = {}
i = 0
genes = set()
for line in open(gffFile):
    line = line.strip()
    if line.startswith("#"):
        pass
    else:
        data = line.split("\t")
        if data[2] == "gene":
            if data[0] in scaffolds:
                s, e = data[3], data[4]
                g = data[8].split("=")[1].split('.')[0]
                g+='_'+taxa
                genes.add(g)
                sc = data[0]
                if data[0] not in chromo:
                    i += 1
                    chromo[data[0]]=name+str(i)
                ch = name+str(i)
                string = ch + "\t" + g + "\t" + s + "\t" + e
                print(string, file=outfile1)
outfile1.close()

outfile2 = open(taxa+'.chr.conversion.txt','w')
for c in chromo:
    print(c+'\t'+chromo[c],file=outfile2)
outfile2.close()

outfile3 = open(taxa+'.mcs.pep','w')
for g in genes:
    GM.print_sequence(g,''.join(pep[g]),outfile3)
outfile3.close()
print('End...')