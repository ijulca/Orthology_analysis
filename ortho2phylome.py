#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 21:07:29 2021

@author: ijulca
"""
import argparse
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM
import orthology_modules as OM
import general_modules as gmo



### main
parser = argparse.ArgumentParser(description="get formated fasta files for analysis")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="orthogroups file")
parser.add_argument("-p", "--pepFile", dest="pepFile", required=True, help="db of the proteins")
parser.add_argument("-c", "--cdsFile", dest="cdsFile", default='', help="db of the cds (optional)")
parser.add_argument("-s", "--size", dest="size", default='2', help="minimum size to be included. Default=2")
parser.add_argument("-o", "--path", dest="path", default='./Data/', help="path where to create the folders")
args = parser.parse_args()

orthoFile = args.inFile
pepFile = args.pepFile
cdsFile = args.cdsFile
size = int(args.size)
outpath = args.path+'/'

ortho2pep = OM.loadOrthofinder(orthoFile)
ortho2pep = OM.flt_orthogroups(ortho2pep, size)
pepDB = GM.load_sequences(pepFile)
if cdsFile != '':
    cdsDB = GM.load_sequences(cdsFile)

families = list(ortho2pep.keys())
print('Number of orthogroups that will be analysed:', len(families))
r = len(families)+1

gmo.create_folder(outpath)
for i in range(1,r,1000):
    
    z = i + 1000
    if z >r:
        z = r
    outdir1 = outpath+str(i)+'-'+str(z)+'/'
    gmo.create_folder(outdir1)
    for j in range(i,z):
        group = families[j]
        genes = ortho2pep[group]
        outdir2 = outpath+group+'/'
        gmo.create_folder(outdir2)
        outpep = open(outdir2+group+'.fa', 'w')
        if cdsFile != '':
            outcds = open(outdir2+group+'.cds', 'w')
        for g in genes:
            GM.print_sequence(g,''.join(pepDB[g]),outpep)
            if cdsFile != '':
                GM.print_sequence(g,''.join(cdsDB[g]),outcds)
        outpep.close()
        if cdsFile !='':
            outcds.close()

print('End...')