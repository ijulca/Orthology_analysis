#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 14:22:38 2020

@author: ijulca
"""
import argparse
import sys, glob, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM
import general_modules as gmo

def get_sp2num(genes):
    species = {}
    for g in genes:
        sp = g.split('-')[-1]
        if sp not in species:
            species[sp] = 0
        species[sp] += 1
    return species

def issingle(genes, nsp):
    species = get_sp2num(genes)
    toprint = False
    if len(species) >= nsp:
        i = 0
        for s in species:
            if species[s] == 1:
                i+=1
        if i >= nsp:
            toprint = True
    return toprint

def islow(genes, nsp):
    species = get_sp2num(genes)
    toprint = False
    if len(species) >= nsp:
        i,j=0,0
        for s in species:
            if species[s] == 1:
                i+=1
            elif species[s] == 2: ### up to 2 copy
                j+=1
        if (i+j) >= nsp:
            per = 0.2*len(species) ### 20% of species with 2 copies
            if j<=per:
                toprint = True
    return toprint
    

def get_orthogroups2single(orthoFile, num, gaps):
    print('loading single copy orthogroups...')
    orthologs = {}
    nsp = num - gaps
    for line in open(orthoFile): 
        line = line.strip() 
        data = line.split(' ')
        genes = data[1:]
        single = issingle(genes, nsp)
        if single == True:
            name = data[0].replace(':','')
            orthologs[name] = genes
    return orthologs
        
        
def get_orthogroups2low(orthoFile, orthologs, num, gaps):
    print('loading low copy orthogroups...')
    nsp = num - gaps
    for line in open(orthoFile): 
        line = line.strip() 
        data = line.split(' ')
        name = data[0].replace(':','')
        if name not in orthologs:
            genes = data[1:]
            tag = islow(genes, nsp)
            if tag == True:
                orthologs[name] = genes
    return orthologs

def get_proteins(files):
    prot = {}
    for f in files:
        seq = GM.load_sequences(f)
        prot.update(seq)
    return prot

def get_largest(genes, proteins):
    species = {}
    for g in genes:
        sp = g.split('-')[-1]
        if sp not in species:
            species[sp] = []
        species[sp].append(g+'++'+str(len(proteins[g])))
    new_genes = []
    for sp in species:
        if len(species[sp]) >1:
            pep = sorted(species[sp], key=lambda x:int(x.split('++')[1]))
            new_genes.append(pep[-1])
        else:
            new_genes.append(species[sp][0])
    if len(new_genes) != len(species.keys()):
        print('ERROR selecting duplicates...')
    return new_genes
    

def get_concat(orthologs, proteins, path, num, gaps):    
    for o in orthologs:
        outpath = path + o+'/'
        outname = outpath + o +'.seqs'
        gmo.create_folder(outpath)
        outfile = open(outname, 'w')
        genes = orthologs[o]
        tag = issingle(genes, num-gaps)
        if tag == False:
            genes = get_largest(genes, proteins)
        for g in genes:
            GM.print_sequence(g,''.join(proteins[g]), outfile)
        outfile.close()
        


### main
parser = argparse.ArgumentParser(description="get single or low copy orthogroups for concatenation")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="orthogroups file")
parser.add_argument("-p", "--pepFile", dest="pepFile", required=True, help="path where the pep files are stored, it will search for .fa")
parser.add_argument("-g", "--gaps", dest="gaps", default=0, help="allow gaps in number of species (1 = max one species missing). Default=0")
parser.add_argument("-e", "--extra", dest="extra", action='store_true', help="get extra-aligments by adding gene families with 2 copy genes")
args = parser.parse_args()

inFile = args.inFile
pepFiles = glob.glob(args.pepFile+'/*.fa')
gaps = int(args.gaps)
extra = args.extra

print('analysing '+str(len(pepFiles))+' species...')
proteins = get_proteins(pepFiles)
proteins = GM.remove_stopCodon(proteins)

path = './orthogroups_aligment/'
gmo.create_folder(path)

orthologs = get_orthogroups2single(inFile, len(pepFiles), gaps)
path = './orthogroups_aligment/single_copy/'
if os.path.exists(path):
    print('single copy orthologs fasta already exists...')
else:
    gmo.create_folder(path)
    get_concat(orthologs, proteins, path, len(pepFiles), gaps)

if len(orthologs) < 20:
    if extra == True:
        print('few single copy orthologs', len(orthologs))
        orthologs = get_orthogroups2low(inFile, orthologs, len(pepFiles), gaps)
        path = './orthogroup_aligment/low_copy/'
        gmo.create_folder(path)
        get_concat(orthologs, proteins, path, len(pepFiles), gaps)
    else:
        print('Warning, less than 20 single copy genes:', len(orthologs))
        print('Try --gaps or --extra options. Is better to try first --gaps option if was not used before')

print('End...')