#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 17:36:38 2021

@author: ijulca
"""
import argparse, glob
import sys,os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import phylome_analysis as PA
import general_modules as gmo

def check_spider(f):
    folder = f.split('model')[0]
    log = folder+'/genetree.log'
    toprint = False
    if os.path.isfile(log) == True:
        with open(log, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
            if 'Date and Time:' in last_line:
                toprint = True
    return toprint


def create_allTrees(treeFiles,outFile):
    print('creating treeFile...')
    outfile = open(outFile, 'w')
    for treeFile in treeFiles:
        print(treeFile)
        totree = check_spider(treeFile.replace('treefile','log'))
        if totree == True:
            group = treeFile.split('/')[-2]
            tree = gmo.load_list(treeFile)
            if len(tree) == 1:
                print(group+'\t'+tree[0],file=outfile)
            else:
                print('ERROR...', treeFile)
    outfile.close()

def get_distFile(genetrees, outname):
    outfile = open(outname, 'w')
    for group in genetrees:
        tree = genetrees[group]
        dist = PA.get_average_branchLen(tree)
        string = group +'\t'+str(dist)
        print(string)
    return outfile

### main
parser = argparse.ArgumentParser(description="get the gene trees and analysis")
parser.add_argument("-p", "--path", dest="path", default='./Data/', help="path where is the data. Default=./Data/")
#parser.add_argument("-s", "--spTree", dest="spTree", required=True, help="rooted species tree")
args = parser.parse_args()

path = args.path+'/'
#spTreeFile = args.spTree

alltreeFile = 'alltree.txt'
ditFile = 'branch_lenght.txt'

treeFiles = glob.glob(path +'*/*/genetree.treefile')
print('Number of trees found:',len(treeFiles))

if os.path.isfile(alltreeFile) == False:
    create_allTrees(treeFiles,alltreeFile)

### analysis
# genetrees = PA.get_trees_from_file(alltreeFile)
# print('Number of tree to be analysed:', len(genetrees))
# get_distFile(genetrees, ditFile)

#spTree,spe2age = PA.load_species_tree(spTreeFile,'no')
#print(spe2age)
print('End...')
    
    
        
        