#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 16:38:29 2020

@author: ijulca
"""
import glob
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import orthology_modules as OM
import pandas as pd



### main
taxaFile = '/home/ijulca/projects/Malaria/taxa.txt'
orthoFile = '/home/ijulca/projects/Malaria/proteomes/OrthoFinder/Results_Dec03/Orthogroups/Orthogroups.txt'
spTree = '/home/ijulca/projects/Malaria/species_tree/species_tree18.txt' 
pepPath = ''

### outputs:
outpath = '/home/ijulca/projects/Malaria/analysis/'
outOrtho2node = outpath+'ortho2node.txt'
outGene2node = outpath+'gene2node.txt'
treeNodes = outpath+'tree_percentageNodes.svg'

infiles = glob.glob('/home/ijulca/projects/Malaria/Data_plasmodium/*')

# species = OM.loadTaxa(taxaFile)
# orthogroups = OM.loadOrthofinder(orthoFile)

# t,node_names = OM.load_tree_nodes(spTree)
# perOrtho = OM.ortho2node(orthogroups, outOrtho2node, t, node_names)

# t = OM.change_leafName(t,species)
# OM.tree_nodes_orthoper(t, perOrtho, treeNodes)


### create gene to node
ortho2node = OM.load_ortho2node(outOrtho2node)
orthogroups = OM.loadOrthofinder(orthoFile)


##########
# #proteins = OM.load_pepFromPath(pepPath)



# for f in infiles:
#     df = pd.read_csv(f,sep='\t',header=0)