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


def get_table(nodes,genes):
    table = {}
    for n in nodes:
        if n not in table:
            table[n] = 0
    for g in genes:
        n = gene2node[g]
        table[n]+=1
    return table

def get_string(string, nodes,table, tot):
    for n in nodes:
        string += '\t'+str(table[n]*100/tot) ## get percentage
    return string

def get_info_list(df,tag_list,tag_name,nodes,sp,outfile):
    df2 = df[df['FINAL LOCALIZATION'].isin(tag_list)]
    genes = set(df2['GeneID'].values)
    genes = [x+'-'+sp for x in genes]
    table = get_table(nodes,genes)
    string = sp + '\t' +tag_name
    string = get_string(string, nodes,table, len(genes))
    print(string,file=outfile)

def get_info_tag(df,tag,nodes,sp,outfile):
    df2 = df[df['FINAL LOCALIZATION'] == tag]
    genes = set(df2['GeneID'].values)
    genes = [x+'-'+sp for x in genes]
    table = get_table(nodes,genes)
    string = sp + '\t' +tag
    string = get_string(string, nodes,table, len(genes))
    print(string,file=outfile)

### main
taxaFile = '/home/ijulca/projects/Malaria/taxa.txt'
orthoFile = '/home/ijulca/projects/Malaria/proteomes/OrthoFinder/Results_Dec03/Orthogroups/Orthogroups.txt'
spTree = '/home/ijulca/projects/Malaria/species_tree/species_tree18.txt' 
pepPath = ''

infiles = glob.glob('/home/ijulca/projects/Malaria/Data_plasmodium/*')

### outputs:
outpath = '/home/ijulca/projects/Malaria/analysis/'
outOrtho2node = outpath+'ortho2node.txt'
outGene2node = outpath+'gene2node.txt'
treeNodes = outpath+'tree_percentageNodes.svg'

table1 = outpath+'Table_1.csv'
# ### create tree and ortho2node file


# species = OM.loadTaxa(taxaFile)
# orthogroups = OM.loadOrthofinder(orthoFile)

# t,node_names = OM.load_tree_nodes(spTree)
# perOrtho = OM.create_ortho2node(orthogroups, outOrtho2node, t, node_names)

# t = OM.change_leafName(t,species)
# OM.tree_nodes_orthoper(t, perOrtho, treeNodes)


# ### create gene to node
# ortho2node = OM.load_ortho2node(outOrtho2node)
# orthogroups = OM.loadOrthofinder(orthoFile)
# OM.create_gene2node(orthogroups, ortho2node, outGene2node)


##########
#proteins = OM.load_pepFromPath(pepPath)
gene2node = OM.load_ortho2node(outGene2node)
nodes = ['NODE_'+str(i) for i in range(1,18)]
plasmo = ['PLAF7','PLABA', '31271','PLAKH','1323249']
nodes+=plasmo

parasite = ['PARASITE','parasite','PPM']
paras = ['PARASITE','parasite']
exported = ["Cleft's","EXPORTED","GHOST","HCC","PV","PVM","PV or PVM","VESICLE"]

outfile = open(table1,'w')
print('species\tTag\t'+'\t'.join(nodes))
for f in infiles:
    sp = OM.get_prefix(f).split('_all')[0]
    df = pd.read_csv(f,sep='\t',header=0)
    get_info_list(df,parasite,'Parasite_all',nodes,sp,outfile)
    get_info_list(df,paras,'Parasite',nodes,sp,outfile)
    get_info_tag(df,'PPM',nodes,sp,outfile)
    get_info_list(df,exported,'Exported_all',nodes,sp,outfile)
    for tag in exported:
        get_info_tag(df,tag,nodes,sp,outfile)
    
outfile.close()
    