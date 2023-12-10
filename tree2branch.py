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
import ete4
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns



def check_spider(log):
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
        totree = check_spider(treeFile.replace('treefile','log'))
        if totree == True:
            group = treeFile.split('/')[-2]
            tree = gmo.load_list(treeFile)
            if len(tree) == 1:
                print(group+'\t'+tree[0],file=outfile)
            else:
                print('ERROR...', treeFile)
        else:
            print('ERROR...bad tree', treeFile)
    outfile.close()

def get_distFile(genetrees, outname):
    outfile = open(outname, 'w')
    i = 0
    for group in genetrees:
        tree = genetrees[group]
        dist0 = PA.get_average_branchLen(tree,0)
        dist1 = PA.get_average_branchLen(tree,95)
        if len(dist0) == 0:
            dist0 = ['NA']
        if len(dist1) == 0:
            dist1 = ['NA']
        dist0 = [str(x) for x in dist0]
        dist1 = [str(x) for x in dist1]
        string = group +'\t'+str(';'.join(dist0))+'\t'+str(';'.join(dist1))
        print(string,file=outfile)
    print(i)
    return outfile


### main
# parser = argparse.ArgumentParser(description="get the gene trees and analysis")
# parser.add_argument("-p", "--path", dest="path", default='./Data/', help="path where is the data. Default=./Data/")
# #parser.add_argument("-s", "--spTree", dest="spTree", required=True, help="rooted species tree")
# args = parser.parse_args()

# path = args.path+'/'
# #spTreeFile = args.spTree

# #####outputs
# alltreeFile = 'alltree.txt'
# ditFile = 'branch_lenght.txt'
# ####


# if os.path.isfile(alltreeFile) == False:
#     treeFiles = glob.glob(path +'*/*/genetree.treefile')
#     print('Number of trees found:',len(treeFiles))
#     create_allTrees(treeFiles,alltreeFile)

# ### analysis
# genetrees = PA.get_trees_from_file(alltreeFile)
# print('Number of tree to be analysed:', len(genetrees))
# get_distFile(genetrees, ditFile)

# #spTree,spe2age = PA.load_species_tree(spTreeFile,'no')
# #print(spe2age)
# print('End...')
    
    
#######################
#### PANTHER TREES #### 
#######################

def get_tree_panther(inFile):
    i = 0
    for line in open(inFile):
        line = line.strip()
        i += 1
        if i == 1:
            tree = line
    return tree

pathPlots = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/results/plots/'
pathTables = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/results/tables/'
treepath = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/panther-18.0/trees/'
speciesFile = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/panther-18.0/species_tree.nhx'
expPlantdata = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/Plant_expression/Samples/'

#### Get number of duplications
# treeFiles = glob.glob(treepath+'*')

# outfile = open(pathTables+'tree_duplications.tsv2','w')
# for f in treeFiles:
#     name = f.split('/')[-1].split('.')[0]
#     evolevnts = {'UNK':'U','0>1':'S','1>0':'D','0>0':'H'}
#     tree = get_tree_panther(f)       
#     t = ete4.Tree(tree)
#     dup = 0
#     for node in t.traverse():
#         if node.is_root:
#             pass
#         else:
#             props = node.props
#             if 'Ev' in props:
#                 ev = evolevnts[props['Ev']]
#                 if ev == 'D':
#                     dup += 1
#     string = name +'\t'+ str(dup)
#     print(string,file=outfile)
# outfile.close()

##### family-specific rate ### Supplementary Figure 1
# inFile = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/results/tables/fitted_family_info.tsv'
# outfigure1 = pathPlots +'histogram_familyrates.svg'
# df = pd.read_csv(inFile, sep='\t', header=0)
# ax = sns.histplot(df,x='rate', bins=200)
# ax.set_title('Histogram of family-specific factor')
# ax.set_xlabel('Family-specific rate')
# ylims = ax.get_ylim()
# plt.vlines(1, ylims[0]-1000, ylims[1]+1000, color='red', linestyle='--')#, ymin=0, ymax=)\n",
# ax.set_ylim(*ylims)
# plt.savefig(outfigure1, bbox_inches='tight')
# plt.show()
# print(df['rate'].mean())


##### Species tree

# lineages = {'Viridiplantae':[], 'Archaea':[], 'Eubacteria':[], 'Amoebozoa':[], 'Fungi':[],
#             'Deuterostomia':[], 'Protostomia':[], 'Alveolata-Stramenopiles':[],
#             'Excavates':[], 'NEMVE':['NEMVE'], 'TRIAD':['TRIAD'],'MONBE':['MONBE']}

# t = ete4.Tree(open(speciesFile))

# taken = set()
# for node in t.traverse('postorder'):
#     props = node.props
#     if 'S' in props:
#         lin = props['S']
#         # print(lin)
#         if lin in lineages:
#             taken.add(node)
#             for leaf in node:
#                 lineages[lin].append(leaf.name)
#                 if leaf.name not in taken:
#                     taken.add(leaf.name)
#                 else:
#                     print('Errror...', leaf.name)

# print(lineages)
# i = 0
# outfile = open(pathTables+'major_lineages_sp.txt','w')
# for e in lineages:
#     string = e+'\t'+str(len(lineages[e]))+'\t'+'; '.join(lineages[e])
#     print(string,file=outfile)
#     i += len(lineages[e])
# print(i)
# outfile.close()


# treeFiles = glob.glob(treepath+'*')

# outfile = open(pathTables+'plants_geneNames.txt','w')
# for f in treeFiles:
#     i = 0
#     for line in open(f):
#         line = line.strip()
#         i += 1
#         if i >1:
#             data = line.split(':')[1].split('|')
#             data = [x.replace(';','') for x in data]
#             if data[0] in lineages['Viridiplantae']:
#                 print('\t'.join(data), file=outfile)
# outfile.close()
    
#### Getting the names of the gene expression matrices and panther
expPlantdata = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/Plant_expression/new_samples/'
key = 'VITVI'
genes = set()
for line in open(expPlantdata+key+'/'+key+'.cds.fa'):
    line = line.strip()
    if '>' in line:
        genes.add(line.split('>')[1])
    # data = line.split('\t')
    # if data[0] != 'gene':
    #     genes.add(data[0].split('_')[0])

# dades = {}
# for line in open(expPlantdata+key+'/AGI2uniprot-Jul2023.txt'):
#     line = line.strip()
#     data = line.split('\t')
#     dades[data[1].split('-')[0]] =data[0].split('.')[0]

table = {}
for line in open(expPlantdata+key+'/'+key+'.genesNames.txt'):
    line = line.strip()
    data = line.split('\t')
    name = '|'.join(data)
    if 'EnsemblGenome' in data[1]:
        n = data[1].split('=')[1]
    # n = data[2].split('=')[1]
    # n = data[1].split('=')[1].replace('CISIN_','orange1.').replace('mg','m.g')
        table[n] = name

# outfile = open(expPlantdata+key+'/'+key+'.conversion.panther.txt','w')
for e in table:
    if e not in genes:
        print(e, table[e])
#         string = e+'\t'+table[e]
#         print(string,file=outfile)
# outfile.close()
        
# print(genes)   
        
       
        
       