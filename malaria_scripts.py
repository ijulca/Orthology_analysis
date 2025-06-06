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
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import random
import numpy as np


def get_table(nodes,genes, gene2node):
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

def get_info_list(df,tag_list,tag_name,nodes,gene2node,sp,outfile):
    df2 = df[df['FINAL LOCALIZATION'].isin(tag_list)]
    genes = set(df2['GeneID'].values) #### new list of genes, based on the localization
    genes = [x+'-'+sp for x in genes]
    table = get_table(nodes,genes, gene2node)
    string = sp + '\t' +tag_name
    string = get_string(string, nodes,table, len(genes))
    print(string,file=outfile)

def get_info_tag(df,tag,nodes,gene2node,sp,outfile):
    df2 = df[df['FINAL LOCALIZATION'] == tag]
    genes = set(df2['GeneID'].values)
    genes = [x+'-'+sp for x in genes]
    table = get_table(nodes,genes, gene2node)
    string = sp + '\t' +tag
    string = get_string(string, nodes,table, len(genes))
    print(string,file=outfile)
    
def get_info_listTot(df,tag_list,tag_name,nodes,gene2node,sp,outfile,tot):
    df2 = df[df['FINAL LOCALIZATION'].isin(tag_list)]
    genes = set(df2['GeneID'].values)
    genes = [x+'-'+sp for x in genes]
    table = get_table(nodes,genes, gene2node)
    string = sp + '\t' +tag_name
    string = get_string(string, nodes,table, tot)
    print(string,file=outfile)

def get_info_tagTot(df,tag,nodes,gene2node,sp,outfile,tot):
    df2 = df[df['FINAL LOCALIZATION'] == tag]
    genes = set(df2['GeneID'].values)
    genes = [x+'-'+sp for x in genes]
    table = get_table(nodes,genes, gene2node)
    string = sp + '\t' +tag
    string = get_string(string, nodes,table, tot)
    print(string,file=outfile)

def get_sp2gene2loc(files):
    table = {}
    for f in files:
        sp = OM.get_prefix(f).split('_all')[0]
        table[sp] = {}
        for line in open(f):
            line = line.strip()
            data = line.split('\t')
            if data[0] != 'GeneID':
                if len(data) == 8:
                    key = data[0]+'-'+sp
                    if data[7] == 'PARASITE' or data[7] == 'parasite':
                        data[7] = 'PARASITE'
                    table[sp][key] = data[7] ## get localization
                else:
                    pass ## not include the genes without localization
                    
        print(sp,len(table[sp]))
    return table

def get_sp2gene2locEnrich(files):
    table = {}
    for f in files:
        sp = OM.get_prefix(f).split('_all')[0]
        table[sp] = {}
        for line in open(f):
            line = line.strip()
            data = line.split('\t')
            if data[0] != 'GeneID':
                if len(data) == 8:
                    key = data[0]+'-'+sp
                    if data[7] != 'EXPORTED' and data[7] != 'PV or PVM': ### remove this loc because are ambiguos
                        if data[7] == 'PARASITE' or data[7] == 'parasite':
                            data[7] = 'PARASITE'
                        elif data[7] == 'VESICLE' or data[7] == "Cleft's":
                            data[7] = "VESICLE+Cleft's"
                        table[sp][key] = data[7] ## get localization
                else:
                    pass ## not include the genes without localization
                    
    return table

def get_sp2gene2node(infiles):
    table = {}
    for f in infiles:
        sp = OM.get_prefix(f).split('_')[0]
        table[sp] = {}
        for line in open(f):
            line = line.strip()
            data = line.split('\t')
            if data[0] != 'GeneID':
                key = data[0] +'-'+sp
                table[sp][key] = data[-1]
    return table


#### enrichment
def load_data_plot(inFile):
    table = {}
    nodes = set([])
    for line in open(inFile):
        line = line.strip()
        data = line.split("\t")
        if data[0] != "Species":
            sp = data[0]
            org = data[1]
            node = data[2]
            nodes.add(node)
            if data[4] == 'enrichment':
                if sp not in table:
                    table[sp] = {}
                if node not in table[sp]:
                    table[sp][node] = {}
                table[sp][node][org]=float(data[3])
    return table, nodes

def get_enrichNodePlot(inFile, sp2nodes, outname, samples):
    table,nodes = load_data_plot(inFile)
    new_nodes = ['NODE_1', 'NODE_2', 'NODE_3', 'NODE_4', 'NODE_5', 'NODE_6', 'NODE_7', 'NODE_8', 'NODE_9',
                 'NODE_10', 'NODE_11', 'NODE_12', 'NODE_13', 'NODE_14', 'NODE_15', 'NODE_16', 'NODE_17', 
                 'PLAF7','PLAKH','31271','1323249', 'PLABA']
    #print(nodes)
    for sp in table:
        new_table = {}
        for n in new_nodes:
            if n in table[sp]:
                if n not in new_table:
                    new_table[n] = []
                for sam in samples:
                    if sam in table[sp][n]:
                        new_table[n].append(table[sp][n][sam])
                    else:
                        new_table[n].append(10)
        df = pd.DataFrame(new_table, index=samples)
        plt.clf()
        fig, ax = plt.subplots(figsize=(10,5))
        _ = sns.heatmap(df, cmap='Blues_r', mask=df>0.05, xticklabels=True, vmax=0.06)
        plt.savefig(outname+'_'+sp+'.svg',bbox_inches = "tight")

def get_sp2gene2node2(genenode):
    table = {}
    for line in open(genenode):
        line = line.strip()
        data = line.split('\t')
        sp = data[0].split('-')[-1]
        if sp not in table:
            table[sp] = {}
        table[sp][data[0]] = data[1]
    return table

def transform_table_dot(inFile):
    outfile = open(inFile.split('.')[0]+'_transform.csv','w')
    for line in open(inFile):
        line = line.strip()
        data = line.split('\t')
        if data[0] == 'Node':
            print(line, file=outfile)
        else:
            values = [float(x) for x in data[1:]]
            string = data[0]
            for v in values:
                if v >=0.8:
                    v = 1
                elif v == 0:
                    v = v
                elif v>=0.6 and v<0.8:
                    v=0.9
                elif v>=0.4 and v<0.6:
                    v = 0.8
                elif v>=0.3 and v<0.4:
                    v=0.7
                elif v>=0.2 and v<0.3:
                    v = 0.6
                elif v>=0.1 and v<0.2:
                    v = 0.5
                elif v>=0.06 and v<0.1:
                    v = 0.4
                elif v>=0.04 and v<0.06:
                    v = 0.3
                elif v>=0.02 and v<0.04:
                    v = 0.2
                elif v<0.02:
                    v = 0.1
                else:
                    print('ERROR', v)
                string += '\t'+str(v)
            print(string,file=outfile)
    outfile.close()
    

### main
taxaFile = '/home/ijulca/projects/Malaria/taxa.txt'
orthoFile = '/home/ijulca/projects/Malaria/proteomes/OrthoFinder/Results_Apr26/Orthogroups/Orthogroups.txt'
#spTree = '/home/ijulca/projects/Malaria/species_tree/species_tree18.txt' 
#spTree = '/home/ijulca/projects/Malaria/species_tree/species_tree18HM.txt'
spTree = '/home/ijulca/projects/Malaria/species_tree/species_tree19.txt'
pepPath = ''

infiles = glob.glob('/home/ijulca/projects/Malaria/Data_plasmodium/*')

nodeFiles = glob.glob('/home/ijulca/projects/Malaria/analysis/*_nodes.csv')

### outputs:
outpath = '/home/ijulca/projects/Malaria/analysis/'
outOrtho2node = outpath+'ortho2node'
outGene2node = outpath+'gene2node.txt'
treeNodes = outpath+'tree_percentageNodes.svg'

table1 = outpath+'Table_1.csv'
table2 = outpath+'Table_2.csv'

tablex2 = outpath+'Table_X2.csv'
tablex3 = outpath+'Table_X3.csv' ## table summary including number of genes, percentage of genes, and p-val
### enrichment
table3 = outpath+'Table_3.csv' ## percentage of genes per node
table4 = outpath+'Table_4.csv' ## p-values per node
outenrichplot = outpath+'Figure_dotsEnrich.svg' ## p-values per node
###

table5 = outpath +'Table_5.csv'
table6 = outpath +'Table_6.csv'



#### 1) create tree and ortho2node file

# species = OM.loadTaxa(taxaFile)
# orthogroups = OM.loadOrthofinder(orthoFile)

# t,node_names = OM.load_tree_nodes(spTree)
# perOrtho = OM.create_ortho2node(orthogroups, outOrtho2node, t, node_names)

# t = OM.change_leafName(t,species)
# OM.tree_nodes_orthoper(t, perOrtho, treeNodes)


# ### 2) create gene to node
# ortho2node = OM.load_ortho2node(outOrtho2node+'.txt')
# orthogroups = OM.loadOrthofinder(orthoFile)
# OM.create_gene2node(orthogroups, ortho2node, outGene2node)


########## 3) Table 1

# gene2node = OM.load_ortho2node(outGene2node)
# nodes = ['NODE_'+str(i) for i in range(1,19)] ### check for the number of nodes (+1)
# plasmo = ['PLAF7','PLABA', '31271','PLAKH','1323249']
# nodes+=plasmo

# parasite = ['PARASITE','parasite','PPM']
# paras = ['PARASITE','parasite']
# exported = ["Cleft's","EXPORTED","GHOST","HCC","PV","PVM","PV or PVM","VESICLE"]

# outfile = open(table1,'w')
# print('species\tTag\t'+'\t'.join(nodes),file=outfile)
# for f in infiles:
#     sp = OM.get_prefix(f).split('_all')[0]
#     df = pd.read_csv(f,sep='\t',header=0)
#     get_info_list(df,parasite,'Parasite_all',nodes,gene2node,sp,outfile)
#     get_info_list(df,paras,'Parasite',nodes,gene2node,sp,outfile)
#     get_info_tag(df,'PPM',nodes,gene2node,sp,outfile)
#     get_info_list(df,exported,'Exported_all',nodes,gene2node,sp,outfile)
#     for tag in exported:
#         get_info_tag(df,tag,nodes,gene2node,sp,outfile)
# outfile.close()

# ####### 4) Table 2

# outfile = open(table2,'w')
# print('species\tTag\t'+'\t'.join(nodes),file=outfile)
# for f in infiles:
#     sp = OM.get_prefix(f).split('_all')[0]
#     df = pd.read_csv(f,sep='\t',header=0)
#     tot = len(set(df['GeneID'].values)) ### Omar Table
#     get_info_listTot(df,parasite,'Parasite_all',nodes,gene2node,sp,outfile,tot)
#     get_info_listTot(df,paras,'Parasite',nodes,gene2node,sp,outfile,tot)
#     get_info_tagTot(df,'PPM',nodes,gene2node,sp,outfile,tot)
#     get_info_listTot(df,exported,'Exported_all',nodes,gene2node,sp,outfile,tot)
#     for tag in exported:
#         get_info_tagTot(df,tag,nodes,gene2node,sp,outfile,tot)
    
# outfile.close()

###  5) create tables of each file with the node feature

# for f in infiles:
#     sp = OM.get_prefix(f).split('_all')[0]
#     outfile = open(outpath+sp+'_nodes.csv', 'w')
#     for line in open(f):
#         line = line.strip()
#         data = line.split('\t')
#         if data[0] == 'GeneID':
#             line += '\t'+'phylostrata'
#         else:
#             if len(data) <8:
#                 line = line +'\t-'*(8-len(data))
#             key = data[0]+'-'+sp
#             line+= '\t'+gene2node[key]
#         print(line,file=outfile)
#     outfile.close()


# ##### 6) Calculate enrichment per loc per node:

species = OM.loadTaxa(taxaFile)
    
sp2gene2node = get_sp2gene2node(nodeFiles) ### genes only from Omar table. Species, genes to node
###sp2gene2node = get_sp2gene2node2(outGene2node) ### all genes of the species

nodes = ['NODE_1', 'NODE_3','NODE_4','NODE_6','NODE_7','NODE_8','NODE_9','NODE_11','NODE_13','NODE_15','NODE_16','NODE_17','NODE_18']
sp2nodes = {'PLAF7':nodes+['PLAF7'],'PLAKH':nodes+['NODE_15','NODE_16','PLAKH'],'31271':nodes+['NODE_15','NODE_17','31271'], 
            '1323249':nodes+['NODE_15','NODE_17','NODE_18','1323249'],'PLABA':nodes+['NODE_15','NODE_17','NODE_18','PLABA']}

node2phylo = {'NODE_1':'Eukaryota','NODE_3':'SAR','NODE_4':'Alveolata','NODE_6':'Myzozoa','NODE_7':'Myzozoa',
              'NODE_8':'Apicomplexa','NODE_9':'Coccidiomorpha','NODE_11':'Haematozoa','NODE_13':'Plasmodium',
              'NODE_15':'Non-falciparum','NODE_16':'P. knowlesi/P. vivax','NODE_17':'Rodent Plasmodium','NODE_18':'P. berghei/P. yoelii',
              'PLAF7':'P. falciparum','PLAKH':'P. knowlesi','1323249':'P. yoelii','PLABA':'P. berghei','31271':'P. chabaudi'}

x=10000  # number of replicates for the enrichment analysis

parasite = ['PARASITE','PPM']
# exported = ["Cleft's","EXPORTED","GHOST","HCC","PV","PVM","PV or PVM","VESICLE"]
peripheral = ["PV","PVM"]
exported = ["HCC","VESICLE+Cleft's","GHOST"]


sp2gene2loc = get_sp2gene2locEnrich(infiles) ## Species, genes to loc

sp2loc2gene = {}
for sp in sp2gene2loc:
    sp2loc2gene[sp] = {}
    for g in sp2gene2loc[sp]:
        loc = sp2gene2loc[sp][g]
        if loc not in sp2loc2gene[sp]:
            sp2loc2gene[sp][loc] = set([])
        sp2loc2gene[sp][loc].add(g)

outfile = open(tablex2,'w')
print('Species\tGene\tSubcellular fraction\tPhylostratum',file=outfile)
dades = {}
for sp in sp2loc2gene:
    ##spgenes = list(sp2gene2node[sp].keys()) ### all proteome
    spgenes = list(sp2gene2loc[sp].keys()) ## only genes with localization
    ##spgenes = list(sp2loc2gene[sp]['PARASITE'])+list(sp2loc2gene[sp]['PPM']) ### using this loc as population
    print(sp,len(spgenes),sep='\t')
    dades[sp] = {}
    for loc in sp2loc2gene[sp]:
        dades[sp][loc] = {}
        genes = sp2loc2gene[sp][loc]
        for g in genes:
            string = 'P. '+species[sp].split(' ')[1]+'\t'+g.split('-')[0]+'\t'+loc+'\t'+node2phylo[sp2gene2node[sp][g]]
            print(string,file=outfile)
        #print(sp,loc,len(genes), sep='\t')
        nodes2genes = OM.get_spe_node(genes,sp2gene2node[sp])
        for n in nodes2genes:
            dades[sp][loc][n] = [len(nodes2genes[n])/len(genes),0.1,0.1,len(nodes2genes[n])]  ### len(genes) percentage relative to the location, len(spgenes) ## relative to species
        values = {}
        for i in range(0, x): 
            new_genes = random.sample(spgenes, len(genes))  
            new_nodes2genes = OM.get_spe_node(new_genes,sp2gene2node[sp])
            for n in nodes2genes:
                if n in new_nodes2genes:
                    if n not in values:
                        values[n] = 0
                    obs = len(nodes2genes[n])
                    calc = len(new_nodes2genes[n])
                    if calc >= obs: ### for enrichment
                        values[n]+=1
        pval, nodenames = OM.get_p(values, x)
        padj = OM.get_p_adj(pval)
        i = 0
        for n in nodenames:
            # if padj[i]  < 0.05: ## include this when creating table3 and table4, remove for tablex3
            dades[sp][loc][n][1] = padj[i]
            dades[sp][loc][n][2] = pval[i]
            i+=1
outfile.close()

locs = parasite+peripheral+exported

plasmo = ['PLAF7','PLAKH','1323249','PLABA','31271']
nodesall = nodes+plasmo
nodesall.reverse()

# outfile1 = open(table3, 'w')
# outfile2 = open(table4, 'w')

# head = 'Node'
# for loc in locs:
#     for sp in plasmo:
#         head+='\t'+loc+'-'+sp
# print(head,file=outfile1)
# print(head,file=outfile2)

# for n in nodesall:
#     string1 = n
#     string2 = n
#     for loc in locs:
#         for sp in plasmo:
#             if n in dades[sp][loc]:
#                 num1,num2=dades[sp][loc][n][0],dades[sp][loc][n][1] #### percentage, p-value (enrichment)
#             else:
#                 num1,num2=0,0.1
#             string1 += '\t'+str(num1)
#             string2 += '\t'+str(num2)
#     print(string1,file=outfile1)
#     print(string2,file=outfile2)
# outfile1.close()
# outfile2.close()

## this is tablex3
outfile = open(tablex3,'w')
print('Phylostratum\tSpecies\tSubcellular fraction\tN genes\t% genes\tp-value\tp-adj',file=outfile)
species = OM.loadTaxa(taxaFile)
for n in nodesall:
    for loc in locs:
        for sp in plasmo:
            name = 'P.'+species[sp].split(' ')[1]
            if n in dades[sp][loc]:
                values = [str(x) for x in dades[sp][loc][n]]
                string = node2phylo[n]+'\t'+name+'\t'+loc+'\t'+values[3]+'\t'+values[0]+'\t'+values[2]+'\t'+values[1]
            else:
                string = node2phylo[n]+'\t'+name+'\t'+loc+'\t'+'0'+'\t'+'0'+'\t'+'-'+'\t'+'-'
            print(string,file=outfile)
outfile.close()

## plot dots
##

#transform_table_dot(table3)
# table = table3.split('.')[0]+'_transform.csv'

# df1 = pd.read_csv(table,sep='\t',header=0,index_col=0)
# df2 = pd.read_csv(table4,sep='\t',header=0,index_col=0)

# ylabels = df1.index.values
# xlabels = df1.columns.values
# M,N = len(xlabels),len(ylabels)
# x, y = np.meshgrid(np.arange(M), np.arange(N))

# fig, ax = plt.subplots(figsize=(21.3,9))

# R = df1.to_numpy()/2

# circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
# col = PatchCollection(circles, array=df2.to_numpy().flatten(), cmap="cividis")
# ax.add_collection(col)

# ax.set(xticks=np.arange(M), yticks=np.arange(N),
#         xticklabels=xlabels, yticklabels=ylabels)
# ax.set_xticks(np.arange(M+1)-0.5, minor=True)
# ax.set_yticks(np.arange(N+1)-0.5, minor=True)
# ax.grid(which='minor')
# plt.xticks(rotation=45, ha='right')

# fig.colorbar(col)
# plt.savefig(outenrichplot,bbox_inches = "tight")
# plt.show()


### Create a supplementary table with the percentage and the p-val

# species = OM.loadTaxa(taxaFile)

# df1 = pd.read_csv(table3,sep='\t',header=0,index_col=0) # % of genes
# df2 = pd.read_csv(table4,sep='\t',header=0,index_col=0) # p-values

# nodes2loc2per = df1.to_dict('index')
# nodes2loc2pval = df2.to_dict('index')

# nodes = list(nodes2loc2per.keys())
# nod1 = [x for x in nodes if 'NODE' in x]
# nod1 = sorted(nod1, key=lambda x:int(x.split('_')[1]))
# nod2 = [x for x in nodes if 'NODE' not in x]
# nodes = nod1+nod2

# for n in nodes:
#     for k in nodes2loc2per[n]:
#         loc,sp = k.split('-')
#         sp = 'P. ' +species[sp].split(' ')[1]
#         string = 







### Plot Anthony table
# inFile = "/home/ijulca/projects/Malaria/analysis/Classeur2.csv"
# #transform_table_dot(inFile)
# table = inFile.split('.')[0]+'_transform.csv'
# outfig = '/home/ijulca/projects/Malaria/analysis/Classeur2.svg'

# df1 = pd.read_csv(table,sep='\t',header=0,index_col=0)

# ylabels = df1.index.values
# xlabels = df1.columns.values
# M,N = len(xlabels),len(ylabels)
# x, y = np.meshgrid(np.arange(M), np.arange(N))

# fig, ax = plt.subplots(figsize=(14,3.5))

# R = df1.to_numpy()/2

# circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
# col = PatchCollection(circles, array=df1.to_numpy().flatten(), cmap="cividis")
# ax.add_collection(col)

# ax.set(xticks=np.arange(M), yticks=np.arange(N),
#         xticklabels=xlabels, yticklabels=ylabels)
# ax.set_xticks(np.arange(M+1)-0.5, minor=True)
# ax.set_yticks(np.arange(N+1)-0.5, minor=True)
# ax.grid(which='minor')
# plt.xticks(rotation=45, ha='right')

# fig.colorbar(col)
# plt.savefig(outfig,bbox_inches = "tight")
# plt.show()


#### Ploat Anthony table2

# inFile = "/home/ijulca/projects/Malaria/analysis/Classeur22.csv"
# transform_table_dot(inFile)
# table = inFile.split('.')[0]+'_transform.csv'
# outfig = '/home/ijulca/projects/Malaria/analysis/Classeur2.svg'

# df1 = pd.read_csv(table,sep='\t',header=0,index_col=0)

# ylabels = df1.index.values
# xlabels = df1.columns.values
# M,N = len(xlabels),len(ylabels)
# x, y = np.meshgrid(np.arange(M), np.arange(N))

# fig, ax = plt.subplots(figsize=(5.5,5))

# R = df1.to_numpy()/2

# circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
# col = PatchCollection(circles, array=df1.to_numpy().flatten(), cmap="cividis")
# ax.add_collection(col)

# ax.set(xticks=np.arange(M), yticks=np.arange(N),
#         xticklabels=xlabels, yticklabels=ylabels)
# ax.set_xticks(np.arange(M+1)-0.5, minor=True)
# ax.set_yticks(np.arange(N+1)-0.5, minor=True)
# ax.grid(which='minor')
# plt.xticks(rotation=45, ha='right')

# fig.colorbar(col)
# plt.savefig(outfig,bbox_inches = "tight")
# plt.show()


###################### 7) Table 5 and 6
# orthogroups = OM.loadOrthofinder(orthoFile)
# sp2gene2node = get_sp2gene2node(nodeFiles) ### genes only from Omar table

# plasmodium = ['PLAF7', 'PLABA', '31271', 'PLAKH', '1323249', 'PLAVI']
# hosts = ['HUMAN', 'MOUSE']

# table = {}
# g = 0
# for group in orthogroups:
#     species = OM.get_species(orthogroups[group])
#     plas = [x for x in species if x in plasmodium]
#     host = [x for x in species if x in hosts]
#     if len(plas)>0 and len(host)>0:
#         g+=1
#         if len(plas) == 1:
#             key = plas[0]
#             pep =[0,0] ## just to continue
#             #pep = [x for x in species[key] if x in list(sp2gene2node[key].keys())] ### check if there are genes in Omar table. Not for PLAVI
#             if len(pep) > 0:
#                 #g+=1
#                 if key not in table:
#                     table[key] = [set([]),set([]),set([])] # Human, mouse, both
#                 if len(host) == 2:
#                     table[key][2].add(group)
#                 elif 'HUMAN' in host:
#                     table[key][0].add(group)
#                 elif 'MOUSE' in host:
#                     table[key].add(group)
#         else:
#             toprint = False
#             for sp in plas:
#                 pep =[0,0] ##
#                 #pep = [x for x in species[sp] if x in list(sp2gene2node[sp].keys())]
#                 if len(pep) > 0: #### check if at least one species has a protein in the list of omar
#                     toprint = True
#             if toprint == True:
#                 #g+=1
#                 key = ''
#                 for sp in plasmodium:
#                     if sp in plas:
#                         key += '-'+sp
#                 if key not in table:
#                     table[key] = [set([]),set([]),set([])]
#                 if len(host) == 2:
#                     table[key][2].add(group)
#                 elif 'HUMAN' in host:
#                     table[key][0].add(group)
#                 elif 'MOUSE' in host:
#                     table[key][1].add(group)
                          
# print('Number of orthogroups analysed:', g)
# outfile = open(table5, 'w')
# print('Plasmodium\tHUMAN\tMOUSE\tBOTH\t%HUMAN\t%MOUSE\t%BOTH', file=outfile)
# for key in table:
#     #values = [str(x) for x in table[key]]
#     values = [str(len(x)) for x in table[key]]
#     pervalues = [str(len(x)*100/g) for x in table[key]]
#     string = key+'\t'+'\t'.join(values)+'\t'+'\t'.join(pervalues)
#     print(string,file=outfile)
# outfile.close()

# fun = {}
# for line in open('/home/ijulca/projects/Malaria/analysis/PLAF7_nodes.csv'):
#     line = line.strip()
#     data = line.split('\t')
#     fun[data[0]] = data[1]

# outfile = open(table6, 'w')
# print('orthogroup\tplasmodium\thost\t'+'\t'.join(plasmodium+hosts)+'\tOther species\tFunction in PLAF7',file=outfile)
# for key in table:
#     for i in range(0,len(table[key])):
#         if i == 0:
#             name = 'HUMAN'
#         elif i == 1:
#             name = 'MOUSE'
#         elif i == 2:
#             name = 'BOTH'
#         else:
#             print('Error number of elements...')
#         for g in table[key][i]:
#             species = OM.get_species(orthogroups[g])
#             string = g +'\t'+key+'\t'+name
#             pep_fun = []
#             for sp in list(plasmodium+hosts):
#                 if sp in species:
#                     sp_pep = [x.split('-')[0] for x in species[sp]]
#                     string += '\t'+';'.join(sp_pep)
#                     if sp == 'PLAF7':
#                         pep_fun = [fun[x] for x in sp_pep if x in fun]
#                 else:
#                     string += '\t-'
#             print(pep_fun)
#             pep = [x.split('-')[0] for x in orthogroups[g] if x.split('-')[-1] not in list(plasmodium+hosts)]
#             string +='\t'+ ';'.join(pep)+'\t'+';'.join(pep_fun)
#             print(string,file=outfile)
# outfile.close()


#### Calculate total percentages
# for f in infiles:
#     sp = OM.get_prefix(f).split('_all')[0]
#     df = pd.read_csv(f,sep='\t',header=0)
#     df2 = df[df['FINAL LOCALIZATION'].isin(parasite)]
#     genes = set(df2['GeneID'].values)
#     print(sp, 'Parasite', len(genes), len(genes)*100/proteins[sp])
#     df2 = df[df['FINAL LOCALIZATION'].isin(exported)]
#     genes = set(df2['GeneID'].values)
#     print(sp, 'Exported', len(genes), len(genes)*100/proteins[sp])

########################################################################################################3
#######################  OLD
###################### Calculate p-value
# sp2gene2node = get_sp2gene2node(nodeFiles) ### genes only from Omar table

# nodes = ['NODE_1', 'NODE_3','NODE_4','NODE_6','NODE_7','NODE_8','NODE_9','NODE_11','NODE_13']
# sp2nodes = {'PLAF7':nodes+['PLAF7'],'PLAKH':nodes+['NODE_15','PLAKH'],'31271':nodes+['NODE_15','NODE_16','31271'], 
#             '1323249':nodes+['NODE_15','NODE_16','NODE_17','1323249'],'PLABA':nodes+['NODE_15','NODE_16','NODE_17','PLABA']}

# gene2node = OM.load_ortho2node(outGene2node)
# sp2gene2node = {} ### all proteome
# for g in gene2node:
#     sp = g.split('-')[-1]
#     if sp not in sp2gene2node:
#         sp2gene2node[sp] = {}
#     sp2gene2node[sp][g] = gene2node[g]

# x=1000

# parasite = ['PARASITE','parasite','PPM']
# paras = ['PARASITE','parasite']
# exported = ["Cleft's","EXPORTED","GHOST","HCC","PV","PVM","PV or PVM","VESICLE"]

# table = get_sp2gene2loc(infiles)
# sp2loc2gene, loctag = {}, set(['Parasite_all','Parasite','Exported_all'])
# for sp in table:
#     sp2loc2gene[sp] = {'Parasite_all':set([]), 'Exported_all':set([]),'Parasite':set([])}
#     for g in table[sp]:
#         loc = table[sp][g]
#         if loc in parasite:
#             sp2loc2gene[sp]['Parasite_all'].add(g)
#             if loc in paras:
#                 sp2loc2gene[sp]['Parasite'].add(g)
#             else:
#                 if loc not in sp2loc2gene[sp]:
#                     sp2loc2gene[sp][loc] = set([])
#                     loctag.add(loc)
#                 sp2loc2gene[sp][loc].add(g)
#         if loc in exported:
#             sp2loc2gene[sp]['Exported_all'].add(g)
#             if loc not in sp2loc2gene[sp]:
#                 sp2loc2gene[sp][loc] = set([])
#                 loctag.add(loc)
#             sp2loc2gene[sp][loc].add(g)
        


# OM.get_enrichment2node(Node2Enrich, sp2loc2gene, sp2gene2node, sp2nodes, x)

# loctag = ['Parasite_all','Parasite','PPM','Exported_all','EXPORTED',"Cleft's", 'GHOST', 'VESICLE','HCC', 'PV', 'PVM', 'PV or PVM']
# get_enrichNodePlot(Node2Enrich, sp2nodes, Node2EnrichPlot, loctag)

# a= 0000