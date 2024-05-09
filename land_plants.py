#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 20:11:54 2024

@author: ijulcach
"""
import sys,os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import ete3
import numpy as np
import glob


def get_species(treeFile):
    t = ete3.Tree(treeFile)
    names = []
    for leaf in t:
        names.append(leaf.name)
    return names

def load_hogs(inFile):
    table = {}
    hogs = {}
    for line in open(inFile):
        line = line.strip()
        if line.startswith('RootHOG'):
            pass
        else:
            data = line.split('\t')
            if data[0] not in table:
                table[data[0]] = {'g':set(),'h':set()}
            table[data[0]]['g'].add(data[1])
            table[data[0]]['h'].add(data[2])
            if data[2] not in hogs:
                hogs[data[2]] = set()
            hogs[data[2]].add(data[0])
    return table, hogs

def get_species_hog(hog):
    species = [x[:5] for x in hog]
    values, counts = np.unique(species, return_counts=True)
    return values, counts

def create_new_format(table, oma, outname1, outname2):
    outfile = open(outname1,'w')
    print('HOG\toma_HOG\tn_genes\tn_species\tgenes', file=outfile)
    for e in table:
        genes = table[e]['g']
        ohog = list(table[e]['h'])
        species, count = get_species_hog(genes)
        string = e+'\t'+ohog[0]+'\t'+str(len(genes))+'\t'+str(len(species))+'\t'+';'.join(genes)
        print(string,file=outfile)
    outfile.close()
    outfile = open(outname2,'w')
    print('omahog\ttag\tn_rhog\trhog\tsize_rhog\tns_rhog',file=outfile)
    for hog in oma:
        rhogs = list(oma[hog])
        if len(rhogs) == 1:
            key = 'ONE'
        else:
            key = 'SPLIT'
        sizes,spec = [],[]
        for e in rhogs:
            genes = table[e]['g']
            sizes.append(len(genes))
            species, counts = get_species_hog(genes)
            spec.append(len(set(species)))
        string = hog+'\t'+key+'\t'+str(len(rhogs))+'\t'+';'.join(rhogs)+'\t'+';'.join([str(x) for x in sizes])
        string +='\t'+';'.join([str(x) for x in spec])
        print(string,file=outfile)
    outfile.close()
        
def get_number_genes_in_hogs(table):
    data = {}
    for e in table:
        genes = table[e]['g']
        species,count = get_species_hog(genes)
        for i,s in enumerate(species):
            if s not in data:
                data[s] = 0
            data[s]+=count[i]
    return data

def general_stat(protpath, table, outname):
    protFile = glob.glob(protpath+'/*.fa')   
    proteins = {}
    for f in protFile:
        seqs = GM.load_sequences(f)
        n = f.split('/')[-1].split('.')[0]
        proteins[n] = len(list(seqs.keys()))
    hog_data = get_number_genes_in_hogs(table)
    print('Number of species', len(proteins))
    outfile = open(outname,'w')
    print('species\tgenes\thog_genes\tper_hog_genes', file=outfile)
    for s in proteins:
        npep,nhog = proteins[s], hog_data[s]
        string = s+'\t'+str(npep)+'\t'+str(nhog)+'\t'+str(nhog*100/npep)
        print(string,file=outfile)
    outfile.close()
                                        

#### plots
def plot_stat_species(inFile, outname):
    df = pd.read_csv(inFile, sep='\t')
    df = df.sort_values(by='genes')
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 10), sharey=True)
    axes = axes.flatten()
    sns.barplot(data=df, y='species',x='genes', ax=axes[0])
    axes[0].set_title('Number of genes per species')
    sns.barplot(data=df, y='species',x='per_hog_genes', ax=axes[1])
    axes[1].vlines(50,0,57, color='red', linestyle='--')
    axes[1].vlines(80,0,57, color='green', linestyle='--')
    axes[1].set_title('Percentage of genes in HOGs')
    plt.savefig(outname, bbox_inches='tight')
    plt.show()

def plot_hog_number(inFile, outfig):
    df = pd.read_csv(inFile, sep='\t')
    df['n_genes'].where(df['n_genes'] <= 100, 100, inplace=True)
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 3))#, sharey=True)
    axes = axes.flatten()
    sns.histplot(data=df, x='n_genes', ax=axes[0])#, bins=100)
    sns.histplot(data=df, x='n_species', ax=axes[1])
    plt.savefig(outfig+'_hist.svg', bbox_inches='tight')
    plt.show()
    ax = sns.scatterplot(data=df,x='n_genes', y='n_species')
    plt.savefig(outfig+'_scatter.svg', bbox_inches='tight')
    plt.show()
    print(df.columns.values)

#####################################
######################################
########## FastOMA ###################
######################################
######################################
path = '/home/ijulcach/projects/Land_Plants/Fastoma_plants/oma_plants_fastoma/'
rhogFile = path +'RootHOGs.tsv'
protpath = '/home/ijulcach/projects/Land_Plants/Fastoma_plants/oma_plants/proteome/'

#### outfiles
hogsFile = path+'RootHOGs.line.tsv'
omahogfile = path +'oma_root_split.tsv'
species_stat = path +'species_stat.tsv'

species_statFig = path+'plots/'+'species_stat.svg'
outhogfig1 = path+'plots/'+'hogs_size'
outhogfig2 = path+'plots/'+'oma_hog_split.svg'


####
# table, oma = load_hogs(rhogFile)
# create_new_format(table, oma, hogsFile, omahogfile)
# general_stat(protpath, table, species_stat)
# plot_stat_species(species_stat, species_statFig)
# plot_hog_number(hogsFile, outhogfig1)

df = pd.read_csv(omahogfile, sep='\t')
df = df[df['tag']=='SPLIT']
df['n_rhog'] = df['n_rhog'].apply(lambda x:(x if x<20 else 20))
print(df.columns.values)
ax = sns.histplot(data=df, x='n_rhog', color='red')#, bin=10)
plt.savefig(outhogfig2, bbox_inches='tight')
plt.show()






# path = '/home/ijulcach/projects/Land_Plants/'
# plotPath = path +'plot_results/'
# numbersFile = '/home/ijulcach/projects/Land_Plants/numbers_genomes_oma.tsv.csv'
# treemajorLin = path+'major_lineages.nw'
# treeangios = path+'angiosperms_tree.nw'
# faoFile = path+'Fao_data/FAOSTAT_data_en_3-22-2024_2000_2022.csv'
# ##### OutputFiles
# outfig1 = plotPath+'major_lineages_species.svg'
# outfig2 = plotPath+'major_lineages_genomes.svg'

# ### Format FAO data
# # i = 0
# # for line in open(faoFile):
# #     line = line.strip()
# #     data = line.split('","')
# #     if 'Domain Code' in data[0]:
# #         print(data)
# #     else:
# #         food, year = data[7], data[9]
# #         print(food,year)
# #         print(data[10:])
    


# ## get the order of species to be plotted:
# majorLin = get_species(treemajorLin)
# angiosperms = get_species(treeangios)
# y = 'major_lin'# 'clades' #'major_lin'
# ## load data
# df = pd.read_csv(numbersFile, sep='\t')
# # print(df.columns.values)
# # print(df['genome_avail'].sum()) ## number of species
# # ## select major lineages
# df2 = df.drop(['clades'], axis=1) 
# df2 = df2.groupby(by=['major_lin']).sum()
# df2 = df2.reindex(index = majorLin) # reindex based in the order of species
# # ## select clades of angiosperms
# # df2 = df.drop(['major_lin'], axis=1)
# # df2 = df2.set_index('clades')
# # noangios = [x for x in list(df['clades'].values) if x not in angiosperms]
# # df2 = df2.drop(noangios)
# # print('number of angiosperms', df2['number_species'].sum())
# # ## calculate percentages
# df2['num_per'] = df2['number_species']*100/df2['number_species'].sum()
# df2['geno_per'] = 100*((df2['num_per']/100)*(df2['genome_avail']/df2['number_species']))
# # print('% of species sequenced',df2['genome_avail'].sum()*100/df2['number_species'].sum())
# # ## plot 1
# plt.figure(figsize=(4,6)) ##4,6 ###for major, 4,8
# # # ax = sns.barplot(df2, x="num_per", y=y, color='#98D8AA')
# # # # ax.set_xlim(0,100)
# # # ax = sns.barplot(df2, x='geno_per', y=y, color='#4D96FF')
# # # ax.set_xlabel('% Species')
# # # ax.set_ylabel("")
# # # plt.savefig(outfig1, bbox_inches='tight')
# # # plt.show()

# # ## plot 2
# # df2.loc['Ceratophyllales','genome_avail'] =0
# # df2.loc['Amborellales','genome_avail'] =0
# df2['geno_tper'] = 100*df2['genome_avail']/df2['number_species']
# df2['geno_tper'] = df2['geno_tper'].where(df2['geno_tper'] <= 3, 3)

# df2['oma_per'] = 100*((df2['geno_tper']/100)*(df2['oma_species']/df2['genome_avail']))
# df2['oma_nper'] = 100*((df2['geno_tper']/100)*(df2['future_oma']/df2['genome_avail']))
# ax = sns.barplot(df2, x='geno_tper', y=y, color='#4D96FF')
# ax = sns.barplot(df2, x='oma_nper', y=y, color='#6BCB77')
# ax = sns.barplot(df2, x='oma_per', y=y, color='#FFD93D')
# ax.set_xticklabels(ax.get_xticks().astype(int))
# ax.set_xlabel('% Sequenced genomes')
# ax.set_ylabel("")
# plt.savefig(outfig2, bbox_inches='tight')
# plt.show()





