#!/usr/bin/env python
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
import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr
from scipy.stats import mannwhitneyu
from scipy.stats import ranksums
from statsmodels.sandbox.stats.multicomp import multipletests
import matplotlib.colors
import math


def check_spider(log):
    toprint = False
    if os.path.isfile(log) == True:
        with open(log, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
            if 'Date and Time:' in last_line:
                toprint = True
    return toprint

def load_taxa(inFile):
    table = {}
    for line in open(inFile):
        line = line.strip()
        data = line.split('\t')
        table[data[0]] = data[1].split(';')[-1].strip()
    return table

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

def get_tree_panther(inFile):
    i = 0
    for line in open(inFile):
        line = line.strip()
        i += 1
        if i == 1:
            tree = line
    return tree

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

def get_conversion(inFile):
    table = {}
    for line in open(inFile):
        line = line.strip()
        data = line.split('\t')
        name = data[1].split('=')[-1]
        table[name] = data[0]
    return table

def get_sample_id(inFile):
    sam = {}
    s = '-'
    for line in open(inFile):
        line = line.strip()
        line = line.replace('"','')
        if '##SampleName ' in line:
            if s == '-':
                s = line.split(' ')[1]
                g,n = '-','-'
            else:
                if n != '-':
                    n = n.replace('-','').strip()
                sam[s] = g+'\t'+n
                s = line.split(' ')[1]
                g,n = '-','-'
        else:
            if '/growth condition' in line:
                g = line.split('=')[1]
            elif '/tissue=' in line:
                n = line.split('=')[1]
            elif '/sampling site=' in line:
                n += ' '+line.split('=')[1]
            elif '/source name=' in line:
                n += ' '+line.split('=')[1]
            elif '/plant structure' in line:
                n += ' '+line.split('=')[1]
            elif '/organism part' in line:
                n += ' '+line.split('=')[1]
            elif '/plant body site' in line:
                n += ' '+line.split('=')[1]
            elif '/cell type' in line:
                n += ' '+line.split('=')[1]
            elif '/sample type' in line:
                n += ' '+line.split('=')[1]
    sam[s] = g+'\t'+n
    return sam

def get_wilcoxon_rank_sum(list1, list2):
    wilran = ranksums(list1, list2)
    # stat, p = mannwhitneyu(list1,list2)
    #mannwhitneyu(list1, list2).pvalue
    return wilran[1] # pvalue of the test

def get_p_adj(pval):
    pval_flt = [x for x in pval if str(x) != 'nan'] ### remove nan
    p_new = multipletests(pval_flt, alpha=0.05, method='fdr_bh') ## bonferroni
    p_adj = [x for x in p_new[1]]
    p_adjusted = []
    for item in pval:
        if str(item) == 'nan':
            p_adjusted.append(item)
        else:
            p_adjusted.append(p_adj.pop(0))    
    return p_adjusted

def plot_box_pval(categories, table, outfigure):
    pval = []
    for n in categories[1:]:
        d1 = table[categories[0]]
        d2 = table[n]
        p = get_wilcoxon_rank_sum(d1,d2)
        pval.append(p)
    p_adjusted = get_p_adj(pval)           
            
    dades = [[]]*len(categories)
    for i,n in enumerate(categories):
        dades[i] = table[n]

    # plot
    ax = sns.boxplot(data=dades, palette=color, fliersize=0)#, saturation=0.7, linecolor=color2)
    # ax = sns.swarmplot(data=dades, color=".25")

    i,j = 0,1.1
    for p in p_adjusted:
        p = float(p)
        if p < 0.01:
            sn = np.format_float_scientific(p, precision=1, exp_digits=1)+"**"
        elif p < 0.05:
            sn = np.format_float_scientific(p, precision=1, exp_digits=1)+"*"
        else:
            sn = str(round(p,3))
        x1,x2 = 0,i+1
        y, h, col, z = j+0.04, 0.05+i/4, 'k', 0.03
        plt.plot([x1, x1, x2, x2], [y-z, y+h-z, y+h-z, y-z], lw=0.5, c=col)
        plt.text((x1+x2)*.5, y+h-z, sn, ha='center', va='bottom', color=col, fontsize=8)
        i+=1
        # ax.set_ylim([0.6,1.05])
    _ = plt.xticks([0,1,2,3,4,5], categories)
    _ = plt.xticks(rotation=45, ha='right')
    ticks_to_keep = ax.get_yticks()[1:6] ### remove some yticks
    ax.set_yticks(ticks_to_keep)
    plt.savefig(outfigure, bbox_inches='tight') #outfig1, outfig2, outfig3
    plt.show()

def plot_box_pval_sp(categories, table, outfigure):
    fig, axes = plt.subplots(nrows=4, ncols=5, figsize=(15, 10), sharex=True, sharey=True)
    axes = axes.flatten()
    species = [x for x in list(table.keys()) if x!='None']
    for x,sp in enumerate(species):
        print(x,sp)
        data = table[sp]
        pval = []
        for n in categories[1:]:
            d1 = data[categories[0]]
            d2 = data[n]
            p = get_wilcoxon_rank_sum(d1,d2)
            pval.append(p)
        p_adjusted = get_p_adj(pval)
        dades = [[]]*len(categories)
        for i,n in enumerate(categories):
            dades[i] = data[n]
        
        # plot
        sns.boxplot(data=dades, palette=color, fliersize=0, ax=axes[x])
        i,j = 0,1.1
        for p in p_adjusted:
            p = float(p)
            if p < 0.01:
                sn = np.format_float_scientific(p, precision=1, exp_digits=1)+"**"
            elif p < 0.05:
                sn = np.format_float_scientific(p, precision=1, exp_digits=1)+"*"
            else:
                sn = str(round(p,3))
            x1,x2 = 0,i+1
            y, h, col, z = j+0.04, 0.05+i/4, 'k', 0.03
            axes[x].plot([x1, x1, x2, x2], [y-z, y+h-z, y+h-z, y-z], lw=0.5, c=col)
            axes[x].text((x1+x2), y+h-z, sn, ha='center', va='bottom', color=col, fontsize=8)
            i+=1
        axes[x].set_title(sp)
        axes[x].tick_params(axis='both', which='both', length=0)
    
    ## plot adjustments
    for ax in axes:
        ax.set_xticks([0,1,2,3,4,5])
        ax.set_xticklabels(categories, rotation=45, ha='right')
        ax.set_yticks([-1,-0.5, 0, 0.5, 1, 1.5, 2, 2.5])
        ticks_to_keep = ax.get_yticks()[0:5] ### remove some yticks
        ax.set_yticks(ticks_to_keep)
    plt.savefig(outfigure, bbox_inches='tight')
    plt.show()

    
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
### infiles and outfiles

pathPlots = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/results/plots/'
pathTables = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/results/tables/'
treepath = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/panther-18.0/trees/'
speciesFile = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/panther-18.0/species_tree.nhx'
expPlantdata = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/Plant_expression/Samples/'
taxaFile = '/home/ijulcach/projects/Land_Plants/nemo2taxa.txt'

inFile = pathTables+'pairwise_tests.tsv'
linFile = pathTables+'major_lineages_sp.txt'
expPath = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/Plant_expression/Exp_matrices/'
# outfile1 = pathTables+'plant_exp_pcc_scaler.tsv'
# outfig1 = pathPlots+'plant_exp_pcc_scaler.svg'
# outfile1 = pathTables+'plant_exp_pcc_noscaler.tsv'
# outfig1 = pathPlots+'plant_exp_pcc_noscaler.svg'
# outfile1 = pathTables+'plant_exp_pcc_zscore.tsv' 
# outfig1 = pathPlots+'plant_exp_pcc_zscore.svg' 
# outfile1 = pathTables+'plant_exp_pcc_zscoremedian.tsv'
# outfig1 = pathPlots+'plant_exp_pcc_zscoremedian.svg'
# outfile1 = pathTables+'plant_exp_pcc_Snames.zscoreLog.tsv'
# outfig1 = pathPlots+'plant_exp_pcc_Snames.zscoreLog.svg'
outfile1 = pathTables+'plant_exp_pcc_zscoreLog.tsv'  ### I use this
outfig1 = pathPlots+'plant_exp_pcc_zscoreLog.svg'    ### I use this
outfig2 = pathPlots+'plant_exp_pcc_zscoreLog_nodesS.svg'
outfig3 = pathPlots+'plant_exp_pcc_zscoreLog_nodesI.svg'
outfig4 = pathPlots +'plant_exp_pcc_zscoreLog_sp.svg'
outfig5 = pathPlots +'plant_exp_pcc_zscoreLog_spS.svg'
outfig6 = pathPlots +'plant_exp_pcc_zscoreLog_spI.svg'

files = glob.glob('/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/Plant_expression/Exp_matrices/*/*.Gnames.sample.txt')
outname1 = pathTables+'plant_sampleSPM.tsv'
outname2 = pathTables + 'plant_sampleSPM_matrix.tsv'
outfig7 = pathPlots+'plant_samplesSPM.svg'
outfig8 = pathPlots+'plant_samplesSPMS.svg'
outfig9 = pathPlots+'plant_samplesSPMI.svg'
outfig10 = pathPlots+'plant_samplesSPM_sp.svg'
outfig11 = pathPlots+'plant_samplesSPM_spS.svg'
outfig12 = pathPlots+'plant_samplesSPM_spI.svg'
outfig13 = pathPlots+'supplementary_spm_plants.svg'
outsampTable = pathTables+'sp2number_samples.tsv'
outfigs = pathPlots+'sp2number_samples.svg'

### load names
sp2name = load_taxa(taxaFile)
#### Get number of duplications
treeFiles = glob.glob(treepath+'*')

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

# lineages = ['Deuterostomia', 'Protostomia', 'Fungi', 'Amoebozoa', 'Excavates',
#             'Alveolata-Stramenopiles','Viridiplantae', 'Archaea', 'Eubacteria']
# indLin = ['NEMVE','TRIAD','MONBE']
# lineages = {x:set() for x in lineages}
# lineages.update({x:set([x]) for x in indLin})

# t = ete4.Tree(open(speciesFile))

# for node in t.traverse('postorder'):
#     props = node.props
#     if 'S' in props:
#         lin = props['S']
#         if lin in lineages:
#             lineages[lin].add(lin)
#             for leaf in node:
#                 lineages[lin].add(leaf.name)
#             for n in node.traverse():
#                 info = n.props
#                 if 'S' in info:
#                     name = info['S']
#                     lineages[lin].add(name)
# i = 0
# outfile = open(pathTables+'major_lineages_sp.txt','w')
# for e in lineages:
#     string = e+'\t'+str(len(lineages[e]))+'\t'+'; '.join(lineages[e])
#     print(string,file=outfile)
#     i += len(lineages[e])
# print(i)
# outfile.close()

# internal_nodes = ['Bilateria', 'Eumetazoa', 'Metazoa-Choanoflagellida',
#                   'Opisthokonts', 'Unikonts', 'Eukaryota', 'Bikonts',
#                   'SAR/HA_supergroup', 'Archaea-Eukaryota', 'LUCA']
## Node 1: 'Bilateria'
## Node 2: 'Eumetazoa'
## Node 3: 'Metazoa-Choanoflagellida'
## Node 4: 'Opisthokonts'
## Node 5: 'Unikonts'
## Node 6: 'Eukaryota'
## Node 7: 'Bikonts'
## Node 8: 'SAR/HA_supergroup'
## Node 9: 'Archaea-Eukaryota'
## Node 10: 'LUCA'


##### Get the categories tables
categories = ['normal-normal', 'short-short', 'long-long', 'normal-short', 
             'normal-long', 'short-long']
# table = {}
# dups = {}
# for line in open(pathTables+'duplication_events_with_stat.tsv'):
#     line = line.strip()
#     data = line.split()
#     if data[0] != 'fam_id':
#         taxa = data[4] ### sp_tree_head
#         fam = data[0]
#         p1,p2 = float(data[14]),float(data[15]) ## all normal normal have a p >sig
#         if fam not in dups:
#             dups[fam] = 0
#         dups[fam] += 1  
#         key1 = data[16]+'-'+data[17]
#         key2 = data[17]+'-'+data[16]
#         if taxa not in table:
#             table[taxa] = {x:0 for x in categories}
#         if key1 in table[taxa]:
#             key = key1
#         elif key2 in table[taxa]:
#             key = key2
#         table[taxa][key] += 1
# print('Total number of duplications:', np.sum([dups[x] for x in dups]))
# print('Average duplications per family:', np.mean([dups[x] for x in dups]))
# outfile = open(pathTables+'nodes2categories.tsv','w')
# print('Taxa\t'+'\t'.join(categories), file=outfile)
# for tax in table:
#     string = tax
#     for e in categories:
#         string += '\t'+str(table[tax][e])
#     print(string, file=outfile)
# outfile.close()

# table = {x:[0.0]*len(categories) for x in lineages}
# table.update({x:[0.0]*len(categories) for x in internal_nodes})
# for line in open(pathTables+'nodes2categories.tsv'):
#     line = line.strip()
#     data = line.split('\t')
#     if data[0] != 'Taxa':
#         nums = [float(x) for x in data[1:]]
#         tax = data[0]
#         if tax in table:
#             table[tax] =[x + y for x, y in zip(table[tax], nums)]
#         else:
#             for e in lineages:
#                 if tax in lineages[e]:
#                     tax = e
#             table[tax] =[x + y for x, y in zip(table[tax], nums)]
# outfile = open(pathTables+'nodes2categories.MajorLin.tsv','w')
# print('Taxa\t'+'\t'.join(categories), file=outfile)
# for tax in table:
#     string = tax + '\t' + '\t'.join([str(x) for x in table[tax]])
#     print(string, file=outfile)
# outfile.close()

##### Plot of the categories
# outfigure = pathPlots + 'bar_categories.svg'
# outfigure2 = pathPlots +'bar_dups.svg'
# df = pd.read_csv(pathTables+'nodes2categories.MajorLin.tsv', sep='\t', header=0)
# df['total'] = df.sum(numeric_only=True, axis=1)
# df['totalper']=df['total']*100/df['total'].sum()
# ax = df.plot(x='Taxa', y='total', kind='barh')
# # ax.set_xlim([0, 100])
# plt.savefig(outfigure2, bbox_inches='tight')
# plt.show()
# df.to_csv(path_or_buf=pathTables+'nodes2categories.MajorLin_per.tsv', sep='\t', header=True, index=False)


# df2 = pd.DataFrame()
# for e in categories:
#     df2[e] = df[e]*100/df['total']
# df2['Taxa'] = df['Taxa']
# color = ['#01befe','#ffdd00','#ff7d00','#ff006d','#adff02','#8f00ff']
# ax = df2.plot(x='Taxa', y=categories, kind='barh', stacked=True, color=color, alpha=.7)
# plt.savefig(outfigure, bbox_inches='tight')
# plt.show()
# df2.to_csv(path_or_buf=pathTables+'nodes2categories.MajorLin_per_cat.tsv', sep='\t', header=True, index=False)

##### Plot of the categories per lineage per node
## Create the file
# outname1 = pathTables+'nodes2categories_numsranksMerge.tsv'

# species = set()
# lineages = {}
# for line in open(pathTables+'major_lineages_sp.txt'):
#     line = line.strip()
#     data = line.split('\t')
#     lineages[data[0]] = data[2].split('; ')
#     for e in data[2].split('; '):
#         if len(e) <=5 and e.isupper():
#             species.add(e)
# print('Number of species', len(species))

# table = {}
# for line in open(pathTables+'nodes2categories.tsv'):
#     line = line.strip()
#     data = line.split('\t')
#     if data[0] != 'Taxa':
#         table[data[0]] = [float(x) for x in data[1:]]
# ranks = {}
# for l in lineages:
#     ranks[l] = {}
#     ranks[l]['internal-'+l[0]] = table[l]
#     for tax in lineages[l]:
#         if tax in species:
#             key = 's-s-'+l[0]
#         else:
#             key = 'internal-'+l[0] ### tax
#         if key not in ranks[l]:
#             ranks[l][key]=[0.0]*len(categories)
#         if tax in table:
#             nums = table[tax]
#             ranks[l][key] = [x + y for x, y in zip(ranks[l][key], nums)]
# outfile = open(outname1,'w')
# head = 'Major_Rank\tRank\t'+'\t'.join(categories)
# print(head, file=outfile)
# for l in ranks:
#     for e in ranks[l]:
#         string = l+'\t'+e+'\t'+'\t'.join([str(x) for x in ranks[l][e]])
#         print(string,file=outfile)
# outfile.close()

# ## Create the plot
# outfigure = pathPlots + 'rank_categoriesMerge.svg'
# outfigure2 = pathPlots +'rank_dupsMerge.svg'
# dfname1 = path_or_buf=pathTables+'rank_dupMerge.txt'
# dfname2 = path_or_buf=pathTables+'rank_dupMerge_per.txt'

# df = pd.read_csv(outname1, sep='\t', header=0)
# color = ['#01befe','#ffdd00','#ff7d00','#ff006d','#adff02','#8f00ff']

# ### sort ranks and create a new dataframe
# lin_sort = ['Protostomia','Deuterostomia', 'TRIAD', 'NEMVE','MONBE','Fungi', 
#             'Amoebozoa', 'Excavates','Viridiplantae',
#             'Alveolata-Stramenopiles', 'Archaea', 'Eubacteria']
# ind = ['TRIAD', 'NEMVE','MONBE']
# lin_sort.reverse()
# df_all = pd.DataFrame()
# for mrank in lin_sort:
#     print(mrank)
#     df2 = df[df.Major_Rank.isin([mrank])]
#     df2 = df2.drop(['Major_Rank'], axis=1)
#     # #### sort when using tax
#     # ranks = [x for x in lineages[mrank] if x not in species]
#     # ranks.append('s-s-'+mrank[0])
#     #### end sort
#     ranks = ['s-s-'+mrank[0],'internal-'+mrank[0]]
#     ranks.reverse()
#     rank_category = pd.CategoricalDtype(categories=ranks, ordered=True) ## to sort in order of list
#     if len(ranks) >1:
#         df2['Rank'] = df2['Rank'].astype(rank_category)
#         df2 = df2.sort_values(by='Rank')
#     if mrank in ind:
#         df2 = df2.drop(df2[df2['Rank'] == 's-s-'+mrank[0]].index, axis=0)
    
#     df2['total'] = df2.sum(numeric_only=True, axis=1)
#     df2['totalper']=df2['total']*100/df2['total'].sum()
#     df3 = pd.DataFrame().reindex_like(df2)
#     for e in categories:
#         df3[e] = df2[e]*100/df2['total']
#     df3['Rank'] = df2['Rank']
#     df3['total'] = df2['total']
#     df3['totalper'] = df2['totalper']
#     if df_all.empty:
#         df_all = df3
#     else:
#         df_all = pd.concat([df_all, df3], axis=0)
# df_all.to_csv(dfname1, sep='\t', header=True, index=False)
# df_all.to_csv(dfname2, sep='\t', header=True, index=False)
# ### end sorting
# df_all.set_index('Rank', inplace=True)
# plt.clf()
# fig, ax = plt.subplots()
# bar_width = 0.35

# df_all.plot(y=categories, kind='barh', stacked=True, color=color, alpha=.7,
#                 figsize=(8, 10), width=bar_width, ax=ax)

# bar_positions = range(len(df_all))
# ax.barh([p + bar_width for p in bar_positions],df_all['totalper'], bar_width, 
#         color='grey', label='duplications')
# plt.savefig(outfigure, bbox_inches='tight')

# plt.show()


#### Getting the names of the gene expression matrices and panther
expPlantdata = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/Plant_expression/Exp_matrices/'
expPlantdata = '/home/ijulcach/projects/ldo_project/data/Animal_expression/'

animals = ['XENTR', 'DANRE', 'CANLF']

key = animals[2]
genes = set()
for line in open(expPlantdata+key+'/'+key+'.newMatrix.txt'): #'.cds.fa'):
    line = line.strip()
    # if '>' in line:
    #     genes.add(line.split('>')[1])
    data = line.split('\t')
    if data[0] != 'genes':
        genes.add(data[0])#.split('_')[0])

# dades = {}
# for line in open(expPlantdata+key+'/ensemble2rat.txt'): #'/AGI2uniprot.txt'):
#     line = line.strip()
#     data = line.split('\t')
#     if data[1] not in dades:
#         dades[data[1]] = []
#     dades[data[1]].append(data[0])
#     # dades[data[1].split('-')[0]] =data[0].split('.')[0]

table = {}
for line in open(expPlantdata+key+'/'+key+'.genesNames.txt'):
    line = line.strip()
    data = line.split('\t')
    name = data[1] #'|'.join(data)
    n = name.split('|')[1]#.split('=')[1]#.split('.')[0]
    if 'Ensembl' in name:
        n = n.split('=')[1].split('.')[0]
    if n not in table:
        table[n] = name
    # if 'EnsemblGenome' in data[1]:
    #     n = data[1].split('=')[1]
    # n = data[2].split('=')[1]
    # n = data[1].split('=')[1].replace('CISIN_','orange1.').replace('mg','m.g')
        # table[n] = name

# outfile = open(expPlantdata+key+'/'+key+'.conversion.panther.txt','w')
i = 0
for e in table:
    # if e in dades:
    #     for g in dades[e]:
    if e not in genes:
        print(e)
        i +=1
                # print(e, table[e])
#                 string = g+'\t'+table[e]
#                 print(string,file=outfile)
# outfile.close()
print(key, len(table), i, i*100/len(table)) 

################################
#### Analysis of expression ####     
################################

# expmat = glob.glob(expPath+'*/')
# expsp = [x.split('/')[-2] for x in expmat]
# print('Species with expression data:',len(expsp))

### get the lineages of plants
# table = {}
# for line in open(linFile):
#     line = line.strip()
#     data = line.split('\t')
#     table[data[0]] = data[2].split('; ')
# viridi = table['Viridiplantae']
# species = [x for x in viridi if len(x) <=5 and x.isupper()]
# inter = [x for x in viridi if x not in species]

# #### get the duplications
# dades = {}
# sp2genes = {}
# for line in open(inFile):
#     line = line.strip()
#     data = line.split('\t')
#     if data[0] != 'fam_id':
#         taxa = data[4] ## node where the duplication happened
#         k1,k2 = data[16:18]
#         sp,g1,g2 = data[19:] ## species and gene names of the duplication
#         if sp in expsp: ## analising only sp with expression
#             if taxa not in dades:
#                 dades[taxa] = {x:[] for x in categories}
#             key1 = k1+'-'+k2
#             key2 = k2+'-'+k1
#             if key1 in dades[taxa]:
#                 key = key1
#             elif key2 in dades[taxa]:
#                 key = key2
#                 g1,g2 = g2,g1
#             dades[taxa][key].append(sp+'--'+g1+'--'+g2)
#             if sp not in sp2genes:
#                 sp2genes[sp] = {}
#             sp2genes[sp][g1+'--'+g2] = 'None'

#### get expression values, pearson correlation
# for f in expmat:
#     sp = f.split('/')[-2]
#     print(sp)
#     gnames = get_conversion(f+sp+'.conversion.panther.txt')
#     df = pd.read_csv(f+sp+'.Gnames.tpm.av', sep='\t', header=0, index_col=0) # .Gnames.tpm.av
#     expgenes = df.index.values
#     ex = 0
#     ### Scaler normalization
#     # scaler = StandardScaler()
#     # data_scaler = scaler.fit_transform(df)
#     # df = pd.DataFrame(data_scaler, columns=df.columns,index=df.index)
#     for g in list(sp2genes[sp]):
#         g1,g2 = g.split('--')
#         if g1 in gnames and g2 in gnames:
#             p1,p2 = gnames[g1], gnames[g2]
#             if p1 in expgenes and p2 in expgenes:
#                 row1, row2 = df.loc[p1].values, df.loc[p2].values ## use values for Zlog
#                 if np.sum(row1) == 0 or np.sum(row2) == 0: ### removing genes with only ceros
#                     pass
#                 else:
#                     # x = (row1-row1.mean())/row1.std() ## z-score ## can use median too
#                     # y = (row2-row2.mean())/row2.std()
#                     # x = list(row1/row1.max())
#                     # y = list(row2/row2.max())
#                     x = [(math.log2(x+1)-np.median([math.log2(y+1) for y in row1]))/np.std([math.log2(y+1) for y in row1]) for x in row1]
#                     y = [(math.log2(x+1)-np.median([math.log2(y+1) for y in row2]))/np.std([math.log2(y+1) for y in row2]) for x in row2]
#                     corr,_ = pearsonr(x, y)
#                     sp2genes[sp][g] = corr
#             else:
#                 ex +=1
#             ### continue with Scaler df
#             # if np.sum(row1) == 0: ## avoid rows of 0
#             #     x = list(row1.values)
#             # else:
#             #     x = list(row1/row1.max())
#             # if np.sum(row2) == 0:
#             #     y = list(row2.values)
#             # else:
#             #     y = list(row2/row2.max())
#             # corr,_ = pearsonr(x, y)
#             # sp2genes[sp][g] = corr
# print('genes not in expression',ex)
# outfile = open(outfile1, 'w')
# print('Taxa\tCAT\tSpecies\tgene1\tgene2\tpcc',file=outfile)
# for tax in dades:
#     for cat in categories:
#         data = dades[tax][cat]
#         if len(data) == 0:
#             string = tax+'\t'+cat+'\t'+'None'+'\tNone\tNone\tNone'
#             print(string,file=outfile)
#         else:
#             for e in data:
#                 sp,g1,g2 = e.split('--')
#                 string = tax+'\t'+cat+'\t'+sp+'\t'+g1+'\t'+g2
#                 g = g1+'--'+g2
#                 if g in sp2genes[sp]:
#                     string += '\t'+str(sp2genes[sp][g])
#                 else:
#                     string += '\tNone'
#                 print(string, file=outfile)

# outfile.close()


##### Analysis of pcc

color = ['#01befe','#ffdd00','#ff7d00','#ff006d','#adff02','#8f00ff']
color2 = [matplotlib.colors.to_rgb(x) for x in color]

## Plot of the taxa together and for node
# table = {}
# for line in open(outfile1):
#     line = line.strip()
#     data = line.split('\t')
#     taxa,key,v = data[0],data[1],data[-1]
#     if taxa != 'Taxa':
#         if key not in table:
#             table[key] = []
#         if v != 'None' and v!='nan':
#             # table[key].append(float(v)) ## outfig1
#             if len(taxa) <=5 and taxa.isupper(): ## species, outfig2
#                 # table[key].append(float(v))
#                 pass
#             else:
#                 table[key].append(float(v)) ## internal, outfig3

# plot_box_pval(categories, table, outfig3) ##outfig1, outfig2,outfig3

##### plot per species
# table = {}
# for line in open(outfile1):
#     line = line.strip()
#     data = line.split('\t')
#     taxa,key,sp,v = data[0],data[1],data[2],data[-1]
#     if taxa != 'Taxa':
#         if sp not in table:
#             table[sp] = {x:[] for x in categories}
#         if v != 'None' and v!='nan':
#             # table[sp][key].append(float(v)) ## outfig4
#             if len(taxa) <=5 and taxa.isupper(): ## species, outfig5
#                 # table[sp][key].append(float(v))
#                 pass
#             else:
#                 table[sp][key].append(float(v)) ## internal, outfig6

# plot_box_pval_sp(categories, table, outfig6) ##outfig4, outfig5,outfig6


####### SPM analysis

# genes2panther = {}
# table = {}
# for f in files:
#     sp = f.split('/')[-2]
#     gnames = get_conversion(f.split('.')[0]+'.conversion.panther.txt')
#     genes2panther.update(gnames)
#     for line in open(f):
#         line = line.strip()
#         data = line.split('\t')
#         table[data[0]] = data[2].split('; ')
#     print(sp)

# outfile1 = open(outname1, 'w')

# for line in open(inFile):
#     line = line.strip()
#     data = line.split('\t')
#     if data[0] != 'fam_id':
#         taxa = data[4] ## node where the duplication happened
#         k1,k2 = data[16:18]
#         sp,g1,g2 = data[19:] ## species and gene names of the duplication
#         if sp in expsp: ## analising only sp with expression
#             key1 = k1+'-'+k2
#             key2 = k2+'-'+k1
#             if key1 in categories:
#                 key = key1
#             elif key2 in categories:
#                 key = key2
#                 g1,g2 = g2,g1
#             if g1 in genes2panther and g2 in genes2panther:
#                 p1,p2 = genes2panther[g1], genes2panther[g2]
#                 if p1 in table and p2 in table:
#                     sam1, sam2 = table[p1], table[p2]
#                     if set(sam1) == set(sam2):
#                         SPMdif = 0
#                     else:
#                         if 'Ubiquitous' in sam1 and 'Ubiquitous' not in sam2:
#                             SPMdif = 3
#                         else:
#                             if len(list(set(sam1) & set(sam2))) >0:
#                                 SPMdif = 2
#                             else:
#                                 SPMdif = 1
#                     string = taxa+'\t'+key+'\t'+g1+'\t'+g2+'\t'+p1+'\t'+p2+'\t'+sp
#                     string += '\t'+'; '.join(sam1)+'\t'+'; '.join(sam2)+'\t'+str(SPMdif)
#                     print(string, file=outfile1)
# outfile1.close()

# table = {}
# for line in open(outname1):
#     line = line.strip()
#     data = line.split('\t')
#     taxa = data[0]
#     if data[1] not in table:
#         table[data[1]] = {}
#     k = data[-1]
#     if k not in table[data[1]]:
#         table[data[1]][k] = 0
#     # table[data[1]][k]+=1 ## figure 7
#     if len(taxa) <=5 and taxa.isupper():
#         # table[data[1]][k]+=1 ## figure 8
#         pass
#     else:
#         table[data[1]][k]+=1 ## figure 9

# outfile = open(outname2, 'w')
# nums = ['0','1','2','3']
# head = 'categories\t'+'\t'.join(nums)
# print(head, file=outfile)
# for e in categories:
#     string = e
#     tot = np.sum([table[e][n] for n in nums if n in table[e]])
#     for n in nums:
#         if n in table[e]:
#             v = table[e][n]*100/tot ## get the percentage
#         else:
#             v = 0
#         string += '\t'+ str(v)
#     print(string,file=outfile)
# outfile.close()

# df = pd.read_csv(outname2, sep='\t', header=0, index_col=0)
# print(df.sum(axis=1, numeric_only=True))
# ax = sns.heatmap(df, annot=True, cmap='YlOrBr')#, vmax=10)
# plt.savefig(outfig9, bbox_inches='tight')
# plt.show()


#### per species
# table = {}
# for line in open(outname1):
#     line = line.strip()
#     data = line.split('\t')
#     taxa, sp, k = data[0], data[-4], data[-1]
#     if len(data)>0: ### figure 10, just to have the if
#     # if len(taxa) <=5 and taxa.isupper(): ## figure 11
#     #     pass
#     # else:  ## figure 12
#         if sp not in table:
#             table[sp] = {}
#         if data[1] not in table[sp]:
#             table[sp][data[1]] = {}
#         if k not in table[sp][data[1]]:
#             table[sp][data[1]][k] = 0
#         table[sp][data[1]][k] +=1

# fig, axes = plt.subplots(nrows=4, ncols=5, figsize=(15, 10), sharex=True, sharey=True)
# axes = axes.flatten()
# for i,sp in enumerate((list(table.keys()))):
#     print(sp)
#     outfile = open(outname2, 'w')
#     nums = ['0','1','2','3']
#     head = 'categories\t'+'\t'.join(nums)
#     print(head, file=outfile)
#     dades = table[sp]
#     for e in categories:
#         if e in dades:
#             tot = np.sum([dades[e][x] for x in nums if x in dades[e]])
#             string = e
#             for n in nums:
#                 if n in dades[e]:
#                     v = dades[e][n]*100/tot
#                 else:
#                     v = 0
#                 string += '\t'+str(v)
#         else:
#             string = e+'\t0\t0\t0\t0'
#         print(string,file=outfile)
#     outfile.close()    
#     df = pd.read_csv(outname2, sep='\t', header=0, index_col=0)
#     sns.heatmap(df, annot=True, cmap='YlOrBr', ax=axes[i])
#     axes[i].set_title(sp)
# plt.savefig(outfig10, bbox_inches='tight') ### Figure 10,11,12
# plt.show()


    
# ##### Changing annotation:
# path = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/Plant_expression/'
# infile = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/Plant_expression/SRA_annot_all.tsv'
# outname = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/Plant_expression/SRA_annot_all.tsv2'

# species = ['AMBTC','BRANA','HELAN','MANES','MEDTR','PHYPA','SETIT','SOLTU','WHEAT','ARATH',
#            'CAPAN','CUCSA','MAIZE','MARPO','ORYSJ','POPTR','SELML','SOLLC','SOYBN','VITVI']

# # organs = ['leaf','root', 'stem','meristem','flower'] 

# key1 = 'phloem'
# key2 = 'fruit'

# k = 'stem'

# i = 0
# outfile = open(outname,'w')
# for line in open(infile):
#     line = line.strip()
#     data = line.split('\t')
#     name = data[4].lower()
#     g = data[5]
#     s = data[6]
#     if s == key1:# and key1 in line:# and s == key2:
#         print(line)
#         data[5:] = [k,s]
#         i += 1
#     line = '\t'.join(data)
#     print(line,file=outfile)
# outfile.close()
# print(i)

#### animals:
    
# def get_unique_species(table,list2):
#     list2 = [x.split('(')[0] if '(' in x else x for x in list2]
#     species = {}
#     for u in table:
#         list1 = table[u]
#         sp = [x for x in list1 if x not in list2]
#         for s in sp:
#             if s not in species:
#                 species[s] = []
#             species[s].append(u)
#     new_sp = [s + '('+','.join(species[s])+')' for s in species]
#     print('unique', len(species.keys()))
#     return new_sp
    
# inFile = '/home/ijulcach/projects/ldo_project/data/animal_entity/animal_uberon.txt'
# outname = inFile+'2'

# key = "bladder"
# k = 'bladder organ'

# juntos, u  = {}, ''
# for line in open(inFile):
#     line = line.strip()
#     data = line.split('\t')
#     n = data[1]
#     if key in n and 'gallbladder' not in n:
#         if n == k:
#             u = data[0]
#         else:
#             juntos[data[0]] = data[3].split('; ')
#     else:
#         if n == k:
#             u = data[0]

# print(len(juntos), juntos)
# print('join',u)
# outfile = open(outname, 'w')
# for line in open(inFile):
#     line = line.strip()
#     data = line.split('\t')
#     if data[0] in juntos:
#         pass
#     else:
#         if data[0] == u:
#             species = data[3].split('; ')
#             new = get_unique_species(juntos,species)
#             if len(new)!=0:
#                 data[2] = str(int(data[2])+len(new))
#                 data[3] += '; '+'; '.join(new)
#                 line = '\t'.join(data)
#                 if len(data) == 4:
#                     line += '\t'+'; '.join(list(juntos.keys()))
#                 else:
#                     line += '; '.join(list(juntos.keys()))
#             else:
#                 if len(data) == 4:
#                     line += '\t'+'; '.join(list(juntos.keys()))
#                 else:
#                     line += '; '.join(list(juntos.keys()))
#         print(line,file=outfile)
# outfile.close()
    


############## Ploting distribution of mean and median of leaves
# files = glob.glob('/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/Plant_expression/Exp_matrices/*/*.fltMatrix.txt')

# for f in files:
#     n = f.split('.')[0]
#     sp = n.split('/')[-1]
#     print(n)
#     dades = [[],[]]
#     for line in open(f):
#         line = line.strip()
#         data = line.split('\t')
#         if data[0] == 'genes' or data[0] == 'Snames':
#             pass
#         elif data[0] == 'Gnames':
#             pos = []
#             for j,e in enumerate(data[1:]):
#                 if e == 'leaf':
#                     pos.append(j)
#         else:
#             tpm = [float(x) for x in data[1:]]
#             values = []
#             for p in pos:
#                 values.append(tpm[p])
#             dades[0].append(round(np.mean(values),2))
#             dades[1].append(round(np.median(values),2))

#     ax = plt.boxplot(dades)
#     plt.title(label=sp)
#     plt.savefig(n+'.leaf.mean_median.svg', bbox_inches='tight')
#     plt.show()

#################### Ploting the number of samples per species

# table, organs = {}, []
# for p in expmat:
#     sp = p.split('/')[-2]
#     print(sp)
#     f = p+sp+'.fltMatrix.txt'
#     for line in open(f):
#         line = line.strip()
#         data = line.split('\t')
#         if data[0] == 'Gnames':
#             samples = data[1:]
#             table[sp] = {x:samples.count(x) for x in set(samples)}
#             organs += list(set(samples))
# organs = list(set(organs))
# print('Number of samples:', len(organs))
# outfile = open(outsampTable, 'w')
# head = 'Species\tMnemonic'+'\t'+'\t'.join(organs)
# print(head,file=outfile)
# for sp in table:
#     string = sp2name[sp]+'\t'+sp
#     for o in organs:
#         if o in table[sp]:
#             string += '\t'+str(table[sp][o])
#         else:
#             string += '\t0'
#     print(string,file=outfile)
# outfile.close()

# df = pd.read_csv(outsampTable, sep='\t', header=0, index_col=0)
# fig, ax = plt.subplots(figsize=(15,10)) 
# ax = sns.heatmap(df, annot=True, cmap='Blues', vmax=10)
# plt.savefig(outfigs, bbox_inches='tight')
# plt.show()


######################## SPM supplementary plot
# fig, axes = plt.subplots(nrows=5, ncols=4, figsize=(10, 14), sharex=True, sharey=True)
# axes = axes.flatten()
# for i,p in enumerate(expmat):
#     sp = p.split('/')[-2]
#     species = sp2name[sp]
#     f = p+sp+'.Gnames.spm'
#     num = []
#     for line in open(f):
#         line = line.strip()
#         data = line.split('\t')
#         if data[0]=='#gene' or data[0] == 'Samples' or data[0]=='Gnames' or data[0]=='Snames': 
#             pass
#         else:
#             for d in data[1:]:
#                 num.append(float(d))
#     num_sort = sorted(num, reverse=True)
#     x = int(float(5*len(num))/float(100))
#     n = num_sort[x]
#     print(sp, n)
#     ### plot
#     df = pd.DataFrame(num_sort)
#     c = 'lightsteelblue'
#     df.hist(bins=99, grid=False, color = c, ec=c, ax=axes[i])
#     axes[i].set_title(species, fontsize=12, style='italic')
#     plt.yscale('log')
#     axes[i].axvline(x=n, ymin=0, ymax=1, color='#ff0000ff', linestyle='dashed', linewidth=1.5)
#     x = round(n,1)
#     ymin, ymax = plt.ylim()
#     axes[i].text(x + 0.08, ymax/4, f'x: {x}', color='red')#, rotation=90)
#     axes[i].tick_params(axis='both', which='both', length=0)
    
# fig.text(0.02, 0.5, 'counts (log-scale)', va='center', rotation='vertical', fontsize=15)
# fig.text(0.5, 0.02, 'SPM', ha='center', fontsize=15)
# plt.savefig(outfig13, bbox_inches='tight')
# plt.show()
    
