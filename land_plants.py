#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 20:11:54 2024

@author: ijulcach
"""

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import ete3

def get_species(treeFile):
    t = ete3.Tree(treeFile)
    names = []
    for leaf in t:
        names.append(leaf.name)
    return names

path = '/home/ijulcach/projects/Land_Plants/'
plotPath = path +'plot_results/'
numbersFile = '/home/ijulcach/projects/Land_Plants/numbers_genomes_oma.tsv.csv'
treemajorLin = path+'major_lineages.nw'
treeangios = path+'angiosperms_tree.nw'
faoFile = path+'Fao_data/FAOSTAT_data_en_3-22-2024_2000_2022.csv'
##### OutputFiles
outfig1 = plotPath+'major_lineages_species.svg'
outfig2 = plotPath+'major_lineages_genomes.svg'

### Format FAO data
# i = 0
# for line in open(faoFile):
#     line = line.strip()
#     data = line.split('","')
#     if 'Domain Code' in data[0]:
#         print(data)
#     else:
#         food, year = data[7], data[9]
#         print(food,year)
#         print(data[10:])
    


## get the order of species to be plotted:
majorLin = get_species(treemajorLin)
angiosperms = get_species(treeangios)
y = 'major_lin'# 'clades' #'major_lin'
## load data
df = pd.read_csv(numbersFile, sep='\t')
# print(df.columns.values)
# print(df['genome_avail'].sum()) ## number of species
# ## select major lineages
df2 = df.drop(['clades'], axis=1) 
df2 = df2.groupby(by=['major_lin']).sum()
df2 = df2.reindex(index = majorLin) # reindex based in the order of species
# ## select clades of angiosperms
# df2 = df.drop(['major_lin'], axis=1)
# df2 = df2.set_index('clades')
# noangios = [x for x in list(df['clades'].values) if x not in angiosperms]
# df2 = df2.drop(noangios)
# print('number of angiosperms', df2['number_species'].sum())
# ## calculate percentages
df2['num_per'] = df2['number_species']*100/df2['number_species'].sum()
df2['geno_per'] = 100*((df2['num_per']/100)*(df2['genome_avail']/df2['number_species']))
# print('% of species sequenced',df2['genome_avail'].sum()*100/df2['number_species'].sum())
# ## plot 1
plt.figure(figsize=(4,6)) ##4,6 ###for major, 4,8
# # ax = sns.barplot(df2, x="num_per", y=y, color='#98D8AA')
# # # ax.set_xlim(0,100)
# # ax = sns.barplot(df2, x='geno_per', y=y, color='#4D96FF')
# # ax.set_xlabel('% Species')
# # ax.set_ylabel("")
# # plt.savefig(outfig1, bbox_inches='tight')
# # plt.show()

# ## plot 2
# df2.loc['Ceratophyllales','genome_avail'] =0
# df2.loc['Amborellales','genome_avail'] =0
df2['geno_tper'] = 100*df2['genome_avail']/df2['number_species']
df2['geno_tper'] = df2['geno_tper'].where(df2['geno_tper'] <= 3, 3)

df2['oma_per'] = 100*((df2['geno_tper']/100)*(df2['oma_species']/df2['genome_avail']))
df2['oma_nper'] = 100*((df2['geno_tper']/100)*(df2['future_oma']/df2['genome_avail']))
ax = sns.barplot(df2, x='geno_tper', y=y, color='#4D96FF')
ax = sns.barplot(df2, x='oma_nper', y=y, color='#6BCB77')
ax = sns.barplot(df2, x='oma_per', y=y, color='#FFD93D')
ax.set_xticklabels(ax.get_xticks().astype(int))
ax.set_xlabel('% Sequenced genomes')
ax.set_ylabel("")
plt.savefig(outfig2, bbox_inches='tight')
plt.show()





