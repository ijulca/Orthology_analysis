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
from scipy.stats import kruskal
from scikit_posthocs import posthoc_dunn
from statsmodels.sandbox.stats.multicomp import multipletests
import matplotlib.colors
import math
import ete3

def get_keys_hf5(inFile):
    with pd.HDFStore(inFile) as hdf:
        # This prints a list of all group names:
        keys = hdf.keys()
    return keys

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
        table[data[0]] = data[1]#.split(';')[-1].strip()
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

def format_pval(p):
    if p < 0.01:
        sn = np.format_float_scientific(p, precision=1, exp_digits=1)+"**"
    elif p < 0.05:
        sn = np.format_float_scientific(p, precision=1, exp_digits=1)+"*"
    else:
        sn = str(round(p,3))
    return sn

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
            p = gmo.get_wilcoxon_rank_sum(d1,d2)
            pval.append(p)
        p_adjusted = gmo.get_p_adj(pval)
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

###########################################
###########################################
######### PCC and TAU functions    ########
###########################################
###########################################

def transform_p(p):
    if p < 0.01:
        sn = np.format_float_scientific(p, precision=1, exp_digits=1)+"**"
    elif p < 0.05:
        sn = np.format_float_scientific(p, precision=1, exp_digits=1)+"*"
    else:
        sn = str(round(p,3))
    return sn

def get_data_h5(inFile, values, tag):
    keys = get_keys_hf5(inFile)
    for k in keys:
        # print(k)
        df = pd.read_hdf(inFile, key=k)
        df['tauV'] = df['mdo_tau'] - df['ldo_tau']
        # print(df.columns.values)
        df2 = df.loc[(df['ldo_category'] == 'normal') & (df['mdo_category'] == 'normal')]
        df3 = df.loc[(df['ldo_category'] == 'normal') & (df['mdo_category'] == 'long')]
        if tag == 'pcc':
            val_n = list(df2['pcc_stat'].values)
            val_l = list(df3['pcc_stat'].values)
        elif tag == 'tau':
            val_n = list(df2['tauV'].values)
            val_l = list(df3['tauV'].values)
        values[0]+= [x for x in val_n if str(x) !='nan']
        values[1]+= [x for x in val_l if str(x) !='nan']
    return values

def species_pcc_analisis(inFile, table, tag):
    keys = get_keys_hf5(inFile)
    for k in keys:
        n = k.replace('/','')
        df = pd.read_hdf(inFile, key=k)
        df2 = df.loc[(df['ldo_category'] == 'normal') & (df['mdo_category'] == 'normal')]
        df3 = df.loc[(df['ldo_category'] == 'normal') & (df['mdo_category'] == 'long')]
        if tag == 'pcc':
            val_n = list(df2['pcc_stat'].values)
            val_l = list(df3['pcc_stat'].values)
        elif tag == 'tau':
            df2['tauV'] = df2['mdo_tau'] - df2['ldo_tau']
            df3['tauV'] = df3['mdo_tau'] - df3['ldo_tau']
            val_n = list(df2['tauV'].values)
            val_l = list(df3['tauV'].values)
        table[n] =[[x for x in val_n if str(x) !='nan']]
        table[n].append([x for x in val_l if str(x) !='nan'])
    return table

def nodes_pcc_analisis(inFile, table, name, tag):
    table.update({'internal_'+name:[[],[]], 'terminal_'+name:[[],[]]})
    keys = get_keys_hf5(inFile)
    for k in keys:
        n = k.replace('/','')
        df = pd.read_hdf(inFile, key=k)
        df['tauV'] = df['mdo_tau'] - df['ldo_tau']
        df2 = df.loc[(df['ldo_category'] == 'normal') & (df['mdo_category'] == 'normal')]
        df3 = df.loc[(df['ldo_category'] == 'normal') & (df['mdo_category'] == 'long')]
        dfnt = df2.loc[(df['sp_tree_head'] == n)] ## select terminal node normal
        dfni = df2.loc[(df['sp_tree_head'] != n)]
        dflt = df3.loc[(df['sp_tree_head'] == n)] ## select terminal node n-long
        dfli = df3.loc[(df['sp_tree_head'] != n)]
        if tag == 'pcc':
            vnt = list(dfnt['pcc_stat'].values)
            vni = list(dfni['pcc_stat'].values)
            vlt = list(dflt['pcc_stat'].values)
            vli = list(dfli['pcc_stat'].values)
        elif tag == 'tau':
            vnt = list(dfnt['tauV'].values)
            vni = list(dfni['tauV'].values)
            vlt = list(dflt['tauV'].values)
            vli = list(dfli['tauV'].values)
        table['terminal_'+name][0] += [x for x in vnt  if str(x) !='nan']
        table['terminal_'+name][1] += [x for x in vlt if str(x) !='nan']
        table['internal_'+name][0] += [x for x in vni if str(x) !='nan']
        table['internal_'+name][1] += [x for x in vli if str(x) !='nan']
    return table

def load_data_corr_out(inFile, outname,name):
    keys = get_keys_hf5(inFile)
    outfile = open(outname, 'w')
    cat = ['n1','n2','n3','l1']
    print('sp1\tsp2\tgene\tout_gene\tlineage\tcat\tval', file=outfile)
    for k in keys:
        n = k.replace('/','')
        df = pd.read_hdf(inFile, key=k)
        df = df.loc[(df['out_ldo_pcc_stat'].notnull()) & (df['out_mdo_pcc_stat'].notnull())] ### remove the nan values
        df2 = df.loc[(df['ldo_category'] == 'normal') & (df['mdo_category'] == 'normal')]
        df3 = df.loc[(df['ldo_category'] == 'normal') & (df['mdo_category'] == 'long')]
        species = set(list(df['outgroup_species'].values))
        for s in species:
            df2s = df2.loc[(df2['outgroup_species'] == s)]
            df3s = df3.loc[(df3['outgroup_species'] == s)]
            pcc_n1 = df2s['out_ldo_pcc_stat'].values
            pcc_n2 = df2s['out_mdo_pcc_stat'].values
            pcc_l1 = df3s['out_ldo_pcc_stat'].values
            pcc_l2 = df3s['out_mdo_pcc_stat'].values
            if len(pcc_n1)>=50 and len(pcc_n2)>=50 and len(pcc_l1)>=50 and len(pcc_l2)>=50: ### filter minimum amount of data in each catergory
                if len(set(pcc_n1))>2 and len(set(pcc_n2))>2 and len(set(pcc_l1))>2 and len(set(pcc_l2))>2: ## filter rare values
                        dades = [pcc_n1, pcc_n2, pcc_l1, pcc_l2]
                        for i,e in enumerate(dades):
                            c = cat[i]
                            for j,v in enumerate(e):
                                if c == 'n1':
                                    string = n+'\t'+s+'\t'+df2s['ldo_gene'].values[j]+'\t'+df2s['out_gene'].values[j]
                                elif c =='n2':
                                    string = n+'\t'+s+'\t'+df2s['mdo_gene'].values[j]+'\t'+df2s['out_gene'].values[j]
                                elif c == 'n3':
                                    string = n+'\t'+s+'\t'+df3s['ldo_gene'].values[j]+'\t'+df3s['out_gene'].values[j]
                                elif c == 'l1':
                                    string = n+'\t'+s+'\t'+df3s['mdo_gene'].values[j]+'\t'+df3s['out_gene'].values[j]
                                string += '\t'+name+'\t'+c+'\t'+str(v) 
                                print(string,file=outfile)
    outfile.close()

def load_data_tau_out(inFile, outname, name):
    keys = get_keys_hf5(inFile)
    outfile = open(outname, 'w')
    cat = ['n1','n2','n3','l1']
    print('sp1\tsp2\tgene\tout_gene\tlineage\tcat\tval', file=outfile)
    for k in keys:
        n = k.replace('/','')
        df = pd.read_hdf(inFile, key=k)
        df['tau_l-o'] = df['ldo_tau'] - df['out_tau'] ## difference of tau
        df['tau_m-o'] = df['mdo_tau'] - df['out_tau']
        df = df.loc[(df['tau_l-o'].notnull()) & (df['tau_m-o'].notnull())] ### remove the nan values
        df2 = df.loc[(df['ldo_category'] == 'normal') & (df['mdo_category'] == 'normal')]
        df3 = df.loc[(df['ldo_category'] == 'normal') & (df['mdo_category'] == 'long')]
        species = set(list(df['outgroup_species'].values))
        for s in species:
            df2s = df2.loc[(df2['outgroup_species'] == s)]
            df3s = df3.loc[(df3['outgroup_species'] == s)]
            tau_n1 = df2s['tau_l-o'].values
            tau_n2 = df2s['tau_m-o'].values
            tau_l1 = df3s['tau_l-o'].values
            tau_l2 = df3s['tau_m-o'].values
            if len(tau_n1)>=50 and len(tau_n2)>=50 and len(tau_l1)>=50 and len(tau_l2)>=50: ### filter minimum amount of data in each catergory
                if len(set(tau_n1))>2 and len(set(tau_n2))>2 and len(set(tau_l1))>2 and len(set(tau_l2))>2: ## filter rare values
                        dades = [tau_n1, tau_n2, tau_l1, tau_l2]
                        for i,e in enumerate(dades):
                            c = cat[i]
                            for j,v in enumerate(e):
                                if c == 'n1':
                                    string = n+'\t'+s+'\t'+df2s['ldo_gene'].values[j]+'\t'+df2s['out_gene'].values[j]
                                elif c =='n2':
                                    string = n+'\t'+s+'\t'+df2s['mdo_gene'].values[j]+'\t'+df2s['out_gene'].values[j]
                                elif c == 'n3':
                                    string = n+'\t'+s+'\t'+df3s['ldo_gene'].values[j]+'\t'+df3s['out_gene'].values[j]
                                elif c == 'l1':
                                    string = n+'\t'+s+'\t'+df3s['mdo_gene'].values[j]+'\t'+df3s['out_gene'].values[j]
                                string += '\t'+name+'\t'+c+'\t'+str(v) 
                                print(string,file=outfile)
    outfile.close()

def add_time_dataframe(outtab1, plantTreetime, pref2names):
    table = {'_'.join(pref2names[x].split(' ')):x for x in pref2names} ##[:2] for plants
    t = ete3.Tree(plantTreetime, format=1)
    for leaf in t:
        leaf.name = table[leaf.name]
    outfile = open(outtab1+'_temp','w')
    for line in open(outtab1):
        line = line.strip()
        data = line.split('\t')
        if data[0] == 'sp1':
            line += '\ttime'
        else:
            n = t.get_common_ancestor([data[0], data[1]])
            node = t&data[0]
            i = 0
            while node != n:
                i += node.dist
                node = node.up
            line += '\t' +str(i)
        print(line, file=outfile)
    outfile.close()
        
def load_major_lin(inFile):
    table = {}
    for line in open(inFile):
        line = line.strip()
        data = line.split('\t')
        if data[0] not in table:
            table[data[0]] = set()
        table[data[0]].add(data[0])
        for e in data[2].split('; '):
            table[data[0]].add(e)
    return table

def change2float(x):
    try:
        return float(x)
    except:
        return np.nan
    
def get_structure_data(dupFile, foFiles, outfigname, majorlinFile):
    majorlin = load_major_lin(majorlinFile)
    categories = ['normal-normal','normal-long']
    df = pd.read_csv(dupFile,sep='\t')
    df['genes'] = df['ldo_gene']+'-'+df['mdo_gene']
    df['lddt'] = np.nan
    for f in foFiles:
        keys = get_keys_hf5(f)
        for k in keys:
            stframe = pd.read_hdf(f, key=k)
            stframe["lddt"] = stframe["lddt"].apply(change2float) ### change to float
            stframe['genes'] = stframe['gene1']+'-'+stframe['gene2']
            df.loc[df.genes.isin(stframe.genes), ['lddt']] = stframe['lddt'].values
    # df2 = df.loc[(df['ldo_category'] == 'normal') & (df['mdo_category'] == 'normal')]
    # df3 = df.loc[(df['ldo_category'] == 'normal') & (df['mdo_category'] == 'long')]
    # table = [[x for x in list(df2['lddt'].values) if str(x)!='nan'], [x for x in list(df3['lddt'].values) if str(x)!='nan']]
    # ax = sns.boxplot(data=table, saturation=0.7, fliersize=0)#, linecolor=color2)
    # plt.xticks([0,1], categories)
    # plt.ylabel('LDDT')
    # plt.savefig(outfigname+'_all.svg', bbox_inches='tight')
    # print(len(table[0]), len(table[1]))
    # print(np.median(table[0]), np.median(table[1]), np.average(table[0]), np.average(table[1]))
    table = {}
    for lin in majorlin:
        df1 = df[df['sp_tree_head'].isin(majorlin[lin])]
        df2 = df1.loc[(df1['ldo_category'] == 'normal') & (df1['mdo_category'] == 'normal')]
        normal = [x for x in df2['lddt'].values if str(x)!='nan']
        df3 = df1.loc[(df1['ldo_category'] == 'normal') & (df1['mdo_category'] == 'long')]
        long = [x for x in df3['lddt'].values if str(x)!='nan']
        dades = [normal, long]
        table[lin] = dades
    plot_subplots_barplot(table, 4, 3, list(majorlin.keys()), list(majorlin.keys()), outfigname+'_lineages.svg')

def load_gene_panther_uniprot(inFile):
    table = {}
    for line in open(inFile):
        line = line.strip()
        data = line.split('\t')[1].split('|')
        sp,g = data[0], data[2].split('ProtKB=')[1]
        table[g] = sp
    return table

def get_string_gene(g1,s,c, comu, panterg, lddt, outfile):
    for p in comu:
        v = lddt[g1][p]
        if v != 'nan' and v!='':
            s2 = panterg[p]
            string = s+'\t'+s2+'\t'+g1+'\t'+p+'\t'+c+'\t'+str(v)
            print(string, file=outfile)

def get_structure_out(dupFile, foFiles, pref2names, panterGenesFile, outTable):
    panterg = load_gene_panther_uniprot(panterGenesFile)
    lddt = {}
    for f in foFiles:
        keys = get_keys_hf5(f)
        for k in keys:
            df = pd.read_hdf(f, key=k)
            for index, row in df.iterrows():
                g1,g2 = row['gene1'], row['gene2']
                v = row['lddt']
                if g1 not in lddt:
                    lddt[g1] = {}
                lddt[g1][g2] = v
    outfile = open(outTable, 'w')
    print('sp1\tsp2\tg1\tg2\tcat\tval', file=outfile)
    df = pd.read_csv(dupFile,sep='\t')
    for index, row in df.iterrows():
        c1,c2 = row['ldo_category'],row['mdo_category']
        g1,g2 = row['ldo_gene'], row['mdo_gene']
        s = row['species']
        if g1 in lddt and g2 in lddt:
            gen1,gen2 = list(lddt[g1].keys()),list(lddt[g2].keys())
            comu = list(set(gen1) & set(gen2))
            if len(comu) != 0:
                if c1 == 'normal' and c2 == 'normal':
                    get_string_gene(g1,s,'n1', comu, panterg, lddt, outfile)
                    get_string_gene(g2,s,'n2', comu, panterg, lddt, outfile)
                elif c1 == 'normal' and c2 == 'long':
                    get_string_gene(g1,s,'n3', comu, panterg, lddt, outfile)
                    get_string_gene(g1,s,'l1', comu, panterg, lddt, outfile)
    outfile.close()
        

def get_tables_exp(inFile, outfile):
    keys = get_keys_hf5(inFile)
    for k in keys:
        df = pd.read_hdf(inFile, key=k)
        values = df.columns.values
        n = k.replace('/','')
        string = n +'\t'+str(len(values))
        print(string,file=outfile)
        # for e in values:
        #     string = n+'\t'+e[1]+'\t'+e[0]
        #     print(string, file=outfile)

#### Plots
def plot_subplots_barplot(table, rows, cols, samples, names, outname):
    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(15, 15), sharex=True, sharey=True) ##Fig1: figsize=(3, 3), (10,4), (15,15)
    axes = axes.flatten()
    categories = ['normal-normal','normal-long']
    i,j = 0,1.1
    for x,tag in enumerate(samples):
        dades = table[tag]
        print(tag, np.median(dades[0]), np.median(dades[1]))
        p = gmo.get_wilcoxon_rank_sum(dades[0], dades[1])
        sn = gmo.format_pval(p)
        sns.boxplot(data=dades, saturation=0.7, ax=axes[x], fliersize=0, palette="PiYG")#, PiYG,PRGn linecolor=color2)
        x1,x2 = 0,i+1
        y, h, col, z = j+0.04, 0.05+i/4, 'k', 0.03
        axes[x].plot([x1, x1, x2, x2], [y-z, y+h-z, y+h-z, y-z], lw=0.5, c=col)
        axes[x].text((x1+x2)*.5, y+h-z, sn, ha='center', va='bottom', color=col, fontsize=8)
        name = names[x]
        axes[x].set_title(name)#, style='italic')#, size=20)
        axes[x].set_ylim(-1, 1.5) ## Fig1 1.6
        ticks_to_keep = axes[x].get_yticks()[:5] ### remove some yticks
        axes[x].set_yticks(ticks_to_keep)
        axes[x].set_ylabel('PCC')
    plt.xticks([0,1], categories)#, rotation=45, ha='right')
    # plt.savefig(outname, bbox_inches='tight')
    plt.show()

def plot_subplots_violinplot(table, rows, cols, samples, names, outname):
    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(15, 15), sharex=True, sharey=True) ##Fig1: figsize=(4, 4),(10,4), (15,15)
    axes = axes.flatten()
    categories = ['normal-normal','normal-long']
    for x,tag in enumerate(samples):
        dades = table[tag]
        p = gmo.get_wilcoxon_rank_sum(dades[0], dades[1])
        sn = gmo.format_pval(p)
        if '*' not in sn:
            print(tag, sn)
        sns.violinplot(data=dades, palette="PRGn", ax=axes[x], split=True, density_norm="width")
        name = names[x]
        axes[x].set_title(name)#, style='italic')#, size=20)
        axes[x].set_ylim(-1, 1.5) ## Fig1 1.6
        ticks_to_keep = axes[x].get_yticks()[:5] ### remove some yticks
        axes[x].set_yticks(ticks_to_keep)
        axes[x].set_ylabel('delta TAU')
    plt.xticks([0,1], categories)#, rotation=45, ha='right')
    plt.savefig(outname, bbox_inches='tight')
    plt.show()

def plot_subplots_line_pcc(inFile, pref2names, outfigname):
    df = pd.read_csv(inFile, sep='\t')
    ordersp = list(set(df[df['lineage']=='Plants']['sp1'].values))+list(set(df[df['lineage']=='Animals']['sp1'].values))
    species = {x:0 for x in set(list(df['sp1'].values))}
    for tag in species:
        df2 = df.loc[(df['sp1'] == tag)]
        dades = set(list(df2['sp2'].values))
        species[tag] = len(dades)
    speciesflt = [x for x in species if species[x]>=3] ### remove species with few outgroups
    print('number of initial species:',len(species), 'filtered:',len(speciesflt), [x for x in species if x not in speciesflt])
    species = [x for x in ordersp if x in speciesflt]
    fig, axes = plt.subplots(nrows=7, ncols=4, figsize=(16,14), sharex=True, sharey=True) ##Fig1: Plants: 4,4 figsize=(15, 10), 3,4
    axes = axes.flatten()
    for x,tag in enumerate(species):
        name = ' '.join(pref2names[tag].split(' ')[:2])
        name = name.split(' ')[0][0]+'. '+name.split(' ')[1]
        df2 = df.loc[(df['sp1'] == tag)]
        sns.lineplot(data=df2, x="time", y="val", hue="cat", markers=True, ax=axes[x], errorbar=('ci', 95),
                      palette=sns.color_palette("PiYG", 4)) ## PCC: PiYG, PRGn
        axes[x].set_title(name, style='italic')
        axes[x].set_ylabel('PCC')  ## PCC  ## delta TAU
        axes[x].set_xlabel('Time (MYA)')
    plt.savefig(outfigname, bbox_inches='tight')
    plt.show()

def plot_subplots_line_pcc2(inFile, pref2names, outfigname):
    df = pd.read_csv(inFile, sep='\t')
    lin = set(df['lineage'].values)
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(5,4), sharex=True, sharey=True) ##Fig1: Plants: 4,4 figsize=(15, 10), 3,4
    axes = axes.flatten()
    for x,tag in enumerate(lin):
        df2 = df.loc[(df['lineage'] == tag)]
        sns.lineplot(data=df2, x="time", y="val", hue="cat", markers=True, ax=axes[x], errorbar=('ci', 95),
                      palette=sns.color_palette("PRGn", 4)) ## PCC: PiYG ##PRGn
        axes[x].set_title(tag)
    plt.savefig(outfigname, bbox_inches='tight')
    plt.show()
    

#################
#### inFiles ####
#################
path = '/home/ijulcach/projects/ldo_project/results/'
pathPlots = path + 'plots/'
pathTables = path + 'tables/'
famrateFile = pathTables+'fitted_family_info.tsv'
speciesFile = '/home/ijulcach/projects/ldo_project/data/panther-18.0/species_tree.nhx'

pccFilePlant = path +'correlation/plant_specific.h5'
pccFileAnimal = path +'correlation/animal_maincat.h5'
prefFile = pathTables+ 'nemo2name.txt'
outFilePlant = path +'outgroup_correlation/plant_specific.h5'
outFileAnimal = path +'outgroup_correlation/animal_maincat.h5'
plantTreetime = pathTables + 'plant_species_tree.nwk'
animalTreetime = pathTables +'animal_species.nwk'
FoldsFiles = '/home/ijulcach/projects/ldo_project/structure/res_foldseek/'
UalgFiles = '/home/ijulcach/projects/ldo_project/structure/res/'
allDupFile = '/home/ijulcach/projects/ldo_project/results/pairwise_tests.tsv'
majorlinfile = pathTables+'major_lineages_sp.txt'
panterTimeFile = pathTables+'species_panther.nwk'
panterGenesFile = '/home/ijulcach/projects/ldo_project/data/panther-18.0/Panther.genesNames.txt'
expDataPath = '/home/ijulcach/projects/ldo_project/results/expr/'

##################
#### outfiles ####
##################

##### Duplication analysis
majorLinFile = pathTables+'major_lineages_sp.txt'
node2catFile = pathTables+'nodes2categories.tsv'
node2catLinFile = pathTables+'nodes2categories.MajorLin.tsv'
node2catLin_perFile = pathTables+'nodes2categories.MajorLin_per.tsv'
node2catLin_perCatFile = pathTables+'nodes2categories.MajorLin_per_cat.tsv'
nodesRanksFile = pathTables+'nodes2categories_numsranksMerge.tsv'

dupFigure1 = pathPlots +'histogram_familyrates.svg'  ## Supplementary Figure 1
dupFigure2 = pathPlots + 'bar_dups_categories.svg' ## Panel Figure 2
dupFigure3 = pathPlots + 'internl_terminal_dup.svg'

##### PCC and TAU
outfig1 = pathPlots+'pcc_all.svg'
outfig2 = pathPlots+'pcc_species.svg'
outfig3 = pathPlots+'pcc_nodes.svg'
outfig4 = pathPlots+'tau_all.svg'
outfig5 = pathPlots+'tau_species.svg'
outfig6 = pathPlots+'tau_nodes.svg'
outfig7 = pathPlots+'pcc_outgroup_plant_animal.svg'
outfig8 = pathPlots+'pcc_outgroup_plant_animal_all.svg'
outfig9 = pathPlots+'tau_outgroup_plant_animal.svg'
outfig10 = pathPlots+'tau_outgroup_plant_animal_all.svg'
outfig11 = pathPlots+'structure_lddt'
outfig12 = pathPlots+'structure_lddt_outgroup.svg'

outtab1 = pathTables+'outgroup_pcc_analysis_plant.tsv'
outtab2 = pathTables+'outgroup_pcc_analysis_animal.tsv'
outtab3 = pathTables+'outgroup_tau_analysis_plant.tsv'
outtab4 = pathTables+'outgroup_tau_analysis_animal.tsv'
outtab5 = pathTables+'outgroup_structure.tsv'
outtab6 = pathTables+'expression_samples.tsv'


#######################
#### PANTHER TREES #### 
#######################
### infiles and outfiles

# treepath = '/home/ijulcach/projects/ldo_project/data/panther-18.0/trees/'
# 
# expPlantdata = '/home/ijulcach/projects/ldo_project/data/Plant_expression/Samples/'
# taxaFile = '/home/ijulcach/projects/Land_Plants/nemo2taxa.txt'

# inFile = pathTables+'pairwise_tests.tsv'
# linFile = pathTables+'major_lineages_sp.txt'
# expPath = '/home/ijulcach/projects/ldo_project/data/Plant_expression/Exp_matrices/'
# # outfile1 = pathTables+'plant_exp_pcc_scaler.tsv'
# # outfig1 = pathPlots+'plant_exp_pcc_scaler.svg'
# # outfile1 = pathTables+'plant_exp_pcc_noscaler.tsv'
# # outfig1 = pathPlots+'plant_exp_pcc_noscaler.svg'
# # outfile1 = pathTables+'plant_exp_pcc_zscore.tsv' 
# # outfig1 = pathPlots+'plant_exp_pcc_zscore.svg' 
# # outfile1 = pathTables+'plant_exp_pcc_zscoremedian.tsv'
# # outfig1 = pathPlots+'plant_exp_pcc_zscoremedian.svg'
# # outfile1 = pathTables+'plant_exp_pcc_Snames.zscoreLog.tsv'
# # outfig1 = pathPlots+'plant_exp_pcc_Snames.zscoreLog.svg'
# outfile1 = pathTables+'plant_exp_pcc_zscoreLog.tsv'  ### I use this
# outfig1 = pathPlots+'plant_exp_pcc_zscoreLog.svg'    ### I use this
# outfig2 = pathPlots+'plant_exp_pcc_zscoreLog_nodesS.svg'
# outfig3 = pathPlots+'plant_exp_pcc_zscoreLog_nodesI.svg'
# outfig4 = pathPlots +'plant_exp_pcc_zscoreLog_sp.svg'
# outfig5 = pathPlots +'plant_exp_pcc_zscoreLog_spS.svg'
# outfig6 = pathPlots +'plant_exp_pcc_zscoreLog_spI.svg'

# files = glob.glob('/home/ijulcach/projects/ldo_project/data/Plant_expression/Exp_matrices/*/*.Gnames.sample.txt')
# outname1 = pathTables+'plant_sampleSPM.tsv'
# outname2 = pathTables + 'plant_sampleSPM_matrix.tsv'
# outfig7 = pathPlots+'plant_samplesSPM.svg'
# outfig8 = pathPlots+'plant_samplesSPMS.svg'
# outfig9 = pathPlots+'plant_samplesSPMI.svg'
# outfig10 = pathPlots+'plant_samplesSPM_sp.svg'
# outfig11 = pathPlots+'plant_samplesSPM_spS.svg'
# outfig12 = pathPlots+'plant_samplesSPM_spI.svg'
# outfig13 = pathPlots+'supplementary_spm_plants.svg'
# outsampTable = pathTables+'sp2number_samples.tsv'
# outfigs = pathPlots+'sp2number_samples.svg'


##############################
#### Duplication analysis ####
##############################
categories = ['normal-normal', 'short-short', 'long-long', 'normal-long','normal-short', 'short-long']

##### family-specific rate ### Supplementary Figure 1
# df = pd.read_csv(famrateFile, sep='\t', header=0)
# ax = sns.histplot(df,x='rate', bins=200)
# ax.set_title('Histogram of family-specific factor')
# ax.set_xlabel('Family-specific rate')
# ylims = ax.get_ylim()
# plt.vlines(1, ylims[0]-1000, ylims[1]+1000, color='red', linestyle='--')#, ymin=0, ymax=)\n",
# ax.set_ylim(*ylims)
# plt.savefig(dupFigure1, bbox_inches='tight')
# plt.show()
# print(df['rate'].mean())

###### Species tree, Creating duplication files
# lineages = ['Protostomia', 'Deuterostomia', 'Fungi', 'Amoebozoa', 'Excavates',
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
# outfile = open(majorLinFile,'w')
# for e in lineages:
#     string = e+'\t'+str(len(lineages[e]))+'\t'+'; '.join(lineages[e])
#     print(string,file=outfile)
#     i += len(lineages[e])
# print(i)
# outfile.close()

# internal_nodes = ['Bilateria', 'Eumetazoa', 'Metazoa-Choanoflagellida',
#                   'Opisthokonts', 'Unikonts', 'Eukaryota', 'Bikonts',
#                   'SAR/HA_supergroup', 'Archaea-Eukaryota', 'LUCA']
# # ## Node 1: 'Bilateria'
# # ## Node 2: 'Eumetazoa'
# # ## Node 3: 'Metazoa-Choanoflagellida'
# # ## Node 4: 'Opisthokonts'
# # ## Node 5: 'Unikonts'
# # ## Node 6: 'Eukaryota'
# # ## Node 7: 'Bikonts'
# # ## Node 8: 'SAR/HA_supergroup'
# # ## Node 9: 'Archaea-Eukaryota'
# # ## Node 10: 'LUCA'

##### Get the categories tables

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
# outfile = open(node2catFile,'w')
# print('Taxa\t'+'\t'.join(categories), file=outfile)
# for tax in table:
#     string = tax
#     for e in categories:
#         string += '\t'+str(table[tax][e])
#     print(string, file=outfile)
# outfile.close()

# table = {x:[0.0]*len(categories) for x in lineages}
# table.update({x:[0.0]*len(categories) for x in internal_nodes})
# for line in open(node2catFile):
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
# outfile = open(node2catLinFile,'w')
# print('Taxa\t'+'\t'.join(categories), file=outfile)
# for tax in table:
#     string = tax + '\t' + '\t'.join([str(x) for x in table[tax]])
#     print(string, file=outfile)
# outfile.close()

###### Plot of the categories ## Figure 2
lineages = ['Protostomia', 'Deuterostomia','TRIAD','NEMVE','MONBE', 'Fungi', 'Amoebozoa', 'Excavates',
            'Viridiplantae','Alveolata-Stramenopiles','Archaea', 'Eubacteria']
lineages.reverse()
df = pd.read_csv(node2catLinFile, sep='\t', header=0)
df['total'] = df.sum(numeric_only=True, axis=1)
df['totalper']=df['total']*100/df['total'].sum()
df.to_csv(path_or_buf=node2catLin_perFile, sep='\t', header=True, index=False)
df = df[df['Taxa'].isin(lineages) == True] ## remove the ones in internal nodes
df.sort_values(by=['Taxa'], key=lambda column: column.map(lambda e: lineages.index(e)), inplace=True)

df2 = pd.DataFrame() ## dataframe of categories
for e in categories:
    df2[e] = df[e]*100/df['total']
df2['Taxa'] = df['Taxa']
df2.to_csv(path_or_buf=node2catLin_perCatFile, sep='\t', header=True, index=False)

palette = ['#3182bdff','#9ecae1ff','#deebf7ff','#e6550dff','#fdae6bff','#fee6cec4']
bar_width = 0.4
bar_positions = range(len(df))
fig, ax = plt.subplots()
ax = df2.plot(x='Taxa', y=categories, kind='barh', stacked=True, color=palette, alpha=.7, ax=ax,
              figsize=(6, 7.5), width=bar_width)

ax.barh([p + bar_width for p in bar_positions],df['totalper'], bar_width, color='#808080ff', label='duplications')
plt.xticks(fontsize=23)
plt.savefig(dupFigure2, bbox_inches='tight')
plt.show()


##### Plot of the categories per lineage per node, Figure 3
## Create the file
species = set()
lineages = {}
for line in open(majorLinFile):
    line = line.strip()
    data = line.split('\t')
    lineages[data[0]] = data[2].split('; ')
    for e in data[2].split('; '):
        if len(e) <=5 and e.isupper():
            species.add(e)
print('Number of species', len(species))

table = {}
for line in open(node2catFile):
    line = line.strip()
    data = line.split('\t')
    if data[0] != 'Taxa':
        table[data[0]] = [float(x) for x in data[1:]]
ranks = {}
for l in lineages:
    ranks[l] = {}
    ranks[l]['Internal'] = table[l]
    for tax in lineages[l]:
        if tax in species:
            key = 'Terminal'#+l[0]
        else:
            key = 'Internal'#+l[0] ### tax
        if key not in ranks[l]:
            ranks[l][key]=[0.0]*len(categories)
        if tax in table:
            nums = table[tax]
            ranks[l][key] = [x + y for x, y in zip(ranks[l][key], nums)]
outfile = open(nodesRanksFile,'w')
head = 'Major_Rank\tRank\t'+'\t'.join(categories)
print(head, file=outfile)
for l in ranks:
    for e in ranks[l]:
        string = l+'\t'+e+'\t'+'\t'.join([str(x) for x in ranks[l][e]])
        print(string,file=outfile)
outfile.close()

#### Create the plot Figure 3

lineages = ['Protostomia', 'Deuterostomia', 'Fungi', 'Amoebozoa', 'Excavates', ### we removed 'TRIAD', 'NEMVE','MONBE'
            'Viridiplantae','Alveolata-Stramenopiles','Archaea', 'Eubacteria']
lineages.reverse()
df = pd.read_csv(nodesRanksFile, sep='\t', header=0) 
df = df[df['Major_Rank'].isin(lineages) == True] ### we removed 'TRIAD', 'NEMVE','MONBE'

## calculate the total percentage of duplications
for index, row in df.iterrows():
    taxa = df['Major_Rank'].values[index]
    df2 = df[df['Major_Rank'] == taxa]
    total = df2[categories].sum()  # getting the total number of duplications per lineage
    num = row.values[2:]
    df.at[index, 'TotalPer'] = np.sum(num)*100/total.sum()
    sym = (row[['normal-normal','short-short','long-long']].sum())*100/np.sum(num)
    asym = (row[['normal-long','normal-short','short-long']].sum())*100/np.sum(num)
    df.at[index,'symmetric'] = sym
    df.at[index,'asymmetric'] = asym + sym

palettedup = ['#bdbdbd','#636363']
palettecata = ['#fdae6bff','#e6550dff']  #d95f0e
palettecats = ['#9ecae1ff','#3182bdff']  
fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, gridspec_kw={'height_ratios': [1, 4]})#, gridspec_kw={'height_ratios': [1, 1, 3, 3]})#, figsize=(16,14), sharey=True) ##Fig1: Plants: 4,4 figsize=(15, 10), 3,4
axes = axes.flatten()
sns.barplot(df, x="TotalPer", y="Major_Rank", hue="Rank", palette=palettedup, ax=axes[2], alpha=.7)
sns.barplot(df, x='asymmetric', y='Major_Rank', hue='Rank', palette=palettecata, legend=False, ax=axes[3], alpha=.7)
sns.barplot(df, x='symmetric', y='Major_Rank', hue='Rank', palette=['white','white'], legend=False, ax=axes[3], alpha=.7)
sns.barplot(df, x='symmetric', y='Major_Rank', hue='Rank', palette=palettecats, legend=False, ax=axes[3], alpha=.7)

dfper = df[['Rank','TotalPer']].groupby(by=['Rank']).mean().T
dfper = dfper[['Terminal','Internal']]
dfper.plot(kind='barh', color=palettedup, ax=axes[0], alpha=.7)
dfcat = df[['Rank','symmetric','asymmetric']].groupby(by=['Rank']).mean()
dfcat = dfcat.reindex(['Terminal','Internal'])
dfcat['asymmetric'] = 100 - dfcat['symmetric']
dfcat.plot(kind='barh', color=['#3182bdff','#e6550dff'], stacked=True, ax=axes[1], alpha=.7)
plt.savefig(dupFigure3, bbox_inches='tight') ##Supplementary Fig 3  

plt.show()
    


#### Getting the names of the gene expression matrices and panther
# expPlantdata = '/home/ijulcach/projects/ldo_project/Alex_analysis/ldo_project_2023/data/Plant_expression/Exp_matrices/'
# expPlantdata = '/home/ijulcach/projects/ldo_project/data/Animal_expression/'

# animals = ['XENTR', 'DANRE', 'CANLF']

# key = animals[2]
# genes = set()
# for line in open(expPlantdata+key+'/'+key+'.newMatrix.txt'): #'.cds.fa'):
#     line = line.strip()
#     # if '>' in line:
#     #     genes.add(line.split('>')[1])
#     data = line.split('\t')
#     if data[0] != 'genes':
#         genes.add(data[0])#.split('_')[0])

# # dades = {}
# # for line in open(expPlantdata+key+'/ensemble2rat.txt'): #'/AGI2uniprot.txt'):
# #     line = line.strip()
# #     data = line.split('\t')
# #     if data[1] not in dades:
# #         dades[data[1]] = []
# #     dades[data[1]].append(data[0])
# #     # dades[data[1].split('-')[0]] =data[0].split('.')[0]

# table = {}
# for line in open(expPlantdata+key+'/'+key+'.genesNames.txt'):
#     line = line.strip()
#     data = line.split('\t')
#     name = data[1] #'|'.join(data)
#     n = name.split('|')[1]#.split('=')[1]#.split('.')[0]
#     if 'Ensembl' in name:
#         n = n.split('=')[1].split('.')[0]
#     if n not in table:
#         table[n] = name
#     # if 'EnsemblGenome' in data[1]:
#     #     n = data[1].split('=')[1]
#     # n = data[2].split('=')[1]
#     # n = data[1].split('=')[1].replace('CISIN_','orange1.').replace('mg','m.g')
#         # table[n] = name

# # outfile = open(expPlantdata+key+'/'+key+'.conversion.panther.txt','w')
# i = 0
# for e in table:
#     # if e in dades:
#     #     for g in dades[e]:
#     if e not in genes:
#         print(e)
#         i +=1
#                 # print(e, table[e])
# #                 string = g+'\t'+table[e]
# #                 print(string,file=outfile)
# # outfile.close()
# print(key, len(table), i, i*100/len(table)) 

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

# color = ['#01befe','#ffdd00','#ff7d00','#ff006d','#adff02','#8f00ff']
# color2 = [matplotlib.colors.to_rgb(x) for x in color]

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






#############################################
#############################################
######### PCC and TAU analysis final ########
#############################################
#############################################

#### All PCC
# table = {'plant':[[],[]], 'animal':[[],[]]}
# table['plant'] = get_data_h5(pccFilePlant, table['plant'],'pcc')
# table['animal'] =get_data_h5(pccFileAnimal, table['animal'], 'pcc')
# samples = ['animal','plant']
# names = ['Animals', 'Plants']
# plot_subplots_barplot(table, 1, 2, samples, names, outfig1)
# table = {'plant':[[],[]], 'animal':[[],[]]}
# table['plant'] = get_data_h5(pccFilePlant, table['plant'],'tau')
# table['animal'] =get_data_h5(pccFileAnimal, table['animal'], 'tau')
# # plot_subplots_barplot(table, 1, 2, samples, names, outfig4)
# plot_subplots_violinplot(table, 1, 2, samples, names, outfig4)

# #### PCC species
# pref2names = load_taxa(prefFile)
# table = {}
# table = species_pcc_analisis(pccFilePlant, table, 'pcc')
# plants = list(table.keys())
# table = species_pcc_analisis(pccFileAnimal, table, 'pcc')
# species = plants+[x for x in list(table.keys()) if x not in plants]
# names = [pref2names[x] for x in species]
# names = [x.split(' ')[0][0]+'. '+x.split(' ')[1] for x in names] ## Cut Sci. names
# # plot_subplots_barplot(table, 6, 6, species, names, outfig2)
# table = {}
# table = species_pcc_analisis(pccFilePlant, table, 'tau')
# table = species_pcc_analisis(pccFileAnimal, table, 'tau')
# plot_subplots_violinplot(table, 6, 6, species, names, outfig5)


#### PCC nodes internal terminal
# table = {}
# table = nodes_pcc_analisis(pccFilePlant, table, 'plant', 'pcc')
# table = nodes_pcc_analisis(pccFileAnimal, table, 'animal', 'pcc')
# samples = list(table.keys())
# # plot_subplots_barplot(table, 1, 4, samples, samples, outfig3)
# table = {}
# table = nodes_pcc_analisis(pccFilePlant, table, 'plant', 'tau')
# table = nodes_pcc_analisis(pccFileAnimal, table, 'animal', 'tau')
# plot_subplots_violinplot(table, 1, 4, samples, samples, outfig6)

#### outgroup PCC
# pref2names = load_taxa(prefFile)
##Plants
# load_data_corr_out(outFilePlant, outtab1, 'Plants')
# add_time_dataframe(outtab1, plantTreetime, pref2names) ### change the temp file
##Animals
# load_data_corr_out(outFileAnimal, outtab2, 'Animals')
# add_time_dataframe(outtab2, animalTreetime, pref2names)
# plot_subplots_line_pcc(pathTables+'outgroup_pcc_analysis_plant_animal.tsv', pref2names, outfig7)
# plot_subplots_line_pcc2(pathTables+'outgroup_pcc_analysis_plant_animal.tsv', pref2names, outfig8)

#### outgroup TAU
# load_data_tau_out(outFilePlant, outtab3, 'Plants')
# add_time_dataframe(outtab3, plantTreetime, pref2names) ### change the temp file
# load_data_tau_out(outFileAnimal, outtab4, 'Animals')
# add_time_dataframe(outtab4, animalTreetime, pref2names) ### change the temp file
# plot_subplots_line_pcc(pathTables+'outgroup_tau_analysis_plant_animal.tsv' , pref2names, outfig9)
# plot_subplots_line_pcc2(pathTables+'outgroup_tau_analysis_plant_animal.tsv' , pref2names, outfig10)

#### Structure USagl
### All together and per major lineages
# usfiles = glob.glob(UalgFiles+'*_1.h5')
# fofiles = glob.glob(FoldsFiles+'*_1.h5')
# get_structure_data(allDupFile, fofiles, outfig11, majorlinfile)

### outgroup
# fofiles = glob.glob(FoldsFiles+'*_2.h5')
# fofiles += glob.glob(FoldsFiles+'*_3.h5')
# pref2names = load_taxa(prefFile)
# # get_structure_out(allDupFile, fofiles, pref2names, panterGenesFile, outtab5)
# # add_time_dataframe(outtab5, panterTimeFile, pref2names)
# # plot_subplots_line_pcc(outtab5, pref2names, outfig12)
# df = pd.read_csv(outtab5+'_lin', sep='\t')
# lineages = ['Viridiplantae', 'Deuterostomia', 'Protostomia', 'Fungi']
# for lin in lineages:
#     df2 = df.loc[df['lineage']==lin]
#     ax = sns.lineplot(data=df2, x="time", y="val", hue="cat", markers=True, errorbar=('ci', 95))
#     plt.show()

#### Tables
# outfile = open(outtab6, 'w')
# get_tables_exp(expDataPath+'animal_maincat_tpm_expr.h5', outfile)
# get_tables_exp(expDataPath+'plant_specific_tpm_expr.h5', outfile)
# outfile.close()