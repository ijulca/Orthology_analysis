#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 12:13:57 2023

@author: ijulca
"""
import glob
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM
import phylome_analysis as pa
import ete4
import random
from collections import Counter
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

def parse_sp_name(node_name):  ## getting sp name
    return node_name.split("_")[-1]

def load_species_name(node):
	return node.split("_")[-1] ### change

def get_sp2age(treeFile, name):
    t = ete3.PhyloTree(treeFile,sp_naming_function=load_species_name)
    leaves = list(leaf.name for leaf in t.get_leaves())
    leaves = list(set([x.split("_")[-1] for x in leaves if x.split("_")[-1] != name]))+[name]
    spe2age = {x:i for i,x in enumerate(leaves, start=1)}
    return spe2age

def filter_blast(blastFiles, protFiles):
    blastFiles = glob.glob(blastFiles)
    genes = set()
    for f in blastFiles:
        for line in open(f):
            line = line.strip()
            data = line.split("\t")
            e = float(data[10])
            if e <0.005:
                genes.add(data[1])
    print("number of genes:", len(genes))
    pepFiles = glob.glob(protFiles)
    outfile = open(path +"ssp_gene_all.fa", "w")
    for f in pepFiles:
        seqs = GM.load_sequences(f)
        for s in seqs:
            if s in genes:
                GM.print_sequence(s,seqs[s],outfile)
    outfile.close()

def get_taxa(inFile):
    table = {}
    for line in open(inFile):
        line = line.strip()
        data = line.split("\t")
        table[data[0]] = data[2]
    return table
   
def collapse_leaves(t,table, species):
    a = len(list(t.iter_leaves()))
    seed = ["Solyc02g061990.3.1_HEINZ", "Solyc02g083520.2.1_HEINZ"]
    for node in t.traverse():
        if node.is_leaf()==False:
            # list names of leaves
            leaf_names=[leaf.name for leaf in node.iter_leaves()]
            names = set([table[x.split("_")[-1]] for x in leaf_names])
            spn = set([x.split(" ")[0] for x in names]) ## get solanum
            # spn = set()
            # for s in names:
            #     if s in species:
            #         spn.add(s)
            if len(spn) == 1 and "S." in spn: #len(names):
                comon = [x for x in leaf_names if x in seed]
                if len(comon) >= 1:
                    node.add_feature("collapse", comon[0]+"-collapsed")
                else:
                    node.add_feature("collapse",random.choice(leaf_names)+"-collapsed")
             	
    for node in t.traverse():
        if "collapse" in node.features:
            if node.is_root():
                t = None
            else:
                parent = node.up
                newName = node.collapse
                node.detach()
                parent.add_child(name=newName)
    print("Number of leaves, start, end",a,len(list(t.iter_leaves())))
    return t

def tree_figure(treeFile, prefixFile, outname):
    table = get_taxa(prefixFile)
    t = ete3.Tree(treeFile)
    leaf_color = {"S. lycopersicum":"red", "S. pimpinellifolium":"blue"}
    t = collapse_leaves(t,table, list(leaf_color.keys()))  ### collapsa branches
    ts = ete3.TreeStyle()
    for node in t.traverse():
        node.img_style['size'] = 0
        if node.is_leaf():
            sp = node.name.split("_")[-1]
            if "collapsed" in sp:
                sp = sp.split("-")[0]
            if table[sp] in leaf_color:
                color = leaf_color[table[sp]]
            else:
                color = "black"
            name_face = ete3.TextFace(node.name, fgcolor=color, fsize=100)
            node.add_face(name_face, column=0, position='branch-right')
        else:
            nstyle = ete3.NodeStyle()
            nstyle["fgcolor"] = "darkred"
            nstyle["size"] = node.support
            node.set_style(nstyle)
    # Set up style for circular tree
    ts.scale = 1000
    # Disable the default tip names config
    ts.show_leaf_name = False
    # Draw Tree
    t.render(outname, dpi=300, w=1000, tree_style=ts)

def get_duplication_nodes(t):
    events = t.get_descendant_evol_events()
    dups = []
    for ev in events:
        if ev.etype == "D":
            dups.append(ev)
    return dups

def get_branchData(node):
    data = []
    i = 0
    for n in node.iter_descendants():
        i +=1
        if i <= 2:
            lname = list(n.get_leaf_names())
            data.append(";".join(lname))
    return data

def get_evol_events(treeFile, outFile):
    outfile = open(outFile, 'w')
    t = ete3.PhyloTree(treeFile)
    t.set_species_naming_function(parse_sp_name)
    dups = get_duplication_nodes(t)
    j = 0
    for ev in dups:
        j+=1
        paralogs = get_branchData(ev.node)
        string = str(j)+"\t"+"\t".join(paralogs)
        print(string, file=outfile)
    outfile.close()

def sra2names(inFile):
    table = {}
    for line in open(inFile):
        line = line.strip()
        data = line.split("\t")
        table[data[1]] = data[3]
    return table
        

### Main
path = "/home/ijulcach/projects/tomato_project/SSP2_analysis/"

#### Blast filter
filter_blast(path+"blast_results/*.blast", path+"protein_DB/*.fa")

##### Tree analysis

# treeFile = path + 'spp_gene.alg.clean.treefile'
# rootTreeFile = path+'spp_gene.root.treefile'
# out = 'ARATH'

# spe2age = get_sp2age(treeFile, out)
# t = pa.load_tree(treeFile,spe2age)
# #.write(format=2, outfile=rootTreeFile)
# ancestor = t.get_common_ancestor(['Solyc02g061990.3.1_HEINZ',"Solyc02g083520.2.1_HEINZ", "AT2G17770.2_ARATH"])
# leaves = list(leaf.name for leaf in ancestor.get_leaves())
# pepFiles = glob.glob(path+"protein_DB/*.fa")
# outfile = open(path +"spp_gene2.fa", "w")
# for f in pepFiles:
#     seqs = GM.load_sequences(f)
#     seqs = {s.replace(":","."):seqs[s] for s in seqs}
#     for s in seqs:
#         if s in leaves:
#             GM.print_sequence(s.replace(":","."),seqs[s],outfile)
# outfile.close()
# print(len(leaves))

#### Tree plots
#treeFile = path +"spp_gene.root.treefile"
treeFile = path +"spp_gene2.root.treefile"
outparaFile = path+"spp_gene2_paralogs.txt"
seed = ["Solyc02g061990.3.1_HEINZ", "Solyc02g083520.2.1_HEINZ"] ## SSP2, SSP

#tree_figure(treeFile, path +"taxa2species.txt", path+"tree_2_3.svg")
#get_evol_events(treeFile, outparaFile)

# table = get_taxa(path +"taxa2species.txt")
# names = [x for x in table if table[x].split(" ")[0]=="S."]
# print(len(names))
# for line in open(outparaFile):
#     line = line.strip()
#     data = line.split("\t")
#     if data[0] == "2":
#         para1, para2 = data[1].split(";"), data[2].split(";") ## SSP, SSP2
#         if seed[0] in para2:
#             print("Yessssssssssssssss")
#         para1 = [x for x in para1 if x.split("_")[-1] in names]
#         para2 = [x for x in para2 if x.split("_")[-1] in names]
#         print("ssp:",len(para1),"ssp2", len(para2))
#         sp1,sp2 = [x.split("_")[-1] for x in para1], [x.split("_")[-1] for x in para2]
#         print("ssp:",len(set(sp1)), "ssp2:",len(set(sp2)))
#         loss2 = [x for x in sp1 if x not in sp2]
#         loss1 = [x for x in sp2 if x not in sp1]
#         print("losses of ssp1",loss1)
#         print("losses of ssp2",loss2)
#         counts = Counter(sp2)
#         for e in counts:
#             if counts[e] >1:
#                 print("duplications",e, counts[e])

###### Expression analysis
# expFile = path +'ssp_expression.evorepro.csv'
# expFile = path+"ssp_ARATH.expression.txt"
# table = sra2names(path+"Evorepro.all_sampleAnnotation.csv")
# organs = ['Apical meristem', 'Root meristem','Flower', 'Stem', 'Female', 'Root', 'Male', 'Leaf', 'Seeds']
# for line in open(expFile):
#     line = line.strip()
#     data = line.split("\t")
#     if data[0] == "ids":
#         names = [x+"--"+table[x] for x in data[1:]]
#         new_names = {x:table[x] for x in data[1:]}
#         dropcol = [x for x in data[1:] if table[x] not in organs[:2]]

# sort_map = {day: next((idx for idx, val in enumerate(organs) if val in day),
#                  len(names)) for day in names}
# cols = sorted(names, key=sort_map.__getitem__)
# cols = [x.split("--")[0] for x in cols]

# df = pd.read_csv(expFile, sep="\t", header=0, index_col=0)
# df = df[cols]
# df = df.drop(dropcol, axis=1)
# # df = df.rename(columns=new_names)
# for index, row in df.iterrows():
#     # df.loc[index] = (row-row.mean())/row.std()
#     df.loc[index] = row/row.max()
# plt.figure(figsize=(15,8))
# sns.set(font_scale=1.5)
# ax = sns.heatmap(df, cmap="Blues", xticklabels=True)
# plt.savefig(path+"expression_meri_ARATH.svg",bbox_inches = "tight")
# plt.show()
# #corr = df.loc["Solyc02g061990.3.1"].corr(df.loc["Solyc02g083520.2.1"])
# corr = df.loc["AT2G17770"].corr(df.loc["AT4G35900"])
# print(corr)