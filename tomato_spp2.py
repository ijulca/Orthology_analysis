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
from ete4 import Tree, PhyloTree
from ete4.treeview import NodeStyle, TreeStyle
from ete4.smartview import TreeLayout, TextFace
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

def filter_blast(pathBlast, protFiles):
    blastFiles = glob.glob(pathBlast+'/*.blast')
    genes = set()
    l = 200 ## 1/3 of seq lenght
    for f in blastFiles:
        for line in open(f):
            line = line.strip()
            data = line.split("\t")
            lon = float(data[3])
            e = float(data[10])
            if e <0.005 and lon >l:
                genes.add(data[1])
    print("number of genes:", len(genes))
    pepFiles = glob.glob(protFiles+'*.fa')
    outfile = open(pathBlast+"genes_all_blast.fa", "w")
    for f in pepFiles:
        seqs = GM.load_sequences(f)
        for s in seqs:
            if s in genes:
                GM.print_sequence(s.replace(':','.'),seqs[s],outfile)
    outfile.close()

def get_taxa(inFile):
    table = {}
    for line in open(inFile):
        line = line.strip()
        data = line.split("\t")
        name = data[1]
        if '_' in name:
            name = ' '.join(data[1].split('_')[:2])
        table[data[0]] = name
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
        

def tree_figure2(treeFile, outgroup, outname):
    t = PhyloTree(open(treeFile),sp_naming_function=load_species_name)
    #t.set_outgroup(t[outgroup])  ### root the tree if necessary
    events = t.get_descendant_evol_events()
    t.write(props=['evoltype'], outfile=outname+'.nw')
    t = Tree(open(outname+'.nw')) ## open as tree to modify figure
    dups = list(t.search_nodes(evoltype='D'))
    ##### To collapse the tree
    # table = {}
    # for line in open(linFile):
    #     line = line.strip()
    #     data = line.split('\t')
    #     table[data[0]] = data[1]
    # leaves = [leaf.name for leaf in t]
    # lintab = {'Solanaceae':[],'Brassicaceae':[], 'Malvaceae':[],'Fabaceae':[]}
    # for leaf in leaves:
    #     s = leaf.split('_')[-1]
    #     lin = table[s]
    #     for e in lintab:
    #         if e in lin:
    #             lintab[e].append(leaf)
    # parentcollapse = [t.common_ancestor(lintab[x]) for x in lintab]
    ## Solanaceae specific tree:
    # t = t.common_ancestor(lintab['Solanaceae'])
    ## end sol tree
    # collapsed = []
    # for n in parentcollapse:
    #     child = n.children
    #     for n in child:
    #         collapsed.append(n)
    # for node in t.traverse():
    #     if node in collapsed:
    #         parent = node.up
    #         newName = 'Collapsed-'+[leaf.name for leaf in node][0]
    #         node.detach()
    #         parent.add_child(name=newName)
    ###### till here
    ###### Graph design
    ts = TreeStyle()
    style = NodeStyle()
    vn = 2 ### width of the lines
    style['vt_line_width'] = vn
    style['hz_line_width'] = vn
    style['size'] = 0
    for n in t.traverse():
        n.set_style(style)
    for n in t.traverse():
        if n in dups:
            nstyle = NodeStyle()
            nstyle['fgcolor'] = 'red'
            nstyle['size'] = 0
            nstyle['vt_line_color'] = 'red'
            nstyle['vt_line_width'] = vn
            nstyle['hz_line_width'] = vn
            n.set_style(nstyle)
            style1 = NodeStyle()
            style1['hz_line_color'] = 'red'
            style1['vt_line_width'] = vn
            style1['hz_line_width'] = vn
            n.children[0].img_style = style1
            n.children[1].img_style = style1
        else:
            n.set_style(style)
    ts.show_leaf_name = True
    ts.show_branch_support = True ### only for the complete Tree
    ts.branch_vertical_margin = 8
    ts.scale = 100 #250 #horizontal
    ts.branch_vertical_margin = 2  #vertical
    ###### Add species name:
    table = get_taxa(taxaFile)
    for leaf in t:
        leaf.name = '_'.join(leaf.name.split('_')[:-1])+'-'+table[leaf.name.split('_')[-1]]
    t.render(outname+'.svg', w=183, units='mm', tree_style=ts) #_sol
    t.write(outfile=outname+'.nw', parser=1) #_sol, collapsed
    # t.show(tree_style=ts)

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
        
###########################################################################
### Main
path = "/home/ijulcach/projects/tomato/"

#### Blast filter
# filter_blast(path+"MADS_box_genes/split_fasta/", path+"protein_DB/")

##### Tree analysis

# treeFile = path + 'ssp_877/ssp_gene_877.alg.clean.treefile.root'
# genes = ["Solyc02g083520.2.1_HEINZ","Solyc02g061990.3.1_HEINZ", 'mRNA.Phygri02g013770.1_PHYGR',
#           'AT2G17770.2_ARATH','AT4G35900.1_ARATH', 'AMBTC02854_AMBTC',
#           'CUCSA06496_CUCSA', 'WHEAT38870_WHEAT']

# remove = ['ARATH12252_ARATH','ARATH12250_ARATH','ARATH12251_ARATH',
#           'ARATH30780_ARATH','ARATH30781_ARATH','SOLLC12257_SOLLC',
#           'SOLLC13779_SOLLC','CAPAN30218_CAPAN', 'CAPAN32113_CAPAN',
#           'SOLTU15986_SOLTU','SOLTU17722_SOLTU']
# t = PhyloTree(open(treeFile))
# ancestor = t.common_ancestor(genes)
# leaves = [leaf.name for leaf in ancestor]

# pepFiles = glob.glob(path+"protein_DB/*.fa")
# outfile = open(path +"ssp_gene2.fa", "w")
# for f in pepFiles:
#     seqs = GM.load_sequences(f)
#     seqs = {s.replace(":","."):seqs[s] for s in seqs}
#     for s in seqs:
#         k = s.replace(':','.')
#         if k in leaves and k not in remove:
#             GM.print_sequence(k,seqs[s],outfile)
# outfile.close()
# print(len(leaves)-len(remove))

#### Tree plots
#treeFile = path +"spp_gene.root.treefile"
treeFile = path +"spp_gene2.root.treefile"
outparaFile = path+"spp_gene2_paralogs.txt"
taxaFile = '/home/ijulcach/projects/tomato/taxa2species.txt'
linFile = '/home/ijulcach/projects/Land_Plants/nemo2taxa.txt'
seed = ["Solyc02g061990.3.1_HEINZ", "Solyc02g083520.2.1_HEINZ"] ## SSP2, SSP
# treeFilec = path +'ssp_179/ssp_gene_179.alg.clean.treefile'
treeFilec = path+'ssp_179/oma_gene_D0228852/oma_ssp.iqtree.treefile'
treeFilec = path +'oma_geneTree_D0228852/oma_all_ssp.alg.clean.treefile'
treeFilec = path +'oma_geneTree_D0228852/oma_all_ssp.alg.metalig.treefile'
###MADBOX:
treeFilec = '/home/ijulcach/projects/tomato/MADS_box_genes/oma_hogs/sepa_oma.phy.treefile.root.nw'
outfigtree = treeFilec+ '_fig'

#tree_figure(treeFile, path +"taxa2species.txt", path+"tree_2_3.svg")
#get_evol_events(treeFile, outparaFile)
#### Tree figure final
outgroup = 'PAPSO30197_PAPSO'  ## 'AMBTC02854_AMBTC', 'AMTR_s00152p00081700_AMBTC'
tree_figure2(treeFilec, outgroup, outfigtree)  

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