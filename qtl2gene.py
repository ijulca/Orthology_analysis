#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 18:38:16 2023

@author: ijulca
"""
import argparse
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def load_gff(gffFile):
    gff = {}
    i = 0
    for line in open(gffFile):
        line = line.strip()
        if not line or line.startswith("#"):
            pass
        else:
            data = line.split("\t")
            if data[2] == "gene":
                i += 1
                p1,p2 = int(data[3]),int(data[4])
                if p1 <p2:
                    pos = str(p1)+'-'+str(p2)
                else:
                    pos = str(p2)+'-'+str(p1)
                gene = data[8].split(';')[1].split("=")[1]
                sc = data[0]
                if sc not in gff:
                    gff[sc] = {}
                gff[sc][pos] = gene
    print("gff loaded... ",i, "genes")
    return gff

def get_genes(pos,gff):
    pos = [int(x) for x in pos.split("-")]
    positions = sorted(list(gff.keys()), key=lambda x:int(x.split("-")[0]))
    genes = set()
    if len(pos) == 2:
        q1,q2 = pos
        for num in positions:
            g1,g2 = [int(x) for x in num.split("-")]
            if g1>= q1 and g2<=q2:
                genes.add(gff[num])
            elif g1<=q1 and g2>=q1 and g2>=q2:
                genes.add(gff[num])
            elif g1<=q2 and g1>=q1 and g2>=q2:
                genes.add(gff[num])
            elif g1<=q1 and g2>=q2:
                genes.add(gff[num])
    else:
        q1 = pos[0]
        for num in positions:
            g1,g2 = [int(x) for x in num.split("-")]
            if q1>=g1 and q1<=g2:
                genes.add(gff[num])
    return genes

def analysis_qtl(gffFile, qFile):
    gff = load_gff(gffFile)
    outfile = open(qFile.split(".")[0]+"_genes.txt","w")
    for line in open(qFile):
        line =line.strip()
        data = line.split("\t")
        if data[0] == "Trait":
            line+="\tGenes"
        else:
            pos = data[3]
            genes = get_genes(pos,gff[data[1]])
            if len(genes)>=1:
                line+="\t"+";".join(genes)
            else:
                line +="\tNone"
        print(line, file=outfile)
    outfile.close()
            

# parser = argparse.ArgumentParser(description="get the genes that are in the mQTL")
# parser.add_argument("-q", "--qFile", dest="qFile", required=True, help="tab separated mqtl file")
# parser.add_argument("-g", "--genoFile", dest="genoFile", required=True, help="gff file of the genome")
# args = parser.parse_args()

# qFile = args.qFile
# gffFile = args.genoFile

# analysis_qtl(gffFile, qFile)

# print("End..")

### Plots
path = "/home/ijulca/projects/QTL/"
files = glob.glob("/home/ijulca/projects/QTL/*/*_genes.txt")

# table = {"name":[], "qtlsize":[], 'genes':[]}
# for f in files:
#     sp = f.split("/")[-2].split("_")[0]
#     for line in open(f):
#         line = line.strip()
#         data = line.split("\t")
#         if data[0] != 'Trait':
#             if "-" in data[3]:
#                 p1,p2 = [int(x) for x in data[3].split("-")]
#                 qtlsize = p2-p1
#                 gen = len(set(data[4].split(";")))
#                 table["name"].append(sp)
#                 table["qtlsize"].append(qtlsize/1000000)
#                 table['genes'].append(gen)

# df = pd.DataFrame.from_dict(table)

# ax = sns.boxplot(data=df, x="name", y="qtlsize")
# ax.set_xlabel('')
# ax.set_ylabel('QTL size (Mb)')
# ax.set_title('QTL size', fontsize=14)
# plt.savefig(path+"qtlsize.png", bbox_inches='tight')
# plt.show()

# ax = sns.boxplot(data=df, x="name", y="genes")
# ax.set_xlabel('')
# ax.set_ylabel('Nº of genes')
# ax.set_title('Number of genes in each QTL', fontsize=14)
# plt.savefig(path+"genenumber.png", bbox_inches='tight')
# plt.show()

table = {}
for f in files:
    sp = f.split("/")[-2].split("_")[0]
    if sp not in table:
        table[sp] = set()
    for line in open(f):
        line = line.strip()
        data = line.split("\t")
        if data[0] != 'Trait':
            if "," in data[0]:
                for n in data[0].split(","):
                    table[sp].add(n.strip())
            else:
                table[sp].add(data[0].strip())
names = {}
for line in open(path+"traits_ontology.txt"):
    line = line.strip()
    data = line.split("\t")
    n = data[1].split(";")[0]
    if n not in names:
        names[n] = set()
    names[n].add(data[0])

species = list(table.keys())
dades = {x:[] for x in species}
dades["Trait"] = []

for n in names:
    pref = names[n]
    dades["Trait"].append(n)
    for sp in species:
        comun = list(set(pref) & set(table[sp]))
        if len(comun)>0:
            dades[sp].append(1)
        else:
            dades[sp].append(0)
df = pd.DataFrame.from_dict(dades)
df = df.set_index('Trait')
fig, ax = plt.subplots(figsize=(5,6))  
ax = sns.heatmap(df, cmap="Blues", yticklabels=True)
ax.set_ylabel('')
ax.set_title('Traits', fontsize=14)
plt.savefig(path+"traits.png", bbox_inches='tight')
plt.show()
            
            
