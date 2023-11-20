# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import argparse, glob
import ete3
import sys,os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import phylome_analysis as pa
import pandas as pd
import general_modules as gmo
from matplotlib import pyplot as plt

##########################
#### Getting the data ####
##########################

def parse_sp_name(node_name):  ## getting sp name
    return node_name.split("-")[-1]

def get_duplication_nodes(t):
    events = t.get_descendant_evol_events()
    dups = []
    for ev in events:
        if ev.etype == "D":
            dups.append(ev)
    return dups

def get_branchData(node, bootstrap):
    support = node.support
    data, values, species = [], [], []
    if support >=bootstrap:
        i = 0
        for n in node.iter_descendants():
            i +=1
            if i <= 2:
                branchl = n.dist
                lname = n.get_leaf_names()
                data.append(";".join(lname))
                values.append(branchl)
                species += [x.split("-")[-1] for x in lname] ## getting sp name
    return data,values, list(set(species))

def plot_hist(data, outname, maxv, bins, xlab):
    data = [maxv if x > maxv else x for x in data] ## removing high values
    plt.clf()
    plt.hist(data, bins=bins, alpha=0.5)
    plt.xlabel(xlab)
    plt.ylabel("Counts")
    plt.savefig(outname, bbox_inches='tight')

def my_layout(node): 
    if not node.is_leaf(): 
        F = ete3.TextFace(node.name, tight_text=True) 
        ete3.add_face_to_node(F, node, column=0, position="branch-top")

def my_layout_num(node): 
    if not node.is_leaf(): 
        F = ete3.TextFace(node.name.split('_')[1], tight_text=True)
        ete3.add_face_to_node(F, node, column=0, position="branch-bottom")

def get_nodeSp(sptreeFile, outName):
    spTree,spe2age = pa.load_species_tree(sptreeFile,"no")
    spTree, nodeName = pa.add_name_nodes(spTree)
    ts = ete3.TreeStyle()
    ts.layout_fn = my_layout
    spTree.show(tree_style=ts)
    spTree.render(outName+".spTree.nodes.svg", w=400, tree_style=ts)
    return spTree, spe2age

def get_ancestor_node(spTree, species):
    if len(species) == 1:
        ancestor = species[0]
    else:
        ancestor = spTree.get_common_ancestor(species)
        ancestor = ancestor.name
    return ancestor

def get_DuplicationsData(alltreeFile, sptreeFile, bootstrap, outName):
    outfile = open(outName+".duplications.txt", "w")
    print("Orthogroup\tParalog1\tParalog2\tSpecies\tNode\tBranch1\tBranch2\tABranch\tDB\tDB1A\tDB2A", file=outfile)
    genetrees = pa.get_trees_from_file(alltreeFile)
    print('Number of tree to be analysed:', len(genetrees))
    spTree,spe2age = get_nodeSp(sptreeFile, outName)
    distances, ndups = [], []
    for group in genetrees:
        tree = genetrees[group]
        t = pa.load_tree(tree,spe2age)
        t = pa.delete_empty_leaves(t)
        distances += pa.get_average_branchLen(tree,bootstrap) ## getting all the distances of the tree
        t.set_species_naming_function(parse_sp_name)
        dups = get_duplication_nodes(t)
        ndups.append(len(dups))
        for ev in dups:
            node = ev.node
            leaves, branch, species = get_branchData(node, bootstrap)
            if len(leaves) > 0:
                ancestor = get_ancestor_node(spTree, species)
                string = group+"\t"+"\t".join(leaves)+"\t"+";".join(species)+"\t"+ancestor
                string += "\t"+"\t".join([str(x) for x in branch])+"\t"+str(node.dist)
                string += "\t"+str(abs(branch[0]-branch[1]))+"\t"+str(branch[0]-node.dist)+"\t"+str(branch[1]-node.dist)
                print(string, file=outfile)
    outfile.close()
    plot_hist(distances, outName+"_total_branch_length.png",1,100, "Branch length")
    plot_hist(ndups, outName+"_number_dups.png", 20, 10, "Number of duplications")


    
##################
#### Analysis ####        
##################

def tree_nodesize(sptreeFile, outName, values):
    t,spe2age = pa.load_species_tree(sptreeFile,"no")
    t, nodeName = pa.add_name_nodes(t)
    nstyle = ete3.NodeStyle()
    nstyle['size']= 0
    for n in t.traverse():
        n.set_style(nstyle)
    ts = ete3.TreeStyle()    
    for n in t.traverse():
        nstyle = ete3.NodeStyle()
        nstyle["fgcolor"] = "red"
        nstyle["size"] = values[n.name]
        n.set_style(nstyle)
    t.img_style["size"] = 20
    t.img_style["fgcolor"] = "blue"
    ts = ete3.TreeStyle()
    ts.layout_fn = my_layout_num
    t.render(outName+".parasize.svg", w=400, tree_style=ts)

def getting_DF_paralogs(outName, sptreeFile, fold):
    dists = []
    nodes = {}
    tot,simpara, simparasp, simpara1, difpara = 0,0,0,0,0
    evolgenes = []
    dif = 0
    for line in open(outName+".duplications.txt"):
        line = line.strip()
        data = line.split("\t")
        if data[0] != "Orthogroup":
            dists += [float(x) for x in data[-2:]]
            n = data[4]
            if n not in nodes:
                nodes[n] = 0
            nodes[n]+=1
            tot += 1
            para1, para2 = data[1].split(";"), data[2].split(";")
            sp1,sp2 = [x.split("-")[-1] for x in para1], [x.split("-")[-1] for x in para2]
            br1,br2 = [float(x) for x in data[-2:]]
            if br1 >br2:
                rat = br1/br2
                fast,slow = 1,2
            else:
                rat = br2/br1
                fast,slow = 2,1
            if rat >= fold:
                dif += 1
                if len(para1) == len(para2):
                    simpara +=1
                    if sp1 == sp2:
                        simparasp+=1
                    if len(para1) == 1:
                        simpara1 +=1
                        evolgenes.append(data[fast]+"\t"+data[slow])
    print("Total number of duplication events:",tot)
    print("Number of duplications with different "+str(fold)+" fold:", dif*100/tot)
    print("Number of paralagos with the same number of leaves",simpara*100/dif)
    print("Number of paralagos with the same number of leaves and species",simparasp*100/dif)
    print("Number of paralagos size 2",simpara1*100/dif)
    nodes = {x:nodes[x]*100/tot for x in nodes} ## percentage of duplications per node
    plot_hist(dists, outName+"_para_branch_length.png",1,100, "Branch length")
    tree_nodesize(sptreeFile, outName, nodes)
    evolgenes = sorted(evolgenes, key=lambda x:x.split("-")[-1])
    outfile = open(outName+".fast.slow.txt", "w")
    print("\n".join(evolgenes),file=outfile)
    outfile.close()

def getting_DF_paralogs2(outName, sptreeFile, fold):
    dists, dife = [], []
    nodes = {}
    tot,simpara, simparasp, simpara1, difpara = 0,0,0,0,0
    evolgenes = []
    dif = 0
    for line in open(outName+".duplications.txt"):
        line = line.strip()
        data = line.split("\t")
        if data[0] != "Orthogroup":
            dists += [float(x) for x in data[5:7]]
            values = [float(x) for x in data[8:]]
            dife.append(values[2])
    #         n = data[4]
    #         if n not in nodes:
    #             nodes[n] = 0
    #         nodes[n]+=1
    #         tot += 1
    #         para1, para2 = data[1].split(";"), data[2].split(";")
    #         sp1,sp2 = [x.split("-")[-1] for x in para1], [x.split("-")[-1] for x in para2]
    #         br1,br2 = [float(x) for x in data[-2:]]
    #         if br1 >br2:
    #             rat = br1/br2
    #             fast,slow = 1,2
    #         else:
    #             rat = br2/br1
    #             fast,slow = 2,1
    #         if rat >= fold:
    #             dif += 1
    #             if len(para1) == len(para2):
    #                 simpara +=1
    #                 if sp1 == sp2:
    #                     simparasp+=1
    #                 if len(para1) == 1:
    #                     simpara1 +=1
    #                     evolgenes.append(data[fast]+"\t"+data[slow])
    # print("Total number of duplication events:",tot)
    # print("Number of duplications with different "+str(fold)+" fold:", dif*100/tot)
    # print("Number of paralagos with the same number of leaves",simpara*100/dif)
    # print("Number of paralagos with the same number of leaves and species",simparasp*100/dif)
    # print("Number of paralagos size 2",simpara1*100/dif)
    # nodes = {x:nodes[x]*100/tot for x in nodes} ## percentage of duplications per node
    plot_hist(dists, outName+"_para_branch_length.png",1,100, "Branch length")
    plot_hist(dife, outName+"_para_Dif_branch_length.png",1,100, "Difference in Branch length")
    dife_sort = sorted(dife, reverse=False) 
    x = int(float(85*len(dife))/float(100)) ### choose the filter for distance neofunc
    print(x) 
    n = dife_sort[x] 
    print(n)

    # tree_nodesize(sptreeFile, outName, nodes)
    # evolgenes = sorted(evolgenes, key=lambda x:x.split("-")[-1])
    # outfile = open(outName+".fast.slow.txt", "w")
    # print("\n".join(evolgenes),file=outfile)
    # outfile.close()
            
def get_organ_para(outName,organFile):
    gene2organ = gmo.load_dic(organFile)
    outfile = open(outName+".fast.slow.organ.txt","w")
    for line in open(outName+".fast.slow.txt"):
        line = line.strip()
        data = line.split("\t")
        names = [x.split("-")[0] for x in data]
        string = line
        for n in names:
            if n in gene2organ:
                string += "\t"+gene2organ[n]
            else:
                string += "\tNONE"
        print(string,file=outfile)
    outfile.close()
    
def get_table(inFile):
    table = {}
    for line in open(inFile):
        line = line.strip()
        data = line.split("\t")
        if "Mnemonic" not in line:
            table[data[0]] = data[2:]
    return table

#### Main
path = "/home/ijulca/projects/ldo_project/"
alltreeFile = '/home/ijulca/projects/ldo_project/alltree.txt'
sptreeFile = "/home/ijulca/projects/ldo_project/Trees_sp/SpeciesTree_all.txt"
outName = "/home/ijulca/projects/ldo_project/plants"
organFile = "/home/ijulca/projects/ldo_project/genes2organs.txt"
taxaFile = "/home/ijulca/projects/ldo_project/table_animals.dupev.tsv"
speciesFile = "/home/ijulca/projects/ldo_project/Alex_data/data/panther-17.0/species_tree.nhx"
animalTreeFile = path+"animal_tree.txt"
picPath = "/home/ijulca/Programs/phylopic-tools/images/"
spTreeanimal = path+"Figures/species_tree_animals.svg"
bootstrap = 60
fold = 1.5


### Get the data
#get_DuplicationsData(alltreeFile, sptreeFile, bootstrap, outName)


### Analyse the data
#getting_DF_paralogs(outName, sptreeFile, fold)
##getting_DF_paralogs2(outName, sptreeFile, fold)
#get_organ_para(outName,organFile)


### Alex data
### Get the image of the species tree
# pics= glob.glob(picPath+"*")
# table = get_table(taxaFile)
# t1 = ete3.Tree(speciesFile, format=1)
# t = t1.get_common_ancestor(list(table.keys()))
# t.prune(list(table.keys()))
# t.write(format=2, outfile=animalTreeFile)

# for node in t.traverse():
#     node.img_style['size'] = 0
#     if node.is_leaf():
#         sp = table[node.name][0]
#         name_face = ete3.TextFace(sp, fstyle="italic", fsize=30)
#         node.add_face(name_face, column=0, position='aligned')
#         pic = ete3.faces.ImgFace(picPath+node.name+".svg", width=60, height=60)
#         node.add_face(pic, column=1, position='aligned')
        
# ts = ete3.TreeStyle()
# ts.draw_guiding_lines = True
# ts.scale = 20
# ts.branch_vertical_margin = 10 
# ts.show_leaf_name = False
# t.render(spTreeanimal, dpi=300, w=1000, tree_style=ts)
#t.show(tree_style=ts)

plt.figure(figsize=(15,4))
df = pd.read_csv(taxaFile, sep="\t", header=0, index_col=None)
ax = df.plot.barh(x="Mnemonic",y=["No Outliers", "One Outlier (MDO)", "Both Outliers"], stacked=True)

# ax = df.plot.bar(stacked=True)
plt.show()


    
