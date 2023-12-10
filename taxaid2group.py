#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 09:51:58 2020

@author: ijulca
"""
import argparse 
from Bio import Entrez
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import general_modules as gmo
from ete4 import NCBITaxa
ncbi = NCBITaxa()

### use this link of ncbi to get the lineages:
### https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi

Entrez.email = 'irene.julcac@gmail.com'

def get_lineage(inFile):
    outfile = open(inFile+'.lineage_class', 'w')
    tax_ids = set(gmo.load_list(inFile, sep='\t'))
    for tax in tax_ids:
        handle = Entrez.efetch('taxonomy', id=tax, rettype='xml')
        response = Entrez.read(handle)

        for entry in response:
            sci_name = entry.get('ScientificName')
            lineage_taxa = entry.get('Lineage').split(';')
            lineage_taxa = [x.strip() for x in lineage_taxa]
            print(tax, sci_name, '\t'.join(lineage_taxa), sep='\t',file=outfile)
    outfile.close()

def get_all_taxaid(inFile):
    taxaid = set()
    for line in open(inFile):
        line = line.strip()
        if not line or "taxname" in line:
            pass
        else:
            data = [x.strip() for x in line.split("|")]
            if len(data) >= 4:
                ids = data[3].split(" ")
                for i in ids:
                    taxaid.add(i)
    outfile = open("list_taxaid1.txt", "w")
    print("\n".join(taxaid), file=outfile)
    outfile.close()

def get_list(inFile):
    ids = set()
    for line in open(inFile):
        line = line.strip()
        ids.add(line)
    print("number of ids:",len(ids))
    return ids

def get_ids2lineages(inFile):
    ids2lineages, ids2name = {}, {}
    for line in open(inFile):
        line = line.strip()
        if not line or "taxname" in line:
            pass
        else:
            data = [x.strip() for x in line.split("|")]
            if len(data) >= 4:
                taxon = data[3].split(" ")
                taxon.reverse()
                ids2lineages[data[1]] = taxon
                ids2name[data[1]] = data[2]
    return ids2lineages, ids2name

def get_lineages(ids, ids2lineages, ids2name):
    outfile = open("lineages.txt","w")
    for i in ids:
        if i in ids2name:
            n = ids2name[i]
            taxalin = ids2lineages[i]
            taxanam = [ids2name[x] for x in taxalin]
            string = i+"\t"+n+"\t"+";".join(taxalin)+"\t"+";".join(taxanam)
            print(string, file=outfile)
    outfile.close()

def get_ncbi_genomes(inFile):
    outfile = open("plant_genomes_ncbi.txt", "w")
    for line in open(inFile):
        line = line.strip()
        data = line.split("\t")
        if data[0] != 'Assembly Accession':
            name = data[2]
            size = int(data[10])
            name2taxaid = ncbi.get_name_translator([name])
            taxa = name2taxaid[name][0]
            lineage = ncbi.get_lineage(taxa)
            if 33090 in lineage: ### Plant lineage
                string = name+"\t"+str(taxa)+"\t"+";".join([str(x) for x in lineage])
                print(string, file=outfile)
    outfile.close()
 
def get_ensembl_genomes(inFile):
    outfile = open("plant_genomes_ensembl.txt", "w")
    for line in open(inFile):
        line = line.strip()
        line = line.replace('"','')
        data = [x for x in line.split(",") if x !='']
        if data[0] != 'Name':
            name = data[0]
            taxa = int(data[2])
            lineage = ncbi.get_lineage(taxa)
            if 33090 in lineage: ### Plant lineage
                string = name+"\t"+str(taxa)+"\t"+";".join([str(x) for x in lineage])
                print(string, file=outfile)
    outfile.close()               

def get_phytozome_genomes(inFile):
    outfile = open("plant_genomes_phyto.txt", "w")
    for line in open(inFile):
        name = line.strip() 
        #print(name)
        name2taxaid = ncbi.get_name_translator([name])
        taxa = name2taxaid[name][0]
        lineage = ncbi.get_lineage(taxa)
        if 33090 in lineage: ### Plant lineage
            string = name+"\t"+str(taxa)+"\t"+";".join([str(x) for x in lineage])
            print(string, file=outfile)
    outfile.close()

def analysis2taxa(inFile, linFile):
    lineages = {}
    for line in open(linFile):
        line = line.strip()
        data = line.split("\t")
        if data[0] not in lineages:
            lineages[data[0]] = set()
        lineages[data[0]].add(data[2])
    species = {}
    for line in open(inFile):
        line = line.strip()
        data = line.split("\t")
        sn = data[0].split(" ")
        sp = data[0]
        if len(sn) >2:
            if "x" not in sn:
                sp = " ".join(sn[:2])
        if sp not in species:
            species[sp] = data[2]
    print("number of species in total:",len(species))
    counts = {x:0 for x in lineages}
    for s in species:
        taxids = species[s].split(";")
        key = set()
        for taxon in lineages:
            linid = lineages[taxon]        
            comon = list(set(linid) & set(taxids))
            if len(comon) != 0:
                key.add(taxon)
        if len(key) == 0:
            print(s)
            
def taxaid2lineage(inFile): ### using ete4
    outfile = open(inFile+'.lineage.txt','w')
    for line in open(inFile):
        taxid = line.strip()
        lineage = ncbi.get_lineage(taxid) 
        names = ncbi.get_taxid_translator(lineage)
        lineage_names = [names[x] for x in lineage]
        lineage2ranks = ncbi.get_rank(names)
        lineage_ranks = [lineage2ranks[x] for x in lineage]
        string = taxid+'\t'+'; '.join(lineage_names)+'\t'+'; '.join(lineage_ranks)
        print(string, file=outfile)
    outfile.close()
      
            
### main
parser = argparse.ArgumentParser(description="get the taxonomic classification of samples")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="result of ncbi: tax_report file")
parser.add_argument("-n", "--taxlist", dest="taxlist", help="taxaid list. Give with tag 'a'")
parser.add_argument("-t", "--tag", dest="tag", required=True, help="task to do: l:list of taxaid, \
                    a: list of taxaid and lineages for the taxid list,\
                    n: ncbi dataset, p: phytozome list of species,\
                    b: analysis, t list of taxaids")
args = parser.parse_args()

inFile = args.inFile
idsFile = args.taxlist
tag = args.tag

if tag == "l":
    print("getting the taxaID list...")
    get_all_taxaid(inFile)
elif tag =='t':
    print('getting the taxaID list from list...')
    taxaid2lineage(inFile)
elif tag == "a":
    print("getting the lineages...")
    ids = get_list(idsFile)
    ids2lineages, ids2name = get_ids2lineages(inFile)
    get_lineages(ids, ids2lineages, ids2name)
elif tag == 'n':
    print("getting ncbi genomes...")
    get_ncbi_genomes(inFile)
elif tag == 'e':
    print("getting ensembl genomes...")
    get_ensembl_genomes(inFile)
elif tag == "p":
    print("getting the lineages from list of names")
    get_phytozome_genomes(inFile)    
elif tag == "b":
    print("analysis started...")
    analysis2taxa(inFile, idsFile)
    
    
print('End...')