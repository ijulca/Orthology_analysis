#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 10:46:24 2020

@author: ijulca
"""
import argparse 
from ete3 import NCBITaxa
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import general_modules as gmo


### load NCBITaxa ###
""" if there is an error like IntegrityError: 
    change in line 802 /home/ijulca/anaconda/lib/python3.7/site-packages/ete3/ncbi_taxonomy/ncbiquery.py:
    db.execute("INSERT INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname)) for
    db.execute("INSERT OR REPLACE INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname)) 
    rm -r /home/ijulca/.etetoolkit/ """
    
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

###

def get_taxaGroup(taxaids, group):
    ranks = ncbi.get_rank(taxaids)
    lineages = ncbi.get_taxid_translator(taxaids)
    table = set([])
    for taxaid, rank in ranks.items():
        if rank == group:
            name = lineages[taxaid]
            table.add(name+'-'+str(taxaid))
    return table

def get_taxaRanks(taxa):
    taxonRanks = {}
    descendants = ncbi.get_descendant_taxa(taxa, intermediate_nodes=True)
    clases = get_taxaGroup(descendants, 'class')
    for taxCl in clases:
        if taxCl not in taxonRanks:
            taxonRanks[taxCl] = {}
        descendants = ncbi.get_descendant_taxa(taxCl.split('-')[0], intermediate_nodes=True)
        orders = get_taxaGroup(descendants, 'order')
        for taxOr in orders:
            descendants = ncbi.get_descendant_taxa(taxOr.split('-')[0], intermediate_nodes=True)
            families = get_taxaGroup(descendants, 'family')
            taxonRanks[taxCl][taxOr] = {x:set([]) for x in families}
    return taxonRanks

def getRank_lineage(taxaid,taxonRanks):
    lineage = ncbi.get_lineage(taxaid)
    lineage_names = ncbi.get_taxid_translator(lineage)
    keys = [y+'-'+str(x) for x,y in lineage_names.items()]
    kc = [x for x in keys if x in taxonRanks][0]
    ko = [x for x in keys if x in taxonRanks[kc]][0]
    kf = [x for x in keys if x in taxonRanks[kc][ko]][0]
    return kc,ko,kf  

def count_ranks(inFile, taxonRanks):
    tax_ids = set(gmo.load_list(inFile, sep='\t', remove='taxid'))
    for taxa in tax_ids:
        kc,ko,kf = getRank_lineage(taxa,taxonRanks)
        taxonRanks[kc][ko][kf].add(taxa)
    return taxonRanks

def print_taxa_numbers(taxonRanks,outname):
    outfile = open(outname+'.rank.count', 'w')
    print('Class\tOrder\tFamily\tnum_sp\tspecies',file=outfile)
    for c in taxonRanks:
        print('Clase', c, len(taxonRanks[c]))
        for o in taxonRanks[c]:
            print('Order', o, len(taxonRanks[c][o]))
            for f in taxonRanks[c][o]:
                print('Family',f,len(taxonRanks[c][o][f]))
                string = c+'\t'+o+'\t'+f+'\t'+str(len(taxonRanks[c][o][f]))+'\t'+';'.join(taxonRanks[c][o][f])
                print(string,file=outfile)
    outfile.close()

### main
parser = argparse.ArgumentParser(description="count the number of species per taxonomica level")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="classification file")
parser.add_argument("-o", "--outFile", dest="outFile", required=True, help="outfile name")
parser.add_argument("-t", "--taxa", dest="taxa", required=True, help="higher taxonomic level. Ex: Viridiplantae")
args = parser.parse_args()

inFile = args.inFile
taxon = args.taxa
outname = args.outFile

print('We will analyse the ranks in ', taxon)

taxonRanks = get_taxaRanks(taxon)
taxonRanks = count_ranks(inFile, taxonRanks)
    
print_taxa_numbers(taxonRanks, outname)
print('End...')

