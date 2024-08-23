#!/usr/bin/env python
# coding: utf-8

## Silvia Prieto
## Irene Julca - 23.08.24

import argparse
from Bio import SeqIO
import numpy as np
import re
import numpy as np
import pandas as pd
import os
import gffutils

#####################
#### NCBI format ####
#####################
## For this method, the .fa file needs to have a GeneId (GeneID=)
def create_splicefile_NCBI(file, formato):
    '''
    Require num. Create splicefile from multifasta file, specific for NCBI fastas.
    '''
    geneids = []
    for seq_record in SeqIO.parse(file, formato):
        string_to_split = seq_record.description
        split = string_to_split.split('GeneID=') [1]
        split = split.split(']')[0]
        split = split.replace(']','')
        geneids.append(split)


    splicedic = {}
    counter = -1

    for index, record in enumerate(SeqIO.parse(file, formato)):
        counter+= 1
        key = geneids[counter]
        splicedic.setdefault(key, [])
        splicedic[key].append(record.id)

    
    temp = file.split('.fa')[0] +'.temp'
    with open(temp, 'w') as outfile:
        for key in splicedic:
        
            for i in splicedic[key]:
                print(i, end =";", file=outfile)
            print('\n', file = outfile)
            
    splicefile = temp.split('.temp')[0] + '.splice'
    with open(temp,'r') as infile:
        lines = infile.readlines()
    
        with open(splicefile, 'w') as outfile:
            for l in lines:
                line = re.sub(";$","", l)
                if len(line.strip('\n'))>1 and (len(re.findall(';', line)))>0:
                    outfile.write(line)
                    
    return 'Done. Check your folder'


########################
#### Ensembl format ####
########################

def create_splicefile_OMA(file, formato):
    '''
    Create splicefile from multifasta file, specific for Ensembl fastas with pipes separated header ids (OMA).
    '''
    #Make array with gene IDs
    geneids = []
    for seq_record in SeqIO.parse(file, formato):
        string_to_split = seq_record.description
        split = string_to_split.split(' | ') [3]
        split = split.split(' | ')[0]
        geneids.append(split)


    # Make dictionary with gene ids from previous array as keys. Iterate through the sequences in fasta and add
    # individual unique record ids as values to each gene id's key.
    splicedic = {}
    counter = -1

    for index, record in enumerate(SeqIO.parse(file, formato)):
        counter+= 1
        key = geneids[counter]
        splicedic.setdefault(key, [])
        splicedic[key].append(record.id)

    
    temp = file.split('.fa')[0] +'.temp'
    with open(temp, 'w') as outfile:
        for key in splicedic:
        
            for i in splicedic[key]:
                print(i, end =";", file=outfile)
            print('\n', file = outfile)
            
    splicefile = temp.split('.temp')[0] + '.splice'
    with open(temp,'r') as infile:
        lines = infile.readlines()
    
        with open(splicefile, 'w') as outfile:
            for l in lines:
                line = re.sub(";$","", l)
                if len(line.strip('\n'))>1 and (len(re.findall(';', line)))>0:
                    outfile.write(line)
                    
    return 'Done. Check your folder'


def create_splicefile_Ensembl(file, formato):
    '''
    Require num and biopython. 
    Create splicefile from multifasta file, specific for Ensembl fastas with gene:id and no pipes.
    formato is fasta (should be)
    '''

    splicedic={}

    for seq_record in SeqIO.parse(file, formato):
        string_to_split = seq_record.description
        gid = string_to_split.split('gene:')[1]
        gid = gid.split(' transcript')[0]
        splicedic.setdefault(gid,[])
        splicedic[gid].append(seq_record.id)

    
    temp = file.split('.fa')[0] +'.splice'
    with open(temp, 'w') as outfile:
        for key in splicedic:
            print(';'.join(splicedic[key]), file=outfile)
            
                    
    return 'Done. Check your folder'

#########################
#### GFF information ####
#########################
## Method 2, no GeneID --> done with gff files instead

def load_gff(gff):
    gene_to_prot= {}
    prot_to_gene = {}
    db = gffutils.create_db(gff, ':memory:', merge_strategy="create_unique", keep_order=True)
    # Loop through all genes
    for t in db.features_of_type('gene', order_by='start'):
        gene = t.id
        gene_list = []
        ordered_child = list(db.children(t, featuretype='CDS', order_by='start'))
        
        # Loop through all children of genes
        for child in ordered_child:
            type_attribute = ['protein_id', 'Name']
            
            # Loop through all proteins of children??
            for att_type in type_attribute:
                protein = child.attributes.get(att_type, [None])[0]
                if protein:
            
                    break
            if not protein:
                print('warning')
                print(child)
                continue
            corr_gene = prot_to_gene.get(protein, None)
            if corr_gene and corr_gene!=gene:
                gene = corr_gene
                for other_prot in gene_list:
                    prot_to_gene[other_prot] = gene
                gene_list = gene_to_prot[gene]+gene_list
            else:
                prot_to_gene[protein] = gene
            if protein not in gene_list:
                gene_list.append(protein)
        if len(gene_list)!=0:
            gene_to_prot[gene] = gene_list
    return gene_to_prot
        
def create_splicefile_gff(file, gff):
    gene_to_prot = load_gff(gff)
    temp = file.split('.fa')[0] +'.splice' 
    with open(temp,'w') as handle_output:
        for val in gene_to_prot.values():
            handle_output.write(";".join(val)+'\n')

###############################################
#### Augustus with alternative transcripts ####
###############################################
def create_splicefile_Augustus(file, formato):
    '''
    Require num and biopython. 
    Create splicefile from multifasta file, specific for Augustus fastas with gene= and sequence
    ids of the form ID.t1/ ID.t2/ ID.t3...
    Formato should be fasta
    '''

    from Bio import SeqIO
    import numpy as np
    import re

    splicedic={}

    for seq_record in SeqIO.parse(file, 'fasta'):
        string_to_split = seq_record.description
        gene = string_to_split.split('gene=')[1]
        splicedic.setdefault(gene,[])
        splicedic[gene].append(seq_record.id)

    temp = file.split('.fa')[0] +'.splice'
    with open(temp, 'w') as outfile:
        for key in splicedic:
            print(';'.join(splicedic[key]), file=outfile)
            
    return 'Done. Check your folder'




def create_splicefile_unspecific(file, formato):
    '''
    Require num and biopython. 
    Create splicefile from multifasta file, for fastas with isoform ids 
    of the form geneid.1/ geneid.2/ geneid.3
    Formato should be fasta
    '''

    from Bio import SeqIO
    import numpy as np
    import re

    splicedic={}

    for seq_record in SeqIO.parse(file, 'fasta'):
        string_to_split = seq_record.id
        gene = string_to_split.split('.')[0]
        splicedic.setdefault(gene,[])
        splicedic[gene].append(seq_record.id)
        
    temp = file.split('.fa')[0] +'.splice'
    with open(temp, 'w') as outfile:
        for key in splicedic:
            print(';'.join(splicedic[key]), file=outfile)
            
    return 'Done. Check your folder'



##########
## main ##
##########
parser = argparse.ArgumentParser(description="Create the splice files for OMA")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="fasta file")
parser.add_argument("-t", "--tag", dest="tag", default='ensembl', help="format file, default=ensembl, options=ncbi, ensembl, gff, augustus, other")
parser.add_argument("-g", "--gffFile", dest="gffFile", default='no', help="gff File, when using format gff")
args = parser.parse_args()

inFile = args.inFile
tag = args.tag
gffFile = args.gffFile

if tag == 'ncbi':
    ## For this method, the .fa file needs to have a GeneId (GeneID=)
    print('Parsing NCBI format...')
    create_splicefile_NCBI(inFile, "fasta")
elif tag == 'gff':
    print('Parsing with gff File...')
    create_splicefile_gff(inFile, gffFile)
elif tag == 'ensembl':
    create_splicefile_Ensembl(inFile, "fasta")
elif tag == 'augustus':
    create_splicefile_Augustus(inFile, "fasta")
else:
    create_splicefile_unspecific(inFile, 'fasta')
