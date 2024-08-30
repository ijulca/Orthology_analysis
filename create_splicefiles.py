#!/usr/bin/env python
# coding: utf-8

## Yanis Nevers
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
######################## modify by Irene

def write_splice_file(splice_data, splice_file):
    with open(splice_file,'w') as handle_output:
        for val in splice_data.values():
            handle_output.write(";".join(val)+'\n')
            
def extract_splice_data_ensembl(fasta_file, tag):
    all_splice = {}
    with open(fasta_file, 'r') as handle, open(tag+'.fa', 'w') as outfile: 
        for seq_record in SeqIO.parse(handle, "fasta"):
            new_name = tag + seq_record.id 
            string_to_split = seq_record.description
            split = string_to_split.split('gene:') [1]
            gene_id = split.split(' ')[0]
            all_splice[gene_id] = all_splice.get(gene_id, [])+[new_name]
            seq_record.id = seq_record.description = new_name
            SeqIO.write(seq_record, outfile, "fasta")   
    return all_splice
           
def prepare_data_ensembl(fasta_file, tag):
    print(fasta_file)
    splice_value = extract_splice_data_ensembl(fasta_file, tag)
    write_splice_file(splice_value, tag+'.splice')

##########################
#### Phytozome format ####
########################## modify by Irene

def extract_splice_data_phytozome(fasta_file, tag):
    all_splice = {}
    with open(fasta_file, 'r') as handle, open(tag+'.fa', 'w') as outfile: 
        for seq_record in SeqIO.parse(handle, "fasta"):
            new_name = tag + seq_record.id 
            string_to_split = seq_record.description
            split = string_to_split.split('locus=')[1]
            gene_id = split.split(' ')[0]
            all_splice[gene_id] = all_splice.get(gene_id, [])+[new_name]
            seq_record.id = seq_record.description = new_name
            SeqIO.write(seq_record, outfile, "fasta")   
    return all_splice

def prepare_data_phytozome(fasta_file, tag):
    print(fasta_file)
    splice_value = extract_splice_data_phytozome(fasta_file, tag)
    write_splice_file(splice_value, tag+'.splice')



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


###################
#### No format ####
###################

def create_splicefile_unspecific(file, tag):
    all_splice = {}
    with open(file, 'r') as handle, open(tag+'.fa', 'w') as outfile: 
        for seq_record in SeqIO.parse(handle, "fasta"):
            new_name = tag + seq_record.id 
            gene_id = new_name
            all_splice[gene_id] = all_splice.get(gene_id, [])+[new_name]
            seq_record.id = seq_record.description = new_name
            SeqIO.write(seq_record, outfile, "fasta")
    with open(tag+'.splice','w') as handle_output:
        for val in all_splice.values():
            handle_output.write(";".join(val)+'\n')


##########
## main ##
##########
parser = argparse.ArgumentParser(description="Create the splice files for OMA")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="fasta file")
parser.add_argument("-g", "--gffFile", dest="gffFile", default='no', help="gff File, when using format gff")
parser.add_argument("-f", "--format", dest="format", default='ensembl', help="format file, default=ensembl, options = ncbi, ensembl, phyto, gff, augustus, other")
parser.add_argument("-t", "--tag", dest="tag", default='', help="tag or species name to be added to the genename")
args = parser.parse_args()

inFile = args.inFile
gffFile = args.gffFile
formato = args.format
tag = args.tag

if formato == 'ncbi':
    ## For this method, the .fa file needs to have a GeneId (GeneID=)
    print('Parsing NCBI format...')
    create_splicefile_NCBI(inFile, "fasta")
elif formato == 'gff':
    print('Parsing with gff File...')
    create_splicefile_gff(inFile, gffFile)
elif formato == 'ensembl': ## Tested
    print('Parsing ensembl format...')
    prepare_data_ensembl(inFile, tag)
elif formato =='phyto': ## Tested
    print('Parsing phytozome format...')
    prepare_data_phytozome(inFile, tag)
elif formato == 'augustus':
    create_splicefile_Augustus(inFile, "fasta")
else:
    print('Parsing no format ...') ## when you have just the protein names == genes
    create_splicefile_unspecific(inFile, tag)
