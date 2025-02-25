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


def write_splice_file(splice_data, splice_file):
    with open(splice_file,'w') as handle_output:
        for g,val in splice_data.items():
            handle_output.write(g+'\t'+";".join(val)+'\n')

#####################
#### NCBI format ####
#####################
## For this method, the .fa file needs to have a GeneId (GeneID=)

def load_gff_ncbi(inFile):
    rna2gene = {}
    pep2rna = {}
    for line in open(inFile):
        line = line.strip()
        if line.startswith('#') or not line:
            pass
        else:
            data = line.split('\t')
            if data[2] == 'mRNA':
                tr = data[8].split('ID=')[1].split(';')[0]
                g = data[8].split('Parent=')[1].split(';')[0]
                rna2gene[tr] = g
            elif data[2] == 'CDS':
                p = data[8].split('ID=')[1].split(';')[0]
                if 'cds' in p:
                    p = p.split('cds-')[1]
                tr = data[8].split('Parent=')[1].split(';')[0]
                pep2rna[p] = tr
    pep2gene = {x:rna2gene[pep2rna[x]] for x in pep2rna}
    print(f"Proteins loaded:{len(pep2gene)}")
    return pep2gene

def extract_splice_data_ncbi(fasta_file, gff_file, tag):
    pep2gene = load_gff_ncbi(gff_file)
    all_splice = {}
    with open(fasta_file, 'r') as handle, open(tag+'.fa', 'w') as outfile: 
        for seq_record in SeqIO.parse(handle, "fasta"):
            new_name = tag + seq_record.id 
            gene_id = pep2gene[seq_record.id]
            all_splice[gene_id] = all_splice.get(gene_id, [])+[new_name]
            seq_record.id = seq_record.description = new_name
            SeqIO.write(seq_record, outfile, "fasta")   
    return all_splice
           
def create_splicefile_ncbi(fasta_file, gff_file, tag):
    print(fasta_file)
    splice_value = extract_splice_data_ncbi(fasta_file, gff_file, tag)
    write_splice_file(splice_value, tag+'.gene.splice')


########################
#### Ensembl format ####
########################
            
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
           
def create_splicefile_ensembl(fasta_file, tag):
    print(fasta_file)
    splice_value = extract_splice_data_ensembl(fasta_file, tag)
    write_splice_file(splice_value, tag+'.gene.splice')

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

def create_splicefile_phytozome(fasta_file, tag):
    print(fasta_file)
    splice_value = extract_splice_data_phytozome(fasta_file, tag)
    write_splice_file(splice_value, tag+'.gene.splice')


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

    temp = file.split('.fa')[0] +'.gene.splice'
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
            gene_id = '.'.join(new_name.split('.')[:-1])  ### P modify this part: .
            all_splice[gene_id] = all_splice.get(gene_id, [])+[new_name]
            seq_record.id = seq_record.description = new_name
            SeqIO.write(seq_record, outfile, "fasta")
    write_splice_file(all_splice, tag+'.gene.splice')

################################################
#### Creating specific GFF file for EdgeHOG ####
################################################

def gff_parsing(inFile, tag):
    with open(tag+'.gff','w') as gfh:
        for line in open(inFile):
            line = line.strip()
            data = line.split('\t')
            if line.startswith('#') or not line:
                pass
            else:
                if data[2] == 'gene':
                    gfh.write(f"{line}\n")
                elif data[2] == 'CDS':
                    info = data[8].split(';')
                    name = [x for x in info if 'ID=' in x][0]
                    name = name.split('=')[1]
                    if '.CDS' in name:
                        name = name.split('.CDS')[0]
                    if 'cds-' in name:
                        name = name.split('cds-')[1]
                    new_id = tag+name
                    info = [f"ID={new_id}" if 'ID=' in x else x for x in info]
                    data[8] = ';'.join(info)
                    if 'protein_id=' not in data[-1]:
                        data[8]+=f";protein_id={new_id}"
                    else:
                        info = data[8].split(';')
                        info = [f"protein_id={new_id}" if 'protein_id=' in x else x for x in info]
                        data[8] = ';'.join(info)
                    line = '\t'.join(data)
                    gfh.write(f"{line}\n")


##########
## main ##
##########
parser = argparse.ArgumentParser(description="Create the splice files for OMA")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="fasta file")
parser.add_argument("-g", "--gffFile", dest="gffFile", default='no', help="gff File, when using format gff")
parser.add_argument("-f", "--format", dest="format", default='ensembl', help="format file, default=ensembl, options = ncbi, ensembl, phyto, augustus, other")
parser.add_argument("-t", "--tag", dest="tag", default='', help="tag or species name to be added to the genename")
parser.add_argument("-gff", "--gff", dest="gff", action='store_true', help="activate when you want to create a gff for edgehog and use -g option")
args = parser.parse_args()

inFile = args.inFile
gffFile = args.gffFile
formato = args.format
taxa = args.tag
gfftag = args.gff

if formato == 'ncbi':
    ## For this method, the .fa file needs to have a GeneId (GeneID=)
    print('Parsing NCBI format...')
    create_splicefile_ncbi(inFile, gffFile, taxa)
    if gfftag:
        print('Creating gff file...')
        gff_parsing(gffFile, taxa)
elif formato == 'ensembl': ## Tested
    print('Parsing ensembl format...')
    create_splicefile_ensembl(inFile, taxa)
elif formato =='phyto': ## Tested
    print('Parsing phytozome format...')
    create_splicefile_phytozome(inFile, taxa)
elif formato == 'augustus':
    create_splicefile_Augustus(inFile, "fasta")
else:
    print('Parsing no format ...') ## when you have just the protein names == genes
    create_splicefile_unspecific(inFile, taxa)
    if gfftag:
        print('Creating gff file...')
        gff_parsing(gffFile, taxa)







# file = '/home/ijulcach/projects/Legumes_nodules/Data/LOTJB/Lotus_japonicus_annos1-cds0-id_typename-nu1-upa1-add_chr0.gid63162.gff'
# outfile = open(file+'3','w')

# for line in open(file):
#     line = line.strip()
#     data = line.split('\t')
#     if line.startswith('#') or not line:
#         pass
#     else:
#         if data[0].startswith('LjG') and 'chr' in data[0]:
#             if data[2] == 'mRNA':
#                 name = data[8].split('ID=')[1].split(';')[0].split('.mRNA')[0]
#                 data[8] = data[8]+';product='+name
#             line = '\t'.join(data)
#             print(line,file=outfile)
#                 # if '.' in name:
#                 #     name = name.split('.')[0]
#                 # if '_' in name:
#                 #     name = name.split('_')[0]
#                 # print(name)
#             # else:
#             #     print(line,file=outfile)
# outfile.close()

