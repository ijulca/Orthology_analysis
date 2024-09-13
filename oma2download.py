#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 19:42:41 2024

@author: ijulcach
"""

import argparse
import sys,os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM
import general_modules as gmo
import pyoma.browser.db
from pyoma.browser.models import ProteinEntry
db = pyoma.browser.db.Database('/work/FAC/FBM/DBC/cdessim2/oma/oma-browser/All.Jul2024/data/OmaServer.h5')
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyham


def hog2fasta(inFile, outpath):
    gmo.create_folder(outpath)
    hogs = gmo.load_list(inFile)
    print('Number of HOGs:', len(hogs))
    if len(hogs)>1000:
        for i in range(0,len(hogs), 1000):
            z = i + 1000
            if z >len(hogs):
                z = len(hogs)
            outdir1 = outpath+str(i+1)+'-'+str(z)+'/'
            gmo.create_folder(outdir1)
            for j in range(i,z):
                hog = hogs[j]
                name = hog.split(':')[1]
                outdir2 = outdir1+name+'/'
                if os.path.isdir(outdir2) == False:
                    gmo.create_folder(outdir2)
                    outFile = outdir2+name+'.fasta'
                    with open(outFile, 'wt') as fout:
                        for p in db.member_of_hog_id(hog):
                            pe = ProteinEntry(db, p)
                            fout.write(f">{pe.omaid} {pe.canonicalid}\n")
                            fout.write(pe.sequence)
                            fout.write("\n\n")
                else:
                    print('Done for..', hog)
    else:
        for hog in hogs:
            with open(f'hog_{hog}.fa', 'wt') as fout:
                for p in db.member_of_hog_id(hog):
                    pe = ProteinEntry(db, p)
                    fout.write(f">{pe.omaid} {pe.canonicalid}\n")
                    fout.write(pe.sequence)
                    fout.write("\n\n")


def mnemonic2fasta(inFile, outpath):
    gmo.create_folder(outpath)
    genomes = [pyoma.browser.models.Genome(db, g) for g in db.db.root.Genome.read()]
    names = gmo.load_list(inFile)
    for gen in genomes:
        taxid  = gen.ncbi_taxon_id
        name = gen.sciname
        taxa = gen.uniprot_species_code
        if taxa in names:
            pepfile = open(outpath+taxa+'.fa','w')
            cdsfile = open(outpath+taxa+'.cds.fa','w')
            main_isos = [ProteinEntry(db, e) for e in db.main_isoforms(taxa)]
            print(taxa +'\t'+str(taxid)+'\t'+name+'\t'+str(len(main_isos)))
            for iso in main_isos:
                seq = iso.sequence
                cds = iso.cdna
                identifier = iso.omaid #iso.canonicalid
                #description = "" #iso.description
                GM.print_sequence(identifier, ''.join(seq), pepfile)
                GM.print_sequence(identifier, ''.join(cds), cdsfile)
            pepfile.close()
            cdsfile.close()
    

def creat_spliceFile(taxa, outname):
    main_iso = [ProteinEntry(db,e) for e in db.main_isoforms(taxa)]
    main_iso = {x.omaid:[e.omaid for e in x.alternative_isoforms] for x in main_iso}
    spfile = open(outname,'w')
    for e,g in main_iso.items():
        isos = set([e]+g)
        print(';'.join(isos),file=spfile)
    spfile.close()


def mnemonic2fasta2splice(inFile, outpath):
    gmo.create_folder(outpath)
    genomes = [pyoma.browser.models.Genome(db, g) for g in db.db.root.Genome.read()]
    names = gmo.load_list(inFile)
    for gen in genomes:
        taxa = gen.uniprot_species_code
        if taxa in names:
            print(taxa, '...')
            pepfile = open(outpath+taxa+'.fa','w')
            splicenameFile = outpath+taxa+'.splice'
            creat_spliceFile(taxa, splicenameFile)
            genes = [ProteinEntry(db,e) for e in db.all_proteins_of_genome(taxa)]
            for g in genes:
                seq = g.sequence
                identifier = g.omaid #iso.canonicalid
                GM.print_sequence(identifier, ''.join(seq), pepfile)
            pepfile.close()


def get_chr_entries(g):
    genes = db.main_isoforms(g.uniprot_species_code)
    chr_genes = [ProteinEntry(db, z) for z in genes]
    return chr_genes

def write_gff3(genes, genome, gff):
    with open(gff, 'wt') as gffh:
        for g in genes:
            source_id = g.omaid #[x['xref'] for x in g.xrefs if x['source'] == 'SourceID'][0]
            source_ac = [x['xref'] for x in g.xrefs if x['source'] == 'SourceAC'][0]
            strand = "+" if g.strand > 0 else "-"
            gffh.write(f"{g.chromosome}\t{genome.release.split(';')[0]}\tgene\t{g.locus_start}\t{g.locus_end}\t\t{strand}\t0\tID={source_id}\n")
            for ex in g.exons.as_list_of_dict():
                gffh.write(f"{g.chromosome}\t{genome.release.split(';')[0]}\tCDS\t{ex['start']}\t{ex['end']}\t\t{strand}\t0\tID=CDS:{source_ac};Parent={source_id};protein_id={source_id}\n")

    
def mnemonic2gff(inFile):
    genomes = [pyoma.browser.models.Genome(db, g) for g in db.db.root.Genome.read()]
    for genome in genomes:
        taxa = genome.uniprot_species_code
        if taxa == inFile:
            main_iso = [ProteinEntry(db,e) for e in db.main_isoforms(taxa)]    
            with open(taxa+'.gff','w') as gffh:
                for g in main_iso:
                    source_id = [x['xref'] for x in g.xrefs if x['source'] == 'SourceID'][0]
                    # source_ac = [x['xref'] for x in g.xrefs if x['source'] == 'SourceAC'][0]
                    isof = [g]+g.alternative_isoforms           
                    strand = "+" if g.strand > 0 else "-"
                    gffh.write(f"{g.chromosome}\t{genome.release.split(';')[0]}\tgene\t{g.locus_start}\t{g.locus_end}\t\t{strand}\t0\tID={source_id}\n")
                    for p in isof:
                        pep_id = p.omaid
                        strand = "+" if p.strand > 0 else "-"
                        for ex in p.exons.as_list_of_dict():
                            gffh.write(f"{p.chromosome}\t{genome.release.split(';')[0]}\tCDS\t{ex['start']}\t{ex['end']}\t\t{strand}\t0\tID=CDS:{pep_id};Parent={source_id};protein_id={pep_id}\n")
 

### main
parser = argparse.ArgumentParser(description="download fasta file of hogs (root hogs)")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="list of hogs or list of mnemonic")
parser.add_argument("-p", "--outpath", dest="outpath", default='no', help="folder where to create the files")
parser.add_argument("-t", "--tag", dest="tag", required=True, help="what to download, h: hogs fasta, p: proteomes mainiso, s:proteomes and splice forms, g: gff")
args = parser.parse_args()

inFile = args.inFile
outpath = args.outpath+'/'


if args.tag == 'h':
    print('Downloading HOGS...')
    hog2fasta(inFile, outpath)
elif args.tag == 'p':
    print('Downloading proteomes...')
    mnemonic2fasta(inFile, outpath)
elif args.tag == 's':
    print('Downloading proteomes and splice files...')
    mnemonic2fasta2splice(inFile, outpath)
elif args.tag == 'g': ## use directly the infile as mnemonic, so you can parallel
    print('Downloading gff3...')
    mnemonic2gff(inFile)


