#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 19:42:41 2024

@author: ijulcach
"""

import argparse
import sys,os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import general_modules as gmo
import pyoma.browser.db
from pyoma.browser.models import ProteinEntry
db = pyoma.browser.db.Database('/work/FAC/FBM/DBC/cdessim2/oma/oma-browser/All.Jul2023/data/OmaServer.h5')
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
            all_seq = list()
            print(taxa +'\t'+str(taxid)+'\t'+name)
            main_isos = [ProteinEntry(db, e) for e in db.main_isoforms(taxa)]
            for iso in main_isos:
                seq = iso.sequence
                identifier = iso.omaid
                description = "" #iso.description
                record = SeqRecord(
                Seq(seq),
                id=identifier,
                description=description)
            all_seq.append(record)
            with open(outpath+ taxa +  ".fa", "w") as output_handle:
                SeqIO.write(all_seq, output_handle, "fasta")
    
    
### main
parser = argparse.ArgumentParser(description="download fasta file of hogs (root hogs)")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="list of hogs or list of mnemonic")
parser.add_argument("-p", "--outpath", dest="outpath", required=True, help="folder where to create the files")
parser.add_argument("-t", "--tag", dest="tag", required=True, help="what to download, h: hogs fasta, p: proteomes")
args = parser.parse_args()

inFile = args.inFile
outpath = args.outpath+'/'

if args.tag == 'h':
    print('Downloading HOGS...')
    hog2fasta(inFile, outpath)
elif args.tag == 'p':
    print('Downloading proteomes...')
    mnemonic2fasta(inFile, outpath)



