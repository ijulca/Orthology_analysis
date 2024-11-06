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
from pyoma.browser.models import ProteinEntry, HOG
db = pyoma.browser.db.Database('/work/FAC/FBM/DBC/cdessim2/oma/oma-browser/All.Jul2024/data/OmaServer.h5')
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyham
import collections

### Get similar hogs with shared orthologs
def get_shared_hogs(hog_id, level=None): ### level of interest or None --> rootlevel for hogid
    if 'HOG:' not in hog_id:
        hog_id = f'HOG:E0{hog_id}'
    hog = HOG(db, db.get_hog(hog_id, level=level))
    hog_memb_enrs = [e.entry_nr for e in hog.members]
    all_pw_orth_enr = set()
    for enr in hog_memb_enrs:
       pw = db.get_vpairs(enr)
       all_pw_orth_enr.update(list(int(row['EntryNr2']) for row in pw))

    all_roothogs = collections.defaultdict(int)
    for enr in all_pw_orth_enr:
        if enr in hog_memb_enrs:
           continue
        pe = ProteinEntry(db, enr)
        if pe.hog_family_nr == 0:
           continue
        all_roothogs[pe.hog_family_nr] += 1
    shared_hogs = sorted([(k,v) for k,v in all_roothogs.items()], key=lambda x: -x[1])
    with open(f"{hog_id.split(':')[1]}.shared_hogs.txt","w") as fout:
        for e,v in shared_hogs:
            print(f"{e}\t{v}",file=fout)
 
### Get the fasta file of hogs
def hog2fasta(hog_id):
    if 'HOG:' in hog_id:
        hog_id = hog_id.split(':')[1]
    hog_id = int(hog_id)
    with open(f'{hog_id}.fa', 'wt') as fout:
        for p in db.member_of_fam(int(hog_id)): #db.member_of_hog_id(hog): ###use this for full hog name
            pe = ProteinEntry(db, p)
            fout.write(f">{pe.omaid}\n") # {pe.canonicalid}\n")
            fout.write(pe.sequence)
            fout.write("\n")

### Get main isoforms
def mnemonic2fasta(mnemonic):
    genomes = [pyoma.browser.models.Genome(db, g) for g in db.db.root.Genome.read()]
    for gen in genomes:
        taxid  = gen.ncbi_taxon_id
        name = gen.sciname
        taxa = gen.uniprot_species_code
        if taxa == mnemonic:
            pepfile = open(f"{taxa}.fa",'w')
            cdsfile = open(f"{taxa}.cds.fa",'w')
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
    

### Get proteomes and splice files
def creat_spliceFile(taxa, outname):
    main_iso = [ProteinEntry(db,e) for e in db.main_isoforms(taxa)]
    main_iso = {x.omaid:[e.omaid for e in x.alternative_isoforms] for x in main_iso}
    spfile = open(outname,'w')
    for e,g in main_iso.items():
        isos = set([e]+g)
        print(';'.join(isos),file=spfile)
    spfile.close()


def mnemonic2fasta2splice(mnemonic):
    genomes = [pyoma.browser.models.Genome(db, g) for g in db.db.root.Genome.read()]
    for gen in genomes:
        taxa = gen.uniprot_species_code
        if taxa == mnemonic:
            pepfile = open(f"{taxa}.fa",'w')
            splicenameFile = f"{taxa}.splice"
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

    
def mnemonic2gff(mnemonic):
    genomes = [pyoma.browser.models.Genome(db, g) for g in db.db.root.Genome.read()]
    for genome in genomes:
        taxa = genome.uniprot_species_code
        if taxa == mnemonic:
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
 
    
def pep2hogID(inFile):
    genes = gmo.load_list(inFile)
    db_search = db.id_mapper['XRef']
    outfile = open(f"{inFile}.hogs",'w')
    for g in genes:
        info = db_search.search_xref(g)
        if len(info) == 0:
            hogs = ['None']
        else:
            ids = [x['EntryNr'] for x in info]
            hogs = set()
            for i in ids:
                try:
                    hog = db.hog_family(i)
                    hogs.add(hog)
                except Exception as e:
                    print(f"An error occurred: {e}")
            if len(hogs) == 0:
                hogs = ['None']
        for hog in hogs:
            print(f"{g}\t{hog}", file=outfile)
    outfile.close()
    

### main
parser = argparse.ArgumentParser(description="download files from OMA browser", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i", "--input", dest="input", required=True, help="hog_id or mnemonic, check tag")
parser.add_argument("-t", "--tag", dest="tag", required=True, help="""select from options:
                    hog2shogs: hog_id to similar hogs with shared orthologs,
                    hog2fasta: hog_id to fasta file,
                    pep2main_iso: mnemonic to proteomes, getting only main isoforms,
                    pep2splice: proteomes and splice forms,
                    getgff: mnemonic to gff file,
                    pep2hog: give a list of protein IDs and get hogs,""")
args = parser.parse_args()

inData = args.input

if args.tag == 'hog2fasta':
    print(f"Downloading fasta file of ...{inData}")
    hog2fasta(inData)
elif args.tag =='hog2shogs':
    print(f"Getting similar hogs...{inData}")
    get_shared_hogs(inData)
elif args.tag == 'pep2main_iso':
    print('Downloading main isoforms...{inData}')
    mnemonic2fasta(inData)
elif args.tag == 'pep2splice':
    print('Downloading proteomes and splice files...{inData}')
    mnemonic2fasta2splice(inData)
elif args.tag == 'getgff': ## use directly the infile as mnemonic, so you can parallel
    print('Downloading gff3...{inData}')
    mnemonic2gff(inData)
elif args.tag == 'pep2hog':
    print('Getting hog ids...')
    pep2hogID(inData)


