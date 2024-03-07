#!/usr/bin/env python
# coding: utf-8



from omadb import Client #https://dessimozlab.github.io/pyomadb/build/html/
c = Client()
from ete4 import Tree
import re
import argparse
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM


def get_xref(omaid):
    '''Get the source id cross-reference from an omaid'''
    xrefs = c.entries[omaid].xref
    sourceid = [x for x in xrefs if x['source']=='SourceID'][0]['xref']
    return sourceid

def change_treeIDs(treeFile):
    print('Changing Ids in the tree and generating a new newick file...')
    # Loads a tree structure from a newick string
    tree = Tree(open(treeFile)) 

    #get new leaf names with cross-reference ids

    pattern = re.compile(r'\b[A-Z]{5}\d{5}\b')

    for leaf in tree:
        geneid = leaf.name.split("_")[0].replace("'","")
        species = leaf.name.split("_")[1].replace("'","")
        if pattern.match(geneid):
            xref = get_xref(geneid)
            leaf.name = xref + "_" + species

    #newick
    tree.write(outfile=treeFile+'.converted.nw', parser=2)

def change_fastaIDs(inFile):
    print('Changing fasta file...')
    seq = GM.load_sequences(inFile)
    outfile = open(inFile+'.converted.fa','w')
    pattern = re.compile(r'\b[A-Z]{5}\d{5}\b')
    for s in seq:
        geneid = s.split("_")[0].replace("'","")
        species = s.split("_")[1].replace("'","")
        name = s
        if pattern.match(geneid):
            name = get_xref(geneid) + '_'+species
        GM.print_sequence(name,''.join(seq[s]),outfile)
    outfile.close()
        
    
### main
parser = argparse.ArgumentParser(description="change oma ids to original")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="tree file")
parser.add_argument("-t", "--tag", dest="tag", required=True, help="tree, fasta")
args = parser.parse_args()

inFile = args.inFile
tag = args.tag

if tag == 'tree':
    change_treeIDs(inFile)
elif tag == 'fasta':
    change_fastaIDs(inFile)

print('End...')

