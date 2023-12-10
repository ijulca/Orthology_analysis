#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 16:01:31 2023

@author: ijulcach
"""

import argparse
# import ete4
# from ete4 import NCBITaxa
# ncbi = NCBITaxa()
# ncbi.update_taxonomy_database()



##############
#### main ####
##############

def main(inFile):
    outfilefasta = open(inFile.replace('.fa','.format.fa'), 'w')
    outfilesp = open(inFile.replace('.fa','.species.txt'), 'w')
    for line in open(inFile):
        line = line.strip()
        if not line:
            pass
        else:
            if ">" in line:
                n = line[1:6]
                key = ">"+line.split(">")[1].split(" ")[0]+"_"+n
                print(key, file=outfilefasta)
                print(n+'\t'+line.split('|')[-1].strip().replace('[','').replace(']',''),file=outfilesp)
            else:
                print(line,file=outfilefasta)
    outfilefasta.close()
    outfilesp.close()
    print("End...")

### main
parser = argparse.ArgumentParser(description="Download papers with a key word")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="inFile oma fasta")
args = parser.parse_args()

main(args.inFile)