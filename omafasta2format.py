#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 16:01:31 2023

@author: ijulcach
"""

import argparse
# import sys,os
# sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
# import genome_modules as GM



##############
#### main ####
##############

def main(inFile):
    outfilefasta = open(inFile.replace('.fa','.format.fa'), 'w')
    outfilesp = open(inFile.replace('.fa','.species.txt'), 'w')
    names = set()
    toprint = False
    for line in open(inFile):
        line = line.strip()
        if not line:
            pass
        else:
            if ">" in line:
                n = line[1:6]
                key = ">"+line.split(">")[1].split(" ")[0]+"_"+n
                if key not in names:
                    toprint = True
                    names.add(key)
                    print(key, file=outfilefasta)
                    print(n+'\t'+line.split('|')[-1].strip().replace('[','').replace(']',''),file=outfilesp)
                else:
                    toprint = False
            else:
                if toprint == True:
                    print(line,file=outfilefasta)
    outfilefasta.close()
    outfilesp.close()
    print("End...")

### main
parser = argparse.ArgumentParser(description="Download papers with a key word")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="inFile oma fasta, remove duplicated entries")
args = parser.parse_args()

main(args.inFile)