#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 15:30:38 2022

@author: ijulca
"""
import argparse

### main
parser = argparse.ArgumentParser(description="convert nexus format to newick (tree formats)")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="nexus file")
args = parser.parse_args()

nexFile = args.inFile

for line in open(nexFile):
    line = line.strip()
    if 'tree TREE1 =' in line:
        tree = ''
        data = line.split(' ')[-1].split('[')
        for d in data:
            if ']' in d:
                e = d.split(']')[1]
            else:
                e = d
            tree+=e
print(tree)
        
    
