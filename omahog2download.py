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
from omadb import Client #https://dessimozlab.github.io/pyomadb/build/html/
c = Client()

def get_list_hogs(inFile):
    hogs = set()
    for line in open(inFile):
        line = line.strip()
        data = line.split('.')[0]
        hogs.add(data)
    print('Number of hogs',len(hogs))
    return list(hogs)

def get_fasta_hogs(inFile, path):
    omaw = 'https://oma-stage.vital-it.ch/oma/hog/'
    gmo.create_folder(outpath)
    hogs = get_list_hogs(inFile)
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
                gmo.create_folder(outdir2)
                outFile = outdir2+name+'.fasta'
                level = c.hogs[hog].level
                page = omaw+ hog+'/'+level+'/fasta/'
                cmd = 'wget '+page+' -O '+outFile
                print(cmd)
    else:
        for g in hogs:
            level = c.hogs[g].level
            page = 'https://oma-stage.vital-it.ch/oma/hog/'
            page += g+'/'+level+'/fasta/'
            cmd = 'wget '+page + ' -O '+g.split(':')[1]
            print(cmd)
            # gmo.run_command(cmd)


### main
parser = argparse.ArgumentParser(description="download fasta file of hogs (root hogs)")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="list of hogs")
parser.add_argument("-p", "--outpath", dest="outpath", required=True, help="folder where to create the files")
args = parser.parse_args()

inFile = args.inFile
outpath = args.outpath

get_fasta_hogs(inFile, outpath)

