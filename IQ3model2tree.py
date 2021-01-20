#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 18:48:24 2021

@author: ijulca
"""
import argparse, glob
#import sys, os
#sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')


iqtree = '/home/irene.julca/Programs/iqtree-2.1.2-Linux/bin/iqtree2'

def get_model(inFile):
    for line in open(inFile):
        line = line.strip()
        if 'Best-fit model according' in line:
            model = line.split(':')[1]
            if ' ' in model:
                model = model.replace(' ','')
    return model

def check_genetree(f):
    folder = f.split('model')[0]
    files = glob.glob(folder+'/*')
    toprint = False
    for e in files:
        if 'genetree' in e:
            toprint = True
    return toprint

### main
parser = argparse.ArgumentParser(description="get the model and create the jobs for the tree reconstruction")
parser.add_argument("-p", "--path", dest="path", default='./Data/', help="path where to create the folders")
args = parser.parse_args()

path = args.path

modelFiles = glob.glob(path+'/*/*/model.iqtree')
g=0
outfile = open('genetree.job', 'w')
for f in modelFiles:
    toprint = check_genetree(f)
    if toprint == 'False':
        model = get_model(f)
        alg = f.split('model')[0]+f.split('/')[-2]+'.alg.clean'
        pref = f.split('model')[0]+'genetree'
        cmd = iqtree + ' -s '+alg +' -m '+model+' --prefix '+pref +' -T 2 -B 1000' ### faster -B 1000, but better -b 100
        print(cmd,file=outfile)
    else:
        g+=1
outfile.close()
print('orthogroups that have genetree:',g)
print('orthogroups to be analysed:',len(modelFiles)-g)