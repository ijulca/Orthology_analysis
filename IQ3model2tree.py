#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 18:48:24 2021

@author: ijulca
"""
import argparse, glob
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM
import orthology_modules as OM
import general_modules as gmo

iqtree = '/home/irene.julca/Programs/iqtree-2.1.2-Linux/bin/iqtree2'

def get_model(inFile):
    for line in open(inFile):
        line = line.strip()
        if 'Best-fit model according' in line:
            model = line.split(':')[1]
            if ' ' in model:
                model = model.replace(' ','')
    return model


### main
parser = argparse.ArgumentParser(description="get the model and create the jobs for the tree reconstruction")
parser.add_argument("-p", "--path", dest="path", default='./Data/', help="path where to create the folders")
args = parser.parse_args()

path = args.path

modelFiles = glob.glob(path+'/*/*/model.iqtree')

outfile = open('tree.jobs', 'w')
for f in modelFiles:
    model = get_model(f)
    alg = f.split('model')[0]+f.split('/')[-2]+'.alg.clean'
    pref = f.split('model')[0]+'genetree'
    cmd = iqtree + ' -s '+alg +' -m '+model+' --prefix '+pref
    print(cmd,file=outfile)
outfile.close()
print('End...')