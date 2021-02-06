#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 15:42:19 2021

@author: ijulca
"""
import argparse, glob
import os


def spider_model(group,num_seq):
    log = group+'/model.log'
    toprint = False
    if os.path.isfile(log) == True:
        with open(log, 'r') as f:
            lines = f.read().splitlines()
            if len(lines) != 0:
                last_line = lines[-1]
            else:
                last_line = 'No'
            if 'Date and Time:' in last_line:
                toprint = True
    return toprint, num_seq

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

groups = glob.glob(path+'/*/*')

outfile = open('models_all.txt', 'w')
num_seq = set([])
for g in groups:
    toprint, num_seq = spider_model(g, num_seq)
    if toprint == True:
        model = get_model(g+'/model.iqtree')
        string = g.split('/')[-1]+'\t'+model
        print(string,file=outfile)
            
outfile.close()
print('End..')