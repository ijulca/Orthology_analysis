#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 18:48:24 2021

@author: ijulca
"""
import argparse, glob
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import general_modules as gmo


iqtree = '/home/irene.julca/Programs/iqtree-1.6.12-Linux/bin/iqtree'

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
    toprint = 0
    for e in files:
        if 'genetree' in e:
            toprint = 1
    return toprint

def check_spider(f, num_seq):
    folder = f.split('model')[0]
    log = folder+'/genetree.log'
    toprint = 0
    if os.path.isfile(log) == True:
        with open(log, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
            if 'Date and Time:' in last_line:
                toprint = 1
            elif 'ERROR: It makes no sense to perform bootstrap with less than 4 sequences.' in last_line:
                toprint = 2
                num_seq.add(log)
        if toprint == 0:
            print('unfinished job:',log)
        if toprint == 0 or toprint == 2:
            cmd = 'rm '+folder+'/genetree.*'
            gmo.run_command(cmd)
            
    return toprint, num_seq

### main
parser = argparse.ArgumentParser(description="get the model and create the jobs for the tree reconstruction")
parser.add_argument("-p", "--path", dest="path", default='./Data/', help="path where to create the folders")
parser.add_argument("-s", "--spider", dest="spider", action='store_true', help="if activated will readd the genetree.log files and capture the unfinished jobs")
parser.add_argument("-t", "--threads", dest="threads", default='1', help="number of threads. Default=1")
args = parser.parse_args()

path = args.path
spider = args.spider
threads = args.threads
print('spider',spider)

modelFiles = glob.glob(path+'/*/*/model.iqtree')
g=0
num_seq = set([])
outfile = open('genetree.job', 'w')
for f in modelFiles:
    if spider == True:
        toprint, num_seq = check_spider(f,num_seq)
    else:
        toprint = check_genetree(f)
    if toprint == 0:
        model = get_model(f)
        alg = f.split('model')[0]+f.split('/')[-2]+'.alg.clean'
        pref = f.split('model')[0]+'genetree'
        cmd = iqtree + ' -s '+alg +' -m '+model+' -pre '+pref +' -nt '+threads+' -bb 1000' ### faster -B 1000, but better -b 100
        print(cmd,file=outfile)
    elif toprint == 2:
        model = get_model(f)
        alg = f.split('model')[0]+f.split('/')[-2]+'.alg.clean'
        pref = f.split('model')[0]+'genetree'
        cmd = iqtree + ' -s '+alg +' -m '+model+' --prefix '+pref  ### no bootstrap for less than 4 seq
        print(cmd,file=outfile)
    else:
        g+=1
outfile.close()
print('orthogroups that have genetree:',g-len(num_seq))
print('orthogroups that have <4 seq:', len(num_seq))#, num_seq)
print('orthogroups to be analysed:',len(modelFiles)-g)