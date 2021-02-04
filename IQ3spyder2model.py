#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 10:07:42 2021

@author: ijulca
"""
import argparse, glob
import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import general_modules as gmo

iqtree = '/home/irene.julca/Programs/iqtree-2.1.2-Linux/bin/iqtree2'
algtool = '/home/irene.julca/Programs/git_repository/phylome_analysis/salva_scripts/12.IndividualStep.Alignments.py'

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
            elif 'ERROR: There must be at least 3 sequences' in last_line:
                toprint = True
        if toprint == False:
            print('unfinished job:',log)
            num_seq.add(log)
            cmd = 'rm '+group+'/model.*'
            gmo.run_command(cmd)
    return toprint, num_seq

def spider_alg(path):
    log = folder +'/'+folder.split('/')[-1]+'.log'
    toprint = False
    if os.path.isfile(log) == True:
        with open(log, 'r') as f:
            lines = f.read().splitlines()
            if 'STEP	Multipple Sequence Alignment	END' in lines[-3]:
                toprint = True
    return toprint
        
    

parser = argparse.ArgumentParser(description="spyder for aligment and model")
parser.add_argument("-p", "--path", dest="path", default='./Data/', help="path where to create the folders")
parser.add_argument("-t", "--tag", dest="tag", required=True, help="tag to get unfinished jobs for m=model, a=aligment")
parser.add_argument("-c", "--config", dest="config", default = 'no',help="give for tag: a=aligment")
args = parser.parse_args()

path = args.path
tag = args.tag
configFile = args.config


groups = glob.glob(path+'/*/*')

if tag == 'm':
    print('Spider for models...')
    outfile = open('model.spyder.job', 'w')
    num_seq = set([])
    for g in groups:
        toprint, num_seq = spider_model(g, num_seq)
        if toprint == False:
            cmd = iqtree+' -s '+g+'/'+g.split('/')[-1]+'.alg.clean -m MF --prefix '+g+'/model'#' -T 4'
            print(cmd, file=outfile)
    outfile.close()

    print('Unfinished jobs:', len(num_seq))
elif tag == 'a':
    print('Spider for aligment...')
    outfile = open('alg.spyder.job', 'w')
    for folder in groups:
        toprint = spider_alg(folder)
        if toprint == False:
            cmd = algtool+' -c '+configFile+' -i '+folder +'/'+folder.split('/')[-1]+'.fa --cds '+folder +'/'+folder.split('/')[-1]+'.cds -p '+folder +'/'+folder.split('/')[-1]
            print(cmd,file=outfile)
    outfile.close()
print('End')
    
    
