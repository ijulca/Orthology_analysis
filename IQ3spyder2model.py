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

def check_spider(group,num_seq):
    log = group+'/model.log'
    toprint = False
    if os.path.isfile(log) == True:
        with open(log, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
            if 'Date and Time:' in last_line:
                toprint = True
            elif 'ERROR: There must be at least 3 sequences' in last_line:
                toprint = True
        if toprint == False:
            print('unfinished job:',log)
            num_seq.add(log)
            #cmd = 'rm '+group+'/model.*'
            #gmo.run_command(cmd)
    return toprint, num_seq


parser = argparse.ArgumentParser(description="get the model and create the jobs for the tree reconstruction")
parser.add_argument("-p", "--path", dest="path", default='./Data/', help="path where to create the folders")
args = parser.parse_args()

path = args.path

groups = glob.glob(path+'/*/*')

outfile = open('model.job', 'w')
num_seq = set([])
for g in groups:
    toprint, num_seq = check_spider(g, num_seq)
    if toprint == False:
        cmd = iqtree+' -s '+g+'/'+g.split('/')[-1]+'.alg.clean -m MF --prefix '+g+'/model'
        print(cmd, file=outfile)
outfile.close()

print('Unfinished jobs:', len(num_seq))
print(num_seq)
print('End')
    
    
