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
codeml_py = "/home/irene.julca/Programs/git_repository/Codeml_auto/run_codeml.py"
codemlPair = '/home/irene.julca/Programs/git_repository/Codeml_auto/run_alg2dnds.py'

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
        
def spider_codemlF(folder):
    log = folder+'/mlcM1'
    toprint = False
    if os.path.isfile(log) == True:
        with open(log, 'r') as f:
            lines = f.read().splitlines()
            if 'Time used:' in lines[-1]:
                for l in lines:
                    if 'w ratios as labels for TreeView' in l:
                        toprint = True
            else:
                print('ERROR last line....', folder)
    return toprint

def change_master_codeml(masterctl,outpref, alg, tree):
	outfile = open(outpref, "w")
	for line in open(masterctl):
		if not line or line.startswith("#"): 
			continue
		else:
			line = line.strip()
			if "seqfile" in line:
				line = line.replace("seq.phy", str(alg))
			elif "treefile" in line:
				line = line.replace("tree.txt", str(tree))
			print(line, file=outfile)
	outfile.close()
    
def spider_codemlP(folder):
    log = folder+'dnds_pairs.txt'
    toprint = False
    if os.path.isfile(log) == True:
        alg = folder+'/'+folder.split('/')[-1]+'.alg.clean_cds'
        with open(alg, 'r') as f:
            lines = f.read().splitlines()
            num = lines[0].split(' ')[1]
        items = [x for x in range(0,num)]
        pairs = [(items[i],items[j]) for i in range(len(items)) for j in range(i+1, len(items))]
        r_pair = 0
        for line in open(log):
            line = line.strip()
            r_pair +=1
        if pairs == r_pair:
            toprint = True
    return toprint

parser = argparse.ArgumentParser(description="spyder for aligment and model")
parser.add_argument("-p", "--path", dest="path", default='./Data/', help="path where to create the folders")
parser.add_argument("-t", "--tag", dest="tag", required=True, help="tag to get unfinished jobs for m=model, a=aligment,c1=codemlFree,c2=codemlPairwise")
parser.add_argument("-c", "--config", dest="config", default = 'no',help="give config file for tag: a=aligment, c1=codemlFree, c2=codemlPairwise")
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
elif tag == 'c1':
    print('Spider for codeml freeratio')
    outfile = open('codemlF.spyder.job', 'w')
    for folder in groups:
        toprint = spider_codemlF(folder)
        if toprint == False:
            if '/' in configFile:
                outpref = folder + '/'+configFile.split('/')[-1]
            else:
                outpref = folder + '/'+configFile
            
            group = folder.split('/')[-1]
            alg = group+'.alg.clean_cds.phy'
            treefile = 'genetree.treefile'
            change_master_codeml(configFile,outpref, alg, treefile)
            cmd = codeml_py+' -p '+folder+'/ -c '+outpref
            print(cmd,file=outfile)
    outfile.close()
elif tag == 'c2':
    print('Spider for codeml pairwise')
    outfile = open('codemlP.spyder.job', 'w')
    for folder in groups:
        toprint = spider_codemlP(folder)
        if toprint == False:
            cmd = codemlPair+' -p '+folder+' -c '+configFile
            print(cmd,file=outfile)
    outfile.close()
    
print('End')
    
    
