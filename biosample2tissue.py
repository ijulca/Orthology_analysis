#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 15:56:52 2021

@author: irene
"""
import argparse, glob

#### you can use the command to get the SAMN file:
## for i in `awk '{print $1}' biosample.txt`; do wget -q -O - https://www.ncbi.nlm.nih.gov/biosample/?term=$i >$i.info; done

def get_pref(inFile):
    if '/' in inFile:
        pref = inFile.split('/')[-1]
    else:
        pref = inFile
    pref = pref.split('.')[0]
    return pref

def get_sampleFromFile(inFile, key):
    sample, source ='','No'
    for line in open(inFile):
        line = line.strip()
        if key in line:
            source = key.replace('<th>','')
            sample = line.split(key)[1].split('</td')[0].split('>')[-1].strip()
            if sample == 'missing':
                    sample = ''
    return sample,source

def get_sample_name(sample):
    info = {'Leaf':['Leaf', 'leaf', 'leaves', 'Leaves','leave', 'LEAF','petiole','rosette','foliage'], 
            'Flower':['flower','tassels', 'inflorescence', 'panical', 'Panicle', 'tassel', 'floral','ear', 'Florets','spike', 'tepals', 'Sepal','Ears','Flower','Inflorescences',
                      'panicles','carpel','petal','Peduncle','Petal','Floral','Carpel'], 
            'Stem':['stem', 'Bark', 'Stem', 'cambium','xylem', 'phloem', 'tuber', 'steam', 'vascular', 'Vascular', 'Tuber','Internode','wood','cambial','cork',
                    'Xylem','bark','internode'], 
            'Thallus':['thallus'], 
            'Meristem':['shoot', 'Shoot','buds','bud','Buds','Bud'], 
            'Fruit':['fruit','furits','Fruit','Berry','berry', 'Endocarp', 'pericarp','mesocarp', 'berries'],
            'Seed':['Seed', 'seed', 'Kernel','nut', 'cotyledon','kernel','Cotyledon','grain','ENDOSPERM','endosperm','Pericarp'], 
            'Male':['stamen', 'pollen','anther','Microspore'], 
            'Female':['ovule','pistil','ovary','style','megagametophytes','Ovule','nucellus'],
            'Root':['root', 'Root','ROOT'],
            'Embryo':['embryo', 'zygote','Embryo','embryos']}
    gnames = ['Leaf','Flower','Meristem','Stem', 'Thallus','Fruit', 'Seed', 'Male', 'Female','Root','Embryo']
    toprint = False
    sname = 'Other'
    for name in gnames:
        if toprint == False:
            for e in info[name]:
                if e in sample:
                    toprint = True
                    sname = name
    return sname

### main
parser = argparse.ArgumentParser(description="get the name of the tissue from the SAMN file (sample name)")
parser.add_argument("-p", "--path", dest="path", required=True, help="path to the folder containing the SAMN files")
args = parser.parse_args()

path = args.path

outfile = open('sample_names.txt','w')

files = glob.glob(path+'/*')
print('number of files to analize:', len(files))

for inFile in files:
    pref = get_pref(inFile)
    sources = ['<th>tissue','<th>source name','<th>sample name','<th>isolation source']
    sample, source ='','No'
    for s in sources:
        if sample =='':
            sample, source = get_sampleFromFile(inFile, s)
            
    if sample == '':
        print(pref)
    else:
        string = pref +'\t'+sample+'\t'+get_sample_name(sample)+'\t'+source
        print(string,file=outfile)

outfile.close()
print('End...')