#!/usr/bin/env python3
# change non-aa letters
# Irene Julca 19/02/18

import argparse

import sys, os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import genome_modules as GM

def remove_stop(seqs):
    new_seq = {}
    i = 0
    for s in seqs:
        if seqs[s][-1] == '*' or seqs[s][-1]=='.':
            new_seq[s] = seqs[s][:-1]
            i+=1
        else:
            new_seq[s] = seqs[s]
    print('stop codon removed: ', i, ' from the total:', len(seqs))
    return new_seq

def remove_noAA(seqs):
    noAA = {"B":0, "J":0, "O":0, "U":0, "Z":0,'*':0,".":0}
    new_seq = {}
    for s in seqs: 
        pep = seqs[s]
        for a in noAA:
            n = pep.count(a)
            noAA[a] += n
            if n != 0:
                pep = pep.replace(a,'X')
        new_seq[s] = pep
    print('Letters no aa found:', noAA)
    return new_seq

def remove_pepX(seqs, outname):
    prot = set()
    outfile = open(outname, 'w')
    for s in seqs:
        pep = set(seqs[s])
        if len(pep) == 1 and 'X' in pep:
            prot.add(s)
        else:
            GM.print_sequence(s,''.join(seqs[s]),outfile)
    outfile.close()
    print(f'number of sequences with X aa: {len(prot)}')
    print(prot)
    return prot
            
def change_splitFile(prot, spliceFile):
    outfile = open(spliceFile+".changed", "w")
    for line in open(spliceFile):
        line = line.strip()
        data = [x for x in line.split(';') if x not in prot]
        if len(data)>0:
            line = ';'.join(data)
            print(line,file=outfile)
    outfile.close()
        
        
        

## main
parser = argparse.ArgumentParser(description="Clean proteomes: remove proteins with X only and change the letters that are not aa and remove stop codon")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="infile")
parser.add_argument("-s", "--spliceFile", dest="spliceFile", default='no', help="splice file, default=no")
args = parser.parse_args()

inFile = args.inFile
spliceFile = args.spliceFile

if __name__=="__main__":
    seqs = GM.load_sequences(inFile)
    seqs = remove_stop(seqs)
    seqs = remove_noAA(seqs)
    prot = remove_pepX(seqs, inFile+".changed")
    if spliceFile != 'no':
        print('Changing spliceFile...')
        change_splitFile(prot, spliceFile)
