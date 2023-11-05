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

## main
parser = argparse.ArgumentParser(description="Change the letters that are not aa and remove stop codon.")
parser.add_argument("-i", "--inFile", dest="inFile", required=True, help="infile")
args = parser.parse_args()

inFile = args.inFile

seqs = GM.load_sequences(inFile)
seqs = remove_stop(seqs)


outfile = open(inFile + ".changed", "w")

noAA = {"B":0, "J":0, "O":0, "U":0, "Z":0,'*':0}

for s in seqs:
    pep = ''
    for b in seqs[s]:
        if b in noAA:
            noAA[b] += 1
            pep += "X"
        else:
            pep += b
    GM.print_sequence(s,pep,outfile)
		
outfile.close()

for a in noAA:
	print (a + ":\t" + str(noAA[a]))

print ("...")

