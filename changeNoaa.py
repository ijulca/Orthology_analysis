#!/usr/bin/env python3
# change non-aa letters
# Irene Julca 19/02/18

import argparse

## main
parser = argparse.ArgumentParser(description="Change the letters that are not aa.")
parser.add_argument("-i", "--inf", dest="inf", required=True, help="infile")
args = parser.parse_args()

inFile = args.inf

outfile = open(inFile + ".changed", "w")

noAA = {"B":0, "J":0, "O":0, "U":0, "Z":0}

for line in open(inFile):
	line = line.strip()
	if ">" in line:
		print(line,file=outfile)
	else:
		line = line.upper()
		string = ""
		for b in line:
			if b in noAA:
				noAA[b] += 1
				string += "X"
			else:
				string += b
		print(string,file=outfile)
		
outfile.close()

for a in noAA:
	print (a + ":\t" + str(noAA[a]))

print ("...")

