#!/bin/bash

f=("$@")
n=${f%%.*}

mafft $f >$n.mafft
trimal -in $n.mafft -out $n.phy -phylip -gt 0.5
iqtree2 -s $n.phy -m MFP -B 1000 -T 12
