#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 23:13:49 2022

@author: ijulca
"""
import argparse
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

## main
parser = argparse.ArgumentParser(description="Ks plot")
parser.add_argument("-p", "--path", dest="path", required=True, help="path to the kaks files")
args = parser.parse_args()

files = glob.glob(args.path+'/*.kaks')

spe, ks = [], []
for f in files:
    name = f.split('/')[-1].split('.')[0]
    print(name)
    for line in open(f):
        line = line.strip()
        if not line or line.startswith('#'):
            pass
        else:
            data = line.split('\t')
            val = float(data[-1])
            if val >0 and val <2:
                ks.append(val)
                spe.append(name)

df = pd.DataFrame(list(zip(spe, ks)), columns =['Spe', 'Ks'])
ax = sns.displot(data=df, kind="kde", x="Ks", hue="Spe")

plt.savefig("Ks_plot.svg")
plt.show()
print("end...")