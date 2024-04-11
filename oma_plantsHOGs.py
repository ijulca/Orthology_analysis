#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 16:07:43 2024

@author: ijulcach
"""

import sys,os
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/modules_py/')
import general_modules as gmo
from omadb import Client #https://dessimozlab.github.io/pyomadb/build/html/
c = Client()


def get_root_hogs(outFile):
    outfile = open(outFile, 'w')
    table = c.hogs.list()
    for e in table:
        hog = e.hog_id
        level = e.level
        if level == 'root':
            print(hog,file=outfile)
    outfile.close()

    
##############
#### Main ####
##############
path = '/home/ijulcach/projects/Land_Plants/HOGS/'
omaDb = '/home/ijulcach/Programs/DataBases/LUCA.h5'
list_hogFile = path+'hog_all.txt'
get_root_hogs(list_hogFile)
