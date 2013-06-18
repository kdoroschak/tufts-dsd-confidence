#!/usr/bin/env python2.7
"""
GetAdjacency.py

"""

import sys
import numpy as np
import re
import collections
import PPIparser

infilename = raw_input("Input file name: ")
matrixfilename = raw_input("Matrix file name: ")
protfilename = raw_input("Protein file name: ")

(ppbAdj, names) = PPIparser.GetAdj(infilename, 0)
allnames = names.keys()

N = len(names)

fout = open(matrixfilename, 'w')
fnameout = open(protfilename, 'w')

for i in xrange(0, N):
    fnameout.write(allnames[i])
    fnameout.write('\n')
    temp = ''
    for j in xrange(0, N):
        if ppbAdj[i,j] == 0:
            temp += '0'
        else:
            temp += '1'
    fout.write(temp+'\n')

fout.close()
            
