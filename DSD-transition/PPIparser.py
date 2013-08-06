#!/usr/sup/bin/python
"""
PPIparser.py -- this module parse input file and obtain an adjacency matrix
                (symmetric, binary, all zeros(at diagnal)

DSD version 0.5, Copyright (C) 2013, Tufts University
@author -- Mengfei Cao, mcao01@cs.tufts.edu
161 College Ave., Medford, MA 02155, USA

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA

"""

import sys
import numpy as np
import re
import collections


def GetAdj(filename, threshold):
    """
    filename - the name of input file to be parsed, which should be
               a file with one PPI at each line.

    threshold - the threshold for the existence of PPIs

    returns ppbAdj, an adjacency matrix represented as a numpy array
    """
    validpattern = re.compile('^[\w _\-.,\t"\':;]+$')
    splitpattern = re.compile('[\t ;,]+')
    numericpattern = re.compile('^[0-9. \t,\-]+')
### collect node names ###
    finfile = open(filename, 'r')
    names = {}
    index = 0
    for temp in finfile:
        temp = temp.strip('\t \n\r')
        if temp == "" or (re.search(validpattern, temp) is None):
            continue
        allwords = re.split(splitpattern, temp)
        for i in xrange(0, 2):
            if re.match(numericpattern, allwords[i]) is not None:
                temp = "Error: file inconsistent"
                temp = temp + "(possible node ID is numeric)"
                temp = temp + "       " + allwords[i]
                print >> sys.stderr, temp
                exit(1)
        if allwords[0] not in names:
            names[allwords[0]] = index
            index = index + 1
        if allwords[1] not in names:
            names[allwords[1]] = index
            index = index + 1
    finfile.close()
    names = collections.OrderedDict(sorted(names.items(), key=lambda x: x[1]))
    N = index
    ppbAdj = np.zeros((N, N))
### collect edges ###
    finfile = open(filename, 'r')
    for temp in finfile:
        temp = temp.strip('\t \n\r')
        if temp == "" or (re.search(validpattern, temp) is None):
            continue
        allwords = re.split(splitpattern, temp)
#        print temp
        i = names[allwords[0]]
        j = names[allwords[1]]
        if i != j:
            ppbAdj[i, j] = 1
            ppbAdj[j, i] = 1
    finfile.close()
 #   print ppbAdj
 #   print names
    return (ppbAdj, names)

def getTransition(infile, names):
    '''
        input:
            infile: name of file of format [protein1] [protein2] [transition probability]
            names : an array of names
        returns:
            a NumPy matrix (sorted as the input names array) of transition probabilities)

        NOTE: ! it is important that the input file (infile) has all redundant edges (A-->B and B-->A).
    '''

    N = len(names)
    tmat = np.zeros((N,N))
    
    ppifile = open(infile)

    for line in ppifile:
        parts = line.split()
        p1    = parts[0]
        p2    = parts[1]
        tprob = parts[2]
        try:
            p1index = names[p1]
            p2index = names[p2]
        except KeyError:
            continue
        tmat[p1index, p2index] = tprob

    return tmat
     
