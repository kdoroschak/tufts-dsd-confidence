#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""

@author: mcao01

Used to conduct k-fold Majority Voting

usage: DSDmv.py [-h] [-l LABEL] [-r RDINDEX] [-k K] [-o OUTFILE] [-d DSDFILE]
                [-t NEIGHBOR] [-m {0,1,2}]
                infile

parses PPIs from infile and do majority voting

positional arguments:
  infile                read PPIs from infile

optional arguments:
  -h, --help            show this help message and exit
  -l LABEL, --label LABEL
                        label list file
  -r RDINDEX, --rdindex RDINDEX
                        file of random indeces for labeled nodes
  -k K                  k-fold cross validation, 2 by default
  -o OUTFILE, --outfile OUTFILE
                        output prediction file name, tab delimited tables
  -d DSDFILE, --dsdfile DSDFILE
                        dsd file name
  -t NEIGHBOR, --neighbor NEIGHBOR
                        number of lowest DSD neighbors used
  -m {0,1,2,3}, --mode {0,1,2,3}
                        majority voting modes: 0 for ordinary majority voting;
                        1 for unweighted DSD; 2 for weighted DSD; 3 for weighted DSD in multiple iterations. Type 0 by
                        default

"""


import re
import numpy as np
import collections
#from optparse import OptionParser
import argparse
import myparser
import mvote
import sys


temp = "parses PPIs from infile and do majority voting"
parser = argparse.ArgumentParser(description=temp)

parser.add_argument("infile", help="read PPIs from infile")
parser.add_argument("-l", "--label",
                    help="label list file")
parser.add_argument("-r", "--rdindex",
                    help="file of random indeces for labeled nodes")                    
parser.add_argument("-k", default=2, help="k-fold cross validation,"
                    + " 2 by default", type=int)
parser.add_argument("-o", "--outfile", default="test",
                    help="output prediction file name,"
                    + " tab delimited tables")
parser.add_argument("-d", "--dsdfile", help="dsd file name")
parser.add_argument("-t", "--neighbor",
                    default=10, help="number of lowest DSD neighbors used",
                    type=int)
parser.add_argument("-m", "--mode", default=0,
                    help="majority voting modes:"
                    + " 0 for ordinary majority voting; 1 for unweighted"
                    + " DSD; 2 for weighted DSD."
                    + " Type 0 by default", type=int, choices=[0, 1, 2, 3])
parser.add_argument("-p", "--completeProteinList", default="test",
                    help="complete list of all proteins to be predicted,"
                    + " one protein per line")

#args = ['-l', "template//firstlevellabel.list", 
#        'template//BakerAdjacency.list',
#        '-k', '2', '-o', 'haha.test',
#        '-d', 'template//ExactDSD.list',
#        '-r', 'template//firstLevelRandIndex.txt',
#        '-m', '2', '-t', '10']
#options = parser.parse_args(args)
options = parser.parse_args()

#### Phase 1: Parse All Input Files
ppbAdj = myparser.parsePPI(options.infile)
N = len(ppbAdj[:,0])
ppbLabel = myparser.parseLabel(options.label)
m = len(ppbLabel[0,:]) - 1
numLabels = m
pnRD = myparser.parseRDIndex(options.rdindex, ppbLabel)
pnFoldIndex = myparser.GetFoldIndex(pnRD, N, options.k)

if options.mode != 0 and options.mode != 3:
    ppfDSD = myparser.parseDSD(options.dsdfile)
    #### print ppfDSD[1,3], ppfDSD[5,11]

#### Phase 2: Conduct Majority Voting

if options.mode == 0:
    prediction = mvote.ordinaryMV(ppbAdj, ppbLabel, pnFoldIndex, pnRD)
elif options.neighbor <= 0 or options.neighbor >= N/2:
    options.neighbor == 10
    print >> sys.stderr, 'the setting for top DSD neighbors is invalid,'
    print >> sys.stderr, 'change to 10 instead by default.\n'
if options.mode == 1:
    prediction = mvote.DSDUnweightMV(ppfDSD, ppbLabel, pnFoldIndex, pnRD, options.neighbor)
elif options.mode == 2:
    prediction = mvote.DSDWeightedMV(ppfDSD, ppbLabel, pnFoldIndex, pnRD, options.neighbor)
elif options.mode == 3:
    (masterPredictionMatrix, mapProtNamesToMasterIdx, masterLabelMatrix) = mvote.DSDWeightedMVIterativeSetup(numLabels, options.rdindex, options.completeProteinList, ppbLabel)
    masterPredictionMatrix, masterLabelMatrix = mvote.DSDWeightedMVIterative(options.dsdfile, masterLabelMatrix, masterPredictionMatrix, options.neighbor, options.completeProteinList, mapProtNamesToMasterIdx)

#### Phase 3: Write Output

if options.mode == 0:
    options.outfile = 'OrdMV' + options.outfile
elif options.mode == 1:
    options.outfile = 'DSDUnweight' + options.outfile
elif options.mode == 2:
    options.outfile = 'DSDWeighted' + options.outfile
elif options.mode == 3:
    options.outfile = 'DSDWeightedIterative' + options.outfile
mvote.writeOutput(masterPredictionMatrix, options.outfile)
