#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-


import sys
import numpy as np
import re
import collections
import PPIparser
import argparse

def main():

    # Set up argument parser
    temp = "Extracts adjacency matrix and ordered"
    temp += " list of proteins from PPI file."
    parser = argparse.ArgumentParser(description=temp)

    parser.add_argument("-p", "--ppi", help="PPI input file.",
                        default=None)
    parser.add_argument("-m", "--mat", help="Adjacency matrix output file.",
                        default=None)
    parser.add_argument("-o", "--oprot", help="Ordered list of proteins output file.",
                        default=None)

    # Parse options
    options = parser.parse_args()
    ppi = options.ppi
    mat = options.mat
    oprot = options.oprot


    (ppbAdj, names) = PPIparser.GetAdj(ppi, 0)
    allnames = names.keys()

    N = len(names)

    fout = open(mat, 'w')
    fnameout = open(oprot, 'w')

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

if __name__ == "__main__":
    main()
