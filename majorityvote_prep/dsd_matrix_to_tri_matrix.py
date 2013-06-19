#!/usr/bin/env python2.7

""" Takes in a DSD matrix and an output file and outputs the upper right triangle of that matrix
    above the upperleft-to-lowerright diagonal.
    The DSD matrix must be square, all DSD values must be separated by tabs
    in the DSD matrix, and are separated by tabs in the triangular matrix """
""" Inbar Fried 6/18/2013 """
""" to run: python this_file.py dsd_matrix output_file """

import sys

if len (sys.argv) >= 4:
    dsd_matrixf = open(sys.argv[1], 'r')
    tri_matrixf = open(sys.argv[2], 'w')
    order_list = open(sys.argv[3], 'r')

    lines = dsd_matrixf.readlines()

    i = j = n = 0

    """ count line number in matrix """
    with open(sys.argv[1], 'r') as f:
        for line in f:
            n += 1

    for i in range(1, n):
        words = lines[i].split('\t')
        for j in range(i+1, n):
            if j == n-1:
                tri_matrixf.write(words[j])
            else:
                tri_matrixf.write(words[j] + '\t')

    tri_matrixf.close()
    dsd_matrixf.close()

else:
    print "bad input"
