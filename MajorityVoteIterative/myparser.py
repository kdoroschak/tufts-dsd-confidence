#!/usr/bin/env python2.7
"""
Created on Tue May 21 08:47:39 2013

@author: mcao01

myparser:
    parsePPI, reads PPI from input file and returns adjacency matrix

    parseLabel, reads labels from input file and returns label matrix

    parseDSD, reads DSD from input file and returns DSD matrix

    parseRD, reads random indeces from input file and returns the array

"""

import re
import numpy as np
import collections
import sys


def parsePPI(filename):
    '''
    it parses PPIs from filename
    '''
    print 'Parsing PPI file...'
    infile = open(filename, 'r')
    index = 0
    for line in infile:
        index += 1
    N = index
    infile.close()
    ppbAdj = np.zeros((N, N))
    print 'there are ', N, ' nodes from PPI file\n'
    
    infile = open(filename, 'r')
    index = 0
    for line in infile:
        for j in xrange(index+1, N):
            if line[j] == '1':
                ppbAdj[index, j] = 1
                ppbAdj[j, index] = 1
        index += 1
    degree = sum(ppbAdj)
####    print degree[0], degree[1], degree[N-1], degree[N-2]
    infile.close()
    return ppbAdj


def parseLabel(filename):
    '''
    it parses labels from filename
    '''
    print 'Parsing label file...'
    infile = open(filename, 'r')    
    index = 0
    m = 1
    for line in infile:
        line = line.strip('\r\t\n ')
        index += 1
        if m < len(line) - 1:
            m = len(line) - 1
    N = index
    print 'there are ', N, ' nodes and ', m, ' labels from label file'
    infile.close()
    ppbLabel = np.zeros((N, m+1))
    infile = open(filename, 'r')
    index = 0
    for line in infile:
        if line[0] == '1':
            ppbLabel[index, 0] = 1
        else:
            for i in xrange(1, m+1):
                if line[i] == '1':
                    ppbLabel[index, i] = 1
        index += 1
    infile.close()
    
    m1 = N - sum(ppbLabel[:,0])
    m2 = sum(sum(ppbLabel[:,1:(m+1)]))
    print 'there are ', m1, ' labeled nodes and ', m2, ' labels in total\n'
    return ppbLabel
    
    
def parseRDIndex(filename, ppbLabel):
    '''
    it parses random indeces from filename and verify their labels
    '''
    print 'Parsing random indices of labeled nodes'
    infile = open(filename, 'r')    
    index = 0
    n1 = len(ppbLabel[:,0]) - sum(ppbLabel[:,0])
    pnRD = np.zeros((n1, 1), dtype=np.int)

    for line in infile:
        line = line.strip('\r\n\t ')
        if line == '':
            continue
        temp = int(line)
        if index == n1:
            print >> sys.stderr, 'there are more entries in Random index'
            print >> sys.stderr, 'file than there should be(not enough labels)'
            exit(0)
        if ppbLabel[temp-1, 0] == 1:
            print >> sys.stderr, 'Random Index File Not Consistent: Possibly'
            print >> sys.stderr, 'there are indeces of which node is not'
            print >> sys.stderr, 'labeled (index starting from 1 in the file)'
            exit(0)
        pnRD[index] = temp-1
        index += 1
    infile.close()
    print 'there are ', n1, ' labeled nodes from Random Index, confirmed\n'
    return pnRD


def parseDSD(filename):
    '''
    it parses DSD from input file
    the DSD file must contain exactly N lines, where N is the number
    of nodes, leaving the last data line empty
    '''
    print 'Parsing DSD file'
    infile = open(filename, 'r')
    index = 0
    for line in infile:
        index += 1
    N = index
    print 'there are ', N, ' nodes from DSD files'
    infile.close()

    infile = open(filename, 'r')
    ppfDSD = np.zeros((N, N))
    index = 0
    splitpattern = re.compile('[\r\t\n ]')
    for line in infile:
        line = line.strip('\r\n\t ')
        allnum = re.split(splitpattern, line)
        for j in xrange(0, N-index-1):
            ppfDSD[index, index+j+1] = float(allnum[j])
            ppfDSD[index+j+1, index] = ppfDSD[index, index+j+1]
        index += 1
    infile.close()
    print 'the average DSD is ', sum(sum(ppfDSD))/float(N*N-N), '\n'
    return ppfDSD  


def parsePrediction(filename):
    '''
    parse the prediction file, writen by mvote.writeOutput()
    return the prediction matrix
    
    '''
    ifile = open(filename, 'r')
    m1 = 0
    for line in ifile:
        m1 += 1
    allwords = line.split('\t')
    m = (len(allwords) - 1)/2
    print 'Parsing the prediction file,'
    print 'there are ', m1, ' nodes predicted with ', m, ' labels\n'
    ifile.close()

    prediction = np.zeros((m1, 2*m+1))

    ifile = open(filename, 'r')
    index = 0
    for line in ifile:
        allwords = line.split('\t')
        prediction[index, 0] = int(allwords[0]) - 1
        for iLabel in xrange(0, m):
            prediction[index, 2*iLabel+1] = int(allwords[2*iLabel+1]) - 1
            prediction[index, 2*iLabel+2] = float(allwords[2*iLabel+2])
        index += 1
    ifile.close()
    return prediction


def GetFoldIndex(pnRD, N, kfold):
    pnFoldIndex = np.zeros((N, 1), dtype=np.int)
    m1 = len(pnRD)
    foldsize = int(m1/kfold)
    for i in xrange(1, kfold+1):
        for j in xrange(0, foldsize):
            pnFoldIndex[pnRD[(i-1)*foldsize+j]] = i
    return pnFoldIndex