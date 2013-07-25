#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:56:13 2013

@author: mcao01

mvote, used to do k-fold majority voting

    ordinaryMV: takes PPI, label, Fold Index as input
                output prediction matrix, the first column is the index
                of nodes, starting from 1; the second column is the index
                of label with top vote, starting from 1, followed by the
                votes as the third column. 
"""


import re
import numpy as np
import collections
import sys


def ordinaryMV(ppbAdj, ppbLabel, pnFoldIndex, pnRD):
    '''
    Ordinary Majority Voting
    input: ppbAdj -- adjacency matrix
           ppbLabel -- label matrix
           pnFoldIndex -- the ith entry is the fold number which the ith
                          node belongs to
           pnRD -- the set of random indices

    output: prediction matrix, first row: the index of node
                              second row: the index of label with top votes
                              third row: the votes of the label for the 2nd row
                              ...
                              ...
                              last but 1 row: index of label with least votes
                              last row: votes of the previous label
    all indeces start from 1
    '''
    N = len(ppbAdj[:,1]) ### number of nodes
    m = len(ppbLabel[1,:]) -1 ### number of labels
    m1 = sum(pnFoldIndex[:,0] != 0) ### number of labeled nodes
    kfold = max(pnFoldIndex[:,0])
    print 'there are ', N, ' nodes and ', m, ' labels; in this ', kfold
    print '-fold cross validation, there are ', m1, ' labeled nodes.'
    prediction = np.zeros((m1, 2*m+1))
    for annoPro_i in xrange(0, m1):
        pro_i = pnRD[annoPro_i]
        prediction[annoPro_i, 0] = pro_i
        fold_i = pnFoldIndex[pro_i]
        prelist = np.zeros((1, m))
        for pro_j in xrange(0, N):
            if ppbAdj[pro_i, pro_j] and fold_i != pnFoldIndex[pro_j] and not ppbLabel[pro_j,0]:
                for label_i in xrange(0, m):
                    prelist[0, label_i] += ppbLabel[pro_j, label_i+1]
        sortedlabel = np.argsort(prelist[0,:])
        for label_i in xrange(0, m):
            prediction[annoPro_i, label_i*2+1] = sortedlabel[m-1-label_i]
            prediction[annoPro_i, label_i*2+2] = prelist[0,sortedlabel[m-1-label_i]]
    return prediction


def DSDUnweightMV(ppfDSD, ppbLabel, pnFoldIndex, pnRD, top):
    '''
    Unweighted DSD Majority Voting
    input: ppfDSD -- DSD matrix
           ppbLabel -- label matrix
           pnFoldIndex -- the ith entry is the fold number which the ith
                          node belongs to
           pnRD -- the set of random indices
           top -- the number of nodes with lowest DSD used for voting
    output: prediction matrix, first row: the index of node
                              second row: the index of label with top votes
                              third row: the votes of the label for the 2nd row
                              ...
                              ...
                              last but 1 row: index of label with least votes
                              last row: votes of the previous label
    all indeces start from 1
    '''
    N = len(ppfDSD[:,1]) ### number of nodes
    m = len(ppbLabel[1,:]) -1 ### number of labels
    m1 = sum(pnFoldIndex[:,0] != 0) ### number of labeled nodes
    kfold = max(pnFoldIndex[:,0])
    print 'there are ', N, ' nodes and ', m, ' labels; in this ', kfold
    print '-fold cross validation, there are ', m1, ' labeled nodes;'
    print 'take ', top, ' nodes with lowest DSD for unweighted voting'
    prediction = np.zeros((m1, 2*m+1))
    top = int(top)
    for annoPro_i in xrange(0, m1):
        pro_i = pnRD[annoPro_i]
        prediction[annoPro_i, 0] = pro_i
        fold_i = pnFoldIndex[pro_i]
        prelist = np.zeros((1, m))
        sortedDSD = np.argsort(ppfDSD[pro_i,:])
        j = 1
        count = 0
        while ((j < N) and (count < top)):
            pro_j = sortedDSD[0,j]
            if fold_i != pnFoldIndex[pro_j] and not ppbLabel[pro_j,0]:
                count += 1
                for label_i in xrange(0, m):
                    prelist[0, label_i] += ppbLabel[pro_j, label_i+1]
            j += 1
        sortedlabel = np.argsort(prelist[0,:])
        for label_i in xrange(0, m):
            prediction[annoPro_i, label_i*2+1] = sortedlabel[m-1-label_i]
            prediction[annoPro_i, label_i*2+2] = prelist[0,sortedlabel[m-1-label_i]]
    return prediction


def DSDWeightedMV(ppfDSD, ppbLabel, pnFoldIndex, pnRD, top):
    '''
    Weighted DSD Majority Voting
    input: ppfDSD -- DSD matrix
           ppbLabel -- label matrix
           pnFoldIndex -- the ith entry is the fold number which the ith
                          node belongs to
           pnRD -- the set of random indices
           top -- the number of nodes with lowest DSD used for voting
    output: prediction matrix, first row: the index of node
                              second row: the index of label with top votes
                              third row: the votes of the label for the 2nd row
                              ...
                              ...
                              last but 1 row: index of label with least votes
                              last row: votes of the previous label
    all indeces start from 1
    '''
    N = len(ppfDSD[:,1]) ### number of nodes
    m = len(ppbLabel[1,:]) -1 ### number of labels
    m1 = sum(pnFoldIndex[:,0] != 0) ### number of labeled nodes
    kfold = max(pnFoldIndex[:,0])
    print 'there are ', N, ' nodes and ', m, ' labels; in this ', kfold
    print '-fold cross validation, there are ', m1, ' labeled nodes;'
    print 'take ', top, ' nodes with lowest DSD for weighted voting'
    prediction = np.zeros((m1, 2*m+1))
    for annoPro_i in xrange(0, m1):
        pro_i = pnRD[annoPro_i]
        prediction[annoPro_i, 0] = pro_i
        fold_i = pnFoldIndex[pro_i]
        prelist = np.zeros((1, m))
        sortedDSD = np.argsort(ppfDSD[pro_i,:])
        j = 1
        count = 0
        while ((j < N) and (count < top)):
            pro_j = sortedDSD[0,j]
            if fold_i != pnFoldIndex[pro_j] and not ppbLabel[pro_j,0]:
                count += 1
                for label_i in xrange(0, m):
                    prelist[0, label_i] += (ppbLabel[pro_j, label_i+1]/ppfDSD[pro_i, pro_j])
            j += 1
        sortedlabel = np.argsort(prelist[0,:])
        #print prelist
        #print sortedlabel
        for label_i in xrange(0, m):
            prediction[annoPro_i, label_i*2+1] = sortedlabel[m-1-label_i]
            prediction[annoPro_i, label_i*2+2] = prelist[0,sortedlabel[m-1-label_i]]
    return prediction
    
    
def DSDWeightedMVIterativeSetup(numLabels, randomIdxFile, completeProteinListFile):
    # Read in the full list of proteins (labels for master matrix)
    # Create map from protein name to index in prediction matrix
    mapProtNamesToMasterIdx = {}
    with open(completeProteinListFile, 'r') as proteinFile:
        for i,protein in enumerate(proteinFile):
            mapProtNamesToMasterIdx[protein] = i
    numProteins = i + 1

    # Read in random index file (global)
    with open(randomIdxFile, 'r') as randomIdxFile:
        randomIndices = []
        for i,line in enumerate(randomIdxFile):
            randomIndices.append(line.strip())
        numRandomIndices = len(randomIndices)
        randomIdxSet = set(randomIndices[0:numRandomIndices/2])
        
    # Read in annotation file (global)
    # Set up an annotation list with half of all labeled nodes "covered"
    #numLabels = 0
    #trainingAnns = []
    #with open(annFile, 'r') as annFile:
        #for i,line in enumerate(annFile):
            #line = line.strip()
            #if i in randomIdxSet:
                #numLabels = len(line)
                #if numLabels <= 1:
                    #print "numLabels <= 1, reading in annotation file"
                #trainingAnns.append('1')
            #else:
                #trainingAnns.append(line)
    #print trainingAnns

    # Create empty prediction matrix (numAllProteins X 2*numLabels+1)    
    masterPredictionMatrix = np.zeros((numProteins, 2*numLabels+1))
    
    return (masterPredictionMatrix, mapProtNamesToMasterIdx)
    
    
### BEGINNING OF ITERATIVE 
def DSDWeightedMVIterative(ppfDSD, ppbLabel, pnRD, top, localProteinList):
    '''
    Weighted DSD Majority Voting
    input: ppfDSD -- DSD matrix
           ppbLabel -- label matrix
           pnFoldIndex -- the ith entry is the fold number which the ith
                          node belongs to
           pnRD -- the set of random indices
           top -- the number of nodes with lowest DSD used for voting
    output: prediction matrix, first row: the index of node
                       second row: the index of label with top votes
                       third row: the votes of the label for the 2nd row
                              ...
                              ...
                         last but 1 row: index of label with least votes
                         last row: votes of the previous label
    all indeces start from 1
    '''
    
    # ITERATIONS:
    # Read in DSD for single iteration
    # Read in the ordered protein names for the DSD
    # Read in the internal annotation file (produced in setup, added to in each iteration) and parse to size and order of ordered labels
    # Read in top (number of nodes used for voting)
    # Pass in prediction matrix (to be appended to for results)
    
    # Run modified version of MV
    #     remove k-fold portion
    #     map results to indices of big prediction matrix
    
    # "Return":
    #     Annotation file
    #     Prediction matrix

    
    
   

##### WRITE
def writeOutput(prediction, filename):
    '''
    input:  prediction matrix, first row: the index of node
                              second row: the index of label with top votes
                              third row: the votes of the label for the 2nd row
                              ...
                              ...
                              last but 1 row: index of label with least votes
                              last row: votes of the previous label
            filename: output file name
    output: write a new file with prediction in it, tab-delimimited
    all indeces start from 1
    '''
    print 'writing prediction into file ', filename, '\n'
    sout = open(filename, 'w')
    m1 = len(prediction[:,0])
    m = len(prediction[0, :]-1)/2
    for annoPro_i in xrange(0, m1):
        sout.write(('%d\t' % (int(1+prediction[annoPro_i, 0]))))
        for label_i in xrange(0, m):
            sout.write('%d\t%f\t' % (int(1+prediction[annoPro_i, label_i*2+1]), prediction[annoPro_i, label_i*2+2]))
        sout.write('\n')
    sout.close()
