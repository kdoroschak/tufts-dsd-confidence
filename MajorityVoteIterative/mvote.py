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
import myparser


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


def DSDWeightedMVIterativeSetup(numLabels, randomIdxFile, completeProteinListFile, masterLabelMatrix):
    # Read in the full list of proteins (labels for master matrix)
    # Create map from protein name to index in prediction matrix
    mapProtNamesToMasterIdx = {}
    with open(completeProteinListFile, 'r') as proteinFile:
        for i,protein in enumerate(proteinFile):
            mapProtNamesToMasterIdx[protein.strip()] = i
    numProteins = i + 1

    # Read in random index file (global)
    with open(randomIdxFile, 'r') as randomIdxFile:
        randomIndices = []
        for i,line in enumerate(randomIdxFile):
            randomIndices.append(int(line.strip()))
        numRandomIndices = len(randomIndices)
        indicesToCover = randomIndices[:(numRandomIndices/2)]

    # Read in annotation file (global)
    # Set up an annotation list with half of all labeled nodes "covered"
    for index in indicesToCover:
        for i, item in enumerate(masterLabelMatrix[index-1]):
            masterLabelMatrix[index-1, i] = 0
        masterLabelMatrix[index-1, 0] = 1

    # Create empty prediction matrix (numAllProteins X 2*numLabels+1)    
    masterPredictionMatrix = np.zeros((numProteins, 2*numLabels+1))

    # Set first column to protein index (in order)
    for (protein, index) in mapProtNamesToMasterIdx.items():
        masterPredictionMatrix[index, 0] = index

    return (masterPredictionMatrix, mapProtNamesToMasterIdx, masterLabelMatrix)


### BEGINNING OF ITERATIVE 
def DSDWeightedMVIterative(dsdFile, masterLabelMatrix, masterPredictionMatrix,
        desiredNumClosest, toPredictFile,
        proteinToMasterIndex, localKeyFile):

    '''
    Weighted DSD Majority Voting
    input: DSDfile -- file containing DSD triangular matrix
           labelMatrix -- label matrix
           predictionMatrix -- matrix of predictions. see "output"
           randomIndexSet -- the set of random indices
           desiredNumClosest -- the number of nodes with lowest DSD used for voting
           localProteinList -- local list of proteins. file format.
           proteinToMasterIndex -- dictionary mapping protein names to master indices.
    output: prediction matrix, first row: the index of node
                       second row: the index of label with top votes
                       third row: the votes of the label for the 2nd row
                              ...
                              ...
                         last but 1 row: index of label with least votes
                         last row: votes of the previous label
    all indices start from 1
    '''

    # NAMING: "local" refers to all proteins in the current dsd
    #         "toPredict" refers to all proteins that we need to generate predictions for (<= "local")
    #         "master" refers to all proteins in the complete protein list

    # ITERATIONS:
    # Read in DSD for single iteration
    #   Send to parser
    # Read in the ordered protein names for the DSD
    # Pass in the internal label matrix (produced in setup,
    #   added to in each iteration) and parse to size
    #   and order of ordered labels
    # Pass in top (number of nodes used for voting)
    # Pass in prediction matrix (to be appended to for results)

    # Run modified version of MV (see below)

    # "Return":
    #     Annotation file
    #     Prediction matrix

    # Set up initial values
    numLabels = masterLabelMatrix.shape[1] - 1

    # Parse DSD file
    dsdMatrix = myparser.parseDSD(dsdFile)

    # Read in proteins from dsd key
    proteinNameToLocalIndex = {}
    localProteins = [] # stores all local protein names, not just to predict
    with open(localKeyFile, 'r') as localKeyFile: # loops through all local proteins
        for i, protein in enumerate(localKeyFile):
            # Make list of proteins in order of DSD file
            proteinNameToLocalIndex[protein.strip()] = i
            localProteins.append(protein.strip())
    numLocalProteins = len(proteinNameToLocalIndex.keys())

    # Read in proteins that need to be predicted
    localIndicesToPredict = [] # Stores local (dsd) indices of proteins to predict
    with open(toPredictFile, 'r') as toPredictFile: # only loops through proteins to be predicted
        for line in toPredictFile:
            protein = line.strip()
            localIndex = proteinNameToLocalIndex[protein]
            localIndicesToPredict.append(localIndex)
    numToPredict = len(localIndicesToPredict)

    # Calculate prediction values for each unlabeled protein
    for localProteinIndex in localIndicesToPredict:
        predictionList = np.zeros(numLabels)

        indicesOfSortedDSD = np.argsort(dsdMatrix[localProteinIndex,:])

        # Counter for the number of closest DSD values extracted so far
        #   compared with desiredNumClosest, max # of DSD values for voting
        numClosestChosen = 0

        # Convert to master index
        proteinName = localProteins[localProteinIndex]
        masterProteinIndex = proteinToMasterIndex[proteinName]

        # Go through DSD values from smallest to largest
        #   Calculate prediction values for closest labeled nodes
        for nextClosestDSDIndex in indicesOfSortedDSD[1:]:

            # Look up master index of node with next closest DSD value
            currentDSDProteinName = localProteins[nextClosestDSDIndex]
            masterIndexOfCurrentDSD = proteinToMasterIndex[currentDSDProteinName]
            #print masterIndexOfCurrentDSD

            if numClosestChosen >= desiredNumClosest:
                break

            # Only calculate if we're looking at a labeled protein
            if not masterLabelMatrix[masterIndexOfCurrentDSD, 0]:
                #print "   ", masterIndexOfCurrentDSD, masterLabelMatrix[masterIndexOfCurrentDSD]
                numClosestChosen += 1
                for labeli in xrange(0, numLabels):
                    predictionVal = (masterLabelMatrix[masterIndexOfCurrentDSD, labeli + 1] / dsdMatrix[localProteinIndex, nextClosestDSDIndex])
                    predictionList[labeli] += predictionVal

        # Get indices of the sorted prediction values from high to low confidence
        indicesOfSortedPredictionValues = np.argsort(predictionList)[::-1]

        # Update prediction matrix
        for labeli in xrange(0, numLabels):
            # Populate prediction matrix with ordered prediction values and labels
            labelNum = indicesOfSortedPredictionValues[labeli]

            masterPredictionMatrix[masterProteinIndex, labeli * 2 + 1] = labelNum
            masterPredictionMatrix[masterProteinIndex, labeli * 2 + 2] = predictionList[labelNum]


    # store in label matrix
    for localProteinIndex in localIndicesToPredict:
        for labeli in xrange(0, numLabels):
            # Extract first prediction to put in label matrix
            firstPrediction = masterPredictionMatrix[masterProteinIndex, 1]

            # Make sure there is a prediction, then store in label matrix
            if firstPrediction != 0:
                masterLabelMatrix[masterProteinIndex, firstPrediction + 1] = 1
                masterLabelMatrix[masterProteinIndex, 0] = 0

    return masterPredictionMatrix, masterLabelMatrix

    # in set of random indices, for each labeled node index
    #   insert label into prediction matrix (first column) **Now done in setup
    #   return indices of DSD sorted by DSD value (DSD does not change)
    #   iterate through nodes (N) in order of sorted indices *explain more?
    #       if the node is labeled, and we haven't already gotten enough labeled nodes
    #           extract labeled node (l)
    #           compute prediction value (label/DSD)
    #           add list of predictions for the current node (N)
    #       get list of indices of predictions, sorted by prediction value
    #           (low -> high in MF's code)
    #       for each prediction in order
    #           find use index to find prediction
    #           output prediction to matrix along with prediction value (confidence?)
    # return prediction
'''
    N = len(ppfDSD[:,1]) ### number of nodes
    m = len(ppbLabel[1,:]) -1 ### number of labels
    m1 = sum(pnFoldIndex[:,0] != 0) ### number of labeled nodes
    prediction = np.zeros((m1, 2*m+1))
    for annoPro_i in xrange(0, m1):
        ***pro_i = pnRD[annoPro_i]
        ***prediction[annoPro_i, 0] = pro_i
        prelist = np.zeros((1, m))
        sortedDSD = np.argsort(ppfDSD[pro_i,:])
        j = 1
        count = 0
        while ((j < N) and (count < top)):
            pro_j = sortedDSD[0,j]
            if not ppbLabel[pro_j,0]:
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
'''

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
