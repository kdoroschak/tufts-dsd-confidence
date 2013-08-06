#!/usr/bin/env python2.7
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
    for iAnnoPro in xrange(0, m1):
        iPro = pnRD[iAnnoPro]
        prediction[iAnnoPro, 0] = iPro
        iFold = pnFoldIndex[iPro]
        prelist = np.zeros((1, m))
        for jPro in xrange(0, N):
            if ppbAdj[iPro, jPro] and iFold != pnFoldIndex[jPro] and not ppbLabel[jPro,0]:
                for iLabel in xrange(0, m):
                    prelist[0, iLabel] += ppbLabel[jPro, iLabel+1]
        sortedlabel = np.argsort(prelist[0,:])
        for iLabel in xrange(0, m):
            prediction[iAnnoPro, iLabel*2+1] = sortedlabel[m-1-iLabel]
            prediction[iAnnoPro, iLabel*2+2] = prelist[0,sortedlabel[m-1-iLabel]]
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
    for iAnnoPro in xrange(0, m1):
        iPro = pnRD[iAnnoPro]
        prediction[iAnnoPro, 0] = iPro
        iFold = pnFoldIndex[iPro]
        prelist = np.zeros((1, m))
        sortedDSD = np.argsort(ppfDSD[iPro,:])
        j = 1
        count = 0
        while ((j < N) and (count < top)):
            jPro = sortedDSD[0,j]
            if iFold != pnFoldIndex[jPro] and not ppbLabel[jPro,0]:
                count += 1
                for iLabel in xrange(0, m):
                    prelist[0, iLabel] += ppbLabel[jPro, iLabel+1]
            j += 1
        sortedlabel = np.argsort(prelist[0,:])
        for iLabel in xrange(0, m):
            prediction[iAnnoPro, iLabel*2+1] = sortedlabel[m-1-iLabel]
            prediction[iAnnoPro, iLabel*2+2] = prelist[0,sortedlabel[m-1-iLabel]]
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
    for iAnnoPro in xrange(0, m1):
        iPro = pnRD[iAnnoPro]
        prediction[iAnnoPro, 0] = iPro
        iFold = pnFoldIndex[iPro]
        prelist = np.zeros((1, m))
        sortedDSD = np.argsort(ppfDSD[iPro,:])
        j = 1
        count = 0
        while ((j < N) and (count < top)):
            jPro = sortedDSD[0,j]
            if iFold != pnFoldIndex[jPro] and not ppbLabel[jPro,0]:
                count += 1
                for iLabel in xrange(0, m):
                    prelist[0, iLabel] += (ppbLabel[jPro, iLabel+1]/ppfDSD[iPro, jPro])
            j += 1
        sortedlabel = np.argsort(prelist[0,:])
        for iLabel in xrange(0, m):
            prediction[iAnnoPro, iLabel*2+1] = sortedlabel[m-1-iLabel]
            prediction[iAnnoPro, iLabel*2+2] = prelist[0,sortedlabel[m-1-iLabel]]
    return prediction


#added by Tony Cannistra, June 2013 
def SeqDSDWeightedMV(ppfDSD, ppfSeq, ppbLabel, pnFoldIndex, pnRD, top, DSDScale, SeqScale, SeqNorm):
    '''
    Weighted DSD Majority Voting
    input: ppfDSD -- DSD matrix
           ppfSeq -- sequence matrix
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
    print "DSDScale: ", DSDScale, ", SeqScale: ", SeqScale, ", SeqNorm: ", SeqNorm
    prediction = np.zeros((m1, 2*m+1))
    for iAnnoPro in xrange(0, m1):
        iPro = pnRD[iAnnoPro]
        prediction[iAnnoPro, 0] = iPro
        iFold = pnFoldIndex[iPro]
        prelist = np.zeros((1, m))
        sortedDSD = np.argsort(ppfDSD[iPro,:])
        j = 1
        count = 0
        while ((j < N) and (count < top)):
            jPro = sortedDSD[0,j]
            if iFold != pnFoldIndex[jPro] and not ppbLabel[jPro,0]:
                count += 1
                ### BEWARE OF DIVIDE BY ZERO

                seqWeight = 1
                if(ppfSeq[iPro, jPro] != 0):
                  #NON RECIPROCAL
                  seqWeight = (SeqNorm)*ppfSeq[iPro, jPro]
                else:
                  print "warning! 0 sequence score."

                weight = (DSDScale/ppfDSD[iPro, jPro]) + (SeqScale*seqWeight)
                for iLabel in xrange(0, m):
                    prelist[0, iLabel] += (ppbLabel[jPro, iLabel+1]*weight)
            j += 1
        sortedlabel = np.argsort(prelist[0,:])
        for iLabel in xrange(0, m):
            prediction[iAnnoPro, iLabel*2+1] = sortedlabel[m-1-iLabel]
            prediction[iAnnoPro, iLabel*2+2] = prelist[0,sortedlabel[m-1-iLabel]]
    return prediction


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
    for iAnnoPro in xrange(0, m1):
        sout.write(('%d\t' % (int(1+prediction[iAnnoPro, 0]))))
        for iLabel in xrange(0, m):
            sout.write('%d\t%f\t' % (int(1+prediction[iAnnoPro, iLabel*2+1]), prediction[iAnnoPro, iLabel*2+2]))
        sout.write('\n')
    sout.close()
