#!/usr/sup/bin/python
'''
calcDSD.py -- This module parse calculates DSD given the adjacency matrix
    and output accoding to options

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

'''

import numpy as np
import re


def calculator(adjacency, nRW, transition, quiet=False):
    """
    adjacency - adjacency matrix represented as a numpy array
                assumes graph is fully connected.

    nRW - the length of random walks used to calculate DSD
          if nRW = -1, then calculate 

    quiet - Whether to print verbose messages.

    returns DSD matrix represented as a numpy array
    """
    #### p for transition matrix
    n = np.size(adjacency[0])
    degree = np.zeros((n, 1))
    for j in xrange(0, n):
        degree[j] = sum(adjacency[j])
    if transition != None:
        p = transition
    else:
        p = np.zeros((n, n))
        # print 'there are {0:.0f} nodes'.format(n)
        for j in xrange(0, n):
            #degree[j] = sum(adjacency[j])
            for i in xrange(0, n):
                if degree[j] != 0:
                    p[j] = adjacency[j]/degree[j]
    

    if nRW >= 0:
        #### c for visit count matrix
        #### for example, c(2,3) is the number of times
        ####     that node 3 is visited via random walks
        ####     starting from node 2
        c = np.eye(n)
        for rw in xrange(0, nRW):
            #print c
            c = np.dot(c, p) + np.eye(n)
            #### this is the c matrix for random walks with
            #### length rw, i.e. if length is 0, then c is
            #### simply the identiy matrix, visiting itself
            #### if length is 1, then c is eye() plus the
            #### the one step transition matrix
    else:
        c = np.eye(n)
        c = c - p
        pi = (degree.conj().T)/sum(degree)
        c = c + np.tile(pi, (n, 1))
        c = np.linalg.inv(c)

    DSD = np.zeros((n, n))
    for i in xrange(0, n):
        for j in xrange(i+1, n):
            if degree[i] and degree[j]:
                DSD[i, j] = np.linalg.norm((c[i, :]-c[j, :]), ord=1)
                DSD[j, i] = DSD[i, j]
            else:
                DSD[i, j] = -1
                DSD[j, i] = -1
        if(not quiet) and ((i % 100 == 0) or (i == n-1)):
            print('    finish calculating DSD for %d/%d nodes' % (i+1, n))

    return DSD


def writeoutMatrix(DSD, names, ofile):
    """
    write DSD matrix into a tab delimited csv file

    DSD -- the DSD matrix, an numpy matrix

    names -- the node IDs, a dict

    ofile -- the output file object

    """
    n = np.size(DSD[0])
    count = 0
    temp = "\t"
    sNames = names.keys()
    while(count < n-1):
        temp = temp + sNames[count] + '\t'
        count = count + 1
    temp = temp + sNames[n-1] + '\n'
    ofile.write(temp)
    for i in xrange(0, n):
        temp = sNames[i]
        for j in xrange(0, n):
            if DSD[i, j] < 0:
                temp = temp + '\tNA'
            else:
                temp = temp + '\t%.4f' % DSD[i, j]
        ofile.write(temp + '\n')
    return True


def writeoutList(DSD, names, infile, ofile):
    """
    write DSD into a file, where each line is an interaction from
    input file, followed by the DSD values

    DSD -- the DSD matrix, an numpy matrix

    names -- the node IDs, a dict

    infile -- the name of input file

    ofile -- the output file object

    """
    splitpattern = re.compile('[\t ;,]+')
    validpattern = re.compile('^[\w _\-.,\t\':;"]+$')
    ifile = open(infile, 'r')
    sNames = names.keys()
    for temp in ifile:
        temp = temp.strip('\t \n\r')
        if temp == "" or re.search(validpattern, temp) is None:
            continue
        else:
            allwords = re.split(splitpattern, temp)
            if(allwords[0] not in names) or (allwords[1] not in names):
                temp = allwords[0] + '\t' + allwords[1] + '\tNotConnected\n'
            else:
                i = names[allwords[0]]
                j = names[allwords[1]]
                temp = sNames[i] + '\t' + sNames[j] + '\t%.4f\n' % DSD[i, j]
            ofile.write(temp)
    ifile.close()
    return True
 #   print adjacency
 #   print p
 #   print degree


def writeoutToplist(DSD, names, ofile, nTop):
    """
    write DSD into a file, where for each node, we put at the first
    column, then followed by nTop nodes with lowest DSD, in increasing
    order

    DSD -- the DSD matrix, an numpy matrix

    names -- the node IDs, a dict

    ofile -- the output file object

    nTop -- the number of nodes with lowest DSD

    """
    n = np.size(DSD[0])
    sNames = names.keys()
    for i in xrange(0, n):
        temp = sNames[i]
        TopKIndex = findTopk(DSD[i], nTop, i)
        for j in xrange(0, nTop):
            if TopKIndex[j] >= 0:
                temp = temp + '\t' + sNames[TopKIndex[j]]
                temp = temp + ('(%.4f)' % DSD[i, TopKIndex[j]])
            else:
                break
        ofile.write(temp + '\n')
    return True


def findTopk(values, k, a):
    '''
    it finds the k smallest items in list "values" excluding the "a"th item

    values -- a list of numerical values

    k -- the number of smallest items

    a -- the item to be excluded

    returns a list of indices with smallest values, increasing order

    '''
    n = np.size(values)
    indicator = [True]*n
    indicator[a] = False
    index = [-1]*k
    m = max(values)
    for i in xrange(0, n):
        if values[i] < 0:
            indicator[i] = False
    for i in xrange(0, k):
        temp = m
        for j in xrange(0, n):
            if(values[j] <= temp) and indicator[j]:
                temp = values[j]
                index[i] = j
        if index[i] >= 0:
            indicator[int(index[i])] = False
        else:
            break
    return index
