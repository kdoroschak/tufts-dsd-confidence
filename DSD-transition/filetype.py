#!/usr/sup/bin/python
'''
filetype.py -- detect the type of input file for PPIconvert

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

Descriptions:
  0) suffixDetector: detect file type by suffix

  type 1    line stored PPIs: contain exactly m by c entries, the
                              first two columns are node idendifications
  e.g.:     line1: Pro1 Pro2 2.23
            line2: Pro3 Pro1 2.3
            line3: Pro1 Pro4 4.2
            line4: Pro5 Pro2 1.2

  type 1 filename must have suffix .list

  type 2    (n+1) by (n+1) matrices: the first colulmn and the first row
                                     contains node ID assume any entry
                                     in the first row must also appear
                                     in the first column

  e.g.:     line1:      Pro1 Pro2 Pro3
            line2: Pro2 2    0    1.2
            line3: Pro1 0    4.2  0.1
            line4: Pro3 3    1.2  0

  type 2 filename must have suffix .tab or .csv


  1) ft: used to support other 3 file formats by detecting file formats

  1.1) input files consider any line starting with letters or digits
        as a valid line (skipping the first spaces/tabs)

  1.2) and the module uses two continuous valid line to detect file type

  1.3) consider each valid line as tab/comma/space delimited table

  filetype  file_description
  0         non-consistent files
  e.g.:     line1: abcd 123 21
            line2: 123 abcd
            line3: ->12

  1         line stored PPIs: contain exactly m by c entries, the
                              first two columns are node idendifications
  e.g.:     line1: Pro1 Pro2 2.23
            line2: Pro3 Pro1 2.3
            line3: Pro1 Pro4 4.2
            line4: Pro5 Pro2 1.2

  2         (n+1) by (n+1) matrices: the first colulmn and the first row
                                     contains node ID assume any entry
                                     in the first row must also appear
                                     in the first column
  e.g.:     line1:      Pro1 Pro2 Pro3
            line2: Pro2 2    0    1.2
            line3: Pro1 0    4.2  0.1
            line4: Pro3 3    1.2  0

  3         (n) by (n+1) matrices: the first colulmn contains node ID,
                                   assume columns and rows match
  e.g.:     line1: Pro1 0    1    2
            line2: Pro2 2    0    1.2
            line3: Pro3 1    4.2  0


  4         (n+1) by (n) matrices: the first row contains node ID
                                     contains node ID
  e.g.:     line1: Pro1 Pro2 Pro3
            line2: 2    2.3  1.2
            line3: 1    4.2  0.1
            line4: 3    1.2  0.3

  5         (n) by (n) matrices: the pure adjacency matrix, IDs are
                                 generated in format "Pro$i", $i=1,2,...,n
  e.g.:     line1: 2.3    2.3  1.2
            line2: 1.6    4.2  0.1
            line3: 3.1    1.2  0.3

 Note: no colons are allowed following each ID

'''

import re
import numpy as np
import collections
import sys


def suffixDetector(filename):
    filetype = 0
    if((re.match('[\w\\\/]*.csv$', filename) is not None) or
      (re.match('[\w\\\/]*.tab$', filename) is not None)):
        filetype = 2
    elif re.match('[\w\\\/]*.list$', filename) is not None:
        filetype = 1
    return filetype


def CheckValid(filename):
    '''
    it checks if the file is a valid PPI file and then detect the type

    filename -- the input file name

    returns filetype number, the number of invalid lines from the biginning
    and the number of valid entries

    '''
    filetype = 0
    nSkip = 0
    N = 0

    file = open(filename, 'rb')

    validpattern = re.compile('^[\w _\-.,\t":]+$')
    flag = True
    nLine = -1
    temp = " "
    while(flag and temp != ""):
        temp = file.readline()
        temp = temp.strip(' \n\t\r')
        flag = re.search(validpattern, temp) is None
        nLine = nLine + 1

    nSkip = nLine
    sFirstValidLine = temp
    temp = file.readline()
    if temp == "":
        print "file is too short"
        return(0, 0, 0)
    sSecondValidLine = temp.strip(' \n\t\r')
    if re.match(validpattern, sSecondValidLine) is None:
        print 'line: ', nLine, ' has abnormal characters'
        return(0, 0, 0)

    splitpattern = re.compile('[\t ,]+')
    numericpattern = re.compile('^[0-9. \t,\-]+$')

    allWordsFromFirst = re.split(splitpattern, sFirstValidLine)
    allWordsFromSecond = re.split(splitpattern, sSecondValidLine)

    n1 = np.size(allWordsFromFirst)
    n2 = np.size(allWordsFromSecond)

##########################################################################
### for format No. 2
##########################################################################
    if n1 != n2:
        if n1 != n2-1:
            print 'sizes of first valid line and second are very different'
            return(0, 0, 0)
        N = n1
        filetype = 2
        index = 1
        temp = file.readline()
        temp = temp.strip('\n\t\r ')

        while(temp != ""):
                #print 'Error: expect more lines'
                #return(0, 0 ,0)
            allwords = re.split(splitpattern, temp)
            if np.size(allwords) != N+1:
                temp = 'Error: line:' + str(index+2)
                temp = temp + ' has inconsistent number of values'
                return(0, 0, 0)
            else:
                index = index + 1
            temp = file.readline()
            temp = temp.strip('\n\t\r ')

        return(filetype, nSkip, N)

##########################################################################
### for format No. 5
##########################################################################
    if re.match(numericpattern, sFirstValidLine) is not None:
        if re.match(numericpattern, sSecondValidLine) is not None:
            filetype = 5
            N = np.size(allWordsFromSecond)
            index = 2
            while(index < N):
                temp = file.readline()
                temp = temp.strip('\n\t\r ')
                if re.match(numericpattern, temp) is None:
                    print 'line: ', index+1, ' has abnormal values'
                    return(0, 0, 0)
                else:
                    allwords = re.split(splitpattern, temp)
                    if np.size(allwords) != N:
                        print 'line: ', index+1, ' has inconsistent values'
                        return(0, 0, 0)
                    else:
                        index = index + 1
            return(filetype, nSkip, N)
        else:
            print 'the first valid line is pure numeric but the second is not'
            return(0, 0, 0)

##########################################################################
### for format No. 1, 3, 4
##########################################################################
# if flag is true, the corresponding entry is pure numeric
#    flag00 = re.search(numericpattern, allWordsFromFirst[0]) is not None
    flag01 = re.search(numericpattern, allWordsFromFirst[1]) is not None
    flag10 = re.search(numericpattern, allWordsFromSecond[0]) is not None
#    flag11 = re.search(numericpattern, allWordsFromSecond[1]) is not None

    if flag01:
        filetype = 3
        N = n1 - 1
    else:
        if flag10:
            filetype = 4
            N = n1
        else:
            filetype = 1
            N = 2
            temp = file.readline()
            while(temp != ""):
                temp = file.readline()
                N = N + 1

    file.close()
    return (filetype, nSkip, N)


def parseInteractionList(filename, nSkip, N):
    '''
    it parses PPIs from filename

    filename -- the input file name

    nSkip -- the number of invalid lines from the beginning

    N -- the number of node

    returns ppbAdj -- adjacency matrix

    '''
    print 'Parsing interaction list file...'
### collect node names ###
    infile = open(filename, 'r')
    index = 0
    while(index < nSkip):
        temp = infile.readline()
        index = index + 1

    splitpattern = re.compile('[\t ,]+')
    numericpattern = re.compile('^[0-9. \t,\-]+')
    index = 0
    temp = infile.readline()
    temp = temp.strip('\t \n\r')
    allwords = re.split(splitpattern, temp)
    names = {allwords[0]: index}
    index = index + 1
    if allwords[1] not in names:
            names[allwords[1]] = index
            index = index + 1
    for i in xrange(0, 2):
        if re.search(numericpattern, allwords[i]) is not None:
            temp = "Error: file inconsistenc (possible node ID is numeric)"
            print >> sys.stderr, temp
            exit(1)
    count = 1
    while(count < N):
        temp = infile.readline()
        temp = temp.strip('\t \n\r')
        allwords = re.split(splitpattern, temp)
        for i in xrange(0, 2):
            if re.match(numericpattern, allwords[i]) is not None:
                temp = "Error: file inconsistent(possible node ID is numeric)"
                temp = temp + "       " + allwords[i]
                print >> sys.stderr, temp
                exit(1)
        if allwords[0] not in names:
            names[allwords[0]] = index
            index = index + 1
        if allwords[1] not in names:
            names[allwords[1]] = index
            index = index + 1
        count = count + 1
    infile.close()
    names = collections.OrderedDict(sorted(names.items(), key=lambda x: x[1]))
    N = index
    ppbAdj = np.zeros((N, N))
### collect edges ###
    infile = open(filename, 'r')
    for i in xrange(0, nSkip):
        temp = infile.readline()
    count = 0
    for temp in infile:
        temp = temp.strip('\t \n\r')
        if temp == "":
            continue
        allwords = re.split(splitpattern, temp)
#        print temp
#        print allwords[0]
#        print allwords[1]
        i = names[allwords[0]]
        j = names[allwords[1]]
        if(len(allwords) > 2):
            if re.search(numericpattern, allwords[2]) is not None and i != j:
                    ppbAdj[i, j] = float(allwords[2])
        else:
                    ppbAdj[i, j] = 1
    infile.close()
 #   print ppbAdj
 #   print names
    return (ppbAdj, names)


def parseFullMatrix(filename, nSkip, N):
    '''
    it parses PPIs from filename

    filename -- the input file name

    nSkip -- the number of invalid lines from the beginning

    N -- the number of node

    returns ppbAdj -- adjacency matrix

    '''
    print 'Parsing interaction matrix file with IDs'
    print '    at both first cols and first rows...'
### collect node names ###
    file = open(filename, 'r')
    index = 0
    while(index < nSkip):
        temp = file.readline()
        index = index + 1
    numericpattern = re.compile('^[0-9. \t,\-]+')
    splitpattern = re.compile('[\t ,]+')
    index = 0
    temp = file.readline()
    temp = temp.strip('\t \n\r')
    allwords = re.split(splitpattern, temp)
    names = {allwords[0]: index}
    if re.match(numericpattern, allwords[0]) is not None:
        temp = "Error: file inconsistency (possible: node ID is numeric)"
        temp = temp + "       " + allwords[0]
        print >> sys.stderr, temp
        exit(1)
    for i in xrange(1, N):
        if re.match(numericpattern, allwords[i]) is not None:
            temp = "Error: file inconsistency (possible: node ID is numeric)"
            temp = "       " + allwords[i]
            print >> sys.stderr, temp
            exit(1)
        if allwords[i] in names:
            temp = 'Error: ' + allwords[i]
            temp = temp + ' appears more than once at first row'
            print >> sys.stderr, temp
            exit(1)
        names[allwords[i]] = i
    temp = file.readline()
    temp = temp.strip('\t \n\r')
    NN = N
    while(temp != ""):
        allwords = re.split(splitpattern, temp)
        if re.match(numericpattern, allwords[0]) is not None:
            temp = "Error: file inconsistency (possible: node ID is numeric)"
            temp = temp + "       " + allwords[i]
            print >> sys.stderr, temp
            exit(1)
        if allwords[0] not in names:
            names[allwords[0]] = N
            NN = NN + 1
        temp = file.readline()
        temp = temp.strip('\t \n\r')

    file.close()
    names = collections.OrderedDict(sorted(names.items(), key=lambda x: x[1]))

    Ntemp = N
    N = NN
    NN = Ntemp
    file = open(filename, 'r')
    index = 0
    while(index < nSkip):
        temp = file.readline()
        index = index + 1
    temp = file.readline()
    ppbAdj = np.zeros((N, N))

    index = 0
    while(index < N and temp != ""):
        temp = file.readline()
        temp = temp.strip('\t \n\r')
        if temp == "":
            break
        # print temp
        allwords = re.split(splitpattern, temp)
        i = names[allwords[0]]
        for j in xrange(0, NN):
            if allwords[j+1] != '0':
                ppbAdj[i, j] = float(allwords[j+1])
        ppbAdj[i, i] = 0
        index = index + 1

    file.close()
    return (ppbAdj, names)


def parseColMatrix(filename, nSkip, N):
    '''
    it parses PPIs from filename

    filename -- the input file name

    nSkip -- the number of invalid lines from the beginning

    N -- the number of node

    returns ppbAdj -- adjacency matrix

    '''
    print 'Parsing interaction matrix file with IDs'
    print '    at first cols...'
### collect node names ###
    ifile = open(filename, 'r')
    index = 0
    while(index < nSkip):
        temp = ifile.readline()
        index = index + 1
    numericpattern = re.compile('^[0-9. \t,\-]+')
    splitpattern = re.compile('[\t ,]+')
    index = 0
    temp = ifile.readline()
    temp = temp.strip('\t \n\r')
    allwords = re.split(splitpattern, temp)
    names = {allwords[0]: index}
    if re.match(numericpattern, allwords[0]) is not None:
        temp = "Error: file inconsistency (possible: node ID is numeric)"
        temp = temp + "       " + allwords[0]
        print >> sys.stderr, temp
        exit(1)
    temp = ifile.readline()
    temp = temp.strip('\t \n\r')
    while(temp != ""):
        index = index + 1
        allwords = re.split(splitpattern, temp)
        if re.match(numericpattern, allwords[0]) is not None:
            temp = "Error: file inconsistency (possible: node ID is numeric)"
            temp = temp + "       " + allwords[0]
            print >> sys.stderr, temp
            exit(1)
        if allwords[0] not in names:
            names[allwords[0]] = index
        else:
            temp = "Error: " + allwords[0]
            temp = temp + " appears more than once in first column"
            print >> sys.stderr, temp
            exit(1)
        temp = ifile.readline()
        temp = temp.strip('\t \n\r')
    if index+1 != N:
        temp = ("Error: inconsistent number of rows (" +
                N + ") and cols (" + index + ")")
        print >> sys.stderr, temp
        exit(1)
    ifile.close()
    names = collections.OrderedDict(sorted(names.items(),
                                           key=lambda x: x[1]))

    ifile = open(filename, 'r')
    index = 0
    while(index < nSkip):
        temp = file.readline()
        index = index + 1
    ppbAdj = np.zeros((N, N))
    index = 0
    temp = " "
    while(index < N and temp != ""):
        temp = ifile.readline()
        temp = temp.strip('\t \n\r')
        allwords = re.split(splitpattern, temp)
        i = names[allwords[0]]
        for j in xrange(0, N):
            if allwords[j+1] != '0':
                ppbAdj[i, j] = float(allwords[j+1])
        ppbAdj[i, i] = 0
        index = index + 1
    #print ppbAdj
    ifile.close()
    return (ppbAdj, names)


def parseRowMatrix(filename, nSkip, N):
    '''
    it parses PPIs from filename

    filename -- the input file name

    nSkip -- the number of invalid lines from the beginning

    N -- the number of node

    returns ppbAdj -- adjacency matrix

    '''
    print 'Parsing interaction matrix file with IDs'
    print '    at first rows...'
### collect node names ###
    file = open(filename, 'r')
    index = 0
    while(index < nSkip):
        temp = file.readline()
        index = index + 1
    numericpattern = re.compile('^[0-9. \t,\-]+')
    splitpattern = re.compile('[\t ,]+')
    index = 0
    temp = file.readline()
    temp = temp.strip('\t \n\r')
    allwords = re.split(splitpattern, temp)
    names = {allwords[0]: index}
    if re.match(numericpattern, allwords[0]) is not None:
        temp = "Error: file inconsistency (possible: node ID is numeric)"
        temp = temp + "       " + allwords[0]
        print >> sys.stderr, temp
        exit(1)
    for i in xrange(1, N):
        if re.match(numericpattern, allwords[i]) is not None:
            temp = "Error: file inconsistency (possible: node ID is numeric)"
            temp = temp + "       " + allwords[i]
            print >> sys.stderr, temp
            exit(1)
        if allwords[i] in names:
            temp = 'Error: ' + allwords[i]
            temp = temp + ' appears more than once at first row'
            print >> sys.stderr, temp
            exit(1)
        names[allwords[i]] = i
    names = collections.OrderedDict(sorted(names.items(),
                                           key=lambda x: x[1]))
    ppbAdj = np.zeros((N, N))
    trackline = nSkip
    for i in xrange(0, N):
        trackline = trackline + 1
        temp = file.readline()
        temp = temp.strip('\t \n\r')
        if temp == "":
            print >> sys.stderr, 'Error: expect more lines'
            exit(1)
        allwords = re.split(splitpattern, temp)
        if np.size(allwords) != N:
            temp = 'Error: line' + trackline
            temp = temp + ' has inconsistent number of entries'
            print >> sys.stderr, temp
            exit(1)
        for j in xrange(0, N):
            if allwords[j] != '0':
                ppbAdj[i, j] = float(allwords[j])
        ppbAdj[i, i] = 0
    file.close()
    return (ppbAdj, names)


def parseNoIDmatrix(filename, nSkip, N):
    '''
    it parses PPIs from filename

    filename -- the input file name

    nSkip -- the number of invalid lines from the beginning

    N -- the number of node

    returns ppbAdj -- adjacency matrix

    '''
    print 'Parsing interaction matrix file with no IDs'
    ### collect node names ###
    file = open(filename, 'r')
    index = 0
    while(index < nSkip):
        temp = file.readline()
        index = index + 1
    splitpattern = re.compile('[\t ,]+')
    names = {}
    for i in xrange(1, N+1):
        names[('Pro%04d' % i)] = i-1
    names = collections.OrderedDict(sorted(names.items(),
                                           key=lambda x: x[1]))
    trackline = nSkip
    ppbAdj = np.zeros((N, N))
    for i in xrange(0, N):
        trackline = trackline + 1
        temp = file.readline()
        temp = temp.strip('\t \n\r')
        if temp == "":
            print >> sys.stderr, 'Error: expect more lines'
            exit(1)
        allwords = re.split(splitpattern, temp)
        if np.size(allwords) != N:
            temp = 'Error: line' + trackline
            temp = temp + ' has inconsistent number of entries'
            exit(1)
        for j in xrange(0, N):
            if allwords[j] != '0':
                ppbAdj[i, j] = float(allwords[j])
        ppbAdj[i, i] = 0
    file.close()
    return (ppbAdj, names)


def writePPIforDSD(ppbAdj, names, ofile):
    '''
    it writes PPI list into ofile

    ppbAdj -- adjacency matrix

    names -- node ID string list

    ofile -- output file object

    returns true if correctly write all PPIs

    '''
    n = len(ppbAdj[0])
    for i in xrange(0, n):
        temp = ""
        for j in xrange(0, n):
            if ppbAdj[i, j] != 0:
                temp = temp + names[i] + '\t' + names[j]
                temp = temp + '\t' + str(ppbAdj[i, j]) + '\n'
        ofile.write(temp)
    return True
