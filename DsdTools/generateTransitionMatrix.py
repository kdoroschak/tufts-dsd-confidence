#!/usr/bin/env python2.7

import os, sys, argparse, numpy as np, re, collections, time

def main(argv):
	# Build command line argument parser
	parser = argparse.ArgumentParser(description='Take the average of multiple DSD files using one of several averaging methods.')
	parser.add_argument('-o', required=True, type=str, help='name of the output transition matrix file (tab-delimited 2D matrix)')
	parser.add_argument('-i', required=True, type=str, help='name of PPI input file. must be tab-delimited, fully connected, and have confidence score as 3rd column')
	args = parser.parse_args()
	ppiFile = args.i
	outputFile = args.o
	
	conftime = time.clock()
	(conf, names) = getConf(ppiFile)
	print "Getting confidence matrix took " + str(time.clock() - conftime) + "s."

	transtime = time.clock()
	transList = getTransList(conf, names)
	print "Getting transition list took " + str(time.clock() - transtime) + "s."
	conf = ""
	names = ""

	#transtime = time.clock()
	#(transMatrix, names) = getTransMatrix(conf, names)
	#print "Getting transition matrix took " + str(time.clock() - transtime) + "s."
		
	writetime = time.clock()
	#np.savetxt(outputFile, trans, delimiter='\t')
	with open(outputFile, 'w') as outFile:
		outFile.write(transList)	
	print "Writing the matrix took " + str(time.clock() - writetime) + "s."
	transList = ""


def getConf(ppiFile):
    """
    ppiFile - the name of input file to be parsed, which should be
               a file with one PPI at each line

    returns ppbConf, an adjacency matrix represented as a numpy array
    """
    validpattern = re.compile('^[\w _\-.,\t"\':;]+$')
    splitpattern = re.compile('[\t ;,]+')
    numericpattern = re.compile('^[0-9. \t,\-]+')
	### collect node names ###
    finfile = open(ppiFile, 'r')
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
    ppbConf = np.zeros((N, N))
	### collect edges ###
    finfile = open(ppiFile, 'r')
    for temp in finfile:
        temp = temp.strip('\t \n\r')
        if temp == "" or (re.search(validpattern, temp) is None):
            continue
        allwords = re.split(splitpattern, temp)
        i = names[allwords[0]]
        j = names[allwords[1]]
        templen = len(temp)
        # allow for tab or space delimiters
        temp = temp.split('\t')
        if templen == len(temp):
        	temp = temp.split(' ')
        	
        if len(temp) > 2:
            confidence = temp[2]
        else:
            confidence = 1
        if i != j:
            ppbConf[i, j] = confidence
            ppbConf[j, i] = confidence
    finfile.close()
    return (ppbConf, names)
    
    
def getTransMatrix(ppbConf, names):
	numNodes = np.size(ppbConf[0])
	trans = np.zeros((numNodes, numNodes))
	
	for row in xrange(0, numNodes):
		sumRow = sum(ppbConf[row])
		if sumRow != 0:
			trans[row] = ppbConf[row]/sumRow
	return (trans, names)
	

def getTransList(ppbConf, names):
	numNodes = np.size(ppbConf[0])
	trans = ""
	accessNames = names.items()
	for row in xrange(0, numNodes):
		sumRow = sum(ppbConf[row])
		if sumRow != 0:
			transRow = ppbConf[row]/sumRow
			print transRow
		for col in xrange(0, numNodes):
			proteinA = accessNames[row][0]
			proteinB = accessNames[col][0]
			transval = str(transRow[col]) + "\n"
			edge = "\t".join([proteinA, proteinB, transval])
			#print edge
			trans += edge
	return trans		

	
if __name__ == "__main__":
	main(sys.argv[1:])
