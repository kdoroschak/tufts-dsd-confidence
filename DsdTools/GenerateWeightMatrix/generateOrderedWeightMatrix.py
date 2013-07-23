#!/usr/bin/env python2.7

import os, sys, argparse, numpy as np, re, collections

def main(argv):
	# Build command line argument parser
	parser = argparse.ArgumentParser(description='Take the average of multiple DSD files using one of several averaging methods.')
	parser.add_argument('-o', required=True, type=str, help='name of the output transition matrix file (tab-delimited 2D matrix)')
	parser.add_argument('-i', required=True, type=str, help='name of PPI input file. must be tab-delimited, fully connected, and have confidence score as 3rd column')
	parser.add_argument('-p', required=True, type=str, help='name of the ordered list of proteins, for ordering the matrix')
	args = parser.parse_args()
	ppiInfile = args.i
	outputFile = args.o
	orderedProteinFile = args.p
	
	names = getNames(orderedProteinFile)	
	(conf, names) = getConf(ppiInfile, names)

	np.savetxt(outputFile, conf, delimiter='\t', fmt='%5f')	

def getNames(proteinListFile):
	names = {}
	index = 0
	with open(proteinListFile) as proteinFile:
		for protein in proteinFile:
			protein = protein.strip()
			names[protein] = index
			index += 1
	names = collections.OrderedDict(sorted(names.items(), key=lambda x: x[1]))
	return names


def getConf(ppiInfile, names):

    validpattern = re.compile('^[\w _\-.,\t"\':;]+$')
    splitpattern = re.compile('[\t ;,]+')
    numericpattern = re.compile('^[0-9. \t,\-]+')
	
    N = len(names)
    ppbConf = np.zeros((N, N))    
    
	### collect edges ###
    ppiFile = open(ppiInfile, 'r')
    for line in ppiFile:
        line = line.strip('\t \n\r')
        if line == "" or (re.search(validpattern, line) is None):
            continue
        
        allwords = re.split(splitpattern, line)
        i = names[allwords[0]]
        j = names[allwords[1]]
        
        confidence = float(allwords[2])
        
        if i != j and i != None and j != None:
            ppbConf[i, j] = confidence
            ppbConf[j, i] = confidence
    ppiFile.close()
    
    return (ppbConf, names)
    
	
	
if __name__ == "__main__":
	main(sys.argv[1:])
