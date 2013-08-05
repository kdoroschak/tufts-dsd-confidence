#!/usr/bin/env python2.7

import sys, argparse

def main(argv):
	# Build command line argument parser
	parser = argparse.ArgumentParser(description='Invert confidence scores by taking 1-score for each edge in the PPI.')
	parser.add_argument('-o', required=True, type=str, help='name of the output file (inverted PPI file)')
	parser.add_argument('-p', required=True, type=str, help='name of PPI input file. must be tab-delimited with confidence score as 3rd column')
	args = parser.parse_args()
	ppiFile = args.p
	outputFile = args.o

	ppiUnflipped = open(ppiFile, 'r')
	ppiFlipped   = open(outputFile, 'w')

	for line in ppiUnflipped:
		cols = line.split('\t')
		cols[2] = str(1 - float(cols[2]))
		ppiFlipped.write("\t".join(cols))
		ppiFlipped.write("\n")

if __name__ == "__main__":
	main(sys.argv[1:])
