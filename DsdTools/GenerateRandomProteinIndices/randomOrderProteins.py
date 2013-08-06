#!/usr/bin/env python2.7

from random import shuffle
import getopt, sys

def main(argv):
	outputFile = ''
	annotatedFile = ''

	# Read in command line args and process
	try:
		opts, args = getopt.getopt(argv,"hi:o:m:")
	except getopt.GetoptError:
		print ''
		print 'USAGE: ./randomOrderProteins <options>'
		print 'Use -h for more information about options.'
	for opt, arg in opts:
		if opt == '-h':
			print '-i  annotated protein file, in the same format as required by M. Cao\'s majority vote.'
			print '-o  name of output file (random list of indices corresponding to rows in the input file)'
			sys.exit(1)
		elif opt == '-i':
			annotatedFile = arg
		elif opt == '-o':
			outputFile = arg

	annotatedFile = open(annotatedFile, 'r')
	index = 1
	annotated_indices = []
	for line in annotatedFile:
		if len(line) > 0 and line.strip() != '1':
			annotated_indices.append(index)
			#print index
		index += 1
	annotatedFile.close()
	shuffle(annotated_indices)

	outputFile = open(outputFile, 'w')
	for item in annotated_indices:
		#print item
		outputFile.write(str(item) + '\n')

	print "Done."

if __name__ == "__main__":
	main(sys.argv[1:])
