#!/usr/bin/env python2.7

# Huge memory requirement, must be run on a server like pulsar or meteor.

import os, sys, getopt, numpy, time

def main(argv):
	dsdFiles = []
	dsdPath = 'Please specify a path to the input files using -i'
	numDsdFiles = 0
	outputFile = 'Please specify an output file using -o'
	masterMatrixFile = 'Please specify a master matrix file using -m'
	threshold = 0.5

	# Read in command line args and process
	try:
		opts, args = getopt.getopt(argv,"hi:o:m:")
	except getopt.GetoptError:
		print ''
		print 'USAGE: ./DsdAverage-matrix-maxtopandbottom.py <options>'
		print 'Use -h for more information about options.'
	for opt, arg in opts:
		if opt == '-h':
			print '-i  directory containing dsd files to be averaged'
			print '-o  name of output (averaged) file'
			print '-m  master matrix file. header contains all possible proteins that show up in the input files'
			sys.exit(1)
		elif opt == '-i':
			try:
				dsdPath = arg
				dsdFiles = os.listdir(arg)
				numDsdFiles = len(dsdFiles)
			except getopt.GetoptError:
				print 'Error opening path to files.'
		elif opt == '-o':
			outputFile = arg
		elif opt == '-m':
			masterMatrixFile = arg

#	if threshold > 1 or threshold < 0:
#		print "Please specify threshold using -t. Value should be from 0 to 1."
#		exit()


	time_entireprogram = time.clock()
	# ===== CREATE MATRIX FRAMEWORK =====
	# Read in header of master file
	with open(masterMatrixFile) as masterMatrix:
		labels = masterMatrix.readline()
	labels = labels.split('\t')[1:]
	labels = [label.strip() for label in labels]
	numLabels = len(labels)
		
	# Create map from label -> master index
	mapLabelToIdx = {}
	for label_idx in xrange(numLabels):
		mapLabelToIdx[labels[label_idx].strip()] = label_idx

	# Create giant matrix using all edges in header, initialized to 0.0
	allDsds = numpy.zeros(shape=(numLabels, numLabels, numDsdFiles), dtype=float)

	# ===== POPULATE GIANT MATRIX =====
	# Add each element from individual files to giant matrix
	print "Processing all files at " + dsdPath + "..."
	dsdFiles.sort()
	currentCount = 0
	for currentDsdFile in dsdFiles:
		time_processfile = time.clock()
		dsdFileWithPath = dsdPath + "/" + currentDsdFile
		print "Processing DSD file " + str(currentCount + 1) + "/" + str(numDsdFiles) + ": " + dsdFileWithPath	
		
		# Read file to be averaged into numpy array
		time_loadfile = time.clock()
		dsdFile = open(dsdFileWithPath)
		localLabels = dsdFile.readline().split('\t')[1:]
		localLabels = [label.strip() for label in localLabels]
		numLocalLabels = len(localLabels)

		dsdMatrix = numpy.loadtxt(dsdFileWithPath, delimiter='\t', usecols=xrange(1,numLocalLabels+1), skiprows=1, dtype=float)

		print "  Time to read file into numpy array: " + str(time.clock() - time_loadfile)

		time_lineupfile = time.clock()
		# Line up the labels and populate the big array
		for i in xrange(numLocalLabels):
			for j in xrange(numLocalLabels):
				if i < j:
					xMasterIdx = mapLabelToIdx[localLabels[i]]
					yMasterIdx = mapLabelToIdx[localLabels[j]]
					score = dsdMatrix[i][j]
					allDsds[xMasterIdx, yMasterIdx, currentCount] = score
					
		print "  Time to line up protein labels:     " + str(time.clock() - time_lineupfile)
		print "  Time to process entire file:        " + str(time.clock() - time_processfile)
		currentCount += 1


	# ===== CALCULATE AVERAGE SCORES =====		
	print "Calculating average scores..."
	time_averagescores = time.clock()
	bins = xrange(25) # TODO un-hardcode this
	matrixAverageScores = numpy.zeros(shape=(numLabels, numLabels))
	
	for row_i in xrange(0, numLabels):
		for col_j in xrange(0, numLabels):
			if col_j == row_i: # diagonal
				matrixAverageScores[row_i, col_j] = 0.0
			if col_j > row_i and row_i < 20: # upper triangle
				dsdScores = filter(None, sorted(allDsds[row_i][col_j][:], reverse=True))
				hist,binEdges = numpy.histogram(dsdScores, bins=bins)
				#print hist
				#print binEdges
				#maxCount = max(hist)
				#print maxCount
				maxCountIdx = numpy.where(hist == hist.max())[0] # TODO decide how to break ties
				upperBinEdge = binEdges[maxCountIdx+1][0] # TODO do NOT use just the first idx
				lowerBinEdge = binEdges[maxCountIdx][0]
				#print upperBinEdge, lowerBinEdge
				dsdScores = [x for x in dsdScores if (lowerBinEdge <= x) and (x < upperBinEdge)]
				#print len(dsdScores)
				average = numpy.mean(dsdScores)
				matrixAverageScores[row_i, col_j] = average
				matrixAverageScores[col_j, row_i] = average	
				
	
	print "  Time to average scores:             " + str(time.clock() - time_averagescores)			
				

	# ===== WRITE RESULTS TO FILE =====
	print "Writing results to file..."
	time_writetofile = time.clock()
	# Create a numLabels x 1 array of labels for the top row
	labels_row = numpy.array((labels), dtype='|S12')[numpy.newaxis]
	# Concatenate the row of labels to the matrix along the x axis
	matrixAverageScores = numpy.concatenate((labels_row, matrixAverageScores), 0)
	# Add the corner blank spot to the labels and rotate it to use as the left column
	labels_col = numpy.insert(labels_row, 0, " ")[numpy.newaxis].T
	# Concatenate the column of labels to the matrix along the y axis
	matrixAverageScores = numpy.concatenate((labels_col, matrixAverageScores), 1)

	with open(outputFile, 'w') as outFile:
		for row in matrixAverageScores:
			outFile.write('\t'.join(row))
			outFile.write('\n')		

	print "  Time to write to file:              " + str(time.clock() - time_writetofile)
	print "Done."
	print "  Total execution time:           " + str(time.clock() - time_entireprogram)

if __name__ == "__main__":
	main(sys.argv[1:])
