#!/usr/bin/env python2.7

# Huge memory requirement, must be run on a server like pulsar or meteor.

import os, sys, getopt, numpy, time

def main(argv):
	dsdFiles = []
	dsdPath = ''
	numDsdFiles = 0
	outputFile = ''
	masterMatrixFile = ''
	threshold = -1

	# Read in command line args and process
	try:
		opts, args = getopt.getopt(argv,"hi:o:m:t:")
	except getopt.GetoptError:
		print ''
		print 'USAGE: ./DsdAverage-matrix-threshold.py <options>'
		print 'Use -h for more information about options.'
	for opt, arg in opts:
		if opt == '-h':
			print '-i  directory containing dsd files to be averaged'
			print '-o  name of output (averaged) file'
			print '-m  master matrix file. header contains all possible proteins that show up in the input files'
			print '-t  threshold value from 0 to 1'
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
		elif opt == '-t':
			threshold = arg

	if threshold == -1:
		print "Please specify threshold using -t. Value should be from 0 to 1."
		exit()


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
		with open(dsdFileWithPath) as dsdFile:
			localLabels = dsdFile.readline().split('\t')[1:]
			localLabels = [label.strip() for label in localLabels]
		numLocalLabels = len(localLabels)
		dsdMatrix = numpy.loadtxt(dsdFileWithPath, delimiter='\t', usecols=xrange(1,numLocalLabels+1), skiprows=1)
		print "  Time to read file into numpy array: " + str(time.clock() - time_loadfile)

		for index, score in numpy.ndenumerate(dsdMatrix):
			if index[1] > index[0]: # for humans: col > row, upper triangle only
				# Look up the position of each element in the master matrix
				xMasterIdx = mapLabelToIdx[localLabels[index[0]]]
				yMasterIdx = mapLabelToIdx[localLabels[index[1]]]

				# Add each element to the giant matrix of all DSDs
				allDsds[xMasterIdx, yMasterIdx, currentCount] = float(score)

		print "  Time to process entire file:        " + str(time.clock() - time_processfile)
		currentCount += 1


	# ===== CALCULATE AVERAGE SCORES =====		
	print "Calculating average scores..."
	time_averagescores = time.clock()
	matrixAverageScores = numpy.zeros(shape=(numLabels, numLabels))
	for row_i in xrange(0, numLabels):
		for col_j in xrange(0, numLabels):
			if col_j == row_i: # diagonal
				matrixAverageScores[row_i, col_j] = 0.0
			if col_j > row_i: # upper triangle
				dsdScores = filter(None, sorted(allDsds[row_i][col_j][:], reverse=True))
				numScores = len(dsdScores)
				
				if numScores > 0:
					cutoffIdx = int(numScores * float(threshold)) + 1 # +1 to get ceiling
					dsdScores = dsdScores[:cutoffIdx]
					average = sum(dsdScores)/cutoffIdx
					matrixAverageScores[row_i, col_j] = average
					matrixAverageScores[col_j, row_i] = average
				else:
					matrixAverageScores[row_i, col_j] = 999.9
					matrixAverageScores[col_j, row_i] = 999.9

	print "  Time to average scores:             " + str(time.clock() - time_averagescores)			
				

	# ===== WRITE RESULTS TO FILE =====
	print "Writing results to file..."
	time_writetofile = time.clock()
	labels_row = numpy.array((labels), dtype='|S12')[numpy.newaxis]
	matrixAverageScores = numpy.concatenate((labels_row, matrixAverageScores), 0)
	labels_col = numpy.insert(labels_row, 0, " ")[numpy.newaxis].T
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