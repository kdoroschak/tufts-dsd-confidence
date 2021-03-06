#!/usr/bin/env python2.7

# Huge memory requirement, must be run on a server like pulsar or meteor.

import os, sys, getopt, numpy, time, bisect

def main(argv):
	dsdFiles = []
	dsdPath = 'Please specify a path to the input files using -i'
	numDsdFiles = 0
	outputFile = 'Please specify an output file using -o'
	masterMatrixFile = 'Please specify a master matrix file using -m'

	# Read in command line args and process
	try:
		opts, args = getopt.getopt(argv,"hi:o:m:b:")
	except getopt.GetoptError:
		print ''
		print 'USAGE: ./DsdAverage-matrix-maxtopandbottom.py <options>'
		print 'Use -h for more information about options.'
	for opt, arg in opts:
		if opt == '-h':
			print '-i  directory containing dsd files to be averaged'
			print '-o  name of output (averaged) file'
			print '-m  master matrix file. header contains all possible proteins that show up in the input files'
			print '-b  bisect mode. 0 to bisect along original DSD. 1 for most extreme 1/4. 2 for 3/4.'
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
		elif opt == '-b':
			bisectMode = int(arg)

	time_entireprogram = time.clock()
	# ===== CREATE MATRIX FRAMEWORK =====
	# Read in header of master file
	with open(masterMatrixFile) as masterMatrix:
		labels = masterMatrix.readline()
	labels = labels.split('\t')[1:]
	labels = [label.strip() for label in labels]
	numLabels = len(labels)
	masterMatrix = numpy.loadtxt(masterMatrixFile, delimiter='\t', usecols=xrange(1,numLabels+1), skiprows=1)
	
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

#		dsdMatrix = numpy.loadtxt(dsdFileWithPath, delimiter='\t', usecols=xrange(1,numLocalLabels+1), skiprows=1)
#		print "  Time to read file into numpy array: " + str(time.clock() - time_loadfile)

#		for index, score in numpy.ndenumerate(dsdMatrix):
#			if index[1] > index[0]: # for humans: col > row, upper triangle only
#				# Look up the position of each element in the master matrix
#				xMasterIdx = mapLabelToIdx[localLabels[index[0]]]
#				yMasterIdx = mapLabelToIdx[localLabels[index[1]]]
#
#				# Add each element to the giant matrix of all DSDs
#				allDsds[xMasterIdx, yMasterIdx, currentCount] = float(score)
#
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
				# dsdScores sorted !!!!!!!!!!!!! LOW TO HIGH !!!!!!!!!!!!!!!!
				dsdScores = filter(None, sorted(allDsds[row_i][col_j][:]))
				numScores = len(dsdScores)
				
				if numScores > 0:
					# Avg of all DSD scores
					fullAverage = sum(dsdScores)/numScores
					originalDsdScore = masterMatrix[row_i][col_j]

					#numHigherScores = bisect.bisect(dsdScores, fullAverage)
					numLowerScores = bisect.bisect(dsdScores, originalDsdScore)
					numHigherScores = numScores - numLowerScores
					#print str(dsdScores[:numLowerScores]) + " < " + str(originalDsdScore) + " < " + str(dsdScores[numLowerScores:])

					if bisectMode == 0:
						upperPartition = numLowerScores
						lowerPartition = numLowerScores
					elif bisectMode == 1:
						upperPartition = numLowerScores + int(numHigherScores / 2)
						lowerPartition = numLowerScores/2 + 1
					elif bisectMode == 2:
						upperPartition = numLowerScores / 2
						lowerPartition = numLowerScores + (numHigherScores / 2) + 1

					# More scores > original
					if numHigherScores > numLowerScores:
						topDsdScores = dsdScores[upperPartition:]
						topAverage = sum(topDsdScores)/len(topDsdScores)
						matrixAverageScores[row_i, col_j] = topAverage
						matrixAverageScores[col_j, row_i] = topAverage

					# More scores < original
					elif numHigherScores < numLowerScores:
						bottomDsdScores = dsdScores[:lowerPartition]
						bottomAverage = sum(bottomDsdScores)/len(bottomDsdScores)
						matrixAverageScores[row_i, col_j] = bottomAverage
						matrixAverageScores[col_j, row_i] = bottomAverage
					
					# Equal scores > and < original, use regular average
					else:
						matrixAverageScores[row_i, col_j] = fullAverage
						matrixAverageScores[col_j, row_i] = fullAverage

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
