#!/usr/bin/env python2.7

import os, sys, getopt, numpy, time

def main(argv):
	dsdFiles = []
	dsdPath = ''
	outputFile = ''
	masterMatrixFile = ''

	# Read in command line args and process
	try:
		opts, args = getopt.getopt(argv,"hi:o:m:")
	except getopt.GetoptError:
		print ''
		print 'USAGE: ./DsdAverage <options>'
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
			except getopt.GetoptError:
				print 'Error opening path to files.'
		elif opt == '-o':
			outputFile = arg
		elif opt == '-m':
			masterMatrixFile = arg


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

	# Create matrices using all edges in header, initialized to 0
	matrixTotalScores = numpy.zeros(shape=(numLabels, numLabels))
	matrixNumScores   = numpy.zeros(shape=(numLabels, numLabels))
	
	# Create giant matrix using all edges in header, initialized to 0.0
	numDsdFiles = len(dsdFiles)	
	allDsds = numpy.zeros(shape=(numLabels, numLabels, numDsdFiles), dtype=float)


	# ===== CALCULATE TOTAL SCORES AND COUNTS FOR ALL MATRICES =====
	# Add each new element to total and count
	dsdFiles.sort()
	print "Processing all files at " + dsdPath + "..."
	currentCount = 1
	for currentDsdFile in dsdFiles:
		dsdFileWithPath = dsdPath + "/" + currentDsdFile
		print "Processing DSD file " + str(currentCount) + "/" + str(numDsdFiles) + ": " + dsdFileWithPath	
		time_processfile = time.clock()
		
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
#		print "  Time to read file into numpy array:\t" + str(time.clock() - time_loadfile)

#		for index, score in numpy.ndenumerate(dsdMatrix):
#			if index[1] > index[0]: # for humans: col > row, upper triangle only
#				# Look up the position of each element in the master matrix
#				xMasterIdx = mapLabelToIdx[localLabels[index[0]].strip()]
#				yMasterIdx = mapLabelToIdx[localLabels[index[1]].strip()]
#
#				# Add each element to the master total and count
#				matrixTotalScores[xMasterIdx, yMasterIdx] += score
#				matrixNumScores[xMasterIdx, yMasterIdx] += 1

		print "  Time to process entire file:       \t" + str(time.clock() - time_processfile)
		currentCount += 1
	

	# ===== CALCULATE AVERAGE SCORES =====		
	print "Calculating average scores..."
	time_averagescores = time.clock()
	for index, score in numpy.ndenumerate(matrixTotalScores):
		if index[1] > index[0]: # col > row, upper triangle only. diag = 0
			index_t = (index[1], index[0])
			numScores = matrixNumScores[index]
			if numScores > 0:
				averageScore = score / matrixNumScores[index]
				matrixTotalScores[index]   = averageScore
				matrixTotalScores[index_t] = averageScore
			else:
				matrixTotalScores[index] = 999.9 # sentinel value, effectively +inf in this case
	print "  Time to average scores:            \t" + str(time.clock() - time_averagescores)
	

	# ===== WRITE RESULTS TO FILE =====
	print "Writing results to file..."
	labels_row = numpy.array((labels), dtype='|S12')[numpy.newaxis]
	matrixTotalScores = numpy.concatenate((labels_row, matrixTotalScores), 0)
	labels_col = numpy.insert(labels_row, 0, " ")[numpy.newaxis].T
	matrixTotalScores = numpy.concatenate((labels_col, matrixTotalScores), 1)

	with open(outputFile, 'w') as outFile:
		for row in matrixTotalScores:
			outFile.write('\t'.join(row))
			outFile.write('\n')		

	print "Done."		

if __name__ == "__main__":
	main(sys.argv[1:])
