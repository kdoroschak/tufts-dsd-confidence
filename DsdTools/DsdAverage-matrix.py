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
		print '{' + labels[0] + '}'
	labels = labels.split('\t')[1:]
	numLabels = len(labels)
	
	# Create map from label -> master index
	mapLabelToIdx = {}
	for label_idx in xrange(numLabels):
		mapLabelToIdx[labels[label_idx].strip()] = label_idx

	# Create matrices using all edges in header, initialized to 0
	matrixTotalScores = numpy.zeros(shape=(numLabels, numLabels))
	matrixNumScores   = numpy.zeros(shape=(numLabels, numLabels))


	# ===== CALCULATE TOTAL SCORES AND COUNTS FOR ALL MATRICES =====
	# Add each new element to total and count
	dsdFiles.sort()
	numDsdFiles = len(dsdFiles)	
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
		numLocalLabels = len(localLabels)
		dsdMatrix = numpy.loadtxt(dsdFileWithPath, delimiter='\t', skiprows=1, usecols=xrange(1,numLocalLabels)) # TODO this may be a bug, cutting off last column!
		print "Read file into numpy array:\t" + str(time.clock() - time_loadfile)

		for index, score in numpy.ndenumerate(dsdMatrix):
			if index[1] > index[0]: # for humans: col > row, upper triangle only
				# Look up the position of each element in the master matrix
				xMasterIdx = mapLabelToIdx[localLabels[index[0]]]
				yMasterIdx = mapLabelToIdx[localLabels[index[1]]]

				# Add each element to the master total and count
				matrixTotalScores[xMasterIdx, yMasterIdx] += score
				matrixNumScores[xMasterIdx, yMasterIdx] += 1

		print "Process entire file:       \t" + str(time.clock() - time_processfile)
	

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
				matrixTotalScores[index] = 999 # sentinel value, effectively +inf in this case
	print "Average scores:            \t" + str(time.clock() - time_averagescores)
	

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


#numpy.savetxt(outputFile, matrixTotalScores, delimiter='\t', newline='\n')
#	print "(jk, not implemented yet)"


	print "Done."		

if __name__ == "__main__":
	main(sys.argv[1:])
