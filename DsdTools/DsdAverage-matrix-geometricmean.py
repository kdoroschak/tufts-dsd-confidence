#!/usr/bin/env python2.7

# Huge memory requirement, must be run on a server like pulsar or meteor.

import os, sys, argparse, numpy, time
import module_geometricmean as gm

def main(argv):
	time_entireprogram = time.clock()
	useMasterMatrix = False
	
	methodOptions = ['geometricmean', 'upperthreshold', 'arithmeticmean']

	# Build command line argument parser
	parser = argparse.ArgumentParser(description='Take the average of multiple DSD files using one of several averaging methods.')
	parser.add_argument('--dsdFiles', '-i', required=True, type=str, help='directory containing only the DSD files to be averaged')
	parser.add_argument('-o', required=True, type=str, help='name of the output (averaged) file')
	parser.add_argument('--originalDsd', '-d', required=True, type=str, help='master matrix file, typically the original DSD file. header contains all possible proteins that show up in the input files')
	parser.add_argument('--threshold', '-t', required=False, type=float, help='threshold value from 0 to 1 of top t%% of values to use in output')
	parser.add_argument('--method', '-m', required=True, choices=methodOptions, help='select the averaging method from the above choices')
	parser.add_argument('-a', '--xtraarg', required=False, type=str, help='additional argument if the specific method requires it. TODO add more info about which args are required')
	
	args = parser.parse_args()
	try:
		dsdPath = args.dsdFiles
		dsdFiles = os.listdir(dsdPath)
		numDsdFiles = len(dsdFiles)
	except:
		print 'Error opening path to DSD files.'
	outputFile = args.o
	masterMatrixFile = args.originalDsd
	threshold = args.threshold
	method = args.method
	arg=args.xtraarg
	
	# ========== IMPORT MODULE AND SET VARIABLES BASED ON METHOD ==========
	# Method: Geometric mean
	if method == "geometricmean":
		print "Using geometric mean module to average DSDs."
		if threshold is not None:
			print "  Warning: Threshold is not used and will be ignored."
		if arg is not None:
			print "  Warning: Adt'l argument is not used and will be ignored."
		useMasterMatrix = False
		import module_geometricmean as averager
	
	# Method: Upper threshold
	elif method == "upperthreshold":
		print "Using upper threshold module to average DSDs."
		if threshold is None:
			print "  Error: Please specify threshold."
			exit()
		if arg is not None:
			print "  Warning: Adt'l argument is not used and will be ignored."
		useMasterMatrix = False
		import module_upperthreshold as averager
		
	# Method: Arithmetic mean - all possible scores
	if method == "arithmeticmean":
		print "Using arithmetic mean module to average DSDs."
		if threshold is not None:
			print "  Warning: Threshold is not used and will be ignored."
		if arg is not None:
			print "  Warning: Adt'l argument is not used and will be ignored."
		useMasterMatrix = False
		import module_arithmeticmean as averager	
	
	# Template method
	elif method == "":
		print "Using ____ module to average DSDs."
		# is additional arg required?
		# is threshold required/not being used?
		# is master matrix needed?
		# import module_X as averager
	
	if (threshold > 1 or threshold < 0) and (threshold is not None):
		print "Error: Value of threshold should be from 0 to 1."
		exit()
	
	
	# ========== CREATE MATRIX FRAMEWORK ==========
	print "Creating the matrix framework...",
	sys.stdout.flush()
	time_framework = time.clock()
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
	print "   took " + str(time.clock() - time_framework) + "s."


	# ========== READ ORIGINAL DSD ==========
	if useMasterMatrix:
		print "Reading in the original DSD...",
		sys.stdout.flush()
		time_loadfile = time.clock()
		masterMatrix = numpy.loadtxt(masterMatrixFile, delimiter='\t', usecols=xrange(1,numLabels+1), skiprows=1)
		print "   took " + str(time.clock() - time_loadfile) + "s."
	else:
		masterMatrix = ''


	# ========== POPULATE GIANT MATRIX ==========
	# Add each element from individual files to giant matrix
	print "Processing all files at " + dsdPath + "."
	dsdFiles.sort()
	currentCount = 0
	for currentDsdFile in dsdFiles:
		time_processfile = time.clock()
		dsdFileWithPath = dsdPath + "/" + currentDsdFile
		print "Working on file " + str(currentCount + 1) + "/" + str(numDsdFiles) + ": " + dsdFileWithPath	
		
		# Read file to be averaged into numpy array
		print "  Reading file into numpy array...",
		sys.stdout.flush()
		time_loadfile = time.clock()
		dsdFile = open(dsdFileWithPath)
		localLabels = dsdFile.readline().split('\t')[1:]
		localLabels = [label.strip() for label in localLabels]
		numLocalLabels = len(localLabels)

		dsdMatrix = numpy.loadtxt(dsdFileWithPath, delimiter='\t', usecols=xrange(1,numLocalLabels+1), skiprows=1, dtype=float)

		print "   took " + str(time.clock() - time_loadfile) + "s."

		print "  Matching indices and populating big matrix...",
		sys.stdout.flush()
		time_lineupfile = time.clock()
		# Line up the labels and populate the big array
		for i in xrange(numLocalLabels):
			for j in xrange(numLocalLabels):
				if i < j:
					xMasterIdx = mapLabelToIdx[localLabels[i]]
					yMasterIdx = mapLabelToIdx[localLabels[j]]
					score = dsdMatrix[i][j]
					allDsds[xMasterIdx, yMasterIdx, currentCount] = score
					
		print "   took " + str(time.clock() - time_lineupfile) + "s."
		print "  Entire file took " + str(time.clock() - time_processfile) + "s."
		currentCount += 1
	dsdMatrix = '' # clean up memory


	# ========== CALCULATE AVERAGE SCORES ==========		
	print "Calculating average scores..."
	time_averagescores = time.clock()
		
	matrixAverageScores = averager.calculateAverage(allDsds, masterMatrix, threshold, arg)
	allDsds = '' # clean up memory
	masterMatrix = '' # clean up memory

	print "  Averaging took " + str(time.clock() - time_averagescores) +  "s."
				

	# ========== WRITE RESULTS TO FILE ==========
	print "Writing results to file...",
	sys.stdout.flush()
	time_writetofile = time.clock()
	labels_row = numpy.array((labels), dtype='|S12')[numpy.newaxis]
	matrixAverageScores = numpy.concatenate((labels_row, matrixAverageScores), 0)
	labels_col = numpy.insert(labels_row, 0, " ")[numpy.newaxis].T
	matrixAverageScores = numpy.concatenate((labels_col, matrixAverageScores), 1)

	with open(outputFile, 'w') as outFile:
		for row in matrixAverageScores:
			outFile.write('\t'.join(row))
			outFile.write('\n')		

	print "   took " + str(time.clock() - time_writetofile) + "s."
	print "Done."
	print "  Total execution time:           " + str(time.clock() - time_entireprogram)

if __name__ == "__main__":
	main(sys.argv[1:])
