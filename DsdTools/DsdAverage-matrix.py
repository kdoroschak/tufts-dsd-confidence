#!/usr/bin/env python2.7

import os, sys, getopt, numpy

def calculateNewScore(oldAvgScore, oldN, scoreToAdd):
	newAvgScore = ((oldAvgScore * int(oldN)) + float(scoreToAdd)) / (int(oldN) + 1)
	return newAvgScore

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
	# Need: masterMatrixFile, dsdPath

	# ===== CREATE MATRIX FRAMEWORK =====
	# Read in header of original file
	# Create matrix using all edges in header, initialized to 0
	# Create map from label name -> index
	masterMatrixFile = open(masterMatrixFile, 'r')
	mapLabelToIdx = {}
	i = 0
	for line in masterMatrixFile:
		if i > 0:
			break
		labels = line.split('\t')
		numLabels = len(labels)
		avgMatrixScores = numpy.zeros(shape=(numLabels, numLabels))
		avgMatrixNumScores = numpy.zeros(shape=(numLabels, numLabels), dtype='int')
		label_idx = 0
		for label in labels:
			mapLabelToIdx[label.strip()] = label_idx
			label_idx += 1
		i += 1
			
	# ===== POPULATE AVERAGE MATRIX =====
	# Read each file to be averaged
	print "Processing all files at " + dsdPath + "..."	
	for currentDsdFile in dsdFiles:
		print "Processing DSD file"
		mapIdxToLabel = {}
		dsdFile = open(dsdPath + "/" + currentDsdFile)
		print dsdPath + "/" + currentDsdFile
		row_i = -1
		for line in dsdFile:
			# Create map from index -> label name
			if row_i < 0:
				labels = line.split('\t')
				label_idx = 0
				for label in labels:
					mapIdxToLabel[label_idx] = label.strip()
					label_idx += 1
				row_i += 1
				continue

			# Split up line into individual positions
			values = line.split('\t')
			col_j = -1
			for value in values:
				if col_j < 0:
					col_j += 1
					continue
				# Look up real coordinates in master matrix (row, col)
				# (index -> label -> index)
				label = mapIdxToLabel.get(row_i)
				masterRow = mapLabelToIdx.get(label)
				label = mapIdxToLabel.get(col_j)
				masterCol = mapLabelToIdx.get(label)
				
				# Calculate new score and update avg matrix
				scoreToAdd = value
				oldScore = avgMatrixScores[masterRow, masterCol]
				oldNumScores = avgMatrixNumScores[masterRow, masterCol]
				newScore = calculateNewScore(oldScore, oldNumScores, scoreToAdd)
				avgMatrixScores[masterRow, masterCol] = newScore
				avgMatrixNumScores[masterRow, masterCol] = oldNumScores + 1
				col_j += 1
			row_i += 1
			
	# ===== WRITE AVERAGE MATRIX TO FILE =====
	# Write matrix out to file in the exact order as master file
	print "Matrix averaging completed. Writing results to file..."
	numpy.savetxt('outputfile.txt', avgMatrixScores, delimiter='\t', newline='\n')
	print "Done."
	

if __name__ == "__main__":
	main(sys.argv[1:])
