#!/usr/bin/env python2.7

import os, sys, argparse, numpy, time

def main(argv):
	parser = argparse.ArgumentParser(description='Compare two individual DSDs to illuminate how and where the DSD changed.')
	parser.add_argument('--originalDsd', required=True, type=argparse.FileType('r'), help='DSD to use as a baseline for comparisons with the new DSD.')
	parser.add_argument('--newDsd', required=True, type=argparse.FileType('r'), help='DSD to compare against the original DSD.')
	parser.add_argument('-o', '--outFile', required=True, type=argparse.FileType('w'), help='Output file to store text-based results of comparisons.')
	parser.add_argument('-a', '--adjMatrix', required=False, type=argparse.FileType('r'), nargs=2, help='Adjacency matrix and its specified order to determine neighbors for more detailed comparisons. If this arg is not provided, these comparisons will not be done. Same format as DSD. Specify matrix as first arg and its corresponding ordered list of proteins as the second.')
	parser.add_argument('-p', '--proteins', required=False, type=str, nargs=2, help='Edge eliminated in the new DSD. Will be used for neighborhood calculations. If not provided, these comparisons will not be done.')
	
	args = parser.parse_args()
	
	# Read in DSDs, labels, and adjacency matrix
	if args.adjMatrix is None:
		# Use original matrix as source for global labels
		(originalDsd, originalLabels) = readDsd(args.originalDsd)
		(newDsd, newLabels) = readDsd(args.newDsd, originalLabels)
	else:
		# Use adjacency matrix labels globally
		(adjMatrix, labels) = readAdjMatrix(args.adjMatrix[0], args.adjMatrix[1])
		(originalDsd, originalLabels) = readDsd(args.originalDsd, labels)
		(newDsd, newLabels) = readDsd(args.newDsd, labels)
	
	subDsd = numpy.subtract(newDsd, originalDsd)
	avgDsdValueOriginal = numpy.mean(originalDsd)
	print "Avg DSD value for original DSD: " + str(avgDsdValueOriginal)
	avgDsdValueNew = numpy.mean(newDsd)
	print "Avg DSD value for new DSD: " + str(avgDsdValueNew)
	avgChangeOverall = numpy.mean(subDsd)
	print "Avg change in DSD overall: " + str(avgChangeOverall)		
	numIncrOverall = subDsd[subDsd > 0].size
	print "Number of increasing values: " + str(numIncrOverall)
	numDecrOverall = subDsd[subDsd < 0].size
	print "Number of decreasing values: " + str(numDecrOverall)
	subDsd = ""
	if args.proteins is not None and args.adjMatrix is not None:
		idxProteinA = [i for i,x in enumerate(labels) if x == args.proteins[0]][0]
		idxProteinB = [i for i,x in enumerate(labels) if x == args.proteins[1]][0]
		# Filter adj matrix to show only specified proteins' neighbors
		# Increase adj neighbors to 2, everywhere else is 1 or 0
		adjMatrix[idxProteinA][:] += 1
		adjMatrix[:][idxProteinA] += 1
		adjMatrix[idxProteinB][:] += 1
		adjMatrix[:][idxProteinB] += 1 
		# Row/col intersection was added twice, so subtract 1
		adjMatrix[idxProteinA][idxProteinA] = 0
		adjMatrix[idxProteinB][idxProteinB] = 0
		# Subtract 1 everywhere so adj neighbors are 1, else is < 1.
		adjMatrix = numpy.subtract(adjMatrix, 1)
		# Replace negative values with 0
		adjMatrix = adjMatrix.clip(min=0)
		
		adjOriginalDsd = numpy.multiply(adjMatrix, originalDsd)
		originalDsd = ""
		adjNewDsd = numpy.multiply(adjMatrix, newDsd)
		newDsd = ""
		print ""
		print "NEIGHBORS: "
		subAdjDsd = numpy.subtract(adjNewDsd, adjOriginalDsd)
		avgDsdValueOriginal = numpy.mean(adjOriginalDsd[adjOriginalDsd>0])
		print "Avg DSD value for neighbors in original DSD: " + str(avgDsdValueOriginal)
		avgDsdValueNew = numpy.mean(adjNewDsd[adjNewDsd>0])
		print "Avg DSD value for neighbors in new DSD: " + str(avgDsdValueNew)
		avgChangeOverall = numpy.mean(subAdjDsd[subAdjDsd>0])
		print "Avg change in DSD neighbors overall: " + str(avgChangeOverall)		
		numIncrOverall = subAdjDsd[subAdjDsd > 0].size
		print "Number of increasing values: " + str(numIncrOverall)
		numDecrOverall = subAdjDsd[subAdjDsd < 0].size
		print "Number of decreasing values: " + str(numDecrOverall)
		subAdjDsd = ''
		adjOriginalDsd = ''
		adjNewDsd = ''
	
	newDsd = ""
	originalDsd = ""
	adjMatrix = ""
	args.adjMatrix[0].close()
	args.adjMatrix[1].close()
	args.outFile.close()
	args.newDsd.close()
	args.originalDsd.close()
	# Do general comparisons
		# Subtract original matrix from new matrix
		# Average change over the entire matrix
		# Count number of changes over the entire matrix
		# Count number of increases and decreases
		# Average the increases and decreases
	# Apply mask to filter out all but adjacent neighbors
	# Compare neighbors
		# Subtract original matrix from new matrix
		# Average change over the entire matrix
		# Count number of changes over the entire matrix
		# Count number of increases and decreases
		# Average the increases and decreases
	# Write output to stdout or file

def readDsd(dsdFile, labels):
	# Read in header to get labels
	dsdFile.seek(0)
	localLabels = dsdFile.readline()
	localLabels = localLabels.split('\t')[1:]
	localLabels = [label.strip() for label in localLabels]
	numLocalLabels = len(localLabels)

	if labels is None:
		# Use DSD's labels if none are passed in
		labels = localLabels
		numLabels = numLocalLabels
	else:
		# Use labels from args
		numLabels = len(labels)
		# Check to see whether we have the same order
		if set(localLabels) == set(labels):
			# Wahoo things will be easy!
			rearrange = False
		else:
			rearrange = True
			print "Need to rearrange the DSD. BOOOO! Exiting so Katie can stop being lazy and fix it."
			exit()

	# Read in dsd matrix
	dsdFile.seek(0)
	dsd = numpy.loadtxt(dsdFile, delimiter='\t', usecols=xrange(1,numLabels+1), skiprows=1)
	
	if rearrange:
		# TODO implement rearranging code
		print 'TODO implement DSD rearranging code'
	
	return (dsd, labels)
	
def readAdjMatrix(adjMatrixFile, adjLabelsFile):
	# Process labels
	adjLabels = []
	adjLabelsFile.seek(0)
	for line in adjLabelsFile:
		adjLabels.append(line.strip())
	numLabels = len(adjLabels)
	# Process matrix
	adjMatrix = numpy.zeros((numLabels, numLabels), dtype=int)
	index = 0
	adjMatrixFile.seek(0)
	for line in adjMatrixFile:
		adjMatrix[index] = list(line.strip())
		index += 1
	return (adjMatrix, adjLabels)

if __name__ == "__main__":
	main(sys.argv[1:])
	
