import numpy #, scipy.stats

def calculateAverage(allDsds, masterMatrix, threshold, arg):
	print "  using binning module"
	numLabels = allDsds.shape[0]
	bins = xrange(25) # TODO un-hardcode this
	matrixAverageScores = numpy.zeros(shape=(numLabels, numLabels))
	
	for row_i in xrange(0, numLabels):
		for col_j in xrange(0, numLabels):
			if col_j == row_i: # diagonal
				matrixAverageScores[row_i, col_j] = 0.0
			if col_j > row_i and row_i < 20: # upper triangle
				dsdScores = filter(None, sorted(allDsds[row_i][col_j][:], reverse=True))
				hist,binEdges = numpy.histogram(dsdScores, bins=bins)
				maxCountIdx = numpy.where(hist == hist.max())[0] # TODO decide how to break ties
				if len(maxCountIdx) > 1:	
					maxCount = np.zeros(len(maxCountIdx))
					i = 0
					for idx in maxCountIdx:
						maxCount[idx] = hist[idx-1] + hist[idx] + hist[idx+1]
						i += 1
					
					
				upperBinEdge = binEdges[maxCountIdx+1][0] # TODO do NOT use just the first idx
				lowerBinEdge = binEdges[maxCountIdx][0]
				#print upperBinEdge, lowerBinEdge
				dsdScores = [x for x in dsdScores if (lowerBinEdge <= x) and (x < upperBinEdge)]
				#print len(dsdScores)
				average = numpy.mean(dsdScores)
				matrixAverageScores[row_i, col_j] = average
				matrixAverageScores[col_j, row_i] = average	
				
				
	return matrixAverageScores
