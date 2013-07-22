import numpy

def calculateAverage(allDsds, masterMatrix, threshold, arg):
	print "  using threshold module"
	numLabels = allDsds.shape[0]
	matrixAverageScores = numpy.zeros(shape=(numLabels, numLabels))
	for row_i in xrange(0, numLabels):
		for col_j in xrange(0, numLabels):
			if col_j == row_i: # diagonal
				matrixAverageScores[row_i, col_j] = 0.0
			if col_j > row_i: # upper triangle
				dsdScores = filter(None, sorted(allDsds[row_i][col_j][:], reverse=True))
				numScores = len(dsdScores)
				
				if numScores > 0:
					cutoffIdx = int(numScores * float(threshold)) 
					if cutoffIdx - int(cutoffIdx) > 0:
						cutoffIdx += 1 # ceiling
					dsdScores = dsdScores[:cutoffIdx]
					average = sum(dsdScores)/cutoffIdx
					matrixAverageScores[row_i, col_j] = average
					matrixAverageScores[col_j, row_i] = average
				else:
					matrixAverageScores[row_i, col_j] = 999.9
					matrixAverageScores[col_j, row_i] = 999.9

	return matrixAverageScores
