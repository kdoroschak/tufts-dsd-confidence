import scipy.stats, numpy

def calculateAverage(allDsds, masterMatrix, threshold, arg):
	print "  using the geometric mean module"
	matrixAverageScores = scipy.stats.gmean(allDsds, axis=2, dtype=float)
	matrixAverageScores = matrixAverageScores + matrixAverageScores.T
	matrixAverageScores[matrixAverageScores == 0] = 999.9
	numpy.fill_diagonal(matrixAverageScores, 0)
	
	return matrixAverageScores
