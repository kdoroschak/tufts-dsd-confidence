import numpy

def calculateAverage(allDsds, masterMatrix, threshold, arg):
	print "  using arithmetic mean module"
	#TODO change this to use only top t% of dsd values.
	# Perhaps change lower 1-t% to 0 and then let it filter the rest out
	allDsds.sort(axis=2)
	
	matrixAverageScores = numpy.ma.average(allDsds, axis=2, weights=(allDsds>0))
	matrixAverageScores = matrixAverageScores + matrixAverageScores.T
	matrixAverageScores[matrixAverageScores == 0] = 999.9
	numpy.fill_diagonal(matrixAverageScores, 0)
	
	return matrixAverageScores
