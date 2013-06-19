#!/usr/bin/env python2.7

import os, sys, getopt

def main(argv):
	edgeDict = {}
	dsdFiles = []
	dsdPath = ''
	outputFile = ''

	# Read in command line args and process
	try:
		opts, args = getopt.getopt(argv,"hi:o:")
	except getopt.GetoptError:
		print ''
		print 'USAGE: ./DsdAverage <options>'
		print 'Use -h for more information about options.'
	for opt, arg in opts:
		if opt == '-h':
			print 'USAGE: ./DsdAverage <options>'
			print '  -h	Display this help message.'
			print '  -i	Input path of DSD files to be averaged. All files in this directory will be averaged.'
			print '  -o	Name + path of output files.'
			print ''
			print 'EXAMPLE: ./DsdAverage -i ~/Mint/MintNetworks -o ~/Mint/MintNetworks-averageDSD.txt'
			sys.exit(1)
		elif opt == '-i':
			try:
				dsdPath = arg
				dsdFiles = os.listdir(arg)
			except getopt.GetoptError:
				print 'Error opening path to files.'
		elif opt == '-o':
			outputFile = arg

	print "Processing all files at " + dsdPath + "..."	
	# Loop through DSD files and start populating with information.
	for currentDsdFile in dsdFiles:
		with open(dsdPath + "/" + currentDsdFile) as dsdFile:
			for line in dsdFile:
				cols = line.split('\t')
				vertices = cols[0:2]
				vertices.sort()
				vertices = ','.join(vertices)
				vertexValue = edgeDict.get(vertices)
				if (vertexValue != None):
					vertexValue.addDsdScore(cols[2])
				else:
					newEdge = Edge(cols[0:2])
					newEdge.addDsdScore(cols[2].strip())
					edgeDict[vertices] = newEdge

	# Once we've gone through all rows of all files, output the averages
	print "Computing averages and writing to file " + outputFile + "..."
	allValues = edgeDict.values()
	outputFile = open(outputFile, 'w')
	for value in allValues:
		outputFile.write('\t'.join(['\t'.join(value.getVertices()), str(value.getAvgDsd())]) + '\n')
		
	outputFile.close()
	print "Done."

class Edge:
	vertices = []
	dsdScores = []
	notConnectedCount = 0
	scoreCount = 0
	
	def getDsdScores(self):
		return self.dsdScores
	
	def getVertices(self):
		return self.vertices
	
	def getAvgDsd(self):
		if (abs((self.notConnectedCount - len(self.dsdScores))) <= 5):
			print self.vertices 
			print self.notConnectedCount
			print len(self.dsdScores)
			print self.dsdScores

		# Case 0: Self edge.
		#         Return 0. Covers case when score is 'NotConnected'.
		if (self.vertices[0] == self.vertices[1]):
			return 0.0

		# Case 1: All scores good, none indicating 'NotConnected', total# > 0.
		#         Take average as normal.
		elif (len(self.dsdScores) > 0) and (len(self.dsdScores) == self.scoreCount):
			return sum(self.dsdScores) / len(self.dsdScores)
		
		# Case 2: All scores indicate 'NotConnected'.
		#         Return 'NotConnected'.
		elif (self.notConnectedCount == self.scoreCount):
			return 'NotConnected'
		
		# Case 3: Mostly 'NotConnected', some scores too.
		#         Return 'NotConnected'. &&&& Handling??
		elif (self.notConnectedCount > len(self.dsdScores)):
			print 'Mostly n/c:    ' + str(self.notConnectedCount) + ' NC, ' + str(len(self.dsdScores)) + ' scores'
			return 'NotConnected'

		# Case 4: Mostly scores, some 'NotConnected'
		#         Return average of scores present. &&&& Handling??
		else: # TODO &&&& Handle # scores in array <= 0
			print 'Mostly scores: '+ str(self.notConnectedCount) + ' NC, ' + str(len(self.dsdScores)) + ' scores'
			return sum(self.dsdScores) / len(self.dsdScores)

	def addDsdScore(self, score):
		self.scoreCount += 1
		try:
			scoreToAdd = float(score)
			self.dsdScores.append(scoreToAdd)
		except:
			self.notConnectedCount += 1		

	def setVertices(self, vert):
		self.vertices = vert
		self.vertices.sort()

	def __init__(self, vert):
		self.vertices = vert
		self.vertices.sort()
		self.dsdScores = []
		self.notConnectedCount = 0
		self.scoreCount = 0

if __name__ == "__main__":
	main(sys.argv[1:])
