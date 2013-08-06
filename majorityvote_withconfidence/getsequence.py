#! /usr/bin/env python2.7 

import numpy as np
import argparse
import sys

def print_stats(np_db):
	print "Mean    : ", str(np_db.mean())
	print "Variance: ", str(np_db.var())
	print "StdDev  : ", str(np_db.std())
	print "Min     : ", str(np_db.min())
	print "Max     : ", str(np_db.max())

parser = argparse.ArgumentParser(description="creates ordered matrix of sequence similarity scores")
parser.add_argument('names_file', type=argparse.FileType('r'))
parser.add_argument('sequence_file', type=argparse.FileType('r'))
parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-f', '--field', help='sequence db field to take, default=10', default=10)
parser.add_argument('-s', dest='stats', help="print statistics", default=False, action="store_true")
options = parser.parse_args()



seq_db = {}
names  = []
for line in options.sequence_file:
	if line[0] == '#': continue
	parts = line.split()
	key   = (parts[0].strip(), parts[1].strip())
	value = parts[int(options.field)].strip()
	if key in seq_db or key[::-1] in seq_db: continue
	seq_db[key] = float(value)
	seq_db[key[::-1]] = float(value)

for line in options.names_file:
	names.append(line.strip())

num_names = len(names)

matrix = np.zeros((num_names, num_names))

options.outfile.write('\t')
for name in names:
	options.outfile.write(name+'\t')
options.outfile.write('\n')

for row in xrange(0, num_names):
	options.outfile.write(names[row]+'\t')
	for col in xrange(0, num_names):
		key    = (names[row], names[col])
		try:
			options.outfile.write(str(seq_db[key]) + '\t')
		except KeyError:
			options.outfile.write(str(matrix[row][col]) + '\t')
	options.outfile.write('\n')

if options.stats:
	print_stats(np.array(seq_db.values()))



#trimat = extract_trimat(matrix, range(0, num_names))
#write_matrix(options.outfile, trimat)

# Extracts triangular matrix from DSD matrix, ordering by indices.



