#!/usr/bin/env python2.7

import sys, os, argparse
from random import shuffle
from collections import OrderedDict

'''
    funannotate.py

	given funcat, a list of names, and (optionally) an annotation level and outfile name,
	encodes all MIPS functions in the following form:

	if proteinA is encoded with categories (1,2,5) it will have output
	011001

    Tony Cannistra, June 2013
    Tufts BCB
'''

def main():
	parser = argparse.ArgumentParser(description='produce annotation matrix for yeast genes.')
	parser.add_argument('funcat', metavar="funcat", help="MIPS funcat file")
	parser.add_argument('names', metavar="names", help="File of gene names")
	parser.add_argument("-o", dest="outfile", metavar="outfile", type=argparse.FileType('w'), default=sys.stdout, nargs='?', help="name of output annotation matrix [default: stdout]")
	parser.add_argument("-l", type=int, dest="level", metavar="level", nargs='?', help="MIPS annotation level, default is 1", default=1)
	parser.add_argument("-p", dest="permute", metavar="permute_file", nargs='?', type=argparse.FileType('w'), const=sys.stdout, help="permute option and optional permute_file [if -p is given with no file, default stdout]")
	args = parser.parse_args()

	print "{ Analyzing level %i }" % int(args.level)
	fun_dict, categories = function_db(args.funcat, args.level)
	print "{ Number of unique functional annotations: %i }" % len(categories)
	funs = get_functions(fun_dict, args.names)
	print "{ Writing annotations to %s }" % str(args.outfile.name)
	write_binary(funs, categories, args.outfile)
	if args.permute:
		print '{ Writing annotated indices to %s }' % args.permute.name
		permute_indices(funs, args.permute)

	

def write_binary(ordered_db, sorted_cats, outfile):
	'''
		takes: an OrderedDict database of names:set(annotations)
			   a set of sorted categories
			   an open, writeable Python file object 
		returns:
			   for each gene in order in ordered_db:
			   		if it is annotated, print a '0' followed by 
			   		a binary-encoded list of the annotations in
			   		sorted_cats (see above)
			   		if it is NOT annotated, print a '1'
	'''
	for name, annotation in ordered_db.iteritems():
		line = "1"
		if annotation:
			line = "0"
			for c in sorted_cats:
				if c in annotation:
					line += "1"
				else:
					line += "0"
		outfile.write(line+'\n')

def permute_indices(ordered_db, outfile):
	'''
		takes: an OrderedDict database of names:set(annotations)
			   an open, writeable Python file object
		returns:
			   a random permutation of the indices that 
			   correspond to funcationally-annotated nodes.
	'''
	indices = []
	#nf = open(name_file)
	for index, annotations in enumerate(ordered_db.itervalues()):
		if annotations:
			indices.append(index+1) #not 0 indexed  :/
	shuffle(indices)

	for i in indices:
		outfile.write(str(i)+'\n')


def get_functions(db, name_file):
	'''
		takes: a dictionary of names:annotations
			   a *string* containing the name of a file which has
			     one gene name per line.
		returns:
			   an ordered dictionary of names:annotations for all of the names 
			     in the given file, in the same order as the file.
	'''
	names_annotated = OrderedDict()
	non_annotated = 0
	for name in open(name_file):
		name = name.strip().upper()
		try: 
			names_annotated[name] = sorted(db[name])
		except KeyError:
			names_annotated[name] = set()
			non_annotated += 1
	print "{ %i non-annotated names in %s.}" % (non_annotated, os.path.basename(name_file))
	return names_annotated




def function_db(funcat, level):
	'''
		takes: a string containing the name of the funcat database
		       an integer containing the funcat level of annotation desired
		returns: a TUPLE containing a dictionary of names:annotations for all names in funcat
				 and a sorted set of the unique annotation categories on that level. 
	'''
	fc = open(funcat)
	fun_dict = {}
	tmp_funs  = []
	for line in fc:
		parts = line.strip().split('|')
		if len(parts) != 5: print line
		name = parts[0].upper() #### NOTICE: MADE ALL NAMES UPPERCASE TO MATCH WHAT"S IN MATRIX
		fun_string = parts[1]
		if fun_string.count('.') >= level-1:
			cat_parts = fun_string.split('.')
			cat = ".".join(cat_parts[:level])
			if name in fun_dict:
				fun_dict[name].append(cat)
			else:
				fun_dict[name] = [cat]

	fun_dict = {k:set(v) for k, v in fun_dict.iteritems()} #uniqify, should do this first. 
	
	cats = set()
	for funs in fun_dict.itervalues():
		cats = cats.union(funs)

	return fun_dict, sorted(cats)

## REMEMBER ALL NAMES ARE UPPERCASE


if __name__ == "__main__":
	main()
