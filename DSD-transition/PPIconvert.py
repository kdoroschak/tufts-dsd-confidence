#!/usr/sup/bin/python
"""
PPIconvert.py -- this module convert input file to an PPI list file so as
                 to be parsed by DSD program; only limited formats can
                 be converted!

DSD version 0.5, Copyright (C) 2013, Tufts University
@author -- Mengfei Cao, mcao01@cs.tufts.edu
161 College Ave., Medford, MA 02155, USA

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA

usage: PPIconvert.py [-h] [-o OUTFILE] infile

convert PPIs from infile to outfile in the format that DSD can parse

positional arguments:
  infile                read PPIs from infile, either a .csv or .tab file that
                        contains a tab/comma/space delimited table with both
                        IDs at first row and first column, or a .list file
                        that contains for each line one interacting pair

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        converted file name, by default to stdout


"""

import filetype as ft
import argparse
import sys

temp = "convert PPIs from infile to outfile in the format that DSD can parse"
parser = argparse.ArgumentParser(description=temp)

parser.add_argument("infile", help="read PPIs from infile, either "
                    + " a .csv or .tab file that contains a tab/comma/space"
                    + " delimited table with both IDs at first row and"
                    + " first column, or a .list file that contains for"
                    + " each line one interacting pair")
parser.add_argument("-o", "--outfile", help="converted file name,"
                    + " by default to stdout")
#args = ['testfiles//toy.example', '-o', 'haha.test']
#options = parser.parse_args(args)
options = parser.parse_args()

filename = options.infile

(filetype, nSkip, N) = ft.CheckValid(filename)
if filetype == 1:
    (ppbAdj, names) = ft.parseInteractionList(filename, nSkip, N)
elif filetype == 2:
    (ppbAdj, names) = ft.parseFullMatrix(filename, nSkip, N)
elif filetype == 3:
    (ppbAdj, names) = ft.parseColMatrix(filename, nSkip, N)
elif filetype == 4:
    (ppbAdj, names) = ft.parseRowMatrix(filename, nSkip, N)
elif filetype == 5:
    (ppbAdj, names) = ft.parseNoIDmatrix(filename, nSkip, N)

if options.outfile is not None:
    ofile = open(options.outfile, 'w')
else:
    ofile = sys.stdout
if ft.writePPIforDSD(ppbAdj, names.keys(), ofile):
    print "Finished writing!"
#print filetype
#print ppbAdj
#print options.outfile
ofile.close()
