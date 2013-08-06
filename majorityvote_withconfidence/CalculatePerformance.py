#!/usr/bin/env python2.7
"""
@author: mcao01

Usage:
    python CalculatePerformance.py $1 $2
    $1 -- label file
    $2 -- prediction file
"""


import re
import numpy as np
import collections
#from optparse import OptionParser
import argparse
import myparser
import sys


if len(sys.argv) != 3:
    print >> sys.stderr, '\nUsage:'
    print >> sys.stderr, 'python CalculatePerformance.py $1 $2'
    print >> sys.stderr, '$1 -- label file'
    print >> sys.stderr, '$2 -- prediction file\n'
    exit(0)

ppbLabel = myparser.parseLabel(sys.argv[1])
m = len(ppbLabel[0,:]) - 1

prediction = myparser.parsePrediction(sys.argv[2])
m1 = len(prediction[:,0])
##m1 = len(prediction[:,0])

#### calculate accuracy
count = 0.0
for i in xrange(0, m1):
    if ppbLabel[int(prediction[i, 0]), int(prediction[i, 1])+1]:
        count = count + 1
accuracy = count/float(m1)
print 'The accuracy is ', accuracy

#### calculate precision and recall, and F1 score
correct = 0
predict = 0
allfun = 0
for i in xrange(0, m1):
    allfun += sum(ppbLabel[int(prediction[i, 0]), 1:m+1])
    if prediction[i, 2] > 10e-6:
        predict += 1
        if ppbLabel[int(prediction[i, 0]), int(prediction[i, 1])+1]:
            correct += 1
    if prediction[i, 4] > 10e-6:
        predict += 1
        if ppbLabel[int(prediction[i, 0]), int(prediction[i, 3])+1]:
            correct += 1
    if prediction[i, 6] > 10e-6:
        predict += 1
        if ppbLabel[int(prediction[i, 0]), int(prediction[i, 5])+1]:
            correct += 1
precision = float(correct)/float(predict)
recall = float(correct)/float(allfun)
F1 = 2*precision*recall/(precision+recall)
print 'The precision is: ', precision,' the recall is: ', recall
print 'The F1-score is: ', F1, '\n'
