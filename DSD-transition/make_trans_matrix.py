#! /usr/bin/env python2.7
import sys
import numpy as np
import re
import collections
import PPIparser
import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument("physical", help="physical PPI file")
parser.add_argument("-g", "--genetic", help="genetic PPI file")
parser.add_argument("-o", "--outfile", help="PPI outfile name",
                    type=argparse.FileType('w'),
                    default=sys.stdout)
#parser.add_argument("P", help="weight of physical interactions")
#parser.add_argument("G", help="weight of genetic interactions")
options = parser.parse_args()

P_PROB = 0.6 #probability of physical edge in coin flip

def main():
    sys.stderr.write("Load physicals...")
    p_map = load_physicals(options.physical)
    sys.stderr.write("done.\n")
    if (options.genetic):
        sys.stderr.write("Load genetics...")
        g_map = load_genetics(options.genetic)
        sys.stderr.write("done.\n")
    else:
        g_map = {}
    names = list(set(p_map.keys() + g_map.keys()))
    sys.stderr.write("Build transition matrix...")
    tmat  = build_tmat(p_map, g_map, names)
    sys.stderr.write("done.\n")
    sys.stderr.write("Write ppi list...",)
    write_PPI(tmat, names, options.outfile)
    sys.stderr.write("done.\n")

def load_physicals(p_file):
    i_map = {}
    f_open = open(p_file)
    for line in f_open:
        parts = line.split()
        name1 = parts[0]
        name2 = parts[1]
        if name1 in i_map:
            i_map[name1].append(name2)
        else:
            i_map[name1] = [name2]
        if name2 in i_map:
            i_map[name2].append(name1)
        else:
            i_map[name2] = [name1]
    return i_map

def load_genetics(g_file):
    i_map = {}
    f_open = open(g_file)
    for line in f_open:
        parts = line.split()
        name1 = parts[0]
        name2 = parts[1]
        value = abs(float(parts[2]))
        
        if name1 in i_map:
            i_map[name1].append((name2, value))
        else:
            i_map[name1] = [(name2, value)]
        if name2 in i_map:
            i_map[name2].append((name1, value))
        else:
            i_map[name2] = [(name1, value)]

    return i_map

def build_tmat(p_map, g_map, names):
    N = len(names) # number of nodes
    tmat = np.zeros((N, N))

    for index, n in enumerate(names):
        if n in p_map:
            phys_degree = len(p_map[n])
            if n in g_map:
                g_sum = sum([x[1] for x in g_map[n]])
                coin = flip(P_PROB)
                if coin == 'P':
                    for neighbor in p_map[n]:
                        tmat[index, names.index(neighbor)] = 1/float(phys_degree)
                else: #if genetic
                    for n_name, n_value in g_map[n]:
                        tmat[index, names.index(n_name)] = n_value/float(g_sum)
            else: #only physical
                for neighbor in p_map[n]:
                    tmat[index, names.index(neighbor)] = 1/float(phys_degree)
        else: #only genetic
            for n_name, n_value in g_map[n]:
                tmat[index, names.index(n_name)] = n_value/float(g_sum)

    return tmat

def write_PPI(tmat, names, outfile):
    for i,n_i in enumerate(names):
        for j,n_j in enumerate(names):
            if tmat[i,j] != 0:
                outfile.write(n_i + "\t" + n_j + "\t" + str(tmat[i, j]) + '\n')

def flip(p):
    return "P" if random.random() < p else "G"

if __name__ == "__main__":
    main()
