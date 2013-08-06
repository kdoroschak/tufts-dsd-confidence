#!/usr/sup/bin/python
"""
mygraph.py -- the program checks whether the graph is connected using BFS
              and calculates the largest connected component

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
"""

import numpy as np
import collections


class Node(object):
    def __init__(self, index, neighbors):
        self.index = index
        # A set of indices (integers, not Nodes)
        self.neighbors = set(neighbors)


class Component(object):
    def __init__(self, nodes):
        self.nodes = nodes
        self.adjacentNodes = set()
        for node in nodes:
            self.adjacentNodes.add(node.index)
            self.adjacentNodes.update(node.neighbors)

    @property
    def size(self):
        return len(self.nodes)

    @property
    def nodeIndices(self):
        return set(node.index for node in self.nodes)

    @property
    def isComplete(self):
        return self.nodeIndices == self.adjacentNodes

    def AddNode(self, node):
        self.nodes.append(node)
        self.adjacentNodes.update(node.neighbors)


def CheckConnect(ppbAdj):

    """
    ppbAdj - adjacency matrix represented as a numpy array

    returns True if the graph is connected, False otherwise
    """

    n = np.size(ppbAdj[0])
    indicators = [True] * n
    garage = {'0': 0}

    garagenodes = garage.keys()
    while(len(garage)):
        a = garagenodes[0]
        indicators[int(a)] = False
        for i in xrange(0, n):
            if ppbAdj[i, int(a)] and indicators[i]:
                garage[str(i)] = i
                indicators[i] = False
        del garage[a]
        garagenodes = garage.keys()
    flag = True
    for i in xrange(0, n):
        if indicators[i]:
            flag = False
    return flag


def CalcLargestComponent(ppbAdj, names):
    """
    ppbAdj - adjacency matrix represented as a numpy array, not connected

    names - node ID stored in a dict

    returns the adjacency matrix and node ID mappings for the largest
            connected component
    """
    n = len(ppbAdj[0])
    nodes = {}
    for index in xrange(0, n):
        neighbors = [neighbor for neighbor,
                     value in enumerate(ppbAdj[index])
                     if value == 1]
        nodes[index] = Node(index, neighbors)
    comps = []  # for all components
    index, node = nodes.popitem()
    comp = Component([node])  # one component
    nComps = 1
    while nodes:
        if not comp.isComplete:
            toBeAddedNodeSet = comp.adjacentNodes.difference(comp.nodeIndices)
            toBeAddedNodeIndex = toBeAddedNodeSet.pop()
            toBeAddedNode = nodes.pop(toBeAddedNodeIndex)
            comp.AddNode(toBeAddedNode)
        else:
            comps.append(comp)
            nComps = nComps + 1
            (index, node) = nodes.popitem()
            comp = Component([node])
    comps.append(comp)
    # LCC = largest connected component
    LCCsize = max([component.size for component in comps])
    #print count
    for index in xrange(0, nComps):
        if len(comps[index].nodeIndices) == LCCsize:
            LCCindex = index
            break
    #print "there are", nComps, "components, and LCC is the", LCCindex
    LCCNodeIndex = comps[LCCindex].nodeIndices
    newAdj = np.zeros((LCCsize, LCCsize))
    newNames = {}
    indicator = [True]*n
    newIndex = 0
    OldNameKeys = names.keys()
    for index in LCCNodeIndex:
        indicator[index] = False
        newNames[OldNameKeys[index]] = newIndex
        newIndex = newIndex + 1
        #print index+1, " is in LCC"
    newNames = collections.OrderedDict(sorted(newNames.items(),
                                              key=lambda x: x[1]))
    for i in xrange(0, n):
        for j in xrange(i+1, n):
            if not indicator[i] and not indicator[j] and ppbAdj[i, j]:
                newI = newNames[OldNameKeys[i]]
                newJ = newNames[OldNameKeys[j]]
                newAdj[newI, newJ] = 1
                newAdj[newJ, newI] = 1
    return(newAdj, newNames, nComps)
