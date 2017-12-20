#!/usr/bin/python

######################################################################
# mdsSemSim.py                                                       #
# Author:  Dario Ghersi                                              #
# Version: 20171211                                                  #
# Goal:    The script reads a pairwise semantic similarity file      #
#          between GO terms, turns into a distance matrix, and       #
#          embeds the terms in X,Y coordinates using multidimensional#
#          scaling                                                   #
#                                                                    #
# Usage:   ./mdsSemSim.py SEMSIM_FILE OUTPUT_FILE                    #
# Note:    This script requires a large amount of memory to run.     #
#          The assumption is that the semantic similarity used is    #
#          the Lin index.                                            #
######################################################################

from scipy.spatial.distance import squareform
from sklearn import manifold
import numpy
import sys

######################################################################
# FUNCTIONS                                                          #
######################################################################

def getDistMat(semSimFileName):

  ## store the distances into a dictionary
  semDistList = []
  semSimFile = open(semSimFileName, "r")
  allGO = []
  go1 = ""
  go2 = ""
  previousGO = ""
  for line in semSimFile:
    go1, go2, semSim = line[:-1].split()
    if previousGO != go1:
      allGO.append(go1)
    previousGO = go1

    semDist = 1.0 - float(semSim)
    semDistList.append(semDist)

  semSimFile.close()
  allGO.append(go2)

  ## build a distance matrix
  distMat = squareform(semDistList)

  return distMat, allGO

######################################################################

def printResults(outFileName, allGO, pos):
  """
  print the 2D embedding of the GO terms
  """

  outFile = open(outFileName, "w")
  for i in range(len(pos)):
    outFile.write("%s\t%.3f\t%.3f\n" % (allGO[i], pos[i][0], pos[i][1])) 
  outFile.close()

######################################################################
# MAIN PROGRAM                                                       #
######################################################################

## parse the parameters
if len(sys.argv) != 3:
  print "Usage: msdSemSim.py SEMSIM_FILE OUTPUT_FILE"
  sys.exit(1)
semSimFileName, outFileName = sys.argv[1:]

## build the distance matrix
distMat, allGO = getDistMat(semSimFileName)

## perform MDS on the distance matrix
seed = numpy.random.RandomState(seed=1)
mds = manifold.MDS(n_components=2, random_state=seed,
                   dissimilarity="precomputed", n_jobs=-1)
pos = mds.fit(distMat).embedding_

## print the results
printResults(outFileName, allGO, pos)
