#!/usr/bin/python

######################################################################
# mdsSpectralCombo.py                                                #
# Author:  Dario Ghersi                                              #
# Version: 20171221                                                  #
# Usage:   ./spectralClustering.py SEMSIMFILE OUTFILE NUM_CLUSTERS   #
######################################################################

import numpy
import math
from sklearn import manifold
from scipy.spatial.distance import squareform
from sklearn.cluster import SpectralClustering
import sys

######################################################################
# CONSTANTS                                                          #
######################################################################

ALPHA = 1 # parameter for the kernel function

######################################################################
# FUNCTIONS                                                          #
######################################################################

def doMDS(semDistFileName):
  """
  perform MDS on the distance matrix
  """

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
  semDistMat = squareform(semDistList)  
  
  ## perform MDS
  seed = numpy.random.RandomState(seed=1)
  mds = manifold.MDS(n_components=2, random_state=seed,
                     dissimilarity="precomputed", n_jobs=-1)
  pos = mds.fit(semDistMat).embedding_

  return pos, allGO
  
######################################################################

def findMedoids(distMat, labels):
  """
  extract the representative term for each cluster
  """

  clusterIndices = list(set(labels))
  clusterIndices.sort()
  medoids = []

  ## process each cluster
  for index in clusterIndices:
      
    # find the members of the cluster
    membInd = numpy.where([x == index for x in labels])

    # calculate the average distance of a member against all members
    avgDist = []
    for item in membInd[0]:
      dist = 0.0
      for item2 in membInd[0]:
        dist += distMat[item][item2]

      dist /= len(distMat)
      avgDist.append(dist)

    medoids.append(membInd[0][numpy.argmin(avgDist)])
    
  return medoids

######################################################################

def getAffMat(xyCoords):
  """
  compute the affinity matrix for spectral clustering
  """

  ## store the XY
  affMat = []
  distMat = []

  ## compute the distance matrix
  numTerms = len(xyCoords)
  for i in range(numTerms - 1):
    x1 = xyCoords[i][0]
    y1 = xyCoords[i][1]
    
    for j in range(i + 1, numTerms):
      x2 = xyCoords[j][0]
      y2 = xyCoords[j][1]

      euclD = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
      d = math.exp(-ALPHA * euclD)
      distMat.append(euclD)
      affMat.append(d)

  distMat = squareform(distMat)
  affMat = squareform(affMat)
   
  return affMat, distMat

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
# MAIN PROGRAM                                                       #
######################################################################

## parse the parameters
if len(sys.argv) != 4:
  print "Usage: ./mdsSpectralCombo.py SEMSIMFILE OUTFILE NUM_CLUSTERS"
  sys.exit(1)
semSimFileName, outFileName, numCl = sys.argv[1:]
numCl = int(numCl)

## perform MDS
mdsXY, allGO = doMDS(semSimFileName)

## get the distance matrix
affMat, distMat = getAffMat(mdsXY)

## perform spectral clustering on the distance matrix
numpy.random.seed(1)
sc = SpectralClustering(n_clusters=numCl, affinity='precomputed')
sc.fit(affMat)
labels = sc.labels_

## find the medoids
medoids = findMedoids(distMat, labels)

## print the results for index in cluster
outfile = open(outFileName, "w")
for i in range(len(allGO)):
  outfile.write(allGO[i] + "\t" + str(labels[i]) + "\t" +
                str(i in medoids) + "\n")


