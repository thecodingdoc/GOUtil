#!/usr/bin/python

######################################################################
# spectralClustering.py                                              #
# Author:  Dario Ghersi                                              #
# Version: 20190307                                                  #
# Usage:   ./spectralClustering.py MDS_XY OUTFILE NUM_CLUSTERS       #
# Note:    if NUM_CLUSTERS is -1 (i.e., not specified by the user),  #
#          the script uses the eigengap heuristic to determine       #
#          the "optimal" number of clusters as proposed in           #
#          Luxburg (tutorial on spectral clustering)                 #
######################################################################

import numpy as np
import math
import scipy
from sklearn.cluster import SpectralClustering
import sys

######################################################################
# CONSTANTS                                                          #
######################################################################

ALPHA = 1 # parameter for the kernel function
KMAX = 10 # maximum number of clusters

######################################################################
# FUNCTIONS                                                          #
######################################################################

def eigenGapHeuristic(affMat, maxK=20):
  """
  compute the eigengap heuristic to estimate the "optimal" number of
  clusters, following von Luxburg
  """

  ## compute the laplacian
  L = np.identity(np.size(affMat, 1)) - affMat
  
  ## compute the eigenvalues of the laplacian
  w, v = np.linalg.eig(L)

  ## sort the eigenvalues
  w.sort()

  ## compute the eigengap
  k = min(maxK, len(w))
  if k == 1:
    return 1
  else:
    gaps = np.diff(w[:k])

    maxG = gaps[0]
    numCl = 2
    for i in range(1, k-1):
      if gaps[i] > maxG:
        maxG = gaps[i]
        numCl = i + 2
      
  return numCl
  
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
    membInd = np.where([x == index for x in labels])

    # calculate the average distance of a member against all members
    avgDist = []
    for item in membInd[0]:
      dist = 0.0
      for item2 in membInd[0]:
        dist += distMat[item][item2]

      dist /= len(distMat)
      avgDist.append(dist)

    medoids.append(membInd[0][np.argmin(avgDist)])
    
  return medoids

######################################################################

def getAffMat(mdsXYFileName, nn=5):
  """
  compute a locally adapted affinity matrix for spectral clustering.
  nn is the parameter for the nearest neighbor
  """

  ## store the XY
  xyCoords = []
  allGO = []
  affMat = []
  distMat = []
  mdsXYFile = open(mdsXYFileName, "r")
  for line in mdsXYFile:
    go, x, y = line[:-1].split()
    x = float(x)
    y = float(y)
    allGO.append(go)
    xyCoords.append([x, y])
    
  mdsXYFile.close()

  ## make sure nn is not larger than the number of points
  nn = min(nn, len(xyCoords))

  ## z-score normalize the data
  xyCoords = np.array(xyCoords)
  xyCoords = scipy.stats.zscore(xyCoords, axis=0)

  ## compute the distance matrix
  distMat = scipy.spatial.distance_matrix(xyCoords, xyCoords)

  ## compute the sigmas
  sigmas = []
  for i in range(distMat.shape[0]):
    temp = np.copy(distMat[i])
    temp.sort()
    sigmas.append(np.median(temp[:nn]))
  
  ## compute the sigma matrix
  sigmaMat = np.matrix(sigmas)
  sigmaMat = np.matmul(np.transpose(sigmaMat), sigmaMat)
    
  ## build the kernel function
  affMat = np.exp(-np.square(distMat) / sigmaMat)

  return affMat, distMat, allGO

######################################################################
# MAIN PROGRAM                                                       #
######################################################################

## parse the parameters
if len(sys.argv) != 4:
  print "Usage: ./spectralClustering.py MDS_XY OUTFILE NUM_CLUSTERS"
  sys.exit(1)
mdsXYFileName, outFileName, numCl = sys.argv[1:]
numCl = int(numCl)

## initialize the random number generator


## get the distance matrix
affMat, distMat, allGO = getAffMat(mdsXYFileName, nn=10)

## if the number of clusters is not specified (-1) use the
## eigengap heuristic to predict the "optimal" number of clusters
if numCl < 1:
  numCl = eigenGapHeuristic(affMat)

## perform spectral clustering on the distance matrix
sc = SpectralClustering(n_clusters=numCl, affinity='precomputed',
                        assign_labels="discretize", random_state=1)
sc.fit(affMat)
labels = sc.labels_

## find the medoids
medoids = findMedoids(distMat, labels)

## print the results for index in cluster
outfile = open(outFileName, "w")
for i in range(len(allGO)):
  outfile.write(allGO[i] + "\t" + str(labels[i]) + "\t" +
                str(i in medoids) + "\n")


