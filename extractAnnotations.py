#!/usr/bin/python

#######################################################################
# extractAnnotations.py                                               #
# Author:  Dario Ghersi                                               #
# Version: 20171121                                                   #
# Goal:    parse a Gene Ontology 'gene_association' file and prints   #
#          1. the gene centric annotations                            #
#          2. the term centric annotations                            #
#                                                                     #
# Usage:   ./extractAnnotations.py GENE_ASSOCIATION NAMESPACE         #
#                                  IEA,ND,RCA IDTYPE OUTNAME          #
#          where GENE_ASSOCIATOIN is the annotation file (as          #
#          downloaded from http://www.geneontology.org/),             #
#          NAMESPACE is one of [P, F, C] for biological process,      #
#          molecular function, and cellular component, respectively   #
#          the third argument is the blacklist for the evidence codes #
#          that should be excluded, and OUTNAME is the output name    #
#######################################################################

import sys

#######################################################################
# CONSTANTS                                                           #
#######################################################################

pos = {"id": 1, "symbol": 2, "qualifier": 3, "term": 4, "evidence": 6,
       "namespace": 8, "name": 10}
allowedIdType = ["id", "symbol"]

#######################################################################
# FUNCTIONS                                                           #
#######################################################################

def skipHeader(infile):
  """
  skip the lines beginning with a '!'
  """
  ## find where the header ends
  counter = 0
  while infile.readline().startswith("!"):
    counter += 1

  ## reposition the file iterator
  infile.seek(0)
  for i in range(0, counter):
    infile.readline()

  return infile

#######################################################################

def parseAnnFile(goAnnFileName, namespace, blacklist, idType):
  """
  parse the annotations file and return two dictionaries, one with the
  annotations assigned to the primary IDs and the other with the gene
  names assigned to the primary IDs
  """
  ann = {}
  termCentric = {}
  dictID = {}
  goAnnFile = open(goAnnFileName, "r")
  goAnnFile = skipHeader(goAnnFile) # skip the header
  
  for line in goAnnFile:
    fields = line.split("\t")
    
    if fields[pos["namespace"]] == namespace and not\
       fields[pos["evidence"]] in blacklist and\
       fields[pos["qualifier"]] != "NOT":

      ## gene centric annotations
      if ann.has_key(fields[pos[idType]]):
        if not fields[pos["term"]] in ann[fields[pos[idType]]]:
          ann[fields[pos[idType]]].append(fields[pos["term"]])
      else:
        ann[fields[pos[idType]]] = [fields[pos["term"]]]

      ## term centric annotations
      if termCentric.has_key(fields[pos["term"]]):
        termCentric[fields[pos["term"]]].append(fields[pos[idType]])
      else:
        termCentric[fields[pos["term"]]] = [fields[pos[idType]]]
        
      if not dictID.has_key(fields[pos[idType]]):
        dictID[fields[pos[idType]]] = fields[pos["name"]]

  ## close the file
  goAnnFile.close()

  return [ann, termCentric, dictID]

#######################################################################

def printResults(ann, termCentric, dictID, outName):
  """
  print the gene centric and term centric results
  """
  
  geneCentricFile = open(outName + ".txt", "w")

  ## print the gene centric annotations
  for id in ann:
    printAnn = False
    genesNames = dictID[id].split("|")
    geneCentricFile.write(id + "\t" + "\t".join(ann[id]) + "\n")

  geneCentricFile.close()

  ## print the term centric annotations
  termCentricFile = open("termCentric" + outName.capitalize() +
                         ".txt", "w")
  for term in termCentric:
    termCentricFile.write(term + "\t" + "\t".join(termCentric[term])
                          + "\n")
  
  termCentricFile.close()

#######################################################################
# MAIN PROGRAM                                                        #
#######################################################################

## parse the parameters
if len(sys.argv) < 5:
  print "Usage: ./extractAnnotations.py gene_association namespace IEA,ND,RCA,IPI idType outName"
  sys.exit(1)
goAnnFileName, namespace, blacklist, idType, outName = sys.argv[1:6]


## make sure idType is one of ["id", "symbol"]
if not idType in allowedIdType:
  print "idType can only be one of: [" + ", ".join(allowedIdType) + "]"
  sys.exit(1)

## process the blacklist
blacklist = blacklist.split(",")

## parse the annotations file
ann, termCentric, dictID = parseAnnFile(goAnnFileName, namespace,
                                        blacklist, idType)

## print the results (with the name of the gene)
printResults(ann, termCentric, dictID, outName)
