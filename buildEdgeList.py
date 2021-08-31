#!/usr/bin/python

######################################################################
# buildEdgeList.py                                                   # 
# Author:  Dario Ghersi                                              #
# Version: 20210831                                                  #
# Goal:    Parses the go.obo file and outputs all the child-parent   #
#          relationships in a given namespace (one of                #
#          biological_process, molecular_function or                 #
#          cellular_component).                                      #
# N.B.:    Types of relationship include: 'is_a', and 'part_of'      #
#          relationships.                                            #
# Usage:   ./buildEdgeList.py go.obo namespace outputFile            #
######################################################################

import sys

## parse the arguments
if len(sys.argv) != 4:
    print("Usage: ./buildEdgeList.py go.obo namespace outputfile\n")
    print("Example: ./buildEdgeList.py go.obo biological_process edgeList.txt")
    sys.exit(1)
oboFileName, namespace, outputFileName = sys.argv[1:4]

## build the edge list as a dictionary
edgeList = {}
termID = ""
altIDs = []
parents = []
nameDict = {}
with open(oboFileName, "r") as oboFile:
    for line in oboFile:

        ## beginning of a new term
        if line[:6] == "[Term]" or line[:9] == "[Typedef]":

            ## check whether to add a new term
            if termID != "" and ns == namespace:
                edgeList[termID] = []
                if  len(altIDs) > 0:
                    for term in altIDs:
                        edgeList[term] = []
                for term in parents:
                    edgeList[termID].append(term)
                    for altID in altIDs:
                        edgeList[altID].append(term)
            
            ## reset the variables
            termID = ""
            altIDs = []
            ns = ""
            name = ""
            parents = []

        ## term id
        if line[:3] == "id:":
            termID = line[:-1].split("id:")[1].strip()

        ## name
        if line[:5] == "name:":
            name = line[:-1].split("name:")[1].strip()
            nameDict[termID] = name
            
        ## alt id
        if line[:7] == "alt_id:":
            altID = line[:-1].split("alt_id:")[1].strip()
            altIDs.append(altID)
            nameDict[altID] = name

        ## namespace
        if line[:10] == "namespace:":
            ns = line[:-1].split("namespace:")[1].strip()

        ## "is_a" relationship
        if line[:5] == "is_a:":
            parents.append(line.split("is_a:")[1].split("!")[0].strip())

        ## "part_of" relationship
        if line[:21] == "relationship: part_of":
            parents.append(line.split("part_of")[1].split("!")[0].strip())

## add the last term
if termID != "" and ns == namespace:
    edgeList[termID] = []
    if  len(altIDs) > 0:
        for term in altIDs:
            edgeList[term] = []
            for term in parents:
                edgeList[termID].append(term)
                for altID in altIDs:
                    edgeList[altID].append(term)

## write all the edges
with open(outputFileName, "w") as outputFile:
    for child in edgeList:
        for parent in edgeList[child]:
            outputFile.write(child + "\t" + nameDict[child] + "\t" +
                                 parent + "\t" + nameDict[parent] + "\n")
