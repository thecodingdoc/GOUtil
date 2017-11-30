#!/usr/bin/python

#######################################################################
# Requires pronto package                                             #
# pip install pronto should work for most systems.                    #
#                                                                     #
# It parses the go.obo file to output all the relations in a given    #
# namespace. (one of biological_process or molecular_function or      #
# cellular_component).                                                #
# Saves 'is_a', 'part_of', and 'has_part' relations.                  #    
#######################################################################

import pronto
import sys

## read the input parameters
if len(sys.argv) != 4:
    print('Usage: ./buildEdgeList.py go.obo namespace outputfile\n')
    print('Example: ./buildEdgeList.py go.obo biological_process edgeList.txt\n')
    sys.exit(1)

obofile, namespace, output = sys.argv[1:4]

## read the obo file
go = pronto.Ontology(obofile)
edgelist = []
count = 0
## iterate over each term and find its parents and save the pairs
for term in go:
    if term.other['namespace'][0] == namespace:
        for parent in term.parents:
            edge = term.id + " " + parent.id
            if edge not in edgelist: edgelist.append(edge)


## write all the pairs into given output file
with open(output,'w') as fw:
    for e in edgelist: fw.write(e + "\n")

## prints the total number of edges
print(str(len(edgelist))+" edges found.")
