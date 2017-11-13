GOUtil v0.1 -- 11/12/17
Copyright (c) 2017  Dario Ghersi.

PURPOSE

The GOUtil suite can be used to perform enrichment analysis,
semantic similarity calculations, and other basic analyses using
Gene Ontology. The suite should work with other controlled
vocabularies or pathway languages (e.g., Reactome).

---------------------------------------------------------------------

INSTALLATION

The program uses the C++ Boost library, which should be installed
on your system before attempting to compile the source code.
The Boost library can be found at:

www.boost.org

To compile the suite with the GCC compiler,
run the following command:

g++ -O3 -o enrich enrich.C utilities.C --std=gnu++11

---------------------------------------------------------------------

USAGE

To perform enrichment analysis calculations, run the "enrich" program
as follows:

./enrich -a data/annHuman20171106.txt -e edgeList.txt
         -t data/target.txt -b data/background.txt -o output.txt

The annHuman20171106.txt contains Gene Ontology annotations
with the following structure:

GENE1 TERM1 TERM2 ...
GENE2 TERM1 TERM3 TERM6 ...
...

edgeList.txt describes the Gene Ontology graph with a list of
child-parent directed edges.

target.txt and background.txt contain the target gene list and
the background gene list, respectively.

output.txt is the name of the file where the results will be stored.
