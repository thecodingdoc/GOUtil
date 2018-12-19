# GOUtil v0.3
Copyright (c) 2018  Dario Ghersi, Matt L. Hale, Ishwor Thapa.

## Purpose

The GOUtil suite can be used to perform enrichment analysis,
semantic similarity calculations, and other basic analyses using
Gene Ontology. The suite should work with other controlled
vocabularies or pathway languages (e.g., Reactome).

## Installation

The program uses the C++ Boost library, which should be installed
on your system before attempting to compile the source code.
The Boost library can be found at:

[http://www.boost.org](http://www.boost.org)

To compile the enrichment program with the GCC compiler,
run the following command:

```bash
g++ -O3 -o enrich enrich.C utilities.C --std=gnu++11
```

The semantic similarity program can be compiled by running
the following command:

```bash
g++ -O3 -o funSim funSim.C utilities.C --std=gnu++11
```


## Usage
To perform enrichment analysis calculations, run the `enrich` program as follows:

```bash
## downloading dataset for input to enrich program##
git clone https://github.com/MLHale/GOUtildata.git

./enrich -a GOUtildata/jan-2018/ann.hsa.bp.txt -e GOUtildata/jan-2018/edgeList.bp.txt\
  -t GOUtildata/jan-2018/target.txt -b GOUtildata/jan-2018/background.hsa.bp.txt \
  -o target_output_2018.txt -p 0.05
```

The annotation file ann.hsa.bp.txt contains Gene Ontology annotations
with the following structure:

```text
GENE1 TERM1 TERM2 ...
GENE2 TERM1 TERM3 TERM6 ...
...
```
`edgeList.txt` describes the Gene Ontology graph with a list of
child-parent directed edges and is obtained by running the
`buildEdgeList.py` script on a Gene Ontology obo file.

`target.txt` and background.txt contain the target gene list and
the background gene list, respectively.

`output.txt` is the name of the file where the results will be stored.
The output contains the following fields:

```text
GO_TERM GO_DEFINITION ADJUSTED_PVALUE ENRICHMENT_SCORE GENES
```

The suite also contains a program to calculate the semantic similarity between
pairs of Gene Ontology terms, using either the Lin or the Resnik index.
To run the calculations, use the following command:

```bash
./funSim -a GOUtildata/jan-2018/ann.hsa.bp.txt -e GOUtildata/jan-2018/edgeList.bp.txt\
  -o output.txt -t INDEX_TYPE
```

where INDEX_TYPE is one of [Lin, Resnik]. If the optional `-f ENRICH_OUTPUT` is passed
to the program, then `funSim` only computes the semantic similarity between terms in the
enrichment output.



## Docker
A containerized version of the runtime is also provided using docker.

### Install with docker
- Install [Docker](https://www.docker.com/get-docker)
- Run the following:

```bash
docker-compose build
docker-compose run enrich g++ -O3 -o enrich enrich.C utilities.C --std=gnu++11
```

### Run with Docker
To run with the default parameters:

```bash
docker-compose run enrich
```

to provide alternate parameters:
```bash
docker-compose run enrich ./enrich -a GOUtildata/jan-2018/ann.hsa.bp.txt -e GOUtildata/jan-2018/edgeList.bp.txt\
  -t GOUtildata/jan-2018/target.txt -b GOUtildata/jan-2018/background.hsa.bp.txt \
  -o target_output_2018.txt -p 0.05
```
