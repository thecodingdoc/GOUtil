# GOUtil v0.1 -- 11/12/17
Copyright (c) 2017  Dario Ghersi.

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

To compile the suite with the GCC compiler,
run the following command:

```bash
g++ -O3 -o enrich enrich.C utilities.C --std=gnu++11
```

## Usage
To perform enrichment analysis calculations, run the `enrich` program as follows:

```bash
./enrich -a data/annHuman20171106.txt -e data/edgeList.txt -t data/target.txt -b data/background.txt -o output.txt
```

The annHuman20171106.txt contains Gene Ontology annotations
with the following structure:

```text
GENE1 TERM1 TERM2 ...
GENE2 TERM1 TERM3 TERM6 ...
...
```
`edgeList.txt` describes the Gene Ontology graph with a list of
child-parent directed edges.

`target.txt` and background.txt contain the target gene list and
the background gene list, respectively.

`output.txt` is the name of the file where the results will be stored.
The output contains the following fields:

```text
GO_TERM ADJUSTED_PVALUE ENRICHMENT_SCORE
```

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
To run with the default parameters

```bash
docker-compose run enrich ./enrich -a data/annHuman20171106.txt -e data/edgeList.txt -t data/target.txt -b data/background.txt -o output.txt
```
