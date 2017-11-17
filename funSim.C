//////////////////////////////////////////////////////////////////////
// funSim.C  Copyright (c) 2017 Dario Ghersi                        //
// Version:  20171214                                               //
// Goal:     semantic similarity functions                          //
//                                                                  //
// This file is part of the GOUtil suite.                           //
// GOUtil is free software: you can redistribute it and/or modify   //
// it under the terms of the GNU General Public License as          //
// published by the Free Software Foundation, either version 3 of   //
// the License, or (at your option) any later version.              //
//                                                                  //
// GOUtil is distributed in the hope that it will be useful,        //
// but WITHOUT ANY WARRANTY; without even the implied warranty of   //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    //
// GNU General Public License for more details.                     //
//                                                                  //
// You should have received a copy of the GNU General Public        //
// License along with GOUtil.                                       //
// If not, see <http://www.gnu.org/licenses/>.                      //
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iterator>
#include <math.h>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
#include "utilities.h"
#include "funSim.h"

//////////////////////////////////////////////////////////////////////
// DEFINITIONS                                                      //
//////////////////////////////////////////////////////////////////////

void calculateFreq(vector<unsigned int> &freq,
			     vector<vector<string> > termCentric,
			     Graph goG)
{

  // initialize the frequency to 0
  for (unsigned int i = 0; i < freq.size(); i++) {
    freq[i] = 0;
  }
    
  // calculate the frequency for all terms
  for (unsigned int it = 0; it < freq.size(); it++) {
    
    // get all children terms
    set<unsigned int>children;
    findChildren(it, children, goG);

    // calculate the background freq
    set<string>withTerm;

    for (set<unsigned int>::iterator cit = children.begin();
	 cit != children.end(); cit++) {
      
      copy(termCentric[*cit].begin(), termCentric[*cit].end(),
	   inserter(withTerm, withTerm.end()));
    }
    freq[it] = withTerm.size();
  }
}

//////////////////////////////////////////////////////////////////////

void calculateIC(vector<double> &IC, vector<unsigned int> &freq,
		 unsigned int totSize)
{

  for (unsigned int i = 0; i < freq.size(); i++) {
    if (freq[i] > 0) {
      IC[i] = -log10(float(freq[i]) / totSize);
    }
    else {
      IC[i] = -1.0;
    }
  }
}

//////////////////////////////////////////////////////////////////////

void checkCommandLineArgs(char **argv, int argc)
{
  // check all the parameters have been provided

  bool err = false;

  if (!cmdOptionExists(argv, argv+argc, "-e")) {
    cerr << "Edge list file missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-a")) {
    cerr << "Annotation file missing\n";
    err = true;
  }

  if (!cmdOptionExists(argv, argv+argc, "-o")) {
    cerr << "Output file missing\n";
    err = true;
  }

  if (!cmdOptionExists(argv, argv+argc, "-t")) {
    cerr << "Index type missing\n";
    err = true;
  }
 
  if (err) {
    cout << USAGE;
    exit(1);
  }
}

//////////////////////////////////////////////////////////////////////

void storeTermCentricAnn(vector<vector<string> > &termCentric,
			 string annFileName,
                         unordered_map<string, unsigned int> nodeHash,
			 unsigned int &totSize)
{

  string line;

  // open the input file
  fstream termFile;
  termFile.open(annFileName, fstream::in);

  // complain if the file doesn't exist
  if (! termFile.good()) {
    cerr << "Can't open " << annFileName << endl;
    exit(1);
  }
  
  // process each gene
  set <string> allGenes;
  while (getline(termFile, line)) {

    vector<string> tokens;
    istringstream iss(line);

    string gene;
    iss >> gene;

    allGenes.insert(gene);

    do {
      string term;
      iss >> term;
      if (term != "") {
        termCentric[nodeHash[term]].push_back(gene);
      }
    }
    while (iss);
  }

  totSize = allGenes.size();
  
  // close the file
  termFile.close();
}

//////////////////////////////////////////////////////////////////////

double semanticSim(vector<double> &IC, Graph &goG, unsigned int i,
		   unsigned int j, string indexType)
{

  double semanticSim = -1;

  // deal with identical terms
  if (i == j) {
    if (indexType == "Lin") {
      return(1.0);
    }
    else if (indexType == "Resnik") {
      return(IC[i]);
    }
  }	

  // find the common ancestors of term i and j
  set<unsigned int> setA; set<unsigned int> setB;
  setA.insert(i); setB.insert(j);
  set<unsigned int> ancA; set<unsigned int> ancB;
  findAllAncestors(setA, ancA, goG); findAllAncestors(setB, ancB, goG);
  set<unsigned int> commonAnc;
  set_intersection(ancA.begin(), ancA.end(), ancB.begin(), ancB.end(),
                   inserter(commonAnc, commonAnc.begin()));

  // find the IC of the Most Informative Common Ancestor
  double ICMICA = -1.0;
  for (set<unsigned int>::iterator it = commonAnc.begin();
       it != commonAnc.end(); it++) {
    if (ICMICA < IC[*it]) {
      ICMICA = IC[*it];
    }
  }

  // return the semantic similarity
  if (indexType == "Resnik") {
    semanticSim = ICMICA;
  }
  else if (indexType == "Lin") {
    if (IC[i] > 0 && IC[j] > 0) {
      semanticSim = 2 * ICMICA / (IC[i] + IC[j]);
    }
    else {
      return -1;
    }
  }
  
  return semanticSim;
}

//////////////////////////////////////////////////////////////////////

Parameters::Parameters(char **argv, int argc)
{
  // parse the command-line arguments

  edgesFileName = getCmdOption(argv, argv + argc, "-e");
  annotationsFileName = getCmdOption(argv, argv + argc, "-a");
  outFileName = getCmdOption(argv, argv + argc, "-o");
  indexType = getCmdOption(argv, argv + argc, "-t");

  if (indexType != "Resnik" && indexType != "Lin") {
    cerr << "The index type must be one of [Resnik, Lin]" << endl;
    exit(1);
  }
}

//////////////////////////////////////////////////////////////////////
// MAIN PROGRAM                                                     //
//////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // check the command-line arguments
  checkCommandLineArgs(argv, argc);

  // get the parameters
  Parameters p(argv, argc);
  
  // build a hash table with term->index relationship
  unordered_map<string, unsigned int> nodeHash;
  unordered_map<unsigned int, string> revNodeHash;
  buildHashTable(p.edgesFileName, nodeHash, revNodeHash);

  // build the ontology graph
  Graph goG(nodeHash.size());
  buildGraph(goG, nodeHash, p.edgesFileName);

  // store the annotations
  vector<vector<string> > termCentricAnn(nodeHash.size());
  unsigned int totSize;
  storeTermCentricAnn(termCentricAnn, p.annotationsFileName, nodeHash,
		      totSize);

  // compute the term frequency
  vector<unsigned int> freq(termCentricAnn.size());
  calculateFreq(freq, termCentricAnn, goG);

  // compute the Information Content (IC) of each term
  vector<double> IC(freq.size());
  calculateIC(IC, freq, totSize);

  // calculate and print the semantic similarity
  // between pairs of terms
  fstream outFile;
  outFile.open(p.outFileName, fstream::out);
  outFile << std::scientific;
  for (unsigned int i = 0; i < IC.size() - 1; i++) {
    cout << i << endl;
    for (unsigned int j = i + 1; j < IC.size(); j++) {
      double semSim = semanticSim(IC, goG, i, j, p.indexType);
      if (semSim > 0) {
        outFile << revNodeHash[i] << "\t" << revNodeHash[j] << "\t" <<
	  semSim << endl;
      }
    }
  }

  outFile.close();
  
  return 0;
}