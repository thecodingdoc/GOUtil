//////////////////////////////////////////////////////////////////////
// GOUtil.C  Copyright (c) 2017 Dario Ghersi                        //
// Version: 20171213                                                //
// Goal: Enrichment Analysis tools                                  //    
// Usage: GOUtil -e EDGE_LIST -a ANNOTATIONS -b BACKGROUND          //
//               -t TARGET -o OUTFILE [-u]                          //
//          See User's Guide for more details                       //
//                                                                  //
// This file is part of the GOUtil suite.                           //
// GOUtil is free software: you can redistribute it and/or          //
// modify it under the terms of the GNU General Public License as   //
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

#include <cstring>
#include <iostream>
#include <fstream>
#include <iterator>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <boost/math/distributions/hypergeometric.hpp>

using namespace std;
#include "utilities.h"
#include "enrich.h"

//////////////////////////////////////////////////////////////////////
// DEFINITIONS                                                      //
//////////////////////////////////////////////////////////////////////

void calculateBackgroundFreq(vector<unsigned int> &backgroundFreq,
			     vector<unsigned int> &targetFreq,
			     vector<vector<string> > termCentric,
                             vector<vector<string> > termCentricTarget,
			     Graph goG,
			     unordered_map<unsigned int, set<string> >
			     &withTerm)
{

  // initialize the frequency to 0
  for (unsigned int i = 0; i < backgroundFreq.size(); i++) {
    backgroundFreq[i] = 0;
    targetFreq[i] = 0;
  }
  

  // find all terms of interest
  set<unsigned int> termsOI;
  for (unsigned int i = 0; i < termCentricTarget.size(); i++) {
    
    if (termCentricTarget[i].size() > 0) {

      termsOI.insert(i);
    }
  }

  // find all ancestral terms
  set<unsigned int> ancestors;
  findAllAncestors(termsOI, ancestors, goG);
    
  // calculate the frequency for all terms
  for (set<unsigned int>::iterator it = ancestors.begin();
       it != ancestors.end(); it++) {
      
    // get all children terms
    set<unsigned int>children;
    findChildren(*it, children, goG);

    // calculate the background freq
    set<string>backWithTerm;
    set<string>targetWithTerm;

    for (set<unsigned int>::iterator cit = children.begin();
	 cit != children.end(); cit++) {
      
      copy(termCentric[*cit].begin(), termCentric[*cit].end(),
	   inserter(backWithTerm, backWithTerm.end()));

      copy(termCentricTarget[*cit].begin(), termCentricTarget[*cit].end(),
	   inserter(targetWithTerm, targetWithTerm.end()));

    }
    backgroundFreq[*it] = backWithTerm.size();
    targetFreq[*it] += targetWithTerm.size();
    withTerm[*it] = targetWithTerm;
  }
}

/////////////////////////////////////////////////////////////////////

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
  if (!cmdOptionExists(argv, argv+argc, "-b")) {
    cerr << "Background set file missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-t")) {
    cerr << "Target set file missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-o")) {
    cerr << "Output file missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-p")) {
    cerr << "FDR threshold missing\n";
    err = true;
  }
 
  if (err) {
    cout << USAGE;
    exit(1);
  }
}

//////////////////////////////////////////////////////////////////////

bool comparator2(const intDouble &pair1, const intDouble &pair2)
{
  // comparison function

  return pair1.second < pair2.second;
}

//////////////////////////////////////////////////////////////////////

void doEnrichment(unsigned int targetSize,
		  unsigned int backgroundSize,
                  Graph goG, vector<vector<string> > termCentric,
		  vector<vector<string> > termCentricTarget,
		  EnrichedTerms &enrichTerms) {
  
  // perform enrichment analysis using the hypergeometric test

  
  // compute the background distribution
  vector<unsigned int> backgroundFreq(termCentric.size());
  vector<unsigned int> targetFreq(termCentric.size());
  calculateBackgroundFreq(backgroundFreq, targetFreq, termCentric,
			  termCentricTarget, goG, enrichTerms.withTerm);

  // perform the hypergeometric calculations for each term
  for (unsigned int i = 0; i < targetFreq.size(); i++) {
    
    if (targetFreq[i] > 0) {
      enrichTerms.termIndex.push_back(i);
      boost::math::hypergeometric_distribution<double>
        hgDist(targetSize, backgroundFreq[i], backgroundSize);
      try {
	// calculate the p-value
        enrichTerms.pvalues.push_back(1.0 - boost::math::cdf<double>(hgDist,
      	   			      targetFreq[i] - 1));

	// calculate the enrichment factor (number of actual terms
	// in the target set over number of expected terms)
	double expected = float(targetSize) * backgroundFreq[i] /
	  backgroundSize;
	enrichTerms.enrichFactor.push_back(targetFreq[i] / expected);
      }
      catch (...) {
	enrichTerms.pvalues.push_back(1);
	enrichTerms.enrichFactor.push_back(1);
      }

    }
  } 
}

//////////////////////////////////////////////////////////////////////

void EnrichedTerms::addID(unordered_map<unsigned int, string> revHash)
{

  for (unsigned int i = 0; i < pvalues.size(); i++) {
    termID.push_back(revHash[termIndex[i]]);
  }
}

//////////////////////////////////////////////////////////////////////

void EnrichedTerms::fdrCorrection()
{
  // adjust the p-values by applying the Benjamini-Hochberg correction

  // combine the p-values wih their index for sorting
  vector<intDouble> pPairs;
  intDouble foo;
  for (unsigned int i = 0; i < pvalues.size(); i++) {
    foo = make_pair(i, pvalues[i]);
    pPairs.push_back(foo);
  }

  // sort the p-values
  std::sort(pPairs.begin(), pPairs.end(), comparator2);

  // apply the Benjamini-Hochberg procedure to correct the p-values
  // for multiple hypothesis testing
  double currMin, value;
  vector<double> minValues;
  unsigned int m = pPairs.size();
  for (unsigned int i = 0; i < m; i++) {
    currMin = pPairs[i].second * m / (i + 1);
    for (unsigned int j = i + 1; j < m; j++) {
      value = pPairs[j].second * m / (j + 1);
      if (value < currMin) {
        currMin = value;
      }
    }
    minValues.push_back(min(currMin, 1.0));
  }

  // put the adjusted values in the original p-value order
  adjustedP.resize(pPairs.size());
  sortedOrder.resize(pPairs.size());
  for (unsigned int i = 0; i < pPairs.size(); i++) {
    adjustedP[pPairs[i].first] = minValues[i];
    sortedOrder[i] = pPairs[i].first;
  }
}

//////////////////////////////////////////////////////////////////////

void EnrichedTerms::printResults(string outFileName, double threshold,
				 unordered_map<unsigned int, string>
				 &definition)
{
  // print the results of the enrichment analysis using the
  // following format:
  // TERM_ID TERM_DEFINITION ADJUSTED_PVALUE

  // print the enrichment results
  fstream outFile;
  outFile.open(outFileName, fstream::out);
  outFile << std::scientific;

  // print the results
  for (unsigned int i = 0; i < sortedOrder.size(); i++) {
    if (adjustedP[sortedOrder[i]] <= threshold) {
      outFile << termID[sortedOrder[i]] << "\t" <<
	definition[termIndex[sortedOrder[i]]] << "\t" <<
	adjustedP[sortedOrder[i]] << "\t" <<
	enrichFactor[sortedOrder[i]] << "\t";

      // print the genes contributing to the enrichment
      bool isFirst = true;
      for (set<string>::iterator it =
	     withTerm[termIndex[sortedOrder[i]]].begin();
	   it != withTerm[termIndex[sortedOrder[i]]].end(); it++) {
	if (isFirst) {
	  isFirst = false;
	  outFile << *it;
	}
	else {
	  outFile << " " << *it;
	}
      }
      outFile << endl;
    }
  }
  
  // close the output file
  outFile.close();
}

//////////////////////////////////////////////////////////////////////

void storeSet(set<string> &genes, string fileName)
{

  string line;

  // open the input file
  fstream infile;
  infile.open(fileName, fstream::in);

  // complain if the file doesn't exist
  if (! infile.good()) {
    cerr << "Can't open " << fileName << endl;
    exit(1);
  }

  // process each term
  while (getline(infile, line)) {
    genes.insert(line);
  }

  // close the file
  infile.close();

}

//////////////////////////////////////////////////////////////////////

Parameters::Parameters(char **argv, int argc)
{
  // parse the command-line arguments

  edgesFileName = getCmdOption(argv, argv + argc, "-e");
  annotationsFileName = getCmdOption(argv, argv + argc, "-a");
  backgroundSetFileName = getCmdOption(argv, argv + argc, "-b");
  targetSetFileName = getCmdOption(argv, argv + argc, "-t");
  outFileName = getCmdOption(argv, argv + argc, "-o");
  threshold = stod(getCmdOption(argv, argv + argc, "-p"));
}

//////////////////////////////////////////////////////////////////////

void storeTermCentricAnn(set<string> &background,
			 vector<vector<string> > &termCentric,
			 vector<vector<string> > &termCentricTarget,
			 string annFileName,
                         unordered_map<string, unsigned int> nodeHash,
			 set<string> &target)
{
  // store the annotations in a term-centric fashion for target and
  // background.
  // Filter out non-annotated genes from target set and background set

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
  set<string> allAnnGenes;
  while (getline(termFile, line)) {
    vector<string> tokens;
    istringstream iss(line);

    string gene;
    iss >> gene;

    // add the gene to the set of annotated genes
    allAnnGenes.insert(gene);

    // check if gene is in the background set
    const bool isInBack = background.find(gene) != background.end();

    if (isInBack) {
      // check if gene is in the target set
      const bool isInTarget = target.find(gene) != target.end();

      do {
	string term;
	iss >> term;
	if (term != "") {
          if (isInTarget) {
	    termCentricTarget[nodeHash[term]].push_back(gene);
	  }
	  termCentric[nodeHash[term]].push_back(gene);
	}
      }
      while (iss);
    }
  }

  // insersect the background set with the annotated set
  set<string> newBackground;
  set_intersection(background.begin(), background.end(),
		   allAnnGenes.begin(), allAnnGenes.end(),
		   inserter(newBackground, newBackground.begin()));

  background = newBackground;

  // insersect the target set with the annotated set
  set<string> newTarget;
  set_intersection(target.begin(), target.end(),
		   allAnnGenes.begin(), allAnnGenes.end(),
		   inserter(newTarget, newTarget.begin()));

  target = newTarget;

  // close the file
  termFile.close();
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

  // store the background set
  set<string> backgroundSet;
  storeSet(backgroundSet, p.backgroundSetFileName);

  // store the target set
  set<string> targetSet;
  storeSet(targetSet, p.targetSetFileName);

  // build a hash table with term->index relationship
  unordered_map<string, unsigned int> nodeHash;
  unordered_map<unsigned int, string> revNodeHash;
  unordered_map<unsigned int, string> definition;
  buildHashTable(p.edgesFileName, nodeHash, revNodeHash, definition);
  
  // build the ontology graph
  Graph goG(nodeHash.size());
  buildGraph(goG, nodeHash, p.edgesFileName);

  // store the annotations
  vector<vector<string> > termCentricAnn(nodeHash.size());
  vector<vector<string> > termCentricAnnTarget(nodeHash.size());

  storeTermCentricAnn(backgroundSet, termCentricAnn,
		      termCentricAnnTarget, p.annotationsFileName,
  		      nodeHash, targetSet);

  // sanity checks
  if (backgroundSet.size() < 1) {
    cerr << "The background set has no genes in it...aborting" << endl;
    exit(1);
  }
  if (targetSet.size() < 1) {
    cerr << "The target set has no genes in it...aborting" << endl;
  }
  if (targetSet.size() > backgroundSet.size()) {
    cerr << "More genes in the target than in the background...aborting" << endl;
    exit(1);
  }

  // perform enrichment analysis
  EnrichedTerms enrichTerms;
  doEnrichment(targetSet.size(), backgroundSet.size(), goG,
	       termCentricAnn, termCentricAnnTarget, enrichTerms);
  
  // assign term ID, definitions, and perform FDR correction
  enrichTerms.addID(revNodeHash);
  enrichTerms.fdrCorrection();

  // print the results
  enrichTerms.printResults(p.outFileName, p.threshold, definition);
 
  return 0;
}
