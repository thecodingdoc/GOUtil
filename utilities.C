//////////////////////////////////////////////////////////////////////
// utilities.C  Copyright (c) 2017 Dario Ghersi                     //
// Version: 20171212                                                //
// Goal:    utility functions                                       //
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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

#include "utilities.h"

Graph::out_edge_iterator out_begin, out_end;
Graph::in_edge_iterator in_begin, in_end;


//////////////////////////////////////////////////////////////////////
// DEFINITIONS                                                      //
//////////////////////////////////////////////////////////////////////

void buildGraph(Graph &goG, unordered_map<string, unsigned int> nodeHash,
		string fileName)
{

  string line, node1, node2;
  
  // open the input file
  fstream infile;
  infile.open(fileName, fstream::in);

  // complain if the file doesn't exist
  if (! infile.good()) {
    cerr << "Can't open " << fileName << endl;
    exit(1);
  }

  // process each edge
  while (getline(infile, line)) {
    istringstream iss(line);
    iss >> node1; iss >> node2;

    boost::add_edge(nodeHash[node1], nodeHash[node2], goG);
  }

  infile.close(); 
}

//////////////////////////////////////////////////////////////////////

bool cmdOptionExists(char **begin, char **end,
                     const string& option)
{
  return find(begin, end, option) != end;
}

//////////////////////////////////////////////////////////////////////

void findChildren(unsigned int v, set<unsigned int> & children,
		  Graph goG) {
  // find all children of a node

  unsigned int node, neighbor;
  vector<unsigned int> stack;
  stack.push_back(v);
  
  while (stack.size() > 0) {
    node = stack.back();
    stack.pop_back();
    children.insert(node);
    
    for (boost::tie(in_begin, in_end) = in_edges(node, goG);
         in_begin != in_end; ++in_begin) {
      neighbor = source(*in_begin, goG);
      children.insert(neighbor);
      stack.push_back(neighbor);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void findAllAncestors(set<unsigned int> &termsOI,
		      set<unsigned int> &ancestors, Graph goG)
{

  // initialize the stack
  vector<unsigned int> stack;
  for (set<unsigned int>::iterator it = termsOI.begin();
       it != termsOI.end(); it++) {
    stack.push_back(*it);
  }

  // do depth-first search
  unsigned int node, neighbor;
  while (stack.size() > 0) {
    node = stack.back();
    stack.pop_back();
    ancestors.insert(node);
    for (boost::tie(out_begin, out_end) = out_edges(node, goG);
         out_begin != out_end; ++out_begin) {
      neighbor = target(*out_begin, goG);
      ancestors.insert(neighbor);
      stack.push_back(neighbor);
    }
  }
}

//////////////////////////////////////////////////////////////////////

char *getCmdOption(char **begin, char **end,
                   const string & option)
{
  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) {
    return *itr;
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////

void storeTermCentricAnn(set<string> background,
			 vector<vector<string> > &termCentric,
			 vector<vector<string> > &termCentricTarget,
			 string annFileName,
                         unordered_map<string, unsigned int> nodeHash,
			 set<string> target)
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
  while (getline(termFile, line)) {
    vector<string> tokens;
    istringstream iss(line);

    string gene;
    iss >> gene;

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
  
  // close the file
  termFile.close();
}

//////////////////////////////////////////////////////////////////////
