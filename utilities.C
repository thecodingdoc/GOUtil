//////////////////////////////////////////////////////////////////////
// utilities.C  Copyright (c) 2017 Dario Ghersi                     //
// Version: 20171209                                                //
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

void buildGraph(Graph &goG, unordered_map<string, unsigned int> &nodeHash,
		string fileName)
{

  string line, node1, node2, temp;
  
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
    stringstream linestream(line);
    getline(linestream, node1, '\t'); getline(linestream, temp, '\t');
    getline(linestream, node2, '\t'); getline(linestream, temp, '\t');

    boost::add_edge(nodeHash[node1], nodeHash[node2], goG);
  }

  infile.close(); 
}

//////////////////////////////////////////////////////////////////////

void buildHashTable(string fileName,
 		    unordered_map<string, unsigned int> &nodeHash,
		    unordered_map<unsigned int, string> &revNodeHash,
		    unordered_map<unsigned int, string> &definition)
{
  // assign a unique integer to each node name
  
  unsigned int value = 0;

  string line, node1, node2, def1, def2, data;

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
    stringstream linestream(line);
    getline(linestream, node1, '\t'); getline(linestream, def1, '\t');
    getline(linestream, node2, '\t'); getline(linestream, def2, '\t');

    if (!nodeHash.count(node1)) {
      nodeHash[node1] = value;
      revNodeHash[value] = node1;
      definition[value] = def1;
      value++;
    }

    if (!nodeHash.count(node2)) {
      nodeHash[node2] = value;
      revNodeHash[value] = node2;
      definition[value] = def2;
      value++;
    }
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
		  Graph & goG) {
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
		      set<unsigned int> &ancestors, Graph & goG)
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

