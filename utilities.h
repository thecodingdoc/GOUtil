//////////////////////////////////////////////////////////////////////
// utilities.h  Copyright (c) 2017 Dario Ghersi                     //
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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <unordered_map>

typedef boost::adjacency_list<boost::vecS, boost::vecS,
			      boost::bidirectionalS> Graph;


//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

void buildGraph(Graph &, unordered_map<string, unsigned int>, string);
bool cmdOptionExists(char **, char **, const std::string &);
void findChildren(unsigned int, set<unsigned int> &, Graph);
char *getCmdOption(char **, char **, const std::string &);
unordered_map<string, unsigned int> getNodeHash(set<string>);
void storeTermCentricAnn(set<string>,
			 vector<vector<string> > &,
			 vector<vector<string> > &,
			 string,
                         unordered_map<string, unsigned int>,
			 set<string> target);

