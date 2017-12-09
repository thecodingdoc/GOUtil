//////////////////////////////////////////////////////////////////////
// GOUtil.h  Copyright (c) 2017 Dario Ghersi                        //
// Version: 20171209                                                //
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


//////////////////////////////////////////////////////////////////////
// CONSTANTS                                                        //
//////////////////////////////////////////////////////////////////////

const string USAGE = "\nUsage:\nGOUtil -e EDGE_LIST -a ANNOTATIONS -b BACKGROUND -t TARGET -o OUTFILE -p FDR_THRESHOLD\n";

//////////////////////////////////////////////////////////////////////
// CLASSES, STRUCTS, AND TYPEDEFS                                   //
//////////////////////////////////////////////////////////////////////

class EnrichedTerms {

 public:
  vector<unsigned int> termIndex;
  vector<string> termID;
  vector<string> definition;
  vector<double> pvalues;
  vector<double> adjustedP;
  vector<double> enrichFactor;
  vector<unsigned int> sortedOrder;
  unordered_map<unsigned int, set<string> > withTerm;

  void addID(unordered_map<unsigned int, string>);
  void fdrCorrection();
  void printResults(string, double,
		    unordered_map<unsigned int, string> &);
};

//////////////////////////////////////////////////////////////////////

class Parameters {

 public:
  char *edgesFileName;
  string annotationsFileName;
  string backgroundSetFileName;
  string targetSetFileName;
  string outFileName;
  double threshold;

  Parameters(char **, int);
};

typedef pair<unsigned int, double> intDouble;



//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

void calculateBackgroundFreq(vector<unsigned int> &,
			     vector<unsigned int> &,
			     vector<vector<string> >,
                             vector<vector<string> >,
			     Graph,
			     unordered_map<unsigned int, set<string> > &);
bool cmdOptionExists(char **, char **, const string &);
char *getCmdOption(char **, char **, const string &);
void doEnrichment(unsigned int, unsigned int, Graph,
		  vector<vector<string> >, vector<vector<string> >,
		  EnrichedTerms &);
void findAllAncestors(set<unsigned int> &, set<unsigned int> &, Graph &);
void storeSet(set<string> &, string);
void storeTermCentricAnn(set<string>,
			 vector<vector<string> > &,
			 vector<vector<string> > &, string,
                         unordered_map<string, unsigned int>,
			 set<string> target);



