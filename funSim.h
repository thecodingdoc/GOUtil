//////////////////////////////////////////////////////////////////////
// funSim.h  Copyright (c) 2017 Dario Ghersi                        //
// Version:  20171213                                               //
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
/////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// CONSTANTS                                                        //
//////////////////////////////////////////////////////////////////////

const string USAGE = "\nUsage:\nfunSim -e EDGE_LIST -a ANNOTATIONS -o OUTFILE -t INDEX_TYPE [-f ENRICH_OUTPUT]\n";


//////////////////////////////////////////////////////////////////////
// CLASSES                                                          //
//////////////////////////////////////////////////////////////////////

class Parameters {

 public:
  char *edgesFileName;
  string annotationsFileName;
  string outFileName;
  string enrichFileName;
  string indexType;

  Parameters(char **, int);
};

//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

void calcAllSemSim(string, string, vector<double> &IC, Graph &,
		   unordered_map<unsigned int, string> &);
void calculateFreq(vector<unsigned int> &, vector<vector<string> > &,
		   Graph );
void calculateFreqTarget(vector<unsigned int> &,
			 set<unsigned int> &,
			 vector<vector<string> > &,
			 Graph);
void calculateIC(vector<double> &, vector<unsigned int> &,
		 unsigned int);
void checkCommandLineArgs(char **, int);
void readEnrich(string,  unordered_map<string, unsigned int> &,
		set<unsigned int> &);
double semanticSim(vector<double> &, Graph &, unsigned int ,
		   unsigned int, string);
void storeTermCentricAnn(vector<vector<string> > &, string,
                         unordered_map<string, unsigned int>,
			 unsigned int &);
