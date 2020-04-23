#ifndef ISC_H
#define ISC_H

#include <string> 
#include <map> 
#include <vector> 
#include <fstream> 
#include <iostream> 
#include <sstream>
#include <Eigen/Dense> 

#include "basic.h"
#include "residue.h"

using namespace std;
using namespace Eigen;

extern map<string, string> INIT_list;

void addInitConsToPossiblePhipsis(vector<vector<TAFeature>>& possiblePhipsis, const vector<Residue>& residuesData_cons);


#endif