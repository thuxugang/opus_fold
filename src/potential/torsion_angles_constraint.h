#ifndef TAC_H
#define TAC_H

#include <string> 
#include <map> 
#include <vector> 
#include <fstream> 
#include <iostream> 
#include <sstream>
#include <Eigen/Dense> 

#include "config.h"
#include "basic.h"
#include "residue.h"
#include "myutils.h"

using namespace std;
using namespace Eigen;

extern map<string, string> TA_list;

//use # or // to exclude the irrelevant
vector<Vector2d> readTAConstraintFiles(const string& name);

void addTAConsToPossiblePhipsis(vector<vector<TAFeature>>& possiblePhipsis, const vector<Vector2d>& phipsi_cons);

double getTAPontial(const vector<Residue>& residuesData, const vector<vector<TAFeature>>& possiblePhipsis);

#endif