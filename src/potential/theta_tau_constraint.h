#ifndef TTC_H
#define TTC_H

#include <string> 
#include <map> 
#include <vector> 
#include <fstream> 
#include <iostream> 
#include <sstream>
#include <Eigen/Dense> 

#include "config.h"
#include "angle.h"
#include "residue.h"
#include "myutils.h"

using namespace std;
using namespace Eigen;

extern map<string, string> TT_list;

//use # or // to exclude the irrelevant
vector<Vector2d> readTTConstraintFiles(const string& name);

double getTTPontial(const vector<Residue>& residuesData, const vector<Vector2d>& thetatau_cons);

#endif