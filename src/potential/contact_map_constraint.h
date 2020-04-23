#ifndef CMC_H
#define CMC_H

#include <string> 
#include <map> 
#include <vector> 
#include <fstream> 
#include <iostream> 
#include <sstream>
#include <Eigen/Dense> 

#include "config.h"
#include "myutils.h"
#include "angle.h"

using namespace std;
using namespace Eigen;

extern map<string, string> CM_list;

struct CMInfo{
	MatrixXd m_cm_bool;
	MatrixXd m_cm_value;
};

//use # or // to exclude the irrelevant, use 0 8
CMInfo readCMConstraintFiles(const string& name, int size);

CMInfo getContactMapInfo(const vector<Residue>& residuesData);
double getCMConsScore(const CMInfo& cm_cons, const CMInfo& cm_predict);

#endif