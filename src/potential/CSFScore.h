#ifndef CSFSCORE_H
#define CSFSCORE_H

#include <vector> 
#include <assert.h>
#include <Eigen/Dense> 
#include "residue.h"
#include "basic.h"

using namespace std;
using namespace Eigen;

double getCSFPotential(const vector<Residue>& residuesData, const vector<map<string, vector<CSFFeature>>>& CSFData);

#endif