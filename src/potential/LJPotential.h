#ifndef LJPOTENTIAL_H
#define LJPOTENTIAL_H

#include <vector> 
#include <Eigen/Dense> 

#include "residue.h"

using namespace std;
using namespace Eigen;

//main-chain & main-chain
double getMMLJPotential(const vector<Residue>& residuesData);

//main-chain & side-chain
double getMSLJPotential(const vector<Residue>& residuesData);

//side-chain & side-chain
double getSSLJPotential(const vector<Residue>& residuesData);

#endif