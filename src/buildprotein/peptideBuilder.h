#ifndef PEPTIDEBUILDER_H
#define PEPTIDEBUILDER_H

#include <iostream> 
#include <Eigen/Dense> 
#include <string> 
#include <vector> 

#include "residue.h"
#include "geometry.h"

using namespace std;
using namespace Eigen;

//vector<Residue> initResidueChain(const vector<char>& resnames);

Geo* getGeo(int resid, char resname, int* totalAtoms);

//reconstruct from residue.geo
void reconstructFromGeo(vector<Residue>& residuesData);
void reconstructSideChainFromGeo(vector<Residue>& residuesData);

//
//vector<Geo*> initGeoChain(vector<Residue>& residuesData);
//
//vector<Atom> rebuildSideChain(vector<Geo*>& geosData, vector<Residue>& residuesData);
//
//vector<Atom> outputAtomsData(const vector<Geo*>& geosData);

#endif