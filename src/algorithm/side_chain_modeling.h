#ifndef SIDECHAINMODELING_H
#define SIDECHAINMODELING_H

#include <iostream> 
#include <Eigen/Dense> 
#include <string> 
#include <vector> 

#include "residue.h"
#include "basic.h"
#include "config.h"
#include "DASFScore.h"

using namespace std;
using namespace Eigen;


//init side chain
void initSideChain(vector<Residue>& residuesData, const map<string, vector<RotamerFeature>>& RotamerData,
	const vector<map<string, vector<DASFFeature>>>& DASFData);

//select Rotamer in ga mutation
void selectOptRotamer(vector<Residue>& residuesData, const vector<int>& mutation_indexs,
	const map<string, vector<RotamerFeature>>& RotamerData);

//
//vector<Geo*> initGeoChain(vector<Residue>& residuesData);
//
//vector<Atom> rebuildSideChain(vector<Geo*>& geosData, vector<Residue>& residuesData);
//
//vector<Atom> outputAtomsData(const vector<Geo*>& geosData);

#endif