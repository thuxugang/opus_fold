#ifndef MULT_THREADS_H
#define MULT_THREADS_H

#include <string>
#include <vector>
#include <assert.h>
#include <Eigen/Dense> 
#include <thread>
#include <omp.h>

#include "genetic_algorithm.h"
#include "basic.h"
#include "residue.h"
#include "config.h"
#include "side_chain_modeling.h"
#include "pdb.h"
#include "other_potentials.h"

using namespace std;
using namespace Eigen;

vector<Residue> optimize_mt(const vector<vector<TAFeature>>& possiblePhipsis, const CMInfo& cm_cons,
	const vector<map<string, vector<CSFFeature>>>& CSFData, double coverage, const vector<Atom>& atomsData_cons,
	const string& fasta, const map<string, vector<RotamerFeature>>& RotamerData,
	const vector<map<string, vector<DASFFeature>>>& DASFData, const vector<DSSPInfo>& dssp_cons, 
	const vector<Vector2d>& theta_tau_cons);

#endif