#ifndef GA_H
#define GA_H

#include <string>
#include <vector>
#include <assert.h>
#include <Eigen/Dense> 

#include "myutils.h"
#include "side_chain_modeling.h"
#include "config.h"
#include "pdb.h"
#include "basic.h"
#include "residue.h"
#include "peptideBuilder.h"
#include "CSFScore.h"
#include "theta_tau_constraint.h"
#include "torsion_angles_constraint.h"
#include "contact_map_constraint.h"
#include "LJPotential.h"
#include "other_potentials.h"

using namespace std;
using namespace Eigen;

class GA{
public:

	int index;

	double current_score = 999999999;
	double m_coverage;

	bool global_refine;
	bool use_sidechain;
	bool use_potentials_in_selection;
	bool use_dssp;

	int res_length;
	double ga_steps;
	int ga_steps_counter;

	double T;
	double delta;

	double mutation_ratio;
	double combination_ratio;
	double mutation_times_outer;
	double mutation_times_inner;

	GA(const vector<Residue>& residuesData, double coverage, RandomThreads& random_thread, int index);

	void init(vector<Residue>& residuesData, const vector<vector<TAFeature>>& possiblePhipsis,
		const CMInfo& cm_cons, const vector<map<string, vector<CSFFeature>>>& CSFData, 
		const vector<DSSPInfo>& dssp_cons, const vector<Vector2d>& theta_tau_cons);

	double mutation_op(vector<Residue>& residuesData, const vector<vector<TAFeature>>& possiblePhipsis,
		const CMInfo& cm_cons, const vector<map<string, vector<CSFFeature>>>& CSFData,
		const vector<int>& mutation_indexs, RandomThreads& random_thread, const map<string, vector<RotamerFeature>>& RotamerData,
		const vector<DSSPInfo>& dssp_cons, const vector<Vector2d>& theta_tau_cons);

};

//find the structure with minimum potential
int choose_op(vector<double>& scores);

//combine the structure with minimum potential with the others
void combination_op(const vector<Residue>& residuesData_min, vector<Residue>& residuesData, const vector<int>& combination_indexs);

double getTotalPotentials(const vector<Residue>& residuesData, const vector<vector<TAFeature>>& possiblePhipsis,
	const CMInfo& cm_cons, const vector<map<string, vector<CSFFeature>>>& CSFData, const GA& ga, 
	const vector<DSSPInfo>& dssp_cons, const vector<Vector2d>& theta_tau_cons);

#endif