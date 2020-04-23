#ifndef MYUTILS_H
#define MYUTILS_H

#include <string>
#include <vector>
#include <assert.h>
#include <Eigen/Dense> 

#include "config.h"
#include "basic.h"
#include "residue.h"

using namespace std;
using namespace Eigen;

vector<double> splitd(stringstream& ss, const string &s, char delim);
vector<string> splits(stringstream& ss, const string &s, char delim);

class RandomThreads{
public:
	RandomThreads(int seed);
	default_random_engine E;
	//N(0,1)
	normal_distribution<double> normal;
	//uniform
	uniform_int_distribution<int> u_int_res;
	uniform_int_distribution<int> u_int_ta;
	uniform_real_distribution<double> u_real;
	uniform_real_distribution<double> u_real_ta;
	
};

vector<int> getRandomIndexs(int num_total, int num_indexs, RandomThreads& random_thread);

vector<vector<TAFeature>> getPossiblePhiPsis(const vector<Residue>& residuesData, const vector<map<string, vector<TAFeature>>>& TAData);

//use first TA means
void initPhiPsiFromTA(vector<Residue>& residuesData, const vector<vector<TAFeature>>& possiblePhipsis, int init_index);

class OldPhiPsi{
public:
	OldPhiPsi(double phi, double psi, double omega, double ta_possibility);
	double m_phi;
	double m_psi;
	double m_omega;
	double m_ta_possibility;

	vector<double> m_rotamer;
	double m_rotamer_score;
};

void sampleSpecificPhiPsiFromTA(vector<Residue>& residuesData, const vector<vector<TAFeature>>& possiblePhipsis,
	const vector<int>& mutation_indexs, RandomThreads& random_thread);
void sampleOmega(vector<Residue>& residuesData, const vector<int>& mutation_indexs, RandomThreads& random_thread);

vector<OldPhiPsi> saveOldPhiPsi(const vector<Residue>& residuesData, const vector<int>& mutation_indexs, bool use_sidechain);
void restoreOldPhiPsi(const vector<OldPhiPsi>& old_phipsi, vector<Residue>& residuesData, const vector<int>& mutation_indexs, 
	bool use_sidechain);

#endif