#ifndef OP_H
#define OP_H

#include <string> 
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <sstream>
#include <fstream> 
#include <iostream> 
#include <map>

#include "config.h"
#include "myutils.h"
#include "residue.h"
#include "mkdssp.h"

using namespace std;
using namespace Eigen;

/*
0-C
1-S
2-T
3-H
4-G
5-I
6-E
7-B
*/
extern map<string, string> DSSP_list;

class DSSPInfo{
public:
	map<char, double> m_dssp_dict;
	double m_asa;
	char m_dssp;

	DSSPInfo();
	DSSPInfo(char dssp);
	DSSPInfo(double asa, char dssp);
	DSSPInfo(double asa, vector<double>& results);

};

vector<DSSPInfo> readDSSPConstraintFiles(const string& name);

vector<DSSPInfo> getDSSPResults(const vector<Atom>& atomsData, bool use_sidechain);

vector<double> getOtherPotentials(int ga_population);

#endif