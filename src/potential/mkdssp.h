#ifndef MKDSSP_H
#define MKDSSP_H

#include <vector>
#include <string>

using namespace std;

class OPUSAtom{
public:
	int m_atom_id;
	string m_atom_type;
	string m_res_type;
	int m_res_id;
	double m_x;
	double m_y;
	double m_z;

	OPUSAtom();
	OPUSAtom(int atom_id, string& atom_type, string& res_type, int res_id, double x, double y, double z);
};

class DSSPOriResults{
public:
	char m_dssp8;
	double m_asa;

	DSSPOriResults();
	DSSPOriResults(char dssp8, double asa);
};

vector<char> getDSSP8Results(vector<OPUSAtom>& opusAtoms);

vector<DSSPOriResults> getDSSP8AndASAResults(vector<OPUSAtom>& opusAtoms);

#endif