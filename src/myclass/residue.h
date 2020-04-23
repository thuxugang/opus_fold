#ifndef RESIDUE_H
#define RESIDUE_H

#include <iostream> 
#include <string> 
#include <vector> 
#include <map>
#include <hash_map>
#include <Eigen/Dense> 

#include "geometry.h"
#include "atom.h"

using namespace std;
using namespace Eigen;

//const string LJ_MAIN[] = { "N", "CA", "C", "O", "CB" };

//class DASFReference{
//public:
//	Vector3d ref;
//	Matrix4d rotation_matrix;
//	DASFReference();
//	Vector3d transCoordinate(Vector3d& position);
//};
//
//
//class Contact{
//public:
//	int m_i;
//	int m_j;
//	double m_bb_distance;
//	Contact();
//	Contact(int i, int j, double bb_distance);
//};
//
//class RotamerLib{
//public:
//	int m_id;
//	double m_prob;
//	double m_RLScores;
//	double m_MS_LJScores;
//	double m_DASFScores;
//	double m_totalScores;
//	vector<double> m_dihedral;
//	vector<Atom> m_atoms;
//	hash_map<int, double> m_rr_sslj;
//
//	RotamerLib(int id, double prob, double x1, double x2, double x3, double x4);
//};


class Residue{
public:
	int m_resid;   
	char m_resname;	

	Geo* m_geo;
	//double m_phi;
	//double m_psi;
	//string m_phipsi;
	//int m_rl_id;

	//vector<RotamerLib> m_rotamerlibs;
	//vector<double> m_dihedrals;
	map<string, Atom> m_atoms;
	//vector<Contact> m_contacts;
	//vector<vector<DASFFeature>> m_standerDASFs;
	//DASFReference m_DASFReference;

	Residue();
	Residue(int param_resid, char param_resname);

	//void getDASFReference();
	//void setDihedrals(vector<double> param_dihedrals);
	//Atom getAtom(const string& atomname);
};

string getTriResname(char resname);
char getResname(const string& resnameTri);

vector<Residue> fastaToResiduesData(const string& fasta);

//PDB atom to residue
vector<Residue> atomToResidue(const vector<Atom>& atomsData, bool main_chain=true);

void getResidueInfo(vector<Residue>& residuesData, bool from_model=false);

//copy information from m_geo to m_atoms
void geoToResidue(vector<Residue>& residuesData);

#endif