#ifndef ATOM_H
#define ATOM_H

#include <iostream> 
#include <string> 
#include <Eigen/Dense> 
#include "atom.h"

using namespace std;
using namespace Eigen;

class Atom{
public:
	int m_id;
	string m_name;
	char m_resname;
	int m_resid;
	char m_restype;
	Vector3d m_position;

	int m_LJType;
	double m_radii;
	double m_well;

	bool m_isMainChain;

	Atom();

	//Geo
	Atom(const string& atom_name);

	//geoToResidue
	Atom(const Atom* other);

	//PDB
	Atom(int atom_id, const string& atom_name, char atom_restype, char atom_resname, int atom_resid, const Vector3d& atom_position);
	
	void setLJParams();
};

#endif