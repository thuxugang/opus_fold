#include "atom.h"

const double RADII[18] = { 2, 2, 2, 2, 2, 1.8, 2, 2, 1.75, 1.75, 1.75, 1.75, 1.55, 1.55, 1.55, 1.55, 1.9, 1.9 };
const double WELL[18] = { 0.0486, 0.14, 0.0486, 0.1142, 0.1811, 0.08, 0.12, 0.1142, 0.2384, 0.2384, 0.2384, 0.2384, 0.1591, 0.1591, 0.21, 0.1591, 0.16, 0.16 };

void Atom::setLJParams(){
	if (m_name.compare("N") == 0){
		if (m_resname == 'P'){
			this->m_LJType = 12;
		}
		else{
			this->m_LJType = 9;
		}
	}
	else if (m_name.compare("CA") == 0){
		this->m_LJType = 1;
	}
	else if (m_name.compare("C") == 0){
		this->m_LJType = 2;
	}
	else if (m_name.compare("O") == 0){
		this->m_LJType = 13;
	}
	else if (m_resname == 'A'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 5;
		}
	}
	else if (m_resname == 'V'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 3;
		}
		else if (m_name.compare("CG1") == 0){
			this->m_LJType = 5;
		}
		else if (m_name.compare("CG2") == 0){
			this->m_LJType = 5;
		}
	}
	else if (m_resname == 'I'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 3;
		}
		else if (m_name.compare("CG1") == 0){
			this->m_LJType = 5;
		}
		else if (m_name.compare("CG2") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CD1") == 0 || m_name.compare("CD") == 0){
			this->m_LJType = 5;
		}
	}
	else if (m_resname == 'L'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 3;
		}
		else if (m_name.compare("CD1") == 0){
			this->m_LJType = 5;
		}
		else if (m_name.compare("CD2") == 0){
			this->m_LJType = 5;
		}
	}
	else if (m_resname == 'S'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("OG") == 0){
			this->m_LJType = 16;
		}
	}
	else if (m_resname == 'T'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 3;
		}
		else if (m_name.compare("OG1") == 0){
			this->m_LJType = 16;
		}
		else if (m_name.compare("CG2") == 0){
			this->m_LJType = 5;
		}
	}
	else if (m_resname == 'D'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 7;
		}
		else if (m_name.compare("OD1") == 0){
			this->m_LJType = 15;
		}
		else if (m_name.compare("OD2") == 0){
			this->m_LJType = 16;
		}
	}
	else if (m_resname == 'N'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 7;
		}
		else if (m_name.compare("OD1") == 0){
			this->m_LJType = 14;
		}
		else if (m_name.compare("ND2") == 0){
			this->m_LJType = 10;
		}
	}
	else if (m_resname == 'E'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CD") == 0){
			this->m_LJType = 7;
		}
		else if (m_name.compare("OE1") == 0){
			this->m_LJType = 15;
		}
		else if (m_name.compare("OE2") == 0){
			this->m_LJType = 16;
		}
	}
	else if (m_resname == 'Q'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CD") == 0){
			this->m_LJType = 7;
		}
		else if (m_name.compare("OE1") == 0){
			this->m_LJType = 14;
		}
		else if (m_name.compare("NE2") == 0){
			this->m_LJType = 10;
		}
	}
	else if (m_resname == 'K'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CD") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CE") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("NZ") == 0){
			this->m_LJType = 10;
		}
	}
	else if (m_resname == 'R'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CD") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("NE") == 0){
			this->m_LJType = 10;
		}
		else if (m_name.compare("CZ") == 0){
			this->m_LJType = 7;
		}
		else if (m_name.compare("NH1") == 0){
			this->m_LJType = 10;
		}
		else if (m_name.compare("NH2") == 0){
			this->m_LJType = 10;
		}
	}
	else if (m_resname == 'C'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 8;
		}
		else if (m_name.compare("SG") == 0){
			this->m_LJType = 17;
		}
	}
	else if (m_resname == 'M'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("SD") == 0){
			this->m_LJType = 18;
		}
		else if (m_name.compare("CE") == 0){
			this->m_LJType = 5;
		}
	}
	else if (m_resname == 'F'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CD1") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CD2") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CE1") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CE2") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CZ") == 0){
			this->m_LJType = 6;
		}
	}
	else if (m_resname == 'Y'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CD1") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CD2") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CE1") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CE2") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CZ") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("OH") == 0){
			this->m_LJType = 16;
		}
	}
	else if (m_resname == 'W'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CD1") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CD2") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("NE1") == 0){
			this->m_LJType = 11;
		}
		else if (m_name.compare("CE2") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CE3") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CZ2") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CZ3") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CH2") == 0){
			this->m_LJType = 6;
		}
	}
	else if (m_resname == 'H'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 4;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("ND1") == 0){
			this->m_LJType = 11;
		}
		else if (m_name.compare("CD2") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CE1") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("NE2") == 0){
			this->m_LJType = 11;
		}
	}
	else if (m_resname == 'P'){
		if (m_name.compare("CB") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CG") == 0){
			this->m_LJType = 6;
		}
		else if (m_name.compare("CD") == 0){
			this->m_LJType = 6;
		}
	}

	if (this->m_LJType == -1){
		if (m_name.substr(0,1).compare("H") != 0){
			cout << "LJ params not found: " + to_string(m_resid) + " " + m_resname + " " + m_name << endl;
		}
	}
	else{
		this->m_radii = RADII[this->m_LJType - 1];
		this->m_well = WELL[this->m_LJType - 1];
	}

}

Atom::Atom(){}

Atom::Atom(const string& atom_name){
	m_id = -1;
	m_name = atom_name;
	m_restype = 'X';
	m_resname = 'X';
	m_resid = -1;
	m_position = Vector3d(0,0,0);
	if (m_name.compare("N") == 0 || m_name.compare("CA") == 0 || m_name.compare("C") == 0 || m_name.compare("O") == 0){
		m_isMainChain = true;
	}
	else{
		m_isMainChain = false;
	}
	m_LJType = -1;
	m_radii = -1;
	m_well = -1;

}

Atom::Atom(const Atom* other){
	this->m_id = other->m_id;
	this->m_name = other->m_name;
	this->m_restype = other->m_restype;
	this->m_resname = other->m_resname;
	this->m_resid = other->m_resid;
	this->m_position = other->m_position;
	this->m_isMainChain = other->m_isMainChain;
	this->m_LJType = other->m_LJType;
	this->m_radii = other->m_radii;
	this->m_well = other->m_well;
}

Atom::Atom(int atom_id, const string& atom_name, char atom_restype, char atom_resname, int atom_resid, const Vector3d& atom_position){
	m_id = atom_id;
	m_name = atom_name;
	m_restype = atom_restype;
	m_resname = atom_resname;
	m_resid = atom_resid;
	m_position = atom_position;
	if(m_name.compare("N") == 0 || m_name.compare("CA") == 0 || m_name.compare("C") == 0 || m_name.compare("O") == 0
		|| m_name.compare("CB") == 0){
		m_isMainChain = true;
	}else{
		m_isMainChain = false;
	}
	m_LJType = -1;
	m_radii = -1;
	m_well = -1;

	setLJParams();
}


/*
int main(){
	Vector3d v(1,1,1);
	Atom a(10,"CB",'A',10,v);

	cout << a.m_id << endl;
	cout << a.m_name << endl;
	cout << a.m_resname << endl;
	cout << a.m_LJType << endl;
	cout << a.m_radii << endl;
	cout << a.m_well << endl;
	cout << a.m_isMainChain << endl;

	a.setLJParams('M',"SD");
	
	cout << RADII[17] << endl;
	cout << a.m_LJType << endl;
	cout << a.m_radii << endl;
	cout << a.m_well << endl;
	return 0;
}*/


