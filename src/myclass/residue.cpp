#include "angle.h"
#include "peptideBuilder.h"

#include "residue.h"

//DASFReference::DASFReference(){
//}
//
//
//
//Vector3d DASFReference::transCoordinate(Vector3d& position){
//	Vector3d position_new = position - ref;
//	Vector4d position_new4d(position_new[0], position_new[1], position_new[2], 1);
//	Vector4d r = position_new4d.transpose()*rotation_matrix;
//	return Vector3d(r[0],r[1],r[2]);
//}
//
//Contact::Contact(int i, int j, double bb_distance){
//	m_i = i;
//	m_j = j;
//	m_bb_distance = bb_distance;
//}
//
//RotamerLib::RotamerLib(int id, double prob, double x1, double x2, double x3, double x4){
//	m_id = id;
//	m_prob = prob;
//	m_dihedral.push_back(x1);
//	m_dihedral.push_back(x2);
//	m_dihedral.push_back(x3);
//	m_dihedral.push_back(x4);
//}
//
//

//======class======
Residue::Residue(){
	//m_rl_id = 0;
}

Residue::Residue(int param_resid, char param_resname){
	//m_rl_id = 0;
	m_resid = param_resid;
	m_resname = param_resname;
}
//======class======

//======getResname======
map<string, char> getResnameMap(){
	map<string, char> resname_map;
	resname_map["GLY"] = 'G';
	resname_map["ALA"] = 'A';
	resname_map["SER"] = 'S';
	resname_map["CYS"] = 'C';
	resname_map["VAL"] = 'V';
	resname_map["ILE"] = 'I';
	resname_map["LEU"] = 'L';
	resname_map["THR"] = 'T';
	resname_map["ARG"] = 'R';
	resname_map["LYS"] = 'K';

	resname_map["ASP"] = 'D';
	resname_map["GLU"] = 'E';
	resname_map["ASN"] = 'N';
	resname_map["GLN"] = 'Q';
	resname_map["MET"] = 'M';
	resname_map["HIS"] = 'H';
	resname_map["PRO"] = 'P';
	resname_map["PHE"] = 'F';
	resname_map["TYR"] = 'Y';
	resname_map["TRP"] = 'W';
	
	return resname_map;
}
const map<string, char> resname_map = getResnameMap();

map<char, string> getTriResnameMap(){

	map<char, string> triresname_map;
	triresname_map['G'] = "GLY";
	triresname_map['A'] = "ALA";
	triresname_map['S'] = "SER";
	triresname_map['C'] = "CYS";
	triresname_map['V'] = "VAL";
	triresname_map['I'] = "ILE";
	triresname_map['L'] = "LEU";
	triresname_map['T'] = "THR";
	triresname_map['R'] = "ARG";
	triresname_map['K'] = "LYS";

	triresname_map['D'] = "ASP";
	triresname_map['E'] = "GLU";
	triresname_map['N'] = "ASN";
	triresname_map['Q'] = "GLN";
	triresname_map['M'] = "MET";
	triresname_map['H'] = "HIS";
	triresname_map['P'] = "PRO";
	triresname_map['F'] = "PHE";
	triresname_map['Y'] = "TYR";
	triresname_map['W'] = "TRP";

	return triresname_map;
}
const map<char, string> triresname_map = getTriResnameMap();

string getTriResname(char resname){

	try{
		return triresname_map.at(resname);
	}catch (...){
		cout << "Uncommon resname: " << resname << endl;
		throw "residue.getTriResname() wrong!";
	}
}

char getResname(const string& resnameTri){

	try{
		return resname_map.at(resnameTri);
	}catch (...){
		cout << "Uncommon resname: " << resnameTri << endl;
		throw "residue.getResname() wrong!";
	}
}
//======getResname======

//======fastaToResiduesData======
vector<Residue> fastaToResiduesData(const string& fasta){
	vector<Residue> residuesData;
	for (int i = 0; i < fasta.size(); ++i){
		residuesData.emplace_back(i + 1, fasta[i]);
	}
	return residuesData;
}
//======fastaToResiduesData======

//void Residue::getDASFReference(){
//	Vector3d ca_ref = m_atoms["CA"].m_position;
//	Vector3d c_ref = m_atoms["C"].m_position;
//	Vector3d o_ref = m_atoms["O"].m_position;
//
//	m_DASFReference.ref = ca_ref;
//	Vector3d c_ref_new = c_ref - m_DASFReference.ref;
//	Vector3d o_ref_new = o_ref - m_DASFReference.ref;
//
//	//c - ca
//	Vector3d x_axis = c_ref_new / c_ref_new.norm();
//	Vector3d c_o = o_ref_new - c_ref_new;
//
//	//o - c perpendicular to x_axis
//	Vector3d y_axis = c_o - (x_axis.dot(c_o) / x_axis.dot(x_axis) * x_axis);
//	y_axis = y_axis / y_axis.norm();
//
//	Vector3d z_axis = x_axis.cross(y_axis);
//
//	m_DASFReference.rotation_matrix << x_axis[0], y_axis[0], z_axis[0], 0, 
//									   x_axis[1], y_axis[1], z_axis[1], 0,
//									   x_axis[2], y_axis[2], z_axis[2], 0,
//									   0, 0, 0, 1;
//
//}
//
//Atom Residue::getAtom(const string& atomname){
//	if(m_atoms.count(atomname)>0){
//		return m_atoms[atomname];
//	}else{
//		throw atomname + "lost at: " + to_string(m_resid) + " " + m_resname;
//	}
//}

//void Residue::setDihedrals(vector<double> param_dihedrals){
//	m_dihedrals = param_dihedrals;
//}

//string getTriResname(char resname){
//	string resnameTri;
//	switch(resname){
//		case('G'):
//			resnameTri = "GLY";
//			break;
//		case('A'):
//			resnameTri = "ALA";
//			break;
//		case('S'):
//			resnameTri = "SER";
//			break;
//		case('C'):
//			resnameTri = "CYS";
//			break;
//		case('V'):
//			resnameTri = "VAL";
//			break;
//		case('I'):
//			resnameTri = "ILE";
//			break;
//		case('L'):
//			resnameTri = "LEU";
//			break;
//		case('T'):
//			resnameTri = "THR";
//			break;
//		case('R'):
//			resnameTri = "ARG";
//			break;
//		case('K'):
//			resnameTri = "LYS";
//			break;
//		case('D'):
//			resnameTri = "ASP";
//			break;
//		case('E'):
//			resnameTri = "GLU";
//			break;
//		case('N'):
//			resnameTri = "ASN";
//			break;
//		case('Q'):
//			resnameTri = "GLN";
//			break;
//		case('M'):
//			resnameTri = "MET";
//			break;
//		case('H'):
//			resnameTri = "HIS";
//			break;
//		case('P'):
//			resnameTri = "PRO";
//			break;
//		case('F'):
//			resnameTri = "PHE";
//			break;		
//		case('Y'):
//			resnameTri = "TYR";
//			break;
//		case('W'):
//			resnameTri = "TRP";
//			break;
//		default:
//			cout << resname << endl;
//			throw "residue.getTriResname() wrong";
//	}
//	return resnameTri;
//}


//======atomToResidue======
void addAtom(vector<Residue>& residuesData, const Atom& atom){

	int current_length = residuesData.size();
	int resid = atom.m_resid;
	char resname = atom.m_resname;

	//first residue
	if(current_length == 0){
		Residue r(resid, resname);
		r.m_atoms.insert(pair<string, Atom>(atom.m_name,atom));
		residuesData.push_back(r);
	}else{
		if(resid == residuesData[current_length-1].m_resid){
			residuesData[current_length-1].m_atoms.insert(pair<string, Atom>(atom.m_name,atom));
		}else{
			Residue r(resid, resname);
			r.m_atoms.insert(pair<string, Atom>(atom.m_name,atom));
			residuesData.push_back(r);		
		}
	} 
}

vector<Residue> atomToResidue(const vector<Atom>& atomsData, bool main_chain){
	vector<Residue> residuesData;
	int length = atomsData.size();
	for(int i=0; i<length; ++i){
		if(main_chain && !atomsData[i].m_isMainChain){
			continue;
		}
		if (atomsData[i].m_LJType == -1 || atomsData[i].m_restype != 'A'){
			cout << "Uncommon atom at: " + to_string(atomsData[i].m_resid) + " " + atomsData[i].m_resname + " " + atomsData[i].m_name  << endl;
			continue;
		}
		addAtom(residuesData, atomsData[i]);
	}
	return residuesData;

}
//======atomToResidue======

//======geoToResidue======
void geoToResidue(vector<Residue>& residuesData){
	//copy information from m_geo to m_atoms
	map<string, Atom*>::iterator iter;
	for (Residue& residue : residuesData){
		for (iter = residue.m_geo->m_geo_atoms.begin(); iter != residue.m_geo->m_geo_atoms.end(); ++iter) {
			residue.m_atoms[iter->first] = Atom(iter->second);
		}
	}
}
//======geoToResidue======



//======getResidueInfo======
void getInitAngleBond(vector<Residue>& residuesData){
	int totalAtoms = 1;
	for (Residue& residue : residuesData){
		residue.m_geo = getGeo(residue.m_resid, residue.m_resname, &totalAtoms);
		//cout << ((GeoContainCB*)residue.m_geo)->CB.m_resname << endl;
	}
}

void getTorsionAngles(vector<Residue>& residuesData){
	int length = residuesData.size();
	for(int i = 0; i < length; ++i){
		if(i == 0){
			residuesData[i].m_geo->phi = -60;
		}else{
			residuesData[i].m_geo->phi = calDihedral(residuesData[i - 1].m_atoms.at("C").m_position, residuesData[i].m_atoms.at("N").m_position,
				residuesData[i].m_atoms.at("CA").m_position, residuesData[i].m_atoms.at("C").m_position);
			residuesData[i].m_geo->omega = calDihedral(residuesData[i - 1].m_atoms.at("CA").m_position, residuesData[i - 1].m_atoms.at("C").m_position,
				residuesData[i].m_atoms.at("N").m_position, residuesData[i].m_atoms.at("CA").m_position);
		}
		if(i == length-1){
			residuesData[i].m_geo->psi = 60;
		}else{
			residuesData[i].m_geo->psi = calDihedral(residuesData[i].m_atoms.at("N").m_position, residuesData[i].m_atoms.at("CA").m_position,
				residuesData[i].m_atoms.at("C").m_position, residuesData[i + 1].m_atoms.at("N").m_position);
		}
	}
}

//missing residue or main-chain atoms is not allowed
void getResidueInfo(vector<Residue>& residuesData, bool from_model){

	//init bond&angle
	getInitAngleBond(residuesData);

	if (from_model){
		//get model's bond&angle
		for (int i = 0; i < residuesData.size(); ++i){
			//first
			if (i == 0){
				residuesData[i].m_geo->CA_N_length = getBondLength(residuesData[i].m_atoms.at("CA").m_position,
					residuesData[i].m_atoms.at("N").m_position);
				residuesData[i].m_geo->CA_C_length = getBondLength(residuesData[i].m_atoms.at("CA").m_position,
					residuesData[i].m_atoms.at("C").m_position);
				residuesData[i].m_geo->N_CA_C_angle = getBondAngle(residuesData[i].m_atoms.at("N").m_position,
					residuesData[i].m_atoms.at("CA").m_position, residuesData[i].m_atoms.at("C").m_position);
			}else{
				//N
				residuesData[i].m_geo->_peptide_bond = getBondLength(residuesData[i - 1].m_atoms.at("C").m_position, 
					residuesData[i].m_atoms.at("N").m_position);
				residuesData[i].m_geo->_CA_C_N_angle = getBondAngle(residuesData[i - 1].m_atoms.at("CA").m_position, 
					residuesData[i - 1].m_atoms.at("C").m_position, residuesData[i].m_atoms.at("N").m_position);
				//CA
				residuesData[i].m_geo->CA_N_length = getBondLength(residuesData[i].m_atoms.at("CA").m_position, 
					residuesData[i].m_atoms.at("N").m_position);
				residuesData[i].m_geo->_C_N_CA_angle = getBondAngle(residuesData[i - 1].m_atoms.at("C").m_position, 
					residuesData[i].m_atoms.at("N").m_position, residuesData[i].m_atoms.at("CA").m_position);
				//C
				residuesData[i].m_geo->CA_C_length = getBondLength(residuesData[i].m_atoms.at("CA").m_position, 
					residuesData[i].m_atoms.at("C").m_position);
				residuesData[i].m_geo->N_CA_C_angle = getBondAngle(residuesData[i].m_atoms.at("N").m_position, 
					residuesData[i].m_atoms.at("CA").m_position, residuesData[i].m_atoms.at("C").m_position);
			}	
		}

		//get model's phipsiomega
		getTorsionAngles(residuesData);
	}
}
//======getResidueInfo======
