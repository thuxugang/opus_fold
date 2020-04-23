#include "geometry.h"
#include "pdb.h"

using namespace std;

vector<Atom> getAtomsData(const vector<Residue>& residuesData, bool use_sidechain){

	vector<Atom> atomsData;
	for (const Residue & residue : residuesData){
		residue.m_geo->m_id_counter = 0;
		for (auto iter = residue.m_geo->m_geo_atoms.begin(); iter != residue.m_geo->m_geo_atoms.end(); ++iter) {
			residue.m_geo->m_geo_atoms[iter->first]->m_id = residue.m_geo->m_atomid;
		}
		residue.m_geo->output(atomsData);
		if (use_sidechain && residue.m_resname != 'G' && residue.m_resname != 'A'){
			((GeoContainSideChain*)residue.m_geo)->output_sidechain(atomsData);
		}
	}

	return atomsData;
}

void outputPDB(const string& outfile, const vector<Residue>& residuesData, bool use_sidechain){

	vector<Atom> atomsData;
	for (const Residue & residue : residuesData){
		residue.m_geo->m_id_counter = 0;
		for (auto iter = residue.m_geo->m_geo_atoms.begin(); iter != residue.m_geo->m_geo_atoms.end(); ++iter) {
			residue.m_geo->m_geo_atoms[iter->first]->m_id = residue.m_geo->m_atomid;
		}
		residue.m_geo->output(atomsData);
		if (use_sidechain && residue.m_resname != 'G' && residue.m_resname != 'A'){
			((GeoContainSideChain*)residue.m_geo)->output_sidechain(atomsData);
		}
	}

	ofstream fout(outfile);
	char data[80];
	int length = atomsData.size();
	for(int i=0; i< length; i++){
		string triresname = getTriResname(atomsData[i].m_resname);
		double x = atomsData[i].m_position[0];
		double y = atomsData[i].m_position[1];
		double z = atomsData[i].m_position[2];
		sprintf(data, "%-6s%5d%2c%-3s%1c%3s%2c%4d%4c%8.3f%8.3f%8.3f", "ATOM", atomsData[i].m_id, ' ', atomsData[i].m_name.c_str(), ' ', triresname.c_str(), ' ', atomsData[i].m_resid, ' ', x, y, z);
		//sprintf_s(data, "%-6s%5d%2c%-3s%1c%3s%2c%4d%4c%8.3f%8.3f%8.3f","ATOM",i+1, ' ', atomsData[i].m_name.c_str(), ' ',triresname.c_str(),' ',atomsData[i].m_resid, ' ',x,y,z);
		fout << data << endl;
	}
	fout << flush;
	fout.close();
}

vector<Atom> readPDB(const string& infile){

	vector<Atom> atomsData;

	ifstream fin(infile);
	stringstream current_val;

	string current_line;

	string label;
	int atom_id;
	string atom_name;
	string atom_restype_string;
	char atom_restype;
	string triresname;
	char atom_resname;
	int atom_resid;
	double x, y, z;

	while (!fin.eof()){

		getline(fin, current_line);
		//cout << current_line << endl;

		label = current_line.substr(0, 4);

		if (label.compare("ATOM") == 0){

			current_val.str(current_line.substr(6, 5));
			current_val >> atom_id;
			current_val.clear();

			current_val.str(current_line.substr(12, 4));
			current_val >> atom_name;
			current_val.clear();

			atom_restype_string = current_line.substr(16, 1);
			if (atom_restype_string.compare(" ") == 0){
				atom_restype = 'A';
			}
			else{
				current_val.str(atom_restype_string);
				current_val >> atom_restype;
				current_val.clear();
			}

			current_val.str(current_line.substr(17, 3));
			current_val >> triresname;
			atom_resname = getResname(triresname);
			current_val.clear();

			current_val.str(current_line.substr(22, 4));
			current_val >> atom_resid;
			current_val.clear();

			current_val.str(current_line.substr(30, 8));
			current_val >> x;
			current_val.clear();
			current_val.str(current_line.substr(38, 8));
			current_val >> y;
			current_val.clear();
			current_val.str(current_line.substr(46, 8));
			current_val >> z;
			current_val.clear();
			Vector3d atom_position(x, y, z);

			atomsData.emplace_back(atom_id, atom_name, atom_restype, atom_resname, atom_resid, atom_position);
		}
	}
	fin.close();
	return atomsData;
}

//int main(){
//	//vector<Atom> aa;
//	//Vector3d v(1,1,1);
//	//Atom a(10, "CB", 'G', 1, v);
//	//aa.push_back(a);
//	//outputPDB("test.pdb", aa);
//
//	vector<Atom> ad = readPDB("data/ori1a7s.pdb");
//
//	outputPDB("data/1a7s_rebuild.pdb", ad);
//	return 0;
//
//}
