#include "theta_tau_constraint.h"

map<string, string> TT_list;

vector<Vector2d> readTTConstraintFiles(const string& name){

	vector<Vector2d> thetatau_cons;

	stringstream ss;

	ifstream ttFile;
	try{
		ttFile.open(TT_list.at(name));
	}catch(...){
		cout << "Theta tau list format error!";
		ttFile.close();
		exit(-1);
	}
	
	string str_line;
	try{
		while (!ttFile.eof()){
			getline(ttFile, str_line);
			str_line = trim(str_line);
			if (str_line.find('#') == 0 || str_line.find('=') == 0 || str_line.size() == 0){
				continue;
			}
			char spliter = ' ';
			if (str_line.find('\t') != str_line.npos){
				spliter = '\t';
			}
			vector<double> items = splitd(ss, str_line, spliter);
			thetatau_cons.push_back(Vector2d(items[stoi(CONFIG["row_theta"])], items[stoi(CONFIG["row_tau"])]));
		}
	}catch (...){
		cout << "Theta tau file format error!";
		ttFile.close();
		exit(-1);
	}
	ttFile.close();

	return thetatau_cons;
}


double getTTPontial(const vector<Residue>& residuesData, const vector<Vector2d>& thetatau_cons){
	int length = residuesData.size();
	double tt_scores = 0;
	//exclude first one and last two
	for (int i = 1; i < (length - 2); ++i){
		double theta = getBondAngle(residuesData[i - 1].m_atoms.at("CA").m_position, residuesData[i].m_atoms.at("CA").m_position,
			residuesData[i + 1].m_atoms.at("CA").m_position);
		double tau = calDihedral(residuesData[i - 1].m_atoms.at("CA").m_position, residuesData[i].m_atoms.at("CA").m_position,
			residuesData[i + 1].m_atoms.at("CA").m_position, residuesData[i + 2].m_atoms.at("CA").m_position);
		tt_scores += (abs(theta - thetatau_cons[i][0]) + abs(tau - thetatau_cons[i][1]));
	}
	return tt_scores;
}
