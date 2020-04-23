#include "torsion_angles_constraint.h"

map<string, string> TA_list;

vector<Vector2d> readTAConstraintFiles(const string& name){

	vector<Vector2d> phipsi_cons;

	stringstream ss;

	ifstream taFile;
	try{
		taFile.open(TA_list.at(name));
	}catch(...){
		cout << "Torsion angles list format error!";
		taFile.close();
		exit(-1);
	}
	
	string str_line;
	try{
		while (!taFile.eof()){
			getline(taFile, str_line);
			str_line = trim(str_line);
			if (str_line.find('#') == 0 || str_line.find('=') == 0 || str_line.size() == 0){
				continue;
			}
			char spliter = ' ';
			if (str_line.find('\t') != str_line.npos){
				spliter = '\t';
			}
			vector<double> items = splitd(ss, str_line, spliter);
			phipsi_cons.push_back(Vector2d(items[stoi(CONFIG["row_phi"])], items[stoi(CONFIG["row_psi"])]));
		}
	}catch (...){
		cout << "Torsion angles file format error!";
		taFile.close();
		exit(-1);
	}
	taFile.close();

	return phipsi_cons;
}

void addTAConsToPossiblePhipsis(vector<vector<TAFeature>>& possiblePhipsis, const vector<Vector2d>& phipsi_cons){

	int length = possiblePhipsis.size();
	for (int i = 0; i < length; ++i){

		vector<double> mean = { phipsi_cons[i][0], phipsi_cons[i][1] };
		vector<double> sd = { 10, 0, 0, 10 };
		double possibility = possiblePhipsis[i][0].m_possibility*1.5;
		int window_len = possiblePhipsis[i][0].m_window_len;

		possiblePhipsis[i].insert(possiblePhipsis[i].begin(), TAFeature(mean, sd, possibility, window_len, (1 + possibility)));

	}
}

double getTAPontial(const vector<Residue>& residuesData, const vector<vector<TAFeature>>& possiblePhipsis){
	int length = residuesData.size();
	double ta_scores = 0;
	for (int i = 0; i < length; ++i){
		ta_scores += (abs(residuesData[i].m_geo->phi - possiblePhipsis[i][0].m_mean[0]) +
			abs(residuesData[i].m_geo->psi - possiblePhipsis[i][0].m_mean[1]));
	}
	return ta_scores;
}
