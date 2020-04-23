#include "contact_map_constraint.h"

map<string, string> CM_list;

//exclude the contact if the distance of their res_id is less than GAP
int GAP;

CMInfo readCMConstraintFiles(const string& name, int size){

	GAP = stoi(CONFIG["cm_gap"]);

	CMInfo cm_info;
	cm_info.m_cm_bool = MatrixXd(size, size);
	cm_info.m_cm_bool.setZero();
	cm_info.m_cm_value = MatrixXd(size, size);
	cm_info.m_cm_value.setZero();

	stringstream ss;

	ifstream cmFile;
	try{
		cmFile.open(CM_list.at(name));
	}catch(...){
		cout << "Contact map list format error!";
		cmFile.close();
		exit(-1);
	}
	
	string str_line;
	try{
		while (!cmFile.eof()){
			getline(cmFile, str_line);
			str_line = trim(str_line);
			if (str_line.find('#') == 0 || str_line.find('=') == 0 || str_line.size() == 0){
				continue;
			}
			char spliter = ' ';
			if (str_line.find('\t') != str_line.npos){
				spliter = '\t';
			}
			vector<double> items = splitd(ss, str_line, spliter);
			//on
			if ((items[1] - items[0]) >= GAP){
				cm_info.m_cm_value(int(items[0] - 1), int(items[1] - 1)) = items[4];
				if (items[4] > 0.5){
					cm_info.m_cm_bool(int(items[0] - 1), int(items[1] - 1)) = 1;
				}
			}		
		}
	}catch (...){
		cout << "Contact map file format error!";
		cmFile.close();
		exit(-1);
	}
	cmFile.close();

	return cm_info;
}

CMInfo getContactMapInfo(const vector<Residue>& residuesData){
	
	int length = residuesData.size();

	CMInfo cm_info;
	cm_info.m_cm_bool = MatrixXd(length, length);
	cm_info.m_cm_bool.setZero();
	cm_info.m_cm_value = MatrixXd(length, length);
	cm_info.m_cm_value.setZero();

	const Vector3d* i_position;
	const Vector3d* j_position;
	double distance = 0;

	for (int i = 0; i < length; ++i){
		if (residuesData[i].m_resname == 'G'){
			i_position = &(residuesData[i].m_atoms.at("CA").m_position);
		}else{
			i_position = &(residuesData[i].m_atoms.at("CB").m_position);
		}
		for (int j = i + GAP; j < length; ++j){
			if (residuesData[j].m_resname == 'G'){
				j_position = &(residuesData[j].m_atoms.at("CA").m_position);
			}
			else{
				j_position = &(residuesData[j].m_atoms.at("CB").m_position);
			}
			distance = getBondLength(*i_position, *j_position);
			if (distance <= 8){
				cm_info.m_cm_bool(i, j) = 1;
			}
			cm_info.m_cm_value(i, j) = distance;
		}
	
	}	

	return cm_info;
}

double getCMConsScore(const CMInfo& cm_cons, const CMInfo& cm_predict){
	double scores = 0;
	MatrixXd diff_bool = (cm_cons.m_cm_bool - cm_predict.m_cm_bool).cwiseAbs();
	for (int i = 0; i < diff_bool.cols(); ++i){
		for (int j = i + GAP; j < diff_bool.rows(); ++j){
			if (diff_bool(i, j) == 1){
				scores += (abs(cm_predict.m_cm_value(i, j) - 8)*(0.5+abs(cm_cons.m_cm_value(i, j) - 0.5)));
			}
		}
	}
	return scores;
}
