#include "side_chain_modeling.h"


//=====================rotamer=====================
//transfer phipsi into RotamerData format
string rotamerFormat(double anlge){
	int angle_i;
	if (anlge > 0){
		angle_i = ((int)(anlge / 10 + 0.5)) * 10;
	}else{
		angle_i = ((int)(anlge / 10 - 0.5)) * 10;
	}
	if (angle_i > 180){
		angle_i -= 360;
	}else if (angle_i < -180){
		angle_i += 360;
	}
	return to_string(angle_i);
}

//add phipsi_rotamerformat to geo
void getPossibleRotamers(Residue& residue){
	residue.m_geo->phipsi_rotamerformat = getTriResname(residue.m_resname) +
		"_" + rotamerFormat(residue.m_geo->phi) + "_" + rotamerFormat(residue.m_geo->psi);	
}

void addAllPossibleRotamers(vector<Residue>& residuesData){
	for (Residue& residue : residuesData){
		if (residue.m_resname == 'G' || residue.m_resname == 'A'){
			continue;
		}
		getPossibleRotamers(residue);
	}
}
//=====================rotamer=====================

//=====================stander DASF=====================
//add standerDASFs to geo
void getStanderDASFs(vector<Residue>& residuesData, const map<string, vector<DASFFeature>>& DASFData, int window_len){

	int ref_id = (window_len - 1) / 2;

	int start_id = 0;
	if (window_len == 7 || window_len == 11){
		start_id = 1;
	}

	int length = residuesData.size();
	for (int i = 0; i <= (length - window_len); ++i){

		//key
		string key;
		//double id_sum = 0;
		for (auto it = residuesData.begin() + i; it < (residuesData.begin() + i + window_len); ++it){
			key += it->m_resname;
			//	id_sum += it->m_resid;
		}

		//assure no skip
		//if (id_sum / window_len != residuesData[i+ref_id].m_resid){
		//	continue;
		//}

		auto it = DASFData.find(key);
		if (it != DASFData.end()){
			//5 & 7 & 9 & 11
			int index = 0;
			for (int j = start_id; j < window_len; j = j + 2){
				//residuesData[i + j].m_geo->standerDASFs.clear();
				residuesData[i + j].m_geo->standerDASFs.push_back(
					vector<DASFFeature>(it->second.begin() + index, 
					it->second.begin() + index + residuesData[i + j].m_geo->m_num_dasfs));
				index += residuesData[i + j].m_geo->m_num_dasfs;
			}
			assert(index == it->second.size());
		}
	}
}

void addAllStanderDASFs(vector<Residue>& residuesData, const vector<map<string, vector<DASFFeature>>>& DASFData){
	getStanderDASFs(residuesData, DASFData[0], 5);
	getStanderDASFs(residuesData, DASFData[1], 7);
	getStanderDASFs(residuesData, DASFData[2], 9);
	getStanderDASFs(residuesData, DASFData[3], 11);
}
//=====================stander DASF=====================


void initSideChain(vector<Residue>& residuesData, const map<string, vector<RotamerFeature>>& RotamerData,
	const vector<map<string, vector<DASFFeature>>>& DASFData){

	//add phipsi_rotamerformat to geo
	addAllPossibleRotamers(residuesData);

	//add standerDASFs to geo
	//TODO:try only keep the longest one
	addAllStanderDASFs(residuesData, DASFData);

	int length = residuesData.size();
	for (int i = 0; i < length; ++i){

		if (residuesData[i].m_resname == 'G' || residuesData[i].m_resname == 'A'){
			continue;
		}
		vector<RotamerFeature> rotamers = RotamerData.at(residuesData[i].m_geo->phipsi_rotamerformat);
		DASFReference reference = getDASFReference(residuesData[i]);
		double score_min = 999999;
		int id_min = 0;
		double dasf_score = 0;
		double rotamer_score = 0;
		double total_scores = 0;
		int length = stoi(CONFIG["top_k"]) < rotamers.size() ? stoi(CONFIG["top_k"]) : rotamers.size();
		for (int k = 0; k < length; ++k){
			residuesData[i].m_geo->setRotamers(rotamers[k].m_dihedral);
			residuesData[i].m_geo->addSideChain();
			dasf_score = getDASFPotential(residuesData[i], reference);
			//rotamer_score = rotamers[k].m_prob;
			rotamer_score = -log(rotamers[k].m_prob / rotamers[0].m_prob);
			rotamer_score = rotamer_score > 5 ? 5 : rotamer_score;
			rotamer_score = rotamer_score < -5 ? -5 : rotamer_score;
			total_scores = stof(CONFIG["w_dasf"])*dasf_score + stof(CONFIG["w_rotamer"])*rotamer_score;
			if (total_scores < score_min){
				score_min = total_scores;
				id_min = k;
			}
			//cout << stof(CONFIG["w_dasf"])*dasf_score << ' ' << stof(CONFIG["w_rotamer"])*rotamer_score << endl;
		}
		residuesData[i].m_geo->dihedral = rotamers[id_min].m_dihedral;
		residuesData[i].m_geo->dihedral_score = score_min;
	}

}

void selectOptRotamer(vector<Residue>& residuesData, const vector<int>& mutation_indexs, 
	const map<string, vector<RotamerFeature>>& RotamerData){
	
	int length = mutation_indexs.size();
	for (int i = 0; i < length; ++i){
		if (residuesData[mutation_indexs[i]].m_resname == 'G' || residuesData[mutation_indexs[i]].m_resname == 'A'){
			continue;
		}

		//add phipsi_rotamerformat to geo->phipsi_rotamerformat
		getPossibleRotamers(residuesData[mutation_indexs[i]]);

		vector<RotamerFeature> rotamers = RotamerData.at(residuesData[mutation_indexs[i]].m_geo->phipsi_rotamerformat);
		DASFReference reference = getDASFReference(residuesData[mutation_indexs[i]]);
		double score_min = 999999;
		int id_min = 0;
		double dasf_score = 0;
		double rotamer_score = 0;
		double total_scores = 0;
		int length = stoi(CONFIG["top_k"]) < rotamers.size() ? stoi(CONFIG["top_k"]) : rotamers.size();
		for (int k = 0; k < length; ++k){
			residuesData[mutation_indexs[i]].m_geo->setRotamers(rotamers[k].m_dihedral);
			residuesData[mutation_indexs[i]].m_geo->addSideChain();
			dasf_score = getDASFPotential(residuesData[mutation_indexs[i]], reference);
			//rotamer_score = rotamers[k].m_prob;
			rotamer_score = -log(rotamers[k].m_prob / rotamers[0].m_prob);
			rotamer_score = rotamer_score > 5 ? 5 : rotamer_score;
			rotamer_score = rotamer_score < -5 ? -5 : rotamer_score;
			total_scores = stof(CONFIG["w_dasf"])*dasf_score + stof(CONFIG["w_rotamer"])*rotamer_score;
			if (total_scores < score_min){
				score_min = total_scores;
				id_min = k;
			}
			//cout << stof(CONFIG["w_dasf"])*dasf_score << ' ' << stof(CONFIG["w_rotamer"])*rotamer_score << endl;
		}
		residuesData[mutation_indexs[i]].m_geo->dihedral = rotamers[id_min].m_dihedral;
		residuesData[mutation_indexs[i]].m_geo->dihedral_score = score_min;
	}
}