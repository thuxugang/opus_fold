#include "myutils.h"

vector<string> splits(stringstream& ss, const string &s, char delim){
	ss << s;
	string item;
	vector<string> elems;
	while (getline(ss, item, delim)){
		elems.push_back(item);
	}
	ss.clear();
	return elems;
}

vector<double> splitd(stringstream& ss, const string &s, char delim){
	ss << s;
	string item;
	vector<double> elems;
	while (getline(ss, item, delim)){
		elems.push_back(atof(item.c_str()));
	}
	ss.clear();
	return elems;
}

RandomThreads::RandomThreads(int seed){
	this->E.seed(seed);
	this->u_real.param(uniform_real_distribution<>::param_type{ 0, 1 });
}

vector<int> getRandomIndexs(int num_total, int num_indexs, RandomThreads& random_thread){

	vector<int> v;
	for (int i = 0; i < num_total; ++i) {
		v.push_back(i);
	}
	shuffle(v.begin(), v.end(), random_thread.E);

	return vector<int>(v.begin(), v.begin() + num_indexs);
}

//==============getPossiblePhipsis==============
void addTAData(vector<vector<TAFeature>>& possiblePhipsis, const vector<Residue>& residuesData, 
	const map<string, vector<TAFeature>>& TAData, int window_len){

	int ref_id = (window_len - 1) / 2;
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
		//if (id_sum / window_len != residuesData[i + ref_id].m_resid){
		//	continue;
		//}

		vector<TAFeature> data;
		try{
			data = TAData.at(key);
		}
		catch (const out_of_range& oor){
			continue;
		}
		catch (const exception &e){
			cout << "myutils.addTAData() wrong: " << e.what() << endl;
		}
		possiblePhipsis[i + ref_id] = data;
	}
}

bool TA_comp(const TAFeature &a, const TAFeature &b){
	return a.m_possibility > b.m_possibility;
}

vector<vector<TAFeature>> getPossiblePhiPsis(const vector<Residue>& residuesData, const vector<map<string, vector<TAFeature>>>& TAData){
	vector<vector<TAFeature>> possiblePhipsis(residuesData.size());
	addTAData(possiblePhipsis, residuesData, TAData[0], 3);
	addTAData(possiblePhipsis, residuesData, TAData[1], 5);
	addTAData(possiblePhipsis, residuesData, TAData[2], 7);
	
	//first and last
	string key;
	key += "G";
	key += residuesData[0].m_resname;
	key += residuesData[1].m_resname;
	possiblePhipsis[0] = TAData[0].at(key);
	key.clear();
	key += residuesData[residuesData.size() - 2].m_resname;
	key += residuesData[residuesData.size() - 1].m_resname;
	key += "G";
	possiblePhipsis[residuesData.size() - 1] = TAData[0].at(key);

	//sort
	for (auto & possiblePhipsi : possiblePhipsis){
		sort(possiblePhipsi.begin(), possiblePhipsi.end(), TA_comp);
	}

	return possiblePhipsis;
}
//==============getPossiblePhipsis==============


//==============setPhiPsiFromTA==============
Vector2d sampleFromGMM(const TAFeature& feature, RandomThreads& random_thread){

	VectorXd x(feature.m_gmm_sampler.m_dimension);
	for (int i = 0; i < feature.m_gmm_sampler.m_dimension; ++i){
		x[i] = random_thread.normal(random_thread.E);
	}
	VectorXd phi_psi = feature.m_gmm_sampler.m_Q * x + feature.m_gmm_sampler.m_mean;
	//cout << phi_psi << feature.m_gmm_sampler.m_mean << feature.m_mean << endl;
	return phi_psi;	

}

int selectGMMFromTA(const vector<TAFeature>& features, RandomThreads& random_thread){

	if (random_thread.u_real(random_thread.E) < 0.6){
		random_thread.u_real_ta.param(uniform_real_distribution<>::param_type{ 0, features[0].m_total_possibility });
		double value = random_thread.u_real_ta(random_thread.E);

		int i = 0;
		double sum = 0;
		for (; i < features.size(); ++i){
			sum += features[i].m_possibility;
			if (sum > value){
				break;
			}
		}
		return i;	
	}else{
		random_thread.u_int_ta.param(uniform_int_distribution<>::param_type{ 0, int(features.size() - 1) });
		return random_thread.u_int_ta(random_thread.E);	
	}
}

void initPhiPsiFromTA(vector<Residue>& residuesData, const vector<vector<TAFeature>>& possiblePhipsis, int init_index){
	assert(residuesData.size() == possiblePhipsis.size());
	int length = residuesData.size();
	for (int i = 0; i < length; ++i){
		residuesData[i].m_geo->ta_possibility = possiblePhipsis[i][init_index].m_possibility;
		Vector2d phipsi = possiblePhipsis[i][init_index].m_mean;
		residuesData[i].m_geo->phi = phipsi[0];
		residuesData[i].m_geo->psi = phipsi[1];
	}
}

void sampleSpecificPhiPsiFromTA(vector<Residue>& residuesData, const vector<vector<TAFeature>>& possiblePhipsis,
	const vector<int>& mutation_indexs, RandomThreads& random_thread){
	int length = mutation_indexs.size();
	for (int i = 0; i < length; ++i){
		int index = selectGMMFromTA(possiblePhipsis[mutation_indexs[i]], random_thread);
		residuesData[mutation_indexs[i]].m_geo->ta_possibility = possiblePhipsis[mutation_indexs[i]][index].m_possibility;
		Vector2d phipsi = sampleFromGMM(possiblePhipsis[mutation_indexs[i]][index], random_thread);
		residuesData[mutation_indexs[i]].m_geo->phi = phipsi[0];
		residuesData[mutation_indexs[i]].m_geo->psi = phipsi[1];
	}
}

void sampleOmega(vector<Residue>& residuesData, const vector<int>& mutation_indexs, RandomThreads& random_thread){
	int length = mutation_indexs.size();
	for (int i = 0; i < length; ++i){
		if (residuesData[mutation_indexs[i]].m_resname == 'P' && random_thread.u_real(random_thread.E) > 0.8){
			residuesData[mutation_indexs[i]].m_geo->omega = 0;
		}else if (residuesData[mutation_indexs[i]].m_resname != 'P' && random_thread.u_real(random_thread.E) > 0.997){
			residuesData[mutation_indexs[i]].m_geo->omega = 0;
		}
	}
}
//==============setPhiPsiFromTA==============

//==============save & restore oldone==============
OldPhiPsi::OldPhiPsi(double phi, double psi, double omega, double ta_possibility){
	this->m_phi = phi;
	this->m_psi = psi;
	this->m_omega = omega;
	this->m_ta_possibility = ta_possibility;
	this->m_rotamer_score = 0;
}

//save old phipsi and rotamer
vector<OldPhiPsi> saveOldPhiPsi(const vector<Residue>& residuesData, const vector<int>& mutation_indexs, bool use_sidechain){
	vector<OldPhiPsi> old_phipsi;
	int length = mutation_indexs.size();
	for (int i = 0; i < length; ++i){
		old_phipsi.emplace_back(residuesData[mutation_indexs[i]].m_geo->phi, residuesData[mutation_indexs[i]].m_geo->psi,
			residuesData[mutation_indexs[i]].m_geo->omega, residuesData[mutation_indexs[i]].m_geo->ta_possibility);
		if (use_sidechain){
			if (residuesData[mutation_indexs[i]].m_geo->m_resname == 'G' || residuesData[mutation_indexs[i]].m_geo->m_resname == 'A'){
				continue;
			}
			old_phipsi[i].m_rotamer = residuesData[mutation_indexs[i]].m_geo->dihedral;
			old_phipsi[i].m_rotamer_score = residuesData[mutation_indexs[i]].m_geo->dihedral_score;
		}
	}
	return old_phipsi;
}

void restoreOldPhiPsi(const vector<OldPhiPsi>& old_phipsi, vector<Residue>& residuesData, const vector<int>& mutation_indexs, 
	bool use_sidechain){
	int length = mutation_indexs.size();
	for (int i = 0; i < length; ++i){
		residuesData[mutation_indexs[i]].m_geo->phi = old_phipsi[i].m_phi;
		residuesData[mutation_indexs[i]].m_geo->psi = old_phipsi[i].m_psi;
		residuesData[mutation_indexs[i]].m_geo->omega = old_phipsi[i].m_omega;
		residuesData[mutation_indexs[i]].m_geo->ta_possibility = old_phipsi[i].m_ta_possibility;
		if (use_sidechain){
			if (residuesData[mutation_indexs[i]].m_geo->m_resname == 'G' || residuesData[mutation_indexs[i]].m_geo->m_resname == 'A'){
				continue;
			}
			residuesData[mutation_indexs[i]].m_geo->dihedral = old_phipsi[i].m_rotamer;
			residuesData[mutation_indexs[i]].m_geo->dihedral_score = old_phipsi[i].m_rotamer_score;
		}
	}
}
//==============save & restore oldone==============
