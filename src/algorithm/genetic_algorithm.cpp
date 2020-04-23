#include "genetic_algorithm.h"

GA::GA(const vector<Residue>& residuesData, double coverage, RandomThreads& random_thread, int index){

	this->index = index;

	random_thread.u_int_res.param(uniform_int_distribution<>::param_type{ 0, int(residuesData.size() - 1) });
	
	this->m_coverage = coverage;
	this->global_refine = false;
	this->use_sidechain = false;
	this->use_potentials_in_selection = false;
	this->use_dssp = false;

	this->res_length = residuesData.size();
	this->ga_steps = stof(CONFIG["ga_steps"]);
	this->ga_steps_counter = 0;

	this->T = stof(CONFIG["Temperature"]);
	this->delta = stof(CONFIG["delta"]);

	this->mutation_ratio = stof(CONFIG["mutation_ratio"]); 
	this->combination_ratio = stof(CONFIG["combination_ratio"]);
	this->mutation_times_outer = stof(CONFIG["mutation_times_outer"]);
	this->mutation_times_inner = stof(CONFIG["mutation_times_inner"]);

}

void GA::init(vector<Residue>& residuesData, const vector<vector<TAFeature>>& possiblePhipsis,
	const CMInfo& cm_cons, const vector<map<string, vector<CSFFeature>>>& CSFData, 
	const vector<DSSPInfo>& dssp_cons, const vector<Vector2d>& theta_tau_cons){

	//use first means
	if (stoi(CONFIG["init_cons_dynamic"])){
		initPhiPsiFromTA(residuesData, possiblePhipsis, 1);
	}else{
		initPhiPsiFromTA(residuesData, possiblePhipsis, 0);
	}
	
	//reconstruct from residue.geo
	reconstructFromGeo(residuesData);
	//copy information from m_geo to m_atoms
	geoToResidue(residuesData);

	double before_opt = getTotalPotentials(residuesData, possiblePhipsis, cm_cons, CSFData, *this, dssp_cons, theta_tau_cons);
	//cout << "Before opt: " << before_opt << endl;
	this->current_score = before_opt;

}

double GA::mutation_op(vector<Residue>& residuesData, const vector<vector<TAFeature>>& possiblePhipsis, 
	const CMInfo& cm_cons, const vector<map<string, vector<CSFFeature>>>& CSFData,
	const vector<int>& mutation_indexs, RandomThreads& random_thread, const map<string, vector<RotamerFeature>>& RotamerData,
	const vector<DSSPInfo>& dssp_cons, const vector<Vector2d>& theta_tau_cons){

	double total_score_old = this->current_score;
	double total_score = 0;
	for (int i = 0; i < this->mutation_times_inner; ++i){

		//save old phipsi and rotamer
		vector<OldPhiPsi> old_phipsi = saveOldPhiPsi(residuesData, mutation_indexs, this->use_sidechain);

		//sampling
		sampleSpecificPhiPsiFromTA(residuesData, possiblePhipsis, mutation_indexs, random_thread);
		sampleOmega(residuesData, mutation_indexs, random_thread);

		//reconstruct from residue.geo
		reconstructFromGeo(residuesData);
		if (this->use_sidechain){
			selectOptRotamer(residuesData, mutation_indexs, RotamerData);
			reconstructSideChainFromGeo(residuesData);
		}

		//copy information from m_geo to m_atoms
		geoToResidue(residuesData);

		total_score = getTotalPotentials(residuesData, possiblePhipsis, cm_cons, CSFData, *this, dssp_cons, theta_tau_cons);
		
		double delta_E = total_score - total_score_old;
		if (delta_E < 0 || (random_thread.u_real(random_thread.E)<exp(-delta_E / this->T))){
			//substitute
			total_score_old = total_score;		
		}else{
			//restore old phipsi
			restoreOldPhiPsi(old_phipsi, residuesData, mutation_indexs, this->use_sidechain);
			total_score = total_score_old;
		}
	}
	return total_score;
}


int choose_op(vector<double>& scores){
	vector<double>::iterator minpos = min_element(scores.begin(), scores.end());
	int min_id = distance(begin(scores), minpos);
	return min_id;
}

void combination_op(const vector<Residue>& residuesData_min, vector<Residue>& residuesData, const vector<int>& combination_indexs){

	int length = combination_indexs.size();
	for (int i = 0; i < length; ++i){
		residuesData[combination_indexs[i]].m_geo->phi = residuesData_min[combination_indexs[i]].m_geo->phi;
		residuesData[combination_indexs[i]].m_geo->psi = residuesData_min[combination_indexs[i]].m_geo->psi;
		residuesData[combination_indexs[i]].m_geo->omega = residuesData_min[combination_indexs[i]].m_geo->omega;
		residuesData[combination_indexs[i]].m_geo->ta_possibility = residuesData_min[combination_indexs[i]].m_geo->ta_possibility;
		residuesData[combination_indexs[i]].m_geo->dihedral = residuesData_min[combination_indexs[i]].m_geo->dihedral;
		residuesData[combination_indexs[i]].m_geo->dihedral_score = residuesData_min[combination_indexs[i]].m_geo->dihedral_score;
	}
}


double getTotalPotentials(const vector<Residue>& residuesData, const vector<vector<TAFeature>>& possiblePhipsis,
	const CMInfo& cm_cons, const vector<map<string, vector<CSFFeature>>>& CSFData, const GA& ga,
	const vector<DSSPInfo>& dssp_cons, const vector<Vector2d>& theta_tau_cons){

	int length = residuesData.size();

	/*
	//output pdb
	string dir_path;
	if (stoi(CONFIG["dssp_cons"])){
		dir_path = CONFIG["tmp_dir"] + "/tmp_" + to_string(ga.index) + ".pdb";
		outputPDB(dir_path, residuesData, ga.use_sidechain);
	}
	*/

	//cal asa and dssp
	double asa_score = 0;
	double dssp8_score = 0;
	if (stoi(CONFIG["dssp_cons"]) && ga.use_dssp){

		//output atomsData
		vector<Atom> atomsData;
		atomsData = getAtomsData(residuesData, ga.use_sidechain);
		
		vector<DSSPInfo> dssp_infos = getDSSPResults(atomsData, ga.use_sidechain);
		assert(dssp_infos.size() == length);
		for (int i = 0; i < length; ++i){
			dssp8_score += dssp_cons[i].m_dssp_dict.at(dssp_infos[i].m_dssp);
			if (ga.use_sidechain){
				asa_score += abs(dssp_infos[i].m_asa - dssp_cons[i].m_asa);
			}
		}
		
	}

	//refine local first
	double cm_cons_score = 0;
	if (stoi(CONFIG["cm_cons"]) && ga.global_refine){
		//cal cm_cons_score
		CMInfo cm_preidict = getContactMapInfo(residuesData);
		cm_cons_score = getCMConsScore(cm_cons, cm_preidict);
	}

	//cal CSFScores
	double csf_score = getCSFPotential(residuesData, CSFData);

	//cal ta_cons_score
	double ta_cons_score = 0;
	if (stoi(CONFIG["ta_cons"])){
		ta_cons_score = getTAPontial(residuesData, possiblePhipsis);
	}

	//cal theta_tau_score
	double tt_cons_score = 0;
	if (stoi(CONFIG["tt_cons"])){
		tt_cons_score = getTTPontial(residuesData, theta_tau_cons);
	}

	//possibility
	double ta_scores = 0;
	for (int i = 0; i < length; ++i){
		/*
		double ta_score = log(residuesData[i].m_geo->ta_possibility / possiblePhipsis[i][0].m_possibility);
		if (ta_score > 5){
			ta_score = 5;
		}else if (ta_score < -5){
			ta_score = -5;
		}
		
		ta_scores += ta_score;
		*/
		ta_scores += residuesData[i].m_geo->ta_possibility;
	}

	//cal LJ
	double lj_scores = 0;
	if (ga.use_sidechain){

		lj_scores = stof(CONFIG["w_mmlj"])*getMMLJPotential(residuesData) + 
			stof(CONFIG["w_mslj"])*ga.m_coverage*getMSLJPotential(residuesData) +
			stof(CONFIG["w_sslj"])*ga.m_coverage*getSSLJPotential(residuesData);
	}

	//cal dasf
	double dasf_scores = 0;
	if (ga.use_sidechain){
		for (int i = 0; i < length; ++i){
			if (residuesData[i].m_resname == 'G' || residuesData[i].m_resname == 'A'){
				continue;
			}
			dasf_scores += residuesData[i].m_geo->dihedral_score;
		}
	}

	double total_score = stof(CONFIG["w_csf"])*ga.m_coverage*csf_score + stof(CONFIG["w_ta"])*(-ta_scores) +
		stof(CONFIG["w_ta_cons"])*ta_cons_score + stof(CONFIG["w_tt_cons"])*tt_cons_score +
		stof(CONFIG["w_cm_cons"])*cm_cons_score + lj_scores + stof(CONFIG["w_side_chain"])*dasf_scores + 
		stof(CONFIG["w_dssp8"])*(-dssp8_score) + stof(CONFIG["w_asa"])*asa_score;

	//cout << stof(CONFIG["w_csf"])*ga.m_coverage*csf_score << ' ' << stof(CONFIG["w_ta"])*(-ta_scores) << ' ' <<
	//	stof(CONFIG["w_ta_cons"])*ta_cons_score << ' ' << stof(CONFIG["w_cm_cons"])*cm_cons_score << ' ' <<
	//	lj_scores << ' ' << stof(CONFIG["w_side_chain"])*dasf_scores << endl;

	return total_score;
}
