#include "multi_threads.h"


vector<Residue> optimize_mt(const vector<vector<TAFeature>>& possiblePhipsis, const CMInfo& cm_cons,
	const vector<map<string, vector<CSFFeature>>>& CSFData, double coverage, const vector<Atom>& atomsData_cons,
	const string& fasta, const map<string, vector<RotamerFeature>>& RotamerData, 
	const vector<map<string, vector<DASFFeature>>>& DASFData, const vector<DSSPInfo>& dssp_cons,
	const vector<Vector2d>& theta_tau_cons){

	int threads_num = stoi(CONFIG["threads_num"]);
	int ga_population = stoi(CONFIG["ga_population"]);
omp_set_num_threads(threads_num);
	
	cout << "Start optimization..." << endl;

	//copy initial residuesData
	//TODO: copy is hard to implement
	vector<vector<Residue>> residuesDatas_threads;
	vector<RandomThreads> random_threads;
#pragma omp parallel for ordered schedule(dynamic)
	for (int i = 0; i < ga_population; ++i){
		//if use initial structure constraint
		if (stoi(CONFIG["init_cons_dynamic"])){
			bool MAIN_CHAIN = true;
			vector<Residue> residuesData_threads = atomToResidue(atomsData_cons, MAIN_CHAIN);
			bool FROM_MODEL = true;
			getResidueInfo(residuesData_threads, FROM_MODEL);
#pragma omp ordered
			residuesDatas_threads.push_back(residuesData_threads);
#pragma omp ordered
			random_threads.push_back(RandomThreads(i));
		}else{
			//construct residuesData from fasta file
			vector<Residue> residuesData_threads = fastaToResiduesData(fasta);
			bool FROM_MODEL = false;
			getResidueInfo(residuesData_threads, FROM_MODEL);
#pragma omp ordered
			residuesDatas_threads.push_back(residuesData_threads);
#pragma omp ordered
			random_threads.push_back(RandomThreads(i));
		}					
	}
	
	//ga init
	vector<GA> gas_threads;	
#pragma omp parallel for ordered schedule(dynamic)
	for (int i = 0; i < ga_population; ++i){
		//genetic algorithm mutation
		GA ga = GA(residuesDatas_threads[i], coverage, random_threads[i], i);
		ga.init(residuesDatas_threads[i], possiblePhipsis, cm_cons, CSFData, dssp_cons, theta_tau_cons);
#pragma omp ordered
		gas_threads.push_back(ga);
	}

	//directly from global refine
	if (stoi(CONFIG["init_cons_dynamic"])){
		for (int i = 0; i < ga_population; ++i){
			gas_threads[i].global_refine = true;		
		}
		//cout << "Start global refine... " << endl;
	}

	//output temp structure list
	ofstream fout(CONFIG["tmp_dir"] + "/tmp.txt");
	char data[80];
	for (int i = 0; i < ga_population; ++i){
		fout << CONFIG["tmp_dir"] << "/tmp_" << i << endl;
	}
	fout << flush;
	fout.close();

	//ga optimization
	int min_id = 0;
	while (gas_threads[0].ga_steps_counter < gas_threads[0].ga_steps){

		//only dasf term will be taken into consideration since we want side-chain conformation to be a driven force, 
		//LJ potential is the medium.
		//side-chain info init
		if (stoi(CONFIG["sc_cons"]) && stoi(CONFIG["sc_start_step"]) == gas_threads[0].ga_steps_counter){			
#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < ga_population; ++i){

				gas_threads[i].use_sidechain = true;
				if (stoi(CONFIG["dssp_cons"])){
					gas_threads[i].use_dssp = true;
				}				

				reconstructFromGeo(residuesDatas_threads[i]);
				
				initSideChain(residuesDatas_threads[i], RotamerData, DASFData);
				
				reconstructSideChainFromGeo(residuesDatas_threads[i]);
				geoToResidue(residuesDatas_threads[i]);
				gas_threads[i].current_score = getTotalPotentials(residuesDatas_threads[i],
					possiblePhipsis, cm_cons, CSFData, gas_threads[i], dssp_cons, theta_tau_cons);
				if (i == 0){
					cout << "Add side-chain info... " << endl;
					if (stoi(CONFIG["dssp_cons"])){
						cout << "Add dssp info... " << endl;
					}
				}

			}
#pragma omp barrier
		}

		//output temp structure for other potentials in the selection step
		if (stoi(CONFIG["use_potentials_in_selection"]) && 
			stoi(CONFIG["other_potentials_start_step"]) == gas_threads[0].ga_steps_counter){
			for (int i = 0; i < ga_population; ++i){
				gas_threads[i].use_potentials_in_selection = true;
				if (i == 0){
					cout << "Use other potentials in the selection step... " << endl;
				}
			}

		}

		//================mutation_op================
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < ga_population; ++i){
			//if (gas_threads[i].ga_steps_counter % 10 == 0){
			//	gas_threads[i].use_other_potential = true;
			//}else{
			//	gas_threads[i].use_other_potential = false;
			//}
			//local refine
			//wrong but works! (!gas_threads[i].global_refine && )
			if (gas_threads[i].ga_steps_counter == int(0.1*gas_threads[0].ga_steps)){

				gas_threads[i].global_refine = true;

				reconstructFromGeo(residuesDatas_threads[i]);
				geoToResidue(residuesDatas_threads[i]);
				
				gas_threads[i].current_score = getTotalPotentials(residuesDatas_threads[i],
					possiblePhipsis, cm_cons, CSFData, gas_threads[i], dssp_cons, theta_tau_cons);
				if (i == 0){
					cout << "Start global refine... " << endl;
				}
			}

			//mutation outer
			int mutation_num = 1;
			if (gas_threads[i].ga_steps_counter < 0.9*gas_threads[i].ga_steps){
				mutation_num = int(gas_threads[i].res_length*gas_threads[i].mutation_ratio
					*(1 - gas_threads[i].ga_steps_counter / (0.9*gas_threads[i].ga_steps)) + 1);
			}
			vector<int> mutation_indexs;
			for (int j = 0; j < gas_threads[i].mutation_times_outer; ++j){

				mutation_indexs = getRandomIndexs(residuesDatas_threads[i].size(), mutation_num, random_threads[i]);

				//mutation inner
				gas_threads[i].current_score = gas_threads[i].mutation_op(residuesDatas_threads[i], possiblePhipsis,
					cm_cons, CSFData, mutation_indexs, random_threads[i], RotamerData, dssp_cons, theta_tau_cons);

				//cout << i << " : Mutation num: " << mutation_num << " Score: " << gas_threads[i].current_score <<
				//	" T: " << gas_threads[i].T << endl;
			}

			//cout << gas_threads[i].ga_steps_counter << " : Mutation num: " << mutation_indexs.size() << " Score: " << gas_threads[i].current_score <<
			//	" T: " << gas_threads[i].T << endl;
		
			++gas_threads[i].ga_steps_counter;
			gas_threads[i].T *= gas_threads[i].delta;

			if (gas_threads[i].use_potentials_in_selection){
				string dir_path = CONFIG["tmp_dir"] + "/tmp_" + to_string(i) + ".pdb";
				//reconstruct from residue.geo
				reconstructFromGeo(residuesDatas_threads[i]);
				if (gas_threads[i].use_sidechain){
					reconstructSideChainFromGeo(residuesDatas_threads[i]);
				}
				//copy information from m_geo to m_atoms
				geoToResidue(residuesDatas_threads[i]);

				outputPDB(dir_path, residuesDatas_threads[i], gas_threads[i].use_sidechain);
			}
		}
#pragma omp barrier

		//================choose_op================
		//find the structure with minimum potential
		vector<double> scores;
		for (int i = 0; i < ga_population; ++i){
			//cout << gas_threads[i].current_score << endl;
			scores.push_back(gas_threads[i].current_score);
		}

		//if use potentials in the selection step
		if (gas_threads[0].use_potentials_in_selection){			
			vector<double> other_potentials = getOtherPotentials(ga_population);			
			assert(other_potentials.size() == ga_population);			
			for (int i = 0; i < ga_population; ++i){
				scores[i] += other_potentials[i];
			}
		}

		min_id = choose_op(scores);
		//cout << gas_threads[0].ga_steps_counter << ' ' << scores[min_id] << endl;
		
		//not last one
		if (gas_threads[0].ga_steps_counter != (gas_threads[0].ga_steps - 1)){
			//================combination_op================
			//combine the structure with minimum potential with the others
#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < ga_population; ++i){
				if (i != min_id){
					int combination_num = (gas_threads[i].combination_ratio * gas_threads[i].res_length + 1);
					vector<int> combination_indexs = getRandomIndexs(residuesDatas_threads[i].size(), combination_num, random_threads[i]);
					combination_op(residuesDatas_threads[min_id], residuesDatas_threads[i], combination_indexs);
				}
			}
#pragma omp barrier			
		}
		
	}

	//reconstruct from residue.geo
	reconstructFromGeo(residuesDatas_threads[min_id]);
	if (gas_threads[min_id].use_sidechain){
		reconstructSideChainFromGeo(residuesDatas_threads[min_id]);
	}
	//copy information from m_geo to m_atoms
	geoToResidue(residuesDatas_threads[min_id]);

	double after_opt = getTotalPotentials(residuesDatas_threads[min_id], possiblePhipsis, cm_cons, CSFData, gas_threads[min_id],
		dssp_cons, theta_tau_cons);
	cout << "After opt: " << after_opt << endl;

	return residuesDatas_threads[min_id];
}