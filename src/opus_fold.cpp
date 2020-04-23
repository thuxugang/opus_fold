#include <iostream> 
#include <string> 
#include <vector> 
#include <map>
#include <time.h>

#include "config.h"
#include "pdb.h"
#include "mongo.h"
#include "residue.h"
#include "CSFScore.h"
#include "theta_tau_constraint.h"
#include "torsion_angles_constraint.h"
#include "contact_map_constraint.h"
#include "init_structure_constraint.h"
#include "myutils.h"
#include "genetic_algorithm.h"
#include "multi_threads.h"
#include "peptideBuilder.h"

using namespace std;

int main(){

	cout << "Start..." << endl;

	//read config
	readConfigFile(string("./opus_fold.ini"));
	cout << "=========Parameters========="<< endl;
	for (auto &conf : CONFIG){
		if (conf.first.compare("") != 0){
			cout << conf.first << ": " << conf.second << endl;
		}
	}
	cout << "=========Parameters=========" << endl;

	//connect mongo
	cout << "Connect database..." << endl;
	Mongo m;
	cout << "Connection established..." << endl;

	//if use torsion angles constraint
	if (stoi(CONFIG["ta_cons"])){
		//init ta_constraint_files path to torsion_angles_constraint::TA_list
		getConstraintList(CONFIG["ta_cons_list"], TA_list);
	}

	//if use theta tau constraint
	if (stoi(CONFIG["tt_cons"])){
		//init tt_constraint_files path to torsion_angles_constraint::TT_list
		getConstraintList(CONFIG["tt_cons_list"], TT_list);
	}

	//if use dssp8 and asa constraints
	if (stoi(CONFIG["dssp_cons"])){
		//init dssp8_constraint_files path to dssp8_constraint::DSSP_list
		getConstraintList(CONFIG["dssp_cons_list"], DSSP_list);
	}

	//if use contact map constraint
	if (stoi(CONFIG["cm_cons"])){
		//init cm_constraint_files path to contact_map_constraint::CM_list
		getConstraintList(CONFIG["cm_cons_list"], CM_list);
	}

	//if use initial structure constraint
	if (stoi(CONFIG["init_cons"])){
		//init initial_constraint_files path to initial_constraint::INIT_list
		getConstraintList(CONFIG["init_cons_list"], INIT_list);
	}

	//if add side chain
	map<string, vector<RotamerFeature>> RotamerData;
	if (stoi(CONFIG["sc_cons"])){
		RotamerData = m.getRotamerData();
	}

	time_t star_time = time(NULL);

	map<string, string> files = getFastaFile(CONFIG["input_list"]);
	for (auto &file : files){

		cout << "Modeling " << file.first <<  endl;
		cout << file.second << endl;

		string output_file = CONFIG["output_dir"] + "/fold_" + file.first;

		//construct residuesData from fasta file
		vector<Residue> residuesData = fastaToResiduesData(file.second);

		//read csf data from mongo
		vector<map<string, vector<CSFFeature>>> CSFData = m.getCSFData(residuesData);

		//cal coverage
		double total = 0.0;
		for_each(CSFData.begin(), CSFData.end(), [&](map<string, vector<CSFFeature>> &x){total += x.size(); });
		double coverage = total / CSFData.size() / residuesData.size();
		cout << coverage << endl;

		if (coverage < 0.8 && stoi(CONFIG["init_cons"])){
			CONFIG["init_cons_dynamic"] = "1";
		}else{
			CONFIG["init_cons_dynamic"] = "0";
		}

		//read TA data from mongo
		vector<map<string, vector<TAFeature>>> TAData = m.getTAData(residuesData);
		vector<vector<TAFeature>> possiblePhipsis = getPossiblePhiPsis(residuesData, TAData);

		//if use torsion angles constraint
		if (stoi(CONFIG["ta_cons"])){
			//use # or // to exclude the irrelevant
			vector<Vector2d> phipsi_cons = readTAConstraintFiles(file.first);
			assert(phipsi_cons.size() == possiblePhipsis.size());
			//phipsi_cons->TA
			addTAConsToPossiblePhipsis(possiblePhipsis, phipsi_cons);
		}

		//if use theta tau constraint
		vector<Vector2d> theta_tau_cons;
		if (stoi(CONFIG["tt_cons"])){
			//use # or // to exclude the irrelevant
			theta_tau_cons = readTTConstraintFiles(file.first);
			assert(theta_tau_cons.size() == possiblePhipsis.size());
		}

		//if use dssp8 and asa constraints
		vector<DSSPInfo> dssp_cons;
		if (stoi(CONFIG["dssp_cons"])){
			//use # or // to exclude the irrelevant
			dssp_cons = readDSSPConstraintFiles(file.first);
			assert(dssp_cons.size() == possiblePhipsis.size());
		}

		//if use contact map constraint
		CMInfo cm_cons;
		if (stoi(CONFIG["cm_cons"])){
			//use # or // to exclude the irrelevant, use 0 8
			cm_cons = readCMConstraintFiles(file.first, residuesData.size());
		}

		//if add side chain
		vector<map<string, vector<DASFFeature>>> DASFData;
		if (stoi(CONFIG["sc_cons"])){
			DASFData = m.getDASFData(residuesData);
		}
		
		//if use initial structure constraint
		vector<Atom> atomsData_cons;
		if (stoi(CONFIG["init_cons_dynamic"])){
			//read pdb
			atomsData_cons = readPDB(INIT_list.at(file.first));
			bool MAIN_CHAIN = true;
			vector<Residue> residuesData_cons = atomToResidue(atomsData_cons, MAIN_CHAIN);
			bool FROM_MODEL = true;
			getResidueInfo(residuesData_cons, FROM_MODEL);

			//phipsi_cons->init_cons->TA
			addInitConsToPossiblePhipsis(possiblePhipsis, residuesData_cons);

			assert(residuesData.size() == residuesData_cons.size());
			residuesData = residuesData_cons;
		}else{
			bool FROM_MODEL = false;
			getResidueInfo(residuesData, FROM_MODEL);
		}
		
		residuesData = optimize_mt(possiblePhipsis, cm_cons, CSFData, coverage, atomsData_cons,
			file.second, RotamerData, DASFData, dssp_cons, theta_tau_cons);

		if (stoi(CONFIG["sc_cons"])){
			outputPDB(output_file, residuesData, true);
		}else{
			outputPDB(output_file, residuesData);
		}
		
	
		time_t end_time = time(NULL);
		cout << "Total time: " << (end_time - star_time) << "s..." << endl;		
		
		cout << file.first << " done..." << endl;
	}

	//TODO: too lazy to delete residue.geo

	time_t end_time = time(NULL);
	cout << "Total time: " << (end_time - star_time) << "s..." << endl;

	m.destroy();

 	cout << "Mission complete..." << endl;

	return 0;
}