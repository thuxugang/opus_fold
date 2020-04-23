#include "other_potentials.h"

map<string, string> DSSP_list;

DSSPInfo::DSSPInfo(){
}

DSSPInfo::DSSPInfo(double asa, vector<double>& results){
	this->m_asa = asa;
	this->m_dssp_dict['C'] = results[0];
	this->m_dssp_dict['S'] = results[1];
	this->m_dssp_dict['T'] = results[2];
	this->m_dssp_dict['H'] = results[3];
	this->m_dssp_dict['G'] = results[4];
	this->m_dssp_dict['I'] = results[5];
	this->m_dssp_dict['E'] = results[6];
	this->m_dssp_dict['B'] = results[7];
}

DSSPInfo::DSSPInfo(char dssp){
	this->m_asa = 0;
	this->m_dssp = dssp;
}

DSSPInfo::DSSPInfo(double asa, char dssp){
	this->m_asa = asa;
	this->m_dssp = dssp;
}

vector<DSSPInfo> readDSSPConstraintFiles(const string& name){

	vector<DSSPInfo> dssp_cons;

	stringstream ss;

	ifstream dsspFile;
	try{
		dsspFile.open(DSSP_list.at(name));
	}catch (...){
		cout << "DSSP list format error!";
		dsspFile.close();
		exit(-1);
	}

	string str_line;
	try{
		while (!dsspFile.eof()){
			getline(dsspFile, str_line);
			str_line = trim(str_line);
			if (str_line.find('#') == 0 || str_line.find('=') == 0 || str_line.size() == 0){
				continue;
			}
			char spliter = ' ';
			if (str_line.find('\t') != str_line.npos){
				spliter = '\t';
			}
			vector<double> items = splitd(ss, str_line, spliter);
			double asa = stod(CONFIG["row_asa"]);
			vector<double> dssp8(items.begin() + stod(CONFIG["row_dssp8_start"]), items.end());
			assert(dssp8.size() == 8);
			dssp_cons.push_back(DSSPInfo(asa, dssp8));
		}
	}catch (...){
		cout << "DSSP file format error!";
		dsspFile.close();
		exit(-1);
	}
	dsspFile.close();

	return dssp_cons;
}

vector<DSSPInfo> getDSSPResults(const vector<Atom>& atomsData, bool use_sidechain){

	vector<OPUSAtom> opusAtoms;
	for (Atom atom : atomsData){
		string triResname = getTriResname(atom.m_resname);
		opusAtoms.emplace_back(atom.m_id, atom.m_name, triResname, atom.m_resid,
			atom.m_position[0], atom.m_position[1], atom.m_position[2]);
	}

	vector<DSSPInfo> dssp_results;
	if (use_sidechain){
		vector<DSSPOriResults> results = getDSSP8AndASAResults(opusAtoms);
		for (DSSPOriResults result : results){
			dssp_results.emplace_back(result.m_asa, result.m_dssp8);
		}
	}else{
		vector<char> results = getDSSP8Results(opusAtoms);
		for (char result : results){
			dssp_results.emplace_back(result);
		}
	}

	return dssp_results;
}

void get_KORPE(vector<double>& other_potentials){

	system(("./other_potentials/korpe/korpe " + CONFIG["tmp_dir"] + "/tmp.txt \
			--score_file ./other_potentials/korpe/korp6Dv1.bin \
			-o ./other_potentials/korpe/korp6D \
			> ./other_potentials/korpe/korp6D.log").c_str());

	stringstream ss;

	ifstream korpeFile;
	try{
		korpeFile.open("./other_potentials/korpe/korp6D_score.txt");
	}catch (...){
		cout << "korp6D_score file error!";
		korpeFile.close();
		exit(-1);
	}

	string str_line;
	try{
		while (!korpeFile.eof()){
			getline(korpeFile, str_line);
			str_line = trim(str_line);
			if (str_line.find('#') == 0 || str_line.find('=') == 0 || str_line.size() == 0){
				continue;
			}
			char spliter = ' ';
			if (str_line.find('\t') != str_line.npos){
				spliter = '\t';
			}
			vector<double> items = splitd(ss, str_line, spliter);
			other_potentials.push_back(stof(CONFIG["w_korpe"])*items[1]);
		}
	}catch (...){
		cout << "korp6D_score format error!";
		korpeFile.close();
		exit(-1);
	}
	korpeFile.close();
}

void get_ITDA(vector<double>& other_potentials, int ga_population){

	system("cd ./other_potentials/ITDA/ITDA_V1.2 && \
		   cp -r ../../../tmp_structures/*.pdb ./ && \
		   ./ITDA -list itda_list.txt > itda_scores.txt");

	ifstream itdaFile;
	try{
		itdaFile.open("./other_potentials/ITDA/ITDA_V1.2/itda_scores.txt");
	}
	catch (...){
		cout << "ITDA file error!";
		itdaFile.close();
		exit(-1);
	}
	stringstream ss;
	string str_line;
	try{
		while (!itdaFile.eof()){
			getline(itdaFile, str_line);
			str_line = trim(str_line);
			if (str_line.find('#') == 0 || str_line.find('=') == 0 || str_line.size() == 0){
				continue;
			}
			char spliter = ':';
			vector<double> items = splitd(ss, str_line, spliter);
			other_potentials.push_back(stof(CONFIG["w_itda"])*items[1]);
		}
	}catch (...){
		cout << "itdaFile_scores format error!";
		itdaFile.close();
		exit(-1);
	}
	itdaFile.close();
}

void get_ANDIS(vector<double>& other_potentials, int ga_population){

	system("cd ./other_potentials/ANDIS/ANDIS && \
		   	./ANDIS ../../../tmp_structures/ andis_list.txt > andis_scores.txt");

	ifstream andisFile;
	try{
		andisFile.open("./other_potentials/ANDIS/ANDIS/andis_scores.txt");
	}catch (...){
		cout << "ANDIS file error!";
		andisFile.close();
		exit(-1);
	}
	stringstream ss;
	string str_line;
	try{
		while (!andisFile.eof()){
			getline(andisFile, str_line);
			str_line = trim(str_line);
			if (str_line.find('#') == 0 || str_line.find('*') == 0 || str_line.size() == 0){
				continue;
			}
			char spliter = ' ';
			if (str_line.find('\t') != str_line.npos){
				spliter = '\t';
			}
			vector<double> items = splitd(ss, str_line, spliter);
			other_potentials.push_back(stof(CONFIG["w_andis"])*items.back());
		}
	}catch (...){
		cout << "andisFile format error!";
		andisFile.close();
		exit(-1);
	}
	andisFile.close();
}

vector<double> getOtherPotentials(int ga_population){

	vector<double> other_potentials;
	other_potentials.assign(ga_population, 0);

	if (stoi(CONFIG["use_korpe"])){
		vector<double> other_potential;
		get_KORPE(other_potential);
		assert(other_potential.size() == ga_population);
		for (int i = 0; i < ga_population; ++i){
			other_potentials[i] += other_potential[i];
		}
	}
	if (stoi(CONFIG["use_itda"])){
		vector<double> other_potential;
		get_ITDA(other_potential, ga_population);
		assert(other_potential.size() == ga_population);
		for (int i = 0; i < ga_population; ++i){
			other_potentials[i] += other_potential[i];
		}
	}
	if (stoi(CONFIG["use_andis"])){
		vector<double> other_potential;
		get_ANDIS(other_potential, ga_population);
		assert(other_potential.size() == ga_population);
		for (int i = 0; i < ga_population; ++i){
			other_potentials[i] += other_potential[i];
		}
	}
	return other_potentials;
}