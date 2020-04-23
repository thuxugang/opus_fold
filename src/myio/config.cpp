#include "config.h"

using namespace std;

map<string, string> CONFIG;

string& trim(string &s){
	if (s.empty()){
		return s;
	}
	s.erase(0, s.find_first_not_of(" \n\r\t"));
	s.erase(s.find_last_not_of(" \n\r\t") + 1);
	return s;
}

void readConfigFile(const string& path){
	ifstream configFile;
	configFile.open(path.c_str());
	string str_line;
	if (configFile.is_open()){
		while (!configFile.eof()){
			getline(configFile, str_line);
			if (str_line.find('#') == 0 || str_line.find('=') == 0 ){
				continue;
			}
			size_t pos = str_line.find('=');
			string str_key = str_line.substr(0, pos);
			string str_value = str_line.substr(pos + 1);
			CONFIG.insert(pair<string, string>(trim(str_key), trim(str_value)));
		}
	}else{
		cout << "Can not open config file!";
		configFile.close();
		exit(-1);
	}
	configFile.close();
}

map<string, string> getFastaFile(const string& path){
	map<string, string> inputfiles;
	ifstream configFile;
	configFile.open(path.c_str());
	string str_line;
	if (configFile.is_open()){
		while (!configFile.eof()){
			getline(configFile, str_line);
			if (str_line.find('>') == 0){
				string name = str_line.substr(1);
				name = trim(name);
				if (name.compare("") == 0){
					name = "TMP_" + to_string(inputfiles.size()) + ".pdb";
				}
				getline(configFile, str_line);
				inputfiles.insert(pair<string, string>(name, trim(str_line)));
			}
		}
	}
	else{
		cout << "Fasta format error!";
		configFile.close();
		exit(-1);
	}
	configFile.close();
	return inputfiles;
}

void getConstraintList(string& path, map<string, string>& list){

	ifstream file;
	file.open(path.c_str());
	string str_line;
	try{
		while (!file.eof()){
			getline(file, str_line);
			str_line = trim(str_line);
			size_t pos = str_line.find(' ');
			string str_key = str_line.substr(0, pos);
			string str_value = str_line.substr(pos + 1);
			list.insert(pair<string, string>(trim(str_key), trim(str_value)));
		}
	}catch (...){
		cout << path << ": constraint list format error!" << endl;
		file.close();
		exit(-1);
	}
	file.close();

}

//vector<string> getInputFiles(const string& path){
//
//	vector<string> files;
//	long long hFile = 0;
//	struct _finddata_t fileinfo;
//	string p;
//	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1){
//		do{
//			if (!(fileinfo.attrib &  _A_SUBDIR)){
//				files.push_back(fileinfo.name);
//			}
//
//		} while (_findnext(hFile, &fileinfo) == 0);
//		_findclose(hFile);
//	}
//	return files;
//}