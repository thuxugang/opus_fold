#include "mongo.h"

Mongo::Mongo(){
}

void Mongo::destroy(){
}

//==============CSFData==============
vector<CSFFeature> Mongo::getCSFResult(const string& collection_name, const string& key, int window_len){

	vector<CSFFeature> fs;
	vector<string> db_col = splits(ss, collection_name, '.');
	try{
		mongocxx::collection coll = c[db_col[0]][db_col[1]];
		bsoncxx::stdx::optional<bsoncxx::document::value> result =
			coll.find_one(bsoncxx::builder::stream::document{} << "key" << key << bsoncxx::builder::stream::finalize);
		if (result){
			bsoncxx::document::view view = (*result).view();
			string result = view["value"].get_utf8().value.to_string();

			vector<string> features = splits(ss, result, ';');
			int length_features = features.size();
			for (int i = 0; i<length_features; ++i){
				vector<string> feature = splits(ss, features[i], '#');
				fs.emplace_back(splitd(ss, feature[0], '_'), splitd(ss, feature[1], '_'), atoi(feature[2].c_str()), window_len);
			}
		}
	}catch (...){
		cout << "CSF query failed!" << endl;
		throw "CSF query failed!";
	}
	
	return fs;

}

void Mongo::getCollectionData(vector<map<string, vector<CSFFeature>>>& CSFResults, 
	const vector<Residue>& residuesData, const string& collection_name, int window_len){

	map <string, vector<CSFFeature>> collectionsResults;
	int length = residuesData.size();
	for (int i = 0; i <= (length - window_len); ++i){
		string key;
		for (int j = 0; j < window_len; ++j){
			key += residuesData[i + j].m_resname;
		}
		vector<CSFFeature> fs = getCSFResult(collection_name, key, window_len);
		if (fs.size() != 0){
			collectionsResults[key] = fs;
		}
	}
	CSFResults.push_back(collectionsResults);
}

vector<map<string, vector<CSFFeature>>> Mongo::getCSFData(const vector<Residue>& residuesData){
	
	vector<map<string, vector<CSFFeature>>> CSFResults;
	getCollectionData(CSFResults, residuesData, CONFIG["csf_collection5"], 5);
	getCollectionData(CSFResults, residuesData, CONFIG["csf_collection7"], 7);
	getCollectionData(CSFResults, residuesData, CONFIG["csf_collection9"], 9);
	getCollectionData(CSFResults, residuesData, CONFIG["csf_collection11"], 11);
	ss.str("");
	return CSFResults;
}
//==============CSFData==============


//==============TAData==============
vector<TAFeature> Mongo::getTAResult(const string& collection_name, const string& key, int window_len){

	vector<TAFeature> fs;

	vector<string> db_col = splits(ss, collection_name, '.');
	try{
		mongocxx::collection coll = c[db_col[0]][db_col[1]];
		bsoncxx::stdx::optional<bsoncxx::document::value> result =
			coll.find_one(bsoncxx::builder::stream::document{} << "key" << key << bsoncxx::builder::stream::finalize);
		if (result){
			bsoncxx::document::view view = (*result).view();
			string result = view["value"].get_utf8().value.to_string();

			vector<string> features = splits(ss, result, ';');
			int length_features = features.size();
			for (int i = 0; i < length_features; ++i){
				vector<string> feature = splits(ss, features[i], '#');
				fs.emplace_back(splitd(ss, feature[0], '_'), splitd(ss, feature[1], '_'), atof(feature[2].c_str()), window_len);
			}
		}
	}catch (...){
		cout << "TA query failed!" << endl;
		throw "TA query failed!";
	}

	return fs;

}

void Mongo::getCollectionData(vector<map<string, vector<TAFeature>>>& TAResults,
	const vector<Residue>& residuesData, const string& collection_name, int window_len){

	map <string, vector<TAFeature>> collectionsResults;
	int length = residuesData.size();
	for (int i = 0; i <= (length - window_len); ++i){
		string key;
		for (int j = 0; j < window_len; ++j){
			key += residuesData[i + j].m_resname;
		}
		vector<TAFeature> fs = getTAResult(collection_name, key, window_len);
		if (fs.size() != 0){
			collectionsResults[key] = fs;
		}
	}

	//first and last
	if (window_len == 3){
		string key;
		key += "G";
		key += residuesData[0].m_resname;
		key += residuesData[1].m_resname;
		collectionsResults[key] = getTAResult(collection_name, key, window_len);

		key.clear();
		key += residuesData[residuesData.size() - 2].m_resname;
		key += residuesData[residuesData.size() - 1].m_resname;
		key += "G";
		collectionsResults[key] = getTAResult(collection_name, key, window_len);

	}
	TAResults.push_back(collectionsResults);
}

vector<map<string, vector<TAFeature>>> Mongo::getTAData(const vector<Residue>& residuesData){

	vector<map<string, vector<TAFeature>>> TAResults;
	getCollectionData(TAResults, residuesData, CONFIG["ta_collection3"], 3);
	getCollectionData(TAResults, residuesData, CONFIG["ta_collection5"], 5);
	getCollectionData(TAResults, residuesData, CONFIG["ta_collection7"], 7);
	ss.str("");
	return TAResults;

}
//==============TAData==============

//==============DASFData==============
vector<DASFFeature> Mongo::getDASFResult(const string& collection_name, const string& key, int window_len){

	vector<DASFFeature> fs;
	vector<string> db_col = splits(ss, collection_name, '.');
	try{
		mongocxx::collection coll = c[db_col[0]][db_col[1]];
		bsoncxx::stdx::optional<bsoncxx::document::value> result =
			coll.find_one(bsoncxx::builder::stream::document{} << "key" << key << bsoncxx::builder::stream::finalize);
		if (result){
			bsoncxx::document::view view = (*result).view();
			string result = view["value"].get_utf8().value.to_string();

			vector<string> features = splits(ss, result, ';');
			int length_features = features.size();
			for (int i = 0; i < length_features; ++i){
				vector<string> feature = splits(ss, features[i], '#');
				fs.emplace_back(splitd(ss, feature[0], '_'), splitd(ss, feature[1], '_'), atoi(feature[2].c_str()), window_len);
			}
		}
	}catch (...){
		cout << "DASF query failed!" << endl;
		throw "DASF query failed!";
	}

	return fs;
}

void Mongo::getCollectionData(vector<map<string, vector<DASFFeature>>>& DASFResults,
	const vector<Residue>& residuesData, const string& collection_name, int window_len){

	map <string, vector<DASFFeature>> collectionsResults;
	int length = residuesData.size();
	for (int i = 0; i <= (length - window_len); ++i){
		string key;
		for (int j = 0; j < window_len; ++j){
			key += residuesData[i + j].m_resname;
		}
		vector<DASFFeature> fs = getDASFResult(collection_name, key, window_len);
		if (fs.size() != 0){
			collectionsResults[key] = fs;
		}
	}
	DASFResults.push_back(collectionsResults);
}

vector<map<string, vector<DASFFeature>>> Mongo::getDASFData(const vector<Residue>& residuesData){

	vector<map<string, vector<DASFFeature>>> DASFResults;
	getCollectionData(DASFResults, residuesData, CONFIG["dasf_collection5"], 5);
	getCollectionData(DASFResults, residuesData, CONFIG["dasf_collection7"], 7);
	getCollectionData(DASFResults, residuesData, CONFIG["dasf_collection9"], 9);
	getCollectionData(DASFResults, residuesData, CONFIG["dasf_collection11"], 11);
	ss.str("");
	return DASFResults;
}
//==============DASFData==============

//==============rotamer==============
void Mongo::getRotamerResult(map<string, vector<RotamerFeature>>& RotamerResults,
	const string& collection_name){

	vector<RotamerFeature> fs;

	vector<string> db_col = splits(ss, collection_name, '.');
	mongocxx::collection coll = c[db_col[0]][db_col[1]];
	mongocxx::cursor cursor = coll.find(bsoncxx::builder::stream::document{} << bsoncxx::builder::stream::finalize);

	for (auto doc : cursor){
		string key = doc["key"].get_utf8().value.to_string();
		string result = doc["value"].get_utf8().value.to_string();
		vector<string> rls = splits(ss, result, ';');
		int length_rl = rls.size();
		for (int i = 0; i < length_rl; i++){
			vector<double> rotamers = splitd(ss, rls[i], '_');
			assert(rotamers.size() == 5);
			fs.emplace_back(rotamers[0], rotamers[1], rotamers[2], rotamers[3], rotamers[4]);
		}
		RotamerResults[key] = fs;
		fs.clear();
	}
}

map<string, vector<RotamerFeature>> Mongo::getRotamerData(){

	map<string, vector<RotamerFeature>> RotamerResults;
	getRotamerResult(RotamerResults, CONFIG["rotamer_collection"]);
	ss.str("");
	return RotamerResults;
}

//==============rotamer==============
