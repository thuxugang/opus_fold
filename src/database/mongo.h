#ifndef MONGO_H
#define MONGO_H

#include <map>  
#include <string>
#include <iostream>
//#include <winsock2.h>
#include <sstream>
#include <vector>
#include <Eigen/Dense> 

#include "config.h"
#include "residue.h"
#include "basic.h"
#include "myutils.h"

#include <bsoncxx/builder/stream/document.hpp>       
#include <bsoncxx/json.hpp>
#include <mongocxx/client.hpp>
#include <mongocxx/stdx.hpp>
#include <mongocxx/uri.hpp>
#include <mongocxx/instance.hpp>

using namespace std;
using namespace Eigen;

class Mongo{
public:

	mongocxx::instance instance{};
	mongocxx::client c{ mongocxx::uri{} };

	stringstream ss;

	Mongo();

	vector<map<string, vector<CSFFeature>>> getCSFData(const vector<Residue>& residuesData);
	vector<map<string, vector<TAFeature>>> getTAData(const vector<Residue>& residuesData);
	vector<map<string, vector<DASFFeature>>> getDASFData(const vector<Residue>& residuesData);
	
	map<string, vector<RotamerFeature>> getRotamerData();

	void destroy();

private:
	vector<CSFFeature> getCSFResult(const string& collection_name, const string& key, int window_len);
	void getCollectionData(vector<map<string, vector<CSFFeature>>>& CSFResults, 
		const vector<Residue>& residuesData, const string& collection_name, int window_len);

	vector<TAFeature> getTAResult(const string& collection_name, const string& key, int window_len);
	void getCollectionData(vector<map<string, vector<TAFeature>>>& TAResults,
		const vector<Residue>& residuesData, const string& collection_name, int window_len);

	vector<DASFFeature> getDASFResult(const string& collection_name, const string& key, int window_len);
	void getCollectionData(vector<map<string, vector<DASFFeature>>>& DASFResults,
		const vector<Residue>& residuesData, const string& collection_name, int window_len);

	void getRotamerResult(map<string, vector<RotamerFeature>>& RotamerResults,
		const string& collection_name);

};


#endif