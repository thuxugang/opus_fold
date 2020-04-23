#ifndef CONFIG_H
#define CONFIG_H

#include <map>  
#include <string> 
#include <random>
#include <vector> 
#include <fstream> 
#include <iostream> 
//#include <direct.h>
//#include <io.h>

using namespace std;

extern map<string, string> CONFIG;

string& trim(string &s);

void readConfigFile(const string& path);

//vector<string> getInputFiles(const string& path);

map<string, string> getFastaFile(const string& path);

//init constraint_files path to list
void getConstraintList(string& path, map<string, string>& list);

#endif

