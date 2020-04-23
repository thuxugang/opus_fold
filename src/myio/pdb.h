#ifndef PDB_H
#define PDB_H

#include <string> 
#include <vector> 
#include <fstream> 
#include <iostream> 
#include <sstream>
#include "residue.h"
#include "atom.h"

using namespace std;

vector<Atom> getAtomsData(const vector<Residue>& residuesData, bool use_sidechain=false);

void outputPDB(const string& outfile, const vector<Residue>& residuesData, bool use_sidechain=false);

vector<Atom> readPDB(const string& input_file);

#endif

