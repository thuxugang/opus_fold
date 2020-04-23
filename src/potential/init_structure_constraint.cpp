#include "init_structure_constraint.h"

map<string, string> INIT_list;

void addInitConsToPossiblePhipsis(vector<vector<TAFeature>>& possiblePhipsis, const vector<Residue>& residuesData_cons){

	int length = possiblePhipsis.size();
	for (int i = 0; i < length; ++i){
		vector<double> mean;
		if (i == 0 || i == (length-1)){
			mean = { possiblePhipsis[i][0].m_mean[0], possiblePhipsis[i][0].m_mean[1] };
		}else{
			mean = { residuesData_cons[i].m_geo->phi, residuesData_cons[i].m_geo->psi };
		}
		
		vector<double> sd = { 10, 0, 0, 10 };
		double possibility = possiblePhipsis[i][0].m_possibility;
		int window_len = possiblePhipsis[i][0].m_window_len;
		possiblePhipsis[i].insert(possiblePhipsis[i].begin()+1, TAFeature(mean, sd, possibility, window_len));
		possiblePhipsis[i][0].m_total_possibility += possibility;
	}
}

