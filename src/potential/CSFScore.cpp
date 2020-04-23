#include "CSFScore.h"

struct RotationMatrix{
	Matrix4d rotation_matrix;
	Vector3d ca_ref;
};

inline RotationMatrix getRotationMatrix(const Residue& residue){

	Vector3d ca_ref = residue.m_atoms.at("CA").m_position;
	Vector3d c_ref = residue.m_atoms.at("C").m_position;
	Vector3d o_ref = residue.m_atoms.at("O").m_position;

	Vector3d c_ref_new = c_ref - ca_ref;
	Vector3d o_ref_new = o_ref - ca_ref;

	//c - ca
	Vector3d x_axis = c_ref_new / c_ref_new.norm();
	Vector3d c_o = o_ref_new - c_ref_new;

	//o - c perpendicular to x_axis
	Vector3d y_axis = c_o - (x_axis.dot(c_o) / x_axis.dot(x_axis) * x_axis);
	y_axis = y_axis / y_axis.norm();

	Vector3d z_axis = x_axis.cross(y_axis);

	Matrix4d rotation_matrix;
	rotation_matrix  << x_axis[0], y_axis[0], z_axis[0], 0,
						x_axis[1], y_axis[1], z_axis[1], 0,
						x_axis[2], y_axis[2], z_axis[2], 0,
						0, 0, 0, 1;

	return RotationMatrix{rotation_matrix, ca_ref};
}

inline Vector3d transCoordinate(const RotationMatrix& rm, const Residue& residue){

	Vector3d position_new = residue.m_atoms.at("C").m_position - rm.ca_ref;
	Vector4d position_new4d(position_new[0], position_new[1], position_new[2], 1);
	Vector4d r = position_new4d.transpose()*rm.rotation_matrix;

	return Vector3d(r[0], r[1], r[2]);
}

inline double getZscore(const CSFFeature& stander, const Vector3d& features){

	return ((features - stander.m_mean).array() / stander.m_sd.array()).cwiseAbs().sum();
}

double getPotential(const vector<Residue>& residuesData, const map<string, vector<CSFFeature>>& CSFData, int window_len){
	
	int ref_id = (window_len-1)/2;

	int start_id = 0;
	if (window_len == 7 || window_len == 11){
		start_id = 1;
	}

	double totals = 0;
	int length = residuesData.size();
	for (int i = 0; i <= (length - window_len); ++i){
		
		//key
		string key;
		//double id_sum = 0;
		for (auto it = residuesData.begin() + i; it < (residuesData.begin() + i + window_len); ++it){
			key += it->m_resname;
		//	id_sum += it->m_resid;
		}

		//assure no skip
		//if (id_sum / window_len != residuesData[i+ref_id].m_resid){
		//	continue;
		//}

		vector<CSFFeature> stander;
		auto it = CSFData.find(key);
		if (it == CSFData.end()){
			continue;
		}else{
			stander = it->second;
		}
		
		RotationMatrix rm = getRotationMatrix(residuesData[i + ref_id]);
		vector<Vector3d> features;
		for (int j = start_id; j < window_len; j+=2){
			//exclue ref_id
			if (j != ref_id){
				features.push_back(transCoordinate(rm, residuesData[i + j]));
			}
		}

		assert(stander.size() == features.size());
		for (int k = 0; k < stander.size(); ++k){
			totals += getZscore(stander[k], features[k]);
		}
	}
	return totals;
}

double getCSFPotential(const vector<Residue>& residuesData, const vector<map<string, vector<CSFFeature>>>& CSFData){
	double totals = 0;
	totals += getPotential(residuesData, CSFData[0], 5);
	totals += getPotential(residuesData, CSFData[1], 7);
	totals += getPotential(residuesData, CSFData[2], 9);
	totals += getPotential(residuesData, CSFData[3], 11);
	return totals;
}
