#include "DASFScore.h"

DASFReference::DASFReference(){
}

DASFReference getDASFReference(const Residue& residue){

	DASFReference reference;

	Vector3d ca_ref = residue.m_geo->CA.m_position;
	Vector3d c_ref = residue.m_geo->C.m_position;
	Vector3d o_ref = residue.m_geo->O.m_position;

	reference.ref = ca_ref;
	Vector3d c_ref_new = c_ref - reference.ref;
	Vector3d o_ref_new = o_ref - reference.ref;

	//c - ca
	Vector3d x_axis = c_ref_new / c_ref_new.norm();
	Vector3d c_o = o_ref_new - c_ref_new;

	//o - c perpendicular to x_axis
	Vector3d y_axis = c_o - (x_axis.dot(c_o) / x_axis.dot(x_axis) * x_axis);
	y_axis = y_axis / y_axis.norm();

	Vector3d z_axis = x_axis.cross(y_axis);

	reference.rotation_matrix << x_axis[0], y_axis[0], z_axis[0], 0,
		x_axis[1], y_axis[1], z_axis[1], 0,
		x_axis[2], y_axis[2], z_axis[2], 0,
		0, 0, 0, 1;

	return reference;
}

Vector3d transCoordinate(const DASFReference& reference, const Vector3d& position){
	Vector3d position_new = position - reference.ref;
	Vector4d position_new4d(position_new[0], position_new[1], position_new[2], 1);
	Vector4d r = position_new4d.transpose()*reference.rotation_matrix;
	return Vector3d(r[0], r[1], r[2]);
}

vector<Vector3d> getFeatures(const Residue& residue){
	vector<Vector3d> fs;
	switch (residue.m_resname){
		case 'V':
			fs.push_back(residue.m_geo->m_geo_atoms["CG1"]->m_position);	//"CG1"
			break;
		case 'I':
			fs.push_back(residue.m_geo->m_geo_atoms["CG1"]->m_position);	//"CG1"
			fs.push_back(residue.m_geo->m_geo_atoms["CD1"]->m_position);	//"CD1"
			break;
		case 'L':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["CD1"]->m_position);	//"CD1"
			break;
		case 'S':
			fs.push_back(residue.m_geo->m_geo_atoms["OG"]->m_position);	//"OG"
			break;
		case 'T':
			fs.push_back(residue.m_geo->m_geo_atoms["OG1"]->m_position);	//"OG1"
			break;                            
		case 'D':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["OD1"]->m_position);	//"OD1"
			break;
		case 'N':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["OD1"]->m_position);	//"OD1"
			break;
		case 'E':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["CD"]->m_position);	//"CD"
			fs.push_back(residue.m_geo->m_geo_atoms["OE1"]->m_position);	//"OE1"
			break;
		case 'Q':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["CD"]->m_position);	//"CD"
			fs.push_back(residue.m_geo->m_geo_atoms["OE1"]->m_position);	//"OE1"
			break;
		case 'K':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["CD"]->m_position);	//"CD"
			fs.push_back(residue.m_geo->m_geo_atoms["CE"]->m_position);	//"CE"
			fs.push_back(residue.m_geo->m_geo_atoms["NZ"]->m_position);	//"NZ"
			break;
		case 'R':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["CD"]->m_position);	//"CD"
			fs.push_back(residue.m_geo->m_geo_atoms["NE"]->m_position);	//"NE"
			fs.push_back(residue.m_geo->m_geo_atoms["CZ"]->m_position);	//"CZ"
			break;
		case 'C':
			fs.push_back(residue.m_geo->m_geo_atoms["SG"]->m_position);	//"SG"
			break;
		case 'M':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["SD"]->m_position);	//"SD"
			fs.push_back(residue.m_geo->m_geo_atoms["CE"]->m_position);	//"CE"
			break;
		case 'F':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["CD1"]->m_position);	//"CD1"
			break;
		case 'Y':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["CD1"]->m_position);	//"CD1"
			break;
		case 'W':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["CD1"]->m_position);	//"CD1"
			break;
		case 'H':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["ND1"]->m_position);	//"ND1"
			break;
		case 'P':
			fs.push_back(residue.m_geo->m_geo_atoms["CG"]->m_position);	//"CG"
			fs.push_back(residue.m_geo->m_geo_atoms["CD"]->m_position);	//"CD"
			break;
	}
	return fs;
}

double getDASFPotential(const Residue& residue, const DASFReference& reference){
	double results = 0;
	if (residue.m_resname != 'G' && residue.m_resname != 'A'){
		if (residue.m_geo->standerDASFs.size() == 0){
			return results;
		}
		vector<Vector3d> fs = getFeatures(residue);

		//5,7,9,11
		for (vector<DASFFeature> standerDASF : residue.m_geo->standerDASFs){
			assert(fs.size() == standerDASF.size());
			int fs_length = fs.size();
			//different atoms
			for (int i = 0; i < fs_length; ++i){
				Vector3d f = transCoordinate(reference, fs[i]);
				double result = standerDASF[0].m_window_len/10.0*
					(((f - standerDASF[i].m_mean).array() / standerDASF[i].m_sd.array()).cwiseAbs().sum());
				results = results + result;				
			}		
		}
		//normal
		//results = results / fs.size() / residue.m_geo->standerDASFs.size();
	}
	return results;
}