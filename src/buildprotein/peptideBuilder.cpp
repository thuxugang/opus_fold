#include "angle.h"
#include "peptideBuilder.h"


Geo* getGeo(int resid, char resname, int* totalAtoms){
	Geo* geo;
	switch (resname){
		case 'G':
			geo = new GlyGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 4;
			break;
		case 'A':
			geo = new AlaGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 5;
			break;
		case 'S':
			geo = new SerGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 6;
			break;
		case 'C':
			geo = new CysGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 6;
			break;
		case 'V':
			geo = new ValGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 7;
			break;
		case 'I':
			geo = new IleGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 8;
			break;
		case 'L':
			geo = new LeuGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 8;
			break;
		case 'T':
			geo = new ThrGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 7;
			break;
		case 'R':
			geo = new ArgGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 11;
			break;
		case 'K':
			geo = new LysGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 9;
			break;
		case 'D':
			geo = new AspGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 8;
			break;
		case 'N':
			geo = new AsnGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 8;
			break;
		case 'E':
			geo = new GluGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 9;
			break;
		case 'Q':
			geo = new GlnGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 9;
			break;
		case 'M':
			geo = new MetGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 8;
			break;
		case 'H':
			geo = new HisGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 10;
			break;
		case 'P':
			geo = new ProGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 7;
			break;
		case 'F':
			geo = new PheGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 11;
			break;
		case 'Y':
			geo = new TyrGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 12;
			break;
		case 'W':
			geo = new TrpGeo(resid, *totalAtoms);
			*totalAtoms = *totalAtoms + 14;
			break;
		default:
			throw "peptideBuilder.getGeo() wrong!";
	}
	return geo;
}

//vector<Residue> initResidueChain(const vector<char>& resnames){
//	vector<Residue> residuesData;
//	int length = resnames.size();
//	for(int i=0; i<length; i++){
//		Residue aa = Residue(i+1,resnames[i]);
//		residuesData.push_back(aa);
//	}
//	return residuesData;
//}


//
//vector<Geo*> initGeoChain(vector<Residue>& residuesData){
//	vector<Geo*> geosData;
//	int totalAtoms = 1;
//	int length = residuesData.size();
//	for(int i=0; i<length; i++){
//		Geo* geo = getGeo(residuesData[i].m_resid,residuesData[i].m_resname, &totalAtoms);
//		geo->copyMainChain(residuesData[i]);
//		geo->addCB(residuesData[i]);
//		geosData.push_back(geo);
//	}
//	return geosData;
//}
//
//vector<Atom> rebuildSideChain(vector<Geo*>& geosData, vector<Residue>& residuesData){
//	vector<Atom> atomsData;
//	int length = geosData.size();
//	for(int i=0; i<length; i++){
//		//cout << residuesData[i].m_resid << "\t" << residuesData[i].m_resname << endl;
//		atomsData.push_back(residuesData[i].getAtom("N"));
//		atomsData.push_back(residuesData[i].getAtom("CA"));
//		atomsData.push_back(residuesData[i].getAtom("C"));
//		atomsData.push_back(residuesData[i].getAtom("O"));
//		if(residuesData[i].m_resname != 'G'){
//			atomsData.push_back(residuesData[i].getAtom("CB"));
//			if(residuesData[i].m_resname != 'A'){
//				int choose_rl = residuesData[i].m_rl_id;
//				int atom_length = residuesData[i].m_rotamerlibs[choose_rl].m_atoms.size();
//				for(int j=0; j < atom_length; j++){
//					atomsData.push_back(residuesData[i].m_rotamerlibs[choose_rl].m_atoms[j]);
//				}
//			}
//		}
//
//	}
//	return atomsData;
//}

//void rebuildSideChain(vector<Geo*>& geosData, const vector<Residue>& residuesData){
//
//	int length = geosData.size();
//	for(int i=0; i<length; i++){
//		vector<double> dihedrals = residuesData[i].m_dihedrals;
//		geosData[i]->setRotamers(dihedrals);
//		geosData[i]->addSideChain();
//	}
//}
//
//vector<Atom> outputAtomsData(const vector<Geo*>& geosData){
//	vector<Atom> atomsData;
//	int length = geosData.size();
//	for(int i=0; i<length; i++){
//		geosData[i]->output(atomsData);
//	}
//	return atomsData;
//}
//
//



//================reconstruct================
//CB new
//void reconstructFromGeo(vector<Residue>& residuesData){
//	//reconstruct from residue.geo
//	int length = residuesData.size();
//	Geo* geo;
//	Geo* _geo;
//	GeoContainCB* geo_cb;
//	for (int i = 0; i < length; ++i){
//		geo = residuesData[i].m_geo;
//		//set first
//		if (i == 0){
//			//N
//			geo->N.m_position = Vector3d(geo->CA_N_length*cos(geo->N_CA_C_angle*(M_PI / 180.0)),
//				geo->CA_N_length*sin(geo->N_CA_C_angle*(M_PI / 180.0)), 0);
//			//CA
//			geo->CA.m_position = Vector3d(0, 0, 0);
//			//C
//			geo->C.m_position = Vector3d(geo->CA_C_length, 0, 0);
//		}else{
//			_geo = residuesData[i - 1].m_geo;
//			//N
//			geo->N.m_position = calCoordinates(_geo->N, _geo->CA, _geo->C,
//				geo->_peptide_bond, geo->_CA_C_N_angle, _geo->psi);
//			//CA
//			geo->CA.m_position = calCoordinates(_geo->CA, _geo->C, geo->N,
//				geo->CA_N_length, geo->_C_N_CA_angle, geo->omega);
//			//C
//			geo->C.m_position = calCoordinates(_geo->C, geo->N, geo->CA,
//				geo->CA_C_length, geo->N_CA_C_angle, geo->phi);
//
//			//last one
//			if (i == length - 1){
//				//O
//				geo->O.m_position = calCoordinates(geo->N, geo->CA, geo->C,
//					geo->C_O_length, geo->CA_C_O_angle, geo->N_CA_C_O_diangle);
//
//				//CB assume 120
//				if (residuesData[i].m_resname != 'G'){
//					geo_cb = (GeoContainCB*)geo;
//					geo_cb->CB.m_position = calOCoordinates(geo_cb->N, geo_cb->CA, geo_cb->C, geo_cb->CA_CB_length);
//				}
//
//			}
//		}
//	}
//
//	//Need next atom
//	Geo* geo_;
//	for (int i = 0; i < length; ++i){
//		geo = residuesData[i].m_geo;
//		//set first
//		if (i == 0){
//			//O
//			geo_ = residuesData[i + 1].m_geo;
//			geo->O.m_position = calOCoordinates(geo->CA, geo->C, geo_->N, geo->C_O_length);
//			//CB backward
//			if (residuesData[i].m_resname != 'G'){
//				geo_cb = (GeoContainCB*)geo;
//				geo_cb->CB.m_position = calCoordinates(geo_->N, geo_cb->C, geo_cb->CA,
//					geo_cb->CA_CB_length, geo_cb->C_CA_CB_angle, geo_cb->N_C_CA_CB_diangle);
//			}
//		}
//		else{
//			if (i != length - 1){
//				//O assume 120
//				geo_ = residuesData[i + 1].m_geo;
//				geo->O.m_position = calOCoordinates(geo->CA, geo->C, geo_->N, geo->C_O_length);
//				//CB backward
//				if (residuesData[i].m_resname != 'G'){
//					geo_cb = (GeoContainCB*)geo;
//					geo_cb->CB.m_position = calCoordinates(geo_->N, geo_cb->C, geo_cb->CA,
//						geo_cb->CA_CB_length, geo_cb->C_CA_CB_angle, geo_cb->N_C_CA_CB_diangle);
//				}
//			}
//		}
//
//	}
//}

//CB old (wrong but work)
void reconstructFromGeo(vector<Residue>& residuesData){
	//reconstruct from residue.geo
	int length = residuesData.size();
	Geo* geo;
	Geo* _geo;
	GeoContainCB* geo_cb;
	for (int i = 0; i < length; ++i){
		geo = residuesData[i].m_geo;
		//set first
		if (i == 0){
			//N
			geo->N.m_position = Vector3d(geo->CA_N_length*cos(geo->N_CA_C_angle*(M_PI / 180.0)),
				geo->CA_N_length*sin(geo->N_CA_C_angle*(M_PI / 180.0)), 0);
			//CA
			geo->CA.m_position = Vector3d(0, 0, 0);
			//C
			geo->C.m_position = Vector3d(geo->CA_C_length, 0, 0);
		}else{
			_geo = residuesData[i - 1].m_geo;
			//N
			geo->N.m_position = calCoordinates(_geo->N, _geo->CA, _geo->C,
				geo->_peptide_bond, geo->_CA_C_N_angle, _geo->psi);
			//CA
			geo->CA.m_position = calCoordinates(_geo->CA, _geo->C, geo->N,
				geo->CA_N_length, geo->_C_N_CA_angle, geo->omega);
			//C
			geo->C.m_position = calCoordinates(_geo->C, geo->N, geo->CA,
				geo->CA_C_length, geo->N_CA_C_angle, geo->phi);

			//last one
			if (i == length - 1){
				//O
				geo->O.m_position = calCoordinates(geo->N, geo->CA, geo->C,
					geo->C_O_length, geo->CA_C_O_angle, geo->N_CA_C_O_diangle);

			}
		}

		//CB
		if (residuesData[i].m_resname != 'G'){
			geo_cb = (GeoContainCB*)geo;
			geo_cb->CB.m_position = calCoordinates(geo_cb->N, geo_cb->C, geo_cb->CA,
				geo_cb->CA_CB_length, geo_cb->C_CA_CB_angle, geo_cb->N_C_CA_CB_diangle);
		}
	}

	//Need next atom
	Geo* geo_;
	for (int i = 0; i < length; ++i){
		geo = residuesData[i].m_geo;
		//O assume 120
		if (i != length - 1){
			//O
			geo_ = residuesData[i + 1].m_geo;
			geo->O.m_position = calOCoordinates(geo->CA, geo->C, geo_->N, geo->C_O_length);
		}
	}
}

void reconstructSideChainFromGeo(vector<Residue>& residuesData){
	int length = residuesData.size();
	for (int i = 0; i < length; ++i){
		if (residuesData[i].m_resname == 'G' || residuesData[i].m_resname == 'A'){
			continue;
		}
		residuesData[i].m_geo->setRotamers(residuesData[i].m_geo->dihedral);
		residuesData[i].m_geo->addSideChain();
	}
}
//================reconstruct================

