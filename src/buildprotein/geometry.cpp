#include "geometry.h"
#include "angle.h"

void Geo::addSideChain(){
}

void Geo::setRotamers(const vector<double>& rotamers){
}

void Geo::initAtomInGeo(){
	map<string, Atom*>::iterator iter;
	for (iter = this->m_geo_atoms.begin(); iter != this->m_geo_atoms.end(); ++iter) {
		//for side-chain atoms
		if (this->m_geo_atoms[iter->first]->m_id != -1){
			continue;
		}
		this->m_geo_atoms[iter->first]->m_resname = this->m_resname;
		this->m_geo_atoms[iter->first]->m_resid = this->m_resid;
		this->m_geo_atoms[iter->first]->m_id = this->m_atomid;
		this->m_geo_atoms[iter->first]->setLJParams();

		if (iter->first.compare("N") == 0 || iter->first.compare("CA") == 0 || iter->first.compare("C") == 0 
			|| iter->first.compare("O") == 0 || iter->first.compare("CB") == 0){
			this->m_geo_atoms[iter->first]->m_isMainChain = true;
		}else{
			this->m_geo_atoms[iter->first]->m_isMainChain = false;
		}
	}
}

void Geo::addMainChainAtomsToMap(){
	this->m_geo_atoms["N"] = &this->N;
	this->m_geo_atoms["CA"] = &this->CA;
	this->m_geo_atoms["C"] = &this->C;
	this->m_geo_atoms["O"] = &this->O;
}

void GeoContainCB::addMainChainAtomsToMap(){
	this->Geo::addMainChainAtomsToMap();

	this->m_geo_atoms["CB"] = &this->CB;
}

inline void outputAtom(vector<Atom>& atomsData, Atom& atom, int& counter){
	//add counter to init geo.m_atomid
	atom.m_id += counter;
	counter++;
	atomsData.push_back(atom);
}

void Geo::output(vector<Atom>& atomsData){
	outputAtom(atomsData, N, this->m_id_counter);
	outputAtom(atomsData, CA, this->m_id_counter);
	outputAtom(atomsData, C, this->m_id_counter);
	outputAtom(atomsData, O, this->m_id_counter);
}

void GeoContainCB::output(vector<Atom>& atomsData){
	this->Geo::output(atomsData);
	outputAtom(atomsData, CB, this->m_id_counter);
}

void GeoContainSideChain::output_sidechain(vector<Atom>& atomsData){
}

GlyGeo::GlyGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.8914;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5117;
	N_CA_C_O_diangle = 180.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	m_resname = 'G';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 2.5;
	m_num_dasfs = 0;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}

AlaGeo::AlaGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.068;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5;
	N_CA_C_O_diangle = -60.5;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6860;

	m_resname = 'A';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.7;
	m_num_dasfs = 0;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();

}

SerGeo::SerGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.2812;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5;
	N_CA_C_O_diangle = -60.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;
    

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6618;

	CB_OG_length = 1.417;
	CA_CB_OG_angle = 110.773;
	N_CA_CB_OG_diangle = -63.3;

	m_resname = 'S';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.7;
	m_num_dasfs = 1;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();

}
void SerGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_OG_diangle = rotamers[0];
	}catch(...){
		cout << "SerGeo rotamers list: not long enough" << endl;
	}
}
void SerGeo::addSideChain(){
	this->OG.m_position = calCoordinates(N, CA, CB, CB_OG_length, CA_CB_OG_angle, N_CA_CB_OG_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["OG"] = &this->OG;
		this->initAtomInGeo();
	}
}
void SerGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, OG, this->m_id_counter);
}

CysGeo::CysGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.8856;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5;
	N_CA_C_O_diangle = -60.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;


	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.5037;

	CB_SG_length = 1.808;
	CA_CB_SG_angle = 113.8169;
	N_CA_CB_SG_diangle = -62.2;

	m_resname = 'C';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.4;
	m_num_dasfs = 1;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();

}
void CysGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_SG_diangle = rotamers[0];
	}catch(...){
		cout << "Geo rotamers list: not long enough" << endl;
	}
}
void CysGeo::addSideChain(){
	this->SG.m_position = calCoordinates(N, CA, CB, CB_SG_length, CA_CB_SG_angle, N_CA_CB_SG_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["SG"] = &this->SG;
		this->initAtomInGeo();
	}
}
void CysGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, SG, this->m_id_counter);
}
                                       
ValGeo::ValGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 109.7698;
        
	C_O_length = 1.23;
	CA_C_O_angle = 120.5686;
	N_CA_C_O_diangle = -60.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 123.2347;

	CB_CG1_length = 1.527;
	CA_CB_CG1_angle = 110.7;
	N_CA_CB_CG1_diangle = 177.2;

	CB_CG2_length = 1.527;
	CG1_CB_CG2_angle = 109.50; 
	CA_CG1_CB_CG2_diangle = -120;

	m_resname = 'V';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.5;
	m_num_dasfs = 1;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void ValGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG1_diangle = rotamers[0];
	}catch(...){
		cout << "ValGeo rotamers list: not long enough" << endl;
	}
}
void ValGeo::addSideChain(){
	this->CG1.m_position = calCoordinates(N, CA, CB, CB_CG1_length, CA_CB_CG1_angle, N_CA_CB_CG1_diangle);
	this->CG2.m_position = calCoordinates(CA, CG1, CB, CB_CG2_length, CG1_CB_CG2_angle, CA_CG1_CB_CG2_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG1"] = &this->CG1;
		this->m_geo_atoms["CG2"] = &this->CG2;
		this->initAtomInGeo();
	}
}
void ValGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG1, this->m_id_counter);
	outputAtom(atomsData, CG2, this->m_id_counter);
}

IleGeo::IleGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 109.7202;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5403;
	N_CA_C_O_diangle = -60.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 123.2347;

	CB_CG1_length = 1.527;
	CA_CB_CG1_angle = 110.7;
	N_CA_CB_CG1_diangle = 59.7;

	CB_CG2_length = 1.527;
	CG1_CB_CG2_angle = 109.50;
	CA_CG1_CB_CG2_diangle = 120;

	CG1_CD1_length = 1.52;
	CB_CG1_CD1_angle = 113.97;
	CA_CB_CG1_CD1_diangle = 169.8;

	m_resname = 'I';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.5;
	m_num_dasfs = 2;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void IleGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG1_diangle = rotamers[0];
		CA_CB_CG1_CD1_diangle = rotamers[1];
	}catch(...){
		cout << "IleGeo rotamers list: not long enough" << endl;
	}
}
void IleGeo::addSideChain(){
	this->CG1.m_position = calCoordinates(N, CA, CB, CB_CG1_length, CA_CB_CG1_angle, N_CA_CB_CG1_diangle);	
	this->CG2.m_position = calCoordinates(CA, CG1, CB, CB_CG2_length, CG1_CB_CG2_angle, CA_CG1_CB_CG2_diangle);	
	this->CD1.m_position = calCoordinates(CA, CB, CG1, CG1_CD1_length, CB_CG1_CD1_angle, CA_CB_CG1_CD1_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG1"] = &this->CG1;
		this->m_geo_atoms["CG2"] = &this->CG2;
		this->m_geo_atoms["CD1"] = &this->CD1;
		this->initAtomInGeo();
	}
}
void IleGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG1, this->m_id_counter);
	outputAtom(atomsData, CG2, this->m_id_counter);
	outputAtom(atomsData, CD1, this->m_id_counter);
}

LeuGeo::LeuGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.8652;

	C_O_length = 1.23;
	CA_C_O_angle = 120.4647;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.4948;

	CB_CG_length = 1.53;
	CA_CB_CG_angle = 116.10;
	N_CA_CB_CG_diangle = -60.1;

	CG_CD1_length = 1.524;
	CB_CG_CD1_angle = 110.27;
	CA_CB_CG_CD1_diangle = 174.9;

	CG_CD2_length = 1.525;
	CD1_CG_CD2_angle = 109.50;
	CB_CD1_CG_CD2_diangle = -120;

	m_resname =  'L';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.6;
	m_num_dasfs = 2;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void LeuGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD1_diangle = rotamers[1];
	}catch(...){
		cout << "LeuGeo rotamers list: not long enough" << endl;
	}
}
void LeuGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->CD1.m_position = calCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle);	
	this->CD2.m_position = calCoordinates(CB, CD1, CG, CG_CD2_length, CD1_CG_CD2_angle, CB_CD1_CG_CD2_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["CD1"] = &this->CD1;
		this->m_geo_atoms["CD2"] = &this->CD2;
		this->initAtomInGeo();
	}
}
void LeuGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, CD1, this->m_id_counter);
	outputAtom(atomsData, CD2, this->m_id_counter);
}

ThrGeo::ThrGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.7014;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5359;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 123.0953;

	CB_OG1_length = 1.43;
	CA_CB_OG1_angle = 109.18;
	N_CA_CB_OG1_diangle = 60.0;

	CB_CG2_length = 1.53;
	OG1_CB_CG2_angle = 109.50;
	CA_OG1_CB_CG2_diangle = 120;

	m_resname =  'T';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.6;
	m_num_dasfs = 1;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void ThrGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_OG1_diangle = rotamers[0];
	}catch(...){
		cout << "ThrGeo rotamers list: not long enough" << endl;
	}
}
void ThrGeo::addSideChain(){
	this->OG1.m_position = calCoordinates(N, CA, CB, CB_OG1_length, CA_CB_OG1_angle, N_CA_CB_OG1_diangle);
	this->CG2.m_position = calCoordinates(CA, OG1, CB, CB_CG2_length, OG1_CB_CG2_angle, CA_OG1_CB_CG2_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["OG1"] = &this->OG1;
		this->m_geo_atoms["CG2"] = &this->CG2;
		this->initAtomInGeo();
	}
}
void ThrGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, OG1, this->m_id_counter);
	outputAtom(atomsData, CG2, this->m_id_counter);
}

ArgGeo::ArgGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.98;

	C_O_length = 1.23;
	CA_C_O_angle = 120.54;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.76;

	CB_CG_length = 1.52;
	CA_CB_CG_angle = 113.83;
	N_CA_CB_CG_diangle = -65.2;

	CG_CD_length = 1.52;
	CB_CG_CD_angle = 111.79;
	CA_CB_CG_CD_diangle = -179.2;

	CD_NE_length = 1.46;
	CG_CD_NE_angle = 111.68;
	CB_CG_CD_NE_diangle = -179.3;

	NE_CZ_length = 1.33;
	CD_NE_CZ_angle = 124.79;
	CG_CD_NE_CZ_diangle = -178.7;

	CZ_NH1_length = 1.33;
	NE_CZ_NH1_angle = 120.64;
	CD_NE_CZ_NH1_diangle = 0.0;

	CZ_NH2_length = 1.33;
	NE_CZ_NH2_angle = 119.63;
	CD_NE_CZ_NH2_diangle = 180.0;

	m_resname =  'R';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 6.2;
	m_num_dasfs = 4;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void ArgGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD_diangle = rotamers[1];
		CB_CG_CD_NE_diangle = rotamers[2];
		CG_CD_NE_CZ_diangle = rotamers[3];
	}catch(...){
		cout << "ArgGeo rotamers list: not long enough" << endl;
	}
}
void ArgGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->CD.m_position = calCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
	this->NE.m_position = calCoordinates(CB, CG, CD, CD_NE_length, CG_CD_NE_angle, CB_CG_CD_NE_diangle);
	this->CZ.m_position = calCoordinates(CG, CD, NE, NE_CZ_length, CD_NE_CZ_angle, CG_CD_NE_CZ_diangle);
	this->NH1.m_position = calCoordinates(CD, NE, CZ, CZ_NH1_length, NE_CZ_NH1_angle, CD_NE_CZ_NH1_diangle);
	this->NH2.m_position = calCoordinates(CD, NE, CZ, CZ_NH2_length, NE_CZ_NH2_angle, CD_NE_CZ_NH2_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["CD"] = &this->CD;
		this->m_geo_atoms["NE"] = &this->NE;
		this->m_geo_atoms["CZ"] = &this->CZ;
		this->m_geo_atoms["NH1"] = &this->NH1;
		this->m_geo_atoms["NH2"] = &this->NH2;
		this->initAtomInGeo();
	}
}
void ArgGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, CD, this->m_id_counter);
	outputAtom(atomsData, NE, this->m_id_counter);
	outputAtom(atomsData, CZ, this->m_id_counter);
	outputAtom(atomsData, NH1, this->m_id_counter);
	outputAtom(atomsData, NH2, this->m_id_counter);
}

LysGeo::LysGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.08;

	C_O_length = 1.23;
	CA_C_O_angle = 120.54;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.76;

	CB_CG_length = 1.52;
	CA_CB_CG_angle = 113.83;
	N_CA_CB_CG_diangle = -64.5;

	CG_CD_length = 1.52;
	CB_CG_CD_angle = 111.79;
	CA_CB_CG_CD_diangle = -178.1;

	CD_CE_length = 1.46;
	CG_CD_CE_angle = 111.68;
	CB_CG_CD_CE_diangle = -179.6;

	CE_NZ_length = 1.33;
	CD_CE_NZ_angle = 124.79;
	CG_CD_CE_NZ_diangle = 179.6;

	m_resname =  'K';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 5;
	m_num_dasfs = 4;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void LysGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD_diangle = rotamers[1];
		CB_CG_CD_CE_diangle = rotamers[2];
		CG_CD_CE_NZ_diangle = rotamers[3];
	}catch(...){
		cout << "LysGeo rotamers list: not long enough" << endl;
	}
}
void LysGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->CD.m_position = calCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
	this->CE.m_position = calCoordinates(CB, CG, CD, CD_CE_length, CG_CD_CE_angle, CB_CG_CD_CE_diangle);
	this->NZ.m_position = calCoordinates(CG, CD, CE, CE_NZ_length, CD_CE_NZ_angle, CG_CD_CE_NZ_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["CD"] = &this->CD;
		this->m_geo_atoms["CE"] = &this->CE;
		this->m_geo_atoms["NZ"] = &this->NZ;
		this->initAtomInGeo();
	}
}
void LysGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, CD, this->m_id_counter);
	outputAtom(atomsData, CE, this->m_id_counter);
	outputAtom(atomsData, NZ, this->m_id_counter);
}

AspGeo::AspGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.03;    

	C_O_length = 1.23;
	CA_C_O_angle = 120.51;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.82;

	CB_CG_length = 1.52;
	CA_CB_CG_angle = 113.06;
	N_CA_CB_CG_diangle = -66.4;

	CG_OD1_length = 1.25;
	CB_CG_OD1_angle = 119.22;
	CA_CB_CG_OD1_diangle = -46.7;

	CG_OD2_length = 1.25;
	CB_CG_OD2_angle = 118.218;
	CA_CB_CG_OD2_diangle = 180+CA_CB_CG_OD1_diangle;

	m_resname =  'D';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.6;
	m_num_dasfs = 2;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void AspGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_OD1_diangle = rotamers[1];
		if (CA_CB_CG_OD1_diangle > 0){
			CA_CB_CG_OD2_diangle = CA_CB_CG_OD1_diangle-180.0;
		}else{
			CA_CB_CG_OD2_diangle = CA_CB_CG_OD1_diangle+180.0;
		}
	}catch(...){
		cout << "AspGeo rotamers list: not long enough" << endl;
	}
}
void AspGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->OD1.m_position = calCoordinates(CA, CB, CG, CG_OD1_length, CB_CG_OD1_angle, CA_CB_CG_OD1_diangle);
	this->OD2.m_position = calCoordinates(CA, CB, CG, CG_OD2_length, CB_CG_OD2_angle, CA_CB_CG_OD2_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["OD1"] = &this->OD1;
		this->m_geo_atoms["OD2"] = &this->OD2;
		this->initAtomInGeo();
	}
}
void AspGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, OD1, this->m_id_counter);
	outputAtom(atomsData, OD2, this->m_id_counter);
}

AsnGeo::AsnGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.5;

	C_O_length = 1.23;
	CA_C_O_angle = 120.4826;
	N_CA_C_O_diangle =  -60.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 123.2254;

	CB_CG_length = 1.52;
	CA_CB_CG_angle = 112.62;
	N_CA_CB_CG_diangle = -65.5;

	CG_OD1_length = 1.23;
	CB_CG_OD1_angle = 120.85;
	CA_CB_CG_OD1_diangle = -58.3;

	CG_ND2_length = 1.33;
	CB_CG_ND2_angle = 116.48;
	CA_CB_CG_ND2_diangle = 180.0+CA_CB_CG_OD1_diangle;

	m_resname =  'N';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.6;
	m_num_dasfs = 2;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void AsnGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_OD1_diangle = rotamers[1];
		if (CA_CB_CG_OD1_diangle > 0){
			CA_CB_CG_ND2_diangle = CA_CB_CG_OD1_diangle-180.0;
		}else{
			CA_CB_CG_ND2_diangle = CA_CB_CG_OD1_diangle+180.0;
		}
	}catch(...){
		cout << "AsnGeo rotamers list: not long enough" << endl;
	}
}
void AsnGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->OD1.m_position = calCoordinates(CA, CB, CG, CG_OD1_length, CB_CG_OD1_angle, CA_CB_CG_OD1_diangle);
	this->ND2.m_position = calCoordinates(CA, CB, CG, CG_ND2_length, CB_CG_ND2_angle, CA_CB_CG_ND2_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["OD1"] = &this->OD1;
		this->m_geo_atoms["ND2"] = &this->ND2;
		this->initAtomInGeo();
	}
}
void AsnGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, OD1, this->m_id_counter);
	outputAtom(atomsData, ND2, this->m_id_counter);
}

GluGeo::GluGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.1703;

	C_O_length = 1.23;
	CA_C_O_angle = 120.511;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.8702;

	CB_CG_length = 1.52;
	CA_CB_CG_angle = 113.82;
	N_CA_CB_CG_diangle = -63.8;

	CG_CD_length = 1.52;
	CB_CG_CD_angle = 113.31;
	CA_CB_CG_CD_diangle = -179.8;

	CD_OE1_length = 1.25;
	CG_CD_OE1_angle = 119.02;
	CB_CG_CD_OE1_diangle = -6.2;

	CD_OE2_length = 1.25;
	CG_CD_OE2_angle = 118.08;
	CB_CG_CD_OE2_diangle = 180.0+CB_CG_CD_OE1_diangle;

	m_resname =  'E';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.7;
	m_num_dasfs = 3;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void GluGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD_diangle = rotamers[1];
		CB_CG_CD_OE1_diangle = rotamers[2];
		if (CB_CG_CD_OE1_diangle > 0){
			CB_CG_CD_OE2_diangle = CB_CG_CD_OE1_diangle-180.0;
		}else{
			CB_CG_CD_OE2_diangle = CB_CG_CD_OE1_diangle+180.0;
		}
	}catch(...){
		cout << "GluGeo rotamers list: not long enough" << endl;
	}
}
void GluGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->CD.m_position = calCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
	this->OE1.m_position = calCoordinates(CB, CG, CD, CD_OE1_length, CG_CD_OE1_angle, CB_CG_CD_OE1_diangle);
	this->OE2.m_position = calCoordinates(CB, CG, CD, CD_OE2_length, CG_CD_OE2_angle, CB_CG_CD_OE2_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["CD"] = &this->CD;
		this->m_geo_atoms["OE1"] = &this->OE1;
		this->m_geo_atoms["OE2"] = &this->OE2;
		this->initAtomInGeo();
	}
}
void GluGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, CD, this->m_id_counter);
	outputAtom(atomsData, OE1, this->m_id_counter);
	outputAtom(atomsData, OE2, this->m_id_counter);
}

GlnGeo::GlnGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.0849;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5029;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.8134;

	CB_CG_length = 1.52;
	CA_CB_CG_angle = 113.75;
	N_CA_CB_CG_diangle = -60.2;

	CG_CD_length = 1.52;
	CB_CG_CD_angle = 112.78;
	CA_CB_CG_CD_diangle = -69.6;

	CD_OE1_length = 1.24;
	CG_CD_OE1_angle = 120.86;
	CB_CG_CD_OE1_diangle = -50.5;

	CD_NE2_length = 1.33;
	CG_CD_NE2_angle = 116.50;
	CB_CG_CD_NE2_diangle = 180+CB_CG_CD_OE1_diangle;

	m_resname =  'Q';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.8;
	m_num_dasfs = 3;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void GlnGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD_diangle = rotamers[1];
		CB_CG_CD_OE1_diangle = rotamers[2];
		if (CB_CG_CD_OE1_diangle > 0){
			CB_CG_CD_NE2_diangle = CB_CG_CD_OE1_diangle-180.0;
		}else{
			CB_CG_CD_NE2_diangle = CB_CG_CD_OE1_diangle+180.0;
		}
	}catch(...){
		cout << "GlnGeo rotamers list: not long enough" << endl;
	}
}
void GlnGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->CD.m_position = calCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
	this->OE1.m_position = calCoordinates(CB, CG, CD, CD_OE1_length, CG_CD_OE1_angle, CB_CG_CD_OE1_diangle);
	this->NE2.m_position = calCoordinates(CB, CG, CD, CD_NE2_length, CG_CD_NE2_angle, CB_CG_CD_NE2_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["CD"] = &this->CD;
		this->m_geo_atoms["OE1"] = &this->OE1;
		this->m_geo_atoms["NE2"] = &this->NE2;
		this->initAtomInGeo();
	}
}
void GlnGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, CD, this->m_id_counter);
	outputAtom(atomsData, OE1, this->m_id_counter);
	outputAtom(atomsData, NE2, this->m_id_counter);
}

MetGeo::MetGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.9416;

	C_O_length = 1.23;
	CA_C_O_angle = 120.4816;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6733;

	CB_CG_length = 1.52;
	CA_CB_CG_angle =  113.68;
	N_CA_CB_CG_diangle = -64.4;

	CG_SD_length = 1.81;
	CB_CG_SD_angle = 112.69;
	CA_CB_CG_SD_diangle = -179.6;

	SD_CE_length = 1.79;
	CG_SD_CE_angle = 100.61;
	CB_CG_SD_CE_diangle = 70.1;

	m_resname = 'M';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 4.2;
	m_num_dasfs = 3;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void MetGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_SD_diangle = rotamers[1];
		CB_CG_SD_CE_diangle = rotamers[2];
	}catch(...){
		cout << "MetGeo rotamers list: not long enough" << endl;
	}
}
void MetGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->SD.m_position = calCoordinates(CA, CB, CG, CG_SD_length, CB_CG_SD_angle, CA_CB_CG_SD_diangle);
	this->CE.m_position = calCoordinates(CB, CG, SD, SD_CE_length, CG_SD_CE_angle, CB_CG_SD_CE_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["SD"] = &this->SD;
		this->m_geo_atoms["CE"] = &this->CE;
		this->initAtomInGeo();
	}
}
void MetGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, SD, this->m_id_counter);
	outputAtom(atomsData, CE, this->m_id_counter);
}

HisGeo::HisGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 111.0859;

	C_O_length = 1.23;
	CA_C_O_angle = 120.4732;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6711;

	CB_CG_length = 1.49;
	CA_CB_CG_angle = 113.74;
	N_CA_CB_CG_diangle = -63.2;

	CG_ND1_length = 1.38;
	CB_CG_ND1_angle = 122.85;
	CA_CB_CG_ND1_diangle = -75.7;    

	CG_CD2_length = 1.35;
	CB_CG_CD2_angle = 130.61;
	CA_CB_CG_CD2_diangle = 180.0+CA_CB_CG_ND1_diangle;

	ND1_CE1_length = 1.32;
	CG_ND1_CE1_angle = 108.5;
	CB_CG_ND1_CE1_diangle = 180.0;

	CD2_NE2_length = 1.35;
	CG_CD2_NE2_angle = 108.5;
	CB_CG_CD2_NE2_diangle = 180.0;

	m_resname =  'H';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.7;
	m_num_dasfs = 2;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void HisGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_ND1_diangle = rotamers[1];
		if (CA_CB_CG_ND1_diangle > 0){
			CA_CB_CG_CD2_diangle = CA_CB_CG_ND1_diangle-180.0;
		}else{
			CA_CB_CG_CD2_diangle = CA_CB_CG_ND1_diangle+180.0;
		}
	}catch(...){
		cout << "HisGeo rotamers list: not long enough" << endl;
	}
}
void HisGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->ND1.m_position = calCoordinates(CA, CB, CG, CG_ND1_length, CB_CG_ND1_angle, CA_CB_CG_ND1_diangle);
	this->CD2.m_position = calCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
	this->CE1.m_position = calCoordinates(CB, CG, ND1, ND1_CE1_length, CG_ND1_CE1_angle, CB_CG_ND1_CE1_diangle);
	this->NE2.m_position = calCoordinates(CB, CG, CD2, CD2_NE2_length, CG_CD2_NE2_angle, CB_CG_CD2_NE2_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["ND1"] = &this->ND1;
		this->m_geo_atoms["CD2"] = &this->CD2;
		this->m_geo_atoms["CE1"] = &this->CE1;
		this->m_geo_atoms["NE2"] = &this->NE2;
		this->initAtomInGeo();
	}
}
void HisGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, ND1, this->m_id_counter);
	outputAtom(atomsData, CD2, this->m_id_counter);
	outputAtom(atomsData, CE1, this->m_id_counter);
	outputAtom(atomsData, NE2, this->m_id_counter);
}

ProGeo::ProGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 112.7499;

	C_O_length = 1.23;
	CA_C_O_angle = 120.2945;
	N_CA_C_O_diangle = -45.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle  = 116.642992978143;
	_C_N_CA_angle =  121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 115.2975;

	CB_CG_length = 1.49;
	CA_CB_CG_angle = 104.21;
	N_CA_CB_CG_diangle = 29.6;

	CG_CD_length = 1.50;
	CB_CG_CD_angle = 105.03;
	CA_CB_CG_CD_diangle = -34.8;

	m_resname =  'P';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 3.5;
	m_num_dasfs = 2;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void ProGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD_diangle = rotamers[1];
	}catch(...){
		cout << "ProGeo rotamers list: not long enough" << endl;
	}
}
void ProGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->CD.m_position = calCoordinates(CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["CD"] = &this->CD;
		this->initAtomInGeo();
	}
}
void ProGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, CD, this->m_id_counter);
}

PheGeo::PheGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.7528;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5316;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6054;

	CB_CG_length = 1.50;
	CA_CB_CG_angle = 113.85;
	N_CA_CB_CG_diangle = -64.7;

	CG_CD1_length = 1.39;
	CB_CG_CD1_angle = 120.0;
	CA_CB_CG_CD1_diangle = 93.3;

	CG_CD2_length = 1.39;
	CB_CG_CD2_angle = 120.0;
	CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle-180.0;

	CD1_CE1_length = 1.39;
	CG_CD1_CE1_angle = 120.0;
	CB_CG_CD1_CE1_diangle = 180.0;

	CD2_CE2_length = 1.39;
	CG_CD2_CE2_angle = 120.0;
	CB_CG_CD2_CE2_diangle = 180.0;

	CE1_CZ_length = 1.39;
	CD1_CE1_CZ_angle = 120.0;
	CG_CD1_CE1_CZ_diangle = 0.0;

	m_resname =  'F';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 4.3;
	m_num_dasfs = 2;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void PheGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD1_diangle = rotamers[1];
		if (CA_CB_CG_CD1_diangle > 0){
			CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle-180.0;
		}else{
			CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle+180.0;
		}
	}catch(...){
		cout << "PheGeo rotamers list: not long enough" << endl;
	}
}
void PheGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->CD1.m_position = calCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle);
	this->CD2.m_position = calCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
	this->CE1.m_position = calCoordinates(CB, CG, CD1, CD1_CE1_length, CG_CD1_CE1_angle, CB_CG_CD1_CE1_diangle);
	this->CE2.m_position = calCoordinates(CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle);
	this->CZ.m_position = calCoordinates(CG, CD1, CE1, CE1_CZ_length, CD1_CE1_CZ_angle, CG_CD1_CE1_CZ_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["CD1"] = &this->CD1;
		this->m_geo_atoms["CD2"] = &this->CD2;
		this->m_geo_atoms["CE1"] = &this->CE1;
		this->m_geo_atoms["CE2"] = &this->CE2;
		this->m_geo_atoms["CZ"] = &this->CZ;
		this->initAtomInGeo();
	}
}
void PheGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, CD1, this->m_id_counter);
	outputAtom(atomsData, CD2, this->m_id_counter);
	outputAtom(atomsData, CE1, this->m_id_counter);
	outputAtom(atomsData, CE2, this->m_id_counter);
	outputAtom(atomsData, CZ, this->m_id_counter);
}

TyrGeo::TyrGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.9288;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5434;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6023;

	CB_CG_length = 1.51;
	CA_CB_CG_angle =  113.8;
	N_CA_CB_CG_diangle = -64.3;

	CG_CD1_length = 1.39;
	CB_CG_CD1_angle = 120.98;
	CA_CB_CG_CD1_diangle = 93.1;

	CG_CD2_length = 1.39;
	CB_CG_CD2_angle = 120.82;
	CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle-180.0;

	CD1_CE1_length = 1.39;
	CG_CD1_CE1_angle = 120.0;
	CB_CG_CD1_CE1_diangle = 180.0;

	CD2_CE2_length = 1.39;
	CG_CD2_CE2_angle = 120.0;
	CB_CG_CD2_CE2_diangle = 180.0;

	CE1_CZ_length = 1.39;
	CD1_CE1_CZ_angle = 120.0;
	CG_CD1_CE1_CZ_diangle = 0.0;

	CZ_OH_length = 1.39;
	CE1_CZ_OH_angle = 119.78;
	CD1_CE1_CZ_OH_diangle = 180.0;

	m_resname =  'Y';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 5.7;
	m_num_dasfs = 2;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void TyrGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD1_diangle = rotamers[1];
		if (CA_CB_CG_CD1_diangle > 0){
			CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle-180.0;
		}else{
			CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle+180.0;
		}
	}catch(...){
		cout << "TyrGeo rotamers list: not long enough" << endl;
	}

}
void TyrGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->CD1.m_position = calCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle);
	this->CD2.m_position = calCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
	this->CE1.m_position = calCoordinates(CB, CG, CD1, CD1_CE1_length, CG_CD1_CE1_angle, CB_CG_CD1_CE1_diangle);
	this->CE2.m_position = calCoordinates(CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle);
	this->CZ.m_position = calCoordinates(CG, CD1, CE1, CE1_CZ_length, CD1_CE1_CZ_angle, CG_CD1_CE1_CZ_diangle);
	this->OH.m_position = calCoordinates(CD1, CE1, CZ, CZ_OH_length, CE1_CZ_OH_angle, CD1_CE1_CZ_OH_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["CD1"] = &this->CD1;
		this->m_geo_atoms["CD2"] = &this->CD2;
		this->m_geo_atoms["CE1"] = &this->CE1;
		this->m_geo_atoms["CE2"] = &this->CE2;
		this->m_geo_atoms["CZ"] = &this->CZ;
		this->m_geo_atoms["OH"] = &this->OH;
		this->initAtomInGeo();
	}
}
void TyrGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, CD1, this->m_id_counter);
	outputAtom(atomsData, CD2, this->m_id_counter);
	outputAtom(atomsData, CE1, this->m_id_counter);
	outputAtom(atomsData, CE2, this->m_id_counter);
	outputAtom(atomsData, CZ, this->m_id_counter);
	outputAtom(atomsData, OH, this->m_id_counter);
}

TrpGeo::TrpGeo(int param_resid, int param_atomid){
	CA_N_length = 1.46;
	CA_C_length = 1.52;
	N_CA_C_angle = 110.8914;

	C_O_length = 1.23;
	CA_C_O_angle = 120.5117;
	N_CA_C_O_diangle = 120.0;

	phi = -120;
	psi = 140;
	omega = 180.0;
	_peptide_bond = 1.33;
	_CA_C_N_angle = 116.642992978143;
	_C_N_CA_angle = 121.382215820277;

	CA_CB_length = 1.52;
	C_CA_CB_angle = 109.5;
	N_C_CA_CB_diangle = 122.6112;

	CB_CG_length = 1.50;
	CA_CB_CG_angle = 114.10;
	N_CA_CB_CG_diangle = -66.4;

	CG_CD1_length = 1.37;
	CB_CG_CD1_angle = 127.07;
	CA_CB_CG_CD1_diangle = 96.3;

	CG_CD2_length = 1.43;
	CB_CG_CD2_angle = 126.66;
	CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle-180.0;

	CD1_NE1_length = 1.38;
	CG_CD1_NE1_angle = 108.5;
	CB_CG_CD1_NE1_diangle = 180.0;

	CD2_CE2_length = 1.40;
	CG_CD2_CE2_angle = 108.5;
	CB_CG_CD2_CE2_diangle = 180.0;

	CD2_CE3_length = 1.40;
	CG_CD2_CE3_angle = 133.83;
	CB_CG_CD2_CE3_diangle = 0.0;

	CE2_CZ2_length = 1.40;
	CD2_CE2_CZ2_angle = 120.0;
	CG_CD2_CE2_CZ2_diangle = 180.0;

	CE3_CZ3_length = 1.40;
	CD2_CE3_CZ3_angle = 120.0;
	CG_CD2_CE3_CZ3_diangle = 180.0;

	CZ2_CH2_length = 1.40;
	CE2_CZ2_CH2_angle = 120.0;
	CD2_CE2_CZ2_CH2_diangle = 0.0;

	m_resname =  'W';
	m_resid = param_resid;
	m_atomid = param_atomid;
	m_maxdis = 5.4;
	m_num_dasfs = 2;

	this->addMainChainAtomsToMap();
	this->initAtomInGeo();
}
void TrpGeo::setRotamers(const vector<double>& rotamers){
	try{
		N_CA_CB_CG_diangle = rotamers[0];
		CA_CB_CG_CD1_diangle = rotamers[1];
		if (CA_CB_CG_CD1_diangle > 0){
			CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle-180.0;
		}else{
			CA_CB_CG_CD2_diangle = CA_CB_CG_CD1_diangle+180.0;
		}
	}catch(...){
		cout << "TrpGeo rotamers list: not long enough" << endl;
	}
}
void TrpGeo::addSideChain(){
	this->CG.m_position = calCoordinates(N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle);
	this->CD1.m_position = calCoordinates(CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle);
	this->CD2.m_position = calCoordinates(CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle);
	this->NE1.m_position = calCoordinates(CB, CG, CD1, CD1_NE1_length, CG_CD1_NE1_angle, CB_CG_CD1_NE1_diangle);
	this->CE2.m_position = calCoordinates(CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle);
	this->CE3.m_position = calCoordinates(CB, CG, CD2, CD2_CE3_length, CG_CD2_CE3_angle, CB_CG_CD2_CE3_diangle);
	this->CZ2.m_position = calCoordinates(CG, CD2, CE2, CE2_CZ2_length, CD2_CE2_CZ2_angle, CG_CD2_CE2_CZ2_diangle);
	this->CZ3.m_position = calCoordinates(CG, CD2, CE3, CE3_CZ3_length, CD2_CE3_CZ3_angle, CG_CD2_CE3_CZ3_diangle);
	this->CH2.m_position = calCoordinates(CD2, CE2, CZ2, CZ2_CH2_length, CE2_CZ2_CH2_angle, CD2_CE2_CZ2_CH2_diangle);
	if (this->m_geo_atoms.size() == 5){
		this->m_geo_atoms["CG"] = &this->CG;
		this->m_geo_atoms["CD1"] = &this->CD1;
		this->m_geo_atoms["CD2"] = &this->CD2;
		this->m_geo_atoms["NE1"] = &this->NE1;
		this->m_geo_atoms["CE2"] = &this->CE2;
		this->m_geo_atoms["CE3"] = &this->CE3;
		this->m_geo_atoms["CZ2"] = &this->CZ2;
		this->m_geo_atoms["CZ3"] = &this->CZ3;
		this->m_geo_atoms["CH2"] = &this->CH2;
		this->initAtomInGeo();
	}
}
void TrpGeo::output_sidechain(vector<Atom>& atomsData){
	outputAtom(atomsData, CG, this->m_id_counter);
	outputAtom(atomsData, CD1, this->m_id_counter);
	outputAtom(atomsData, CD2, this->m_id_counter);
	outputAtom(atomsData, NE1, this->m_id_counter);
	outputAtom(atomsData, CE2, this->m_id_counter);
	outputAtom(atomsData, CE3, this->m_id_counter);
	outputAtom(atomsData, CZ2, this->m_id_counter);
	outputAtom(atomsData, CZ3, this->m_id_counter);
	outputAtom(atomsData, CH2, this->m_id_counter);
}

