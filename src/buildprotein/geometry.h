#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream> 
#include <vector> 
#include <string> 
#include <map> 
#include "atom.h"
#include "basic.h"

using namespace std;

class Geo{
public:
	double CA_N_length;
	double CA_C_length;
	double N_CA_C_angle;

	double C_O_length;
	double CA_C_O_angle;
	double N_CA_C_O_diangle;

	double phi;
	double psi;
	double omega;
	double ta_possibility;

	vector<double> dihedral;
	double dihedral_score;
	vector<vector<DASFFeature>> standerDASFs;

	double _peptide_bond;
	double _CA_C_N_angle;
	double _C_N_CA_angle;

	char m_resname;
	int m_resid;
	int m_atomid;

	//for side-chain modeling
	string phipsi_rotamerformat;
	double m_maxdis;
	int m_num_dasfs;

	Atom N = Atom("N");
	Atom CA = Atom("CA");
	Atom C = Atom("C");
	Atom O = Atom("O");

	map<string, Atom*> m_geo_atoms;

	virtual void addSideChain();
	virtual void setRotamers(const vector<double>& rotamers);
	
	int m_id_counter = 0;
	void initAtomInGeo();

	virtual void addMainChainAtomsToMap();
	virtual void output(vector<Atom>& atomsData);

};

class GeoContainCB: public Geo{
public:
	double CA_CB_length;
	double C_CA_CB_angle;
	double N_C_CA_CB_diangle;
	
	Atom CB = Atom("CB");

	void addMainChainAtomsToMap();
	void output(vector<Atom>& atomsData);
};

class GeoContainSideChain: public GeoContainCB{
public:
	virtual void output_sidechain(vector<Atom>& atomsData);
};

class GlyGeo: public Geo{
public:

	GlyGeo(int param_resid, int param_atomid);

};

class AlaGeo: public GeoContainCB{
public:

	AlaGeo(int param_resid, int param_atomid);

};

class SerGeo: public GeoContainSideChain{
public:

	double CB_OG_length;
	double CA_CB_OG_angle;
	double N_CA_CB_OG_diangle;

	Atom OG = Atom("OG");

	SerGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);

};

class CysGeo: public GeoContainSideChain{
public:

	double CB_SG_length;
	double CA_CB_SG_angle;
	double N_CA_CB_SG_diangle;

	Atom SG = Atom("SG");

	CysGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);

};                                        

class ValGeo: public GeoContainSideChain{				
public:

	double CB_CG1_length;
	double CA_CB_CG1_angle;
	double N_CA_CB_CG1_diangle;

	double CB_CG2_length;
	double CG1_CB_CG2_angle; 
	double CA_CG1_CB_CG2_diangle;

	Atom CG1 = Atom("CG1");
	Atom CG2 = Atom("CG2");

	ValGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};                        

class IleGeo: public GeoContainSideChain{
public:

	double CB_CG1_length;
	double CA_CB_CG1_angle;
	double N_CA_CB_CG1_diangle;

	double CB_CG2_length;
	double CG1_CB_CG2_angle;
	double CA_CG1_CB_CG2_diangle;

	double CG1_CD1_length;
	double CB_CG1_CD1_angle;
	double CA_CB_CG1_CD1_diangle;

	Atom CG1 = Atom("CG1");
	Atom CG2 = Atom("CG2");
	Atom CD1 = Atom("CD1");

	IleGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
}; 

class LeuGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_CD1_length;
	double CB_CG_CD1_angle;
	double CA_CB_CG_CD1_diangle;

	double CG_CD2_length;
	double CD1_CG_CD2_angle;
	double CB_CD1_CG_CD2_diangle;

	Atom CG = Atom("CG");
	Atom CD1 = Atom("CD1");
	Atom CD2 = Atom("CD2");

	LeuGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};                              

class ThrGeo: public GeoContainSideChain{
public:

	double CB_OG1_length;
	double CA_CB_OG1_angle;
	double N_CA_CB_OG1_diangle;

	double CB_CG2_length;
	double OG1_CB_CG2_angle;
	double CA_OG1_CB_CG2_diangle;

	Atom OG1 = Atom("OG1");
	Atom CG2 = Atom("CG2");

	ThrGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
	
};

class ArgGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_CD_length;
	double CB_CG_CD_angle;
	double CA_CB_CG_CD_diangle;

	double CD_NE_length;
	double CG_CD_NE_angle;
	double CB_CG_CD_NE_diangle;

	double NE_CZ_length;
	double CD_NE_CZ_angle;
	double CG_CD_NE_CZ_diangle;

	double CZ_NH1_length;
	double NE_CZ_NH1_angle;
	double CD_NE_CZ_NH1_diangle;

	double CZ_NH2_length;
	double NE_CZ_NH2_angle;
	double CD_NE_CZ_NH2_diangle;

	Atom CG = Atom("CG");
	Atom CD = Atom("CD");
	Atom NE = Atom("NE");
	Atom CZ = Atom("CZ");
	Atom NH1 = Atom("NH1");
	Atom NH2 = Atom("NH2");

	ArgGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};

class LysGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_CD_length;
	double CB_CG_CD_angle;
	double CA_CB_CG_CD_diangle;

	double CD_CE_length;
	double CG_CD_CE_angle;
	double CB_CG_CD_CE_diangle;

	double CE_NZ_length;
	double CD_CE_NZ_angle;
	double CG_CD_CE_NZ_diangle;

	Atom CG = Atom("CG");
	Atom CD = Atom("CD");
	Atom CE = Atom("CE");
	Atom NZ = Atom("NZ");

	LysGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};

class AspGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_OD1_length;
	double CB_CG_OD1_angle;
	double CA_CB_CG_OD1_diangle;

	double CG_OD2_length;
	double CB_CG_OD2_angle;
	double CA_CB_CG_OD2_diangle;

	Atom CG = Atom("CG");
	Atom OD1 = Atom("OD1");
	Atom OD2 = Atom("OD2");

	AspGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};

class AsnGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_OD1_length;
	double CB_CG_OD1_angle;
	double CA_CB_CG_OD1_diangle;

	double CG_ND2_length;
	double CB_CG_ND2_angle;
	double CA_CB_CG_ND2_diangle;

	Atom CG = Atom("CG");
	Atom OD1 = Atom("OD1");
	Atom ND2 = Atom("ND2");

	AsnGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};

class GluGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_CD_length;
	double CB_CG_CD_angle;
	double CA_CB_CG_CD_diangle;

	double CD_OE1_length;
	double CG_CD_OE1_angle;
	double CB_CG_CD_OE1_diangle;

	double CD_OE2_length;
	double CG_CD_OE2_angle;
	double CB_CG_CD_OE2_diangle;

	Atom CG = Atom("CG");
	Atom CD = Atom("CD");
	Atom OE1 = Atom("OE1");
	Atom OE2 = Atom("OE2");

	GluGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};

class GlnGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_CD_length;
	double CB_CG_CD_angle;
	double CA_CB_CG_CD_diangle;

	double CD_OE1_length;
	double CG_CD_OE1_angle;
	double CB_CG_CD_OE1_diangle;

	double CD_NE2_length;
	double CG_CD_NE2_angle;
	double CB_CG_CD_NE2_diangle;

	Atom CG = Atom("CG");
	Atom CD = Atom("CD");
	Atom OE1 = Atom("OE1");
	Atom NE2 = Atom("NE2");

	GlnGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};              

class MetGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_SD_length;
	double CB_CG_SD_angle;
	double CA_CB_CG_SD_diangle;

	double SD_CE_length;
	double CG_SD_CE_angle;
	double CB_CG_SD_CE_diangle;

	Atom CG = Atom("CG");
	Atom SD = Atom("SD");
	Atom CE = Atom("CE");

	MetGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};           

class HisGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_ND1_length;
	double CB_CG_ND1_angle;
	double CA_CB_CG_ND1_diangle;    

	double CG_CD2_length;
	double CB_CG_CD2_angle;
	double CA_CB_CG_CD2_diangle;

	double ND1_CE1_length;
	double CG_ND1_CE1_angle;
	double CB_CG_ND1_CE1_diangle;

	double CD2_NE2_length;
	double CG_CD2_NE2_angle;
	double CB_CG_CD2_NE2_diangle;

	Atom CG = Atom("CG");
	Atom ND1 = Atom("ND1");
	Atom CD2 = Atom("CD2");
	Atom CE1 = Atom("CE1");
	Atom NE2 = Atom("NE2");

	HisGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};              

class ProGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_CD_length;
	double CB_CG_CD_angle;
	double CA_CB_CG_CD_diangle;

	Atom CG = Atom("CG");
	Atom CD = Atom("CD");

	ProGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};

class PheGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_CD1_length;
	double CB_CG_CD1_angle;
	double CA_CB_CG_CD1_diangle;

	double CG_CD2_length;
	double CB_CG_CD2_angle;
	double CA_CB_CG_CD2_diangle;

	double CD1_CE1_length;
	double CG_CD1_CE1_angle;
	double CB_CG_CD1_CE1_diangle;

	double CD2_CE2_length;
	double CG_CD2_CE2_angle;
	double CB_CG_CD2_CE2_diangle;

	double CE1_CZ_length;
	double CD1_CE1_CZ_angle;
	double CG_CD1_CE1_CZ_diangle;

	Atom CG = Atom("CG");
	Atom CD1 = Atom("CD1");
	Atom CD2 = Atom("CD2");
	Atom CE1 = Atom("CE1");
	Atom CE2 = Atom("CE2");
	Atom CZ = Atom("CZ");

	PheGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};

class TyrGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_CD1_length;
	double CB_CG_CD1_angle;
	double CA_CB_CG_CD1_diangle;

	double CG_CD2_length;
	double CB_CG_CD2_angle;
	double CA_CB_CG_CD2_diangle;

	double CD1_CE1_length;
	double CG_CD1_CE1_angle;
	double CB_CG_CD1_CE1_diangle;

	double CD2_CE2_length;
	double CG_CD2_CE2_angle;
	double CB_CG_CD2_CE2_diangle;

	double CE1_CZ_length;
	double CD1_CE1_CZ_angle;
	double CG_CD1_CE1_CZ_diangle;

	double CZ_OH_length;
	double CE1_CZ_OH_angle;
	double CD1_CE1_CZ_OH_diangle;

	Atom CG = Atom("CG");
	Atom CD1 = Atom("CD1");
	Atom CD2 = Atom("CD2");
	Atom CE1 = Atom("CE1");
	Atom CE2 = Atom("CE2");
	Atom CZ = Atom("CZ");
	Atom OH = Atom("OH");

	TyrGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
};

class TrpGeo: public GeoContainSideChain{
public:

	double CB_CG_length;
	double CA_CB_CG_angle;
	double N_CA_CB_CG_diangle;

	double CG_CD1_length;
	double CB_CG_CD1_angle;
	double CA_CB_CG_CD1_diangle;

	double CG_CD2_length;
	double CB_CG_CD2_angle;
	double CA_CB_CG_CD2_diangle;

	double CD1_NE1_length;
	double CG_CD1_NE1_angle;
	double CB_CG_CD1_NE1_diangle;

	double CD2_CE2_length;
	double CG_CD2_CE2_angle;
	double CB_CG_CD2_CE2_diangle;

	double CD2_CE3_length;
	double CG_CD2_CE3_angle;
	double CB_CG_CD2_CE3_diangle;

	double CE2_CZ2_length;
	double CD2_CE2_CZ2_angle;
	double CG_CD2_CE2_CZ2_diangle;

	double CE3_CZ3_length;
	double CD2_CE3_CZ3_angle;
	double CG_CD2_CE3_CZ3_diangle;

	double CZ2_CH2_length;
	double CE2_CZ2_CH2_angle;
	double CD2_CE2_CZ2_CH2_diangle;

	Atom CG = Atom("CG");
	Atom CD1 = Atom("CD1");
	Atom CD2 = Atom("CD2");
	Atom NE1 = Atom("NE1");
	Atom CE2 = Atom("CE2");
	Atom CE3 = Atom("CE3");
	Atom CZ2 = Atom("CZ2");
	Atom CZ3 = Atom("CZ3");
	Atom CH2 = Atom("CH2");

	TrpGeo(int param_resid, int param_atomid);
	
	void setRotamers(const vector<double>& rotamers);
	void addSideChain();
	void output_sidechain(vector<Atom>& atomsData);
}; 

#endif