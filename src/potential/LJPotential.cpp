#include "LJPotential.h"

double calLJPotential(double d_star2, double eij, double lambd){
	double lj_potential = 0;
	if (d_star2 >= 0 && d_star2 <= 0.565){
		lj_potential = lambd*(49.69 - 40.06*sqrt(d_star2));
	}else if (d_star2 > 0.565 && d_star2 <= 0.797){
		lj_potential = lambd*eij*(pow(d_star2, -6) - 2 * pow(d_star2, -3));
	}else if (d_star2 > 0.797 && d_star2 < 6.25){
		lj_potential = eij*(pow(d_star2, -6) - 2 * pow(d_star2, -3));
	}else{
		lj_potential = 0;
	}

	return lj_potential;
}

double getAij2(const Atom& a, const Atom& b){
	int type1, type2;
	if(a.m_LJType > b.m_LJType){
		type1 = b.m_LJType;
		type2 = a.m_LJType;
	}else{
		type1 = a.m_LJType;
		type2 = b.m_LJType;	
	}
	double aij2 = -1;
	//Special cases for some atoms
	if(type1==9 || type1==10 || type1==11){
		if(type2==11 || type2==12){
			aij2=9.61;
		}else if(type2==13 || type2==14 || type2==15 || type2==16){
			aij2=8.41;
		}
	}else if(type1==12){
		if(type2==15 || type2==16){
			aij2=8.41;	
		}
	}else if(type1==13){
		if(type2==15 || type2==16){
			aij2=7.84;	
		}
	}else if(type1==15 || type1==16){
		if(type2==16){
			aij2=7.84;	
		}
	}else if(type1==17 && type2==17){
		aij2=4;
	}else if(type1==8 && type2==17){
		aij2=9;
	}

	if(aij2 == -1){
		aij2 = pow((a.m_radii + b.m_radii), 2);
	}

	return aij2;
}

double getMSLJPotential(const vector<Residue>& residuesData){

	double potentials = 0;

	int length = residuesData.size();
	for (int i = 0; i < length; ++i){
		//a atom side chain
		for (auto atom_a_iter = residuesData[i].m_atoms.begin(); atom_a_iter != residuesData[i].m_atoms.end(); ++atom_a_iter) {
			if (atom_a_iter->second.m_isMainChain){
				continue;
			}
			//b atom main chain
			for (int j = 0; j < length; ++j){
				//exclude same id
				if (i == j){
					continue;
				}
				for (auto atom_b_iter = residuesData[j].m_atoms.begin(); atom_b_iter != residuesData[j].m_atoms.end(); ++atom_b_iter) {
					if (!atom_b_iter->second.m_isMainChain){
						continue;
					}
					//atom mightbe contact
					double dij2 = (atom_a_iter->second.m_position - atom_b_iter->second.m_position).squaredNorm();
					double aij2 = getAij2(atom_a_iter->second, atom_b_iter->second);
					double d_star2 = dij2 / aij2;
					if (d_star2 < 6.25){
						double lambd = 1.6;
						if (atom_a_iter->second.m_LJType == 6 && atom_b_iter->second.m_LJType == 6){
							lambd = 1;
						}
						double eij = sqrt(atom_a_iter->second.m_well*atom_b_iter->second.m_well);
						potentials += calLJPotential(d_star2, eij, lambd);
					}
				}
			
			}
		}	
	}
	return potentials;
}

double getMMLJPotential(const vector<Residue>& residuesData){

	double potentials = 0;

	int length = residuesData.size();
	for (int i = 0; i < length; ++i){
		//a atom main chain
		for (auto atom_a_iter = residuesData[i].m_atoms.begin(); atom_a_iter != residuesData[i].m_atoms.end(); ++atom_a_iter) {
			if (!atom_a_iter->second.m_isMainChain){
				continue;
			}
			//b atom main chain
			for (int j = i + 2; j < length; ++j){
				for (auto atom_b_iter = residuesData[j].m_atoms.begin(); atom_b_iter != residuesData[j].m_atoms.end(); ++atom_b_iter) {
					if (!atom_b_iter->second.m_isMainChain){
						continue;
					}
					//atom mightbe contact
					double dij2 = (atom_a_iter->second.m_position - atom_b_iter->second.m_position).squaredNorm();
					double aij2 = getAij2(atom_a_iter->second, atom_b_iter->second);
					double d_star2 = dij2 / aij2;
					if (d_star2 < 6.25){
						double lambd = 1.6;
						if (atom_a_iter->second.m_LJType == 6 && atom_b_iter->second.m_LJType == 6){
							lambd = 1;
						}
						double eij = sqrt(atom_a_iter->second.m_well*atom_b_iter->second.m_well);
						potentials += calLJPotential(d_star2, eij, lambd);
					}
				}

			}
		}
	}
	return potentials;
}

double getSSLJPotential(const vector<Residue>& residuesData){

	double potentials = 0;

	int length = residuesData.size();
	for (int i = 0; i < length; ++i){
		//a atom main chain
		for (auto atom_a_iter = residuesData[i].m_atoms.begin(); atom_a_iter != residuesData[i].m_atoms.end(); ++atom_a_iter) {
			if (atom_a_iter->second.m_isMainChain){
				continue;
			}
			//b atom side chain
			for (int j = i + 1; j < length; ++j){
				for (auto atom_b_iter = residuesData[j].m_atoms.begin(); atom_b_iter != residuesData[j].m_atoms.end(); ++atom_b_iter) {
					if (atom_b_iter->second.m_isMainChain){
						continue;
					}
					//atom mightbe contact
					double dij2 = (atom_a_iter->second.m_position - atom_b_iter->second.m_position).squaredNorm();
					double aij2 = getAij2(atom_a_iter->second, atom_b_iter->second);
					double d_star2 = dij2 / aij2;
					if (d_star2 < 6.25){
						double lambd = 1.6;
						if (atom_a_iter->second.m_LJType == 6 && atom_b_iter->second.m_LJType == 6){
							lambd = 1;
						}
						double eij = sqrt(atom_a_iter->second.m_well*atom_b_iter->second.m_well);
						potentials += calLJPotential(d_star2, eij, lambd);
					}
				}

			}
		}
	}
	return potentials;
}