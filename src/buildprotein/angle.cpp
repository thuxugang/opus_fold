#include "angle.h"

double getNorm(const Vector3d& v){
	//Return vector norm
	return sqrt(v.dot(v));
}

double getAngle(const Vector3d& v1, const Vector3d& v2){
	//Return angle between two vectors
	double n1 = getNorm(v1);
	double n2 = getNorm(v2);
	double c = (v1.dot(v2)) / (n1 * n2);
	c = c>1 ? 1 : c;
	c = c<-1 ? -1 : c;
	return acos(c) / M_PI * 180;
}

double getBondLength(const Vector3d& v1, const Vector3d& v2){
	return getNorm(v1 - v2);
}

double getBondAngle(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3){
	//Calculate the angle between 3 vectors
	Vector3d v4 = v1 - v2;
	Vector3d v5 = v3 - v2; 

	return getAngle(v4, v5);
} 

double calDihedral(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3, const Vector3d& v4){
	//Return dihedral between four vectors
	Vector3d ab = v1 - v2;
	Vector3d cb = v3 - v2;
	Vector3d db = v4 - v3;
	Vector3d u = ab.cross(cb);
	Vector3d v = db.cross(cb);
	Vector3d w = u.cross(v);

	double angle = getAngle(u, v);

	//Determine sign of angle
	if (getAngle(cb, w) > 0.001){
		angle = -angle;
	}

	return angle;
}


//
//double calDihedral(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3, const Vector3d& v4){
//	//Return dihedral between four vectors
//	Vector3d ab = v1 - v2;
//	Vector3d cb = v3 - v2;
//	Vector3d db = v4 - v3;
//	Vector3d u = ab.cross(cb);
//	Vector3d v = db.cross(cb);
//	Vector3d w = u.cross(v);
//
//	double angle = getAngle(u, v);
//
//	//Determine sign of angle
//	if (getAngle(cb, w) > 0.001){
//		angle = -angle;
//	}
//
//	return angle;
//}
//
//double calAngle(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3){
//	//Calculate the angle between 3 vectors
//	Vector3d v4 = v1 - v2;
//	Vector3d v5 = v3 - v2; 
//
//	return getAngle(v4, v5);
//} 

MatrixXd rotaxis(double theta, const Vector3d& v){
	
	Vector3d v_norm = v/getNorm(v);
	double c = cos(theta);
	double s = sin(theta); 
	double t = 1 - c;
	double x = v_norm[0];
	double y = v_norm[1];
	double z = v_norm[2];

	MatrixXd rot(3,3);
	rot(0,0) = t * x * x + c;
	rot(0,1) = t * x * y - s * z;
	rot(0,2) = t * x * z + s * y; 
	rot(1,0) = t * x * y + s * z; 
	rot(1,1) = t * y * y + c; 
	rot(1,2) = t * y * z - s * x; 
	rot(2,0) = t * x * z - s * y; 
	rot(2,1) = t * y * z + s * x; 
	rot(2,2) = t * z * z + c; 
	return rot;
}

Vector3d calCoordinates(const Atom& refA, const Atom& refB, const Atom& refC, double L, double ang, double di){

	Vector3d AV = refA.m_position;
	Vector3d BV = refB.m_position;
	Vector3d CV = refC.m_position;

	Vector3d CA = AV - CV;
	Vector3d CB = BV - CV;

	//CA vector
	double AX = CA[0];
	double AY = CA[1];
	double AZ = CA[2];

	//CB vector
	double BX = CB[0];
	double BY = CB[1];
	double BZ = CB[2];

	//Plane Parameters
	double A = (AY*BZ)-(AZ*BY);
	double B = (AZ*BX)-(AX*BZ);
	double G = (AX*BY)-(AY*BX);

	//Dot Product Constant
	double F = sqrt(BX*BX + BY*BY + BZ*BZ) * L * cos(ang*(M_PI/180.0));

	//Constants
	double cons = sqrt(pow((B*BZ-BY*G),2) *(-(F*F)*(A*A+B*B+G*G)+(B*B*(BX*BX+BZ*BZ) + A*A*(BY*BY+BZ*BZ)- (2*A*BX*BZ*G) + (BX*BX+ BY*BY)*G*G - (2*B*BY)*(A*BX+BZ*G))*L*L));
	double denom = (B*B)*(BX*BX+BZ*BZ)+ (A*A)*(BY*BY+BZ*BZ) - (2*A*BX*BZ*G) + (BX*BX+BY*BY)*(G*G) - (2*B*BY)*(A*BX+BZ*G);

	double X = ((B*B*BX*F)-(A*B*BY*F)+(F*G)*(-A*BZ+BX*G)+cons)/denom;

	double Y,Z;
	if((B==0 || BZ==0) && (BY==0 || G==0)){
		double cons1 = sqrt(G*G*(-A*A*X*X+(B*B+G*G)*(L-X)*(L+X)));
		Y = ((-A*B*X)+cons1)/(B*B+G*G);
		Z = -(A*G*G*X+B*cons1)/(G*(B*B+G*G));
	}else{
		Y = ((A*A*BY*F)*(B*BZ-BY*G)+ G*( -F*pow(B*BZ-BY*G,2) + BX*cons) - A*( B*B*BX*BZ*F- B*BX*BY*F*G + BZ*cons)) / ((B*BZ-BY*G)*denom);
		Z = ((A*A*BZ*F)*(B*BZ-BY*G) + (B*F)*pow(B*BZ-BY*G,2) + (A*BX*F*G)*(-B*BZ+BY*G) - B*BX*cons + A*BY*cons) / ((B*BZ-BY*G)*denom);

	}

	//get the new vector from the origin
	Vector3d D(X+CV[0],Y+CV[1],Z+CV[2]);

	double temp = calDihedral(AV, BV, CV, D);
	di = di-temp;
	MatrixXd rot = rotaxis(M_PI*(di/180.0), CV-BV);
	D = rot*(D-BV)+BV;

	return D;
}

Vector3d calOCoordinates(const Atom& refA, const Atom& refB, const Atom& refC, double L){

	Vector3d ca = refA.m_position;
	Vector3d c = refB.m_position;
	Vector3d n1 = refC.m_position;

	Vector3d c_ca = ca - c;
	Vector3d c_n1 = n1 - c;

	//assume 120
	Vector3d c_o_dir = -(c_ca / getNorm(c_ca) + c_n1 / getNorm(c_n1));

	Vector3d c_o_u = c_o_dir / getNorm(c_o_dir);

	return c_o_u * L + c;
}

//
//string rotamerFormat(double anlge){
//	int angle_i;
//	if(anlge>0){
//		angle_i = ((int)(anlge/10+0.5))*10;
//	}else{
//		angle_i = ((int)(anlge/10-0.5))*10;
//	}
//	return to_string(angle_i);
//}
//
//void addPhiPsi(vector<Residue>& residuesData){
//	int length = residuesData.size();
//	for(int i=0; i<length; i++){
//		Residue* res = &residuesData[i];
//		if(i == 0){
//			res->m_phi = -60;
//		}else{
//			Residue* _res = &residuesData[i-1];
//			res->m_phi = calDihedral(_res->getAtom("C").m_position, res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("C").m_position);
//		}
//		if(i == length-1){
//			res->m_psi = 60;
//		}else{
//			Residue* res_ = &residuesData[i+1];
//			res->m_psi = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("C").m_position, res_->getAtom("N").m_position);
//		}
//		string phi = rotamerFormat(res->m_phi);
//		string psi = rotamerFormat(res->m_psi);
//		res->m_phipsi = getTriResname(res->m_resname)+"_"+phi+"_"+psi;
//	}
//};
//
//void addDihedral(vector<Residue>& residuesData){
//	int length = residuesData.size();
//	for(int i=0; i<length; i++){
//		Residue* res = &residuesData[i];
//		double x1;
//		double x2;
//		double x3;
//		double x4;
//		try{
//			switch(res->m_resname){
//				case('G'):
//					break;
//				case('A'):
//					break;
//				case('S'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("OG").m_position);
//					res->m_dihedrals.push_back(x1);
//					break;
//				case('C'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("SG").m_position);
//					res->m_dihedrals.push_back(x1);
//					break;
//				case('V'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG1").m_position);
//					res->m_dihedrals.push_back(x1);
//					break;
//				case('I'):
// 					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG1").m_position);
//					res->m_dihedrals.push_back(x1);
//					//CD1
//					try{
//						x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG1").m_position, res->getAtom("CD").m_position);
//						res->m_dihedrals.push_back(x2);
//					}catch(...){
//						x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG1").m_position, res->getAtom("CD1").m_position);
//						res->m_dihedrals.push_back(x2);
//					}
//					break;
//				case('L'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD1").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2);
//					break;
//				case('T'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("OG1").m_position);
//					res->m_dihedrals.push_back(x1);
//					break;
//				case('R'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD").m_position);
//					x3 = calDihedral(res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD").m_position, res->getAtom("NE").m_position);
//					x4 = calDihedral(res->getAtom("CG").m_position, res->getAtom("CD").m_position, res->getAtom("NE").m_position, res->getAtom("CZ").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2);
//					res->m_dihedrals.push_back(x3);
//					res->m_dihedrals.push_back(x4);
//					break;
//				case('K'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD").m_position);
//					x3 = calDihedral(res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD").m_position, res->getAtom("CE").m_position);
//					x4 = calDihedral(res->getAtom("CG").m_position, res->getAtom("CD").m_position, res->getAtom("CE").m_position, res->getAtom("NZ").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2);
//					res->m_dihedrals.push_back(x3);
//					res->m_dihedrals.push_back(x4);
//					break;
//				case('D'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("OD1").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2);    
//					break;
//				case('E'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD").m_position);
//					x3 = calDihedral(res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD").m_position, res->getAtom("OE1").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2);
//					res->m_dihedrals.push_back(x3);
//					break;
//				case('N'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("OD1").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2); 
//					break;
//				case('Q'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD").m_position);
//					x3 = calDihedral(res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD").m_position, res->getAtom("OE1").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2);
//					res->m_dihedrals.push_back(x3);
//					break;
//				case('M'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("SD").m_position);
//					x3 = calDihedral(res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("SD").m_position, res->getAtom("CE").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2);
//					res->m_dihedrals.push_back(x3);
//					break;
//				case('H'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("ND1").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2); 
//					break;
//				case('P'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2); 
//					break;
//				case('F'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD1").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2); 
//					break;		
//				case('Y'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD1").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2); 
//					break;
//				case('W'):
//					x1 = calDihedral(res->getAtom("N").m_position, res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position);
//					x2 = calDihedral(res->getAtom("CA").m_position, res->getAtom("CB").m_position, res->getAtom("CG").m_position, res->getAtom("CD1").m_position);
//					res->m_dihedrals.push_back(x1);
//					res->m_dihedrals.push_back(x2); 
//					break;
//				default:
//					cout << res->m_resname << endl;
//					throw "angle.calDihedral() wrong";
//			}
//		}catch(...){
//			cout << "calDihedral wrong at: " + to_string(res->m_resid) + " " + res->m_resname  << endl;
//		}
//	}
//
//};

/*
int main(){

	Vector3d v1(23.729,  30.406,  37.362);
	Vector3d v2(25.075,  29.663,  37.446);
	Vector3d v3(25.526,  29.691,  38.860);
	Vector3d v4(26.474,  28.873,  39.300);
	
	double angle1 = calDihedral(v1, v2, v3, v4);
	//double angle2 = calAngle(v1, v2, v3);
	//MatrixXd rot = rotaxis(45, v1);
	//Vector3d v5 = rot*v4;
	cout << angle1 << endl;
	//cout << angle2 << endl;
	//cout << rot << endl;
	//cout << v5[0] <<  v5[1] << v5[2]<< endl;

	

	return 0;
}
*/

/*
int main(){
	Vector3d p1(-5, -1, 1);
	Vector3d p2(0, 6, 0);
	Vector3d p3(1, 0, 7);
	Atom a(0, "N", 'G', 0, p1);
	Atom b(1, "C", 'G', 0, p2);
	Atom c(2, "CG", 'G', 0, p3);

	Vector3d v = calCoordinates(a, b, c, 1, 2, 3);

	cout << v << endl;
	
	

	return 0;

}*/