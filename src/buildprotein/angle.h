#ifndef ANGLE_H
#define ANGLE_H

#include <iostream> 
#include <string> 
#include <vector>
#include <Eigen/Dense> 

#include "residue.h"

using namespace std;
using namespace Eigen;

//const double M_PI = 3.14159265358979323846;

//double getNorm(const Vector3d& v);
//double getAngle(const Vector3d& v1, const Vector3d& v2);
//double calDihedral(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3, const Vector3d& v4);
//double calAngle(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3);
//MatrixXd rotaxis(double theta, const Vector3d& v);
//
//void addDihedral(vector<Residue>& residuesData);

double getBondLength(const Vector3d& v1, const Vector3d& v2);

double getBondAngle(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3);

double calDihedral(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3, const Vector3d& v4);

Vector3d calCoordinates(const Atom& refA, const Atom& refB, const Atom& refC, double L, double ang, double di);
Vector3d calOCoordinates(const Atom& refA, const Atom& refB, const Atom& refC, double L);

#endif