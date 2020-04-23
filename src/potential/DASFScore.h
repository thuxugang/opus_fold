#ifndef DASFSCORE_H
#define DASFSCORE_H

#include <Eigen/Dense> 
#include <vector> 
#include "residue.h"

using namespace std;
using namespace Eigen;

class DASFReference{
public:
	DASFReference();

	Vector3d ref;
	Matrix4d rotation_matrix;
};

DASFReference getDASFReference(const Residue& residue);

double getDASFPotential(const Residue& residue, const DASFReference& reference);

#endif