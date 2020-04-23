#ifndef BASIC_H
#define BASIC_H

#include <vector>
#include <Eigen/Dense> 

using namespace std;
using namespace Eigen;

class CSFFeature{
public:
	Vector3d m_mean;
	Vector3d m_sd;
	int m_nums;
	int m_window_len;

	CSFFeature(const vector<double>& mean, const vector<double>& sd, int nums, int window_len);
};

class MultivariateNorm{
public:
	MultivariateNorm();
	MultivariateNorm(VectorXd mean, MatrixXd sigma);
	MatrixXd m_Q;

	VectorXd m_mean;
	MatrixXd m_cov;
	int m_dimension;
};

class TAFeature{
public:
	Vector2d m_mean;
	Matrix2d m_covar;
	double m_possibility;
	int m_window_len;

	double m_total_possibility;

	MultivariateNorm m_gmm_sampler;

	TAFeature(const vector<double>& mean, const vector<double>& sd, double nums, int window_len, double total_possibility = 1);

};

class DASFFeature{
public:
	Vector3d m_mean;
	Vector3d m_sd;
	int m_nums;
	int m_window_len;

	DASFFeature(const vector<double>& mean, const vector<double>& sd, int nums, int window_len);
};

class RotamerFeature{
public:

	double m_prob;
	vector<double> m_dihedral;

	RotamerFeature();
	RotamerFeature(double prob, double x1, double x2, double x3, double x4);
};

#endif