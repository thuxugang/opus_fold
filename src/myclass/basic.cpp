#include "basic.h"

MultivariateNorm::MultivariateNorm(){
}

MultivariateNorm::MultivariateNorm(VectorXd mean, MatrixXd cov){

	m_mean = mean;
	m_cov = cov;
	m_dimension = m_mean.size();

	// Find the eigen vectors of the covariance matrix
	SelfAdjointEigenSolver<MatrixXd> eigen_solver(m_cov);
	MatrixXd eigenvectors = eigen_solver.eigenvectors().real();

	// Find the eigenvalues of the covariance matrix
	MatrixXd eigenvalues = eigen_solver.eigenvalues().real().asDiagonal();

	// Find the transformation matrix
	SelfAdjointEigenSolver<MatrixXd> es(eigenvalues);
	MatrixXd sqrt_eigenvalues = es.operatorSqrt();
	m_Q = eigenvectors * sqrt_eigenvalues;
}

CSFFeature::CSFFeature(const vector<double>& mean, const vector<double>& sd, int nums, int window_len){
	m_mean = Vector3d(mean[0], mean[1], mean[2]);
	m_sd = Vector3d(sd[0], sd[1], sd[2]);
	m_nums = nums;
	m_window_len = window_len;
}

TAFeature::TAFeature(const vector<double>& mean, const vector<double>& sd, double possibility, int window_len, double total_possibility){
	m_mean = Vector2d(mean[0], mean[1]);
	m_covar << sd[0], sd[1],
		sd[2], sd[3];
	m_possibility = possibility;
	m_window_len = window_len;
	m_total_possibility = total_possibility;
	m_gmm_sampler = MultivariateNorm(m_mean, m_covar);
}

DASFFeature::DASFFeature(const vector<double>& mean, const vector<double>& sd, int nums, int window_len){
	m_mean = Vector3d(mean[0], mean[1], mean[2]);
	m_sd = Vector3d(sd[0], sd[1], sd[2]);
	m_nums = nums;
	m_window_len = window_len;
}

RotamerFeature::RotamerFeature(){
}

RotamerFeature::RotamerFeature(double prob, double x1, double x2, double x3, double x4){
	m_prob = prob;
	m_dihedral.push_back(x1);
	m_dihedral.push_back(x2);
	m_dihedral.push_back(x3);
	m_dihedral.push_back(x4);
}

