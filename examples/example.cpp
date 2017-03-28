#include<iostream>
#include"./../header/SS.hpp"
#include<Eigen/Dense>

int main(int argc, char* argv[]) {
	int N		=	atoi(argv[1]);
	int p		=	atoi(argv[2]);
	Eigen::MatrixXd U	=	Eigen::MatrixXd::Random(N,p);
	Eigen::MatrixXd V	=	Eigen::MatrixXd::Random(N,p);
	Eigen::VectorXd d	=	Eigen::VectorXd::Random(N);
	SS* matrix 			=	new SS(N,p,U,V,d);
	Eigen::VectorXd x	=	matrix->Eigen::VectorXd::Random(N);
	Eigen::VectorXd b	=	matrix->mat_vec_prod(x);
	matrix->obtain_Cholesky();
	
	delete matrix;
}