//	SS.hpp
//	Created by Sivaram Ambikasaran on March 28th, 2017

#ifndef __SS_HPP__
#define __SS_HPP__

#include<vector>
#include<Eigen/Dense>
#include<Eigen/Sparse>

class SS {
	int N;				//	Number of unknowns
	int p;				//	Rank of the semi-separable part
	//	The semi-sparable matrix is of the form diag(d) + tril(U*V',-1) + triu(V*U',1)
	Eigen::MatrixXd U;
	Eigen::MatrixXd V;
	Eigen::VectorXd d;	//	Diagonal of the matrix
	//	The Cholesky factorization is of the form diag(d_Chol) + tril(U X_chol,-1)
	bool cholesky_obtained;
	Eigen::MatrixXd X_chol;
	Eigen::VectorXd d_chol;
public:
	SS(const int N, const int p, const Eigen::MatrixXd U, const Eigen::MatrixXd V, const Eigen::VectorXd d);
	void change_Diagonal(Eigen::VectorXd d);

	void Cholesky_compute();	//	Returns a lower-triangular matrix L such that LL' = A, where L = diag(d_Chol) + tril(U X_chol,-1)

	Eigen::VectorXd mat_vec_prod(Eigen::VectorXd x);	//	Returns b	=	Ax
	Eigen::VectorXd Cholesky_solve(Eigen::VectorXd b);			//	Returns x 	= 	A^{-1}b
	Eigen::VectorXd lower_tri_solve(Eigen::VectorXd d, Eigen::MatrixXd U, Eigen::MatrixXd V, Eigen::VectorXd b); //	Returns x 	= 	A^{-1}b, where A = diag(d) + tril(U*V',-1)
	Eigen::VectorXd upper_tri_solve(Eigen::VectorXd d, Eigen::MatrixXd U, Eigen::MatrixXd V, Eigen::VectorXd b); //	Returns x 	= 	A^{-1}b, where A = diag(d) + triu(V*U',-1)
	Eigen::VectorXd apply_factor(Eigen::VectorXd x);	//	Returns b 	=	Lx
	Eigen::VectorXd generate_Realization();				//	Returns a realization having the SS matrix as its covariance

	Eigen::MatrixXd mat_mat_prod(Eigen::MatrixXd x);	//	Returns B	=	AX
	Eigen::MatrixXd solve(Eigen::MatrixXd b);			//	Returns X 	= 	A^{-1}B
	Eigen::MatrixXd lower_tri_solve(Eigen::VectorXd d, Eigen::MatrixXd U, Eigen::MatrixXd V, Eigen::MatrixXd b); //	Returns X 	= 	A^{-1}B, where A = diag(d) + tril(U*V',-1)
	Eigen::MatrixXd upper_tri_solve(Eigen::VectorXd d, Eigen::MatrixXd U, Eigen::MatrixXd V, Eigen::MatrixXd b); //	Returns X 	= 	A^{-1}B, where A = diag(d) + triu(V*U',-1)
	Eigen::MatrixXd apply_factor(Eigen::MatrixXd x);	//	Returns B 	=	LX
	Eigen::MatrixXd generate_Realization(int n);		//	Returns 'n' realizations having the SS matrix as its covariance
};

SS::SS(const int N, const int p, const Eigen::MatrixXd U, const Eigen::MatrixXd V, const Eigen::VectorXd d) {
	/************************************************/
	/*	Assign all the variables inside the class	*/
	/************************************************/
	this->N	=	N;
	this->p	=	p;
	this->U	=	U;
	this->V	=	V;
	change_Diagonal(d);
	d_chol	=	Eigen::VectorXd::Zero(p);
	X_chol	=	Eigen::MatrixXd::Zero(N,p);
	cholesky_obtained	=	false;
}

void SS::change_Diagonal(Eigen::VectorXd d) {
	this->d	=	d;
}

Eigen::VectorXd SS::mat_vec_prod(Eigen::VectorXd x) {
	Eigen::MatrixXd b	=	d.asDiagonal()*x;
	Eigen::MatrixXd v	=	Eigen::MatrixXd::Zero(p,N);
	Eigen::MatrixXd y	=	Eigen::MatrixXd::Zero(p,N);
	for (int k=0; k<p; ++k) {
		v(k,0)				=	0;
		y(k,N-1)			=	0;
		for (int i=1; i<N; ++i) {
			v(k,i)		=	v(k,i-1)+V(i-1,k)*x(i-1);
		}
		for (int i=N-2; i>-1; --i) {
			y(k,i)		=	y(k,i+1)+U(i+1,k)*x(i+1);
		}
	}
	for (int i=0; i<N; ++i) {
		for (int k=0; k<p; ++k) {
			b(i)+=(U(i,k)*v(k,i)+V(i,k)*y(k,i));
		}
	}
	return b;
}

void SS::Cholesky_compute() {
	Eigen::MatrixXd S	=	Eigen::MatrixXd::Zero(p,p);
	for (int i=0; i<N; ++i) {
		temp			=	U.row(i)*S;
		d_chol(i)		=	sqrt(d(i)-temp*U.row(i),transpose());
		X.row(i)		=	(V.row(i)-temp)/d_Chol(i);
		S+=(X.row(i).transpose()*X.row(i));
	}
	cholesky_obtained	=	true;
}

Eigen::VectorXd SS:lower_tri_solve(Eigen::VectorXd d, Eigen::MatrixXd U, Eigen::MatrixXd V, Eigen::VectorXd b) {
	Eigen::VectorXd x	=	Eigen::VectorXd::Zero(b.size());
	x(0)				=	b(0)/d(0);
	Eigen::VectorXd z	=	Eigen::VectorXd::Zero(V.cols());
	for (int i=1; i<N; ++i) {
		z+=x(i-1)*V.row(i-1).transpose();
		x(i)			=	(b(i)-U.row(i)*z)/d(i);
	}
	return x;
}

Eigen::VectorXd SS::Cholesky_solve(Eigen::VectorXd b) {
	if (cholesky_obtained == false) {
		obtain_Cholesky();
	}
	Eigen::VectorXd x_temp	=	lower_tri_solve(d_chol, U, X, b);
	Eigen::VectorXd x		=	upper_tri_solve(d_chol, U, V, x_temp);
	return x;
}

#endif /*(defined __SS_HPP__)*/