#include "JacobiSolver.hh"
#include "VectorMath.hh"

#include <iostream>
#include <cassert>

// vector<double> operator+(const vector<double> &a, const vector<double> &b) {

// 	int N = a.size(); 

// 	// assert same size 
// 	assert(N == b.size()); 

// 	vector<double> c(N); 

// 	for (int i=0; i<N; i++) {

// 		c[i] = a[i] + b[i]; 

// 	}

// 	return c; 

// }

// vector<double> operator-(const vector<double> &a, const vector<double> &b) {

// 	int N = a.size(); 

// 	// assert same size 
// 	assert(N == b.size()); 

// 	vector<double> c(N); 

// 	for (int i=0; i<N; i++) {

// 		c[i] = a[i] - b[i]; 

// 	}

// 	return c; 

// }

// vector<double> operator*(const vector<double> &a, const double alpha) {

// 	vector<double> b(a.size()); 

// 	for (int i=0; i<a.size(); i++) {

// 		b[i] = a[i] * alpha; 

// 	}

// 	return b; 

// }

double max(const vector<double>& a_a) {

	int N = a_a.size(); 

	double max = 0.; 

	for (int i=0; i<N; i++) {

		if (abs(a_a[i]) > max) max = abs(a_a[i]); 

	}

	return max; 

}

double JacobiSolver::solve(
	vector<double>& a_phi, 
	const SparseMatrix& a_A, 
	const vector<double>& a_rhs, 
	const double& a_tolerance, 
	int a_iter) 
{

	double alpha = .9; 

	// find largest diagonal 
	double max_diag = 0; 
	array<int, 2> loc = {0, 0}; 
	for (int i=0; i<a_rhs.size(); i++) {

		loc[0] = i; 
		loc[1] = i; 

		if (a_A[loc] > max_diag) max_diag = a_A[loc]; 

	}

	double relax = alpha / max_diag; 

	vector<double> phi_old(a_rhs.size(), 0); 
	a_phi.resize(a_rhs.size()); 

	double R; 

	const double max_rhs = max(a_rhs); 

	cout << "norm(rhs) = " << max_rhs << endl; 

	int i; 
	for (i=0; i<a_iter; i++) {

		a_phi = phi_old + (a_rhs - (a_A*phi_old))*relax; 

		vector<double> resid = (a_A*a_phi) - a_rhs; 

		R = max(resid)/max_rhs;

		cout << R << "  " << i << " \r";
		cout.flush(); 

		if (R < a_tolerance) break; 

		phi_old = a_phi; 

	}
	
	cout << "final tolerance: " << R << endl; 
	cout << "iter to converge: " << i << endl; 
	return R; 

}

