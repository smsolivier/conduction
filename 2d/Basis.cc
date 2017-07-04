#include "Basis.hh"
#include "helper.hh"
#include "LinearAlgebra.hh"
#include "LinearSolver.hh"

#include <iostream>
#include <fstream>

Basis::Basis(vector<vector<double>> &pts) : pts(pts) {

	p = pts.size(); // number of points 

	MatrixResize(coef, p); 

	MatrixResize(box, p, 2);
	box[0] = {-1, -1}; 
	box[1] = {1, -1};
	box[2] = {-1, 1};  
	box[3] = {1, 1}; 

	// solve for coefficients 
	for (int i=0; i<p; i++) {

		vector<vector<double>> A; 
		MatrixResize(A, p); 

		for (int j=0; j<p; j++) {

			A[j][0] = 1; 
			A[j][1] = box[j][0]; 
			A[j][2] = box[j][1]; 
			A[j][3] = box[j][0] * box[j][1];

		}

		vector<double> b(p, 0); 
		b[i] = 1; 

		int status = gauss_elim(p, A, coef[i], b); 

	}	

}

double Basis::B(int i, double xi, double eta) {

	return coef[i][0]  + coef[i][1]*xi + coef[i][2]*eta + coef[i][3]*xi*eta; 

}

double Basis::dBxi(int i, double xi, double eta) {

	return coef[i][1] + coef[i][3]*eta; 

}

double Basis::dBeta(int i, double xi, double eta) {

	return coef[i][2] + coef[i][3]*xi; 

}

double Basis::dBx(int i, double xi, double eta) {

	return dBxi(i, xi, eta) * J_in(0, 0, xi, eta) + dBeta(i, xi, eta) * J_in(0, 1, xi, eta); 

}

double Basis::dBy(int i, double xi, double eta) {

	return dBxi(i, xi, eta) * J_in(1, 0, xi, eta) + dBeta(i, xi, eta) * J_in(1, 1, xi, eta);
	// return dBeta(i, xi, eta) * J_in(1, 1, xi, eta); 

}

double Basis::Jacobian(int i, int j, double xi, double eta) {

	double sum = 0; 

	// use xi derivative 
	if (j == 0) {

		for (int ii=0; ii<p; ii++) {

			sum += dBxi(ii, xi, eta) * pts[ii][i]; 

		}

	}

	// use eta derivative 
	else {

		for (int ii=0; ii<p; ii++) {

			sum += dBeta(ii, xi, eta) * pts[ii][i]; 

		}

	}

	return sum; 

}

double Basis::detJ(double xi, double eta) {

	return Jacobian(0, 0, xi, eta) * Jacobian(1, 1, xi, eta) - 
		Jacobian(0, 1, xi, eta) * Jacobian(1, 0, xi, eta);

}

double Basis::J_in(int i, int j, double xi, double eta) {

	double val; 
	if (i == 0 && j == 0) {

		val = Jacobian(1, 1, xi, eta); 

	} 

	else if (i == 0 && j == 1) {

		val = -1*Jacobian(0, 1, xi, eta); 

	}

	else if (i == 1 && j == 0) {

		val = -1*Jacobian(1, 0, xi, eta); 

	}

	else if (i == 1 && j == 1) {

		val = Jacobian(0, 0, xi, eta); 

	}

	return 1/detJ(xi, eta) * val; 

}

// int main() {

// 	vector<vector<double>> pts = {{0, 0}, {1, 0}, {0, 1}, {1,1}}; 

// 	Basis basis(pts); 

// 	double xi = -.5; 
// 	double eta = -.5; 

// 	vector<vector<double>> J, J_in;
// 	MatrixResize(J, 2); 
// 	MatrixResize(J_in, 2);  
// 	for (int i=0; i<2; i++) {

// 		for (int j=0; j<2; j++) {

// 			J[i][j] = basis.Jacobian(i, j, xi, eta); 
// 			J_in[i][j] = basis.J_in(i, j, xi, eta); 

// 		}

// 	}

// 	printVector(J * J_in); 

// }