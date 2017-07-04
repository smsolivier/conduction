#include "Basis.hh"
#include "helper.hh"
#include "LinearAlgebra.hh"
#include "LinearSolver.hh"

#include <iostream>
#include <fstream>
#include <cmath>

Basis::Basis(vector<vector<double>> &pts) : pts(pts) {

	p = pow(pts.size(), .5) - 1; 

	N = pts.size(); 

	x = linspace(-1, 1, p+1); 
	y = linspace(-1, 1, p+1); 

	xind.resize(p+1); 
	yind.resize(p+1); 
	for (int i=0; i<p+1; i++) {

		xind[i] = i; 
		yind[i] = i; 

	}

	MatrixResize(box, N, 2);
	ind.resize(N); 
	for (int i=0; i<p+1; i++) {

		ind[i].resize(2); 

		for (int j=0; j<p+1; j++) {

			box[i*(p+1)+j] = {x[j], y[i]}; 

			ind[i*(p+1)+j] = {xind[j], yind[i]}; 

		}

	}

	MatrixResize(coef, p+1); 

	// solve for coefficients 
	for (int i=0; i<p+1; i++) {

		vector<vector<double>> A; 
		MatrixResize(A, p+1); 

		// build rows of A
		for (int j=0; j<p+1; j++) {

			for (int k=0; k<p+1; k++) {

				A[j][k] = pow(x[j], k); 

			}

		}

		vector<double> b(p+1, 0); 
		b[i] = 1; 

		int status = gauss_elim(p+1, A, coef[i], b); 

	}

	for (int i=0; i<p+1; i++) {

		basis.push_back(Polynomial(coef[i])); 
		dbasis.push_back(basis[i].derivative()); 

	}

	// make sure basis is correct 
	// for (int i=0; i<N; i++) {

	// 	for (int j=0; j<N; j++) {

	// 		cout << B(i, box[j][0], box[j][1]) << endl; 
	// 	}

	// 	cout << endl; 

	// }

}

double Basis::B(int i, double xi, double eta) {

	return basis[ind[i][0]].evaluate(xi) * basis[ind[i][1]].evaluate(eta); 

}

double Basis::dBxi(int i, double xi, double eta) {

	return dbasis[ind[i][0]].evaluate(xi) * basis[ind[i][1]].evaluate(eta); 

}

double Basis::dBeta(int i, double xi, double eta) {

	return basis[ind[i][0]].evaluate(xi) * dbasis[ind[i][1]].evaluate(eta); 

}

double Basis::dBx(int i, double xi, double eta) {

	return dBxi(i, xi, eta) * J_in(0, 0, xi, eta) + dBeta(i, xi, eta) * J_in(0, 1, xi, eta); 

}

double Basis::dBy(int i, double xi, double eta) {

	return dBxi(i, xi, eta) * J_in(1, 0, xi, eta) + dBeta(i, xi, eta) * J_in(1, 1, xi, eta);

}

double Basis::Jacobian(int i, int j, double xi, double eta) {

	double sum = 0; 

	// use xi derivative 
	if (j == 0) {

		for (int ii=0; ii<N; ii++) {

			sum += dBxi(ii, xi, eta) * pts[ii][i]; 

		}

	}

	// use eta derivative 
	else {

		for (int ii=0; ii<N; ii++) {

			sum += dBeta(ii, xi, eta) * pts[ii][i]; 

		}

	}

	return sum; 

}

double Basis::detJ(double xi, double eta) {

	double det = Jacobian(0, 0, xi, eta) * Jacobian(1, 1, xi, eta) - 
		Jacobian(0, 1, xi, eta) * Jacobian(1, 0, xi, eta); 

	if (det == 0) {

		cout << "BIG TROUBLE: determinant of jacobian = 0" << endl; 
		exit(0); 

	}

	return det; 

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

void plotBasis() {

	vector<vector<double>> pts = {{0,0}, {.5, 0}, {1,0}, 
		{0,.5}, {.5, .5}, {.5, 1}, 
		{0,1}, {.5, 1}, {1, 1}}; 

	Basis basis(pts); 

	int N = 250; 

	vector<double> x = linspace(-1, 1, N); 

	ofstream file; 
	file.open("basis"); 
	ofstream X; 
	X.open("X"); 
	ofstream Y; 
	Y.open("Y"); 
	for (int i=0; i<N; i++) {

		for (int j=0; j<N; j++) {

			file << basis.B(8, x[i], x[j]) << " "; 
			X << x[i] << " "; 
			Y << x[j] << " "; 

		}

		file << endl; 
		X << endl; 
		Y << endl; 

	}

	file.close(); 
	X.close();
	Y.close(); 

}

void testJacobian() {

	vector<vector<double>> pts = {{0,0}, {1,0}, {0,1}, {1,1}}; 

	Basis basis(pts); 

	vector<vector<double>> J, J_in; 
	MatrixResize(J, 2); 
	MatrixResize(J_in, 2); 

	double xi = 0; 
	double eta = 0; 

	cout << basis.detJ(xi, eta) << endl; 

	for (int i=0; i<2; i++) {

		for (int j=0; j<2; j++) {

			J[i][j] = basis.Jacobian(i, j, xi, eta); 
			J_in[i][j] = basis.J_in(i, j, xi, eta); 

		}
	}

	printVector(J_in*J); 

}

// int main() {

// 	testJacobian(); 

// }