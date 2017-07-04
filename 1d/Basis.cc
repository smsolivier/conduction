#include "Basis.hh"
#include "helper.hh"
#include "LinearSolver.hh"

#include <cmath>

Basis::Basis(vector<double> &xloc) : p(xloc.size()) {

	MatrixResize(coef, p); 

	for (int k=0; k<p; k++) {

		vector<vector<double>> A(p, vector<double>(p)); // store system for each basis function 

		for (int i=0; i<p; i++) {

			for (int j=0; j<p; j++) {

				A[i][j] = pow(xloc[i], j); 

			}
		}

		vector<double> b(p); // rhs 

		b[k] = 1; 

		int err = gauss_elim(p, A, coef[k], b); 

	}

}

double Basis::B(int i, double xi) {

	double sum = 0; 
	for (int j=0; j<p; j++) {

		sum += coef[i][j] * pow(xi, j); 

	}

	return sum; 

}

double Basis::dB(int i, double xi) {

	double sum = 0; 
	for (int j=1; j<p; j++) {

		sum += j*coef[i][j] * pow(xi, j-1); 
	}

	return sum; 

}

// int main() {

// 	vector<double> xloc = {-1, 0, 1}; 

// 	Basis basis(xloc); 

// 	cout << basis.dB(0, 1) << endl; 

// }