#include "Element.hh"

#include "GaussQuad.hh"
#include "helper.hh"
#include "LinearAlgebra.hh"

#include <cmath>
#include <iostream>

Element::Element(vector<double> &box, vector<int> &globalNodes, int p) :
	globalNodes(globalNodes), p(p) {

	N = pow(p+1, 2); // number of nodes 

	int int_order = 10; 

	xglob = linspace(box[0], box[1], p+1); 
	yglob = linspace(box[2], box[3], p+1); 

	xloc = linspace(-1, 1, p+1); 
	yloc = linspace(-1, 1, p+1); 

	MatrixResize(pts, N, 2); 
	MatrixResize(pts_global, N, 2); 
	for (int i=0; i<yloc.size(); i++) {

		for (int j=0; j<xloc.size(); j++) {

			pts[i*yloc.size()+j] = {xloc[j], yloc[i]}; 
			pts_global[i*yloc.size()+j] = {xglob[j], yglob[i]}; 

		}

	}

	basis = new Basis(pts_global); // basis object 

	// generate A matrix (B_i * B_j * |J| )
	MatrixResize(A, N); 
	for (int i=0; i<N; i++ ) {

		for (int j=0; j<N; j++) {

			auto func = [this, i, j] (double xi, double eta) {

				return basis->B(i, xi, eta) * basis->B(j, xi, eta) * basis->detJ(xi, eta); 

			};

			A[i][j] = GaussQuad2D(func, int_order); 

		}

	}

	// generate B matrix (dB_i/dx * dB_j/dx * |J|)
	MatrixResize(B, N); 
	for (int i=0; i<N; i++) {

		for (int j=0; j<N; j++) {

			auto func = [this, i, j] (double xi, double eta) {

				return basis->dBx(i, xi, eta) * basis->dBx(j, xi, eta) * basis->detJ(xi, eta); 

			}; 

			B[i][j] = GaussQuad2D(func, int_order); 

			if (isnan(B[i][j])) cout << "integration issue in B" << endl; 

		}

	}

	// generate C matrix (dB_i/dy * dB_j/dy * |J|)
	MatrixResize(C, N); 
	for (int i=0; i<N; i++) {

		for (int j=0; j<N; j++) {

			auto func = [this, i, j] (double xi, double eta) {

				return basis->dBy(i, xi, eta) * basis->dBy(j, xi, eta) * basis->detJ(xi, eta); 

			}; 

			C[i][j] = GaussQuad2D(func, int_order); 

			if (isnan(C[i][j])) cout << "integration issue in C" << endl; 

		}

	}

	// generate D matrix (B_i * dB_j/dx)
	MatrixResize(D, N); 
	for (int i=0; i<N; i++) {

		for (int j=0; j<N; j++) {

			auto func = [this, i, j] (double xi, double eta) {

				return basis->dBx(j, xi, eta) * basis->B(i, xi, eta) * basis->detJ(xi, eta);  

			}; 

			D[i][j] = GaussQuad2D(func, int_order); 

		}

	}

	// matrix E (B_i * dB_j/dy)
	MatrixResize(E, N); 
	for (int i=0; i<N; i++) {

		for (int j=0; j<N; j++) {

			auto func = [this, i, j] (double xi, double eta) {

				return basis->dBy(j, xi, eta) * basis->B(i, xi, eta) * basis->detJ(xi, eta); 

			}; 

			E[i][j] = GaussQuad2D(func, int_order); 

		}

	}

}

double Element::evaluate(vector<double> &fin, double xi, double eta) {

	double sum = 0; 

	for (int i=0; i<fin.size(); i++) {

		sum += fin[i] * basis->B(i, xi, eta); 

	}

	return sum; 

}

vector<double> Element::local2global(double xi, double eta) {

	vector<double> out(2); // x, y pair 

	for (int i=0; i<N; i++) {

		out[0] += pts_global[i][0] * basis->B(i, xi, eta); 
		out[1] += pts_global[i][1] * basis->B(i, xi, eta); 

	}

	return out; 

}

double Element::interpolate(vector<double> &fin, vector<double> &out) {

	out = local2global(0, 0); 

	double fout = evaluate(fin, 0, 0); 

	return fout; 

}
