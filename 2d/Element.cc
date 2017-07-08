#include "Element.hh"

#include "GaussQuad.hh"
#include "helper.hh"

#include <cmath>
#include <iostream>

Element::Element(vector<double> &box, vector<int> &neighbors, int p) :
	neighbors(neighbors), p(p) {

	N = pow(p+1, 2); // number of nodes 

	globalNodes.resize(N, -1); 

	// int int_order = (p+1)/2 + 1; 
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

void Element::setFace(int ind, vector<int> &nodes) {

	// set face 0 
	if (ind == 0) {

		for (int i=0; i<nodes.size(); i++) {

			int index = i; 

			if (globalNodes[index] == -1) {

				globalNodes[index] = nodes[i]; 

			}

		}

	}

	// set face 1 
	else if (ind == 1) {

		for (int i=0; i<p+1; i++) {

			int index = i*(p+1); 

			if (globalNodes[index] == -1) {

				globalNodes[index] = nodes[i]; 

			}

		}

	}

	// set face 2 
	else if (ind == 2) {

		for (int i=0; i<p+1; i++) {

			int index = p*(p+1)+i;

			if (globalNodes[index] == -1) {

				globalNodes[index] = nodes[i]; 

			}

		}

	}

	// set face 3 
	else if (ind == 3) {

		for (int i=0; i<p+1; i++) {

			int index = p+i*(p+1);

			if (globalNodes[index] == -1) {

				globalNodes[index] = nodes[i]; 

			}

		}

	}

	else cout << "face not defined" << endl; 

}

void Element::setFace(int ind, int &count) {

	// set face 0 
	if (ind == 0) {

		for (int i=0; i<p+1; i++) {

			int index = i; 

			if (globalNodes[index] == -1) {

				globalNodes[index] = count;
				count++; 

			}

		}

	}

	// set face 1 
	else if (ind == 1) {

		for (int i=0; i<p+1; i++) {

			int index = i*(p+1); 

			if (globalNodes[index] == -1) {

				globalNodes[index] = count;
				count++; 

			}

		}

	}

	// set face 2 
	else if (ind == 2) {

		for (int i=0; i<p+1; i++) {

			int index = p*(p+1)+i;

			if (globalNodes[index] == -1) {

				globalNodes[index] = count;
				count++; 

			}

		}

	}

	// set face 3 
	else if (ind == 3) {

		for (int i=0; i<p+1; i++) {

			int index = p+i*(p+1);

			if (globalNodes[index] == -1) {

				globalNodes[index] = count;
				count++; 

			}

		}

	}

	else cout << "face not defined" << endl; 

}

vector<int> Element::getFace(int ind) {

	vector<int> mNode(p+1); 

	// get face 0 
	if (ind == 0) {

		for (int i=0; i<p+1; i++) {

			mNode[i] = globalNodes[i]; 

		}

	}

	// get face 1 
	else if (ind == 1) {

		for (int i=0; i<p+1; i++) {

			int index = i*(p+1);
			mNode[i] = globalNodes[index]; 

		}

	}

	// get face 2 
	else if (ind == 2) {

		for (int i=0; i<p+1; i++) {

			int index = p*(p+1)+i; 
			mNode[i] = globalNodes[index]; 

		}

	}

	// get face 3 
	else if (ind == 3) {

		for (int i=0; i<p+1; i++) {

			int index = p+i*(p+1); 
			mNode[i] = globalNodes[index]; 

		}

	}

	else cout << "face not defined" << endl; 

	return mNode; 

}

void Element::fillMiddle(int &count) {

	for (int i=0; i<N; i++) { 

		if (globalNodes[i] == -1) {

			globalNodes[i] = count; 

			count++; 

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

vector<vector<double>> Element::interpolate(vector<double> &fin, vector<vector<double>> &xout,
	vector<vector<double>> &yout) {

	vector<double> xmid(p), ymid(p); 
	vector<vector<double>> fout; 
	MatrixResize(fout, p);
	MatrixResize(xout, p);
	MatrixResize(yout, p);

	for (int i=0; i<p; i++) {

		xmid[i] = (xloc[i+1] + xloc[i])/2; 
		ymid[i] = (yloc[i+1] + yloc[i])/2; 

	}

	for (int i=0; i<p; i++) {

		for (int j=0; j<p; j++) {

			vector<double> glob = local2global(xmid[i], ymid[j]); 

			fout[j][i] = evaluate(fin, xmid[i], ymid[j]); 
			xout[j][i] = glob[0]; 
			yout[j][i] = glob[1]; 

		}

	}

	return fout; 

}

double Element::interpolateSingle(vector<double> &fin, double &xout, double &yout) {

	double val = evaluate(fin, 0, 0); 

	vector<double> glob = local2global(0, 0); 

	xout = glob[0]; 

	yout = glob[1]; 

	return val; 
	
}