#include "Element.hh" 

#include "GaussQuad.hh" 
#include "helper.hh"

#include <cmath>
#include <iostream>

Element::Element(vector<double> &box, vector<int> &nodes) : 
	node(nodes) {

	p = node.size(); 

	xglob = linspace(box[0], box[1], p); // generate global points, p evenly spaced points 

	xloc = linspace(-1, 1, p); // local points, p evenly spaced points 

	basis = new Basis(xloc); 

	// generate X matrix 
	MatrixResize(X, p); 
	for (int i=0; i<p; i++) {

		for (int j=0; j<p; j++) {

			// lambda function to integrate 
			auto func = [this, i, j] (double xi) {

				// return dB[i].evaluate(xi) * dB[j].evaluate(xi) * 1/Jacobian(xi); 
				return basis->dB(i, xi) * basis->dB(j, xi) * 1/Jacobian(xi); 

			};

			// pass function to integrate with gauss quadrature 
			X[i][j] = GaussQuad(func, p+1); 

		}

	}

	// generate Y matrix 
	// Y.resize(p);
	// for (int i=0; i<p; i++) {

	// 	Y[i].resize(p); 		

	// 	for (int j=0; j<p; j++) {

	// 		auto func = [this, i, j] (double xi) {

	// 			return B[i].evaluate(xi) * dB[j].evaluate(xi); 

	// 		};

	// 		Y[i][j] = GaussQuad(func, p+1); 

	// 	}
	// }

	// generate Z matrix 
	MatrixResize(Z, p); 
	for (int i=0; i<p; i++) {

		for (int j=0; j<p; j++) {

			auto func = [this, i, j] (double xi) {

				// return B[i].evaluate(xi) * B[j].evaluate(xi) * Jacobian(xi); 
				return basis->B(i, xi) * basis->B(j, xi) * Jacobian(xi); 

			};

			Z[i][j] = GaussQuad(func, p+1); 

		}
	}

}

double Element::Jacobian(double xi) {
	/* Compute Jacobian (dx/dxi) */ 

	double sum = 0; 
	for (int i=0; i<p; i++) {

		// sum += xglob[i] * dB[i].evaluate(xi); 
		sum += xglob[i] * basis->dB(i, xi); 

	}

	return sum; 

}

double Element::xi2x(double xi) {

	double sum = 0; 

	for (int i=0; i<p; i++) {

		// sum += xglob[i] * B[i].evaluate(xi); 
		sum += xglob[i] * basis->B(i, xi); 

	}

	return sum; 

}

double Element::evaluate(vector<double> &fin, double xi) {

	double sum = 0; 
	for (int i=0; i<p; i++) {

		// sum += fin[i] * B[i].evaluate(xi); 
		sum += fin[i] * basis->B(i, xi); 

	}

	return sum; 

}

vector<double> Element::interpolate(vector<double> &fin, vector<double> &xout) {

	vector<double> fout(p-1); 
	xout.resize(p-1); 

	for (int i=0; i<p-1; i++) {

		double mid = (xloc[i+1] + xloc[i])/2; 

		fout[i] = evaluate(fin, mid); 

		xout[i] = xi2x(mid); 

	}

	return fout; 

}
// int main() {

// 	Element el(0, 1, 4); 

// 	printVector(el.D); 

// 	cout << el.xi2x(-1) << endl; 

// }

