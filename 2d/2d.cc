#include "Element.hh"
#include "helper.hh"
#include "Timer.hh"

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std; 

#define OMP_NUM_THREADS 4

void applyBC(vector<vector<double>> &A, vector<double> &rhs, vector<int> &loc, double val) {

	for (int i=0; i<loc.size(); i++) {

		int ind = loc[i]; 

		// remove equation 
		for (int j=0; j<A.size(); j++) {

			A[ind][j] = 0; 

		}

		// subtract first column 
		for (int j=0; j<A.size(); j++) {

			rhs[j] -= val*A[j][ind];
			A[j][ind] = 0; // remove column 

		}

		A[ind][ind] = 1; 
		rhs[ind] = val; 

	}

}

int main() {

	Timer timer; 
	timer.start(); 
	// --- start program --- 


	const int Nx = 20; 
	const int Ny = 20; 
	const int N = Nx * Ny; // number of elements 
	const int Nt = 150; // number of time steps 
	const double T0 = 0; 
	const double kappax = 1;
	const double kappay = 1;  

	const double a = 1; 
	const double b = 1; 
	const double c = 1; 
	const double alpha = 1; 

	int p = 2; 

	int nNodes = (Ny+1)*(Nx+1); 

	cout << "nNodes = " << nNodes << endl; 

	double xb = .1; 
	double yb = .1; 
	double tend = .005; 

	vector<double> t = linspace(0, tend, Nt+1); 

	auto Q = [kappax, kappay, a, b, c, xb, yb] (double x, double y) {

		// conduction MMS
		// return kappax*pow(M_PI/xb, 2)*sin(M_PI*x/xb)*sin(M_PI*y/yb) + 
		// 	kappay*pow(M_PI/yb, 2)*sin(M_PI*x/xb)*sin(M_PI*y/yb); 

		return 0; 
	};

	auto kappa = [xb, yb] (double x, double y) {

		double val = 1.0; 
		// if (x > 5*xb/10 && x < 8*xb/10 && y > 5*yb/10 && y < 8*yb/10) val = 10000; 

		return val; 
	};

	vector<double> x = linspace(0, xb, Nx+1);
	vector<double> y = linspace(0, yb, Ny+1);  

	vector<Element *> el(N); 

	vector<int> bL, bR, bB, bT; 

	Timer meshTime("Mesh Time = "); 
	meshTime.start(); 

	// build mesh of elements 
	// --- start omp parallel --- 
	#pragma omp parallel num_threads(OMP_NUM_THREADS)
	{

		#pragma omp for schedule(dynamic, 1) nowait
		for (int i=0; i<Nx; i++ ) {

			for (int j=0; j<Ny; j++) {

				int elNum = i*Ny+j; 

				// vertices of quad 
				vector<double> box = {x[i], x[i+1], y[j], y[j+1]}; 

				// global node numbers 
				vector<int> nodes = {(Ny+1)*i+j, (Ny+1)*(i+1)+j, (Ny+1)*i+j+1, (Ny+1)*(i+1)+j+1}; 

				el[elNum] = new Element(box, nodes, p); 

			}

		}

	}
	// --- end omp parallel --- 
	meshTime.stop(); 

	// determine boundary 

	for (int i=0; i<Nx+1; i++ ) {

		for (int j=0; j<Ny+1; j++) {

			int elNum = i*(Ny+1) + j; 

			if (i == 0) bL.push_back(elNum); 
			if (j == 0) bB.push_back(elNum); 
			if (i == Nx) bR.push_back(elNum); 
			if (j == Ny) bT.push_back(elNum); 

		}

	}

	// store solution 
	vector<double> T(nNodes); 

	// set initial conditions 
	// for (int i=0; i<N; i++) {

	// 	Element mEl = *el[i]; 

	// 	for (int j=0; j<4; j++) {

	// 		if (mEl.xglob[j] > xb/3 && mEl.xglob[j] < 2*xb/3 && 
	// 			mEl.yglob[j] > yb/3 && mEl.yglob[j] < 2*yb/3) {

	// 			T[mEl.globalNodes[j]] = 1; 

	// 		}

	// 	}

	// }

	// loop through time 
	for (int l=1; l<t.size(); l++) {

		cout << double(l)/t.size() << " \r";   
		cout.flush(); 

		double dt = t[l] - t[l-1]; 

		// global system 
		vector<vector<double>> A; 
		MatrixResize(A, nNodes); 
		vector<double> rhs(nNodes); 

		// loop through elements 
		for (int i=0; i<N; i++) {

			Element mEl = *el[i]; // element i 

			// row of local system 
			for (int j=0; j<4; j++) {

				// column of local system 
				for (int k=0; k<4; k++) {

					double tc = kappa(mEl.pts_global[k][0], mEl.pts_global[k][1]); 

					int row = mEl.globalNodes[j]; 
					int col = mEl.globalNodes[k]; 

					A[row][col] += a/dt*mEl.A[j][k]; 
					A[row][col] += alpha*tc*mEl.B[j][k]; 
					A[row][col] += alpha*tc*mEl.C[j][k]; 

					// A[row][col] += a*mEl.Z[j][k];
					// A[row][col] += b*mEl.W[j][k]; 
					// A[row][col] += c*mEl.X[j][k]; 

					rhs[row] += (a/dt*mEl.A[j][k] - 
						(1-alpha)*(tc*mEl.B[j][k] + tc*mEl.C[j][k]))*T[row]; 
					rhs[row] += Q(mEl.pts_global[k][0], mEl.pts_global[k][1])*mEl.A[j][k]; 

				}

			}

		}

		// boundary conditions 
		applyBC(A, rhs, bL, T0); 
		// applyBC(A, rhs, bT, T0); 
		applyBC(A, rhs, bB, T0); 
		// applyBC(A, rhs, bR, T0); 

		for (int i=0; i<bT.size(); i++) {

			int ind = bT[i]; 

			rhs[ind] += 50; 

		}

		for (int i=0; i<bR.size(); i++) {

			int ind = bR[i]; 

			rhs[ind] += 50; 

		}

		// solve system 
		int status = gauss_elim(A.size(), A, T, rhs); 

		if (status != 0) cout << "linear solver problem" << endl; 

		// print solution to file 
		vector<vector<double>> sol, X, Y; 
		MatrixResize(sol, Ny, Nx); 
		MatrixResize(X, Ny, Nx); 
		MatrixResize(Y, Ny, Nx); 

		for (int i=0; i<Nx; i++) {

			for (int j=0; j<Ny; j++) {

				Element mEl = *el[i*Ny+j]; 

				vector<double> fin(4); 
				for (int k=0; k<4; k++) {

					fin[k] = T[mEl.globalNodes[k]]; 

				}

				vector<double> out; 
				double val = mEl.interpolate(fin, out); 

				sol[j][i] = val; 

				X[j][i] = out[0]; 
				Y[j][i] = out[1]; 

			}

		}

		ofstream file; 
		file.open("data/out"+to_string(l-1)); 

		ofstream gridX;
		ofstream gridY;  
		gridX.open("data/X"+to_string(l-1)); 
		gridY.open("data/Y"+to_string(l-1)); 
		for (int i=0; i<Nx; i++) {

			for (int j=0; j<Ny; j++) {

				file << sol[j][i] << " "; 
				gridX << X[j][i] << " "; 
				gridY << Y[j][i] << " "; 

			}

			file << endl; 
			gridX << endl; 
			gridY << endl; 

		}

		file.close(); 
		gridX.close(); 
		gridY.close(); 

	}


	// --- end program --- 
	timer.stop(); 

}