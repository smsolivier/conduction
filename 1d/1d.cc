#include <omp.h>

#include "table.hh"
#include "Element.hh"
#include "helper.hh"
#include "Timer.hh"

#include <vector>
#include <iostream>

#define OMP_NUM_THREADS 4

using namespace std; 

int main() {

	Timer timer("Wall Time = "); 
	timer.start(); 

	// number of elements 
	int N = 30; 

	// domain boundary 
	double xb = 1; 

	// time end 
	double tend = .3;

	// number of time steps 
	int Nt = 100; 
	double dt = tend/Nt; // time step size 

	// boundary value 
	double TL = 10;
	double TR = 10;  

	double a = 1; 
	double Q = 0; 
	double alpha = .5; 

	auto kappa = [xb] (double x) {

		if (x < xb/2) return 1.0; 
		else if (x >= xb/2) return 1.0; 

	};

	// time locations 
	vector<double> t = linspace(0, tend, Nt+1); 

	vector<double> x = linspace(0, xb, N+1); 

	// polynomial order 
	int p = 4; 

	// number of nodes 
	int nNodes = (p-1) * N + 1; 

	cout << "nNodes = " << nNodes << endl; 

	// build mesh of elements 
	vector<Element *> el(N);

	// --- start omp parallel ---
	#pragma omp parallel num_threads(OMP_NUM_THREADS)
	{

		#pragma omp for schedule(dynamic, 1) nowait 
		for (int i=0; i<N; i++) {

			vector<double> box = {x[i], x[i+1]}; 

			vector<int> nodes(p); 
			for (int j=0; j<p; j++) {

				nodes[j] = i*(p-1) + j; 

			}
	
			el[i] = new Element(box, nodes); 
	
		} 

	}
	// -- end omp parallel ---

	vector<double> T(nNodes); 

	for (int i=0; i<N; i++) {

		double tval; 

		if (x[i] >= 4*xb/5) tval = TR; 
		else if (x[i] >= 3*xb/5) tval = 0; 
		else if (x[i] >= 2*xb/5) tval = TL; 
		else if (x[i] >= 1*xb/5) tval = 0; 
		else if (x[i] >= 0) tval = TL; 

		for (int j=0; j<p; j++) {

			T[i*(p-1)+j] = tval; 

		}
	}

	for (int l=1; l<Nt+1; l++) {

		// build global system 
		vector<vector<double>> A; 
		MatrixResize(A, nNodes); 

		vector<double> rhs(nNodes); 

		// loop through elements 
		for (int i=0; i<N; i++) {

			Element mEl = *el[i]; // element i 

			// row of local system 
			for (int j=0; j<p; j++) {

				// column of local system 
				for (int k=0; k<p; k++) {

					A[mEl.node[j]][mEl.node[k]] += a/dt*mEl.Z[j][k] + kappa(x[i])*alpha*mEl.X[j][k]; 

					rhs[mEl.node[j]] += a/dt*mEl.Z[j][k]*T[mEl.node[j]] - 
						kappa(x[i])*(1-alpha)*T[mEl.node[j]]*mEl.X[j][k] + Q*mEl.Z[j][k];

				}

			}

		}

		// boundary conditions 

		// remove first equation 
		for (int i=0; i<nNodes; i++) {

			A[0][i] = 0; 

		}

		// subtract first column 
		for (int i=0; i<nNodes; i++) {

			rhs[i] -= TL * A[i][0]; 
			A[i][0] = 0; // remove first column 

		}

		A[0][0] = 1; // set first equation to boundary value 
		rhs[0] = TL; // set boundary value 

		// remove last equation 
		for (int i=0; i<nNodes; i++) {

			A[nNodes-1][i] = 0; 

		}

		// subtract first column 
		for (int i=0; i<nNodes; i++) {

			rhs[i] -= TR * A[i][nNodes-1]; 
			A[i][nNodes-1] = 0; // remove first column 

		}

		A[nNodes-1][nNodes-1] = 1; // set first equation to boundary value 
		rhs[nNodes-1] = TR; // set boundary value 

		// solve global system 
		int status = gauss_elim(A.size(), A, T, rhs); 

		if (status == -1) cout << "linear solver failed" << endl; 

		vector<double> T_int; 
		vector<double> x_int; 

		for (int i=0; i<N; i++) {

			vector<double> fin(p); 
			for (int j=0; j<p; j++) {

				fin[j] = T[i*(p-1)+j]; 

			}

			vector<double> xout; 

			vector<double> fout = el[i]->interpolate(fin, xout); 

			for (int j=0; j<p-1; j++) {

				T_int.push_back(fout[j]); 
				x_int.push_back(xout[j]); 

			}

		}

		// write solution to file 
		Table table("data/out_"+to_string(l-1)); 

		table.addColumn(x_int); 
		table.addColumn(T_int); 

		table.write(); 

	}

	timer.stop(); 

}