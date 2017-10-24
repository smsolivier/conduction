#include <omp.h>
#include "VectorMath.hh"
#include <iostream>

// #define CUDA

// vector vector addition 
vector<double> operator+(const vector<double> &A, const vector<double> &B) {

	vector<double> C(A.size()); 

	#pragma omp parallel 
	{

		int i; 
		#pragma omp for private(i) schedule(dynamic, 1) nowait 
		for (i=0; i<A.size(); i++) {

			C[i] = A[i] + B[i]; 

		}

	}

	return C; 

}

// matrix matrix addition 
vector<vector<double>> operator+(vector<vector<double>> &A, vector<vector<double>> &B) {

	vector<vector<double>> C(A.size(), vector<double>(A.size())); 

	#pragma omp parallel 
	{

		#pragma omp for schedule(dynamic, 1) nowait 
		for (int i=0; i<A.size(); i++) {

			for (int j=0; j<A.size(); j++) {

				C[i][j] = A[i][j] + B[i][j]; 

			}

		}

	}

	return C; 

} 

// matrix matrix product 
vector<vector<double>> operator*(vector<vector<double>> &A, vector<vector<double>> &B) {

	int N = A.size(); 

	vector<vector<double>> C(N, vector<double>(N)); 

	#pragma omp parallel 
	{

		#pragma omp for schedule(dynamic, 1) nowait
		for (int i=0; i<N; i++) {

			for (int j=0; j<N; j++) {

				for (int k=0; k<N; k++) {

					C[i][j] += A[i][k] * B[k][j]; 

				}

			}

		}

	}

	return C; 
}

// matrix vector product 
vector<double> operator*(vector<vector<double>> &A, vector<double> &x) {

	int N = x.size(); 

	vector<double> b(N); 

	// --- OMP parallel --- 
	#ifndef CUDA

	#pragma omp parallel 
	{

		#pragma omp for schedule(dynamic, 1) nowait 
		for (int i=0; i<N; i++) {

			for (int j=0; j<N; j++) {

				b[i] += A[i][j] * x[j]; 

			}

		}

	}

	#endif 
	// --- end OMP parallel --- 


	// --- CUDA parallel --- 
	#ifdef CUDA

	b = MVP_GPU(A, x); 

	#endif 
	// --- end CUDA parallel --- 


	return b; 

}

// vector scalar product 
vector<double> operator*(const vector<double> &A, const double x) {

	int N = A.size(); 

	vector<double> b(N); 

	#pragma omp parallel 
	{

		#pragma omp for schedule(dynamic, 1) nowait 
		for (int i=0; i<N; i++) {

			b[i] = A[i] * x; 

		}

	}

	return b; 

}

// matrix scalar product 
vector<vector<double>> operator*(vector<vector<double>> &A, double x) {

	int N = A.size(); 

	vector<vector<double>> b(N, vector<double>(N)); 

	#pragma omp parallel 
	{

		#pragma omp for schedule(dynamic, 1) nowait 
		for (int i=0; i<N; i++) {

			for (int j=0; j<N; j++) {

				b[i][j] = A[i][j] * x; 

			}

		}

	}

	return b; 

}

// vector vector subtraction 
vector<double> operator-(const vector<double> &A, const vector<double> &B) {

	vector<double> C(A.size()); 

	#pragma omp parallel 
	{

		int i; 
		#pragma omp for private(i) schedule(dynamic, 1) nowait 
		for (i=0; i<A.size(); i++) {

			C[i] = A[i] - B[i]; 

		}

	}

	return C; 

}

// matrix matrix subtraction  
vector<vector<double>> operator-(vector<vector<double>> &A, vector<vector<double>> &B) {

	vector<vector<double>> C(A.size(), vector<double>(A.size())); 

	#pragma omp parallel 
	{

		#pragma omp for schedule(dynamic, 1) nowait 
		for (int i=0; i<A.size(); i++) {

			for (int j=0; j<A.size(); j++) {

				C[i][j] = A[i][j] - B[i][j]; 

			}

		}

	}

	return C; 

} 

// return the identity matrix 
vector<vector<double>> Identity(int N) {


	vector<vector<double>> I(N, vector<double>(N)); 

	for (int i=0; i<N; i++) {

		I[i][i] = 1; 

	}

	return I; 
}