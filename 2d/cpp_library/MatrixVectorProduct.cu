#include "VectorMath.hh"

#include "helper.hh"
#include <vector>

#include <iostream>

using namespace std; 

#define THREADS_PER_BLOCK 1024

__global__ void MVP_kernel(double * A, double * x, double * b, int N) {

	int tid = blockIdx.x * blockDim.x + threadIdx.x; 

	double sum = 0; 

	for (int i=0; i<N; i++) {

		int ind = tid*N + i; 

		sum += A[ind] * x[i]; 

	}

	b[tid] = sum; 

}

vector<double> MVP_GPU(vector<vector<double>> &A, vector<double> &x) {

	int N = A.size(); 

	size_t dSize = sizeof(double); 

	// device pointers 
	double *d_A, *d_x, *d_b;
	double *b; 

	b = new double[N]; 

	cudaMalloc(&d_A, N*N*dSize); 
	cudaMalloc(&d_x, N*dSize); 
	cudaMalloc(&d_b, N*dSize); 

	double * array = new double[N*N]; 

	// convert from vector to array 
	for (int i=0; i<N; i++) {

		for (int j=0; j<N; j++) {

			int ind = N*i + j; 

			array[ind] = A[i][j]; 

		}

	}

	cudaMemcpy(d_A, array, N*N*dSize, cudaMemcpyHostToDevice); 
	cudaMemcpy(d_x, &x[0], N*dSize, cudaMemcpyHostToDevice); 

	MVP_kernel<<<N/THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(d_A, d_x, d_b, N); 

	cudaMemcpy(b, d_b, N*dSize, cudaMemcpyDeviceToHost);

	// convert to vector 
	vector<double> b_vec(N); 
	for (int i=0; i<N; i++) {

		b_vec[i] = b[i]; 

	}

	cudaFree(d_A); cudaFree(d_x); cudaFree(d_b); 

	delete b; delete array; 

	return b_vec; 

}