#ifndef __MATRIX_VECTOR_PRODUCT_HH__
#define __MATRIX_VECTOR_PRODUCT_HH__

#include <vector>

using namespace std; 

// CPU

// --- addition --- 
// vector vector 
vector<double> operator+(const vector<double> &A, const vector<double> &B);
// matrix matrix 
vector<vector<double>> operator+(vector<vector<double>> &A, vector<vector<double>> &B);

// --- subtraction ---
// vector vector 
vector<double> operator-(const vector<double> &A, const vector<double> &B); 
// matrix matrix 
vector<vector<double>> operator-(vector<vector<double>> &A, vector<vector<double>> &B);

// --- multiplication ---
// matrix matrix 
vector<vector<double>> operator*(vector<vector<double>> &A, vector<vector<double>> &B); 
// matrix vector 
vector<double> operator*(vector<vector<double>> &A, vector<double> &x);
// vector scalar 
vector<double> operator*(const vector<double> &A, const double x); 
// matrix scalar 
vector<vector<double>> operator*(vector<vector<double>> &A, double x);

// identity 
vector<vector<double>> Identity(int N); 

// GPU
vector<double> MVP_GPU(vector<vector<double>> &A, vector<double> &x);  

#endif 