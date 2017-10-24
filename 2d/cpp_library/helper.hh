#ifndef helper_hh
#define helper_hh

#include <vector>
#include <string> 

using namespace std; 

vector<double> linspace(double start, double end, int N); // generate evenly spaced points 

void printVector(vector<double> x); // print 1D vector 
void printVector(vector<vector<double>> x); // print 2D vector 
void printVector(vector<int> x); 
void printVector(vector<vector<int>> x); 

void writeTable(string name, vector<double> &x, vector<double> &y); 

void MatrixResize(vector<vector<double>> &A, int N); // resize a matrix 
void MatrixResize(vector<vector<double>> &A, int rows, int cols); 

#endif 