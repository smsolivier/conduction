#ifndef LinearSolver_hh
#define LinearSolver_hh

#include <vector>

using namespace std; 

int gauss_elim(int dim, vector<vector<double>> &a, vector<double> &x, vector<double> &b); 

#endif 