#ifndef __BASIS_HH__
#define __BASIS_HH__

#include <vector> 

using namespace std; 

class Basis {

public: 

	Basis(vector<vector<double>> &pts); 
	double B(int i, double xi, double eta); 
	double dBx(int i, double xi, double eta); 
	double dBy(int i, double xi, double eta); 
	double Jacobian(int i, int j, double xi, double eta); 
	double detJ(double xi, double eta); 
	double J_in(int i, int j, double xi, double eta); 

private:

	double dBxi(int i, double xi, double eta); 
	double dBeta(int i, double xi, double eta); 

	vector<vector<double>> coef, box; 

	vector<vector<double>> pts; 

	int p; 
};

#endif 