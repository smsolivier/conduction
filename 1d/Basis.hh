#ifndef __BASIS_HH__
#define __BASIS_HH__

#include <vector> 

using namespace std; 

class Basis {

public: 

	Basis(vector<double> &xloc); 
	double B(int i, double xi); 
	double dB(int i, double xi); 

private:

	vector<vector<double>> coef;

	int p; 
};

#endif 