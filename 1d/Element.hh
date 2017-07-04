#ifndef __ELEMENT_HH__
#define __ELEMENT_HH__

#include <vector>
#include "Basis.hh"
#include "LinearSolver.hh"

using namespace std; 

class Element {

public:

	Element(vector<double> &box, vector<int> &nodes); 

	vector<double> interpolate(vector<double> &fin, vector<double> &xout); // interpolate answer

	// basis combination matrices 
	vector<vector<double>> X, Y, Z; 

	vector<int> node; 

private:

	double Jacobian(double xi); // computes jacobian at xi 
	double xi2x(double xi); // convert local to global coordinates 
	double evaluate(vector<double> &fin, double xi); // evaluate answer at xi 

	vector<double> xglob; // global x values within element 
	vector<double> xloc; // local x values within element 

	int p; // polynomial order 

	Basis * basis; 

};

#endif 