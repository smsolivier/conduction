#ifndef __ELEMENT_HH__
#define __ELEMENT_HH__

#include "Basis.hh"
#include <vector>
#include "Polynomial.hh"
#include "LinearSolver.hh"

using namespace std; 

class Element {

public:

	Element(vector<double> &box, vector<int> &globalNodes, int p); 
	/*  box = global positions of quad 
		globalNodes = global node numbers 
		p = polynomial order 
	*/

	double interpolate(vector<double> &fin, vector<double> &xout); 

	vector<int> globalNodes; 

	// basis combination matrices 
	vector<vector<double>> A, B, C, D, E; 

	vector<vector<double>> pts; // local (xi, eta) points 
	vector<vector<double>> pts_global; // global (xi, eta) points 

	vector<double> xglob; // global x locations
	vector<double> yglob; // global y locations 

private:

	double evaluate(vector<double> &fin, double xi, double eta); 
	vector<double> local2global(double xi, double eta); 

	int p; // polynomial order 
	int N; // number of nodes 

	Basis * basis; 

	vector<double> xloc; // local x locations 
	vector<double> yloc; // local y locations

}; 

#endif 