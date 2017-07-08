#ifndef __ELEMENT_HH__
#define __ELEMENT_HH__

#include "Basis.hh"
#include <vector>
#include "Polynomial.hh"
#include "LinearSolver.hh"

using namespace std; 

class Element {

public:

	Element(vector<double> &box, vector<int> &neighbors, int p); 
	/*  box = global positions of quad 
		globalNodes = global node numbers 
		p = polynomial order 
	*/

	vector<vector<double>> interpolate(vector<double> &fin, vector<vector<double>> &xout, 
		vector<vector<double>> &yout); 
	double interpolateSingle(vector<double> &fin, double &xout, double &yout); 
	void setFace(int i, vector<int> &nodes); 
	void setFace(int ind, int &count); 
	vector<int> getFace(int ind); 
	void fillMiddle(int &count); 

	vector<int> globalNodes; 

	// basis combination matrices 
	vector<vector<double>> A, B, C, D, E; 

	vector<vector<double>> pts; // local (xi, eta) points 
	vector<vector<double>> pts_global; // global (xi, eta) points 

	vector<double> xglob; // global x locations
	vector<double> yglob; // global y locations 

	vector<int> neighbors; 

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