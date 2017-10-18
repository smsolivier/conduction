#include "WriteCurve.hh"

#include <fstream>

WriteCurve::WriteCurve(string dir, string name) : dir(dir), name(name) {

	fcount = 0; // initialize to zero 

}

void WriteCurve::print(vector<double> &x, vector<double> &f) {

	ofstream file; 

	file.open(dir + name + to_string(fcount) + ".curve"); 

	file << "# " << name << endl; 

	for (int i=0; i<x.size(); i++) {

		file << x[i] << " " << f[i] << endl; 

	}

	file.close(); 

	fcount += 1; 

}