#ifndef WRITECURVE_HH
#define WRITECURVE_HH

#include <string> 
#include <vector>

using namespace std; 

class WriteCurve {

public:

	string dir; // name of ouput directory 
	string name; // name of curve 
	WriteCurve(string dir, string name); 

	void print(vector<double> &x, vector<double> &f); 

	int fcount; // store how many files written 

};

#endif 