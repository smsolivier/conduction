#ifndef __TABLE_HH__
#define __TABLE_HH__

#include <string> 
#include <vector>

using namespace std; 

class Table {

public: 

	Table(string name); 

	void addColumn(vector<double> &x); 

	void write(); 

private:

	int nCol; // store number of columns 
	int nRow; // store number of rows 

	string name; // output file name 

	vector<vector<double>> table; 

};

#endif 