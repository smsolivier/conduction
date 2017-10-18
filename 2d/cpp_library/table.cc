#include "table.hh"

#include <fstream>

Table::Table(string name) : name(name) {

	nCol = 0;
	nRow = 0; 

}

void Table::addColumn(vector<double> &x) {

	nCol++; 

	nRow = x.size(); 

	table.push_back(x); 

}

void Table::write() {

	ofstream file; 
	file.open(name); 

	for (int i=0; i<nRow; i++) {

		for (int j=0; j<nCol; j++) {

			file << table[j][i] << " "; 

		}

		file << endl; 

	}

	file.close(); 

}