#include "SparseMatrix.hh"
#include <iostream>
#include <fstream>
#include <cassert>

SparseMatrix::SparseMatrix() {

	m_m = 0; 
	m_n = 0; 

	m_zero = 0.0; 

}

SparseMatrix::SparseMatrix(int a_M, int a_N) {

	// store matrix size 
	m_m = a_M; // rows 
	m_n = a_N; // cols 

	// set zero 
	m_zero = 0.0; 

	// resize the vectors to m x 0 vectors 
	m_colIndex.resize(m_m); 
	m_data.resize(m_m);

}

vector<double> SparseMatrix::operator*(const vector<double>& a_v) const {

	int N = a_v.size(); 

	// assert(N == m_n); 
	if (N != m_n) {
		cout << N << " " << m_n << endl; 
	}

	vector<double> sol(m_m, 0); 

	#pragma omp parallel
	{

		#pragma omp for schedule(static)
		for (int i=0; i<m_m; i++) {

			for (int j=0; j<m_colIndex[i].size(); j++) {

				int col = m_colIndex[i].at(j); 

				sol[i] += m_data[i][j] * a_v.at(col); 

			}

		}

	}

	return sol; 
}

double& SparseMatrix::operator[](array<int, 2>& a_index) {

	// extract the row of indices from m_colIndex 
	vector<int> row = m_colIndex[a_index[0]]; 

	int col = -1; 
	for (int i=0; i<row.size(); i++) {

		if (row[i] == a_index[1]) col = i; 

	}

	// if not found 
	if (col == -1) {

		// make it zero in data 
		m_data[a_index[0]].push_back(m_zero); 

		// add the index into colindex 
		m_colIndex[a_index[0]].push_back(a_index[1]); 

		col = row.size(); 

	}

	return m_data[a_index[0]][col]; 

}

const double& SparseMatrix::operator[](array<int, 2>& a_index) const {

		// extract the row of indices from m_colIndex 
	vector<int> row = m_colIndex[a_index[0]]; 

	int col = -1; 
	for (int i=0; i<row.size(); i++) {

		if (row[i] == a_index[1]) col = i; 

	}

	// if not found 
	if (col == -1) return m_zero; 
	else return m_data[a_index[0]][col]; 

}

void SparseMatrix::zero() {

	for (int i=0; i<m_m; i++) {

		for (int j=0; j<m_data[i].size(); j++) {

			m_data[i][j] = m_zero; 

		}

	}

}

SparseMatrix SparseMatrix::transpose() const {

	// create new sparse matrix 
	SparseMatrix sp(m_m, m_n); 

	for (int i=0; i<m_m; i++) {

		for (int j=0; j<m_colIndex[i].size(); j++) {

			// flip coordinates 
			array<int, 2> p = {m_colIndex[i][j], i}; 

			sp[p] = m_data[i][j]; 

		}

	}

	return sp; 

}

bool SparseMatrix::symmetric() const {

	const SparseMatrix T = this->transpose(); 

	// check this == T 
	for (int i=0; i<m_m; i++) {

		for (int j=0; j<m_colIndex[i].size(); j++) {

			array<int, 2> p = {i, m_colIndex[i][j]}; 

			if (T[p] != m_data[i][j]) return false; 

		}

	}

	return true; 

}

void SparseMatrix::print() const {

	ofstream out; 
	out.open("matrix.csv"); 

	// loop through rows 
	for (int i=0; i<m_m; i++) {

		// build a vector of length m_n 
		vector<double> row(m_n, 0); 

		// fill vector with nonzero elements 
		for (int j=0; j<m_colIndex[i].size(); j++) {

			row[m_colIndex[i][j]] = m_data[i][j]; 

		}

		// print the vector 
		for (int j=0; j<m_n; j++) {

			cout << row[j] << " ";

			out << row[j]; 

			if (j != m_n-1) out << ",";  

		}

		// skip a line 
		cout << endl; 
		out << endl;

	}

	out.close(); 

}

unsigned int SparseMatrix::M() const {

	return m_m; 

}

unsigned int SparseMatrix::N() const {

	return m_n; 

}