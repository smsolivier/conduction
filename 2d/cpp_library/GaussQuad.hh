#include <vector>
#include <math.h>
#include <iostream>

using namespace std; 

// generate weights, abcisae 
vector<vector<double>> createRef1dQuad(double x1, double x2, int order);

double GaussQuad(auto &f, int p) {
	/* do gauss integration using weights, absicae from createRef1dQuad */ 

	vector<vector<double>> quad = createRef1dQuad(-1, 1, p);  

	double sum = 0; 

	for (int i=0; i<quad.size(); i++) {

		sum += f(quad[i][0])*quad[i][1]; 

	}

	return sum; 

}

double GaussQuad2D(auto &f, int p) {

	vector<vector<double>> quad = createRef1dQuad(-1, 1, p); 

	double sum = 0; 

	for (int i=0; i<quad.size(); i++) {

		for (int j=0; j<quad.size(); j++) {

			sum += f(quad[i][0], quad[j][0]) * quad[i][1] * quad[j][1]; 

		}

	}

	return sum; 

}

// double GaussQuad2D(auto &f) {

// 	double diff = 1000; 
// 	double tol = 1e-8; 

// 	double int1 = 0; 
// 	double int2 = 0; 

// 	int p = 1; 

// 	while (diff > tol) {

// 		int1 = GaussQuad2D(f, p); 

// 		while (isnan(int1)) {

// 			p++; 

// 			int1 = GaussQuad2D(f, p); 

// 		}

// 		p++; 
// 		int2 = GaussQuad2D(f, p); 

// 		while (isnan(int2)) {

// 			p++; 

// 			int2 = GaussQuad2D(f, p); 

// 		}

// 		diff = abs(int1 - int2)/int2; 

// 	}

// 	cout << "p used = " << p << endl; 

// 	return int2; 

// }