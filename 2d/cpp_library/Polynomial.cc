#include <iostream>
#include <cmath> 

#include "Polynomial.hh"

using namespace std; 

Polynomial::Polynomial(vector<double> ccoef) {

	coef = ccoef; 

	p = coef.size(); // number of coefficients 
}

double Polynomial::evaluate(double x) {

	double sum = 0; 

	for (int i=0; i<p; i++) {

		sum += coef[i]*pow(x, i); 
	}

	return sum; 
}

vector<double> Polynomial::evaluate(vector<double> x) {

	int N = x.size(); // number of points to evaluate 

	vector<double> sum(N); 

	for (int i=0; i<N; i++) {

		double isum = evaluate(x[i]); 

		sum[i] = isum; 

	}

	return sum; 
}

Polynomial Polynomial::derivative() {

	vector<double> dcoef; // store derivatives coefficients 

	for (int i=1; i<p; i++) {

		dcoef.push_back(coef[i]*i); 
	}

	Polynomial deriv(dcoef); 

	return deriv; 
}

void Polynomial::print() {

	for (int i=0; i<p; i++) {

		cout << coef[i] << " x^" << i << " "; 

	}

	cout << endl; 
}

// int main() {

// 	Polynomial p(vector<double>{1, 1, 1}); 

// 	Polynomial deriv = p.derivative(); 

// 	cout << deriv.evaluate(1) << endl; 

// }