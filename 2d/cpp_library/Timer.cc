#include "Timer.hh"

#include <omp.h>

#include <iostream>

Timer::Timer(string message) : message(message) {} 

void Timer::start() {

	begin = omp_get_wtime(); 

}

double Timer::stop() {

	end = omp_get_wtime(); 

	double wtime = double(end - begin); 

	double elapsed = wtime; 

	string units = "seconds"; 

	if (elapsed > 3600) {

		elapsed /= 3600; 
		units = "hours"; 

	}

	else if (elapsed > 60) {

		elapsed /= 60; 
		units = "minutes"; 

	}

	cout << message << elapsed << " " << units << endl; 

	return wtime; 

}