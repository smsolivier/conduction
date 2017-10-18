#ifndef __TIMER_HH__
#define __TIMER_HH__

#include <ctime> 
#include <string>

using namespace std; 

class Timer {

public:

	Timer(string message="Wall Time = ");

	void start(); 
	double stop(); 

private:

	string message; 

	double begin; 
	double end; 

}; 

#endif 