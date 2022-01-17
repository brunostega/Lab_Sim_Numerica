
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "lib.h"
#include "/home/bruno/Scrivania/universita/magistrale1/lab_simulaz/esercizi/random.h"
#include <string>

using namespace std;

double function(double x){         //function to integrate
	return M_PI/2*cos(x*M_PI/2);
}

double retta(double x){            //function for the importance sampling
	return 2*(1-x);
}

double cumulative_retta(double x){ //function used to sample "retta" with the inverse of the cumulative function
	return 1-sqrt(1-x);
}




