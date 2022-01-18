#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "lib.h"
#include <string>

using namespace std;

double modulo_quadro(int dim, double *x){
	double appo=0;
	for(int i=0;i<dim; i++){
		appo+=x[i]*x[i];
	}
	return appo;
}

double mean(int N, double * vett){
	double appo=0;
	for(int i=0;i<N;i++){
		appo+=vett[i];
	}
return appo/N;
}



