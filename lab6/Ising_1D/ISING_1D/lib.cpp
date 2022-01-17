#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "lib.h"
#include <string>

using namespace std;



//BOLTZMANN PROBABILITY//////////////////////////////////////
double Boltzmann_prob::eval(double DeltaE){
	return exp(-beta*DeltaE)
}


//METROPOLIS//////////////////////////////////

acceptance::acceptance(){
	Random rnd;
	rnd.Initialize();
	_rnd=rnd;
}
acceptance::~acceptance(){}



Metropolis::Metropolis(int passi, probability *prob, acceptance *T){   //initializing 
	_passi=passi;
	_rnd.Initialize();
	_probability= prob;
	_acceptance=T;
}

Metropolis::~Metropolis(){}




posizione Metropolis::Acceptance(posizione x_old, posizione x_new){
	double r=_rnd.Rannyu();
	if(1<_probability->eval(x_new)/_probability->eval(x_old)){
		return x_new;
	}
	else{
		if(_probability->eval(x_new)/_probability->eval(x_old)<r){
			return x_old;
		}
		else{
			return x_new;
		}
	}
}

posizione Metropolis::MR_alghoritm(posizione x_start){
	posizione x_new(0,0,0);  
	for(int i=0;i<_passi;i++){
		x_new=this->Acceptance(x_start,_acceptance->eval(x_start));
		x_start=x_new;
	}
	return x_new;
}
