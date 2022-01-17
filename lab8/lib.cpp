#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "lib.h"
#include <string>

using namespace std;

double a_0=1;      //valore effettivo=0.0529;

//POSIZIONE////////////////////////////////

posizione::posizione(){
	m_x=0;
}
posizione::posizione(double x){
	m_x=x;
}

posizione::~posizione(){}

double posizione::getX() const{
return m_x;
}



double posizione::getDistanza( const posizione & b) const{
return (sqrt(pow(getX()-b.getX(),2)));
}




//WAVE FUNCTION SQUARED//////////////////////////////////////
probability::probability(){}
probability::~probability(){}

Wave_func::Wave_func(double mu, double sigma){
	_sigma=sigma;
	_mu=mu;
}

void Wave_func::reset(double mu, double sigma){
	_sigma=sigma;
	_mu=mu;
}
double Wave_func::getMu(){return _mu;}
double Wave_func::getSigma(){return _sigma;}

Wave_func::~Wave_func(){}

double Wave_func::eval(posizione x){
	double exp1=x.getX()+_mu;
	double exp2=x.getX()-_mu;
	return pow(   exp(-pow(exp2,2)/(_sigma*_sigma*2)) + exp( -pow(exp1,2)/(_sigma*_sigma*2))  ,2);
	//return exp(-pow(x.getX(),2)/(_sigma*_sigma*2));
}




//METROPOLIS//////////////////////////////////

acceptance::acceptance(){
	Random rnd;
	rnd.Initialize();
	_rnd=rnd;
}
acceptance::~acceptance(){}

T_uniform::T_uniform(){}
T_uniform::T_uniform(double delta){
	_delta=delta;
}
T_uniform::~T_uniform(){}

T_gaussian::T_gaussian(){}
T_gaussian::T_gaussian(double delta){
	_delta=delta;
}
T_gaussian::~T_gaussian(){}

//acceptance
posizione T_gaussian::eval(posizione x_old){    //generate new point starting from x_old with gaussian transition probability
	//cout<<_f<<endl;
	double x=_rnd.Gauss(x_old.getX(), _delta);

	posizione x_new(x);
	return x_new;
}

posizione T_uniform::eval(posizione x_old){    //generate new point starting from x_old with uniform transition probability

	double x=x_old.getX()+_rnd.Rannyu(-0.5,0.5)*_delta;
	posizione x_new=posizione(x);
	return x_new;
}

void Metropolis::SetProbabillity(probability *prob){
	_probability= prob;

}

Metropolis::Metropolis(int passi, probability *prob, acceptance *T){   //initializing 
	_passi=passi;
	_rnd.Initialize();
	_probability= prob;
	_acceptance=T;

}

Metropolis::~Metropolis(){}

void Metropolis::setX_start(posizione x){
	_x_start=x;
}


void Metropolis::Find_delta(){
	int N=10000;
	double perc=0;
	posizione x_new;
	posizione x_start=posizione();
	for(int i=0;i<N;i++){
		x_new=this->Acceptance(x_start,_acceptance->eval(x_start));
	
		if( (x_start.getX()!=x_new.getX())){
			perc+=1;
		}
		x_start=x_new;
	}
	cout<<"controllo che la percentuale di acceptance dell'algoritmo Metropolis sia circa 50%"<<endl;
	cout<<"percentuale(%)="<<perc/N*100<<endl;
	cout<<endl;
	return;
}


posizione Metropolis::Acceptance(posizione x_old, posizione x_new){
	double r=_rnd.Rannyu();
	if(1<_probability->eval(x_new)/_probability->eval(x_old)){
		return x_new;
	}
	else{
		if(_probability->eval(x_new)/_probability->eval(x_old)>r){
			return x_new;
		}
		else{
			return x_old;
		}
	}
}

posizione Metropolis::MR_alghoritm(){
	posizione x_new(0);
	//double x=x_start.getX()+_rnd.Rannyu(-0.5,0.5)*5;
	//posizione x_test=posizione(x); 
	x_new=this->Acceptance(_x_start,_acceptance->eval(_x_start));
	_x_start=x_new;

	return x_new;
}
