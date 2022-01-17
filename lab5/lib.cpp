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
	m_y=0;
	m_z=0;
}

posizione::posizione(double x, double y, double z){
	m_x=x;
	m_y=y;
	m_z=z;
}

posizione::~posizione(){}

double posizione::getX() const{
return m_x;
}

double posizione::getY() const{
return m_y;
}

double posizione::getZ() const{
return m_z;
}

double posizione::getR() const{            //coordinate sferiche
return (sqrt(pow(m_x,2)+pow(m_y,2)+pow(m_z,2)));
}

double posizione::getTheta() const{
return (acos(m_z/getR()));
}

double posizione::getPhi() const{
return (atan2(m_y,m_x));
}

double posizione::getRho() const{         //coordinate cilindriche
return (sqrt(m_x*m_x+m_y*m_y));
}

double posizione::getDistanza( const posizione & b) const{
return (sqrt(pow(getX()-b.getX(),2)+pow(getY()-b.getY(),2)+pow(getZ()-b.getZ(),2)));
}




//HYDROGEN//////////////////////////////////////

double Hydrogen100::eval(posizione x){
	double H100_wave_func=pow( a_0, -1.5 ) / sqrt( M_PI ) * exp( -x.getR() / a_0 );
	return pow( H100_wave_func , 2 );
}

double Hydrogen210::eval(posizione x){
	double H210_wave_func=pow( a_0, -2.5 )/8 * sqrt(2/ M_PI ) * x.getR() * exp( -x.getR() / ( 2 * a_0 ) ) * cos( x.getTheta() );
	return pow( H210_wave_func , 2 );
}

double Hydrogen200::eval(posizione x){
	double H200_wave_func=pow( a_0, -1.5 )/4 * sqrt(1./( 2 * M_PI )) * ( 2 - x.getR()/a_0 ) * exp( -x.getR() / ( 2 * a_0 ) );
	return pow( H200_wave_func , 2 );
}




//METROPOLIS//////////////////////////////////

acceptance::acceptance(){
	Random rnd;
	rnd.Initialize();
	_rnd=rnd;
}
double acceptance::getDelta(){return _delta;} 
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
	double y=_rnd.Gauss(x_old.getY(), _delta);
	double z=_rnd.Gauss(x_old.getZ(), _delta);

	posizione x_new(x,y,z);
	return x_new;
}

posizione T_uniform::eval(posizione x_old){    //generate new point starting from x_old with uniform transition probability

	double x=x_old.getX()+_rnd.Rannyu(-1,1)*_delta;
	double y=x_old.getY()+_rnd.Rannyu(-1,1)*_delta;
	double z=x_old.getZ()+_rnd.Rannyu(-1,1)*_delta;

	posizione x_new=posizione(x,y,z);
	return x_new;
}


Metropolis::Metropolis(int N_equilib,probability *prob, acceptance *T, posizione x, string nome){   //initializing 
	_rnd.Initialize();
	_probability= prob;
	_acceptance=T;
	_xstart=x;
	_N_equilib=N_equilib;

	ofstream equilib_out;
	equilib_out.open("output/equilibration"+nome+".dat");
	for(int i=0;i<N_equilib;i++){
		equilib_out<<this->MR_alghoritm().getR()<<endl;
	}
	equilib_out.close();
}

Metropolis::~Metropolis(){}

double Metropolis::Find_delta(){
	int N=10000;
	double perc=0;
	posizione x_new;
	posizione x_start=posizione(0.01,0,0.003);
	for(int i=0;i<N;i++){
		x_new=this->Acceptance(x_start,_acceptance->eval(x_start));
		//cout<<_acceptance->eval(x_start).getR()<<endl;
		if( (x_start.getX()!=x_new.getX()) & (x_start.getY()!=x_new.getY()) & (x_start.getZ()!=x_new.getZ()) ){
			perc+=1;
		}
		x_start=x_new;
	}
	cout<<"controllo che la percentuale di acceptance dell'algoritmo Metropolis sia circa 50%"<<endl;
	cout<<"percentuale(%)="<<perc/N*100<<endl;
	cout<<endl;

	return perc/N ;
}


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

posizione Metropolis::MR_alghoritm(){
	posizione x_new(0,0,0);  
	x_new=this->Acceptance(_xstart,_acceptance->eval(_xstart));
	_xstart=x_new;
	return x_new;
}
