#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "lib.h"
#include <string>

using namespace std;


////////////////////////BLACKSHOLES////////////////////////////////////
BlackScholes::BlackScholes(double T, double K, double r, double mu, double sigma, double S0, Random & rnd){
        _T=T;
        _K=K;
        _r=r;
        _mu=mu;
        _sigma=sigma;
        _S0=S0;
        _rnd=rnd;
}

BlackScholes::~BlackScholes(){}
double BlackScholes::Eur_call_direct(double t){    //estimation of call price
        double s=this->S(t);
        return exp(-t*_r)*max(0.,s-_K);
        //return s * this->N( this->d1(t,s) ) - _K * exp(-_r * (_T-t) ) * this->N( this->d2(t,s) );
}
double BlackScholes::Eur_put_direct(double t){     //estimation of put price
        double s=this->S(t);
        return exp(-t*_r)*max(0.,_K-s);
        //return s * (this->N( this->d1(t,s) ) - 1 ) - _K * exp(-_r * (_T-t) ) * (this->N( this->d2(t,s) ) -1 );
}
double BlackScholes::Eur_call_discrete(double t, int N){
        double s=this->S_discrete(t,N);
        return exp(-t*_r)*max(0.,s-_K);
        //return s * this->N( this->d1(t, s) ) - _K * exp(-_r * (_T-t) ) * this->N( this->d2(t, s) );
}
double BlackScholes::Eur_put_discrete(double t, int N){
        double s=this->S_discrete(t,N);
        return exp(-t*_r)*max(0.,_K-s);
        //return s * (this->N( this->d1(t, s) ) - 1 ) - _K * exp(-_r * (_T-t) ) * (this->N( this->d2(t, s) ) -1 );
}
double BlackScholes::W(double t){     
        return _rnd.Gauss(0, sqrt(t) );
}
double BlackScholes::S(double t){
        return _S0 * exp( ( _mu-(0.5 * _sigma*_sigma ) ) * t + _sigma * this->W(t) );

}
double BlackScholes::S_discrete(double t, int N){
        double appo=_S0;
        double dt=t/N;
        for(int i=0;i<N; i++){
                appo=appo * exp( ( _mu-0.5 * _sigma*_sigma ) * dt + _sigma * this->W(1)*sqrt(dt) );
        }
        return appo;
}

/*
double BlackScholes::d1(double t,double s){
        return 1./( _sigma * sqrt(_T-t) )*( log(s / _K) + _r + _sigma*_sigma/2 * (_T-t) );
}
double BlackScholes::d2(double t,double s){
        return this->d1(t,s)-_sigma*sqrt(_T-t);
}
double BlackScholes::N(double x){
        return 0.5*( 1 + erf(x/sqrt(2)) );
}


*/





