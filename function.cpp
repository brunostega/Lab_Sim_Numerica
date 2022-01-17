
#include "function.h"
#include <iostream>
#include <cmath>


using namespace std;



/////////////////////////FUNZIONEBASE///////////////////


////////////////////////SENO/////////////////////////////


double seno::eval(double x) const{

return sin(x);
}


double coseno::eval(double x) const{

return cos(x);
}

////////////////////////PARABOLA////////////////////////

parabola::parabola(double a, double b, double c){
        m_a=a;
        m_b=b;
        m_c=c;

}

parabola::~parabola(){}

double parabola::getA(){
        return m_a;
}

double parabola::getB(){
        return m_b;
}

double parabola::getC(){
        return m_c;
}

double parabola::eval(double x) const{

return x*x*m_a+x*m_b+m_c;
}



/////////////////////////GAUSSIANA////////////////////////////////////                                                                                                                        
gaussiana::gaussiana(double mean, double sigma){
        _mean=mean;
        _sigma=sigma;
}
gaussiana::~gaussiana(){}
double gaussiana::getMean(){return _mean;}
double gaussiana::getSigma(){return _sigma;}
double gaussiana::eval(double x) const{
        return (1./(sqrt(2*M_PI)*_sigma))*exp(-(pow((x-_mean),2)/(2*pow(_sigma,2)))) ;
}



