
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "lib.h"
#include "/home/bruno/Scrivania/universita/magistrale1/lab_simulaz/esercizi/random.h"
#include <string>

using namespace std;

void exp(double *Sn,int N, int n, double lambda, Random rnd){  //calculating the sum of random variables with exponential distribution
	                                                           //for n realizations with increasing number (N) of random variables 
															   
	double sum=0;
	for(int i=0; i<N;i++){
		for(int j=0;j<n;j++){
			sum+=rnd.Exponential(lambda);
		}
		Sn[i]=1./n*sum;
		sum=0;
	}	
}

void cauchy(double *Sn,int N, int n, double mean,double gamma, Random rnd){  //calculating the sum of random variables with cauchy distribution
																			 //for n realizations with increasing number (N) of random variables
	double sum=0;
	for(int i=0; i<N;i++){
		for(int j=0;j<n;j++){
			sum+=rnd.Cauchy(mean, gamma);
		}
		Sn[i]=1./n*sum;
		sum=0;
	}	
}

void dado(double *Sn,int N, int n, int facce ,Random rnd){   //calculating the sum of random variables with uniform distribution (dado a 6 facce)
															 //for n realizations with increasing number (N) of random variables
	double sum=0;  
	for(int i=0; i<N;i++){
		for(int j=0;j<n;j++){
			sum+=rnd.Dado(facce);
		}
		Sn[i]=1./n*sum;
		sum=0;
	}
}



