
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "lib.h"
#include "/home/bruno/Scrivania/universita/magistrale1/lab_simulaz/esercizi/random.h"
#include <string>

using namespace std;




double pi(int M,double L, double d, Random & rnd){    //stima di pigreco con esperimento di Buffon

	double pos=0;
	double x=0;
	double y=0;
	double X=0;
	int count=0;
	double prob=0;

	for(int i=0;i<M;i++){
		pos=rnd.Rannyu(0.,d);
		x=rnd.Rannyu();
		y=rnd.Rannyu();
		while(pow(x,2)+pow(y,2)>1){
			x=rnd.Rannyu();
			y=rnd.Rannyu();
		}
		X=x/pow( pow(x,2) + pow(y,2) , 0.5 );
		if( (pos-X*L/2>0) & (pos+X*L/2<d) ){
			count+=1;
		}
	}

	prob=float(M-count)/M;
return 2*L/(prob*d);
}



