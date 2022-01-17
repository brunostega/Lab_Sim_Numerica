
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
		pos=rnd.Rannyu(0.,d);    //genero la posizione del centro dell'ago
		x=rnd.Rannyu();   
		y=rnd.Rannyu();
		while(pow(x,2)+pow(y,2)>1){   //genero due punti all'interno di una sfera di raggio unitario (necessario affinchè la distribuzione dell'angolo sia uniforme su tutta la sfera)
			x=rnd.Rannyu();
			y=rnd.Rannyu();
		}
		costeta=x/pow( pow(x,2) + pow(y,2) , 0.5 );    //il coseno dell'angolo associato è dato da x normalizzato opportunamente
		if( (pos-costeta*L/2>0) & (pos+costeta*L/2<d) ){  //se la proiezione sull'asse contenente le barre non interseca le barre, aumento il conteggio
			count+=1;
		}
	}

	prob=float(M-count)/M;
return 2*L/(prob*d);
}



