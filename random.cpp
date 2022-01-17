#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "function.h"
#include <string>

using namespace std;

const string path="/home/bruno/Scrivania/universita/magistrale1/lab_simulaz/esercizi/";

Random :: Random(){}         //coustruttore (posso eventualmente dargli una costruzione predefinita ma generalemte dovrò all'interno del codice settare la classe)
 
Random :: ~Random(){}        //disreuttore (necessario quando nella classe ci sono vettori dinamici che vanno deallocati dalla memoria)

void Random :: SaveSeed(){   //scrive su un file output i valori l_i del private
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Exponential(double lambda){           //generates random variable with exponential distribution
	double y=Rannyu();
	return -1./lambda *log(1.-y);
}

double Random :: Cauchy(double mean, double gamma) {   //generates random variable with cauchy distribution
   double s=Rannyu();
   
   double x=tan(M_PI*(s-1./2.));
   return mean + x * gamma;
}

double Random :: Dado(int facce){                      //generates random variable with integer uniform distribution in [1,6]
   double s=Rannyu();
   return int(s*facce)+1;
}

double Random :: Rannyu(double min, double max){       //genera un numero casuale tra [min, max]
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){                         //genera un numero casuale tra [0,1]
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){  //Set parametri della classe (4 seed e p1,p2 nel nostro caso)
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

void Random :: Initialize(){     //inizializzo la variabile random

   int seed[4];                  // definisco un array di dimensione 4
   int p1, p2;
   ifstream Primes(path+"Primes");
   if (Primes.is_open()){        //apro il file Primes
      Primes >> p1 >> p2 ;       //associo i primi numeri di Primes a p1,p2
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input(path+"seed.in");     //apro il file seed.in
   string property;               //definisco una stringa di nome proprerty
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;                 //nel file la prima riga è una stringa
         if( property == "RANDOMSEED" ){    //se il nome della stringa è questa allora..
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];         //imput i 4 dati del file
            this->SetRandom(seed,p1,p2);                               //alloca set random
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   this->SaveSeed(); //salva i nuovi seed in un file output
}

double Random::acceptReject(funzioneBase *F, double max, double a, double b){   //Accept reject method for a function F with maximum max in interval [a,b]
        double x=this->Rannyu(a,b);       //generate a random point in [a,b] 
        double y=max*this->Rannyu();      //generate a random point in [0,max]

        while(y>=F->eval(x)){             //continues to generate until (x,y) is accepted 
                x=this->Rannyu(a,b);
                y=max*this->Rannyu();
        }
return x;
}


void Random :: RW_discrete(double passo, double *x, int dim ,int N_passi){    //Discrete Random Walk
                                                                              //makes "N_passi" step of lenght "passo" in a casual (possible) direction of the different dimensions (dim)
	for(int i=0;i<N_passi;i++){
		int comp=int(this->Rannyu(0,dim));    //component in wich make the step
		double appo=this->Rannyu();           //random number to decide to step forward or backwards
		if(appo<0.5){
			x[comp]+=-passo;		
		}
		else{
			x[comp]+=passo;
		}
	}
}


void Random :: RW_continum3D(double passo, double *x, int dim ,int N_passi){   //Continuum Random Walk 3D
                                                                               //makes "N_passi" step of lenght "passo" in a casual direction 
	if(dim!=3){   
		cout<<"dimensionalità sbagliata. Questo metodo va bene solo per 3 dim"<<endl;
	}
	else{
		for(int i=0;i<N_passi;i++){
			double theta=this->Rannyu(0,M_PI);     //generate a casual zenital angle (theta in [0,pi])
			double phi=this->Rannyu(0,2*M_PI);     //generate a casual azimutal angle (phi in [0,2*pi])
			x[2]+=passo*cos(theta);                //makes the step
			x[1]+=passo*sin(theta)*sin(phi);
			x[0]+=passo*sin(theta)*cos(phi);
		}
	}
}


int Random::UpDown(){
   double r=this->Rannyu();
   if(r<0.5){return -1;}
   else{return +1;}
}
