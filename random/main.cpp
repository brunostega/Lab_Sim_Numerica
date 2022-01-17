/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "random.h"

using namespace std;

 
int main (int argc, char *argv[]){

int M=10000;

   Random rnd;
   int seed[4];                  // definisco un array di dimensione 4
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){        //apro il file Primes
      Primes >> p1 >> p2 ;       //associo i primi numeri di Primes a p1,p2
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");     //apro il file seed.in
   string property;               //definisco una stringa di nome proprerty
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;                 //nel file la prima riga è una stringa
         if( property == "RANDOMSEED" ){    //se il nome della stringa è questa allora..
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];         //imput i 4 dati del file
            rnd.SetRandom(seed,p1,p2);                                 //alloca set random
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   for(int i=0; i<20; i++){
      cout << rnd.Rannyu() << endl;       //crea i numeri casuali in [0,1]
   }

   rnd.SaveSeed();                          //salva i nuovi seed in un file output


ofstream output;
output.open("random.txt");
if(output.is_open()){
	for(int i=0;i<M;i++){
		output<<rnd.Rannyu()<<endl;
	}
}
else{
	cerr<<"PROBEL, couldn't open the file"<<endl;
}
output.close();

   return 0;
}



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
