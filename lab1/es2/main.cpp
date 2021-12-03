#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "/home/bruno/Scrivania/universita/magistrale1/lab_simulaz/esercizi/random.h"
#include "lib.h"
#include <string>

using namespace std;

int main (int argc, char *argv[]){

//INITIALISING SOME VARIABLES                                                   

Random rnd;         //random generator     
rnd.Initialize();   //initializing random generator

int N=10000;        //# of realization
double mean=0;      //mean of the cauchy distribution
double lambda=1.;   //second parameter of the cauchy distribution
double gamma=1.;    //parameter for the exponential distribution
int facce=6;        //numero di facce del dado 

int N1=1;           //different number of throws
int N2=2;
int N10=10;
int N100=100;

double S1[N]{};     //arrays for different number of trhows
double S2[N]{};
double S10[N]{};
double S100[N]{};


//make all "experiments" for the different probability distribution to test the central limit theorem

//DADO A 6 FACCE
cout<<"testing central limit theorem for a dice"<<endl;

dado(S1,N,N1,facce, rnd);
dado(S2,N,N2,facce, rnd);
dado(S10,N,N10,facce, rnd);
dado(S100,N,N100,facce, rnd);

//PRINTING OUTPUT
const string file_out_dado="dado.csv";
cout<<"printing the result in "<<file_out_dado<<endl<<endl;
ofstream output;
output.open(file_out_dado);
if(output.is_open()){
        for(int i=0;i<N;i++){
                output<<S1[i]<<","<<S2[i]<<","<<S10[i]<<","<<S100[i]<<endl;
        }
}
else{
        cerr<<"PROBEL, couldn't open the file"<<endl;
}
output.close();


//EXPONENTIAL DISTRIBUTION
cout<<"testing central limit theorem for a exponential distribution"<<endl;

exp(S1,N,N1,lambda, rnd);
exp(S2,N,N2,lambda, rnd);
exp(S10,N,N10,lambda, rnd);
exp(S100,N,N100,lambda, rnd);

//PRINTING THE OUTPUT
const string file_out_exp="exp.csv";
cout<<"printing the result in "<<file_out_exp<<endl<<endl;
output.open(file_out_exp);
if(output.is_open()){
        for(int i=0;i<N;i++){
                output<<S1[i]<<","<<S2[i]<<","<<S10[i]<<","<<S100[i]<<endl;
        }
}
else{
        cerr<<"PROBEL, couldn't open the file"<<endl;
}
output.close();


//CAUCHY DISTRIBUTION
cout<<"testing central limit theorem for a cauchy distribution"<<endl;

cauchy(S1,N,N1,mean,gamma, rnd);
cauchy(S2,N,N2,mean,gamma, rnd);
cauchy(S10,N,N10,mean,gamma, rnd);
cauchy(S100,N,N100,mean,gamma, rnd);

//PRINTING OUTPUT
const string file_out_cauchy="cauchy.csv";
cout<<"printing the result in "<<file_out_cauchy<<endl<<endl;
output.open(file_out_cauchy);
if(output.is_open()){
        for(int i=0;i<N;i++){
                output<<S1[i]<<","<<S2[i]<<","<<S10[i]<<","<<S100[i]<<endl;
        }
}
else{
        cerr<<"PROBEL, couldn't open the file"<<endl;
}
output.close();


return 0;
}

