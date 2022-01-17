#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <cmath>
#include "lib.h"
#include "../random.h"
#include "../block_stat.h"
#include "../function.h"


using namespace std;


int main (int argc, char *argv[]){

	//INITIALIZING SOME VARIABLES
    Random rnd;               
    rnd.Initialize();     //initialize random variabile

	double T=1;              //delivery time
	double K=100;            //strike price  
	double r=0.1;            //risk free intrest rate
	double mu=r;             //intrest rate (>=r)
	double sigma=0.25;       //volatility 
	double S0=100;           //asset price at t=0

	int N_throws=70000;      //# throws
	int N_blocks=100;         //# blocks
	int L=N_throws/N_blocks;  //# throws per block
	int N=100;                //# of temporal steps

    const string  E_C_direct="Eur_call_direct.csv";         //output file for direct European call price
    const string  E_P_direct="Eur_put_direct.csv";          //output file for direct European put price
    const string  E_C_discrete="Eur_call_discrete.csv";     //output file for direct European call price
    const string  E_P_discrete="Eur_put_discrete.csv";      //output file for direct European put price

	double ave_C_direct[N_blocks]={};      //arrays for average
	double ave_P_direct[N_blocks]={};
	double ave_C_discrete[N_blocks]={};
	double ave_P_discrete[N_blocks]={};

	BlockStat block_stat(N_blocks,N_throws);      //construct the class for block statistics
	BlackScholes bh(T,K,r,mu,sigma,S0,rnd);       //construct the class BlackScholes

	for(int i=0; i<N_blocks;i++){

		for(int j=0;j<L;j++){

			ave_C_direct[i]+=bh.Eur_call_direct(T);         
			ave_P_direct[i]+=bh.Eur_put_direct(T);
			ave_C_discrete[i]+=bh.Eur_call_discrete(T, N);
			ave_P_discrete[i]+=bh.Eur_put_discrete(T, N);

		}
		ave_C_direct[i]/=L;
		ave_P_direct[i]/=L;
		ave_C_discrete[i]/=L;
		ave_P_discrete[i]/=L;
		
	}

	block_stat.block_stat_all(E_C_direct, ave_C_direct);             //print statistical block progression of European direct call price in file 
	block_stat.block_stat_all(E_P_direct, ave_P_direct);             //print statistical block progression of European direct put price in file
	block_stat.block_stat_all(E_C_discrete, ave_C_discrete);         //print statistical block progression of European discrete call price in file 
	block_stat.block_stat_all(E_P_discrete, ave_P_discrete);         //print statistical block progression of European discrete put price in file

return 0;
}
