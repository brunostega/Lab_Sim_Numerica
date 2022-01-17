#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "/home/bruno/Scrivania/universita/magistrale1/lab_simulaz/esercizi/random.h"
#include "/home/bruno/Scrivania/universita/magistrale1/lab_simulaz/esercizi/block_stat.h"
#include "lib.h"
//#include "vettore.h"
#include <string>
#include <cmath>


using namespace std;


int main (int argc, char *argv[]){

	//INITIALIZING SOME VARIABLES                                                   

    Random rnd;                 //random generator     
    rnd.Initialize();           //initialize random generator

	int N=100;                  //# of rw steps
	double passo=1;             //step lenght of rw
	int dim=3;                	//dimensionality of random walk
	double *abs_distance=new double[N];	//array of mean absolute distance as a function of number of steps
	double *err=new double[N];  //error as a function of number of steps

	int N_throws=10000;         //# of diffent "throws"
	int N_blocks=100;           //#of blocks
	int L=N_throws/N_blocks;    //# of differnt "throws" per block
	BlockStat block_stat=BlockStat(N_blocks, N_throws);     //initializing class for block statistics


//DISCRETE RANDOM WALK
	//estimating (with block statistics) the distance of a discrete random walk as a function of the number of steps
	//for every step generates N_throws random walks and estimate the mean squared distance and thenerror with block statistics


	double ave[N_blocks]{};   //array of mean in each block
	double av2[N_blocks]{};   //array of mean squared in each block

	//cicle over the number of steps
	for(int i=0;i<N;i++){      

		//cicle over the number of blocks
		for(int j=0;j<N_blocks;j++){   
			double appo=0;

			//cicle over the number of "throws" per block
			for(int l=0;l<L;l++){   
				double x[dim]{};        //position array variable   
				rnd.RW_discrete(passo, x , dim, i+1);       //discrete "i" random walks on "dim" dimensional lattice of passo "passo"
				appo+=modulo_quadro(dim,x);		            //sum module squared distance of the random walk		
			}	
			ave[j]=appo/L;            //estimating mean for each block
			av2[j]=pow(ave[j],2);     //estimating mean squared for each block

		}


		//BLOCK STATISTICS
		block_stat.reset();
		block_stat.progressive_mean(ave);           //progressive mean
		block_stat.progressive_mean_squared(av2);   //progressive mean squared
		block_stat.progressive_error();             //progressive error         //for the first step this will obviusly give error related to a zero variance. it's not a problem

		abs_distance[i]=sqrt(block_stat.get_final_mean());          //estimating sqrt of the mean distance (mean using all blocks)      
		err[i]=1/(2*abs_distance[i])*block_stat.get_final_error();  //estimating error (error using all blocks) associated with mean distance (using error propagation)
		cout<<block_stat.get_final_error()<<endl;
	}

	//OUTPUT
	string file="discrete_RW.csv";
	cout<<"printing the output (array of distance as a function on the number of steps and associated error) in "<<file<<endl;
	ofstream output;
	output.open(file);
	if(output.is_open()){
        	for(int i=0;i<N;i++){
        	        output<<abs_distance[i]<<","<<err[i]<<endl;
        	}
	}
	else{
        	cerr<<"PROBEL, couldn't open the file"<<endl;
	}
	output.close();



//COUNTINUUM RANDOM WALK
	//estimating (with block statistics) the distance of a continuum 3D random walk as a function of the number of steps
	//for every step generates N_throws random walks and estimate the mean squared distance and thenerror with block statistics

	//cicle over the number of steps
	for(int i=0;i<N;i++){

		//cicle over the number of blocks
		for(int j=0;j<N_blocks;j++){
			double appo=0;

			//cicle over the number of "throws" per block
			for(int l=0;l<L;l++){
				double x[dim]{};          //position array variable   
				rnd.RW_continum3D(passo, x , dim, i+1);     //continuum random walksin 3D "passo"
				appo+=modulo_quadro(dim,x);		            //sum module squared distance of the random walk		
			}	
			ave[j]=appo/L;            //estimating mean for each block
			av2[j]=pow(ave[j],2);     //estimating mean squared for each block   	
		}

        //BLOCK STATISTICS

		block_stat.reset();
		block_stat.progressive_mean(ave);           //progressive mean
		block_stat.progressive_mean_squared(av2);   //progressive mean squared
		block_stat.progressive_error();             //progressive error

		abs_distance[i]=sqrt(block_stat.get_final_mean());          //estimating sqrt of the mean distance (mean using all blocks)      
		err[i]=1/(2*abs_distance[i])*block_stat.get_final_error();  //estimating error (error using all blocks) associated with mean distance (using error propagation)
	}

	//OUTPUT
	file="continum_RW.csv";
	cout<<"printing the output (array of distance as a function on the number of steps and associated error) in "<<file<<endl;
	output.open(file);
	if(output.is_open()){
        	for(int i=0;i<N;i++){
        	        output<<abs_distance[i]<<","<<err[i]<<endl;
        	}
	}
	else{
        	cerr<<"PROBEL, couldn't open the file"<<endl;
	}
	output.close();

return 0;
}
