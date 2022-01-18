#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "../../random.h"
#include "../../block_stat.h"
#include "../../function.h"
#include "lib.h"
#include <string>
#include <cmath>


using namespace std;


int main (int argc, char *argv[]){

	//INITIALIZING SOME VARIABLES                                                   

    Random rnd;          //random generator     
    rnd.Initialize();    //initializing random generator

	int M=50000;         //# of points sampled to evaluate the integral 

	int N_blocks=100;         //# of blocks
    int N_throws=10000;       //# of throws
    int L=N_throws/N_blocks;  //# of throws in each block
	BlockStat block_stat=BlockStat(N_blocks, N_throws);     //initializing class for block statistics
	

    double ave[N_blocks]{};   //array of mean in each block
    double av2[N_blocks]{};   //array of mean squared in each block


//UNIFORM SAMPLING
	//evaluating integral of "function" sampled by a uniform distribution

	//EVALUATE AVERAGE IN EACH BLOCK
    cout<<"Evaluating the integral with uniform sampling"<<endl;

	for(int i=0; i<N_blocks; i++){
		double sum=0.;
    	for(int j=0;j<L;j++){ 
			double appo=0.;
        	for(int i=0;i<M;i++){

				appo+=function(rnd.Rannyu());     //function point generated with a uniform distribution
		
			}         
			appo=appo/M;         //evaluating the mean of the function 
			sum+=appo;           //cumulative sum in each block
		}
		ave[i]=sum/L;            //estimating mean for each block
		av2[i]=pow(ave[i],2);    //estimating mean squared for each block
	}

	//BLOCK STATISTICS
	cout<<"Doing block statistics"<<endl;

	block_stat.progressive_mean(ave);           //progressive mean
	block_stat.progressive_mean_squared(av2);   //progressive mean squared
	block_stat.progressive_error();             //progressive standard deviation

	//PRINTING OUPUT
	const string  file_out="I_unif_samp.csv";   //output file
	cout<<"Printing output in file "<<file_out<<endl<<endl;
	block_stat.print_output(file_out);          //printing the outpu: progressive mean and progressive error 


//IMPORTANCE SAMPLING (RETTA)
	//evaluating integral of "function" using importance sampling with a probability distribution given by 
	//2*(1-x)    ("retta")

	//EVALUATE AVERAGE IN EACH BLOCK
    cout<<"Evaluating the integral with importance sampling"<<endl;
	

	for(int i=0; i<N_blocks; i++){
		double sum=0;
    	for(int j=0;j<L;j++){ 
			double x_retta=0.;
			double integ_new=0.;
        	for(int i=0;i<M;i++){
				x_retta=cumulative_retta(rnd.Rannyu());        //generate random variable distributed as "retta" starting from a uniform distribution
				integ_new+=function(x_retta)/retta(x_retta);   //sum of the new function to integrate (starting function)/(new probability distribution)
			}         
			integ_new/=M;       //evaluating the mean of the function
			sum+=integ_new;     //cumulative sum in each block
       	}
		ave[i]=sum/L;           //estimating mean for each block
		av2[i]=pow(ave[i],2);   //estimating mean squared for each block
    }

	//BLOCK STATISTICS
	cout<<"Doing block statistics"<<endl;

	block_stat.progressive_mean(ave);           //progressive mean
	block_stat.progressive_mean_squared(av2);   //progressive mean squared
	block_stat.progressive_error();             //progressive standard deviation

	//PRINTING OUPUT
	const string  file_out2="I_imp_samp.csv";    //output file
	cout<<"Printing output in file "<<file_out2<<endl<<endl;
	block_stat.print_output(file_out2);          //printing the outpu: progressive mean and progressive error 


//IMPORTANCE SAMPLING ACCEPT REJECT

	//EVALUATE AVERAGE IN EACH BLOCK
    cout<<"Evaluating the integral with accept reject"<<endl;
	
	double a=-3./2;   //definisco la seconda funzione di importance sampling (parabola) ax^2+bx+c
	double b=0.;
	double c=3./2;
	parabola * imp_samp2=new parabola(a,b,c);


	for(int i=0; i<N_blocks; i++){
		double sum=0;
    	for(int j=0;j<L;j++){ 
			double x_acc_rej=0.;
			double integ_new=0.;
        	for(int i=0;i<M;i++){
				x_acc_rej=rnd.acceptReject(imp_samp2, 2., 0., 1.);          //generate point distributed with parabola probability function using accept reject method
				integ_new+=function(x_acc_rej)/imp_samp2->eval(x_acc_rej);  //sum of the new function to integrate (starting function)/(new probability distribution)
			}         
			integ_new/=M;         //evaluating the mean of the function
			sum+=integ_new;       //cumulative sum in each block
       	}
		ave[i]=sum/L;             //estimating mean for each block
		av2[i]=pow(ave[i],2);     //estimating mean squared for each block
	}

	//BLOCK STATISTICS
	cout<<"Doing block statistics"<<endl;

	block_stat.progressive_mean(ave);           //progressive mean
	block_stat.progressive_mean_squared(av2);   //progressive mean squared
	block_stat.progressive_error();             //progressive standard deviation

	//PRINTING OUPUT
	const string  file_out3="I_imp_samp2.csv";   //output file
	cout<<"Printing output in file "<<file_out3<<endl<<endl;
	block_stat.print_output(file_out3);          //printing the outpu: progressive mean and progressive error 

return 0;
}


